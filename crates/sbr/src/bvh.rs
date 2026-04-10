//! AABB BVH (Bounding Volume Hierarchy) for fast ray–triangle intersection.
//!
//! Uses Surface Area Heuristic (SAH) splitting with a flat node array for
//! cache-friendly traversal. Compatible with WASM (no unsafe, no C FFI).

use crate::ray::{ray_triangle, RayHit, sub3, add3, scale3};
use rem_mom::surface_mesh::SurfaceMesh;
use std::sync::Arc;

// ---------------------------------------------------------------------------
// AABB
// ---------------------------------------------------------------------------

/// Axis-aligned bounding box.
#[derive(Debug, Clone, Copy)]
pub struct Aabb {
    pub min: [f64; 3],
    pub max: [f64; 3],
}

impl Aabb {
    pub fn empty() -> Self {
        Self {
            min: [f64::INFINITY;  3],
            max: [f64::NEG_INFINITY; 3],
        }
    }

    pub fn from_triangle(p0: &[f64; 3], p1: &[f64; 3], p2: &[f64; 3]) -> Self {
        Self {
            min: [
                p0[0].min(p1[0]).min(p2[0]),
                p0[1].min(p1[1]).min(p2[1]),
                p0[2].min(p1[2]).min(p2[2]),
            ],
            max: [
                p0[0].max(p1[0]).max(p2[0]),
                p0[1].max(p1[1]).max(p2[1]),
                p0[2].max(p1[2]).max(p2[2]),
            ],
        }
    }

    pub fn expand(&mut self, other: &Aabb) {
        for i in 0..3 {
            self.min[i] = self.min[i].min(other.min[i]);
            self.max[i] = self.max[i].max(other.max[i]);
        }
    }

    pub fn surface_area(&self) -> f64 {
        let d = sub3(&self.max, &self.min);
        2.0 * (d[0]*d[1] + d[1]*d[2] + d[0]*d[2])
    }

    pub fn centroid(&self) -> [f64; 3] {
        scale3(&add3(&self.min, &self.max), 0.5)
    }

    /// Slab method ray–AABB intersection. Returns minimum positive t, or None.
    pub fn intersect_ray(&self, origin: &[f64; 3], inv_dir: &[f64; 3]) -> Option<f64> {
        let mut t_min = f64::NEG_INFINITY;
        let mut t_max = f64::INFINITY;

        for i in 0..3 {
            let t0 = (self.min[i] - origin[i]) * inv_dir[i];
            let t1 = (self.max[i] - origin[i]) * inv_dir[i];
            let (lo, hi) = if t0 < t1 { (t0, t1) } else { (t1, t0) };
            t_min = t_min.max(lo);
            t_max = t_max.min(hi);
        }

        if t_max >= t_min.max(0.0) { Some(t_min) } else { None }
    }
}

// ---------------------------------------------------------------------------
// BVH node (flat array, index-based children)
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
enum BvhNode {
    Leaf {
        bounds: Aabb,
        face_start: u32,
        face_count: u32,
    },
    Interior {
        bounds: Aabb,
        left: u32,   // index into Bvh::nodes
        right: u32,
        #[allow(dead_code)]
        split_axis: u8,
    },
}

impl BvhNode {
    fn bounds(&self) -> &Aabb {
        match self {
            BvhNode::Leaf { bounds, .. } => bounds,
            BvhNode::Interior { bounds, .. } => bounds,
        }
    }
}

// ---------------------------------------------------------------------------
// BVH tree
// ---------------------------------------------------------------------------

const MAX_LEAF_FACES: u32 = 4;
const SAH_BUCKETS:    usize = 12;

/// BVH acceleration structure over a `SurfaceMesh`.
pub struct Bvh {
    nodes: Vec<BvhNode>,
    /// Permuted face indices (leaves reference contiguous slices)
    face_indices: Vec<usize>,
    surf: Arc<SurfaceMesh>,
}

impl Bvh {
    /// Build BVH from a surface mesh in O(N log N).
    pub fn build(surf: Arc<SurfaceMesh>) -> Self {
        let n = surf.faces.len();
        let mut face_indices: Vec<usize> = (0..n).collect();

        // Precompute per-face AABB and centroid
        let face_bounds: Vec<Aabb> = surf.faces.iter().map(|f| {
            let p0 = &surf.nodes[f.nodes[0]];
            let p1 = &surf.nodes[f.nodes[1]];
            let p2 = &surf.nodes[f.nodes[2]];
            Aabb::from_triangle(p0, p1, p2)
        }).collect();

        let mut nodes = Vec::with_capacity(2 * n);
        build_recursive(&face_bounds, &mut face_indices, 0, n, &mut nodes);

        Bvh { nodes, face_indices, surf }
    }

    /// Nearest hit query. Returns the closest `RayHit` along the ray, or `None`.
    pub fn intersect(&self, origin: &[f64; 3], dir: &[f64; 3]) -> Option<RayHit> {
        let inv_dir = [1.0/dir[0], 1.0/dir[1], 1.0/dir[2]];
        let mut best: Option<RayHit> = None;
        let mut stack = [0u32; 64];
        let mut sp = 0usize;
        stack[sp] = 0; sp += 1;

        while sp > 0 {
            sp -= 1;
            let node_idx = stack[sp] as usize;
            let node = &self.nodes[node_idx];

            let t_box = node.bounds().intersect_ray(origin, &inv_dir);
            let t_cur = best.as_ref().map_or(f64::INFINITY, |h| h.t);
            if t_box.map_or(true, |t| t >= t_cur) {
                continue;
            }

            match node {
                BvhNode::Leaf { face_start, face_count, .. } => {
                    let start = *face_start as usize;
                    let end   = start + *face_count as usize;
                    for &fi in &self.face_indices[start..end] {
                        let face = &self.surf.faces[fi];
                        let p0 = &self.surf.nodes[face.nodes[0]];
                        let p1 = &self.surf.nodes[face.nodes[1]];
                        let p2 = &self.surf.nodes[face.nodes[2]];
                        let t_min = best.as_ref().map_or(1e-6, |h| h.t);
                        if let Some((t, u, v)) = ray_triangle(origin, dir, p0, p1, p2, 1e-6) {
                            if best.as_ref().map_or(true, |h| t < h.t) && t > t_min {
                                let w = 1.0 - u - v;
                                let pt = [
                                    w*p0[0] + u*p1[0] + v*p2[0],
                                    w*p0[1] + u*p1[1] + v*p2[1],
                                    w*p0[2] + u*p1[2] + v*p2[2],
                                ];
                                best = Some(RayHit {
                                    t,
                                    face_idx: fi,
                                    bary: [u, v, w],
                                    point: pt,
                                    normal: face.normal,
                                });
                            }
                        }
                    }
                }
                BvhNode::Interior { left, right, .. } => {
                    // Push both children; closer one last (visited first)
                    if sp + 2 <= stack.len() {
                        stack[sp] = *right; sp += 1;
                        stack[sp] = *left;  sp += 1;
                    }
                }
            }
        }
        best
    }

    /// Shadow query: does *any* face block the ray within distance `max_t`?
    pub fn any_hit(&self, origin: &[f64; 3], dir: &[f64; 3], max_t: f64) -> bool {
        let inv_dir = [1.0/dir[0], 1.0/dir[1], 1.0/dir[2]];
        let mut stack = [0u32; 64];
        let mut sp = 0usize;
        stack[sp] = 0; sp += 1;

        while sp > 0 {
            sp -= 1;
            let node_idx = stack[sp] as usize;
            let node = &self.nodes[node_idx];

            if node.bounds().intersect_ray(origin, &inv_dir)
                .map_or(true, |t| t >= max_t)
            {
                continue;
            }

            match node {
                BvhNode::Leaf { face_start, face_count, .. } => {
                    let start = *face_start as usize;
                    let end   = start + *face_count as usize;
                    for &fi in &self.face_indices[start..end] {
                        let face = &self.surf.faces[fi];
                        let p0 = &self.surf.nodes[face.nodes[0]];
                        let p1 = &self.surf.nodes[face.nodes[1]];
                        let p2 = &self.surf.nodes[face.nodes[2]];
                        if let Some((t, ..)) = ray_triangle(origin, dir, p0, p1, p2, 1e-6) {
                            if t < max_t { return true; }
                        }
                    }
                }
                BvhNode::Interior { left, right, .. } => {
                    if sp + 2 <= stack.len() {
                        stack[sp] = *right; sp += 1;
                        stack[sp] = *left;  sp += 1;
                    }
                }
            }
        }
        false
    }

    pub fn surf(&self) -> &SurfaceMesh { &self.surf }
}

// ---------------------------------------------------------------------------
// Recursive SAH build
// ---------------------------------------------------------------------------

/// Returns the index of the newly created node in `nodes`.
fn build_recursive(
    face_bounds: &[Aabb],
    face_indices: &mut [usize],
    start: usize,
    end: usize,
    nodes: &mut Vec<BvhNode>,
) -> u32 {
    let count = end - start;
    let node_idx = nodes.len() as u32;
    nodes.push(BvhNode::Leaf { bounds: Aabb::empty(), face_start: 0, face_count: 0 }); // placeholder

    // Compute bounding box of all faces in [start, end)
    let mut bounds = Aabb::empty();
    for &fi in &face_indices[start..end] {
        bounds.expand(&face_bounds[fi]);
    }

    if count <= MAX_LEAF_FACES as usize {
        nodes[node_idx as usize] = BvhNode::Leaf {
            bounds,
            face_start: start as u32,
            face_count: count as u32,
        };
        return node_idx;
    }

    // SAH: find best split across 3 axes
    let (split_axis, split_pos) = sah_split(face_bounds, &face_indices[start..end], &bounds);

    // Partition face_indices in-place around split_pos
    let mid = partition_faces(face_bounds, &mut face_indices[start..end], split_axis, split_pos) + start;

    // Degenerate partition: force equal halves
    let mid = if mid == start || mid == end { (start + end) / 2 } else { mid };

    let left  = build_recursive(face_bounds, face_indices, start, mid, nodes);
    let right = build_recursive(face_bounds, face_indices, mid,   end, nodes);

    nodes[node_idx as usize] = BvhNode::Interior {
        bounds,
        left,
        right,
        split_axis: split_axis as u8,
    };
    node_idx
}

/// SAH bucket split: returns (axis, split centroid value).
fn sah_split(face_bounds: &[Aabb], faces: &[usize], node_bounds: &Aabb) -> (usize, f64) {
    let node_sa = node_bounds.surface_area().max(1e-30);
    let mut best_cost = f64::INFINITY;
    let mut best_axis = 0usize;
    let mut best_split = 0.0f64;

    for axis in 0..3usize {
        let lo = node_bounds.min[axis];
        let hi = node_bounds.max[axis];
        if (hi - lo).abs() < 1e-14 { continue; }

        // Build SAH buckets
        let mut bucket_bounds = [Aabb::empty(); SAH_BUCKETS];
        let mut bucket_count  = [0u32; SAH_BUCKETS];

        for &fi in faces {
            let c = face_bounds[fi].centroid()[axis];
            let b = ((c - lo) / (hi - lo) * SAH_BUCKETS as f64) as usize;
            let b = b.min(SAH_BUCKETS - 1);
            bucket_bounds[b].expand(&face_bounds[fi]);
            bucket_count[b] += 1;
        }

        // Evaluate cost for each split plane
        for split in 1..SAH_BUCKETS {
            let mut left_b = Aabb::empty(); let mut left_n = 0u32;
            let mut right_b = Aabb::empty(); let mut right_n = 0u32;
            for b in 0..split      { left_b.expand(&bucket_bounds[b]);  left_n  += bucket_count[b]; }
            for b in split..SAH_BUCKETS { right_b.expand(&bucket_bounds[b]); right_n += bucket_count[b]; }

            let cost = 0.125 + (left_n as f64 * left_b.surface_area()
                + right_n as f64 * right_b.surface_area()) / node_sa;

            if cost < best_cost {
                best_cost  = cost;
                best_axis  = axis;
                best_split = lo + (split as f64 / SAH_BUCKETS as f64) * (hi - lo);
            }
        }
    }
    (best_axis, best_split)
}

/// Partition `faces` in-place so faces with centroid[axis] < split come first.
/// Returns the index (relative to slice start) of the first right-side face.
fn partition_faces(face_bounds: &[Aabb], faces: &mut [usize], axis: usize, split: f64) -> usize {
    let mut lo = 0;
    let mut hi = faces.len();
    while lo < hi {
        if face_bounds[faces[lo]].centroid()[axis] < split {
            lo += 1;
        } else {
            hi -= 1;
            faces.swap(lo, hi);
        }
    }
    lo
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn single_tri_bvh() -> Bvh {
        use rem_mom::surface_mesh::{SurfaceMesh, TriFace};
        use rem_mom::surface_mesh::tri_geometry;

        let nodes = vec![[0.0f64,0.0,0.0],[1.0,0.0,0.0],[0.0,1.0,0.0]];
        let (c, n, a) = tri_geometry(&nodes[0], &nodes[1], &nodes[2]);
        let surf = SurfaceMesh {
            nodes,
            faces: vec![TriFace { nodes:[0,1,2], centroid:c, normal:n, area:a }],
            edges: vec![],
            boundary_edges: vec![],
            face_attrs: vec![0],
        };
        Bvh::build(Arc::new(surf))
    }

    #[test]
    fn bvh_hit_single_triangle() {
        let bvh = single_tri_bvh();
        let hit = bvh.intersect(&[0.1, 0.1, 2.0], &[0.0, 0.0, -1.0]);
        assert!(hit.is_some(), "should hit triangle");
        assert!((hit.unwrap().t - 2.0).abs() < 1e-10);
    }

    #[test]
    fn bvh_miss_outside() {
        let bvh = single_tri_bvh();
        let hit = bvh.intersect(&[2.0, 2.0, 1.0], &[0.0, 0.0, -1.0]);
        assert!(hit.is_none());
    }

    #[test]
    fn bvh_any_hit() {
        let bvh = single_tri_bvh();
        assert!(bvh.any_hit(&[0.1, 0.1, 2.0], &[0.0, 0.0, -1.0], 10.0));
        assert!(!bvh.any_hit(&[2.0, 2.0, 2.0], &[0.0, 0.0, -1.0], 10.0));
    }

    #[test]
    fn aabb_ray_hit() {
        let b = Aabb { min: [-1.0, -1.0, -1.0], max: [1.0, 1.0, 1.0] };
        let o = [0.0, 0.0, 5.0];
        let inv = [f64::INFINITY, f64::INFINITY, -1.0];
        assert!(b.intersect_ray(&o, &inv).is_some());
    }
}
