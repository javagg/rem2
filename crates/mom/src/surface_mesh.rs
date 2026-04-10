//! Surface mesh extraction and edge-face topology for MoM/BEM.
//!
//! Extracts the triangular surface mesh from `RemMesh` boundary elements
//! and builds the shared-edge data structure needed for RWG basis functions.

use rem_core::{RemError, RemResult};
use rem_mesh::{RemMesh, ElementKind, Element};
use std::collections::HashMap;

// ---------------------------------------------------------------------------
// Data structures
// ---------------------------------------------------------------------------

/// A single triangular surface face.
#[derive(Debug, Clone)]
pub struct TriFace {
    /// Global node indices (into `SurfaceMesh::nodes`)
    pub nodes: [usize; 3],
    /// Face centroid [m]
    pub centroid: [f64; 3],
    /// Outward unit normal
    pub normal: [f64; 3],
    /// Face area [m²]
    pub area: f64,
}

/// A shared interior edge shared by exactly two triangles (carrier of one RWG basis function).
#[derive(Debug, Clone)]
pub struct SharedEdge {
    /// Sorted global node indices of this edge
    pub nodes: [usize; 2],
    /// Index of T⁺ face in `SurfaceMesh::faces`
    pub plus_face: usize,
    /// Index of T⁻ face in `SurfaceMesh::faces`
    pub minus_face: usize,
    /// Edge length [m]
    pub length: f64,
}

/// Extracted triangular surface mesh with edge-face topology.
pub struct SurfaceMesh {
    /// Node coordinates [m]
    pub nodes: Vec<[f64; 3]>,
    /// All triangular faces
    pub faces: Vec<TriFace>,
    /// Interior shared edges (one RWG basis function per edge)
    pub edges: Vec<SharedEdge>,
    /// Boundary (open) edges belonging to exactly one face
    pub boundary_edges: Vec<[usize; 2]>,
    /// Per-face physical-group attribute tag (from the mesh boundary tags).
    /// Zero means "untagged / not extracted from a named boundary".
    pub face_attrs: Vec<u32>,
    /// Global node IDs in the parent RemMesh (index i → global node ID).
    /// Allows FE-BI coupling to map surface local DOFs back to volume DOFs.
    pub global_node_ids: Vec<usize>,
}

// ---------------------------------------------------------------------------
// Implementation
// ---------------------------------------------------------------------------

impl SurfaceMesh {
    /// Extract triangular faces from `RemMesh` whose physical tag is in `pec_attrs`.
    ///
    /// Only `Tri3` elements are supported at this stage; `Tri6` nodes are silently
    /// reduced to their first three (corner) nodes.
    ///
    /// For surface-only meshes (no volume elements, e.g. SBR+ sphere meshes) the
    /// Tri3 elements appear in `volume_elements` rather than `boundary_elements`.
    /// This function searches `boundary_elements` first, then falls back to
    /// `volume_elements` so both use-cases work transparently.
    pub fn extract(rem_mesh: &RemMesh, pec_attrs: &[u32]) -> RemResult<Self> {
        let attr_set: std::collections::HashSet<u32> = pec_attrs.iter().copied().collect();

        let is_tri = |e: &&Element| {
            attr_set.contains(&e.tag) &&
            matches!(e.kind, ElementKind::Tri3 | ElementKind::Tri6)
        };

        // ── 1. Filter triangular elements (boundary first, then volume fallback) ─
        let tri_elems: Vec<_> = {
            let from_boundary: Vec<_> = rem_mesh.boundary_elements.iter().filter(is_tri).collect();
            if !from_boundary.is_empty() {
                from_boundary
            } else {
                // Surface-only mesh: Tri3 are in volume_elements
                rem_mesh.volume_elements.iter().filter(is_tri).collect()
            }
        };

        if tri_elems.is_empty() {
            return Err(RemError::Mesh(format!(
                "No triangular boundary elements found for PEC attributes {:?}", pec_attrs
            )));
        }

        // ── 2. First pass: collect all node IDs used by the triangles ────────
        let mut global_ids: Vec<usize> = tri_elems.iter()
            .flat_map(|e| e.node_ids[..3].iter().copied())
            .collect();
        global_ids.sort_unstable();
        global_ids.dedup();

        // Build global→local map and node list
        let global_to_local: std::collections::HashMap<usize, usize> = global_ids.iter()
            .enumerate()
            .map(|(li, &gi)| (gi, li))
            .collect();

        let nodes: Vec<[f64; 3]> = global_ids.iter()
            .map(|&gi| {
                let n = &rem_mesh.nodes[gi];
                [n.x, n.y, n.z]
            })
            .collect();

        // ── 3. Build faces ───────────────────────────────────────────────────
        let mut faces: Vec<TriFace> = Vec::with_capacity(tri_elems.len());
        let mut face_attrs: Vec<u32> = Vec::with_capacity(tri_elems.len());
        for elem in &tri_elems {
            let l0 = global_to_local[&elem.node_ids[0]];
            let l1 = global_to_local[&elem.node_ids[1]];
            let l2 = global_to_local[&elem.node_ids[2]];

            let p0 = nodes[l0];
            let p1 = nodes[l1];
            let p2 = nodes[l2];

            let (centroid, normal, area) = tri_geometry(&p0, &p1, &p2);
            faces.push(TriFace { nodes: [l0, l1, l2], centroid, normal, area });
            face_attrs.push(elem.tag as u32);
        }

        // ── 4. Build edge-face topology ──────────────────────────────────────
        let (mut edges, boundary_edges) = build_edge_topology(&faces);
        patch_edge_lengths(&mut edges, &nodes);

        Ok(SurfaceMesh { nodes, faces, edges, boundary_edges, face_attrs, global_node_ids: global_ids })
    }

    /// Number of RWG basis functions = number of shared interior edges.
    pub fn n_rwg(&self) -> usize { self.edges.len() }

    /// Number of pulse basis functions = number of faces.
    pub fn n_pulse(&self) -> usize { self.faces.len() }
}

// ---------------------------------------------------------------------------
// Geometry helpers
// ---------------------------------------------------------------------------

/// Compute centroid, outward unit normal, and area of a triangle.
pub fn tri_geometry(p0: &[f64; 3], p1: &[f64; 3], p2: &[f64; 3])
    -> ([f64; 3], [f64; 3], f64)
{
    let e1 = [p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2]];
    let e2 = [p2[0]-p0[0], p2[1]-p0[1], p2[2]-p0[2]];

    // Cross product e1 × e2
    let cx = e1[1]*e2[2] - e1[2]*e2[1];
    let cy = e1[2]*e2[0] - e1[0]*e2[2];
    let cz = e1[0]*e2[1] - e1[1]*e2[0];
    let cross_len = (cx*cx + cy*cy + cz*cz).sqrt();

    let area = 0.5 * cross_len;
    let normal = if cross_len > 1e-30 {
        [cx / cross_len, cy / cross_len, cz / cross_len]
    } else {
        [0.0, 0.0, 1.0]   // degenerate fallback
    };

    let centroid = [
        (p0[0] + p1[0] + p2[0]) / 3.0,
        (p0[1] + p1[1] + p2[1]) / 3.0,
        (p0[2] + p1[2] + p2[2]) / 3.0,
    ];

    (centroid, normal, area)
}

/// Compute edge length between two nodes.
pub fn edge_length(nodes: &[[f64; 3]], n0: usize, n1: usize) -> f64 {
    let p0 = &nodes[n0];
    let p1 = &nodes[n1];
    let dx = p1[0]-p0[0];
    let dy = p1[1]-p0[1];
    let dz = p1[2]-p0[2];
    (dx*dx + dy*dy + dz*dz).sqrt()
}

// ---------------------------------------------------------------------------
// Edge topology
// ---------------------------------------------------------------------------

/// Build `SharedEdge` list and boundary edge list from face connectivity.
///
/// Algorithm: for every edge (pair of nodes), track which faces contain it.
/// Edges shared by exactly 2 faces → `SharedEdge`. Edges owned by 1 face → boundary.
fn build_edge_topology(faces: &[TriFace]) -> (Vec<SharedEdge>, Vec<[usize; 2]>) {
    // Map: sorted edge (n0, n1) → list of face indices
    let mut edge_map: HashMap<(usize, usize), Vec<usize>> = HashMap::new();

    for (fi, face) in faces.iter().enumerate() {
        let [a, b, c] = face.nodes;
        for &(u, v) in &[(a,b), (b,c), (c,a)] {
            let key = if u < v { (u, v) } else { (v, u) };
            edge_map.entry(key).or_default().push(fi);
        }
    }

    // Identify node coordinates length (faces reference indices, we need actual coords
    // but `TriFace` doesn't carry them — compute length lazily via `SurfaceMesh::nodes`).
    // We'll fill `length` after constructing the SurfaceMesh using a dummy 0.0 here;
    // the caller (SurfaceMesh::extract) will patch it in the returned value.
    // Actually: we don't have nodes here. We store the node pair and compute length later.

    let mut shared: Vec<SharedEdge> = Vec::new();
    let mut boundary: Vec<[usize; 2]> = Vec::new();

    for ((n0, n1), face_list) in &edge_map {
        match face_list.len() {
            1 => { boundary.push([*n0, *n1]); }
            2 => {
                shared.push(SharedEdge {
                    nodes: [*n0, *n1],
                    plus_face:  face_list[0],
                    minus_face: face_list[1],
                    length: 0.0,  // patched below
                });
            }
            n => {
                log::warn!("Edge ({},{}) shared by {} faces — non-manifold surface", n0, n1, n);
            }
        }
    }

    (shared, boundary)
}

/// Patch `SharedEdge::length` using the node coordinate array.
/// Called by `SurfaceMesh::extract` after both `nodes` and `edges` are built.
pub fn patch_edge_lengths(edges: &mut [SharedEdge], nodes: &[[f64; 3]]) {
    for e in edges.iter_mut() {
        e.length = edge_length(nodes, e.nodes[0], e.nodes[1]);
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn unit_tri() -> ([f64; 3], [f64; 3], [f64; 3]) {
        ([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0])
    }

    #[test]
    fn test_tri_geometry_unit() {
        let (p0, p1, p2) = unit_tri();
        let (centroid, normal, area) = tri_geometry(&p0, &p1, &p2);
        assert!((area - 0.5).abs() < 1e-14, "area = {}", area);
        assert!((normal[2] - 1.0).abs() < 1e-14, "normal = {:?}", normal);
        assert!((centroid[0] - 1.0/3.0).abs() < 1e-14);
    }

    #[test]
    fn test_edge_topology_two_triangles() {
        // Two triangles sharing edge (1,2)
        //   0---1
        //   |\ /|
        //   | 2 |
        //   |/ \|
        //   (not actually that but two triangles sharing an edge)
        let faces = vec![
            TriFace { nodes: [0,1,2], centroid: [0.0;3], normal: [0.0,0.0,1.0], area: 0.5 },
            TriFace { nodes: [1,2,3], centroid: [0.0;3], normal: [0.0,0.0,1.0], area: 0.5 },
        ];
        let (shared, boundary) = build_edge_topology(&faces);
        // Edge (1,2) is shared; edges (0,1),(0,2),(1,3),(2,3) are boundary
        assert_eq!(shared.len(), 1, "should have 1 shared edge");
        assert_eq!(boundary.len(), 4, "should have 4 boundary edges");
        let e = &shared[0];
        assert!(e.nodes == [1,2] || e.nodes == [1,2]);
    }
}
