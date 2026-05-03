//! MoM adaptive mesh refinement (AMR).
//!
//! Implements a surface-current-based error indicator, Dörfler marking,
//! and uniform midpoint-subdivision refinement for triangular surface meshes.
//!
//! # Workflow
//! ```text
//! solve MoM → currents  →  AmrtIndicator::compute(surf, bases, currents)
//!                                ↓
//!                         dorfer_mark(indicators, θ=0.5)  →  Vec<bool>
//!                                ↓
//!                         refine_surface(surf, marked)  →  SurfaceMesh
//! ```

use crate::basis::rwg::{RwgBasis, generate_rwg_bases};
use crate::surface_mesh::{SurfaceMesh, TriFace, SharedEdge, tri_geometry, patch_edge_lengths, build_edge_topology};
use num_complex::Complex64;
use std::collections::HashMap;

// ─── Error indicator ─────────────────────────────────────────────────────────

/// Per-face surface-current error indicator.
///
/// η_m = √(Σ_{n∈m} |I_n|² · l_n²) × √(A_m)
///
/// where the sum runs over all RWG basis functions n whose shared edge
/// belongs to face m (either T⁺ or T⁻).
pub struct AmrtIndicator {
    /// One indicator per face in the surface mesh.
    pub face_errors: Vec<f64>,
}

impl AmrtIndicator {
    /// Compute the error indicator from solved RWG current coefficients.
    pub fn compute(
        surf: &SurfaceMesh,
        bases: &[RwgBasis],
        currents: &[Complex64],
    ) -> Self {
        let n_faces = surf.faces.len();
        // Sum |I_n|² * l_n² into each adjacent face
        let mut face_sum = vec![0.0_f64; n_faces];
        for (n, basis) in bases.iter().enumerate() {
            let amp = currents[n].norm_sqr() * basis.length * basis.length;
            face_sum[basis.plus_face]  += amp;
            face_sum[basis.minus_face] += amp;
        }
        let face_errors = surf.faces.iter().zip(face_sum.iter())
            .map(|(f, &s)| s.sqrt() * f.area.sqrt())
            .collect();
        Self { face_errors }
    }
}

// ─── Dörfler marking ─────────────────────────────────────────────────────────

/// Dörfler (bulk) marking strategy.
///
/// Marks the minimal subset of faces such that the sum of their squared error
/// indicators is ≥ `theta` × (total sum of squared indicators).
///
/// `theta` should be in (0, 1]; a typical value is 0.5.
pub fn dorfer_mark(indicators: &[f64], theta: f64) -> Vec<bool> {
    let total: f64 = indicators.iter().map(|e| e * e).sum();
    if total == 0.0 {
        return vec![false; indicators.len()];
    }
    // Sort face indices by error² descending
    let mut order: Vec<usize> = (0..indicators.len()).collect();
    order.sort_unstable_by(|&a, &b| {
        indicators[b].partial_cmp(&indicators[a]).unwrap_or(std::cmp::Ordering::Equal)
    });

    let target = theta * total;
    let mut accumulated = 0.0;
    let mut marked = vec![false; indicators.len()];
    for idx in order {
        if accumulated >= target {
            break;
        }
        marked[idx] = true;
        accumulated += indicators[idx] * indicators[idx];
    }
    marked
}

// ─── Mesh refinement ─────────────────────────────────────────────────────────

/// Refine marked faces by 1→4 uniform midpoint subdivision.
///
/// Each marked face is split into four sub-triangles by inserting midpoints on
/// all three edges.  Unmarked faces are kept as-is.  The resulting mesh
/// preserves the attribute tags of the original faces.
pub fn refine_surface(surf: &SurfaceMesh, marked: &[bool]) -> SurfaceMesh {
    assert_eq!(marked.len(), surf.faces.len());

    let mut new_nodes: Vec<[f64; 3]> = surf.nodes.clone();
    let mut new_face_nodes: Vec<[usize; 3]> = Vec::new();
    let mut new_face_attrs: Vec<u32> = Vec::new();
    let mut new_global_ids: Vec<usize> = surf.global_node_ids.clone();

    // Midpoint cache: edge (min, max) → new node index
    let mut midpoint_cache: HashMap<(usize, usize), usize> = HashMap::new();

    let mut get_or_create_midpoint =
        |n0: usize, n1: usize,
         nodes: &mut Vec<[f64; 3]>,
         gids: &mut Vec<usize>,
         cache: &mut HashMap<(usize, usize), usize>| -> usize {
            let key = (n0.min(n1), n0.max(n1));
            if let Some(&idx) = cache.get(&key) {
                return idx;
            }
            let p0 = nodes[n0];
            let p1 = nodes[n1];
            let mid = [
                (p0[0] + p1[0]) * 0.5,
                (p0[1] + p1[1]) * 0.5,
                (p0[2] + p1[2]) * 0.5,
            ];
            let idx = nodes.len();
            nodes.push(mid);
            // global id: 0 (new interior node, no parent mesh node)
            gids.push(0);
            cache.insert(key, idx);
            idx
        };

    for (fi, face) in surf.faces.iter().enumerate() {
        let [n0, n1, n2] = face.nodes;
        let attr = surf.face_attrs.get(fi).copied().unwrap_or(0);

        if !marked[fi] {
            new_face_nodes.push([n0, n1, n2]);
            new_face_attrs.push(attr);
        } else {
            // Insert edge midpoints
            let m01 = get_or_create_midpoint(n0, n1, &mut new_nodes, &mut new_global_ids, &mut midpoint_cache);
            let m12 = get_or_create_midpoint(n1, n2, &mut new_nodes, &mut new_global_ids, &mut midpoint_cache);
            let m20 = get_or_create_midpoint(n2, n0, &mut new_nodes, &mut new_global_ids, &mut midpoint_cache);
            // 4 child triangles
            new_face_nodes.push([n0,  m01, m20]);
            new_face_nodes.push([m01, n1,  m12]);
            new_face_nodes.push([m20, m12, n2 ]);
            new_face_nodes.push([m01, m12, m20]);
            for _ in 0..4 { new_face_attrs.push(attr); }
        }
    }

    // Rebuild TriFace array and topology
    let new_faces: Vec<TriFace> = new_face_nodes.iter().map(|&ns| {
        let (centroid, normal, area) = tri_geometry(
            &new_nodes[ns[0]], &new_nodes[ns[1]], &new_nodes[ns[2]],
        );
        TriFace { nodes: ns, centroid, normal, area }
    }).collect();

    SurfaceMesh::from_parts(new_nodes, new_faces, new_face_attrs, new_global_ids)
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::surface_mesh::SurfaceMesh;
    use crate::basis::rwg::generate_rwg_bases;
    use num_complex::Complex64;

    /// Build a flat 2-triangle surface (one shared edge) at z=0.
    fn two_tri_surf() -> SurfaceMesh {
        // Two triangles sharing edge [1,2]:
        //   T+: nodes 0,1,2
        //   T-: nodes 1,3,2
        let nodes = vec![
            [0.0_f64, 0.0, 0.0],  // 0
            [1.0,     0.0, 0.0],  // 1
            [0.5,     1.0, 0.0],  // 2
            [-0.5,    1.0, 0.0],  // 3
        ];
        let (c0, n0, a0) = tri_geometry(&nodes[0], &nodes[1], &nodes[2]);
        let (c1, n1, a1) = tri_geometry(&nodes[1], &nodes[3], &nodes[2]);
        let faces = vec![
            TriFace { nodes: [0, 1, 2], centroid: c0, normal: n0, area: a0 },
            TriFace { nodes: [1, 3, 2], centroid: c1, normal: n1, area: a1 },
        ];
        SurfaceMesh::from_parts(nodes, faces, vec![1, 1], vec![])
    }

    #[test]
    fn error_indicator_nonzero() {
        let surf = two_tri_surf();
        let bases = generate_rwg_bases(&surf);
        assert_eq!(bases.len(), 1, "expected 1 RWG basis");
        let currents = vec![Complex64::new(1.0, 0.5)];
        let ind = AmrtIndicator::compute(&surf, &bases, &currents);
        assert_eq!(ind.face_errors.len(), 2);
        assert!(ind.face_errors[0] > 0.0);
        assert!(ind.face_errors[1] > 0.0);
    }

    #[test]
    fn dorfer_mark_marks_highest() {
        let errors = vec![0.1, 0.9, 0.3, 0.05];
        let marked = dorfer_mark(&errors, 0.5);
        // Total² = 0.01+0.81+0.09+0.0025 = 0.9125
        // Target  = 0.45625; largest is 0.81 → face 1 marked → 0.81 >= 0.45625
        assert!(marked[1], "face 1 (highest error) must be marked");
        // Face 0 should NOT be marked (0.81 alone covers >50%)
        assert!(!marked[0]);
        assert!(!marked[2]);
    }

    #[test]
    fn refine_marks_all_makes_four_times_faces() {
        let surf = two_tri_surf();
        let marked = vec![true, true];
        let refined = refine_surface(&surf, &marked);
        assert_eq!(refined.faces.len(), 8, "2 faces × 4 = 8 after full refinement");
        // Node count: 4 original + 3 new midpoints (one per edge of each face
        // but T+ and T- share edge [1,2] → only 1 midpoint there)
        // Edges: 3 (T+: 01,12,20) + 2 new unique (T- has 13,32; shares 12 midpoint)
        // Actually midpoints on 01, 12, 20, 13, 32 = 5 unique edges → 5 midpoints
        // but edge 12 is shared → 4+4+1 = wait let me count:
        // T+ edges: (0,1), (1,2), (2,0)  = 3 new midpoints
        // T- edges: (1,3), (3,2), (2,1)  = (2,1) == (1,2) shared → 2 new midpoints
        // Total new nodes = 3 + 2 = 5, total nodes = 9
        assert!(refined.nodes.len() >= 7, "need at least 7 nodes after refinement");
    }

    #[test]
    fn refine_unmarked_face_unchanged() {
        let surf = two_tri_surf();
        let marked = vec![true, false];
        let refined = refine_surface(&surf, &marked);
        // T+ is refined (→4), T- is kept (→1): total 5 faces
        assert_eq!(refined.faces.len(), 5);
    }
}
