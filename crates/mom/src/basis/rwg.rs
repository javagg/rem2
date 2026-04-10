//! RWG (Rao-Wilton-Glisson) vector basis functions.
//!
//! Reference: S.M. Rao, D.R. Wilton, A.W. Glisson,
//! "Electromagnetic scattering by surfaces of arbitrary shape",
//! IEEE TAP 30(3) 1982, 409-418.

use crate::surface_mesh::{SurfaceMesh, SharedEdge};

/// One RWG basis function, defined on a shared interior edge.
#[derive(Debug, Clone)]
pub struct RwgBasis {
    /// Index into `SurfaceMesh::edges`
    pub edge_idx: usize,
    /// T⁺ face index
    pub plus_face: usize,
    /// T⁻ face index
    pub minus_face: usize,
    /// Free vertex of T⁺ (the vertex not on the shared edge)
    pub free_node_plus: usize,
    /// Free vertex of T⁻
    pub free_node_minus: usize,
    /// Edge length lₙ [m]
    pub length: f64,
}

impl RwgBasis {
    /// Evaluate f_n(r) on T⁺ (`in_plus = true`) or T⁻ (`in_plus = false`).
    ///
    /// ```text
    /// f_n(r) = ±lₙ/(2Aᵢ) * (r - rᵢ_free)
    /// ```
    pub fn eval(&self, r: &[f64; 3], surf: &SurfaceMesh, in_plus: bool) -> [f64; 3] {
        let (face_idx, free_node, sign) = if in_plus {
            (self.plus_face, self.free_node_plus, 1.0_f64)
        } else {
            (self.minus_face, self.free_node_minus, -1.0_f64)
        };
        let area = surf.faces[face_idx].area;
        let free = &surf.nodes[free_node];
        let scale = sign * self.length / (2.0 * area);
        [
            scale * (r[0] - free[0]),
            scale * (r[1] - free[1]),
            scale * (r[2] - free[2]),
        ]
    }

    /// Surface divergence: ∇_s · f_n = ±lₙ/Aᵢ (constant on each triangle).
    pub fn divergence(&self, surf: &SurfaceMesh, in_plus: bool) -> f64 {
        let (face_idx, sign) = if in_plus {
            (self.plus_face, 1.0_f64)
        } else {
            (self.minus_face, -1.0_f64)
        };
        sign * self.length / surf.faces[face_idx].area
    }
}

/// Generate one `RwgBasis` for each shared interior edge.
pub fn generate_rwg_bases(surf: &SurfaceMesh) -> Vec<RwgBasis> {
    surf.edges.iter().enumerate().map(|(ei, edge)| {
        let (free_plus, free_minus) = find_free_nodes(edge, surf);
        RwgBasis {
            edge_idx:        ei,
            plus_face:       edge.plus_face,
            minus_face:      edge.minus_face,
            free_node_plus:  free_plus,
            free_node_minus: free_minus,
            length:          edge.length,
        }
    }).collect()
}

/// Find the free vertex (not on the shared edge) in each adjacent triangle.
fn find_free_nodes(edge: &SharedEdge, surf: &SurfaceMesh) -> (usize, usize) {
    let en = &edge.nodes;
    let free = |fi: usize| -> usize {
        surf.faces[fi].nodes.iter()
            .copied()
            .find(|&n| n != en[0] && n != en[1])
            .expect("triangle must have a vertex not on the shared edge")
    };
    (free(edge.plus_face), free(edge.minus_face))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::surface_mesh::SurfaceMesh;

    /// Build a minimal two-triangle surface mesh for testing.
    fn two_tri_mesh() -> SurfaceMesh {
        // Vertices:  0=(0,0,0)  1=(1,0,0)  2=(0,1,0)  3=(1,1,0)
        // T0: [0,1,2],  T1: [1,3,2]  — share edge (1,2)
        use rem_mesh::RemMesh;
        // We can't build RemMesh easily here; instead exercise via SurfaceMesh directly.
        // Build manually.
        use crate::surface_mesh::{TriFace, SharedEdge, tri_geometry};

        let nodes = vec![
            [0.0,0.0,0.0],
            [1.0,0.0,0.0],
            [0.0,1.0,0.0],
            [1.0,1.0,0.0],
        ];
        let (c0,n0,a0) = tri_geometry(&nodes[0], &nodes[1], &nodes[2]);
        let (c1,n1,a1) = tri_geometry(&nodes[1], &nodes[3], &nodes[2]);
        let faces = vec![
            TriFace { nodes:[0,1,2], centroid:c0, normal:n0, area:a0 },
            TriFace { nodes:[1,3,2], centroid:c1, normal:n1, area:a1 },
        ];
        let edges = vec![SharedEdge {
            nodes: [1,2],
            plus_face: 0,
            minus_face: 1,
            length: (1.0f64 + 1.0).sqrt(),
        }];
        SurfaceMesh { nodes, faces, edges, boundary_edges: vec![], face_attrs: vec![0, 0] }
    }

    #[test]
    fn rwg_divergence_signs() {
        let surf = two_tri_mesh();
        let bases = generate_rwg_bases(&surf);
        assert_eq!(bases.len(), 1);
        let b = &bases[0];
        // Divergence on T+ must be positive, T- negative
        assert!(b.divergence(&surf, true)  > 0.0);
        assert!(b.divergence(&surf, false) < 0.0);
    }

    #[test]
    fn rwg_normal_continuity() {
        // The normal (to the shared edge) component of f_n must be continuous
        // across the edge. Evaluate at the edge midpoint.
        let surf = two_tri_mesh();
        let bases = generate_rwg_bases(&surf);
        let b = &bases[0];
        let edge = &surf.edges[b.edge_idx];
        let mid = [
            (surf.nodes[edge.nodes[0]][0] + surf.nodes[edge.nodes[1]][0]) / 2.0,
            (surf.nodes[edge.nodes[0]][1] + surf.nodes[edge.nodes[1]][1]) / 2.0,
            (surf.nodes[edge.nodes[0]][2] + surf.nodes[edge.nodes[1]][2]) / 2.0,
        ];
        let fp = b.eval(&mid, &surf, true);
        let fm = b.eval(&mid, &surf, false);
        // Edge direction vector
        let e_dir = {
            let n0 = &surf.nodes[edge.nodes[0]];
            let n1 = &surf.nodes[edge.nodes[1]];
            let len = edge.length;
            [(n1[0]-n0[0])/len, (n1[1]-n0[1])/len, (n1[2]-n0[2])/len]
        };
        // Face normal (same for both since they're coplanar here)
        let face_normal = surf.faces[b.plus_face].normal;
        // Edge outward normal = face_normal × e_dir
        let en_normal = cross(&face_normal, &e_dir);
        let dot_p = dot(&fp, &en_normal);
        let dot_m = dot(&fm, &en_normal);
        assert!((dot_p - dot_m).abs() < 1e-10,
            "Normal continuity failed: f⁺·n = {}, f⁻·n = {}", dot_p, dot_m);
    }

    fn cross(a: &[f64;3], b: &[f64;3]) -> [f64;3] {
        [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
    }
    fn dot(a: &[f64;3], b: &[f64;3]) -> f64 {
        a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
    }
}
