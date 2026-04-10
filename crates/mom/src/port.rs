//! MoM lumped-port model for S-parameter extraction.
//!
//! A lumped port is defined over a set of surface triangles (identified by
//! boundary attribute IDs).  The port produces:
//! - An excitation RHS vector (one solve per active port).
//! - V/I extraction from the solved surface current coefficients.

use crate::surface_mesh::SurfaceMesh;
use crate::basis::rwg::RwgBasis;
use num_complex::Complex64;
use rem_core::{RemResult, RemError};

/// Lumped port: a set of RWG indices + excitation direction + reference Z₀.
#[derive(Debug, Clone)]
pub struct MomLumpedPort {
    /// 1-based port index (matches Touchstone ordering).
    pub index: u32,
    /// Indices into the RWG basis array that are on this port's surface.
    pub rwg_indices: Vec<usize>,
    /// Dominant field direction unit vector [x, y, z].
    pub direction: [f64; 3],
    /// Reference impedance Z₀ [Ω].
    pub z0: f64,
}

impl MomLumpedPort {
    /// Find the RWG basis functions whose shared edge lies on a face that
    /// belongs to one of `port_attrs`.  A face "belongs to a port" if its
    /// attribute tag (stored in `SurfaceMesh::face_attrs`) matches.
    ///
    /// `direction_str` — "x", "y", or "z"; mapped to unit vector.
    pub fn from_surface(
        surf: &SurfaceMesh,
        bases: &[RwgBasis],
        port_attrs: &[u32],
        index: u32,
        direction_str: &str,
        z0: f64,
    ) -> RemResult<Self> {
        let port_attr_set: std::collections::HashSet<u32> =
            port_attrs.iter().copied().collect();

        let rwg_indices: Vec<usize> = bases.iter().enumerate()
            .filter(|(_, b)| {
                let plus_attr  = surf.face_attrs.get(b.plus_face).copied().unwrap_or(0);
                let minus_attr = surf.face_attrs.get(b.minus_face).copied().unwrap_or(0);
                port_attr_set.contains(&plus_attr) && port_attr_set.contains(&minus_attr)
            })
            .map(|(i, _)| i)
            .collect();

        Ok(Self {
            index,
            rwg_indices,
            direction: direction_vec(direction_str),
            z0,
        })
    }

    /// Build the N-element excitation RHS for this port (one entry per RWG).
    ///
    /// For RWG basis, the port contribution at basis m is:
    ///   V_m = -∫_{port} f_m(r) · d̂ dS   (d̂ = direction unit vector)
    ///
    /// Only the `rwg_indices` of this port get non-zero entries; all others are 0.
    /// The excitation amplitude `v0` is nominally 1 V.
    pub fn excitation_rhs(
        &self,
        surf: &SurfaceMesh,
        bases: &[RwgBasis],
        n_total: usize,
        v0: Complex64,
    ) -> Vec<Complex64> {
        let mut rhs = vec![Complex64::ZERO; n_total];
        let d = self.direction;
        for &mi in &self.rwg_indices {
            let b = &bases[mi];
            let mut val = Complex64::ZERO;
            // Integrate f_m · d̂ over both support triangles using centroid quadrature
            for &(face_idx, in_plus) in &[(b.plus_face, true), (b.minus_face, false)] {
                let face = &surf.faces[face_idx];
                let fm = b.eval(&face.centroid, surf, in_plus);
                let dot = d[0]*fm[0] + d[1]*fm[1] + d[2]*fm[2];
                val += Complex64::new(dot * face.area, 0.0);
            }
            rhs[mi] = -v0 * val;
        }
        rhs
    }

    /// Extract port current I from solved surface current coefficients.
    ///
    /// I ≈ Σ_{m ∈ port} a_m * (div_m_plus * A_plus + div_m_minus * A_minus)
    ///
    /// where a_m are the solved RWG coefficients.
    pub fn extract_current(
        &self,
        surf: &SurfaceMesh,
        bases: &[RwgBasis],
        coeffs: &[Complex64],
    ) -> Complex64 {
        let mut i_port = Complex64::ZERO;
        for &mi in &self.rwg_indices {
            if mi >= coeffs.len() { continue; }
            let b = &bases[mi];
            let div_p = b.divergence(surf, true);
            let div_m = b.divergence(surf, false);
            let area_p = surf.faces[b.plus_face].area;
            let area_m = surf.faces[b.minus_face].area;
            let contrib = div_p * area_p + div_m * area_m;
            i_port += coeffs[mi] * contrib;
        }
        i_port
    }
}

fn direction_vec(s: &str) -> [f64; 3] {
    match s.to_lowercase().as_str() {
        "y" => [0.0, 1.0, 0.0],
        "z" => [0.0, 0.0, 1.0],
        _   => [1.0, 0.0, 0.0],  // default "x"
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::surface_mesh::{SurfaceMesh, TriFace, SharedEdge, tri_geometry};
    use crate::basis::rwg::generate_rwg_bases;

    /// Build a two-triangle surface with attribute tags.
    fn two_tri_port_surf() -> SurfaceMesh {
        let nodes = vec![
            [0.0_f64, 0.0, 0.0],
            [1.0,     0.0, 0.0],
            [0.5,     1.0, 0.0],
            [-0.5,    1.0, 0.0],
        ];
        let (c0, n0, a0) = tri_geometry(&nodes[0], &nodes[1], &nodes[2]);
        let (c1, n1, a1) = tri_geometry(&nodes[0], &nodes[2], &nodes[3]);
        let faces = vec![
            TriFace { nodes: [0,1,2], centroid: c0, normal: n0, area: a0 },
            TriFace { nodes: [0,2,3], centroid: c1, normal: n1, area: a1 },
        ];
        let edges = vec![SharedEdge {
            nodes: [0, 2],
            plus_face: 0,
            minus_face: 1,
            length: (nodes[2][0].powi(2) + nodes[2][1].powi(2)).sqrt(),
        }];
        SurfaceMesh {
            nodes,
            faces,
            edges,
            boundary_edges: vec![[0,1],[1,2],[2,3],[3,0]],
            face_attrs: vec![1, 1],   // both faces tagged as attr 1
        }
    }

    #[test]
    fn from_surface_finds_rwg_on_attr() {
        let surf = two_tri_port_surf();
        let bases = generate_rwg_bases(&surf);
        assert_eq!(bases.len(), 1, "should have 1 shared edge");
        let port = MomLumpedPort::from_surface(
            &surf, &bases, &[1], 1, "x", 50.0
        ).expect("port construction failed");
        assert_eq!(port.rwg_indices.len(), 1, "all RWG on attr-1 surface should be found");
        assert_eq!(port.rwg_indices[0], 0);
        assert_eq!(port.z0, 50.0);
    }

    #[test]
    fn from_surface_empty_when_no_match() {
        let surf = two_tri_port_surf();
        let bases = generate_rwg_bases(&surf);
        let port = MomLumpedPort::from_surface(
            &surf, &bases, &[99], 1, "x", 50.0
        ).expect("should succeed even with no match");
        assert!(port.rwg_indices.is_empty());
    }

    #[test]
    fn excitation_rhs_zero_outside_port() {
        let surf = two_tri_port_surf();
        let bases = generate_rwg_bases(&surf);
        let port = MomLumpedPort::from_surface(&surf, &bases, &[1], 1, "x", 50.0).unwrap();
        let rhs = port.excitation_rhs(&surf, &bases, bases.len(), Complex64::new(1.0, 0.0));
        assert_eq!(rhs.len(), bases.len());
        // All entries on the port should be non-zero (there's 1 RWG on port attr 1)
        assert!(rhs[0].norm() > 0.0, "port RWG should have non-zero excitation");
    }

    #[test]
    fn extract_current_finite() {
        let surf = two_tri_port_surf();
        let bases = generate_rwg_bases(&surf);
        let port = MomLumpedPort::from_surface(&surf, &bases, &[1], 1, "x", 50.0).unwrap();
        let coeffs = vec![Complex64::new(1.0, 0.5); bases.len()];
        let i = port.extract_current(&surf, &bases, &coeffs);
        assert!(i.re.is_finite() && i.im.is_finite());
    }
}
