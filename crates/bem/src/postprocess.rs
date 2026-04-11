//! BEM post-processing: capacitance, surface potential, far field.

use rem_mom::surface_mesh::SurfaceMesh;
use rem_core::EPS0;

/// Compute the total charge Q = ∫_S σ(r) dS from P0 coefficients.
pub fn total_charge(sigma: &[f64], surf: &SurfaceMesh) -> f64 {
    sigma.iter().zip(surf.faces.iter())
        .map(|(&s, f)| s * f.area)
        .sum()
}

/// Compute capacitance C = Q / V from total charge and applied potential.
///
/// For a grounded sphere at potential V₀, C = Q * ε₀ / V₀.
pub fn capacitance(sigma: &[f64], surf: &SurfaceMesh, v0: f64) -> f64 {
    EPS0 * total_charge(sigma, surf) / v0
}

/// Evaluate the potential at an external observation point r_obs.
///
/// φ(r_obs) = ∫_S G(r_obs, r') σ(r') dS'
///          ≈ Σ_n G(r_obs, centroid_n) * σ_n * area_n
pub fn eval_potential(
    r_obs: &[f64; 3],
    sigma: &[f64],
    surf: &SurfaceMesh,
) -> f64 {
    use crate::kernel::laplace_G;
    sigma.iter().zip(surf.faces.iter())
        .map(|(&s, f)| laplace_G(r_obs, &f.centroid) * s * f.area)
        .sum()
}

#[cfg(test)]
mod tests {
    use super::*;
    use rem_mom::surface_mesh::{SurfaceMesh, TriFace, SharedEdge, tri_geometry, patch_edge_lengths};

    fn single_tri_mesh() -> SurfaceMesh {
        let nodes = vec![[0.0,0.0,0.0],[1.0,0.0,0.0],[0.0,1.0,0.0]];
        let (c,n,a) = tri_geometry(&nodes[0],&nodes[1],&nodes[2]);
        SurfaceMesh {
            nodes,
            faces: vec![TriFace { nodes:[0,1,2], centroid:c, normal:n, area:a }],
            edges: vec![],
            boundary_edges: vec![[0,1],[1,2],[2,0]],
            face_attrs: vec![0],
            global_node_ids: vec![],
        }
    }

    #[test]
    fn total_charge_unit_sigma() {
        let surf = single_tri_mesh();
        let sigma = vec![1.0];
        let q = total_charge(&sigma, &surf);
        assert!((q - surf.faces[0].area).abs() < 1e-14);
    }

    #[test]
    fn eval_potential_positive_at_distance() {
        let surf = single_tri_mesh();
        let sigma = vec![1.0];
        let r_far = [100.0, 0.0, 0.0];
        let phi = eval_potential(&r_far, &sigma, &surf);
        assert!(phi > 0.0, "Potential should be positive for positive charge");
        assert!(phi.is_finite());
    }
}
