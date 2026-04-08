//! Physical Theory of Diffraction (PTD) edge correction for SBR+.
//!
//! Implements the Ufimtsev / Mitzner line-integral edge correction that
//! supplements the PO surface currents with fringe-wave contributions from
//! the open (silhouette) edges of the illuminated surface.  PTD corrections
//! are most significant at large bistatic angles, shadow-boundaries, and for
//! small-kl edges (kl ~ 1–10).
//!
//! # Algorithm (scalar Ufimtsev fringe current)
//!
//! For each boundary edge of length L with mid-point r_e, unit tangent t̂,
//! inward-surface-normal n̂, and outward edge-normal m̂ = t̂ × n̂:
//!
//! 1. Geometric illumination test:  `cos_α = dot(-k̂_inc, t̂)` is the grazing
//!    angle cosine.  Skip edge if it is not in the illuminated shadow boundary.
//!
//! 2. Local TE/TM decomposition of the incident field at r_e:
//!    - E_te = (E_inc · t̂) t̂   (parallel to edge)
//!    - E_tm = E_inc − E_te      (perpendicular to edge)
//!
//! 3. Ufimtsev fringe diffraction coefficients (hard/soft wedge, exterior
//!    angle β_e = π − interior wedge angle; we use β_e = π for flat faces):
//!
//!    D_s(ϕ, ϕ') = −1/(2√(2πk)) · [
//!        1/(sin ϕ · sin ϕ') · (cot((π + Δ⁻)/2) + cot((π − Δ⁻)/2))
//!    ]
//!    where Δ± = ϕ ± ϕ'  and ϕ, ϕ' are the observation/incidence angles in
//!    the edge-fixed coordinate system.
//!
//!    We use the simplified scalar PTD fringe coefficient for a half-plane
//!    (β_e = π):
//!    f_TE = (E_te · t̂) · (−1/2) · g(ϕ, ϕ')   [fringe current density]
//!    f_TM = (E_tm · m̂) · (−1/2) · g(ϕ, ϕ')
//!
//!    where g(ϕ, ϕ') = the half-plane edge diffraction factor below.
//!
//! 4. Far-field integration:
//!    Each edge of length L contributes to the far-field pattern via a
//!    line-integral of the fringe current weighted by exp(+jk r̂·r_e):
//!    ΔN = f_edge · L · exp(jk r̂·r_e)  (discrete mid-point approximation)
//!
//! # Reference
//! Ufimtsev, P.Ya., "Fundamentals of the Physical Theory of Diffraction",
//! IEEE Press, 2007, Chapter 3.

use num_complex::Complex64;
use rem_mom::surface_mesh::SurfaceMesh;
use crate::excitation::PlaneWave;
use crate::ray::{dot3, cross3, normalize3};

const MIN_EDGE_LEN: f64 = 1e-15;

/// A single boundary edge with pre-computed geometry.
#[derive(Debug, Clone)]
pub struct BoundaryEdge {
    /// Mid-point [m]
    pub midpoint: [f64; 3],
    /// Unit tangent (from node0 toward node1)
    pub tangent: [f64; 3],
    /// Edge length [m]
    pub length: f64,
    /// Outward edge normal m̂ = t̂ × n̂_face  (points away from material)
    pub edge_normal: [f64; 3],
    /// Inward surface normal of the adjacent face
    pub face_normal: [f64; 3],
}

/// Extract boundary edges from the SurfaceMesh.
///
/// `boundary_edges` are edges that belong to exactly one face (open boundary).
/// We annotate each with its adjacent face normal so we can compute m̂.
pub fn extract_boundary_edges(surf: &SurfaceMesh) -> Vec<BoundaryEdge> {
    // For each boundary edge [n0, n1] find the adjacent face
    let mut out = Vec::with_capacity(surf.boundary_edges.len());

    'edge: for &[n0, n1] in &surf.boundary_edges {
        let p0 = &surf.nodes[n0];
        let p1 = &surf.nodes[n1];
        let edge_vec = [p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2]];
        let length = (edge_vec[0]*edge_vec[0]+edge_vec[1]*edge_vec[1]+edge_vec[2]*edge_vec[2]).sqrt();
        if length < MIN_EDGE_LEN { continue; }
        let tangent = [edge_vec[0]/length, edge_vec[1]/length, edge_vec[2]/length];
        let midpoint = [
            (p0[0]+p1[0])*0.5,
            (p0[1]+p1[1])*0.5,
            (p0[2]+p1[2])*0.5,
        ];

        // Find adjacent face (contains both n0 and n1)
        for face in &surf.faces {
            let ns = face.nodes;
            if (ns[0]==n0 || ns[1]==n0 || ns[2]==n0) &&
               (ns[0]==n1 || ns[1]==n1 || ns[2]==n1) {
                let n_hat = face.normal;
                // Edge normal = t̂ × n̂ (points outward from material)
                let edge_normal = normalize3(cross3(&tangent, &n_hat));
                out.push(BoundaryEdge {
                    midpoint,
                    tangent,
                    length,
                    edge_normal,
                    face_normal: n_hat,
                });
                continue 'edge;
            }
        }
        // If no adjacent face found, still add with a zero normal (will be skipped)
        out.push(BoundaryEdge {
            midpoint,
            tangent,
            length,
            edge_normal: [0.0; 3],
            face_normal: [0.0; 3],
        });
    }

    out
}

/// Ufimtsev scalar half-plane diffraction coefficient for observation angle ϕ
/// and incidence angle ϕ' (measured from the edge's shadow boundary, in [0, 2π]).
///
/// Uses the PTD fringe function for the PEC half-plane:
///   D(ϕ, ϕ') = D_PO(ϕ, ϕ') − D_exact_Sommerfeld(ϕ, ϕ')  [fringe correction]
///
/// For the simplified mid-edge formulation we use the leading-order half-plane
/// fringe coefficient (Ufimtsev eq. 3.31, β = π):
///   g(ϕ, ϕ') = 1 / (cos(ϕ/2 - ϕ'/2)) − 1 / (cos(ϕ/2 + ϕ'/2))
///
/// (This is the non-uniform fringe current coefficient, exact at ϕ = ϕ'.)
fn half_plane_fringe(phi_obs: f64, phi_inc: f64) -> f64 {
    let half_diff = 0.5 * (phi_obs - phi_inc);
    let half_sum  = 0.5 * (phi_obs + phi_inc);

    let cos_d = half_diff.cos();
    let cos_s = half_sum.cos();

    // Avoid division by zero at shadow/reflection boundaries
    let term1 = if cos_d.abs() > 1e-8 { 1.0 / cos_d } else { 0.0 };
    let term2 = if cos_s.abs() > 1e-8 { 1.0 / cos_s } else { 0.0 };

    term1 - term2
}

/// Compute the PTD fringe-current far-field contribution N_ptd for a single
/// observation direction r̂_obs.
///
/// Returns the complex vector `ΔN = ΔN_θ θ̂ + ΔN_φ φ̂` packed as [Nx, Ny, Nz].
pub fn ptd_far_field_contribution(
    edges: &[BoundaryEdge],
    wave: &PlaneWave,
    k: f64,
    r_hat: &[f64; 3],    // observation unit vector
    e_inc_at: &dyn Fn(&[f64; 3]) -> [Complex64; 3],  // incident E field at a point
) -> [Complex64; 3] {
    let k_hat = wave.k_hat();
    let neg_kh = [-k_hat[0], -k_hat[1], -k_hat[2]];

    let mut n_ptd = [Complex64::ZERO; 3];

    for edge in edges {
        if edge.length < MIN_EDGE_LEN { continue; }
        if edge.face_normal == [0.0f64; 3] { continue; }

        let t = &edge.tangent;
        let m = &edge.edge_normal;

        // ── incidence angle ϕ' in the edge-normal plane ──────────────────────
        // ϕ' = angle between -k̂_inc and the edge-normal m̂, measured in [0, π]
        let cos_phi_inc = dot3(&neg_kh, m);
        // Only illuminate edges where ϕ' ∈ (0, π) (not grazing / parallel)
        if cos_phi_inc.abs() < 1e-4 { continue; }
        let phi_inc = cos_phi_inc.clamp(-1.0, 1.0).acos();

        // ── observation angle ϕ in the edge-normal plane ─────────────────────
        let neg_rhat = [-r_hat[0], -r_hat[1], -r_hat[2]];
        let cos_phi_obs = dot3(&neg_rhat, m);
        let phi_obs = cos_phi_obs.clamp(-1.0, 1.0).acos();

        // ── Fringe coefficient ────────────────────────────────────────────────
        let d_fringe = half_plane_fringe(phi_obs, phi_inc);
        if d_fringe.abs() < 1e-12 { continue; }

        // ── Incident field at edge midpoint ───────────────────────────────────
        let e_mid = e_inc_at(&edge.midpoint);

        // ── TE component: E parallel to edge (along t̂) ───────────────────────
        let e_te: Complex64 = e_mid[0] * t[0] + e_mid[1] * t[1] + e_mid[2] * t[2];

        // ── TM component: E along m̂ ──────────────────────────────────────────
        let e_tm: Complex64 = e_mid[0] * m[0] + e_mid[1] * m[1] + e_mid[2] * m[2];

        // ── Fringe current vector ─────────────────────────────────────────────
        // J_fringe ≈ d_fringe · (TE part along t̂ + TM part along m̂)
        // ΔN = J_fringe · L · exp(jk r̂·r_e)  [mid-point rule]
        let phase_exp = {
            let phase = k * (r_hat[0]*edge.midpoint[0]
                           + r_hat[1]*edge.midpoint[1]
                           + r_hat[2]*edge.midpoint[2]);
            Complex64::new(0.0, phase).exp()
        };

        let scale = Complex64::new(d_fringe * edge.length, 0.0) * phase_exp;

        // TE contribution: (E_te) * t̂ component → adds to N_t̂
        let n_te = e_te * scale;
        // TM contribution: (E_tm) * m̂ component → adds to N_m̂
        let n_tm = e_tm * scale;

        for i in 0..3 {
            n_ptd[i] += n_te * t[i] + n_tm * m[i];
        }
    }

    n_ptd
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use rem_mom::surface_mesh::{SurfaceMesh, TriFace};

    /// Single boundary edge of length 1 along x-axis.
    fn single_edge_surf() -> SurfaceMesh {
        let nodes = vec![
            [0.0_f64, 0.0, 0.0],
            [1.0,     0.0, 0.0],
            [0.5,     1.0, 0.0], // off-edge vertex → forms a face
        ];
        let faces = vec![TriFace {
            nodes: [0, 1, 2],
            centroid: [0.5, 1.0/3.0, 0.0],
            normal: [0.0, 0.0, 1.0],
            area: 0.5,
        }];
        // boundary edge: [0, 1] (the bottom edge)
        let boundary_edges = vec![[0, 1]];
        SurfaceMesh { nodes, faces, edges: vec![], boundary_edges }
    }

    #[test]
    fn extract_boundary_edges_length() {
        let surf = single_edge_surf();
        let bedges = extract_boundary_edges(&surf);
        assert_eq!(bedges.len(), 1);
        assert!((bedges[0].length - 1.0).abs() < 1e-12);
    }

    #[test]
    fn half_plane_fringe_symmetry() {
        // g(ϕ, ϕ') should equal g(ϕ', ϕ) by reciprocity
        let phi1 = 1.0_f64;
        let phi2 = 0.8_f64;
        let g12 = half_plane_fringe(phi1, phi2);
        let g21 = half_plane_fringe(phi2, phi1);
        assert!((g12 - g21).abs() < 1e-12, "g12={g12}, g21={g21}");
    }

    #[test]
    fn half_plane_fringe_backscatter() {
        // At exact backscatter ϕ = ϕ', g = 1/1 - 1/cos(ϕ) → non-zero for ϕ ≠ 0, π
        let phi = std::f64::consts::PI * 0.4;
        let g = half_plane_fringe(phi, phi);
        // Should be finite and non-NaN
        assert!(g.is_finite());
    }

    #[test]
    fn ptd_contribution_finite() {
        let surf = single_edge_surf();
        let edges = extract_boundary_edges(&surf);
        let wave = PlaneWave {
            theta_inc: 0.0,
            phi_inc: 0.0,
            pol: "theta".to_string(),
        };
        let k = 10.0;
        let r_hat = [1.0_f64, 0.0, 0.0];
        let e_fn = |_: &[f64; 3]| -> [Complex64; 3] {
            [Complex64::new(1.0, 0.0), Complex64::ZERO, Complex64::ZERO]
        };
        let n = ptd_far_field_contribution(&edges, &wave, k, &r_hat, &e_fn);
        for i in 0..3 {
            assert!(n[i].re.is_finite() && n[i].im.is_finite(), "n[{i}] not finite");
        }
    }
}
