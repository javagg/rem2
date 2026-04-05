//! Far-field Physical Optics integral and RCS computation for SBR+.
//!
//! Accumulates PO surface currents from all ray hits and integrates them
//! to obtain the scattered far-field and RCS pattern.

use num_complex::Complex64;
use rem_core::{ETA0, RemResult};
use rem_mom::surface_mesh::SurfaceMesh;
use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Face current storage
// ---------------------------------------------------------------------------

/// Accumulated PO currents on a single triangular face.
#[derive(Debug, Clone, Default)]
pub struct FaceCurrent {
    /// Electric surface current density J [A/m] (sum over all bounces)
    pub j: [Complex64; 3],
    /// Equivalent magnetic current density M [V/m] (zero for PEC)
    pub m: [Complex64; 3],
}

/// Per-face current map indexed like `SurfaceMesh::faces`.
pub type CurrentMap = Vec<FaceCurrent>;

/// Allocate a zeroed current map for a surface mesh.
pub fn zero_currents(surf: &SurfaceMesh) -> CurrentMap {
    vec![FaceCurrent::default(); surf.faces.len()]
}

// ---------------------------------------------------------------------------
// Far-field PO integral
// ---------------------------------------------------------------------------

/// Compute the bistatic RCS [m²] at observation direction r̂(θ_s, φ_s).
///
/// ```text
/// N(r̂) = Σ_m J_m * A_m * exp(+jk r̂·r_m)
/// L(r̂) = Σ_m M_m * A_m * exp(+jk r̂·r_m)
///
/// E_scat = -jkη₀/(4π) [ r̂×(r̂×N) + r̂×L/η₀ ]
///
/// σ(r̂) = 4π |E_scat|² / |E_inc|²   [m²]
/// ```
///
/// For PEC targets M=0, and |E_inc| = 1 (unit amplitude plane wave).
pub fn rcs_at(
    currents: &CurrentMap,
    surf: &SurfaceMesh,
    k: f64,
    theta_s: f64,
    phi_s: f64,
) -> f64 {
    let (st, ct) = (theta_s.sin(), theta_s.cos());
    let (sp, cp) = (phi_s.sin(), phi_s.cos());
    let r_hat = [st * cp, st * sp, ct];

    // Radiation vectors N and L
    let mut nx = Complex64::ZERO;
    let mut ny = Complex64::ZERO;
    let mut nz = Complex64::ZERO;
    let mut lx = Complex64::ZERO;
    let mut ly = Complex64::ZERO;
    let mut lz = Complex64::ZERO;

    for (fc, face) in currents.iter().zip(surf.faces.iter()) {
        let r = &face.centroid;
        let phase = k * (r_hat[0]*r[0] + r_hat[1]*r[1] + r_hat[2]*r[2]);
        let phasor = Complex64::new(0.0, phase).exp();
        let a = face.area;

        nx += fc.j[0] * phasor * a;
        ny += fc.j[1] * phasor * a;
        nz += fc.j[2] * phasor * a;
        lx += fc.m[0] * phasor * a;
        ly += fc.m[1] * phasor * a;
        lz += fc.m[2] * phasor * a;
    }

    // r̂ × (r̂ × N)  =  N − (r̂·N) r̂
    let rn = r_hat[0]*nx + r_hat[1]*ny + r_hat[2]*nz;
    let rr_n = [rn * r_hat[0], rn * r_hat[1], rn * r_hat[2]];
    let ex = nx - rr_n[0];
    let ey = ny - rr_n[1];
    let ez = nz - rr_n[2];

    // r̂ × L (cross product)
    let rl_x = Complex64::new(r_hat[1], 0.0)*lz - Complex64::new(r_hat[2], 0.0)*ly;
    let rl_y = Complex64::new(r_hat[2], 0.0)*lx - Complex64::new(r_hat[0], 0.0)*lz;
    let rl_z = Complex64::new(r_hat[0], 0.0)*ly - Complex64::new(r_hat[1], 0.0)*lx;

    let eta0 = ETA0;
    let prefac = k * eta0 / (4.0 * PI); // |−jkη₀/(4π)| = kη₀/(4π)

    // E_scat components (magnitude factor, phase cancels in |·|²)
    let esx = prefac * (ex + rl_x / eta0);
    let esy = prefac * (ey + rl_y / eta0);
    let esz = prefac * (ez + rl_z / eta0);

    let e_sq = esx.norm_sqr() + esy.norm_sqr() + esz.norm_sqr();

    // σ = 4π |E_scat|² / |E_inc|²   (|E_inc|² = 1 for unit amplitude)
    4.0 * PI * e_sq
}

/// Compute RCS pattern over all (θ, φ) observation angles.
/// Returns 2-D array `result[i_theta][i_phi]` in [m²].
pub fn rcs_pattern(
    currents: &CurrentMap,
    surf: &SurfaceMesh,
    k: f64,
    theta_deg: &[f64],
    phi_deg: &[f64],
) -> Vec<Vec<f64>> {
    theta_deg.iter().map(|&th| {
        phi_deg.iter().map(|&ph| {
            rcs_at(currents, surf, k, th.to_radians(), ph.to_radians())
        }).collect()
    }).collect()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use rem_mom::surface_mesh::{SurfaceMesh, TriFace};
    use rem_mom::surface_mesh::tri_geometry;

    fn flat_surf() -> SurfaceMesh {
        let nodes = vec![[0.0f64,0.0,0.0],[1.0,0.0,0.0],[0.0,1.0,0.0]];
        let (c,n,a) = tri_geometry(&nodes[0],&nodes[1],&nodes[2]);
        SurfaceMesh {
            nodes, faces: vec![TriFace{nodes:[0,1,2],centroid:c,normal:n,area:a}],
            edges: vec![], boundary_edges: vec![],
        }
    }

    #[test]
    fn zero_current_gives_zero_rcs() {
        let surf = flat_surf();
        let cur = zero_currents(&surf);
        let sigma = rcs_at(&cur, &surf, 10.0, 0.0, 0.0);
        assert!(sigma < 1e-30);
    }

    #[test]
    fn nonzero_current_gives_positive_rcs() {
        let surf = flat_surf();
        let mut cur = zero_currents(&surf);
        cur[0].j[0] = Complex64::new(1.0, 0.0);
        let sigma = rcs_at(&cur, &surf, 10.0, 0.0, 0.0);
        assert!(sigma > 0.0, "RCS should be positive, got {}", sigma);
    }

    #[test]
    fn pattern_shape() {
        let surf = flat_surf();
        let cur = zero_currents(&surf);
        let pat = rcs_pattern(&cur, &surf, 10.0, &[0.0, 90.0], &[0.0, 90.0]);
        assert_eq!(pat.len(), 2);
        assert_eq!(pat[0].len(), 2);
    }
}
