//! Surface impedance boundary condition (SIBC) support for MoM RWG systems.
//!
//! Adds the Leontovich impedance term
//!   Z_s * \int_\Gamma f_m · f_n dS
//! to the PEC EFIE/CFIE impedance matrix, which models conductor loss through
//! a surface impedance rather than a perfect electric conductor.
//!
//! Supports roughness-corrected effective conductivity via
//! `rem_materials::Material::effective_conductivity(freq)`.

use crate::basis::rwg::RwgBasis;
use crate::quadrature::TriQuad;
use crate::surface_mesh::SurfaceMesh;
use nalgebra::DMatrix;
use num_complex::Complex64;
use rem_core::MU0;
use std::f64::consts::PI;

/// Skin-depth surface impedance for a good conductor.
///
/// Zs = (1 + j) / (sigma * delta_s),   delta_s = sqrt(2 / (omega * mu0 * sigma))
///
/// `sigma` should be the effective conductivity, possibly roughness-corrected
/// (use [`rem_materials::Material::effective_conductivity`] before calling).
pub fn surface_impedance_from_conductivity(sigma: f64, freq: f64) -> Complex64 {
    if sigma <= 0.0 || freq <= 0.0 {
        return Complex64::ZERO;
    }
    let omega = 2.0 * PI * freq;
    let delta_s = (2.0 / (omega * MU0 * sigma)).sqrt();
    Complex64::new(1.0, 1.0) / (sigma * delta_s)
}

/// Apply the SIBC correction Z += Zs * M_surf where
/// M_surf[m,n] = \int_\Gamma f_m · f_n dS.
///
/// `sigma` should be the conductor's **effective** conductivity including
/// surface roughness correction if applicable. Use
/// `rem_materials::Material::effective_conductivity(freq)` to compute it
/// from bulk conductivity and roughness parameters.
pub fn apply_sibc_rwg(
    z_mat: &mut DMatrix<Complex64>,
    surf: &SurfaceMesh,
    bases: &[RwgBasis],
    freq: f64,
    sigma: f64,
    quad: &TriQuad,
) {
    let z_s = surface_impedance_from_conductivity(sigma, freq);
    if z_s.norm() < 1e-30 {
        return;
    }

    let n = bases.len();
    for i in 0..n {
        for j in i..n {
            let overlap = rwg_surface_overlap(&bases[i], &bases[j], surf, quad);
            if overlap.abs() < 1e-30 {
                continue;
            }
            let delta = z_s * overlap;
            z_mat[(i, j)] += delta;
            if i != j {
                z_mat[(j, i)] += delta;
            }
        }
    }
}

fn rwg_surface_overlap(
    bm: &RwgBasis,
    bn: &RwgBasis,
    surf: &SurfaceMesh,
    quad: &TriQuad,
) -> f64 {
    let mut sum = 0.0;

    for &(m_face, m_plus) in &[(bm.plus_face, true), (bm.minus_face, false)] {
        for &(n_face, n_plus) in &[(bn.plus_face, true), (bn.minus_face, false)] {
            if m_face != n_face {
                continue;
            }
            let face = &surf.faces[m_face];
            for (bary, &w) in quad.bary.iter().zip(quad.weights.iter()) {
                let r = TriQuad::global_point(bary, face, &surf.nodes);
                let fm = bm.eval(&r, surf, m_plus);
                let fn_ = bn.eval(&r, surf, n_plus);
                let dot = fm[0] * fn_[0] + fm[1] * fn_[1] + fm[2] * fn_[2];
                sum += dot * (w * 2.0 * face.area);
            }
        }
    }

    sum
}

/// Apply the SIBC correction to a **pulse-basis** EFIE matrix.
///
/// For pulse basis functions `f_m = 1` on face `m` and 0 elsewhere, the
/// surface Gram-matrix overlap integral simplifies to a diagonal:
///
///   ∫_Γ f_m · f_n dS = A_m · δ_{mn}
///
/// where `A_m` is the area of face `m`.  The SIBC correction therefore only
/// adds to the diagonal:
///
///   Z[m,m] += Z_s · A_m
pub fn apply_sibc_pulse(
    z_mat: &mut DMatrix<Complex64>,
    surf: &SurfaceMesh,
    freq: f64,
    sigma: f64,
) {
    let z_s = surface_impedance_from_conductivity(sigma, freq);
    if z_s.norm() < 1e-30 {
        return;
    }
    for (i, face) in surf.faces.iter().enumerate() {
        z_mat[(i, i)] += z_s * face.area;
    }
}

/// Apply the SIBC correction to a **pulse-basis** EFIE matrix with optional roughness + superconductor.
pub fn apply_sibc_pulse_with_config(
    z_mat: &mut DMatrix<Complex64>,
    surf: &SurfaceMesh,
    freq: f64,
    sigma: f64,
    roughness_model: Option<&str>,
    rms_roughness_m: Option<f64>,
    superconductor_ls: f64,
    superconductor_rdc: f64,
    superconductor_rrf: f64,
    superconductor_xdc: f64,
) {
    if sigma <= 0.0 && superconductor_ls <= 0.0 {
        return;
    }
    let z_s = if superconductor_ls > 0.0 || superconductor_rdc > 0.0 || superconductor_rrf > 0.0 || superconductor_xdc > 0.0 {
        surface_impedance_superconductor(superconductor_ls, superconductor_rdc, superconductor_rrf, superconductor_xdc, freq)
    } else {
        match (roughness_model, rms_roughness_m) {
            (Some(model), Some(rms)) if rms > 0.0 => {
                surface_impedance_with_roughness(sigma, freq, rms, model)
            }
            _ => surface_impedance_from_conductivity(sigma, freq),
        }
    };
    if z_s.norm() < 1e-30 {
        return;
    }
    for (i, face) in surf.faces.iter().enumerate() {
        z_mat[(i, i)] += z_s * face.area;
    }
}

// ── Conductor surface roughness models ─────────────────────────────────────

/// Hammerstad–Jensen roughness correction factor K_r.
///
/// Models the effective increase in surface resistance due to conductor
/// surface roughness (rms roughness Δ) relative to the skin depth δ_s.
///
///   K_r = 1 + (1/2)·[1 + erf(−1.4·(Δ/δ_s − 1.9))]
///
/// * `rms_roughness_m` — RMS surface roughness Δ [m]
/// * `freq`            — frequency [Hz]
/// * `sigma`           — conductor conductivity [S/m]
///
/// Returns a multiplicative factor ≥ 1.
pub fn roughness_factor_hammerstad(rms_roughness_m: f64, freq: f64, sigma: f64) -> f64 {
    if freq <= 0.0 || sigma <= 0.0 || rms_roughness_m <= 0.0 {
        return 1.0;
    }
    let omega = 2.0 * PI * freq;
    let delta_s = (2.0 / (omega * MU0 * sigma)).sqrt();
    let ratio = rms_roughness_m / delta_s;
    // erf approximation via Horner / series — use std's f64::erf via libm if available,
    // otherwise a fast rational approximation (max error < 1.5e-7).
    let x = -1.4 * (ratio - 1.9);
    let erf_x = erf_approx(x);
    1.0 + 0.5 * (1.0 + erf_x)
}

/// Groisse roughness correction factor K_r.
///
///   K_r = 1 + (2/π)·arctan[1.4·(Δ/δ_s)²]
///
/// * `rms_roughness_m` — RMS surface roughness Δ [m]
/// * `freq`            — frequency [Hz]
/// * `sigma`           — conductor conductivity [S/m]
pub fn roughness_factor_groisse(rms_roughness_m: f64, freq: f64, sigma: f64) -> f64 {
    if freq <= 0.0 || sigma <= 0.0 || rms_roughness_m <= 0.0 {
        return 1.0;
    }
    let omega = 2.0 * PI * freq;
    let delta_s = (2.0 / (omega * MU0 * sigma)).sqrt();
    let ratio = rms_roughness_m / delta_s;
    1.0 + (2.0 / PI) * (1.4 * ratio * ratio).atan()
}

/// Compute skin-depth surface impedance with optional roughness correction.
///
/// * `model` — `"hammerstad"` or `"groisse"` (case-insensitive); any other
///             value disables roughness correction (K_r = 1).
pub fn surface_impedance_with_roughness(
    sigma: f64,
    freq: f64,
    rms_roughness_m: f64,
    model: &str,
) -> Complex64 {
    let z_smooth = surface_impedance_from_conductivity(sigma, freq);
    if rms_roughness_m <= 0.0 {
        return z_smooth;
    }
    let k_r = match model.to_lowercase().as_str() {
        "hammerstad" | "hammerstad-jensen" => {
            roughness_factor_hammerstad(rms_roughness_m, freq, sigma)
        }
        "groisse" => roughness_factor_groisse(rms_roughness_m, freq, sigma),
        _ => 1.0,
    };
    z_smooth * k_r
}

/// Fast rational approximation of the error function (max error < 1.5e-7).
///
/// Based on Abramowitz & Stegun formula 7.1.26.
fn erf_approx(x: f64) -> f64 {
    const P: f64 = 0.3275911;
    let t = 1.0 / (1.0 + P * x.abs());
    let poly = t * (0.254829592
        + t * (-0.284496736
            + t * (1.421413741
                + t * (-1.453152027 + t * 1.061405429))));
    let result = 1.0 - poly * (-x * x).exp();
    if x >= 0.0 { result } else { -result }
}

/// Superconducting surface impedance: Zs = Rdc + Rrf·√f + j·(2π·f·Ls + Xdc).
/// Parameters: Ls [H], Rdc [Ω], Rrf [Ω/√Hz], Xdc [Ω].
pub fn surface_impedance_superconductor(
    ls: f64, rdc: f64, rrf: f64, xdc: f64, freq: f64,
) -> Complex64 {
    if ls <= 0.0 && rdc <= 0.0 && rrf <= 0.0 && xdc <= 0.0 {
        return Complex64::ZERO;
    }
    let omega = 2.0 * PI * freq;
    let re = rdc + rrf * freq.sqrt();
    let im = omega * ls + xdc;
    Complex64::new(re, im)
}

/// Apply SIBC correction to RWG-basis matrix with optional roughness + superconductor.
pub fn apply_sibc_rwg_with_config(
    z_mat: &mut DMatrix<Complex64>,
    surf: &SurfaceMesh,
    bases: &[RwgBasis],
    freq: f64,
    sigma: f64,
    roughness_model: Option<&str>,
    rms_roughness_m: Option<f64>,
    superconductor_ls: f64,
    superconductor_rdc: f64,
    superconductor_rrf: f64,
    superconductor_xdc: f64,
    quad: &TriQuad,
) {
    if sigma <= 0.0 && superconductor_ls <= 0.0 {
        return;
    }
    let z_s = if superconductor_ls > 0.0 || superconductor_rdc > 0.0 || superconductor_rrf > 0.0 || superconductor_xdc > 0.0 {
        surface_impedance_superconductor(superconductor_ls, superconductor_rdc, superconductor_rrf, superconductor_xdc, freq)
    } else {
        match (roughness_model, rms_roughness_m) {
            (Some(model), Some(rms)) if rms > 0.0 => {
                surface_impedance_with_roughness(sigma, freq, rms, model)
            }
            _ => surface_impedance_from_conductivity(sigma, freq),
        }
    };
    if z_s.norm() < 1e-30 {
        return;
    }
    let n = bases.len();
    for i in 0..n {
        for j in i..n {
            let overlap = rwg_surface_overlap(&bases[i], &bases[j], surf, quad);
            if overlap.abs() < 1e-30 {
                continue;
            }
            let delta = z_s * overlap;
            z_mat[(i, j)] += delta;
            if i != j {
                z_mat[(j, i)] += delta;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::surface_mesh::{SurfaceMesh, TriFace, SharedEdge, tri_geometry};
    use nalgebra::DMatrix;
    use num_complex::Complex64;

    fn two_face_surf() -> SurfaceMesh {
        let nodes = vec![
            [0.0_f64, 0.0, 0.0],
            [1.0,     0.0, 0.0],
            [0.0,     1.0, 0.0],
            [1.0,     1.0, 0.0],
        ];
        let (c0, n0, a0) = tri_geometry(&nodes[0], &nodes[1], &nodes[2]);
        let (c1, n1, a1) = tri_geometry(&nodes[1], &nodes[3], &nodes[2]);
        let faces = vec![
            TriFace { nodes: [0,1,2], centroid: c0, normal: n0, area: a0 },
            TriFace { nodes: [1,3,2], centroid: c1, normal: n1, area: a1 },
        ];
        let elen = {
            let d = [nodes[1][0]-nodes[2][0], nodes[1][1]-nodes[2][1], nodes[1][2]-nodes[2][2]];
            (d[0]*d[0]+d[1]*d[1]+d[2]*d[2]).sqrt()
        };
        let edges = vec![SharedEdge { nodes: [1, 2], plus_face: 0, minus_face: 1, length: elen }];
        SurfaceMesh { nodes, faces, edges, boundary_edges: vec![], face_attrs: vec![0, 0], global_node_ids: vec![0, 1, 2, 3] }
    }

    #[test]
    fn sibc_pulse_diagonal_only() {
        let surf = two_face_surf();
        let n = surf.faces.len();
        let mut z = DMatrix::<Complex64>::zeros(n, n);
        apply_sibc_pulse(&mut z, &surf, 1e9, 5.96e7);
        // Off-diagonal must remain zero
        for i in 0..n {
            for j in 0..n {
                if i != j {
                    assert!(z[(i, j)].norm() < 1e-30, "off-diagonal non-zero");
                }
            }
        }
        // Diagonal must be non-zero (positive real part)
        for i in 0..n {
            assert!(z[(i, i)].re > 0.0, "diagonal real part must be positive");
        }
    }

    #[test]
    fn roughness_factor_hammerstad_limits() {
        // K_r = 1 when roughness = 0
        let kr = roughness_factor_hammerstad(0.0, 1e9, 5.96e7);
        assert!((kr - 1.0).abs() < 1e-10);
        // K_r >= 1 for nonzero roughness
        let kr2 = roughness_factor_hammerstad(1e-6, 1e9, 5.96e7);
        assert!(kr2 >= 1.0);
        // K_r ≤ 2 (physical upper bound)
        assert!(kr2 <= 2.0 + 1e-10);
    }

    #[test]
    fn roughness_factor_groisse_limits() {
        let kr0 = roughness_factor_groisse(0.0, 1e9, 5.96e7);
        assert!((kr0 - 1.0).abs() < 1e-10);
        let kr2 = roughness_factor_groisse(5e-6, 10e9, 5.96e7);
        assert!(kr2 >= 1.0);
        assert!(kr2 < 2.0 + 1e-10);
    }

    #[test]
    fn erf_approx_basic() {
        assert!(erf_approx(0.0).abs() < 1e-6);
        assert!((erf_approx(1.0) - 0.8427).abs() < 2e-4);
        assert!((erf_approx(-1.0) + 0.8427).abs() < 2e-4);
    }

    #[test]
    fn superconductor_zs_niobium_10ghz() {
        // Nb: Ls=0.11 pH, Rdc=Rrf=Xdc=0 at 10 GHz → Zs = j·6.91 mΩ
        let zs = surface_impedance_superconductor(0.11e-12, 0.0, 0.0, 0.0, 10.0e9);
        assert!((zs.re - 0.0).abs() < 1e-15);
        assert!((zs.im - 6.912e-3).abs() < 1e-6);
    }

    #[test]
    fn superconductor_zs_with_rdc() {
        let zs = surface_impedance_superconductor(0.5e-12, 0.1, 0.0, 0.0, 1.0e9);
        assert!((zs.re - 0.1).abs() < 1e-15);
        assert!((zs.im - 3.1416e-3).abs() < 1e-5);
    }

    #[test]
    fn sibc_pulse_with_config_no_sigma_returns_early() {
        let surf = two_face_surf();
        let mut z = DMatrix::<Complex64>::zeros(2, 2);
        apply_sibc_pulse_with_config(&mut z, &surf, 1e9, 0.0, None, None, 0.0, 0.0, 0.0, 0.0);
        assert_eq!(z[(0, 0)], Complex64::ZERO);
    }
}
