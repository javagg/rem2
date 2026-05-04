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
