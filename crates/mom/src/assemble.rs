//! Impedance matrix assembly for MoM (pulse and RWG basis functions).

use crate::surface_mesh::SurfaceMesh;
use crate::quadrature::TriQuad;
use crate::green::green3d;
use crate::singular::zmn_self_duffy_pulse;
use crate::basis::rwg::RwgBasis;
use num_complex::Complex64;
use rem_core::{RemError, RemResult, EPS0, MU0, C0};
use std::f64::consts::PI;
use rayon::prelude::*;
use faer::Mat;

// faer::c64 == num_complex::Complex64 in faer 0.21
type C64 = faer::c64;

/// Assemble N×N impedance matrix for scalar EFIE with pulse (constant) basis functions.
///
/// Z[m,n] = -jωμ₀ ∫_Tm ∫_Tn G(r,r') dS' dS
///
/// Diagonal blocks use Duffy self-integral; off-diagonal use standard Gauss quadrature.
pub fn assemble_efie_pulse(
    surf: &SurfaceMesh,
    freq: f64,
    quad: &TriQuad,
    _singular_tol: f64,
) -> RemResult<Mat<C64>> {
    let n = surf.faces.len();
    let omega = 2.0 * PI * freq;
    let k     = omega / C0;
    let omega_mu0 = omega * MU0;

    let mut z = Mat::<C64>::zeros(n, n);

    // Parallel assembly: compute each column independently
    let cols: Vec<Vec<C64>> = (0..n).into_par_iter().map(|ni| {
        let face_n = &surf.faces[ni];
        let mut col = vec![C64::new(0.0, 0.0); n];
        for mi in 0..n {
            let val = if mi == ni {
                let zself = zmn_self_duffy_pulse(face_n, &surf.nodes, k, omega_mu0, 4);
                to_c64(zself)
            } else {
                let face_m = &surf.faces[mi];
                let zoff = zmn_regular_pulse(face_m, face_n, &surf.nodes, k, omega_mu0, quad);
                to_c64(zoff)
            };
            col[mi] = val;
        }
        col
    }).collect();

    for (ni, col) in cols.into_iter().enumerate() {
        for mi in 0..n {
            z[(mi, ni)] = col[mi];
        }
    }

    Ok(z)
}

/// Off-diagonal element Z[m,n] using standard Gauss quadrature on Tm × Tn.
fn zmn_regular_pulse(
    face_m: &crate::surface_mesh::TriFace,
    face_n: &crate::surface_mesh::TriFace,
    nodes: &[[f64; 3]],
    k: f64,
    omega_mu0: f64,
    quad: &TriQuad,
) -> Complex64 {
    let mut val = Complex64::ZERO;
    for (bm, &wm) in quad.bary.iter().zip(quad.weights.iter()) {
        let rm = TriQuad::global_point(bm, face_m, nodes);
        for (bn, &wn) in quad.bary.iter().zip(quad.weights.iter()) {
            let rn = TriQuad::global_point(bn, face_n, nodes);
            let g  = green3d(&rm, &rn, k);
            val += g * (wm * wn * 4.0 * face_m.area * face_n.area);
        }
    }
    Complex64::new(0.0, -omega_mu0) * val
}

// ---------------------------------------------------------------------------
// RWG CFIE
// ---------------------------------------------------------------------------

/// Assemble N×N impedance matrix using RWG basis functions and CFIE.
///
/// Z_CFIE = alpha * Z_EFIE + (1 - alpha) * eta0 * Z_MFIE
pub fn assemble_cfie_rwg(
    surf: &SurfaceMesh,
    bases: &[RwgBasis],
    freq: f64,
    alpha: f64,
    quad: &TriQuad,
    _singular_tol: f64,
) -> RemResult<Mat<C64>> {
    let n = bases.len();
    if n == 0 {
        return Err(RemError::Mesh("No RWG bases found — check surface mesh".to_string()));
    }

    let omega   = 2.0 * PI * freq;
    let k       = omega / C0;
    let eta0    = (MU0 / EPS0).sqrt();

    let z_efie = assemble_efie_rwg(surf, bases, k, omega, quad)?;

    if alpha >= 1.0 - 1e-9 {
        return Ok(z_efie);
    }

    let z_mfie = assemble_mfie_rwg(surf, bases, k, quad)?;

    let mut z = Mat::<C64>::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            let v = z_efie[(i, j)] * C64::new(alpha, 0.0)
                  + z_mfie[(i, j)] * C64::new((1.0 - alpha) * eta0, 0.0);
            z[(i, j)] = v;
        }
    }
    Ok(z)
}

/// Assemble RWG EFIE impedance matrix.
fn assemble_efie_rwg(
    surf: &SurfaceMesh,
    bases: &[RwgBasis],
    k: f64,
    omega: f64,
    quad: &TriQuad,
) -> RemResult<Mat<C64>> {
    let n = bases.len();
    let omega_mu0 = omega * MU0;
    let inv_omega_eps0 = 1.0 / (omega * EPS0);

    let cols: Vec<Vec<C64>> = (0..n).into_par_iter().map(|ni| {
        let bn = &bases[ni];
        let mut col = vec![C64::new(0.0, 0.0); n];
        for mi in 0..n {
            let bm = &bases[mi];
            let val = zmn_efie_rwg(bm, bn, surf, k, omega_mu0, inv_omega_eps0, quad);
            col[mi] = to_c64(val);
        }
        col
    }).collect();

    let mut z = Mat::<C64>::zeros(n, n);
    for (ni, col) in cols.into_iter().enumerate() {
        for mi in 0..n {
            z[(mi, ni)] = col[mi];
        }
    }
    Ok(z)
}

/// Single EFIE matrix element Z_EFIE[m,n] for RWG bases.
fn zmn_efie_rwg(
    bm: &RwgBasis,
    bn: &RwgBasis,
    surf: &SurfaceMesh,
    k: f64,
    omega_mu0: f64,
    inv_omega_eps0: f64,
    quad: &TriQuad,
) -> Complex64 {
    // Integrate over T_m (observation) and T_n (source), both ± halves
    let mut val = Complex64::ZERO;

    for &(m_face, m_plus) in &[(bm.plus_face, true), (bm.minus_face, false)] {
        for &(n_face, n_plus) in &[(bn.plus_face, true), (bn.minus_face, false)] {
            let face_m = &surf.faces[m_face];
            let face_n = &surf.faces[n_face];
            let div_n  = bn.divergence(surf, n_plus);

            for (bm_pt, &wm) in quad.bary.iter().zip(quad.weights.iter()) {
                let rm = TriQuad::global_point(bm_pt, face_m, &surf.nodes);
                let fm = bm.eval(&rm, surf, m_plus);

                for (bn_pt, &wn) in quad.bary.iter().zip(quad.weights.iter()) {
                    let rn = TriQuad::global_point(bn_pt, face_n, &surf.nodes);
                    let fn_ = bn.eval(&rn, surf, n_plus);
                    let g  = green3d(&rm, &rn, k);

                    let dot_ff = fm[0]*fn_[0] + fm[1]*fn_[1] + fm[2]*fn_[2];
                    let div_m  = bm.divergence(surf, m_plus);

                    let integrand = g * (dot_ff - inv_omega_eps0 / omega_mu0 * div_m * div_n);
                    val += integrand * (wm * wn * 4.0 * face_m.area * face_n.area);
                }
            }
        }
    }

    Complex64::new(0.0, -omega_mu0) * val
}

/// Assemble MFIE matrix (identity + curl-Green term).
fn assemble_mfie_rwg(
    surf: &SurfaceMesh,
    bases: &[RwgBasis],
    k: f64,
    quad: &TriQuad,
) -> RemResult<Mat<C64>> {
    let n = bases.len();

    let cols: Vec<Vec<C64>> = (0..n).into_par_iter().map(|ni| {
        let bn = &bases[ni];
        let mut col = vec![C64::new(0.0, 0.0); n];
        for mi in 0..n {
            let bm = &bases[mi];
            let val = zmn_mfie_rwg(bm, bn, surf, k, quad);
            col[mi] = to_c64(val);
        }
        col
    }).collect();

    let mut z = Mat::<C64>::zeros(n, n);
    for (ni, col) in cols.into_iter().enumerate() {
        for mi in 0..n {
            z[(mi, ni)] = col[mi];
        }
    }
    Ok(z)
}

fn zmn_mfie_rwg(
    bm: &RwgBasis,
    bn: &RwgBasis,
    surf: &SurfaceMesh,
    k: f64,
    quad: &TriQuad,
) -> Complex64 {
    // δ_{mn}/2 term (identity)
    let identity_term = if bm.edge_idx == bn.edge_idx {
        // ∫_T f_m · f_n dS (overlap integral)
        let mut overlap = 0.0f64;
        for &(face_idx, in_plus) in &[
            (bm.plus_face, true), (bm.minus_face, false)
        ] {
            let face = &surf.faces[face_idx];
            for (b_pt, &w) in quad.bary.iter().zip(quad.weights.iter()) {
                let r = TriQuad::global_point(b_pt, face, &surf.nodes);
                let f = bm.eval(&r, surf, in_plus);
                overlap += (f[0]*f[0] + f[1]*f[1] + f[2]*f[2]) * (w * 2.0 * face.area);
            }
        }
        Complex64::new(0.5 * overlap, 0.0)
    } else {
        Complex64::ZERO
    };

    // curl-Green integral term
    let mut curl_term = Complex64::ZERO;
    for &(m_face, m_plus) in &[(bm.plus_face, true), (bm.minus_face, false)] {
        for &(n_face, n_plus) in &[(bn.plus_face, true), (bn.minus_face, false)] {
            let face_m = &surf.faces[m_face];
            let face_n = &surf.faces[n_face];
            let nm = face_m.normal;

            for (bm_pt, &wm) in quad.bary.iter().zip(quad.weights.iter()) {
                let rm = TriQuad::global_point(bm_pt, face_m, &surf.nodes);
                let fm = bm.eval(&rm, surf, m_plus);

                for (bn_pt, &wn) in quad.bary.iter().zip(quad.weights.iter()) {
                    let rn = TriQuad::global_point(bn_pt, face_n, &surf.nodes);
                    let fn_ = bn.eval(&rn, surf, n_plus);
                    let grad_g = green_gradient(&rm, &rn, k);

                    // ∇G × f_n (fn_ is real, grad_g is complex)
                    let fn_c = [
                        Complex64::new(fn_[0], 0.0),
                        Complex64::new(fn_[1], 0.0),
                        Complex64::new(fn_[2], 0.0),
                    ];
                    let curl_gfn = cross_c(&grad_g, &fn_c);

                    // n̂_m × (∇G × f_n) → dot with f_m (nm is real)
                    let nm_c = [
                        Complex64::new(nm[0], 0.0),
                        Complex64::new(nm[1], 0.0),
                        Complex64::new(nm[2], 0.0),
                    ];
                    let n_x_curl = cross_c(&nm_c, &curl_gfn);
                    let dot_val =   fm[0]*n_x_curl[0].re + fm[1]*n_x_curl[1].re + fm[2]*n_x_curl[2].re;

                    curl_term += Complex64::new(dot_val, 0.0)
                               * (wm * wn * 4.0 * face_m.area * face_n.area);
                }
            }
        }
    }

    identity_term + curl_term
}

/// ∇G(r,r') = G(r,r') * (jk + 1/R) * (r - r') / R
fn green_gradient(r: &[f64; 3], rp: &[f64; 3], k: f64) -> [Complex64; 3] {
    let dr = [r[0]-rp[0], r[1]-rp[1], r[2]-rp[2]];
    let dist = (dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]).sqrt();
    if dist < 1e-14 { return [Complex64::ZERO; 3]; }
    let g = green3d(r, rp, k);
    let factor = g * Complex64::new(1.0/dist, k) / dist;
    [factor * dr[0], factor * dr[1], factor * dr[2]]
}

fn cross_c(a: &[Complex64; 3], b: &[Complex64; 3]) -> [Complex64; 3] {
    [
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0],
    ]
}

// ---------------------------------------------------------------------------
// LU solver
// ---------------------------------------------------------------------------

/// Solve Z·I = V using LU decomposition (faer).
pub fn lu_solve(z: &Mat<C64>, rhs: &[Complex64]) -> RemResult<Vec<Complex64>> {
    use faer::linalg::solvers::Solve;

    let n = rhs.len();
    if z.nrows() != n || z.ncols() != n {
        return Err(RemError::Config(format!(
            "lu_solve: matrix {}×{} but rhs length {}", z.nrows(), z.ncols(), n
        )));
    }

    // Build RHS column vector
    let mut b = Mat::<C64>::zeros(n, 1);
    for i in 0..n {
        b[(i, 0)] = C64::new(rhs[i].re, rhs[i].im);
    }

    let lu = z.as_ref().partial_piv_lu();
    let x  = lu.solve(b.as_ref());

    Ok((0..n).map(|i| {
        let v = x[(i, 0)];
        Complex64::new(v.re, v.im)
    }).collect())
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

#[inline]
fn to_c64(v: Complex64) -> C64 { C64::new(v.re, v.im) }
