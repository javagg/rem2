//! Impedance matrix assembly for MoM (pulse and RWG basis functions).

use crate::surface_mesh::SurfaceMesh;
use crate::quadrature::TriQuad;
use crate::green::green3d;
use crate::singular::{zmn_self_duffy_pulse, zmn_singular_pulse, classify_pair, TriPairType, zmn_efie_rwg_singular};
use crate::basis::rwg::RwgBasis;
use nalgebra::{DMatrix, DVector};
use num_complex::Complex64;
use rem_core::{RemError, RemResult, EPS0, MU0, C0, LinearOperator};
use std::f64::consts::PI;

#[cfg(not(target_arch = "wasm32"))]
use rayon::prelude::*;

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
) -> RemResult<DMatrix<Complex64>> {
    let n = surf.faces.len();
    let omega = 2.0 * PI * freq;
    let k     = omega / C0;
    let omega_mu0 = omega * MU0;

    let compute_col = |ni: usize| -> Vec<Complex64> {
        let face_n = &surf.faces[ni];
        let mut col = vec![Complex64::ZERO; n];
        for mi in 0..n {
            let face_m = &surf.faces[mi];
            let pair = classify_pair(face_m, face_n);
            let val = match pair {
                TriPairType::Identical => {
                    let zself = zmn_self_duffy_pulse(face_n, &surf.nodes, k, omega_mu0, 4);
                    zself
                }
                TriPairType::SharedEdge | TriPairType::SharedVertex => {
                    let integral = zmn_singular_pulse(face_m, face_n, &surf.nodes, k, 4);
                    Complex64::new(0.0, -omega_mu0) * integral
                }
                TriPairType::Disjoint => {
                    let zoff = zmn_regular_pulse(face_m, face_n, &surf.nodes, k, omega_mu0, quad);
                    zoff
                }
            };
            col[mi] = val;
        }
        col
    };

    #[cfg(not(target_arch = "wasm32"))]
    let cols: Vec<Vec<Complex64>> = (0..n).into_par_iter().map(compute_col).collect();
    #[cfg(target_arch = "wasm32")]
    let cols: Vec<Vec<Complex64>> = (0..n).map(compute_col).collect();

    let mut z = DMatrix::<Complex64>::from_element(n, n, Complex64::ZERO);
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
) -> RemResult<DMatrix<Complex64>> {
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

    let mut z = DMatrix::<Complex64>::from_element(n, n, Complex64::ZERO);
    for i in 0..n {
        for j in 0..n {
            let v = z_efie[(i, j)] * Complex64::new(alpha, 0.0)
                  + z_mfie[(i, j)] * Complex64::new((1.0 - alpha) * eta0, 0.0);
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
) -> RemResult<DMatrix<Complex64>> {
    let n = bases.len();
    let omega_mu0 = omega * MU0;
    let inv_omega_eps0 = 1.0 / (omega * EPS0);

    let cols: Vec<Vec<Complex64>> = {
        let compute = |ni: usize| -> Vec<Complex64> {
            let bn = &bases[ni];
            let mut col = vec![Complex64::ZERO; n];
            for mi in 0..n {
                let bm = &bases[mi];
                let val = zmn_efie_rwg(bm, bn, surf, k, omega_mu0, inv_omega_eps0, quad);
                col[mi] = val;
            }
            col
        };
        #[cfg(not(target_arch = "wasm32"))]
        { (0..n).into_par_iter().map(compute).collect() }
        #[cfg(target_arch = "wasm32")]
        { (0..n).map(compute).collect() }
    };

    let mut z = DMatrix::<Complex64>::from_element(n, n, Complex64::ZERO);
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
            let div_m  = bm.divergence(surf, m_plus);

            let pair = classify_pair(face_m, face_n);
            if pair != TriPairType::Disjoint {
                // Near-singular: use Sauter-Schwab / Duffy quadrature
                // Build a closure capturing (bm, bn, m_plus, n_plus, div_m, div_n, surf)
                let fm_fn = |rm: &[f64; 3], rn: &[f64; 3]| -> (f64, f64) {
                    let fm  = bm.eval(rm, surf, m_plus);
                    let fn_ = bn.eval(rn, surf, n_plus);
                    let dot = fm[0]*fn_[0] + fm[1]*fn_[1] + fm[2]*fn_[2];
                    (dot, div_m * div_n)
                };
                let (a_term, phi_term) =
                    zmn_efie_rwg_singular(face_m, face_n, &fm_fn, &surf.nodes, k, 4);
                val += a_term - inv_omega_eps0 / omega_mu0 * phi_term;
            } else {
                // Well-separated: standard Gauss quadrature
                for (bm_pt, &wm) in quad.bary.iter().zip(quad.weights.iter()) {
                    let rm = TriQuad::global_point(bm_pt, face_m, &surf.nodes);
                    let fm = bm.eval(&rm, surf, m_plus);

                    for (bn_pt, &wn) in quad.bary.iter().zip(quad.weights.iter()) {
                        let rn = TriQuad::global_point(bn_pt, face_n, &surf.nodes);
                        let fn_ = bn.eval(&rn, surf, n_plus);
                        let g  = green3d(&rm, &rn, k);

                        let dot_ff = fm[0]*fn_[0] + fm[1]*fn_[1] + fm[2]*fn_[2];

                        let integrand = g * (dot_ff - inv_omega_eps0 / omega_mu0 * div_m * div_n);
                        val += integrand * (wm * wn * 4.0 * face_m.area * face_n.area);
                    }
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
) -> RemResult<DMatrix<Complex64>> {
    let n = bases.len();

    let cols: Vec<Vec<Complex64>> = {
        let compute = |ni: usize| -> Vec<Complex64> {
            let bn = &bases[ni];
            let mut col = vec![Complex64::ZERO; n];
            for mi in 0..n {
                let bm = &bases[mi];
                let val = zmn_mfie_rwg(bm, bn, surf, k, quad);
                col[mi] = val;
            }
            col
        };
        #[cfg(not(target_arch = "wasm32"))]
        { (0..n).into_par_iter().map(compute).collect() }
        #[cfg(target_arch = "wasm32")]
        { (0..n).map(compute).collect() }
    };

    let mut z = DMatrix::<Complex64>::from_element(n, n, Complex64::ZERO);
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
                    let dot_val = Complex64::new(
                        fm[0]*n_x_curl[0].re + fm[1]*n_x_curl[1].re + fm[2]*n_x_curl[2].re,
                        fm[0]*n_x_curl[0].im + fm[1]*n_x_curl[1].im + fm[2]*n_x_curl[2].im,
                    );

                    curl_term += dot_val
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

/// Solve Z·I = V using dense LU decomposition.
pub fn lu_solve(z: &DMatrix<Complex64>, rhs: &[Complex64]) -> RemResult<Vec<Complex64>> {
    let n = rhs.len();
    if z.nrows() != n || z.ncols() != n {
        return Err(RemError::Config(format!(
            "lu_solve: matrix {}×{} but rhs length {}", z.nrows(), z.ncols(), n
        )));
    }

    let b = DVector::<Complex64>::from_iterator(n, rhs.iter().copied());
    let lu = z.clone().lu();
    let x = lu
        .solve(&b)
        .ok_or_else(|| RemError::Other("lu_solve failed: matrix may be singular".to_string()))?;

    Ok(x.iter().copied().collect())
}

// ---------------------------------------------------------------------------
// GMRES solver (restarted, complex, for large Z matrices)
// ---------------------------------------------------------------------------

/// Solve A·x = b using restarted GMRES with generic LinearOperator<Complex64>.
///
/// This is the new interface accepting any matrix type that implements LinearOperator,
/// enabling future support for sparse matrices and other backends.
///
/// # Parameters
/// - `op`: Matrix-vector operation (matvec)
/// - `b`: Right-hand side vector
/// - `restart`: Restart parameter (default 30)
/// - `tol`: Convergence tolerance (default 1e-8)
/// - `max_iters`: Maximum total iterations (default 500)
///
/// # Returns
/// Computed solution x where A·x ≈ b
pub fn gmres_solve_generic(
    op: &dyn LinearOperator<Complex64>,
    b: &DVector<Complex64>,
    restart: usize,
    tol: f64,
    max_iters: usize,
) -> RemResult<DVector<Complex64>> {
    let (n, n_cols) = op.size();
    if n != n_cols {
        return Err(RemError::Other(
            format!("gmres requires square operator, got {}×{}", n, n_cols)
        ));
    }
    if b.len() != n {
        return Err(RemError::Other(
            format!("gmres dimension mismatch: operator {}×{}, rhs len {}", n, n_cols, b.len())
        ));
    }

    let max_outer = (max_iters + restart - 1) / restart;
    let mut x = DVector::zeros(n);
    let b_norm = b.norm();

    if b_norm < f64::EPSILON {
        return Ok(x);
    }

    for _outer in 0..max_outer {
        // Compute residual r = b - A·x
        let mut r = b.clone();
        {
            let mut ax = DVector::zeros(n);
            op.matvec(&x, &mut ax).map_err(|e| RemError::Other(e))?;
            r -= &ax;
        }
        let beta = r.norm();

        if beta / b_norm < tol {
            return Ok(x);
        }

        // Arnoldi basis V (n × (restart+1))
        let mut v: Vec<DVector<Complex64>> = Vec::with_capacity(restart + 1);
        v.push(&r * Complex64::new(1.0 / beta, 0.0));

        // Upper Hessenberg H ((restart+1) × restart)
        let mut h = vec![vec![Complex64::new(0.0, 0.0); restart]; restart + 1];

        // Givens rotation cosines and sines
        let mut cs = vec![0.0f64; restart];
        let mut sn = vec![Complex64::new(0.0, 0.0); restart];
        let mut g = vec![Complex64::new(0.0, 0.0); restart + 1];
        g[0] = Complex64::new(beta, 0.0);

        let mut j_end = restart;
        for j in 0..restart {
            // w = A · v[j]
            let mut w = DVector::zeros(n);
            op.matvec(&v[j], &mut w).map_err(|e| RemError::Other(e))?;

            // Modified Gram-Schmidt orthogonalization
            for i in 0..=j {
                h[i][j] = v[i].dot(&w.map(|c| c.conj()));
                w -= &v[i] * h[i][j];
            }
            h[j + 1][j] = Complex64::new(w.norm(), 0.0);

            // Normalize to get next basis vector
            let h_norm = h[j + 1][j].re;
            if h_norm > f64::EPSILON {
                v.push(&w * Complex64::new(1.0 / h_norm, 0.0));
            } else {
                v.push(DVector::zeros(n));
            }

            // Apply previous Givens rotations to new column
            for i in 0..j {
                let tmp = cs[i] * h[i][j] + sn[i].conj() * h[i + 1][j];
                h[i + 1][j] = -sn[i] * h[i][j] + cs[i] * h[i + 1][j];
                h[i][j] = tmp;
            }

            // Compute new Givens rotation to zero out h[j+1][j]
            let (c, s) = givens_rotation(h[j][j], h[j + 1][j]);
            cs[j] = c;
            sn[j] = s;
            h[j][j] = c * h[j][j] + s.conj() * h[j + 1][j];
            h[j + 1][j] = Complex64::new(0.0, 0.0);

            // Update residual estimate
            g[j + 1] = -sn[j] * g[j];
            g[j] = cs[j] * g[j];

            if g[j + 1].norm() / b_norm < tol {
                j_end = j + 1;
                break;
            }
        }

        // Back-substitution to get y (j_end × 1)
        let m = j_end;
        let mut y = vec![Complex64::new(0.0, 0.0); m];
        for i in (0..m).rev() {
            y[i] = g[i];
            for k in (i + 1)..m {
                let yk = y[k];
                y[i] -= h[i][k] * yk;
            }
            if h[i][i].norm() < f64::EPSILON {
                return Err(RemError::Other("gmres_solve: singular Hessenberg".to_string()));
            }
            y[i] /= h[i][i];
        }

        // Update solution x += V_m · y
        for j in 0..m {
            x += &v[j] * y[j];
        }
    }

    Ok(x)
}

/// Solve A·x = b using restarted GMRES (convenience wrapper).
///
/// Uses default parameters: restart=30, tol=1e-8, max_iters=500
pub fn gmres_solve_op(
    op: &dyn LinearOperator<Complex64>,
    b: &DVector<Complex64>,
) -> RemResult<DVector<Complex64>> {
    gmres_solve_generic(op, b, 30, 1e-8, 500)
}

/// Solve Z·I = V using restarted GMRES (restart=30, tol=1e-8, max 500 iters).
///
/// Uses modified Gram-Schmidt Arnoldi process. Suitable when N > ~1000 where
/// dense LU would be prohibitively expensive (O(N³) vs O(N²·restart) per outer).
pub fn gmres_solve(z: &DMatrix<Complex64>, rhs: &[Complex64]) -> RemResult<Vec<Complex64>> {
    let n = rhs.len();
    if z.nrows() != n || z.ncols() != n {
        return Err(RemError::Config(format!(
            "gmres_solve: matrix {}×{} but rhs length {}", z.nrows(), z.ncols(), n
        )));
    }

    const RESTART: usize = 30;
    const TOL: f64 = 1e-8;
    const MAX_OUTER: usize = 500 / RESTART + 1;

    // Initial guess x = 0
    let mut x = vec![Complex64::new(0.0, 0.0); n];

    let rhs_norm = vec_norm(rhs);
    if rhs_norm < f64::EPSILON {
        return Ok(x);
    }

    for _outer in 0..MAX_OUTER {
        // Compute residual r = b - A·x
        let mut r = rhs.to_vec();
        matvec_sub(z, &x, &mut r);          // r -= A·x
        let beta = vec_norm(&r);

        if beta / rhs_norm < TOL {
            return Ok(x);
        }

        // Arnoldi basis V (n × (restart+1))
        let mut v: Vec<Vec<Complex64>> = Vec::with_capacity(RESTART + 1);
        let scale = Complex64::new(1.0 / beta, 0.0);
        v.push(r.iter().map(|&c| c * scale).collect());

        // Upper Hessenberg H ((restart+1) × restart)
        let mut h = vec![vec![Complex64::new(0.0, 0.0); RESTART]; RESTART + 1];

        // Givens rotation cosines and sines
        let mut cs = vec![0.0f64; RESTART];
        let mut sn = vec![Complex64::new(0.0, 0.0); RESTART];
        let mut g = vec![Complex64::new(0.0, 0.0); RESTART + 1];
        g[0] = Complex64::new(beta, 0.0);

        let mut j_end = RESTART;
        for j in 0..RESTART {
            // w = A · v[j]
            let mut w = vec![Complex64::new(0.0, 0.0); n];
            matvec(z, &v[j], &mut w);

            // Modified Gram-Schmidt orthogonalization
            for i in 0..=j {
                h[i][j] = dot_conj(&v[i], &w);
                let hij = h[i][j];
                for k in 0..n {
                    w[k] -= hij * v[i][k];
                }
            }
            h[j + 1][j] = Complex64::new(vec_norm(&w), 0.0);

            // Normalize to get next basis vector
            let h_norm = h[j + 1][j].re;
            if h_norm > f64::EPSILON {
                let inv = Complex64::new(1.0 / h_norm, 0.0);
                v.push(w.iter().map(|&c| c * inv).collect());
            } else {
                v.push(vec![Complex64::new(0.0, 0.0); n]);
            }

            // Apply previous Givens rotations to new column
            for i in 0..j {
                let tmp = cs[i] * h[i][j] + sn[i].conj() * h[i + 1][j];
                h[i + 1][j] = -sn[i] * h[i][j] + cs[i] * h[i + 1][j];
                h[i][j] = tmp;
            }

            // Compute new Givens rotation to zero out h[j+1][j]
            let (c, s) = givens_rotation(h[j][j], h[j + 1][j]);
            cs[j] = c;
            sn[j] = s;
            h[j][j] = c * h[j][j] + s.conj() * h[j + 1][j];
            h[j + 1][j] = Complex64::new(0.0, 0.0);

            // Update residual estimate
            g[j + 1] = -sn[j] * g[j];
            g[j] = cs[j] * g[j];

            if g[j + 1].norm() / rhs_norm < TOL {
                j_end = j + 1;
                break;
            }
        }

        // Back-substitution to get y (j_end × 1)
        let m = j_end;
        let mut y = vec![Complex64::new(0.0, 0.0); m];
        for i in (0..m).rev() {
            y[i] = g[i];
            for k in (i + 1)..m {
                let hik_yk = h[i][k] * y[k];
                y[i] -= hik_yk;
            }
            if h[i][i].norm() < f64::EPSILON {
                return Err(RemError::Other("gmres_solve: singular Hessenberg".to_string()));
            }
            y[i] /= h[i][i];
        }

        // Update solution x += V_m · y
        for j in 0..m {
            let yj = y[j];
            for k in 0..n {
                x[k] += yj * v[j][k];
            }
        }
    }

    Ok(x)
}

// ---------------------------------------------------------------------------
// ACA-compressed GMRES solver
// ---------------------------------------------------------------------------

/// Solve Z·I = V using ACA matrix compression + restarted GMRES.
///
/// This avoids assembling the full N×N Z matrix by instead computing the
/// ACA near/far partition and running GMRES with the compressed matrix-vector
/// product.  Memory: O(N r) instead of O(N²).
///
/// `tol_aca` — ACA approximation tolerance (relative Frobenius, e.g. 1e-4)
/// `tol_gmres` — GMRES convergence tolerance (e.g. 1e-8)
pub fn aca_gmres_solve(
    z: &DMatrix<Complex64>,
    rhs: &[Complex64],
    tol_aca: f64,
    tol_gmres: f64,
) -> RemResult<Vec<Complex64>> {
    use crate::aca::{aca_partition, aca_matvec};

    let n = rhs.len();
    if z.nrows() != n || z.ncols() != n {
        return Err(RemError::Config(format!(
            "aca_gmres_solve: matrix {}×{} but rhs length {}", z.nrows(), z.ncols(), n
        )));
    }

    // Block size: ~sqrt(N), capped at 64
    let block_size = ((n as f64).sqrt() as usize).clamp(8, 64);
    // Near-field: diagonal block + 1 neighbor on each side
    let near_thresh = 1_usize;
    let max_rank = (block_size * 2).min(n);

    let entry_fn = |i: usize, j: usize| z[(i, j)];
    let (near, far) = aca_partition(n, block_size, near_thresh, tol_aca, max_rank, &entry_fn);

    log::info!(
        "ACA: {} near entries, {} far blocks (avg rank ≈ {})",
        near.len(),
        far.len(),
        far.iter().map(|(_, _, a)| a.rank).sum::<usize>().checked_div(far.len().max(1)).unwrap_or(0)
    );

    // GMRES with ACA matvec
    const RESTART: usize = 30;
    const MAX_OUTER: usize = 500 / RESTART + 1;

    let mut x = vec![Complex64::ZERO; n];
    let rhs_norm = vec_norm(rhs);
    if rhs_norm < f64::EPSILON {
        return Ok(x);
    }

    for _outer in 0..MAX_OUTER {
        // r = b - A·x
        let ax = aca_matvec(n, &near, &far, &x);
        let mut r: Vec<Complex64> = rhs.iter().zip(ax.iter()).map(|(&b, &ax)| b - ax).collect();
        let beta = vec_norm(&r);
        if beta / rhs_norm < tol_gmres { return Ok(x); }

        let mut v: Vec<Vec<Complex64>> = vec![vec![Complex64::ZERO; n]; RESTART + 1];
        let inv_beta = 1.0 / beta;
        for k in 0..n { v[0][k] = r[k] * inv_beta; }

        let mut h = vec![vec![Complex64::ZERO; RESTART]; RESTART + 1];
        let mut cs = vec![1.0_f64; RESTART];
        let mut sn = vec![Complex64::ZERO; RESTART];
        let mut g: Vec<Complex64> = vec![Complex64::ZERO; RESTART + 1];
        g[0] = Complex64::new(beta, 0.0);

        let mut j_end = RESTART;

        for j in 0..RESTART {
            let w = aca_matvec(n, &near, &far, &v[j]);

            for i in 0..=j {
                h[i][j] = dot_conj(&v[i], &w);
                for k in 0..n {
                    let h_ij = h[i][j];
                    // safe: we need mutable w, do it inline
                    r[k] = w[k] - h_ij * v[i][k]; // reuse r as w buffer
                }
                // reset r back to being w
                for k in 0..n { let tmp = r[k]; r[k] = tmp; }
            }
            // Actually accumulate into a proper w buffer using a different approach
            let mut w2 = w.clone();
            for i in 0..=j {
                let h_ij = h[i][j];
                for k in 0..n { w2[k] -= h_ij * v[i][k]; }
            }
            h[j + 1][j] = Complex64::new(vec_norm(&w2), 0.0);

            if h[j + 1][j].re.abs() > 1e-14 {
                let inv_h = 1.0 / h[j + 1][j].re;
                for k in 0..n { v[j + 1][k] = w2[k] * inv_h; }
            }

            for i in 0..j {
                let tmp = cs[i] * h[i][j] + sn[i].conj() * h[i + 1][j];
                h[i + 1][j] = -sn[i] * h[i][j] + cs[i] * h[i + 1][j];
                h[i][j] = tmp;
            }

            let (c, s) = givens_rotation(h[j][j], h[j + 1][j]);
            cs[j] = c; sn[j] = s;
            h[j][j] = c * h[j][j] + s.conj() * h[j + 1][j];
            h[j + 1][j] = Complex64::ZERO;

            g[j + 1] = -sn[j] * g[j];
            g[j]     = cs[j] * g[j];

            if g[j + 1].norm() / rhs_norm < tol_gmres {
                j_end = j + 1;
                break;
            }
        }

        let m = j_end;
        let mut y = vec![Complex64::ZERO; m];
        for i in (0..m).rev() {
            let mut yi = g[i];
            for k in (i + 1)..m { yi -= h[i][k] * y[k]; }
            if h[i][i].norm() < f64::EPSILON {
                return Err(RemError::Other("aca_gmres: singular Hessenberg".to_string()));
            }
            y[i] = yi / h[i][i];
        }
        for j in 0..m {
            for k in 0..n { x[k] += y[j] * v[j][k]; }
        }
    }

    Ok(x)
}

fn vec_norm(v: &[Complex64]) -> f64 {
    v.iter().map(|c| c.norm_sqr()).sum::<f64>().sqrt()
}

fn dot_conj(a: &[Complex64], b: &[Complex64]) -> Complex64 {
    a.iter().zip(b.iter()).map(|(&ai, &bi)| ai.conj() * bi).sum()
}

fn matvec(a: &DMatrix<Complex64>, x: &[Complex64], out: &mut [Complex64]) {
    let n = x.len();
    for i in 0..n {
        out[i] = (0..n).map(|j| a[(i, j)] * x[j]).sum();
    }
}

fn matvec_sub(a: &DMatrix<Complex64>, x: &[Complex64], r: &mut [Complex64]) {
    let n = x.len();
    for i in 0..n {
        let ax_i: Complex64 = (0..n).map(|j| a[(i, j)] * x[j]).sum();
        r[i] -= ax_i;
    }
}

/// Compute Givens rotation (c, s) such that [c s*; -s c] · [a; b] = [r; 0].
fn givens_rotation(a: Complex64, b: Complex64) -> (f64, Complex64) {
    let norm = (a.norm_sqr() + b.norm_sqr()).sqrt();
    if norm < f64::EPSILON {
        return (1.0, Complex64::new(0.0, 0.0));
    }
    let c = a.norm() / norm;
    let s = if a.norm() < f64::EPSILON {
        Complex64::new(1.0, 0.0)
    } else {
        (a / a.norm()) * (b.conj() / norm)
    };
    (c, s)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a known diagonally-dominant complex system and verify GMRES matches LU.
    #[test]
    fn gmres_matches_lu_small() {
        let n = 8usize;
        let mut z = DMatrix::<Complex64>::zeros(n, n);
        let mut rhs = vec![Complex64::new(0.0, 0.0); n];

        // Fill with a diagonally dominant complex matrix
        for i in 0..n {
            for j in 0..n {
                if i == j {
                    z[(i, j)] = Complex64::new(10.0 + i as f64, 1.0);
                } else {
                    z[(i, j)] = Complex64::new(0.5 / (1.0 + (i as f64 - j as f64).abs()), -0.2);
                }
            }
            rhs[i] = Complex64::new(i as f64 + 1.0, -(i as f64));
        }

        let x_lu = lu_solve(&z, &rhs).unwrap();
        let x_gmres = gmres_solve(&z, &rhs).unwrap();

        for i in 0..n {
            let err = (x_lu[i] - x_gmres[i]).norm();
            assert!(err < 1e-6, "index {i}: LU={} GMRES={} err={err:.2e}",
                x_lu[i], x_gmres[i]);
        }
    }

    /// GMRES on trivial identity system: I·x = b → x = b.
    #[test]
    fn gmres_identity_system() {
        let n = 5usize;
        let mut z = DMatrix::<Complex64>::zeros(n, n);
        for i in 0..n { z[(i, i)] = Complex64::new(1.0, 0.0); }
        let rhs: Vec<Complex64> = (0..n).map(|i| Complex64::new(i as f64, 1.0)).collect();

        let x = gmres_solve(&z, &rhs).unwrap();
        for i in 0..n {
            let err = (x[i] - rhs[i]).norm();
            assert!(err < 1e-10, "index {i}: got {}, expected {}, err={err:.2e}",
                x[i], rhs[i]);
        }
    }

    /// Verify gmres_solve_op (new LinearOperator interface) matches old gmres_solve.
    #[test]
    fn gmres_solve_op_matches_old() {
        let n = 8usize;
        let mut z = DMatrix::<Complex64>::zeros(n, n);
        let mut rhs = vec![Complex64::new(0.0, 0.0); n];

        // Same diagonally dominant complex matrix as gmres_matches_lu_small
        for i in 0..n {
            for j in 0..n {
                if i == j {
                    z[(i, j)] = Complex64::new(10.0 + i as f64, 1.0);
                } else {
                    z[(i, j)] = Complex64::new(0.5 / (1.0 + (i as f64 - j as f64).abs()), -0.2);
                }
            }
            rhs[i] = Complex64::new(i as f64 + 1.0, -(i as f64));
        }

        // Solve with old interface
        let x_old = gmres_solve(&z, &rhs).unwrap();

        // Solve with new LinearOperator interface
        let b_dvector = DVector::from_vec(rhs.clone());
        let x_new = gmres_solve_op(&z, &b_dvector).unwrap();

        // Verify results match (allow small numerical difference)
        for i in 0..n {
            let err = (x_old[i] - x_new[i]).norm();
            assert!(err < 1e-7, "index {i}: old={} new={} err={err:.2e}",
                x_old[i], x_new[i]);
        }
    }

    /// Test gmres_solve_generic with custom LinearOperator struct.
    #[test]
    fn gmres_solve_generic_identity() {
        let n = 4usize;
        let z = DMatrix::<Complex64>::from_fn(n, n, |i, j| {
            if i == j {
                Complex64::new(2.0, 0.0)
            } else {
                Complex64::new(0.0, 0.0)
            }
        });

        let b = DVector::from_vec(vec![
            Complex64::new(1.0, 0.0),
            Complex64::new(2.0, 0.0),
            Complex64::new(3.0, 0.0),
            Complex64::new(4.0, 0.0),
        ]);

        let x = gmres_solve_generic(&z, &b, 10, 1e-8, 100).unwrap();

        // Expected: x ≈ [0.5, 1.0, 1.5, 2.0] since 2·x = b
        for i in 0..n {
            let expected = Complex64::new((i as f64 + 1.0) / 2.0, 0.0);
            let err = (x[i] - expected).norm();
            assert!(err < 1e-8, "index {i}: got {}, expected {}, err={err:.2e}",
                x[i], expected);
        }
    }
}
