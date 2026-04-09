//! Vector Fitting (VF) for rational approximation of S-parameter data.
//!
//! Given complex frequency-response samples H(jω_k), fits a rational function:
//!
//!   H(s) ≈ Σ_{k=1}^{n} R_k / (s - p_k)  +  D
//!
//! where p_k are poles (stable: Re(p_k) < 0) and R_k are residues.
//! Complex-conjugate poles are stored as a single entry; the conjugate pair
//! R*/(s-p*) is implied and added during evaluation.
//!
//! # Algorithm
//! Gustavsen & Semlyen (1999) "Rational Approximation of Frequency Domain
//! Responses by Vector Fitting".  Key steps:
//! 1. Place initial poles logarithmically over the frequency band.
//! 2. Build a real-valued LS matrix (2N×2n_poles+2) in the partial-fraction basis.
//! 3. Solve via QR → denominator coefficients → extract new poles from Schur
//!    decomposition of the companion matrix.
//! 4. Flip unstable poles (positive real part) to stable.
//! 5. Repeat until pole movement < tol (typically 3–10 iterations).
//! 6. Final LS to fit residues given the converged poles.
//!
//! # Outputs
//! - Touchstone `.s1p` (RI format, 50 Ω)
//! - `circuit_model.csv` (pole-residue table)
//! - `equivalent_circuit.cir` (SPICE Laplace-controlled-source netlist)

use nalgebra::{DMatrix, DVector};
use num_complex::Complex64;
use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// Public data structures
// ─────────────────────────────────────────────────────────────────────────────

/// Pole-residue model H(s) ≈ Σ R_k/(s-p_k) + D.
/// Complex conjugate pairs are stored as a single entry:
/// the term is  R/(s-p) + R*/(s-p*)  evaluated by `eval_vf`.
#[derive(Debug, Clone)]
pub struct VfModel {
    /// Poles p_k.  One entry per logical pole; conjugate pairs stored once.
    pub poles: Vec<Complex64>,
    /// Residues R_k matching `poles`.
    pub residues: Vec<Complex64>,
    /// Constant (direct) term D.
    pub d_term: f64,
    /// Minimum fit frequency [Hz].
    pub f_min_hz: f64,
    /// Maximum fit frequency [Hz].
    pub f_max_hz: f64,
    /// RMS fitting error ‖H_fit − H_data‖ / N.
    pub rms_error: f64,
}

// ─────────────────────────────────────────────────────────────────────────────
// Public API
// ─────────────────────────────────────────────────────────────────────────────

/// Run Vector Fitting on complex S-parameter samples.
///
/// * `freqs_hz` — sample frequencies [Hz], length N ≥ 4
/// * `h_data`   — complex H(jω) values, length N
/// * `n_poles`  — number of poles (logical count; recommended 4–16)
/// * `max_iter` — VF outer iterations (default 10)
/// * `tol`      — convergence: max pole movement per iteration (1e-6 recommended)
///
/// Returns `None` only if the LS solve fails in every iteration.
pub fn vector_fit(
    freqs_hz: &[f64],
    h_data: &[Complex64],
    n_poles: usize,
    max_iter: usize,
    tol: f64,
) -> Option<VfModel> {
    let n = freqs_hz.len();
    if n < 4 || n_poles == 0 { return None; }

    let omegas: Vec<f64> = freqs_hz.iter().map(|&f| 2.0 * PI * f).collect();
    let mut poles = initial_poles(*freqs_hz.first()?, *freqs_hz.last()?, n_poles);

    for _iter in 0..max_iter {
        let (a_mat, b_vec) = build_vf_ls(&omegas, h_data, &poles);
        let svd = a_mat.svd(true, true);
        let x = match svd.solve(&b_vec, 1e-10) {
            Ok(v) => v,
            Err(_) => continue,
        };
        // Split x into numerator half and denominator half.
        // Each complex pair contributes 2 columns, each real pole 1 column.
        let n_num_cols = num_scalar_cols(&poles) + 1; // +1 for D_num constant
        let den_start  = n_num_cols;
        let n_den_cols = num_scalar_cols(&poles);     // denominator coefficients (no D_den constant)
        let den_coeffs: Vec<f64> = x.iter().skip(den_start).take(n_den_cols).copied().collect();

        let new_poles = match extract_poles(&poles, &den_coeffs) {
            Some(p) => p,
            None => continue,
        };
        // Flip unstable poles
        let new_poles: Vec<Complex64> = new_poles.into_iter()
            .map(|p| if p.re > 0.0 { Complex64::new(-p.re, p.im) } else { p })
            .collect();

        // Convergence check
        let max_move = poles.iter().zip(new_poles.iter())
            .map(|(po, pn)| (pn - po).norm())
            .fold(0.0_f64, f64::max);

        poles = new_poles;
        if max_move < tol { break; }
    }

    let (residues, d_term) = fit_residues(&omegas, h_data, &poles)?;
    let rms_error = compute_rms(&omegas, h_data, &poles, &residues, d_term);

    Some(VfModel {
        poles,
        residues,
        d_term,
        f_min_hz: *freqs_hz.first()?,
        f_max_hz: *freqs_hz.last()?,
        rms_error,
    })
}

/// Evaluate the VfModel at the given frequencies.
pub fn eval_vf(model: &VfModel, freqs_hz: &[f64]) -> Vec<Complex64> {
    freqs_hz.iter().map(|&f| {
        let s = Complex64::new(0.0, 2.0 * PI * f);
        let sum: Complex64 = model.poles.iter().zip(model.residues.iter())
            .map(|(&p, &r)| {
                if p.im.abs() > 1e-10 {
                    // Complex conjugate pair: R/(s-p) + R*/(s-p*)
                    r / (s - p) + r.conj() / (s - p.conj())
                } else {
                    r / (s - p)
                }
            })
            .sum();
        sum + Complex64::new(model.d_term, 0.0)
    }).collect()
}

// ─────────────────────────────────────────────────────────────────────────────
// Output format functions
// ─────────────────────────────────────────────────────────────────────────────

/// Write single-port Touchstone (.s1p) in RI format.
pub fn write_touchstone_s1p(freqs_hz: &[f64], s11: &[Complex64], z0_ohm: f64) -> String {
    let mut out = String::new();
    out.push_str("! Touchstone S1P generated by rem2 EM solver\n");
    out.push_str(&format!("# GHz S RI R {:.1}\n", z0_ohm));
    out.push_str("! Freq(GHz)   Re(S11)   Im(S11)\n");
    for (f, s) in freqs_hz.iter().zip(s11.iter()) {
        out.push_str(&format!(
            "{:.9e}   {:.9e}   {:.9e}\n",
            f / 1e9, s.re, s.im
        ));
    }
    out
}

/// Write pole-residue model as CSV.
/// Complex pairs are written as two rows (pole and conjugate).
/// Metadata rows at the end (d_term, rms_error, frequency range).
pub fn write_circuit_model_csv(model: &VfModel) -> String {
    let mut out = String::new();
    out.push_str("pole_re,pole_im,residue_re,residue_im,pole_type\n");
    for (p, r) in model.poles.iter().zip(model.residues.iter()) {
        if p.im.abs() > 1e-10 {
            out.push_str(&format!(
                "{:.9e},{:.9e},{:.9e},{:.9e},complex_pair\n",
                p.re, p.im, r.re, r.im
            ));
            out.push_str(&format!(
                "{:.9e},{:.9e},{:.9e},{:.9e},complex_pair\n",
                p.re, -p.im, r.re, -r.im
            ));
        } else {
            out.push_str(&format!(
                "{:.9e},{:.9e},{:.9e},{:.9e},real\n",
                p.re, 0.0_f64, r.re, r.im
            ));
        }
    }
    out.push_str(&format!("d_term,{:.9e},,,\n", model.d_term));
    out.push_str(&format!("rms_error,{:.6e},,,\n", model.rms_error));
    out.push_str(&format!("f_min_hz,{:.6e},,,\n", model.f_min_hz));
    out.push_str(&format!("f_max_hz,{:.6e},,,\n", model.f_max_hz));
    out
}

/// Write SPICE netlist using Laplace-controlled-source syntax (ngspice compatible).
/// Each pole-residue term is expressed in partial-fraction form inside LAPLACE.
pub fn write_spice_netlist(model: &VfModel, z0_ohm: f64) -> String {
    // Build Laplace expression as sum of partial fractions.
    // Complex pair (σ, ω), residue (a+jb):
    //   R/(s-p) + R*/(s-p*) = [2a(s-σ) - 2bω] / [(s-σ)² + ω²]
    //                       = [2a*(s+|σ|) - 2b*ω] / [(s+|σ|)² + ω²]   (σ < 0)
    // Real pole σ_r, residue a:
    //   a / (s - σ_r) = a / (s + |σ_r|)   (σ_r < 0)
    let mut terms: Vec<String> = Vec::new();

    for (p, r) in model.poles.iter().zip(model.residues.iter()) {
        if p.im.abs() > 1e-10 {
            let sigma = -p.re;  // positive damping
            let omega = p.im.abs();
            let a = 2.0 * r.re;
            let b = if p.im > 0.0 { -2.0 * r.im } else { 2.0 * r.im };
            // Numerator: a*(s+sigma) + b*omega
            // Denominator: (s+sigma)^2 + omega^2
            terms.push(format!(
                "({a:.6e}*(s+{sigma:.6e})+{b:.6e}*{omega:.6e})/((s+{sigma:.6e})^2+({omega:.6e})^2)",
                a = a, sigma = sigma, b = b, omega = omega
            ));
        } else {
            let sigma = -p.re;  // positive
            let a = r.re;
            terms.push(format!("{a:.6e}/(s+{sigma:.6e})", a = a, sigma = sigma));
        }
    }
    if model.d_term.abs() > 1e-15 {
        terms.push(format!("{:.6e}", model.d_term));
    }

    let laplace_expr = if terms.is_empty() {
        "0".to_string()
    } else {
        terms.join("+")
    };

    format!(
        "* S11 equivalent circuit generated by rem2 EM solver\n\
         * H(s) = sum_k R_k/(s-p_k) + D  (partial-fraction form)\n\
         * Fit range: {f_min:.3e} Hz to {f_max:.3e} Hz\n\
         * RMS fitting error: {rms:.4e}\n\
         *\n\
         .SUBCKT S11_MODEL IN GND\n\
         R_port IN GND {z0:.1}\n\
         * S11(s) output on node SOUT (re: GND)\n\
         E_S11 SOUT GND LAPLACE {{V(IN,GND)}} = {{{expr}}}\n\
         .ENDS S11_MODEL\n",
        f_min = model.f_min_hz,
        f_max = model.f_max_hz,
        rms   = model.rms_error,
        z0    = z0_ohm,
        expr  = laplace_expr,
    )
}

// ─────────────────────────────────────────────────────────────────────────────
// Internal helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Count scalar columns needed for a set of poles.
/// Complex pair → 2 columns; real pole → 1 column.
fn num_scalar_cols(poles: &[Complex64]) -> usize {
    poles.iter().map(|p| if p.im.abs() > 1e-10 { 2 } else { 1 }).sum()
}

/// Initial pole placement: n_complex = n/2 complex pairs + n_real = n%2 real poles.
fn initial_poles(f_min: f64, f_max: f64, n: usize) -> Vec<Complex64> {
    let omega_min = 2.0 * PI * f_min;
    let omega_max = 2.0 * PI * f_max;
    let n_pairs = n / 2;
    let n_real  = n % 2;
    let mut poles = Vec::with_capacity(n_pairs + n_real);

    for k in 0..n_pairs {
        let t = if n_pairs > 1 { k as f64 / (n_pairs - 1) as f64 } else { 0.5 };
        let omega = if f_max > f_min * 10.0 {
            // log spacing
            omega_min * (omega_max / omega_min).powf(t)
        } else {
            omega_min + t * (omega_max - omega_min)
        };
        let sigma = -0.01 * omega;
        poles.push(Complex64::new(sigma, omega));
    }
    if n_real > 0 {
        let omega_mid = (omega_min * omega_max).sqrt();
        poles.push(Complex64::new(-0.01 * omega_mid, 0.0));
    }
    poles
}

/// Build the real-valued VF least-squares matrix (2N × 2*n_scalar_cols+2).
///
/// Column layout (n_sc = num_scalar_cols):
///   [0 .. n_sc)           : numerator basis Φ_num
///   [n_sc]                : constant column (D_num)
///   [n_sc+1 .. 2*n_sc+1)  : denominator basis -Re/Im(H·Φ_den) interleaved
///
/// Row layout (row 2k = Re equation at ω_k, row 2k+1 = Im equation):
///   A[2k, col]   = Re( Φ_col(jω_k) )
///   A[2k+1, col] = Im( Φ_col(jω_k) )
///
/// The full equation is:
///   Φ_num·c_num + D_num - H·(Φ_den·c_den) = H  (complex)
/// Split into Re/Im → real LS system, RHS = [Re(H); Im(H)] interleaved.
fn build_vf_ls(
    omegas: &[f64],
    h_data: &[Complex64],
    poles: &[Complex64],
) -> (DMatrix<f64>, DVector<f64>) {
    let n_sc = num_scalar_cols(poles);
    let n_rows = 2 * omegas.len();
    let n_cols = 2 * n_sc + 2; // num_basis + D_num + den_basis  (no D_den: absorbed into D_num)
    // Actually standard VF includes D_den too → n_cols = 2*n_sc + 2
    // We keep: [Φ_num (n_sc) | D_num (1) | -H·Φ_den (n_sc) | -H·D_den? no]
    // Simplified: drop the D_den term (set denominator constant = 1 implicitly).
    // This is the standard "fixed denominator constant" form.

    let mut a = DMatrix::<f64>::zeros(n_rows, n_cols);
    let mut b = DVector::<f64>::zeros(n_rows);

    for (k, (&omega_k, &hk)) in omegas.iter().zip(h_data.iter()).enumerate() {
        let row_re = 2 * k;
        let row_im = 2 * k + 1;
        b[row_re] = hk.re;
        b[row_im] = hk.im;

        let mut col = 0usize;
        for pole in poles.iter() {
            if pole.im.abs() > 1e-10 {
                // Complex pair p = σ ± jω_p
                let sigma = pole.re;
                let omega_p = pole.im;
                let denom = sigma * sigma + (omega_k - omega_p) * (omega_k - omega_p);
                let b_re = 2.0 * sigma / denom;
                let b_im = -2.0 * (omega_k - omega_p) / denom;

                // Numerator columns
                a[(row_re, col)]     = b_re;
                a[(row_im, col)]     = b_im;
                a[(row_re, col + 1)] = b_im;   // note: Im part of B swapped for Re row
                a[(row_im, col + 1)] = -b_re;

                // Denominator columns: -Re(H * B), -Im(H * B)
                // H * B_complex = (h_re + j*h_im)(b_re + j*b_im)
                //               = (h_re*b_re - h_im*b_im) + j(h_re*b_im + h_im*b_re)
                let hb_re = hk.re * b_re - hk.im * b_im;
                let hb_im = hk.re * b_im + hk.im * b_re;
                let den_col = n_sc + 1 + col;
                a[(row_re, den_col)]     = -hb_re;
                a[(row_im, den_col)]     = -hb_im;
                a[(row_re, den_col + 1)] = -hb_im;   // cross terms for Im equation
                a[(row_im, den_col + 1)] =  hb_re;

                col += 2;
            } else {
                // Real pole p = σ_r (purely real, σ_r < 0)
                let sigma_r = pole.re;
                let denom = sigma_r * sigma_r + omega_k * omega_k;
                let phi_re = -sigma_r / denom;   // Re(1/(jΩ - σ_r))
                let phi_im = -omega_k / denom;   // Im(1/(jΩ - σ_r))

                a[(row_re, col)] = phi_re;
                a[(row_im, col)] = phi_im;

                let hp_re = hk.re * phi_re - hk.im * phi_im;
                let hp_im = hk.re * phi_im + hk.im * phi_re;
                let den_col = n_sc + 1 + col;
                a[(row_re, den_col)] = -hp_re;
                a[(row_im, den_col)] = -hp_im;

                col += 1;
            }
        }
        // D_num constant column (index n_sc)
        a[(row_re, n_sc)] = 1.0;
        a[(row_im, n_sc)] = 0.0;
    }

    (a, b)
}

/// Extract new poles from the denominator coefficients via Schur decomposition.
///
/// New poles = eigenvalues of  H_companion = diag(poles) - 1·c_den^T
/// where c_den are the denominator coefficients from the LS solution.
fn extract_poles(
    poles: &[Complex64],
    den_coeffs: &[f64],
) -> Option<Vec<Complex64>> {
    let n_sc = num_scalar_cols(poles);
    if den_coeffs.len() != n_sc { return None; }

    // Build real companion matrix of size n_sc × n_sc.
    // Diagonal: for complex pair (σ, ω) → 2×2 block [[σ, ω], [-ω, σ]]
    //           for real pole σ_r → scalar σ_r
    // Off-diagonal: -b_i * c_j  where b_i = 1 for real, and for complex pair block:
    //   each column j gets  -b * [c_j; c_{j+1}] for the pair
    // Simplified: use b_i = 1 for all scalar columns.

    let mut companion = DMatrix::<f64>::zeros(n_sc, n_sc);

    // Fill diagonal blocks
    let mut col = 0usize;
    for pole in poles.iter() {
        if pole.im.abs() > 1e-10 {
            companion[(col,     col)]     =  pole.re;
            companion[(col,     col + 1)] =  pole.im;
            companion[(col + 1, col)]     = -pole.im;
            companion[(col + 1, col + 1)] =  pole.re;
            col += 2;
        } else {
            companion[(col, col)] = pole.re;
            col += 1;
        }
    }

    // Fill -b * c^T rank-1 update (b = ones vector)
    for i in 0..n_sc {
        for j in 0..n_sc {
            companion[(i, j)] -= den_coeffs[j]; // b_i = 1
        }
    }

    // Schur decomposition → real Schur form (quasi-upper-triangular)
    let schur = nalgebra::linalg::Schur::new(companion);
    let (t, _q) = schur.unpack();

    extract_eigenvalues_from_schur(&t, n_sc)
}

/// Read eigenvalues from real Schur form T (quasi-upper-triangular).
/// 1×1 diagonal blocks → real eigenvalue.
/// 2×2 diagonal blocks with sub-diagonal entry → complex conjugate pair.
fn extract_eigenvalues_from_schur(t: &DMatrix<f64>, n: usize) -> Option<Vec<Complex64>> {
    let mut poles: Vec<Complex64> = Vec::with_capacity(n);
    let mut i = 0usize;
    while i < n {
        if i + 1 < n && t[(i + 1, i)].abs() > 1e-10 {
            // 2×2 block: complex conjugate pair
            let a = t[(i,     i)];
            let b = t[(i,     i + 1)];
            let c = t[(i + 1, i)];
            let d = t[(i + 1, i + 1)];
            let tr   = a + d;
            let det  = a * d - b * c;
            let disc = tr * tr - 4.0 * det;
            if disc < 0.0 {
                let sigma = tr / 2.0;
                let omega = (-disc).sqrt() / 2.0;
                poles.push(Complex64::new(sigma, omega)); // store only one of the pair
            } else {
                // Two real eigenvalues from the block
                let sq = disc.sqrt();
                poles.push(Complex64::new((tr + sq) / 2.0, 0.0));
                poles.push(Complex64::new((tr - sq) / 2.0, 0.0));
            }
            i += 2;
        } else {
            // 1×1 block: real eigenvalue
            poles.push(Complex64::new(t[(i, i)], 0.0));
            i += 1;
        }
    }
    Some(poles)
}

/// Final residue fit: given converged poles, solve for residues and D term.
/// Returns (residues, d_term).
fn fit_residues(
    omegas: &[f64],
    h_data: &[Complex64],
    poles: &[Complex64],
) -> Option<(Vec<Complex64>, f64)> {
    let n_sc = num_scalar_cols(poles);
    let n_rows = 2 * omegas.len();
    let n_cols = n_sc + 1; // residue columns + D column

    let mut a = DMatrix::<f64>::zeros(n_rows, n_cols);
    let mut b = DVector::<f64>::zeros(n_rows);

    for (k, (&omega_k, &hk)) in omegas.iter().zip(h_data.iter()).enumerate() {
        let row_re = 2 * k;
        let row_im = 2 * k + 1;
        b[row_re] = hk.re;
        b[row_im] = hk.im;

        let mut col = 0usize;
        for pole in poles.iter() {
            if pole.im.abs() > 1e-10 {
                let sigma  = pole.re;
                let omega_p = pole.im;
                let denom = sigma * sigma + (omega_k - omega_p) * (omega_k - omega_p);
                let b_re = 2.0 * sigma / denom;
                let b_im = -2.0 * (omega_k - omega_p) / denom;

                // Re row: Re(R) * b_re - Im(R) * b_im
                // Im row: Re(R) * b_im + Im(R) * b_re
                a[(row_re, col)]     = b_re;
                a[(row_re, col + 1)] = -b_im;
                a[(row_im, col)]     = b_im;
                a[(row_im, col + 1)] = b_re;
                col += 2;
            } else {
                let sigma_r = pole.re;
                let denom = sigma_r * sigma_r + omega_k * omega_k;
                let phi_re = -sigma_r / denom;
                let phi_im = -omega_k / denom;
                a[(row_re, col)] = phi_re;
                a[(row_im, col)] = phi_im;
                col += 1;
            }
        }
        // D column
        a[(row_re, n_sc)] = 1.0;
    }

    let svd = a.svd(true, true);
    let x = svd.solve(&b, 1e-10).ok()?;

    // Reconstruct complex residues
    let mut residues: Vec<Complex64> = Vec::with_capacity(poles.len());
    let mut col = 0usize;
    for pole in poles.iter() {
        if pole.im.abs() > 1e-10 {
            residues.push(Complex64::new(x[col], x[col + 1]));
            col += 2;
        } else {
            residues.push(Complex64::new(x[col], 0.0));
            col += 1;
        }
    }
    let d_term = x[n_sc];
    Some((residues, d_term))
}

/// Compute RMS fitting error ‖H_fit − H_data‖ / N.
fn compute_rms(
    omegas: &[f64],
    h_data: &[Complex64],
    poles: &[Complex64],
    residues: &[Complex64],
    d_term: f64,
) -> f64 {
    let freqs: Vec<f64> = omegas.iter().map(|&w| w / (2.0 * PI)).collect();
    let model = VfModel {
        poles: poles.to_vec(),
        residues: residues.to_vec(),
        d_term,
        f_min_hz: freqs.first().copied().unwrap_or(0.0),
        f_max_hz: freqs.last().copied().unwrap_or(0.0),
        rms_error: 0.0,
    };
    let h_fit = eval_vf(&model, &freqs);
    let n = h_data.len() as f64;
    let sse: f64 = h_data.iter().zip(h_fit.iter())
        .map(|(a, b)| (a - b).norm_sqr())
        .sum();
    sse.sqrt() / n
}

// ─────────────────────────────────────────────────────────────────────────────
// Unit tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Fit a known single complex conjugate pair and verify the RMS is low.
    #[test]
    fn test_vf_single_pair_recovery() {
        let f_min = 0.5e9;
        let f_max = 2.0e9;
        let n_pts = 60;
        let freqs: Vec<f64> = (0..n_pts)
            .map(|i| f_min + i as f64 * (f_max - f_min) / (n_pts - 1) as f64)
            .collect();

        // True model: one complex pair p = -1e8 ± j*2π*1e9, residue R = 1
        let p0 = Complex64::new(-1e8, 2.0 * PI * 1e9);
        let r0 = Complex64::new(1.0, 0.0);
        let h_data: Vec<Complex64> = freqs.iter().map(|&f| {
            let s = Complex64::new(0.0, 2.0 * PI * f);
            r0 / (s - p0) + r0.conj() / (s - p0.conj())
        }).collect();

        let model = vector_fit(&freqs, &h_data, 2, 15, 1e-8)
            .expect("VF should converge on a single complex pair");
        assert!(
            model.rms_error < 1e-3,
            "RMS error too large: {:.4e} (expected < 1e-3)",
            model.rms_error
        );
    }

    /// Touchstone output must contain the standard option line.
    #[test]
    fn test_touchstone_format() {
        let freqs = vec![1e9, 2e9, 3e9];
        let s11 = vec![
            Complex64::new(0.8, -0.2),
            Complex64::new(0.5, -0.5),
            Complex64::new(0.3, -0.7),
        ];
        let ts = write_touchstone_s1p(&freqs, &s11, 50.0);
        assert!(ts.contains("# GHz S RI R 50.0"), "Missing Touchstone option line");
        assert!(ts.contains("1.000000000e0") || ts.contains("1.0e0") || ts.contains("1.000"),
                "Missing first frequency entry");
    }

    /// CSV output must have header and metadata footer.
    #[test]
    fn test_circuit_csv_format() {
        let model = VfModel {
            poles:    vec![Complex64::new(-1e8, 6.28e9)],
            residues: vec![Complex64::new(0.1,  0.05)],
            d_term:   0.01,
            f_min_hz: 1e9,
            f_max_hz: 5e9,
            rms_error: 0.001,
        };
        let csv = write_circuit_model_csv(&model);
        assert!(csv.starts_with("pole_re,pole_im"), "Missing CSV header");
        assert!(csv.contains("d_term"),    "Missing d_term footer");
        assert!(csv.contains("complex_pair"), "Missing pole_type label");
        assert!(csv.contains("rms_error"), "Missing rms_error footer");
    }

    /// eval_vf must reconstruct reasonable values after fitting smooth data.
    #[test]
    fn test_eval_vf_consistency() {
        // Use a simple smooth response (no sharp resonances)
        let freqs: Vec<f64> = (0..30)
            .map(|i| 1e9 + i as f64 * 1e8)
            .collect();
        // Artificial smooth S11 decreasing from 0.9 to 0.1
        let h_data: Vec<Complex64> = freqs.iter().enumerate().map(|(i, _)| {
            let t = i as f64 / 29.0;
            Complex64::new(0.9 - 0.8 * t, -0.3 * (PI * t).sin())
        }).collect();

        let model = vector_fit(&freqs, &h_data, 4, 15, 1e-6)
            .expect("VF should converge on smooth data");
        let h_fit = eval_vf(&model, &freqs);

        // Each point should be reasonably close (not exact for 4 poles on 30 samples)
        let max_err = h_data.iter().zip(h_fit.iter())
            .map(|(a, b)| (a - b).norm())
            .fold(0.0_f64, f64::max);
        assert!(max_err < 0.6, "Max pointwise error too large: {:.4}", max_err);
    }
}
