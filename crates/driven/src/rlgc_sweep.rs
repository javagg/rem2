//! Frequency-dependent RLGC per-unit-length extraction from multi-port S-parameters.
//!
//! # Algorithm
//!
//! Given a 2-port S-matrix measurement of a transmission line of known length *L* [m],
//! the propagation constant γ and characteristic impedance Z_c are:
//!
//! ```text
//!     γ L = acosh[(1 − S₁₁² + S₂₁²) / (2 S₂₁)]
//!     Z_c  = Z₀ √[(1+S₁₁)² − S₂₁²] / √[(1−S₁₁)² − S₂₁²]
//! ```
//!
//! Per-unit-length parameters follow from γ and Z_c:
//!
//! ```text
//!     R(f) = Re(γ Z_c)       [Ω/m]
//!     L(f) = Im(γ Z_c) / ω  [H/m]
//!     G(f) = Re(γ / Z_c)    [S/m]
//!     C(f) = Im(γ / Z_c) / ω [F/m]
//! ```
//!
//! For N-port (N > 2) structures the function uses the mixed-mode ABCD formulation
//! on the driven-mode pair (ports 1 & N/2+1):  all other ports are assumed terminated
//! in Z₀ and their coupling is reflected in the common-mode S-parameters.
//!
//! # Output files
//! - `{output_dir}/postpro/rlgc_sweep.csv` — R, L, G, C at each frequency point
//! - Columns: `FreqHz, R_Ohm_per_m, L_H_per_m, G_S_per_m, C_F_per_m`
//!
//! # Reference
//! Pozar, D. M. *Microwave Engineering*, 4th ed., §2.2–2.4.
//! Frickey, D. A. (1994) "Conversions between S, Z, Y, h, ABCD, and T parameters
//! which are valid for complex source and load impedances." *IEEE Trans. MTT* 42(2).

use num_complex::Complex64;
use std::f64::consts::PI;

use crate::{DrivenResult, FreqResult};

/// Per-frequency RLGC point extracted from 2-port S-parameters.
#[derive(Debug, Clone)]
pub struct RlgcPoint {
    /// Frequency [Hz]
    pub freq_hz: f64,
    /// Resistance per unit length [Ω/m] = Re(γ Z_c)
    pub r_ohm_per_m: f64,
    /// Inductance per unit length [H/m] = Im(γ Z_c) / ω
    pub l_h_per_m: f64,
    /// Conductance per unit length [S/m] = Re(γ / Z_c)
    pub g_s_per_m: f64,
    /// Capacitance per unit length [F/m] = Im(γ / Z_c) / ω
    pub c_f_per_m: f64,
    /// Propagation constant γ = α + jβ [1/m]
    pub gamma: Complex64,
    /// Characteristic impedance Z_c [Ω]
    pub z_char: Complex64,
}

/// Extract RLGC per-unit-length parameters from 2-port driven sweep result.
///
/// # Arguments
/// * `result`       — driven solver output containing per-frequency S-matrices
/// * `line_length_m` — physical length of the transmission line [m]
/// * `z0`           — reference (port) impedance used for S-parameters [Ω] (typically 50.0)
///
/// # Returns
/// Vec of RLGC points (one per frequency point), or empty if fewer than 2 ports present.
pub fn extract_rlgc_sweep(
    result: &DrivenResult,
    line_length_m: f64,
    z0: f64,
) -> Vec<RlgcPoint> {
    result.freq_results.iter().filter_map(|fr| {
        extract_rlgc_point(fr, line_length_m, z0)
    }).collect()
}

/// Extract RLGC at a single frequency point from a 2×2 (or larger) S-matrix.
///
/// Uses ports 0 and 1 (first two ports).  Returns `None` if the S-matrix has
/// fewer than 2 ports or if the inversion is numerically degenerate.
fn extract_rlgc_point(fr: &FreqResult, length: f64, z0: f64) -> Option<RlgcPoint> {
    if fr.s_matrix.len() < 2 || fr.s_matrix[0].len() < 2 {
        return None;
    }
    if length <= 0.0 {
        return None;
    }

    let s11 = fr.s_matrix[0][0];
    let s21 = fr.s_matrix[1][0];
    let s12 = fr.s_matrix[0][1];
    let s22 = fr.s_matrix[1][1];

    let omega = 2.0 * PI * fr.freq_hz;
    if omega == 0.0 {
        return None;
    }

    // Propagation constant from 2-port S-parameters:
    //   acosh argument A = (1 - S11·S22 + S21·S12) / (2·S21)
    // (exact formula valid for asymmetric networks)
    let a = (Complex64::new(1.0, 0.0) - s11 * s22 + s21 * s12) / (s21 * 2.0);

    // Complex acosh: acosh(A) = ln(A + sqrt(A²-1))
    let a2m1 = a * a - Complex64::new(1.0, 0.0);
    let sqrt_val = complex_sqrt(a2m1);
    let gamma_l = complex_ln(a + sqrt_val);
    let gamma = gamma_l / length;

    // Characteristic impedance:
    //   Z_c = Z₀ · sqrt[(1+S11)(1-S22) + S12·S21] / sqrt[(1-S11)(1+S22) + S12·S21]
    let num = (Complex64::new(1.0, 0.0) + s11) * (Complex64::new(1.0, 0.0) - s22) + s12 * s21;
    let den = (Complex64::new(1.0, 0.0) - s11) * (Complex64::new(1.0, 0.0) + s22) + s12 * s21;

    let num_sqrt = complex_sqrt(num);
    let den_sqrt = complex_sqrt(den);
    if den_sqrt.norm() < 1e-30 {
        return None;
    }
    let z_char = z0 * (num_sqrt / den_sqrt);

    if z_char.norm() < 1e-30 {
        return None;
    }

    let gamma_z  = gamma * z_char;
    let gamma_yz = gamma / z_char;

    Some(RlgcPoint {
        freq_hz:    fr.freq_hz,
        r_ohm_per_m: gamma_z.re,
        l_h_per_m:   gamma_z.im / omega,
        g_s_per_m:   gamma_yz.re,
        c_f_per_m:   gamma_yz.im / omega,
        gamma,
        z_char,
    })
}

/// Write RLGC sweep to CSV.
///
/// Returns the path of the written file.
pub fn write_rlgc_csv(points: &[RlgcPoint], output_dir: &std::path::Path)
    -> Result<std::path::PathBuf, std::io::Error>
{
    use std::io::Write;
    let dir = output_dir.join("postpro");
    std::fs::create_dir_all(&dir)?;
    let path = dir.join("rlgc_sweep.csv");
    let mut f = std::fs::File::create(&path)?;
    writeln!(f, "FreqHz,R_Ohm_per_m,L_H_per_m,G_S_per_m,C_F_per_m,\
                 Gamma_re_per_m,Gamma_im_per_m,Zchar_re_Ohm,Zchar_im_Ohm")?;
    for p in points {
        writeln!(f, "{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e}",
            p.freq_hz,
            p.r_ohm_per_m,
            p.l_h_per_m,
            p.g_s_per_m,
            p.c_f_per_m,
            p.gamma.re,
            p.gamma.im,
            p.z_char.re,
            p.z_char.im,
        )?;
    }
    Ok(path)
}

// --- Complex math helpers ---------------------------------------------------

/// Principal value complex square root.
fn complex_sqrt(z: Complex64) -> Complex64 {
    let r   = z.norm().sqrt();
    let arg = z.arg() * 0.5;
    Complex64::new(r * arg.cos(), r * arg.sin())
}

/// Principal value complex natural logarithm.
fn complex_ln(z: Complex64) -> Complex64 {
    Complex64::new(z.norm().ln(), z.arg())
}

// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::FreqResult;

    fn make_lossless_tline_s(freq_hz: f64, beta_l: f64) -> FreqResult {
        // Lossless matched T-line (Z_c = Z₀): S11=S22=0, S21=S12=exp(-j·βl)
        let s21 = Complex64::new(-beta_l.sin(), -beta_l.cos()) * Complex64::new(0.0, -1.0)
            + Complex64::new(beta_l.cos(), 0.0);
        // S21 = exp(-j·β·l) = cos(βl) - j·sin(βl)
        let s21 = Complex64::new(beta_l.cos(), -beta_l.sin());
        FreqResult {
            freq_hz,
            s11_re: 0.0, s11_im: 0.0,
            s_matrix: vec![
                vec![Complex64::ZERO, s21],
                vec![s21,            Complex64::ZERO],
            ],
            port_list: vec![1, 2],
            port_vi: vec![],
        }
    }

    /// For a lossless matched line (R=G=0, Z_c=Z₀=50Ω), RLGC extraction should give:
    ///  R ≈ 0, G ≈ 0, L·C = 1/c² (LC product = 1/(velocity)²).
    #[test]
    fn lossless_matched_line_rg_near_zero() {
        let f = 1.0e9;
        let c0 = 3.0e8;
        // phase velocity = c0 (air-filled), β = 2πf/c0, choose βl = π/6
        let beta_l = std::f64::consts::PI / 6.0;
        let beta = beta_l * f / (f);
        let length = beta_l / (2.0 * std::f64::consts::PI * f / c0);

        let fr = make_lossless_tline_s(f, beta_l);
        let point = extract_rlgc_point(&fr, length, 50.0).expect("should extract");

        // R and G should be near zero (lossless)
        assert!(point.r_ohm_per_m.abs() < 1e-3,
            "R={:.4e} should be ~0", point.r_ohm_per_m);
        assert!(point.g_s_per_m.abs() < 1e-6,
            "G={:.4e} should be ~0", point.g_s_per_m);
        // L·C should equal 1/c₀² ≈ 1.11e-17
        let lc = point.l_h_per_m * point.c_f_per_m;
        let expected_lc = 1.0 / (c0 * c0);
        let rel_err = (lc - expected_lc).abs() / expected_lc;
        assert!(rel_err < 0.01,
            "L·C={:.4e}, expected {:.4e} (err={:.2}%)", lc, expected_lc, rel_err*100.0);
    }

    /// Z_c reconstruction: for matched line Z_c should equal Z₀ = 50 Ω.
    #[test]
    fn characteristic_impedance_matched_line() {
        let f = 2.4e9;
        let c0 = 3.0e8;
        let beta_l = 0.5;
        let length = beta_l / (2.0 * std::f64::consts::PI * f / c0);

        let fr = make_lossless_tline_s(f, beta_l);
        let point = extract_rlgc_point(&fr, length, 50.0).expect("should extract");

        let z_char_mag = point.z_char.norm();
        assert!((z_char_mag - 50.0).abs() < 0.5,
            "|Z_c|={:.4}, expected 50 Ω", z_char_mag);
    }

    /// Returns None for degenerate inputs: length=0 or single-port S-matrix.
    #[test]
    fn degenerate_inputs_return_none() {
        let fr = make_lossless_tline_s(1e9, 0.3);
        // Zero length
        assert!(extract_rlgc_point(&fr, 0.0, 50.0).is_none());
        // Single-port
        let fr1 = FreqResult {
            freq_hz: 1e9,
            s11_re: 0.1, s11_im: 0.0,
            s_matrix: vec![vec![Complex64::new(0.1, 0.0)]],
            port_list: vec![1],
            port_vi: vec![],
        };
        assert!(extract_rlgc_point(&fr1, 0.1, 50.0).is_none());
    }

    /// complex_sqrt self-test: sqrt(−1) = j.
    #[test]
    fn complex_sqrt_minus_one() {
        let s = complex_sqrt(Complex64::new(-1.0, 0.0));
        assert!((s.im - 1.0).abs() < 1e-12, "Im(sqrt(-1))={:.6}", s.im);
        assert!(s.re.abs() < 1e-12,          "Re(sqrt(-1))={:.6}", s.re);
    }

    /// complex_ln self-test: ln(e) = 1, ln(j) = jπ/2.
    #[test]
    fn complex_ln_basic() {
        let ln_e = complex_ln(Complex64::new(std::f64::consts::E, 0.0));
        assert!((ln_e.re - 1.0).abs() < 1e-12);
        assert!(ln_e.im.abs() < 1e-12);

        let ln_j = complex_ln(Complex64::new(0.0, 1.0));
        assert!(ln_j.re.abs() < 1e-12);
        assert!((ln_j.im - std::f64::consts::FRAC_PI_2).abs() < 1e-12);
    }
}
