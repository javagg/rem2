//! S-parameter analysis utilities: group delay, stability, gain metrics.
//!
//! # Group Delay
//! The group delay τ_g is the negative derivative of the S-parameter phase
//! with respect to angular frequency ω:
//!
//! ```text
//!     τ_g(f) = −dφ/dω = −dφ/(2π df)
//! ```
//!
//! Computed by central finite differences on the unwrapped phase.
//!
//! # Rollett Stability Factor K
//!
//! The Rollett stability factor (Rollett 1962) determines unconditional
//! stability of a 2-port microwave network:
//!
//! ```text
//!     K = (1 − |S11|² − |S22|² + |Δ|²) / (2 |S12 S21|)
//!
//!     Δ = S11 S22 − S12 S21   (determinant of S-matrix)
//! ```
//!
//! Unconditional stability requires **K > 1** AND **|Δ| < 1** simultaneously.
//!
//! # Maximum Stable Gain / Maximum Available Gain
//!
//! ```text
//!     MSG = |S21| / |S12|          (potentially unstable: K < 1)
//!     MAG = MSG · (K − √(K²−1))   (unconditionally stable: K ≥ 1)
//! ```
//!
//! Both MSG and MAG are in linear (power ratio) units; convert to dB with 10·log₁₀.
//!
//! # Transducer Power Gain
//!
//! ```text
//!     G_T = |S21|²   [linear, reference 50 Ω system]
//! ```
//!
//! # References
//! Rollett, J. M. (1962) "Stability and Power-Gain Invariants of Linear Twoports."
//! *IRE Trans. Circuit Theory* CT-9(1):29–32.
//!
//! Pozar, D. M. *Microwave Engineering*, 4th ed., Ch. 11.

use crate::{DrivenResult, FreqResult};
use num_complex::Complex64;
use std::f64::consts::PI;

/// Per-frequency analysis result for a 2-port network.
#[derive(Debug, Clone)]
pub struct TwoPortAnalysis {
    pub freq_hz: f64,
    /// Group delay of S21 [s] — NaN at endpoints (finite diff needs neighbours)
    pub group_delay_s21_s: f64,
    /// Group delay of S11 [s]
    pub group_delay_s11_s: f64,
    /// Rollett K-factor (unconditionally stable when K > 1 AND delta_mag < 1)
    pub rollett_k: f64,
    /// |Δ| = |S11·S22 − S12·S21|
    pub delta_mag: f64,
    /// Unconditionally stable (K > 1 AND |Δ| < 1)
    pub unconditionally_stable: bool,
    /// Maximum Stable Gain [linear] — valid when K < 1
    pub msg_linear: f64,
    /// Maximum Available Gain [linear] — valid when K ≥ 1
    pub mag_linear: f64,
    /// Transducer power gain |S21|² [linear]
    pub gain_t_linear: f64,
    /// |S11|² (input return loss power)
    pub s11_mag_sq: f64,
    /// |S22|² (output return loss power)
    pub s22_mag_sq: f64,
}

impl TwoPortAnalysis {
    /// MSG in dB (10·log₁₀(MSG)), or f64::NEG_INFINITY if S12 = 0.
    pub fn msg_db(&self) -> f64 {
        if self.msg_linear <= 0.0 { f64::NEG_INFINITY } else { 10.0 * self.msg_linear.log10() }
    }
    /// MAG in dB.
    pub fn mag_db(&self) -> f64 {
        if self.mag_linear <= 0.0 { f64::NEG_INFINITY } else { 10.0 * self.mag_linear.log10() }
    }
    /// Transducer gain in dB.
    pub fn gain_t_db(&self) -> f64 {
        if self.gain_t_linear <= 0.0 { f64::NEG_INFINITY } else { 10.0 * self.gain_t_linear.log10() }
    }
    /// Input return loss in dB (positive = good matching).
    pub fn s11_return_loss_db(&self) -> f64 {
        if self.s11_mag_sq <= 0.0 { f64::INFINITY } else { -10.0 * self.s11_mag_sq.log10() }
    }
}

/// Analyse a 2-port driven sweep.
///
/// Computes group delay (central FD), Rollett K, MSG/MAG, and gain at each
/// frequency point.  The `z0` argument is the reference impedance [Ω] (typically 50.0);
/// it is used only for documentation — S-matrix is already normalised.
///
/// Returns an empty `Vec` if the S-matrix has fewer than 2 ports.
pub fn analyse_two_port(result: &DrivenResult, _z0: f64) -> Vec<TwoPortAnalysis> {
    let freq_results = &result.freq_results;
    let n = freq_results.len();
    if n == 0 { return vec![]; }

    // Check 2-port data is present
    if freq_results[0].s_matrix.len() < 2 { return vec![]; }

    // Compute stability and gain metrics at each point
    let mut out: Vec<TwoPortAnalysis> = freq_results.iter().map(|fr| {
        let s11 = fr.s_matrix[0][0];
        let s12 = fr.s_matrix[0][1];
        let s21 = fr.s_matrix[1][0];
        let s22 = fr.s_matrix[1][1];

        let delta = s11 * s22 - s12 * s21;
        let delta_mag = delta.norm();

        let s11sq = s11.norm_sqr();
        let s22sq = s22.norm_sqr();
        let s12s21 = (s12 * s21).norm();

        let k_denom = 2.0 * s12s21;
        let rollett_k = if k_denom > 1e-30 {
            (1.0 - s11sq - s22sq + delta_mag * delta_mag) / k_denom
        } else {
            f64::INFINITY
        };

        let s21sq = s21.norm_sqr();
        let s12mag = s12.norm();
        let s21mag = s21.norm();

        let msg_linear = if s12mag > 1e-30 { s21mag / s12mag } else { f64::INFINITY };
        let mag_linear = if rollett_k >= 1.0 {
            let k = rollett_k;
            msg_linear * (k - (k * k - 1.0).sqrt())
        } else {
            f64::NAN
        };

        TwoPortAnalysis {
            freq_hz: fr.freq_hz,
            group_delay_s21_s: f64::NAN, // filled in next pass
            group_delay_s11_s: f64::NAN,
            rollett_k,
            delta_mag,
            unconditionally_stable: rollett_k > 1.0 && delta_mag < 1.0,
            msg_linear,
            mag_linear,
            gain_t_linear: s21sq,
            s11_mag_sq: s11sq,
            s22_mag_sq: s22sq,
        }
    }).collect();

    // Group delay: −dφ/dω via central finite differences on unwrapped phase
    let phases_s21: Vec<f64> = freq_results.iter()
        .map(|fr| fr.s_matrix[1][0].arg())
        .collect();
    let phases_s11: Vec<f64> = freq_results.iter()
        .map(|fr| fr.s_matrix[0][0].arg())
        .collect();

    let unwrapped_s21 = unwrap_phase(&phases_s21);
    let unwrapped_s11 = unwrap_phase(&phases_s11);

    for i in 0..n {
        let gd_s21 = if i == 0 || i == n - 1 {
            f64::NAN
        } else {
            let df = freq_results[i + 1].freq_hz - freq_results[i - 1].freq_hz;
            if df.abs() < 1e-6 { f64::NAN }
            else {
                -(unwrapped_s21[i + 1] - unwrapped_s21[i - 1]) / (2.0 * PI * df)
            }
        };
        let gd_s11 = if i == 0 || i == n - 1 {
            f64::NAN
        } else {
            let df = freq_results[i + 1].freq_hz - freq_results[i - 1].freq_hz;
            if df.abs() < 1e-6 { f64::NAN }
            else {
                -(unwrapped_s11[i + 1] - unwrapped_s11[i - 1]) / (2.0 * PI * df)
            }
        };
        out[i].group_delay_s21_s = gd_s21;
        out[i].group_delay_s11_s = gd_s11;
    }

    out
}

/// Phase unwrapping: adjust each successive sample to minimise discontinuities.
///
/// Handles the ±π wrapping by adding multiples of 2π as needed.
fn unwrap_phase(phases: &[f64]) -> Vec<f64> {
    let mut out = phases.to_vec();
    for i in 1..out.len() {
        let mut diff = out[i] - out[i - 1];
        // Wrap diff into (−π, π]
        while diff > PI  { diff -= 2.0 * PI; }
        while diff < -PI { diff += 2.0 * PI; }
        out[i] = out[i - 1] + diff;
    }
    out
}

/// Write 2-port analysis results to CSV.
pub fn write_two_port_analysis_csv(
    data: &[TwoPortAnalysis],
    output_dir: &std::path::Path,
) -> Result<std::path::PathBuf, std::io::Error> {
    use std::io::Write;
    let dir = output_dir.join("postpro");
    std::fs::create_dir_all(&dir)?;
    let path = dir.join("two_port_analysis.csv");
    let mut f = std::fs::File::create(&path)?;
    writeln!(f, "FreqHz,GroupDelay_S21_ps,GroupDelay_S11_ps,\
                 Rollett_K,Delta_mag,UnconditionallyStable,\
                 MSG_dB,MAG_dB,GainT_dB,S11_ReturnLoss_dB,S22_mag_sq")?;
    for p in data {
        let gd_s21_ps = if p.group_delay_s21_s.is_nan() { f64::NAN }
                        else { p.group_delay_s21_s * 1e12 };
        let gd_s11_ps = if p.group_delay_s11_s.is_nan() { f64::NAN }
                        else { p.group_delay_s11_s * 1e12 };
        writeln!(f, "{:.9e},{},{},{:.6e},{:.6e},{},{:.4e},{:.4e},{:.4e},{:.4e},{:.6e}",
            p.freq_hz,
            if gd_s21_ps.is_nan() { "NaN".to_string() } else { format!("{gd_s21_ps:.4e}") },
            if gd_s11_ps.is_nan() { "NaN".to_string() } else { format!("{gd_s11_ps:.4e}") },
            p.rollett_k,
            p.delta_mag,
            if p.unconditionally_stable { 1 } else { 0 },
            p.msg_db(),
            if p.mag_linear.is_nan() { f64::NAN } else { p.mag_db() },
            p.gain_t_db(),
            p.s11_return_loss_db(),
            p.s22_mag_sq,
        )?;
    }
    Ok(path)
}

// ---------------------------------------------------------------------------

/// Compute the unilateral transducer power gain |S21|² at each frequency [linear].
///
/// Also returns |S11|² and |S22|² as a convenience triple.
pub fn gain_sweep(freq_results: &[FreqResult]) -> Vec<(f64, f64, f64, f64)> {
    freq_results.iter().filter_map(|fr| {
        if fr.s_matrix.len() < 2 { return None; }
        let s11sq = fr.s_matrix[0][0].norm_sqr();
        let s21sq = fr.s_matrix[1][0].norm_sqr();
        let s22sq = fr.s_matrix[1][1].norm_sqr();
        Some((fr.freq_hz, s21sq, s11sq, s22sq))
    }).collect()
}

// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{DrivenResult, FreqResult};

    fn make_driven_result(freq_results: Vec<FreqResult>) -> DrivenResult {
        DrivenResult {
            freq_results,
            peak_phi: vec![],
            peak_freq_hz: 0.0,
            far_field_pattern: vec![],
            circuit_model: None,
        }
    }

    fn s2x2(s11: Complex64, s12: Complex64, s21: Complex64, s22: Complex64, f: f64) -> FreqResult {
        FreqResult {
            freq_hz: f,
            s11_re: s11.re, s11_im: s11.im,
            s_matrix: vec![vec![s11, s12], vec![s21, s22]],
            port_list: vec![1, 2],
            port_vi: vec![],
        }
    }

    /// Single-frequency: matched lossless 2-port with S21=S12=exp(−jπ/4), S11=S22=0.
    /// Δ = S11·S22 − S12·S21 = −S21² → |Δ| = |S21|² = 1.
    /// K = (1 − 0 − 0 + 1)/(2·|S12·S21|) = 2/2 = 1.0 (marginally stable).
    #[test]
    fn stability_lossless_matched() {
        let phase = std::f64::consts::FRAC_PI_4;
        let s21 = Complex64::new(phase.cos(), -phase.sin());
        let fr = s2x2(Complex64::ZERO, s21, s21, Complex64::ZERO, 1.0e9);
        let result = make_driven_result(vec![fr]);
        let analysis = analyse_two_port(&result, 50.0);
        let a = &analysis[0];
        // K = 1.0 exactly (marginally stable)
        assert!((a.rollett_k - 1.0).abs() < 1e-10, "K={:.6}", a.rollett_k);
        // |Δ| = |−S21²| = 1.0
        assert!((a.delta_mag - 1.0).abs() < 1e-10, "|Δ|={:.6}", a.delta_mag);
        // K > 1 is false, |Δ| < 1 is false → unconditionally_stable = false
        assert!(!a.unconditionally_stable);
        // MSG = |S21|/|S12| = 1 → 0 dB
        assert!((a.msg_linear - 1.0).abs() < 1e-10);
        assert!(a.msg_db().abs() < 1e-8, "MSG_dB={:.6}", a.msg_db());
    }

    /// Unconditionally stable network: verify K > 1 and |Δ| < 1.
    #[test]
    fn stability_unconditionally_stable() {
        // Small S11, S22, very small S21, large loss
        let s11 = Complex64::new(0.1, 0.0);
        let s22 = Complex64::new(0.1, 0.0);
        let s21 = Complex64::new(0.01, 0.0);
        let s12 = Complex64::new(0.001, 0.0); // very asymmetric
        let fr = s2x2(s11, s12, s21, s22, 1.0e9);
        let result = make_driven_result(vec![fr]);
        let analysis = analyse_two_port(&result, 50.0);
        let a = &analysis[0];
        assert!(a.rollett_k > 1.0, "K={:.4}", a.rollett_k);
        assert!(a.delta_mag < 1.0, "|Δ|={:.4}", a.delta_mag);
        assert!(a.unconditionally_stable);
        // MAG should be finite and > 0
        assert!(a.mag_linear > 0.0 && a.mag_linear.is_finite());
    }

    /// Group delay: for a linear-phase S21 = exp(−j·2πf·τ), group delay = τ.
    #[test]
    fn group_delay_constant_phase_slope() {
        let tau = 100.0e-12_f64; // 100 ps
        let freqs = [1.0e9, 2.0e9, 3.0e9, 4.0e9, 5.0e9];
        let frs: Vec<FreqResult> = freqs.iter().map(|&f| {
            let phi = -2.0 * PI * f * tau;
            let s21 = Complex64::new(phi.cos(), phi.sin());
            s2x2(Complex64::ZERO, s21, s21, Complex64::ZERO, f)
        }).collect();
        let result = make_driven_result(frs);
        let analysis = analyse_two_port(&result, 50.0);
        // Interior points should give ≈ τ
        for a in &analysis[1..analysis.len()-1] {
            if !a.group_delay_s21_s.is_nan() {
                let err = (a.group_delay_s21_s - tau).abs();
                assert!(err < 1e-12, "τ_g={:.4e} vs expected {tau:.4e}", a.group_delay_s21_s);
            }
        }
    }

    /// Phase unwrap: a step from π−ε to −π+ε should be unwrapped to a small positive step.
    #[test]
    fn unwrap_phase_handles_pi_crossing() {
        let eps = 0.01;
        let phases = vec![PI - eps, -(PI - eps)];
        let unwrapped = unwrap_phase(&phases);
        // Second point should be close to π + ε, not flipped
        let expected = PI + eps;
        assert!((unwrapped[1] - expected).abs() < 0.1,
            "unwrapped[1]={:.6}, expected ≈{:.6}", unwrapped[1], expected);
    }

    /// Returns empty Vec for single-port data.
    #[test]
    fn analyse_single_port_returns_empty() {
        let fr = FreqResult {
            freq_hz: 1e9,
            s11_re: 0.1, s11_im: 0.0,
            s_matrix: vec![vec![Complex64::new(0.1, 0.0)]],
            port_list: vec![1],
            port_vi: vec![],
        };
        let result = make_driven_result(vec![fr]);
        let analysis = analyse_two_port(&result, 50.0);
        assert!(analysis.is_empty());
    }

    /// Gain sweep extracts correct |S21|² values.
    #[test]
    fn gain_sweep_values() {
        let s21 = Complex64::new(0.6, 0.8); // |S21|² = 0.36+0.64 = 1.0
        let fr = s2x2(Complex64::ZERO, s21, s21, Complex64::ZERO, 1e9);
        let gains = gain_sweep(&[fr]);
        assert_eq!(gains.len(), 1);
        let (_, g_t, _, _) = gains[0];
        assert!((g_t - 1.0).abs() < 1e-12, "|S21|²={g_t:.6}");
    }
}
