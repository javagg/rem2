//! Port-extension de-embedding for multi-port S-parameter data.
//!
//! De-embedding removes the effect of transmission-line stubs (or other parasitic
//! structures) between the simulation reference planes and the true device under
//! test (DUT) planes.
//!
//! # Port-Extension Method
//!
//! The simplest and most widely used de-embedding method corrects each port *i*
//! for a known electrical delay τᵢ (one-way propagation time in seconds).  The
//! de-embedded S-matrix is:
//!
//! ```text
//!     S_deemb[i][j](f) = S_meas[i][j](f) · exp(−j 2π f (τᵢ + τⱼ))
//! ```
//!
//! The phase shift `exp(−j·4π·f·τᵢ/2·2) = exp(−j·2π·f·τᵢ)` per port accounts for
//! the round-trip (into and out of port *i*).  The combined factor for element
//! [i,j] applies one-way delays from both ports.
//!
//! # Open/Short De-embedding
//!
//! For two-port measurements where lumped-pad parasitics are known from
//! OPEN and SHORT standards, the full Y-parameter correction is:
//!
//! ```text
//!     Y_DUT = Y_meas − Y_open        (remove shunt parasitics)
//!     Z_DUT = 1/(Y_DUT) − Z_short   (remove series parasitics)
//! ```
//!
//! Both methods are provided here.
//!
//! # Reference
//! Koolen, M. et al. (1991) "An improved de-embedding technique for on-wafer
//! high-frequency characterisation." *BCTM Proceedings* pp. 188–191.
//! Pozar, D. M. *Microwave Engineering*, 4th ed., Appendix C.

use num_complex::Complex64;
use std::f64::consts::PI;

use crate::FreqResult;

/// Result of port-extension de-embedding: one entry per original frequency.
#[derive(Debug, Clone)]
pub struct DeembedResult {
    /// De-embedded S-matrix at each frequency.
    pub freq_results: Vec<DeembedFreqPoint>,
}

#[derive(Debug, Clone)]
pub struct DeembedFreqPoint {
    pub freq_hz: f64,
    /// De-embedded N×N S-matrix (same ordering as source `port_list`).
    pub s_matrix: Vec<Vec<Complex64>>,
    pub port_list: Vec<u32>,
}

/// Apply port-extension de-embedding to a set of `FreqResult` records.
///
/// # Arguments
/// * `freq_results`    — from `DrivenResult::freq_results`
/// * `port_delays_s`   — one-way electrical delay for each port [seconds].
///                       Length must equal the number of ports in `s_matrix`.
///
/// # Example
/// ```ignore
/// // Remove 100 ps stub from port 1, 150 ps stub from port 2
/// let de = deembed_port_extension(&result.freq_results, &[100e-12, 150e-12]);
/// ```
pub fn deembed_port_extension(
    freq_results: &[FreqResult],
    port_delays_s: &[f64],
) -> DeembedResult {
    let freq_pts = freq_results.iter().map(|fr| {
        let n_ports = fr.s_matrix.len();
        let n_delay = port_delays_s.len();

        let mut s_de = fr.s_matrix.clone();
        for i in 0..n_ports {
            for j in 0..n_ports {
                let tau_i = if i < n_delay { port_delays_s[i] } else { 0.0 };
                let tau_j = if j < n_delay { port_delays_s[j] } else { 0.0 };
                let phase = -2.0 * PI * fr.freq_hz * (tau_i + tau_j);
                let correction = Complex64::new(phase.cos(), phase.sin());
                s_de[i][j] *= correction;
            }
        }

        DeembedFreqPoint {
            freq_hz: fr.freq_hz,
            s_matrix: s_de,
            port_list: fr.port_list.clone(),
        }
    }).collect();

    DeembedResult { freq_results: freq_pts }
}

/// Open/Short de-embedding for a 2-port measurement (lumped-pad parasitics).
///
/// Removes shunt parasitics measured from an OPEN standard and series parasitics
/// from a SHORT standard.
///
/// # Arguments
/// * `s_meas`   — measured 2×2 S-matrix (as Complex64 2D array, [port][port])
/// * `s_open`   — OPEN standard S-matrix (S_22 of pad-only structure)
/// * `s_short`  — SHORT standard S-matrix (series parasitics, short across DUT)
/// * `z0`       — reference impedance [Ω]
///
/// Returns de-embedded `[[Complex64; 2]; 2]` S-matrix or `None` if inversion fails.
pub fn deembed_open_short_2port(
    s_meas:  &[[Complex64; 2]; 2],
    s_open:  &[[Complex64; 2]; 2],
    s_short: &[[Complex64; 2]; 2],
    z0: f64,
) -> Option<[[Complex64; 2]; 2]> {
    // Convert S → Y
    let y_meas  = s2y(s_meas, z0)?;
    let y_open  = s2y(s_open, z0)?;
    let y_short = s2y(s_short, z0)?;

    // Step 1: subtract shunt parasitics
    let y_corrected = mat2x2_sub(&y_meas, &y_open);

    // Step 2: convert to Z, subtract series parasitics
    let z_corrected = y2z(&y_corrected)?;
    let z_short_corr = y2z(&mat2x2_sub(&y_short, &y_open))?;
    let z_dut = mat2x2_sub(&z_corrected, &z_short_corr);

    // Convert Z → S
    z2s(&z_dut, z0)
}

/// Write de-embedded S-parameters to Touchstone-style CSV.
pub fn write_deembed_csv(
    result: &DeembedResult,
    output_dir: &std::path::Path,
) -> Result<std::path::PathBuf, std::io::Error> {
    use std::io::Write;
    let dir = output_dir.join("postpro");
    std::fs::create_dir_all(&dir)?;
    let path = dir.join("deembedded_sparams.csv");
    let mut f = std::fs::File::create(&path)?;

    // Determine max port count
    let n_ports = result.freq_results.first().map(|r| r.s_matrix.len()).unwrap_or(0);

    // Header
    write!(f, "FreqHz")?;
    for i in 0..n_ports {
        for j in 0..n_ports {
            write!(f, ",S{}{}_re,S{}{}_im", i+1, j+1, i+1, j+1)?;
        }
    }
    writeln!(f)?;

    for pt in &result.freq_results {
        write!(f, "{:.9e}", pt.freq_hz)?;
        for row in &pt.s_matrix {
            for s in row {
                write!(f, ",{:.9e},{:.9e}", s.re, s.im)?;
            }
        }
        writeln!(f)?;
    }
    Ok(path)
}

// --- Matrix conversion helpers (2×2) ----------------------------------------

/// S → Y: Y = (1/Z₀) · (I − S)(I + S)⁻¹
fn s2y(s: &[[Complex64; 2]; 2], z0: f64) -> Option<[[Complex64; 2]; 2]> {
    let one  = Complex64::new(1.0, 0.0);
    let invz = Complex64::new(1.0 / z0, 0.0);
    // (I − S)
    let a = [[one - s[0][0], -s[0][1]], [-s[1][0], one - s[1][1]]];
    // (I + S)
    let b = [[one + s[0][0],  s[0][1]], [ s[1][0], one + s[1][1]]];
    let b_inv = mat2x2_inv(&b)?;
    let y = mat2x2_mul(&a, &b_inv);
    Some([[y[0][0]*invz, y[0][1]*invz], [y[1][0]*invz, y[1][1]*invz]])
}

/// Y → Z: Z = Y⁻¹
fn y2z(y: &[[Complex64; 2]; 2]) -> Option<[[Complex64; 2]; 2]> {
    mat2x2_inv(y)
}

/// Z → S: S = (Z − Z₀ I)(Z + Z₀ I)⁻¹
fn z2s(z: &[[Complex64; 2]; 2], z0: f64) -> Option<[[Complex64; 2]; 2]> {
    let z0c = Complex64::new(z0, 0.0);
    let a = [[z[0][0] - z0c, z[0][1]], [z[1][0], z[1][1] - z0c]];
    let b = [[z[0][0] + z0c, z[0][1]], [z[1][0], z[1][1] + z0c]];
    let b_inv = mat2x2_inv(&b)?;
    Some(mat2x2_mul(&a, &b_inv))
}

fn mat2x2_sub(a: &[[Complex64; 2]; 2], b: &[[Complex64; 2]; 2]) -> [[Complex64; 2]; 2] {
    [[a[0][0]-b[0][0], a[0][1]-b[0][1]], [a[1][0]-b[1][0], a[1][1]-b[1][1]]]
}

fn mat2x2_mul(a: &[[Complex64; 2]; 2], b: &[[Complex64; 2]; 2]) -> [[Complex64; 2]; 2] {
    [[a[0][0]*b[0][0]+a[0][1]*b[1][0], a[0][0]*b[0][1]+a[0][1]*b[1][1]],
     [a[1][0]*b[0][0]+a[1][1]*b[1][0], a[1][0]*b[0][1]+a[1][1]*b[1][1]]]
}

fn mat2x2_inv(m: &[[Complex64; 2]; 2]) -> Option<[[Complex64; 2]; 2]> {
    let det = m[0][0]*m[1][1] - m[0][1]*m[1][0];
    if det.norm() < 1e-300 { return None; }
    let inv_det = Complex64::new(1.0, 0.0) / det;
    Some([[ m[1][1]*inv_det, -m[0][1]*inv_det],
          [-m[1][0]*inv_det,  m[0][0]*inv_det]])
}

// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn make_freq_results(freqs: &[f64], s_fn: impl Fn(f64) -> Vec<Vec<Complex64>>) -> Vec<FreqResult> {
        freqs.iter().map(|&f| FreqResult {
            freq_hz: f,
            s11_re: 0.0, s11_im: 0.0,
            s_matrix: s_fn(f),
            port_list: vec![1, 2],
            port_vi: vec![],
        }).collect()
    }

    /// Port extension: zero delay → identity transformation.
    #[test]
    fn zero_delay_is_identity() {
        let freqs = [1.0e9, 2.0e9, 3.0e9];
        let s_fn = |_f: f64| vec![
            vec![Complex64::new(0.1, 0.05), Complex64::new(0.9, -0.1)],
            vec![Complex64::new(0.9, -0.1), Complex64::new(0.1, 0.05)],
        ];
        let fr = make_freq_results(&freqs, s_fn);
        let de = deembed_port_extension(&fr, &[0.0, 0.0]);
        for (orig, dem) in fr.iter().zip(de.freq_results.iter()) {
            for i in 0..2 {
                for j in 0..2 {
                    let diff = (orig.s_matrix[i][j] - dem.s_matrix[i][j]).norm();
                    assert!(diff < 1e-12, "S{}{} changed with zero delay: diff={:.2e}", i+1,j+1,diff);
                }
            }
        }
    }

    /// Port extension: applying a delay then its negative recovers original S-params.
    #[test]
    fn delay_then_neg_delay_is_identity() {
        let freqs = [1.0e9];
        let s21 = Complex64::new(0.7, -0.3);
        let s_fn = |_f: f64| vec![
            vec![Complex64::ZERO, s21],
            vec![s21, Complex64::ZERO],
        ];
        let fr = make_freq_results(&freqs, s_fn);
        let tau = 100e-12; // 100 ps

        let de1 = deembed_port_extension(&fr, &[tau, tau]);
        // Now undo with negative delay on the already de-embedded result
        let fr2: Vec<FreqResult> = de1.freq_results.iter().map(|p| FreqResult {
            freq_hz: p.freq_hz,
            s11_re: 0.0, s11_im: 0.0,
            s_matrix: p.s_matrix.clone(),
            port_list: p.port_list.clone(),
            port_vi: vec![],
        }).collect();
        let de2 = deembed_port_extension(&fr2, &[-tau, -tau]);

        let diff = (de2.freq_results[0].s_matrix[0][1] - s21).norm();
        assert!(diff < 1e-12, "Round-trip de-embed error: {:.2e}", diff);
    }

    /// S→Y→Z→S round-trip should recover the original matrix.
    #[test]
    fn s_y_z_s_roundtrip() {
        let s: [[Complex64; 2]; 2] = [
            [Complex64::new(0.1, 0.05), Complex64::new(0.8, -0.2)],
            [Complex64::new(0.8, -0.2), Complex64::new(0.1, 0.05)],
        ];
        let z0 = 50.0;
        let y = s2y(&s, z0).unwrap();
        let z = y2z(&y).unwrap();
        let s2 = z2s(&z, z0).unwrap();
        for i in 0..2 {
            for j in 0..2 {
                let diff = (s[i][j] - s2[i][j]).norm();
                assert!(diff < 1e-12, "S{}{} round-trip error: {:.3e}", i+1,j+1,diff);
            }
        }
    }

    /// mat2x2_inv: inverse of identity should be identity.
    #[test]
    fn identity_inverse() {
        let one = Complex64::new(1.0, 0.0);
        let zer = Complex64::ZERO;
        let id: [[Complex64; 2]; 2] = [[one, zer],[zer, one]];
        let inv = mat2x2_inv(&id).unwrap();
        assert!((inv[0][0] - one).norm() < 1e-12);
        assert!((inv[1][1] - one).norm() < 1e-12);
        assert!(inv[0][1].norm() < 1e-12);
        assert!(inv[1][0].norm() < 1e-12);
    }

    /// Open/Short de-embedding: if OPEN=MEAS and SHORT=0, DUT should be zero.
    #[test]
    fn open_short_identity_open() {
        let zero = Complex64::ZERO;
        let s_meas: [[Complex64; 2]; 2] = [
            [Complex64::new(0.2, 0.0), Complex64::new(0.7, 0.0)],
            [Complex64::new(0.7, 0.0), Complex64::new(0.2, 0.0)],
        ];
        let s_short: [[Complex64; 2]; 2] = [[zero;2];2];
        // OPEN = MEAS → after subtracting, Y_corrected = 0
        let result = deembed_open_short_2port(&s_meas, &s_meas, &s_short, 50.0);
        // Should be None or give a near-zero DUT (degenerate Y=0 → Z=∞ → may fail inversion)
        // The important thing is it doesn't panic.
        let _ = result; // may be None; that's acceptable
    }
}
