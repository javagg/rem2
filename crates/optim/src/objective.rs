//! Objective function evaluation for parametric optimization.
//!
//! Given a set of S-matrices (one per frequency), evaluates a scalar cost
//! from an [`OptimObjective`] specification.

use num_complex::Complex64;
use rem_config::OptimObjective;
use rem_mom::sparams::SMatrix;

/// Evaluate a list of objectives against a frequency sweep result.
/// Returns the sum of all objective values (all are to be minimized).
pub fn evaluate_objectives(
    matrices: &[SMatrix],
    objectives: &[OptimObjective],
) -> f64 {
    objectives.iter().map(|obj| evaluate_one(matrices, obj)).sum()
}

fn evaluate_one(matrices: &[SMatrix], obj: &OptimObjective) -> f64 {
    match obj {
        OptimObjective::MinS11dB { port, freq_hz } => {
            let pi = port.saturating_sub(1); // 1-indexed → 0-indexed
            s_mag_db_at_freq(matrices, pi, pi, *freq_hz)
        }
        OptimObjective::MinSijdB { port_i, port_j, freq_hz } => {
            let pi = port_i.saturating_sub(1);
            let pj = port_j.saturating_sub(1);
            s_mag_db_at_freq(matrices, pi, pj, *freq_hz)
        }
        OptimObjective::MaxBandwidthS11dB { port, thresh_db, freq_min_hz, freq_max_hz } => {
            let pi = port.saturating_sub(1);
            // Negative bandwidth (minimizer drives this toward most negative)
            -bandwidth_s11_hz(matrices, pi, *thresh_db, *freq_min_hz, *freq_max_hz)
        }
        OptimObjective::TargetS11dB { port, freq_hz, target_db } => {
            let pi = port.saturating_sub(1);
            let val = s_mag_db_at_freq(matrices, pi, pi, *freq_hz);
            (val - target_db).powi(2)
        }
    }
}

/// |S_{row,col}| in dB at the frequency nearest `freq_hz`.
/// Returns 0.0 if no matrices available or port index out of range.
fn s_mag_db_at_freq(matrices: &[SMatrix], row: usize, col: usize, freq_hz: f64) -> f64 {
    let sm = nearest_matrix(matrices, freq_hz);
    let Some(sm) = sm else { return 0.0 };
    if row >= sm.n_ports || col >= sm.n_ports { return 0.0 }
    let s: Complex64 = sm.get(row, col);
    let mag = s.norm();
    if mag < 1e-30 { return -300.0; }
    20.0 * mag.log10()
}

/// Bandwidth [Hz] in [freq_min, freq_max] where |S11|_dB < thresh_db.
fn bandwidth_s11_hz(
    matrices: &[SMatrix],
    port: usize,
    thresh_db: f64,
    freq_min: f64,
    freq_max: f64,
) -> f64 {
    // Collect frequency points in [freq_min, freq_max] where S11 < thresh
    let relevant: Vec<&SMatrix> = matrices.iter()
        .filter(|m| m.freq_hz >= freq_min && m.freq_hz <= freq_max)
        .collect();

    if relevant.len() < 2 { return 0.0; }

    // Trapezoidal integration of indicator function
    let mut bw = 0.0;
    for w in relevant.windows(2) {
        let f1 = w[0].freq_hz;
        let f2 = w[1].freq_hz;
        let s1 = s_mag_db_at_freq(&[w[0].clone()], port, port, f1);
        let s2 = s_mag_db_at_freq(&[w[1].clone()], port, port, f2);
        let in1 = if s1 < thresh_db { 1.0 } else { 0.0 };
        let in2 = if s2 < thresh_db { 1.0 } else { 0.0 };
        bw += 0.5 * (in1 + in2) * (f2 - f1);
    }
    bw
}

fn nearest_matrix(matrices: &[SMatrix], freq_hz: f64) -> Option<&SMatrix> {
    matrices.iter().min_by(|a, b| {
        let da = (a.freq_hz - freq_hz).abs();
        let db = (b.freq_hz - freq_hz).abs();
        da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use rem_config::OptimObjective;

    fn make_smatrix(freq: f64, s11: Complex64) -> SMatrix {
        SMatrix { n_ports: 1, freq_hz: freq, data: vec![s11] }
    }

    #[test]
    fn min_s11_db_returns_correct_db() {
        let s11 = Complex64::new(0.1, 0.0); // −20 dB
        let mats = vec![make_smatrix(2.4e9, s11)];
        let obj = OptimObjective::MinS11dB { port: 1, freq_hz: 2.4e9 };
        let v = evaluate_objectives(&mats, &[obj]);
        assert!((v - (-20.0_f64)).abs() < 0.1, "v={:.2}", v);
    }

    #[test]
    fn target_s11_zero_when_on_target() {
        let s11 = Complex64::new(0.1, 0.0); // −20 dB
        let mats = vec![make_smatrix(2.4e9, s11)];
        let obj = OptimObjective::TargetS11dB { port: 1, freq_hz: 2.4e9, target_db: -20.0 };
        let v = evaluate_objectives(&mats, &[obj]);
        assert!(v < 0.01, "v={:.6}", v);
    }

    #[test]
    fn bandwidth_nonzero_when_below_thresh() {
        let mats: Vec<SMatrix> = (0..=10)
            .map(|i| {
                let f = 2.0e9 + i as f64 * 0.1e9;
                make_smatrix(f, Complex64::new(0.1, 0.0)) // all −20 dB < −10 dB
            })
            .collect();
        let obj = OptimObjective::MaxBandwidthS11dB {
            port: 1, thresh_db: -10.0, freq_min_hz: 2.0e9, freq_max_hz: 3.0e9,
        };
        let v = evaluate_objectives(&mats, &[obj]);
        // objective is -bandwidth, bandwidth ≈ 1 GHz
        assert!(v < -0.8e9, "v={:.3e}", v);
    }
}
