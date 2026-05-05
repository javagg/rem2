//! Monte Carlo yield analysis for parametric design objectives.
//!
//! Samples each design parameter from an independent Gaussian distribution
//! centred on its nominal value with a user-specified relative standard
//! deviation `σᵢ = McSigmaRel · |pᵢ_nom|`.  For each trial the objective
//! function is evaluated and statistics are accumulated.
//!
//! # Output files
//! - `postpro/monte_carlo_samples.csv`  — one row per trial with sampled
//!   values and the resulting objective
//! - `postpro/monte_carlo_stats.csv`    — summary: mean, std-dev, min, max,
//!   yield at ±1σ / ±3σ windows around the nominal objective
//!
//! # JSON config example
//! ```json
//! "Parametric": {
//!   "Mode": "MonteCarlo",
//!   "Parameters": [
//!     {"Name": "eps_r", "Target": {"Type": "SubstratePermittivity", "Layer": 0},
//!      "Initial": 4.0, "Bounds": [3.0, 5.0]},
//!     {"Name": "h_mm",  "Target": {"Type": "SubstrateThickness",    "Layer": 0},
//!      "Initial": 0.254e-3, "Bounds": [0.1e-3, 0.5e-3]}
//!   ],
//!   "Objectives": [{"Type": "MinS11dB", "Port": 1, "FreqHz": 2.4e9}],
//!   "McSamples":  500,
//!   "McSigmaRel": 0.05,
//!   "McSeed":     42
//! }
//! ```

use rem_config::{PalaceConfig, ParametricConfig};
use rem_core::RemResult;
use rem_parallel::NoComm;

use crate::objective::evaluate_objectives;
use crate::param_apply::apply_params;

// ---------------------------------------------------------------------------
// Minimal portable PRNG + Box-Muller for Gaussian samples.
// Uses xorshift64 — no external crate dependency needed.
// ---------------------------------------------------------------------------

struct Xorshift64 {
    state: u64,
}

impl Xorshift64 {
    fn new(seed: u64) -> Self {
        // Ensure non-zero state (xorshift requires it)
        Self { state: if seed == 0 { 6364136223846793005 } else { seed } }
    }

    /// Returns next pseudo-random u64.
    fn next_u64(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.state = x;
        x
    }

    /// Returns U(0,1) in (0, 1].
    fn next_f64(&mut self) -> f64 {
        // Map to (0, 1] to avoid log(0) in Box-Muller
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64) + f64::EPSILON
    }

    /// Returns a standard-normal sample via Box-Muller transform.
    fn next_normal(&mut self) -> f64 {
        let u1 = self.next_f64();
        let u2 = self.next_f64();
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }
}

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

/// Run Monte Carlo yield analysis and write statistics CSVs.
pub fn run_monte_carlo(config: &PalaceConfig, par_cfg: &ParametricConfig) -> RemResult<()> {
    use std::io::Write;

    if par_cfg.objectives.is_empty() {
        return Err(rem_core::RemError::Config(
            "MonteCarlo: at least one Objectives entry required".to_string(),
        ));
    }
    if par_cfg.parameters.is_empty() {
        return Err(rem_core::RemError::Config(
            "MonteCarlo: at least one Parameters entry required".to_string(),
        ));
    }

    let n_samples   = par_cfg.mc_samples.max(1);
    let sigma_rel   = par_cfg.mc_sigma_rel.abs().max(1e-9);
    let seed        = par_cfg.mc_seed.unwrap_or(0xDEAD_BEEF_1234_5678);

    // Nominal parameter values
    let p_nom: Vec<f64> = par_cfg.parameters.iter().map(|p| {
        p.initial.unwrap_or_else(|| {
            if let Some([lo, hi]) = p.bounds { (lo + hi) * 0.5 } else { 1.0 }
        })
    }).collect();
    let params: Vec<&rem_config::SweepParam> = par_cfg.parameters.iter().collect();
    let bounds: Vec<Option<[f64; 2]>> = par_cfg.parameters.iter().map(|p| p.bounds).collect();
    let n_params = p_nom.len();

    // σᵢ for each parameter
    let sigmas: Vec<f64> = p_nom.iter().map(|&v| sigma_rel * v.abs().max(1e-30)).collect();

    // Pre-load mesh once (it doesn't change between trials)
    let mesh = rem_mesh::load_mesh(config, &NoComm)?;

    // Evaluate objective for a given parameter vector
    let eval = |x: &[f64]| -> RemResult<f64> {
        let clamped: Vec<f64> = x.iter().enumerate().map(|(i, &v)| {
            if let Some([lo, hi]) = bounds[i] { v.clamp(lo, hi) } else { v }
        }).collect();
        let cfg = apply_params(config, &params, &clamped)?;
        let mom_cfg = cfg.solver.mom.as_ref().ok_or_else(|| rem_core::RemError::Config(
            "MonteCarlo: Solver.MoM section required".to_string(),
        ))?;
        let matrices = rem_mom::compute_s_param_sweep_for_optim(&cfg, mom_cfg, &mesh)?;
        Ok(evaluate_objectives(&matrices, &par_cfg.objectives))
    };

    // Nominal objective (used as yield reference)
    let f_nom = eval(&p_nom)?;
    log::info!("MonteCarlo: nominal objective = {f_nom:.6e}");

    // Output directory
    let out_dir = std::path::Path::new(config.problem.output_dir()).join("postpro");
    std::fs::create_dir_all(&out_dir)?;

    let samples_path = out_dir.join("monte_carlo_samples.csv");
    let mut csv = std::fs::File::create(&samples_path)?;

    // Header
    let param_names: Vec<&str> = par_cfg.parameters.iter().map(|p| p.name.as_str()).collect();
    write!(csv, "Trial").map_err(rem_core::RemError::Io)?;
    for name in &param_names { write!(csv, ",{name}").map_err(rem_core::RemError::Io)?; }
    writeln!(csv, ",Objective").map_err(rem_core::RemError::Io)?;

    let mut rng = Xorshift64::new(seed);
    let mut objectives: Vec<f64> = Vec::with_capacity(n_samples);

    for trial in 0..n_samples {
        // Sample perturbed parameter vector
        let x: Vec<f64> = (0..n_params).map(|i| {
            p_nom[i] + sigmas[i] * rng.next_normal()
        }).collect();

        let f = eval(&x)?;
        objectives.push(f);

        // Write sample row
        write!(csv, "{trial}").map_err(rem_core::RemError::Io)?;
        for &v in &x { write!(csv, ",{v:.9e}").map_err(rem_core::RemError::Io)?; }
        writeln!(csv, ",{f:.9e}").map_err(rem_core::RemError::Io)?;

        if (trial + 1) % 10 == 0 || trial == 0 {
            log::info!("  MonteCarlo trial {}/{}: f = {f:.6e}", trial + 1, n_samples);
        }
    }

    // Compute statistics
    let mean   = objectives.iter().sum::<f64>() / n_samples as f64;
    let var    = objectives.iter().map(|&v| (v - mean).powi(2)).sum::<f64>() / n_samples as f64;
    let stddev = var.sqrt();
    let min    = objectives.iter().cloned().fold(f64::INFINITY, f64::min);
    let max    = objectives.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

    // Yield: fraction of samples within ±k·σ of f_nom
    let yield_1sigma = {
        let tol = stddev;
        objectives.iter().filter(|&&v| (v - f_nom).abs() <= tol).count() as f64
            / n_samples as f64
    };
    let yield_3sigma = {
        let tol = 3.0 * stddev;
        objectives.iter().filter(|&&v| (v - f_nom).abs() <= tol).count() as f64
            / n_samples as f64
    };

    log::info!(
        "MonteCarlo stats: mean={mean:.4e}  std={stddev:.4e}  min={min:.4e}  max={max:.4e}"
    );
    log::info!(
        "  Yield ±1σ={:.1}%  ±3σ={:.1}%",
        yield_1sigma * 100.0,
        yield_3sigma * 100.0,
    );

    // Write stats CSV
    let stats_path = out_dir.join("monte_carlo_stats.csv");
    let mut sf = std::fs::File::create(&stats_path)?;
    writeln!(sf, "Statistic,Value").map_err(rem_core::RemError::Io)?;
    writeln!(sf, "NominalObjective,{f_nom:.9e}").map_err(rem_core::RemError::Io)?;
    writeln!(sf, "NSamples,{n_samples}").map_err(rem_core::RemError::Io)?;
    writeln!(sf, "McSigmaRel,{sigma_rel:.6e}").map_err(rem_core::RemError::Io)?;
    writeln!(sf, "Mean,{mean:.9e}").map_err(rem_core::RemError::Io)?;
    writeln!(sf, "StdDev,{stddev:.9e}").map_err(rem_core::RemError::Io)?;
    writeln!(sf, "Min,{min:.9e}").map_err(rem_core::RemError::Io)?;
    writeln!(sf, "Max,{max:.9e}").map_err(rem_core::RemError::Io)?;
    writeln!(sf, "Yield_1sigma_pct,{:.4}", yield_1sigma * 100.0).map_err(rem_core::RemError::Io)?;
    writeln!(sf, "Yield_3sigma_pct,{:.4}", yield_3sigma * 100.0).map_err(rem_core::RemError::Io)?;

    log::info!("MonteCarlo samples → {}", samples_path.display());
    log::info!("MonteCarlo stats   → {}", stats_path.display());
    Ok(())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::Xorshift64;

    #[test]
    fn xorshift_non_zero_and_varies() {
        let mut rng = Xorshift64::new(42);
        let samples: Vec<u64> = (0..100).map(|_| rng.next_u64()).collect();
        assert!(samples.iter().all(|&v| v != 0));
        // At least 90 distinct values expected in 100 draws
        let unique: std::collections::HashSet<u64> = samples.iter().cloned().collect();
        assert!(unique.len() >= 90, "only {} unique values", unique.len());
    }

    #[test]
    fn box_muller_mean_and_std() {
        let mut rng = Xorshift64::new(0xABCD_1234);
        let n = 10_000usize;
        let samples: Vec<f64> = (0..n).map(|_| rng.next_normal()).collect();
        let mean = samples.iter().sum::<f64>() / n as f64;
        let var  = samples.iter().map(|&v| (v - mean).powi(2)).sum::<f64>() / n as f64;
        let std  = var.sqrt();
        // Mean ≈ 0 ± 0.05, std ≈ 1 ± 0.05
        assert!(mean.abs() < 0.05, "mean = {mean:.4}");
        assert!((std - 1.0).abs() < 0.05, "std = {std:.4}");
    }

    #[test]
    fn box_muller_f64_in_unit_interval() {
        let mut rng = Xorshift64::new(99);
        for _ in 0..1000 {
            let v = rng.next_f64();
            assert!(v > 0.0 && v <= 1.0, "out of (0,1]: {v}");
        }
    }

    #[test]
    fn yield_threshold_logic() {
        // Check that the yield formula gives sensible results.
        // 100% of samples within ±∞ of nominal should be 1.0.
        let objectives = vec![1.0_f64, 2.0, 3.0, 4.0, 5.0];
        let n = objectives.len() as f64;
        let mean = objectives.iter().sum::<f64>() / n;
        let var  = objectives.iter().map(|&v| (v - mean).powi(2)).sum::<f64>() / n;
        let std  = var.sqrt();
        let f_nom = mean;
        let y3 = objectives.iter().filter(|&&v| (v - f_nom).abs() <= 3.0 * std).count() as f64 / n;
        // All 5 samples within 3σ = 100%
        assert!((y3 - 1.0).abs() < 1e-9, "expected y3=1.0, got {y3}");
    }
}
