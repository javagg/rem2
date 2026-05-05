//! Finite-difference gradient (sensitivity) analysis for parametric objectives.
//!
//! Given a nominal design point, computes the gradient ∂f/∂pᵢ of the objective
//! with respect to each design parameter using central finite differences:
//!
//!   ∂f/∂pᵢ ≈ [f(p + h·eᵢ) − f(p − h·eᵢ)] / (2h)
//!
//! where h is a relative step size (default 1 %).
//!
//! # Output files
//! - `postpro/sensitivity.csv` — one row per parameter with the gradient and
//!   its normalised form `|∂f/∂pᵢ| · |pᵢ| / |f|`.
//!
//! # Usage (JSON config)
//! ```json
//! "Parametric": {
//!   "Mode": "Sensitivity",
//!   "Parameters": [
//!     {"Name": "eps_r", "Target": {"Type": "SubstratePermittivity", "Layer": 0},
//!      "Initial": 4.0, "Bounds": [3.0, 5.0]}
//!   ],
//!   "Objectives": [{"Type": "MinS11dB", "Port": 1, "FreqHz": 2.4e9}],
//!   "SensRelStep": 0.01
//! }
//! ```

use rem_config::{PalaceConfig, ParametricConfig};
use rem_core::RemResult;
use rem_parallel::NoComm;

use crate::objective::evaluate_objectives;
use crate::param_apply::apply_params;

/// Run finite-difference sensitivity analysis at the nominal design point.
///
/// For each parameter pᵢ the function computes:
/// - `grad[i]` = ∂f/∂pᵢ  (central differences)
/// - `norm_sens[i]` = `|grad[i]| * |p_nominal[i]| / max(|f_nominal|, 1e-30))`
///
/// Results are written to `{output_dir}/postpro/sensitivity.csv`.
pub fn run_sensitivity(config: &PalaceConfig, par_cfg: &ParametricConfig) -> RemResult<()> {
    use std::io::Write;

    if par_cfg.objectives.is_empty() {
        return Err(rem_core::RemError::Config(
            "Sensitivity: at least one Objectives entry required".to_string(),
        ));
    }
    if par_cfg.parameters.is_empty() {
        return Err(rem_core::RemError::Config(
            "Sensitivity: at least one Parameters entry required".to_string(),
        ));
    }

    let rel_step = par_cfg.sens_rel_step.unwrap_or(0.01).max(1e-8);

    // Nominal parameter values
    let p_nom: Vec<f64> = par_cfg.parameters.iter().map(|p| {
        p.initial.unwrap_or_else(|| {
            if let Some([lo, hi]) = p.bounds { (lo + hi) * 0.5 } else { 1.0 }
        })
    }).collect();
    let params: Vec<&rem_config::SweepParam> = par_cfg.parameters.iter().collect();
    let bounds: Vec<Option<[f64; 2]>> = par_cfg.parameters.iter().map(|p| p.bounds).collect();

    // Helper: clamp to bounds, build config, evaluate objective.
    let mesh = rem_mesh::load_mesh(config, &NoComm)?;
    let mut eval = |x: &[f64]| -> RemResult<f64> {
        let clamped: Vec<f64> = x.iter().enumerate().map(|(i, &v)| {
            if let Some([lo, hi]) = bounds[i] { v.clamp(lo, hi) } else { v }
        }).collect();
        let cfg = apply_params(config, &params, &clamped)?;
        let mom_cfg = cfg.solver.mom.as_ref().ok_or_else(|| rem_core::RemError::Config(
            "Sensitivity: Solver.MoM section required".to_string(),
        ))?;
        let matrices = rem_mom::compute_s_param_sweep_for_optim(&cfg, mom_cfg, &mesh)?;
        Ok(evaluate_objectives(&matrices, &par_cfg.objectives))
    };

    // Nominal objective value
    let f_nom = eval(&p_nom)?;
    log::info!("Sensitivity: nominal objective = {f_nom:.6e}");

    let n = par_cfg.parameters.len();
    let mut gradients: Vec<f64>   = vec![0.0; n];
    let mut norm_sens: Vec<f64>   = vec![0.0; n];

    for i in 0..n {
        let pi = p_nom[i];
        let h  = if pi.abs() > 1e-30 { rel_step * pi.abs() } else { rel_step };

        let mut x_plus  = p_nom.clone();
        let mut x_minus = p_nom.clone();
        x_plus[i]  = pi + h;
        x_minus[i] = pi - h;

        let f_plus  = eval(&x_plus)?;
        let f_minus = eval(&x_minus)?;
        let grad_i  = (f_plus - f_minus) / (2.0 * h);

        gradients[i] = grad_i;
        norm_sens[i] = grad_i.abs() * pi.abs() / f_nom.abs().max(1e-30);

        log::info!(
            "  ∂f/∂{} = {:.6e}  (norm.sens={:.4})",
            par_cfg.parameters[i].name,
            grad_i,
            norm_sens[i],
        );
    }

    // Write CSV
    let out_dir = std::path::Path::new(config.problem.output_dir()).join("postpro");
    std::fs::create_dir_all(&out_dir)?;
    let path = out_dir.join("sensitivity.csv");
    let mut f = std::fs::File::create(&path)?;
    writeln!(f, "Parameter,NominalValue,Gradient (df/dp),NormSensitivity").map_err(rem_core::RemError::Io)?;
    for i in 0..n {
        writeln!(f, "{},{:.9e},{:.9e},{:.6e}",
            par_cfg.parameters[i].name,
            p_nom[i],
            gradients[i],
            norm_sens[i],
        ).map_err(rem_core::RemError::Io)?;
    }
    log::info!("Sensitivity results written to {}", path.display());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Gradient of f(x) = x² at x=3 via central FD: should be ≈ 6.
    #[test]
    fn central_fd_quadratic() {
        let x = 3.0_f64;
        let rel_step = 0.01_f64;
        let h = rel_step * x.abs();
        let f = |v: f64| v * v;
        let grad = (f(x + h) - f(x - h)) / (2.0 * h);
        assert!((grad - 6.0).abs() < 1e-6, "grad={grad:.8}, expected ≈ 6.0");
    }

    /// Normalised sensitivity of f(x)=x² at x=3: |grad|*|x|/|f| = 6*3/9 = 2.
    #[test]
    fn normalised_sensitivity_quadratic() {
        let x = 3.0_f64;
        let f_nom = x * x;  // 9
        let grad  = 2.0 * x; // 6  (exact derivative)
        let norm = grad.abs() * x.abs() / f_nom.abs();
        assert!((norm - 2.0).abs() < 1e-12, "norm_sens={norm:.6}, expected 2.0");
    }
}
