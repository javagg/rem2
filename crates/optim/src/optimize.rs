//! Derivative-free Nelder-Mead simplex optimizer for REM parametric optimization.
//!
//! Minimizes a scalar objective (sum of [`OptimObjective`] values) over the
//! named design parameters defined in [`ParametricConfig`].
//!
//! # Algorithm
//! Classic Nelder-Mead simplex with reflection/expansion/contraction/shrink
//! operations (Nelder & Mead 1965).  Convergence criterion: relative simplex
//! diameter < `tolerance`.
//!
//! # Output files
//! - `optimization_trace.csv` — iteration, objective value, best parameter values
//! - `optimization_result.json` — best parameters + objective at termination

use rem_config::{PalaceConfig, ParametricConfig};
use rem_core::RemResult;
use rem_parallel::NoComm;
use crate::objective::evaluate_objectives;
use crate::param_apply::apply_params;

// Nelder-Mead coefficients (standard values)
const NM_REFLECT:  f64 = 1.0;
const NM_EXPAND:   f64 = 2.0;
const NM_CONTRACT: f64 = 0.5;
const NM_SHRINK:   f64 = 0.5;

/// Run Nelder-Mead optimization and write output files.
pub fn run_optimize(config: &PalaceConfig, par_cfg: &ParametricConfig) -> RemResult<()> {
    if par_cfg.objectives.is_empty() {
        return Err(rem_core::RemError::Config(
            "Parametric optimizer: at least one Objectives entry required".to_string(),
        ));
    }
    if par_cfg.parameters.is_empty() {
        return Err(rem_core::RemError::Config(
            "Parametric optimizer: at least one Parameters entry required".to_string(),
        ));
    }

    let output_dir = std::path::Path::new(config.problem.output_dir());
    std::fs::create_dir_all(output_dir)?;

    let _n = par_cfg.parameters.len();
    let bounds: Vec<Option<[f64; 2]>> = par_cfg.parameters.iter().map(|p| p.bounds).collect();

    // Build initial simplex from `initial` values (falling back to midpoint of bounds or 1.0).
    let x0: Vec<f64> = par_cfg.parameters.iter().map(|p| {
        p.initial.unwrap_or_else(|| {
            if let Some([lo, hi]) = p.bounds { (lo + hi) * 0.5 } else { 1.0 }
        })
    }).collect();

    let mesh = rem_mesh::load_mesh(config, &NoComm)?;
    let mut n_evals: usize = 0;

    // Objective closure: clamp to bounds, patch config, run sweep, evaluate.
    let mut eval = |x: &[f64]| -> RemResult<f64> {
        let clamped: Vec<f64> = x.iter().enumerate().map(|(i, &v)| {
            if let Some([lo, hi]) = bounds[i] { v.clamp(lo, hi) } else { v }
        }).collect();
        let params: Vec<&rem_config::SweepParam> = par_cfg.parameters.iter().collect();
        let cfg = apply_params(config, &params, &clamped)?;
        let mom_cfg = cfg.solver.mom.as_ref().ok_or_else(|| rem_core::RemError::Config(
            "Optimizer: Solver.MoM section required".to_string(),
        ))?;
        let matrices = rem_mom::compute_s_param_sweep_for_optim(&cfg, mom_cfg, &mesh)?;
        n_evals += 1;
        Ok(evaluate_objectives(&matrices, &par_cfg.objectives))
    };

    // Nelder-Mead
    let (best_x, best_f, history) = nelder_mead(
        &x0, &mut eval, par_cfg.max_iter, par_cfg.tolerance,
    )?;

    log::info!(
        "Nelder-Mead converged after {} evaluations: objective={:.6e}",
        n_evals, best_f
    );
    for (i, p) in par_cfg.parameters.iter().enumerate() {
        log::info!("  {} = {:.6e}", p.name, best_x[i]);
    }

    // Write trace CSV
    write_trace_csv(output_dir, &history, par_cfg)?;

    // Write result JSON
    write_result_json(output_dir, par_cfg, &best_x, best_f, n_evals)?;

    println!(
        "[rem-optim] Optimization complete: objective = {:.6e}\n  {}",
        best_f,
        best_x.iter().zip(par_cfg.parameters.iter())
            .map(|(v, p)| format!("{}={:.4e}", p.name, v))
            .collect::<Vec<_>>().join(", ")
    );
    Ok(())
}

/// One Nelder-Mead history record.
struct HistoryEntry {
    iter: usize,
    best_f: f64,
    best_x: Vec<f64>,
}

/// Nelder-Mead simplex minimization.
///
/// Returns `(best_x, best_f, history)`.
fn nelder_mead<F>(
    x0: &[f64],
    f: &mut F,
    max_iter: usize,
    tol: f64,
) -> RemResult<(Vec<f64>, f64, Vec<HistoryEntry>)>
where
    F: FnMut(&[f64]) -> RemResult<f64>,
{
    let n = x0.len();

    // Initial simplex: x0, x0 + h·eᵢ for i=0..n-1
    let mut simplex: Vec<Vec<f64>> = Vec::with_capacity(n + 1);
    simplex.push(x0.to_vec());
    for i in 0..n {
        let mut v = x0.to_vec();
        let step = if v[i].abs() > 1e-10 { 0.05 * v[i].abs() } else { 0.00025 };
        v[i] += step;
        simplex.push(v);
    }

    // Evaluate at all simplex vertices.
    let mut fvals: Vec<f64> = simplex.iter().map(|x| f(x)).collect::<RemResult<_>>()?;

    let mut history: Vec<HistoryEntry> = Vec::new();

    for iter in 0..max_iter {
        // Sort by function value.
        let mut order: Vec<usize> = (0..=n).collect();
        order.sort_by(|&a, &b| fvals[a].partial_cmp(&fvals[b]).unwrap_or(std::cmp::Ordering::Equal));
        let best  = order[0];
        let worst = order[n];
        let second_worst = order[n - 1];

        history.push(HistoryEntry {
            iter,
            best_f: fvals[best],
            best_x: simplex[best].clone(),
        });

        // Convergence check: relative simplex diameter
        let diameter: f64 = (1..=n).map(|k| {
            let i = order[k];
            simplex[i].iter().zip(simplex[best].iter())
                .map(|(a, b)| (a - b).powi(2)).sum::<f64>().sqrt()
        }).fold(0.0_f64, f64::max);
        let scale = simplex[best].iter().map(|x| x.abs()).fold(1.0_f64, f64::max);
        if diameter / scale < tol {
            break;
        }

        // Centroid of all but worst.
        let centroid: Vec<f64> = (0..n).map(|d| {
            order[..n].iter().map(|&i| simplex[i][d]).sum::<f64>() / n as f64
        }).collect();

        // Reflection.
        let reflected: Vec<f64> = centroid.iter().zip(simplex[worst].iter())
            .map(|(c, w)| c + NM_REFLECT * (c - w)).collect();
        let f_reflected = f(&reflected)?;

        if f_reflected < fvals[best] {
            // Expansion.
            let expanded: Vec<f64> = centroid.iter().zip(reflected.iter())
                .map(|(c, r)| c + NM_EXPAND * (r - c)).collect();
            let f_expanded = f(&expanded)?;
            if f_expanded < f_reflected {
                simplex[worst] = expanded;
                fvals[worst] = f_expanded;
            } else {
                simplex[worst] = reflected;
                fvals[worst] = f_reflected;
            }
        } else if f_reflected < fvals[second_worst] {
            simplex[worst] = reflected;
            fvals[worst] = f_reflected;
        } else {
            // Contraction.
            let contracted: Vec<f64> = centroid.iter().zip(simplex[worst].iter())
                .map(|(c, w)| c + NM_CONTRACT * (w - c)).collect();
            let f_contracted = f(&contracted)?;
            if f_contracted < fvals[worst] {
                simplex[worst] = contracted;
                fvals[worst] = f_contracted;
            } else {
                // Shrink.
                for k in 1..=n {
                    let i = order[k];
                    let shrunk: Vec<f64> = simplex[best].iter().zip(simplex[i].iter())
                        .map(|(b, v)| b + NM_SHRINK * (v - b)).collect();
                    fvals[i] = f(&shrunk)?;
                    simplex[i] = shrunk;
                }
            }
        }
    }

    // Final sort.
    let mut order: Vec<usize> = (0..=n).collect();
    order.sort_by(|&a, &b| fvals[a].partial_cmp(&fvals[b]).unwrap_or(std::cmp::Ordering::Equal));
    let best_x = simplex[order[0]].clone();
    let best_f = fvals[order[0]];

    Ok((best_x, best_f, history))
}

fn write_trace_csv(
    output_dir: &std::path::Path,
    history: &[HistoryEntry],
    par_cfg: &ParametricConfig,
) -> RemResult<()> {
    use std::io::Write;
    let path = output_dir.join("optimization_trace.csv");
    let mut f = std::fs::File::create(&path)?;

    let mut header = vec!["iter".to_string(), "best_f".to_string()];
    for p in &par_cfg.parameters { header.push(p.name.clone()); }
    writeln!(f, "{}", header.join(","))?;

    for h in history {
        let mut parts = vec![h.iter.to_string(), format!("{:.6e}", h.best_f)];
        for &v in &h.best_x { parts.push(format!("{:.6e}", v)); }
        writeln!(f, "{}", parts.join(","))?;
    }
    Ok(())
}

fn write_result_json(
    output_dir: &std::path::Path,
    par_cfg: &ParametricConfig,
    best_x: &[f64],
    best_f: f64,
    n_evals: usize,
) -> RemResult<()> {
    let mut obj = serde_json::Map::new();
    obj.insert("objective".to_string(), serde_json::json!(best_f));
    obj.insert("n_evaluations".to_string(), serde_json::json!(n_evals));
    let mut params = serde_json::Map::new();
    for (i, p) in par_cfg.parameters.iter().enumerate() {
        params.insert(p.name.clone(), serde_json::json!(best_x[i]));
    }
    obj.insert("best_parameters".to_string(), serde_json::Value::Object(params));
    let json = serde_json::to_string_pretty(&serde_json::Value::Object(obj))
        .map_err(|e| rem_core::RemError::Io(std::io::Error::other(e.to_string())))?;
    std::fs::write(output_dir.join("optimization_result.json"), json)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nelder_mead_minimizes_quadratic() {
        // f(x,y) = (x-2)² + (y+3)² with minimum at (2,-3).
        let x0 = vec![0.0f64, 0.0];
        let mut eval = |x: &[f64]| -> RemResult<f64> {
            Ok((x[0] - 2.0).powi(2) + (x[1] + 3.0).powi(2))
        };
        let (best, fval, _) = nelder_mead(&x0, &mut eval, 500, 1e-8).unwrap();
        assert!((best[0] - 2.0).abs() < 1e-3, "x={:.6}", best[0]);
        assert!((best[1] + 3.0).abs() < 1e-3, "y={:.6}", best[1]);
        assert!(fval < 1e-6, "f={:.2e}", fval);
    }

    #[test]
    fn nelder_mead_1d_minimizes() {
        let mut eval = |x: &[f64]| -> RemResult<f64> { Ok((x[0] - 5.0).powi(2)) };
        let (best, fval, _) = nelder_mead(&[0.0], &mut eval, 200, 1e-8).unwrap();
        assert!((best[0] - 5.0).abs() < 1e-3, "x={:.6}", best[0]);
        assert!(fval < 1e-5, "f={:.2e}", fval);
    }
}
