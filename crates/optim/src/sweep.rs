//! Parametric grid sweep: Cartesian product over all parameter value lists.
//!
//! For each combination of parameter values, patches the base config, runs
//! the MoM S-parameter sweep, and records the result in a CSV table.
//!
//! Output file: `<output_dir>/parametric_sweep.csv`
//!
//! CSV columns: `param_name_1, param_name_2, ..., freq_hz, S11_dB, S21_dB, ...`

use rem_config::{PalaceConfig, ParametricConfig};
use rem_core::RemResult;
use rem_mesh::RemMesh;
use rem_parallel::NoComm;
use crate::param_apply::apply_params;
use crate::objective::evaluate_objectives;

/// Run an exhaustive grid sweep and write `parametric_sweep.csv`.
pub fn run_sweep(config: &PalaceConfig, par_cfg: &ParametricConfig) -> RemResult<()> {
    let output_dir = std::path::Path::new(config.problem.output_dir());
    std::fs::create_dir_all(output_dir)?;

    let base_mesh = rem_mesh::load_mesh(config, &NoComm)?;

    let mut records: Vec<SweepRecord> = Vec::new();
    run_sweep_with_mesh(config, par_cfg, &base_mesh, &mut records)?;

    write_sweep_csv(&records, par_cfg, output_dir)?;
    log::info!(
        "Parametric sweep: {} combinations → {}",
        records.len(),
        output_dir.join("parametric_sweep.csv").display()
    );
    Ok(())
}

/// Record for one grid point: parameter values + per-frequency S-matrix data.
pub struct SweepRecord {
    /// Values of each parameter at this grid point.
    pub param_values: Vec<f64>,
    /// Objective value (sum of objectives, NaN if objectives not configured).
    pub objective: f64,
    /// Per-frequency rows: (freq_hz, S-matrix data in magnitude dB).
    pub freq_rows: Vec<FreqRow>,
}

pub struct FreqRow {
    pub freq_hz: f64,
    /// S-mag dB entries row-major: [S11, S12, ..., S21, S22, ...].
    pub s_mag_db: Vec<f64>,
}

fn run_sweep_with_mesh(
    config: &PalaceConfig,
    par_cfg: &ParametricConfig,
    mesh: &RemMesh,
    records: &mut Vec<SweepRecord>,
) -> RemResult<()> {
    // Build value lists per parameter.
    let all_values: Vec<Vec<f64>> = par_cfg.parameters.iter()
        .map(|p| p.resolved_values())
        .collect();

    // Verify all parameters have at least one value.
    for (i, vals) in all_values.iter().enumerate() {
        if vals.is_empty() {
            return Err(rem_core::RemError::Config(format!(
                "Parametric sweep: parameter '{}' has no values (set Values or Min/Max/Steps)",
                par_cfg.parameters[i].name
            )));
        }
    }

    // Enumerate Cartesian product via index vector.
    let n_params = par_cfg.parameters.len();
    let mut indices = vec![0usize; n_params];
    let counts: Vec<usize> = all_values.iter().map(|v| v.len()).collect();

    loop {
        // Collect current parameter values.
        let values: Vec<f64> = indices.iter().enumerate().map(|(i, &idx)| all_values[i][idx]).collect();
        let params: Vec<&rem_config::SweepParam> = par_cfg.parameters.iter().collect();

        // Patch config and run MoM sweep.
        let cfg = apply_params(config, &params, &values)?;
        let mom_cfg = cfg.solver.mom.as_ref().ok_or_else(|| rem_core::RemError::Config(
            "Parametric sweep: Solver.MoM section required".to_string(),
        ))?;

        log::info!(
            "Parametric sweep [{}]: {}",
            indices.iter().map(|i| i.to_string()).collect::<Vec<_>>().join(","),
            values.iter().zip(par_cfg.parameters.iter()).map(|(v, p)| format!("{}={:.4}", p.name, v)).collect::<Vec<_>>().join(", ")
        );

        let matrices = rem_mom::compute_s_param_sweep_for_optim(&cfg, mom_cfg, mesh)?;

        // Compute objective if configured.
        let objective = if !par_cfg.objectives.is_empty() {
            evaluate_objectives(&matrices, &par_cfg.objectives)
        } else {
            f64::NAN
        };

        // Build frequency rows.
        let freq_rows: Vec<FreqRow> = matrices.iter().map(|sm| {
            let n = sm.n_ports;
            let s_mag_db: Vec<f64> = (0..n).flat_map(|row| {
                (0..n).map(move |col| {
                    let s = sm.get(row, col);
                    let mag = s.norm();
                    if mag < 1e-30 { -300.0 } else { 20.0 * mag.log10() }
                })
            }).collect();
            FreqRow { freq_hz: sm.freq_hz, s_mag_db }
        }).collect();

        records.push(SweepRecord { param_values: values, objective, freq_rows });

        // Advance index vector (little-endian).
        let mut carry = true;
        for i in (0..n_params).rev() {
            if carry {
                indices[i] += 1;
                if indices[i] >= counts[i] {
                    indices[i] = 0;
                } else {
                    carry = false;
                }
            }
        }
        if carry { break; } // All combinations exhausted
    }
    Ok(())
}

fn write_sweep_csv(
    records: &[SweepRecord],
    par_cfg: &ParametricConfig,
    output_dir: &std::path::Path,
) -> RemResult<()> {
    use std::io::Write;

    let path = output_dir.join("parametric_sweep.csv");
    let mut f = std::fs::File::create(&path)?;

    // Determine n_ports from first record.
    let n_ports = records.first()
        .and_then(|r| r.freq_rows.first())
        .map(|row| {
            let n2 = row.s_mag_db.len();
            (n2 as f64).sqrt().round() as usize
        })
        .unwrap_or(0);

    // Header
    let param_headers: Vec<String> = par_cfg.parameters.iter().map(|p| p.name.clone()).collect();
    let mut header_parts = param_headers.clone();
    header_parts.push("freq_hz".to_string());
    for i in 0..n_ports {
        for j in 0..n_ports {
            header_parts.push(format!("S{}{}_dB", i + 1, j + 1));
        }
    }
    if !par_cfg.objectives.is_empty() {
        header_parts.push("objective".to_string());
    }
    writeln!(f, "{}", header_parts.join(","))?;

    // Data rows
    for record in records {
        for row in &record.freq_rows {
            let mut parts: Vec<String> = record.param_values.iter().map(|v| format!("{:.6e}", v)).collect();
            parts.push(format!("{:.6e}", row.freq_hz));
            for &db in &row.s_mag_db {
                parts.push(format!("{:.4}", db));
            }
            if !par_cfg.objectives.is_empty() {
                parts.push(format!("{:.6e}", record.objective));
            }
            writeln!(f, "{}", parts.join(","))?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use rem_config::{SweepParam, ParamTarget, ParametricMode};

    fn make_sweep_param(name: &str, values: Vec<f64>) -> SweepParam {
        SweepParam {
            name: name.to_string(),
            target: ParamTarget::FreqMin,
            values,
            min: None, max: None, steps: None,
            initial: None, bounds: None,
        }
    }

    #[test]
    fn resolved_values_explicit() {
        let p = make_sweep_param("a", vec![1.0, 2.0, 3.0]);
        assert_eq!(p.resolved_values(), vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn resolved_values_range() {
        let p = SweepParam {
            name: "a".to_string(),
            target: ParamTarget::FreqMin,
            values: vec![],
            min: Some(0.0), max: Some(1.0), steps: Some(3),
            initial: None, bounds: None,
        };
        let v = p.resolved_values();
        assert_eq!(v.len(), 3);
        assert!((v[0] - 0.0).abs() < 1e-12);
        assert!((v[1] - 0.5).abs() < 1e-12);
        assert!((v[2] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn cartesian_product_count() {
        // Two params: 3 × 4 values → 12 combos.
        // We test just the counting logic here without running a real solver.
        let a_vals = vec![1.0f64, 2.0, 3.0];
        let b_vals = vec![10.0f64, 20.0, 30.0, 40.0];
        let mut count = 0usize;
        let counts = [a_vals.len(), b_vals.len()];
        let mut indices = [0usize; 2];
        loop {
            count += 1;
            let mut carry = true;
            for i in (0..2).rev() {
                if carry {
                    indices[i] += 1;
                    if indices[i] >= counts[i] { indices[i] = 0; } else { carry = false; }
                }
            }
            if carry { break; }
        }
        assert_eq!(count, 12);
    }
}
