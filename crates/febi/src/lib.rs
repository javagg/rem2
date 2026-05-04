//! rem-febi — Hybrid Finite Element – Boundary Integral (FE-BI) solver

pub mod hybrid_mesh;
pub mod calderon;
pub mod coupling;
pub mod solver;
pub mod postprocess;

use rem_config::{PalaceConfig, FeBiSolverConfig};
use rem_core::RemResult;
use rem_mesh::RemMesh;
use rem_parallel::NoComm;

/// FE-BI 求解结果
#[derive(Debug, Clone)]
pub struct FebiResult {
    pub sparams: Vec<(f64, Vec<num_complex::Complex64>)>,
}

/// CLI 入口
pub fn run(config: &PalaceConfig) -> RemResult<()> {
    let febi_cfg = config.solver.febi.as_ref()
        .ok_or_else(|| rem_core::RemError::Config(
            "Problem.Type = \"FEBI\" requires a Solver.FEBI section".to_string()
        ))?;
    let mesh = rem_mesh::load_mesh(config, &NoComm)?;
    run_with_mesh(config, febi_cfg, &mesh).map(|_| ())
}

/// 在已加载网格上运行 FE-BI 求解
pub fn run_with_mesh(
    _config: &PalaceConfig,
    febi_cfg: &FeBiSolverConfig,
    mesh: &RemMesh,
) -> RemResult<FebiResult> {
    use std::f64::consts::PI;

    log::info!("\n=== FE-BI (Finite Element - Boundary Integral) solver ===\n");

    let freqs = build_freq_list(febi_cfg);
    log::info!("Frequency sweep:");
    log::info!("  {} frequencies", freqs.len());
    log::info!("");

    let surf = hybrid_mesh::extract_radiation_boundary(mesh, &febi_cfg.radiation_boundary)?;
    log::info!("Radiation boundary:");
    log::info!("  {} triangular faces", surf.faces.len());
    log::info!("  {} RWG basis functions", surf.edges.len());
    log::info!("");

    let mut sparams_out = Vec::new();

    for &freq in &freqs {
        let omega = 2.0 * PI * freq;
        log::info!("  Computing f = {:.4e} Hz (ω = {:.4e} rad/s)...", freq, omega);

        // 1. 组装 Calderón BI 矩阵
        let z_bi = calderon::assemble_calderon(&surf, freq, febi_cfg.aca_tol)?;

        // 2. 组装混合 FEM-BI 系统（返回 FebiSystem）
        let system = coupling::assemble_febi_system(
            febi_cfg, mesh, &surf, &z_bi, freq,
        )?;

        // 3. 求解（传入 &system.z_bi 和 &system.rhs）
        let solution = solver::solve_febi(&system.z_bi, &system.rhs)?;

        // 4. S 参数后处理（传入 solution 切片、端口数量、freq）
        let n_ports = febi_cfg.ports.len().max(1);
        let sp = postprocess::extract_sparams(solution.as_slice(), n_ports, freq)?;
        sparams_out.push((freq, sp));
    }

    Ok(FebiResult { sparams: sparams_out })
}

fn build_freq_list(cfg: &FeBiSolverConfig) -> Vec<f64> {
    if cfg.freq_step <= 0.0 || cfg.freq_min >= cfg.freq_max {
        return vec![cfg.freq_min];
    }
    let mut freqs = Vec::new();
    let mut f = cfg.freq_min;
    while f <= cfg.freq_max + cfg.freq_step * 1e-6 {
        freqs.push(f);
        f += cfg.freq_step;
    }
    freqs
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_build_freq_list_single() {
        let cfg = FeBiSolverConfig {
            freq_min: 1e9,
            freq_max: 1e9,
            freq_step: 0.0,
            radiation_boundary: vec![],
            equation: "CFIE".to_string(),
            alpha: 0.5,
            aca_tol: 1e-3,
            gmres_tol: 1e-6,
            gmres_max_iter: 500,
            ports: vec![],
            ref_impedance: 50.0,
            exterior_eps_r: 1.0,
            exterior_mu_r: 1.0,
            output_dir: "output/febi".to_string(),
        };
        let freqs = build_freq_list(&cfg);
        assert_eq!(freqs.len(), 1);
    }

    #[test]
    fn test_build_freq_list_sweep() {
        let cfg = FeBiSolverConfig {
            freq_min: 1e9,
            freq_max: 3e9,
            freq_step: 1e9,
            radiation_boundary: vec![],
            equation: "CFIE".to_string(),
            alpha: 0.5,
            aca_tol: 1e-3,
            gmres_tol: 1e-6,
            gmres_max_iter: 500,
            ports: vec![],
            ref_impedance: 50.0,
            exterior_eps_r: 1.0,
            exterior_mu_r: 1.0,
            output_dir: "output/febi".to_string(),
        };
        let freqs = build_freq_list(&cfg);
        assert_eq!(freqs.len(), 3);
    }
}
