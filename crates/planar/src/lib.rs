//! `rem-planar` — 平面分层介质矩量法（Planar MoM + FFT）求解器
//!
//! 对标 Sonnet Suite，支持：
//! - 分层媒质格林函数（谱域传递矩阵法）
//! - 均匀网格 2D FFT 卷积加速 O(N log N)
//! - 直接求解 / 最速下降迭代

pub mod layered_green;
pub mod grid;
pub mod fft_conv;
pub mod impedance;
pub mod solver;

pub use layered_green::{LayeredMedium, Layer, SpectralGreen};

use rem_config::PalaceConfig;
use rem_core::{RemResult, RemError};

/// Run the Planar MoM solver from a Palace config.
///
/// Expects a `Solver.Planar` section with grid dimensions (Lx, Ly, Nx, Ny),
/// frequency range (FreqMin, FreqMax, FreqStep), and optional substrate layers.
pub fn run(config: &PalaceConfig) -> RemResult<()> {
    let planar_cfg = config.solver.planar.as_ref()
        .ok_or_else(|| RemError::Config(
            "Problem.Type = \"Planar\" requires a Solver.Planar section".to_string()
        ))?;

    if planar_cfg.lx <= 0.0 || planar_cfg.ly <= 0.0 {
        return Err(RemError::Config("Planar domain Lx/Ly must be positive".into()));
    }

    log::info!("\n=== Planar MoM solver ===\n");
    log::info!("Domain: {} × {} m, grid: {}×{}", planar_cfg.lx, planar_cfg.ly, planar_cfg.nx, planar_cfg.ny);

    // Build uniform planar grid
    let grid = grid::PlanarGrid::new(planar_cfg.lx, planar_cfg.ly, planar_cfg.nx, planar_cfg.ny);
    let n_basis = grid.edges.len();
    log::info!("RWG basis functions (edges): {}", n_basis);

    if n_basis == 0 {
        return Err(RemError::Mesh("Planar grid produced no edges".into()));
    }

    let output_dir = std::path::Path::new(config.problem.output_dir());
    #[cfg(not(target_arch = "wasm32"))]
    std::fs::create_dir_all(output_dir)?;

    // Solver config
    let solver_cfg = solver::SolverConfig {
        use_fft: false,
        max_iter: 0, // direct solve
        tol: 1e-6,
    };
    let mom_solver = solver::PlanarMomSolver::new(grid, solver_cfg);

    // Frequency sweep
    let f_min = planar_cfg.freq_min;
    let f_max = planar_cfg.freq_max;
    let f_step = planar_cfg.freq_step;

    if f_step <= 0.0 || f_max < f_min {
        return Err(RemError::Config("Invalid Planar frequency range".into()));
    }

    let n_freqs = ((f_max - f_min) / f_step).ceil() as usize + 1;

    #[allow(unused_imports)]
    use num_complex::Complex64;

    // Unit excitation vector
    let excitation = nalgebra::DVector::from_element(n_basis, Complex64::new(1.0, 0.0));

    #[cfg(not(target_arch = "wasm32"))]
    let mut wtr = csv::Writer::from_path(output_dir.join("planar-s.csv"))
        .map_err(|e| RemError::Io(e.into()))?;

    #[cfg(not(target_arch = "wasm32"))]
    wtr.write_record(&["freq_hz", "s11_re", "s11_im"])
        .map_err(|e| RemError::Io(e.into()))?;

    for i in 0..n_freqs {
        let freq = f_min + i as f64 * f_step;
        if freq > f_max {
            break;
        }

        log::info!("  freq = {:.3e} Hz", freq);

        let solution = mom_solver.solve(freq, &excitation);

        // Estimate S11 from the first current coefficient
        let i0 = solution.coefficients[0];
        let z0 = planar_cfg.ref_impedance;
        let z_in = if i0.norm() > 1e-30 {
            Complex64::new(1.0, 0.0) / i0
        } else {
            Complex64::new(1e12, 0.0)
        };
        let s11 = (z_in - Complex64::new(z0, 0.0)) / (z_in + Complex64::new(z0, 0.0));

        log::info!("    S11 = {:.4}∠{:.1}°", s11.norm(), s11.arg().to_degrees());

        #[cfg(not(target_arch = "wasm32"))]
        wtr.write_record(&[
            format!("{:.6e}", freq),
            format!("{:.6e}", s11.re),
            format!("{:.6e}", s11.im),
        ]).map_err(|e| RemError::Io(e.into()))?;
    }

    #[cfg(not(target_arch = "wasm32"))]
    wtr.flush().map_err(|e| RemError::Io(e.into()))?;

    log::info!("\nPlanar MoM solve completed ({} frequency points).\n", n_freqs);
    Ok(())
}
