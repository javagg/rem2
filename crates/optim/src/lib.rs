//! rem-optim — Parametric sweep and gradient optimization for REM solvers.
//!
//! # Overview
//!
//! This crate wraps the MoM S-parameter solver with two outer-loop strategies:
//!
//! - **Sweep** (`Mode: "Sweep"`): Exhaustive Cartesian-product grid over named
//!   design parameters.  Writes `parametric_sweep.csv`.
//!
//! - **Optimize** (`Mode: "Optimize"`): Derivative-free Nelder-Mead simplex
//!   minimization of a user-defined objective (e.g., `MinS11dB` at a target
//!   frequency).  Writes `optimization_trace.csv` and
//!   `optimization_result.json`.
//!
//! # JSON config example (Sweep)
//! ```json
//! "Parametric": {
//!   "Mode": "Sweep",
//!   "Parameters": [
//!     {"Name": "eps_r", "Target": {"Type": "SubstratePermittivity", "Layer": 0},
//!      "Min": 3.0, "Max": 5.0, "Steps": 5},
//!     {"Name": "h_mm",  "Target": {"Type": "SubstrateThickness",    "Layer": 0},
//!      "Values": [0.2e-3, 0.3e-3, 0.5e-3]}
//!   ]
//! }
//! ```
//!
//! # JSON config example (Optimize)
//! ```json
//! "Parametric": {
//!   "Mode": "Optimize",
//!   "Parameters": [
//!     {"Name": "eps_r", "Target": {"Type": "SubstratePermittivity", "Layer": 0},
//!      "Initial": 4.0, "Bounds": [3.0, 5.0]}
//!   ],
//!   "Objectives": [{"Type": "MinS11dB", "Port": 1, "FreqHz": 2.4e9}],
//!   "MaxIter": 200,
//!   "Tolerance": 1e-4
//! }
//! ```

pub mod sweep;
pub mod optimize;
pub mod objective;
pub mod param_apply;
pub mod sensitivity;
pub mod monte_carlo;

use rem_config::{PalaceConfig, ParametricMode};
use rem_core::RemResult;

/// Entry point called from the CLI when `Solver.Parametric` is present.
///
/// Dispatches to [`sweep::run_sweep`] or [`optimize::run_optimize`] based on
/// `config.solver.parametric.mode`.
pub fn run(config: &PalaceConfig) -> RemResult<()> {
    let par_cfg = config
        .solver
        .parametric
        .as_ref()
        .ok_or_else(|| rem_core::RemError::Config(
            "Parametric solver: Solver.Parametric section is missing".to_string(),
        ))?;

    match par_cfg.mode {
        ParametricMode::Sweep       => sweep::run_sweep(config, par_cfg),
        ParametricMode::Optimize    => optimize::run_optimize(config, par_cfg),
        ParametricMode::Sensitivity => sensitivity::run_sensitivity(config, par_cfg),
        ParametricMode::MonteCarlo  => monte_carlo::run_monte_carlo(config, par_cfg),
    }
}
