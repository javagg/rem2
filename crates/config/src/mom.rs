//! Standalone MoM config format for `rem-suite run mom <config.json>`.
//!
//! JSON is flat — no `Solver.MoM` nesting.  Example:
//!
//! ```json
//! {
//!   "Output": "./output",
//!   "Mesh": "mesh.msh",
//!   "Equation": "CFIE",
//!   "Basis": "RWG",
//!   "FreqMin": 1e9,
//!   "FreqMax": 10e9,
//!   "FreqStep": 0.1e9,
//!   "Box": { "Width": 0.01, "Height": 0.005, "CellsX": 80, "CellsY": 40 },
//!   "Ports": [ {"Index": 1, "Attributes": [1001], "Direction": "x"} ]
//! }
//! ```

use serde::Deserialize;
use std::ops::Deref;

/// Top-level MoM configuration — flat JSON with no `Solver.MoM` wrapper.
///
/// All MoM-specific solver parameters live directly at the top level via
/// `#[serde(flatten)]` + `Deref`, so `cfg.equation` works directly.
#[derive(Debug, Clone, Deserialize)]
pub struct MomConfig {
    /// Output directory.
    #[serde(rename = "Output", default)]
    pub output: Option<String>,

    /// Path to the Gmsh mesh file (`.msh`).
    #[serde(rename = "Mesh", default)]
    pub mesh: String,

    /// Reference length scale [m] (default 1.0).
    #[serde(rename = "L0", default = "default_l0")]
    pub l0: f64,

    /// Flattened MoM solver parameters (all `Solver.MoM.*` fields).
    #[serde(flatten)]
    pub solver: super::MomSolverConfig,
}

/// All MoM solver fields are directly accessible: `cfg.equation`, `cfg.box_config`, etc.
impl Deref for MomConfig {
    type Target = super::MomSolverConfig;
    fn deref(&self) -> &Self::Target { &self.solver }
}

fn default_l0() -> f64 { 1.0 }

/// Load a MoM config from a JSON file path.
pub fn load_mom_config(path: &std::path::Path) -> Result<MomConfig, crate::RemError> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| crate::RemError::Config(format!("reading mom config: {e}")))?;
    let cfg: MomConfig = serde_json::from_str(&text)
        .map_err(|e| crate::RemError::Config(format!("parsing mom config: {e}")))?;
    Ok(cfg)
}
