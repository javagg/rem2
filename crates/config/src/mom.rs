//! Standalone MoM config format for `rem-suite run mom <config.json>`.
//!
//! JSON is flat — no `Solver.MoM` nesting.  Example:
//!
//! ```json
//! {
//!   "Output": "./output",
//!   "Mesh": "mesh.msh",
//!   "Preset": "Sonnet19",
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
    let mut cfg: MomConfig = serde_json::from_str(&text)
        .map_err(|e| crate::RemError::Config(format!("parsing mom config: {e}")))?;
    cfg.solver.apply_preset();
    Ok(cfg)
}

/// Apply the named preset to override all solver defaults.
///
/// When the user picks a preset, it enforces a coherent set of modeling
/// parameters — any explicit values the user set in the JSON will be
/// overwritten.  Use `Preset = "Custom"` for per-field control.
impl super::MomSolverConfig {
    pub fn apply_preset(&mut self) {
        match self.preset.to_ascii_lowercase().as_str() {
            "sonnet19" | "sonnet" => {
                self.equation = "EFIE".into();
                self.basis = "Rooftop".into();
                self.alpha = 1.0;
                self.singular_tol = 1e-12;
                self.fast_solver = "UFFT".into();
                self.mom_type = "Boxed".into();
                self.kernel = "Cavity".into();
                self.mesh_format = "RectGrid".into();
                self.adaptive_sweep = true;
                self.adaptive_target = 100;
            }
            "ads" | "ads_mom" | "momentum" => {
                // ADS Momentum: planar MoM, same core as Sonnet19.
                // Reserve for future ADS-specific port/calibration defaults.
                self.equation = "EFIE".into();
                self.basis = "Rooftop".into();
                self.alpha = 1.0;
                self.singular_tol = 1e-12;
                self.fast_solver = "UFFT".into();
                self.mom_type = "Boxed".into();
                self.kernel = "Cavity".into();
                self.mesh_format = "RectGrid".into();
                self.adaptive_sweep = true;
                self.adaptive_target = 100;
            }
            "q3d" | "q3d_extractor" => {
                // Q3D Extractor: electrostatic BEM, capacitance matrix.
                self.equation = "EFIE".into();
                self.basis = "Pulse".into();
                self.alpha = 1.0;
                self.singular_tol = 1e-12;
                self.fast_solver = "Direct".into();
                self.mom_type = "Capacitance".into();
                self.kernel = "Laplace".into();
                self.mesh_format = "TriSurface".into();
                self.adaptive_sweep = false;
                self.adaptive_target = 0;
            }
            _ => {} // "Custom" or unknown: keep explicit field values
        }
    }
}
