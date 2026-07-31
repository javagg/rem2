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
use crate::schema::MaterialSpec;

/// Top-level MoM configuration — flat JSON with no `Solver.MoM` wrapper.
///
/// All MoM-specific solver parameters live directly at the top level via
/// `#[serde(flatten)]` + `Deref`, so `cfg.equation` works directly.
#[derive(Debug, Clone, Deserialize)]
pub struct MomConfig {
    /// Output directory.
    #[serde(rename = "Output", default = "default_output_dir")]
    pub output_dir: String,

    /// Path to the Gmsh mesh file (`.msh`).
    #[serde(rename = "Mesh", default)]
    pub mesh: String,

    /// Reference length scale [m] (default 1.0).
    #[serde(rename = "L0", default = "default_l0")]
    pub l0: f64,

    /// PEC physical group attributes (surface metal).
    #[serde(rename = "PEC", default)]
    pub pec: Vec<u32>,

    /// Material definitions for volume domains.
    #[serde(rename = "Materials", default)]
    pub materials: Vec<MaterialSpec>,

    /// Conductor terminals for capacitance extraction.
    /// Each entry: `{ "Index": 1, "Attributes": [1001] }`.
    #[serde(rename = "Terminals", default)]
    pub terminals: Vec<TerminalSpec>,

    /// Ground (reference) conductor attributes for capacitance extraction.
    #[serde(rename = "Ground", default)]
    pub ground: Vec<u32>,

    /// Adaptive mesh refinement config.
    #[serde(rename = "Refinement", default)]
    pub refinement: RefinementConfig,

    /// Converter/extractor metadata (PlanarTechLayers, ViaTechLayers, thick metal, etc.).
    #[serde(rename = "Metadata", default)]
    pub metadata: Option<serde_json::Value>,

    /// Flattened MoM solver parameters (all `Solver.MoM.*` fields).
    #[serde(flatten)]
    pub solver: super::MomSolverConfig,
}

/// A single conductor terminal for capacitance extraction.
#[derive(Debug, Clone, Deserialize)]
pub struct TerminalSpec {
    /// Terminal index (0 = ground reference).
    #[serde(rename = "Index")]
    pub index: u32,
    /// Mesh physical-group attributes that belong to this terminal.
    #[serde(rename = "Attributes", deserialize_with = "crate::schema::deserialize_attributes")]
    pub attributes: Vec<u32>,
}

/// Adaptive mesh refinement settings for capacitance extraction.
#[derive(Debug, Clone, Deserialize, Default)]
pub struct RefinementConfig {
    /// Maximum AMR iterations (0 = no refinement).
    #[serde(rename = "MaxIter", default)]
    pub max_iter: usize,
    /// Relative error tolerance for refinement (0 = disable).
    #[serde(rename = "Tol", default = "default_ref_tol")]
    pub tol: f64,
}

fn default_ref_tol() -> f64 { 0.5 }

/// All MoM solver fields are directly accessible: `cfg.equation`, `cfg.box_config`, etc.
impl Deref for MomConfig {
    type Target = super::MomSolverConfig;
    fn deref(&self) -> &Self::Target { &self.solver }
}

fn default_l0() -> f64 { 1.0 }
fn default_output_dir() -> String { "output".to_string() }

/// Load a MoM config from a JSON file path.
pub fn load_mom_config(path: &std::path::Path) -> Result<MomConfig, crate::RemError> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| crate::RemError::Config(format!("reading mom config: {e}")))?;
    let mut cfg: MomConfig = serde_json::from_str(&text)
        .map_err(|e| crate::RemError::Config(format!("parsing mom config: {e}")))?;
    cfg.solver.apply_preset();
    cfg.solver.resolve_calibration();
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
                self.basis = "Rooftop".into();
                self.singular_tol = 1e-12;
                self.fast_solver = "UFFT".into();
                self.kernel = "Cavity".into();
                self.mesh_format = "RectGrid".into();
                self.adaptive_sweep = true;
                self.adaptive_target = 100;
            }
            "ads" => {
                // ADS Momentum: 2.5D planar MoM with conformal triangular mesh,
                // RWG basis, layered media Green's function, open boundaries.
                // ACA + GMRES is preferred over FMM/MLFMA (which target 3D free-space).
                self.basis = "RWG".into();
                self.singular_tol = 1e-12;
                self.fast_solver = "ACA".into();
                self.kernel = "Layered".into();
                self.mesh_format = "TriSurface".into();
                self.adaptive_sweep = false;
                self.adaptive_target = 0;
            }
            "q3d" | "q3d_extractor" => {
                // Q3D Extractor: electrostatic BEM, capacitance matrix.
                // Uses FMM-accelerated GMRES (not Direct LU) for large problems.
                self.basis = "Pulse".into();
                self.singular_tol = 1e-12;
                self.fast_solver = "FMM".into();
                self.kernel = "Laplace".into();
                self.mesh_format = "TriSurface".into();
                self.adaptive_sweep = false;
                self.adaptive_target = 0;
            }
            _ => {} // "Custom" or unknown: keep explicit field values
        }
    }

    /// Resolve per-port calibration — currently disabled (no-op).
    ///
    /// CalPort metadata is preserved but does NOT auto-populate global
    /// TRL/SOLT kits.  Calibration requires fixture-simulated standard
    /// responses which are not available in the current architecture.
    /// See `apply_per_port_calibration()` in calib.rs for details.
    pub fn resolve_calibration(&mut self) {
        let n_cal = self.ports.iter().filter(|p| p.calibration.is_some()).count();
        if n_cal > 0 {
            log::warn!(
                "resolve_calibration: {} CalPort(s) found but calibration is \
                 disabled — no fixture-simulated standard responses available. \
                 CalPort metadata is preserved; raw S-parameters unchanged.",
                n_cal
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_cfg(preset: &str, extra: &str) -> MomConfig {
        let json = format!(
            r#"{{"Mesh":"test.msh","Preset":"{}",{}"FreqMin":1e9,"FreqMax":10e9,"FreqStep":0.1e9}}"#,
            preset, extra
        );
        let mut cfg: MomConfig = serde_json::from_str(&json).unwrap();
        cfg.solver.apply_preset();
        cfg
    }

    fn assert_eq_eps(a: f64, b: f64) {
        assert!((a - b).abs() < 1e-14, "expected {b}, got {a}");
    }

    #[test]
    fn preset_default_is_sonnet19() {
        let cfg = make_cfg("Sonnet19", "");
        assert_eq!(cfg.equation, "EFIE");
        assert_eq!(cfg.basis, "Rooftop");
        assert_eq_eps(cfg.alpha, 1.0);
        assert_eq_eps(cfg.singular_tol, 1e-12);
        assert_eq!(cfg.fast_solver, "UFFT");
        assert_eq!(cfg.kernel, "Cavity");
        assert_eq!(cfg.mesh_format, "RectGrid");
        assert!(cfg.adaptive_sweep);
        assert_eq!(cfg.adaptive_target, 100);
    }

    #[test]
    fn preset_sonnet19_no_preset_field() {
        // When Preset is not set, it defaults to "Sonnet19"
        let json = r#"{"Mesh":"test.msh","FreqMin":1e9,"FreqMax":10e9,"FreqStep":0.1e9}"#;
        let mut cfg: MomConfig = serde_json::from_str(json).unwrap();
        assert_eq!(cfg.preset, "Sonnet19");
        cfg.solver.apply_preset();
        assert_eq!(cfg.equation, "EFIE");
        assert_eq!(cfg.basis, "Rooftop");
    }

    #[test]
    fn preset_ads_sets_open_defaults() {
        let cfg = make_cfg("ADS", "");
        assert_eq!(cfg.equation, "EFIE");
        assert_eq!(cfg.basis, "RWG");
        assert_eq!(cfg.fast_solver, "ACA");
        assert_eq!(cfg.kernel, "Layered");
        assert_eq!(cfg.mesh_format, "TriSurface");
        assert!(!cfg.adaptive_sweep);
        assert_eq!(cfg.adaptive_target, 0);
    }

    #[test]
    fn preset_q3d_sets_electrostatic_defaults() {
        let cfg = make_cfg("Q3D", "");
        assert_eq!(cfg.equation, "EFIE");
        assert_eq!(cfg.basis, "Pulse");
        assert_eq_eps(cfg.alpha, 1.0);
        assert_eq!(cfg.fast_solver, "FMM");
        assert_eq!(cfg.kernel, "Laplace");
        assert_eq!(cfg.mesh_format, "TriSurface");
        assert!(!cfg.adaptive_sweep);
        assert_eq!(cfg.adaptive_target, 0);
    }

    #[test]
    fn preset_custom_preserves_explicit_values() {
        let json = r#"{
            "Mesh":"test.msh",
            "Preset":"Custom",
            "FreqMin":1e9,"FreqMax":10e9,"FreqStep":0.1e9,
            "Equation":"CFIE","Basis":"RWG","Alpha":0.5,
            "FastSolver":"GMRES","AdaptiveSweep":false
        }"#;
        let mut cfg: MomConfig = serde_json::from_str(json).unwrap();
        cfg.solver.apply_preset();
        assert_eq!(cfg.equation, "CFIE");
        assert_eq!(cfg.basis, "RWG");
        assert_eq_eps(cfg.alpha, 0.5);
        assert_eq!(cfg.fast_solver, "GMRES");
        assert!(!cfg.adaptive_sweep);
    }

    #[test]
    fn preset_unknown_treated_as_custom() {
        let cfg = make_cfg("SomeRandomName", "");
        assert_eq!(cfg.equation, "EFIE"); // serde default
        assert_eq!(cfg.basis, "Rooftop"); // serde default
    }
}

    #[test]
    fn mom_config_serialize_roundtrip_keeps_box_and_bricks() {
        // Round-trip a full MomSolverConfig (with Box) through JSON: this is
        // the contract the Sonnet19 converter relies on when it serializes
        // sonnet19_build_rem_config's output and the CLI loads it back.
        let json = r#"{
            "Mesh":"test.msh","Preset":"Custom","Equation":"EFIE","Basis":"Rooftop",
            "FreqMin":1e8,"FreqMax":2e9,"FreqStep":1e7,"FastSolver":"UFFT",
            "Box":{
                "Width":6.4e-4,"Height":6.4e-4,"CellsX":64,"CellsY":64,
                "TopCover":true,"BottomCover":true,"Symmetry":"X",
                "DielectricBricks":[{
                    "IxStart":29,"IxEnd":35,"IyStart":28,"IyEnd":36,
                    "NzLayers":1,"Permittivity":7.5,"Thickness":2e-7
                }]
            }
        }"#;
        let cfg: MomConfig = serde_json::from_str(json).unwrap();
        let box_cfg = cfg.box_config.as_ref().unwrap();
        assert_eq!(box_cfg.dielectric_bricks.len(), 1);
        assert!((box_cfg.dielectric_bricks[0].eps_r - 7.5).abs() < 1e-9);
        assert!((box_cfg.dielectric_bricks[0].thickness - 2e-7).abs() < 1e-12);
        assert_eq!(box_cfg.symmetry, Some(crate::SymmetryAxis::X));

        // Serialize and re-deserialize.
        let text = serde_json::to_string(&cfg.solver).unwrap();
        let back: crate::MomSolverConfig = serde_json::from_str(&text).unwrap();
        let b2 = back.box_config.unwrap();
        assert_eq!(b2.dielectric_bricks.len(), 1);
        assert!((b2.dielectric_bricks[0].eps_r - 7.5).abs() < 1e-9);
        assert!((b2.dielectric_bricks[0].thickness - 2e-7).abs() < 1e-12);
        assert_eq!(b2.symmetry, Some(crate::SymmetryAxis::X));
    }
