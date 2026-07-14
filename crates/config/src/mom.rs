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

    /// Resolve per-port calibration into solver-level TRL/SOLT kit.
    ///
    /// When ports carry `Calibration` data from a Sonnet Cal Port but no
    /// explicit `TrlKit`/`SoltKit` is configured, auto-populate the
    /// solver-level kit from the first CalPort's parameters.  This bridges
    /// the per-port CalPort model into the existing global calibration path.
    ///
    /// Only processes 2-port S-parameter solves (TRL/SOLT require exactly
    /// 2 ports).  When more than one CalPort with conflicting types exists,
    /// the first one wins with a warning.
    pub fn resolve_calibration(&mut self) {
        use super::schema::{CalibrationConfig, TrlKitConfig, SoltKitConfig};

        // Scan ports for Calibration data
        let cal: Option<&CalibrationConfig> = self.ports.iter()
            .find_map(|p| p.calibration.as_ref());

        let Some(cal) = cal else { return };
        if self.trl_kit.is_some() || self.solt_kit.is_some() {
            return; // explicit kit already set
        }

        match cal.cal_type.to_ascii_uppercase().as_str() {
            "TRL" => {
                let line_len = cal.line_length;
                let thru_len = cal.thru_length.max(0.0);
                if line_len <= thru_len {
                    log::warn!("CalPort TRL: LineLength ({:.3e}) <= ThruLength ({:.3e}), skipping",
                        line_len, thru_len);
                    return;
                }
                let delta_l = line_len - thru_len;
                // Estimate ε_eff from line impedance (assume air-line for Z0=50,
                // otherwise back-compute from Z0 characteristic)
                let eps_eff = if (cal.line_impedance - 50.0).abs() < 1.0 {
                    1.0 // air-like
                } else {
                    // Zo = 60/√ε · ln(4h/d) — rough estimate: ε ≈ (60/Zo)²
                    (60.0 / cal.line_impedance.max(1.0)).powi(2)
                };
                self.trl_kit = Some(TrlKitConfig {
                    thru_length: thru_len,
                    line_length: line_len,
                    line_impedance: cal.line_impedance,
                    epsilon_eff: eps_eff,
                    reflect_type: if cal.reflect_type.is_empty() {
                        "SHORT".to_string()
                    } else { cal.reflect_type.clone() },
                    reflect_magnitude: 1.0,
                    solve_side: false,
                });
                log::info!("CalPort: auto-populated TRL kit (Thru={:.3e}, Line={:.3e}, Δl={:.3e})",
                    thru_len, line_len, delta_l);
            }
            "SOLT" => {
                self.solt_kit = Some(SoltKitConfig {
                    short_inductance: 0.0,
                    open_capacitance: 0.0,
                    load_resistance: cal.line_impedance,
                    load_inductance: 0.0,
                    load_capacitance: 0.0,
                    ref_impedance: cal.line_impedance,
                    epsilon_eff: 1.0,
                });
                log::info!("CalPort: auto-populated SOLT kit (Z0={:.1})", cal.line_impedance);
            }
            "LRM" => {
                // LRM maps to TRL with zero-length THRU (ideal match is
                // equivalent to a zero-length thru for the reference plane).
                let line_len = cal.line_length;
                if line_len <= 0.0 {
                    log::warn!("CalPort LRM: LineLength must be positive, skipping");
                    return;
                }
                let eps_eff = if (cal.line_impedance - 50.0).abs() < 1.0 {
                    1.0
                } else {
                    (60.0 / cal.line_impedance.max(1.0)).powi(2)
                };
                self.trl_kit = Some(TrlKitConfig {
                    thru_length: 0.0,
                    line_length: line_len,
                    line_impedance: cal.line_impedance,
                    epsilon_eff: eps_eff,
                    reflect_type: if cal.reflect_type.is_empty() {
                        "SHORT".to_string()
                    } else { cal.reflect_type.clone() },
                    reflect_magnitude: 1.0,
                    solve_side: false,
                });
                log::info!("CalPort: auto-populated LRM kit (Line={:.3e}, mapped to TRL)",
                    line_len);
            }
            _ => {} // NONE / unknown: not supported
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
        assert_eq!(cfg.mom_type, "Boxed");
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
    fn preset_ads_matches_sonnet19() {
        let cfg = make_cfg("ADS", "");
        assert_eq!(cfg.equation, "EFIE");
        assert_eq!(cfg.basis, "Rooftop");
        assert_eq!(cfg.fast_solver, "UFFT");
        assert_eq!(cfg.mom_type, "Boxed");
        assert_eq!(cfg.kernel, "Cavity");
        assert_eq!(cfg.mesh_format, "RectGrid");
        assert!(cfg.adaptive_sweep);
    }

    #[test]
    fn preset_q3d_sets_electrostatic_defaults() {
        let cfg = make_cfg("Q3D", "");
        assert_eq!(cfg.equation, "EFIE");
        assert_eq!(cfg.basis, "Pulse");
        assert_eq_eps(cfg.alpha, 1.0);
        assert_eq!(cfg.fast_solver, "Direct");
        assert_eq!(cfg.mom_type, "Capacitance");
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
