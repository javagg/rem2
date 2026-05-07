//! Runtime validation for Palace config files.
//!
//! Two phases:
//!   1. Structural: parse to `serde_json::Value` and check required keys / known enum values
//!      before full serde deserialization.
//!   2. Semantic: validate the deserialized `PalaceConfig` for logical consistency.

use rem_core::{RemError, RemResult};

fn is_supported_formulation(v: &str) -> bool {
    matches!(v.to_lowercase().as_str(), "" | "auto" | "h1" | "hcurl" | "nedelec")
}

// ---------------------------------------------------------------------------
// Known problem types
// ---------------------------------------------------------------------------

const KNOWN_PROBLEM_TYPES: &[&str] = &[
    "Electrostatic", "Magnetostatic", "Eigenmode", "Driven", "Transient",
    "MoM", "SBR",
];

// ---------------------------------------------------------------------------
// Phase 1: Structural pre-validation (JSON)
// ---------------------------------------------------------------------------

/// Validate a JSON config string for required top-level keys and known enum values.
/// Returns an error with a clear, human-readable message on the first violation.
pub(crate) fn validate_json_structure(json: &str) -> RemResult<()> {
    let v: serde_json::Value = serde_json::from_str(json)
        .map_err(|e| RemError::Config(format!("JSON parse error: {}", e)))?;

    let root = match v.as_object() {
        Some(o) => o,
        None => return Err(RemError::Config(
            "Config must be a JSON object at the top level".to_string()
        )),
    };

    // Required: "Problem"
    let problem = match root.get("Problem") {
        Some(p) => p,
        None => return Err(RemError::Config(
            "Missing required key \"Problem\" in config. \
             Expected: {\"Problem\": {\"Type\": \"Electrostatic\"}, ...}".to_string()
        )),
    };

    // Required: "Problem.Type"
    let ptype = problem.get("Type")
        .and_then(|t| t.as_str())
        .ok_or_else(|| RemError::Config(
            "Missing required key \"Problem.Type\". \
             Must be one of: Electrostatic, Magnetostatic, Eigenmode, Driven, Transient, MoM, SBR".to_string()
        ))?;

    if !KNOWN_PROBLEM_TYPES.contains(&ptype) {
        return Err(RemError::Config(format!(
            "Unknown Problem.Type = \"{}\". \
             Supported types: {}",
            ptype,
            KNOWN_PROBLEM_TYPES.join(", ")
        )));
    }

    // Required: "Model"
    let model = match root.get("Model") {
        Some(m) => m,
        None => return Err(RemError::Config(
            "Missing required key \"Model\" in config. \
             Expected: {\"Model\": {\"Mesh\": \"path/to/mesh.msh\"}}".to_string()
        )),
    };

    // Required: "Model.Mesh"
    if model.get("Mesh").and_then(|m| m.as_str()).is_none() {
        return Err(RemError::Config(
            "Missing required key \"Model.Mesh\". \
             Must be a path to a .msh or .mesh file.".to_string()
        ));
    }

    // Warn (not error) if "Domains" has no materials
    if let Some(domains) = root.get("Domains") {
        let n_mats = domains.get("Materials")
            .and_then(|m| m.as_array())
            .map(|a| a.len())
            .unwrap_or(0);
        if n_mats == 0 {
            log::info!(
                "[REM] Config validation: Domains.Materials is empty — vacuum (ε₀, μ₀) will be used everywhere."
            );
        }
    }

    Ok(())
}

/// Validate a YAML config string (same checks, via serde_yaml → serde_json::Value).
pub(crate) fn validate_yaml_structure(yaml: &str) -> RemResult<()> {
    let v: serde_json::Value = serde_yaml::from_str(yaml)
        .map_err(|e| RemError::Config(format!("YAML parse error: {}", e)))?;
    let json = serde_json::to_string(&v)
        .map_err(|e| RemError::Config(format!("Internal re-encode error: {}", e)))?;
    validate_json_structure(&json)
}

// ---------------------------------------------------------------------------
// Phase 2: Semantic post-deserialization validation
// ---------------------------------------------------------------------------

/// Validate semantic constraints on the fully-deserialized config.
pub(crate) fn validate_config_semantics(cfg: &super::PalaceConfig) -> RemResult<()> {
    use super::ProblemType;

    // Check for overlapping material attribute assignments
    let mut seen_attrs: std::collections::HashMap<u32, usize> = std::collections::HashMap::new();
    for (idx, mat) in cfg.domains.materials.iter().enumerate() {
        for &attr in &mat.attributes {
            if let Some(prev) = seen_attrs.insert(attr, idx) {
                log::warn!(
                    "[REM] Config: physical group attribute {} is assigned to both material {} and {} — \
                     last definition wins. This may cause unexpected results.",
                    attr, prev, idx
                );
            }
        }
    }

    // Driven: check frequency range
    if cfg.problem.problem_type == ProblemType::Driven {
        if let Some(drv) = &cfg.solver.driven {
            let min_f = drv.min_freq;
            let max_f = drv.max_freq;
            if max_f > 0.0 && min_f > max_f {
                return Err(RemError::Config(format!(
                    "Solver.Driven.MinFreq ({:.3e} Hz) > MaxFreq ({:.3e} Hz). \
                     Frequency sweep requires MinFreq ≤ MaxFreq.",
                    min_f, max_f
                )));
            }
            if min_f < 0.0 || max_f < 0.0 {
                return Err(RemError::Config(
                    "Solver.Driven.MinFreq and MaxFreq must be non-negative.".to_string()
                ));
            }
            if !is_supported_formulation(&drv.formulation) {
                return Err(RemError::Config(format!(
                    "Solver.Driven.Formulation = \"{}\" is not supported. Use Auto, HCurl/Nedelec, or H1.",
                    drv.formulation
                )));
            }
            if let Some(order) = drv.hcurl_order {
                if !(1..=2).contains(&order) {
                    return Err(RemError::Config(format!(
                        "Solver.Driven.HCurlOrder = {} is invalid. Supported values: 1 or 2.",
                        order
                    )));
                }
            }
        } else {
            log::warn!(
                "[REM] Config: Problem.Type = \"Driven\" but no Solver.Driven section found. \
                 MinFreq/MaxFreq/FreqStep must be set before running the driven solver."
            );
        }
    }

    // Eigenmode: n_modes must be positive
    if cfg.problem.problem_type == ProblemType::Eigenmode {
        if let Some(eig) = &cfg.solver.eigenmode {
            if eig.n == 0 {
                return Err(RemError::Config(
                    "Solver.Eigenmode.N must be ≥ 1.".to_string()
                ));
            }
            if !is_supported_formulation(&eig.formulation) {
                return Err(RemError::Config(format!(
                    "Solver.Eigenmode.Formulation = \"{}\" is not supported. Use Auto, HCurl/Nedelec, or H1.",
                    eig.formulation
                )));
            }
            if let Some(order) = eig.hcurl_order {
                if !(1..=2).contains(&order) {
                    return Err(RemError::Config(format!(
                        "Solver.Eigenmode.HCurlOrder = {} is invalid. Supported values: 1 or 2.",
                        order
                    )));
                }
            }
        } else {
            log::warn!(
                "[REM] Config: Problem.Type = \"Eigenmode\" but no Solver.Eigenmode section found. \
                 N (number of modes) must be set before running the eigenmode solver."
            );
        }
    }

    // Ports: check for duplicate port indices
    let mut port_indices: std::collections::HashSet<u32> = std::collections::HashSet::new();
    for port in &cfg.boundaries.lumped_port {
        if !port_indices.insert(port.index) {
            log::warn!(
                "[REM] Config: duplicate LumpedPort index {} — only the last definition will be used.",
                port.index
            );
        }
    }
    let mut wp_indices: std::collections::HashSet<u32> = std::collections::HashSet::new();
    for port in &cfg.boundaries.wave_port {
        if !wp_indices.insert(port.index) {
            log::warn!(
                "[REM] Config: duplicate WavePort index {} — only the last definition will be used.",
                port.index
            );
        }
    }

    log::info!("[REM] Config validation passed.");
    Ok(())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn valid_minimal_json() {
        let json = r#"{"Problem":{"Type":"Electrostatic"},"Model":{"Mesh":"x.msh"}}"#;
        assert!(validate_json_structure(json).is_ok());
    }

    #[test]
    fn missing_problem() {
        let json = r#"{"Model":{"Mesh":"x.msh"}}"#;
        let err = validate_json_structure(json).unwrap_err();
        assert!(err.to_string().contains("Problem"));
    }

    #[test]
    fn missing_model_mesh() {
        let json = r#"{"Problem":{"Type":"Electrostatic"},"Model":{}}"#;
        let err = validate_json_structure(json).unwrap_err();
        assert!(err.to_string().contains("Mesh"));
    }

    #[test]
    fn unknown_problem_type() {
        let json = r#"{"Problem":{"Type":"QuantumField"},"Model":{"Mesh":"x.msh"}}"#;
        let err = validate_json_structure(json).unwrap_err();
        assert!(err.to_string().contains("Unknown Problem.Type"));
    }

    #[test]
    fn driven_frequency_order_error() {
        use crate::{load_config_from_str, ConfigFormat};
        let json = r#"{"Problem":{"Type":"Driven"},"Model":{"Mesh":"x.msh"},
            "Solver":{"Driven":{"MinFreq":10e9,"MaxFreq":1e9,"FreqStep":1e8}}}"#;
        let err = load_config_from_str(json, ConfigFormat::Json).unwrap_err();
        assert!(err.to_string().contains("MinFreq"), "got: {}", err);
    }

    #[test]
    fn driven_formulation_override_h1() {
        use crate::{load_config_from_str, ConfigFormat};
        let json = r#"{
            "Problem": {"Type": "Driven"},
            "Model": {"Mesh": "x.msh"},
            "Solver": {
                "Discretization": "HCurl",
                "Driven": {"MinFreq": 1e9, "MaxFreq": 1e9, "FreqStep": 1e8, "Formulation": "H1"}
            }
        }"#;
        let cfg = load_config_from_str(json, ConfigFormat::Json).expect("config should parse");
        assert!(cfg.solver.uses_hcurl());
        assert!(!cfg.solver.uses_hcurl_for_driven());
    }

    #[test]
    fn eigen_formulation_override_hcurl() {
        use crate::{load_config_from_str, ConfigFormat};
        let json = r#"{
            "Problem": {"Type": "Eigenmode"},
            "Model": {"Mesh": "x.msh"},
            "Solver": {
                "Discretization": "H1",
                "Eigenmode": {"N": 1, "Formulation": "HCurl"}
            }
        }"#;
        let cfg = load_config_from_str(json, ConfigFormat::Json).expect("config should parse");
        assert!(!cfg.solver.uses_hcurl());
        assert!(cfg.solver.uses_hcurl_for_eigenmode());
    }

    #[test]
    fn invalid_driven_formulation_rejected() {
        use crate::{load_config_from_str, ConfigFormat};
        let json = r#"{
            "Problem": {"Type": "Driven"},
            "Model": {"Mesh": "x.msh"},
            "Solver": {
                "Driven": {"MinFreq": 1e9, "MaxFreq": 1e9, "FreqStep": 1e8, "Formulation": "Foo"}
            }
        }"#;
        let err = load_config_from_str(json, ConfigFormat::Json).unwrap_err();
        assert!(err.to_string().contains("Formulation"), "got: {}", err);
    }

    #[test]
    fn driven_hcurl_order_override_and_fallback() {
        use crate::{load_config_from_str, ConfigFormat};
        let json = r#"{
            "Problem": {"Type": "Driven"},
            "Model": {"Mesh": "x.msh"},
            "Solver": {
                "Order": 1,
                "Driven": {
                    "MinFreq": 1e9,
                    "MaxFreq": 1e9,
                    "FreqStep": 1e8,
                    "HCurlOrder": 2
                }
            }
        }"#;
        let cfg = load_config_from_str(json, ConfigFormat::Json).expect("config should parse");
        assert_eq!(cfg.solver.order, 1);
        assert_eq!(cfg.solver.driven_hcurl_order(), 2);
    }

    #[test]
    fn eigen_hcurl_order_falls_back_to_solver_order() {
        use crate::{load_config_from_str, ConfigFormat};
        let json = r#"{
            "Problem": {"Type": "Eigenmode"},
            "Model": {"Mesh": "x.msh"},
            "Solver": {
                "Order": 2,
                "Eigenmode": {"N": 1}
            }
        }"#;
        let cfg = load_config_from_str(json, ConfigFormat::Json).expect("config should parse");
        assert_eq!(cfg.solver.eigenmode_hcurl_order(), 2);
    }

    #[test]
    fn invalid_hcurl_order_rejected() {
        use crate::{load_config_from_str, ConfigFormat};
        let json = r#"{
            "Problem": {"Type": "Driven"},
            "Model": {"Mesh": "x.msh"},
            "Solver": {
                "Driven": {
                    "MinFreq": 1e9,
                    "MaxFreq": 1e9,
                    "FreqStep": 1e8,
                    "HCurlOrder": 3
                }
            }
        }"#;
        let err = load_config_from_str(json, ConfigFormat::Json).unwrap_err();
        assert!(err.to_string().contains("HCurlOrder"), "got: {}", err);
    }
}
