//! Standalone MoM config format for `rem-suite run mom <config.json>`.

use serde::Deserialize;

/// Top-level MoM configuration file.
#[derive(Debug, Clone, Deserialize)]
pub struct MomConfig {
    #[serde(rename = "Problem", default)]
    pub problem: MomProblem,
    #[serde(rename = "Model", default)]
    pub model: MomModel,
    #[serde(rename = "Solver")]
    pub solver: super::MomSolverConfig,
}

#[derive(Debug, Clone, Deserialize, Default)]
pub struct MomProblem {
    #[serde(rename = "Output", default)]
    pub output: Option<String>,
}

#[derive(Debug, Clone, Deserialize, Default)]
pub struct MomModel {
    #[serde(rename = "Mesh", default)]
    pub mesh: String,
    #[serde(rename = "L0", default = "default_l0")]
    pub l0: f64,
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
