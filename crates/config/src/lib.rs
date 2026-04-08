pub mod schema;
mod preprocess;

pub use schema::*;
use rem_core::{RemError, RemResult};
use std::path::Path;

fn json_err(e: serde_json::Error) -> RemError { RemError::Config(e.to_string()) }
fn yaml_err(e: serde_yaml::Error) -> RemError { RemError::Config(e.to_string()) }

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConfigFormat {
    Json,
    Yaml,
}

/// Load a Palace-compatible config file (JSON or YAML).
/// JSON supports C++-style `//` and `/* */` comments and attribute ranges like `"1,3-5"`.
pub fn load_config(path: &Path) -> RemResult<PalaceConfig> {
    let content = std::fs::read_to_string(path).map_err(RemError::Io)?;
    let fmt = match path.extension().and_then(|e| e.to_str()) {
        Some("json") => ConfigFormat::Json,
        Some("yaml") | Some("yml") => ConfigFormat::Yaml,
        _ => return Err(RemError::UnknownFormat),
    };
    load_config_from_str(&content, fmt)
}

/// Parse a Palace config from a string with the given format.
pub fn load_config_from_str(content: &str, fmt: ConfigFormat) -> RemResult<PalaceConfig> {
    let cfg = match fmt {
        ConfigFormat::Json => {
            let stripped = preprocess::strip_comments(content);
            serde_json::from_str(&stripped).map_err(json_err)
        }
        ConfigFormat::Yaml => {
            serde_yaml::from_str(content).map_err(yaml_err)
        }
    }?;
    // Warn on Palace fields that are accepted but not fully implemented.
    schema::validate_palace_compat(&cfg);
    Ok(cfg)
}
