use anyhow::bail;
use std::path::Path;

pub fn convert_project_to_rem(
    _project_path: &Path,
    _out_config: &Path,
    _out_msh: &Path,
) -> anyhow::Result<()> {
    bail!("ADS conversion is not implemented yet")
}
