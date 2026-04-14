pub mod ads;
pub mod ansys;
pub mod sonnet19;

use std::path::Path;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProjectFormat {
	Sonnet19,
	Ansys,
	Ads,
}

#[derive(Debug, Clone, Copy, Default)]
pub struct Sonnet19Overrides {
	pub freq_min: Option<f64>,
	pub freq_max: Option<f64>,
	pub freq_step: Option<f64>,
}

pub fn convert_project_to_rem(
	format: ProjectFormat,
	project_path: &Path,
	out_config: &Path,
	out_msh: &Path,
	sonnet19_overrides: Sonnet19Overrides,
) -> anyhow::Result<()> {
	match format {
		ProjectFormat::Sonnet19 => sonnet19::convert_xml_to_rem(
			project_path,
			out_config,
			out_msh,
			sonnet19_overrides.freq_min,
			sonnet19_overrides.freq_max,
			sonnet19_overrides.freq_step,
		),
		ProjectFormat::Ansys => ansys::convert_project_to_rem(project_path, out_config, out_msh),
		ProjectFormat::Ads => ads::convert_project_to_rem(project_path, out_config, out_msh),
	}
}
