mod material;
mod domain_map;

pub use material::Material;
pub use domain_map::DomainMap;

use rem_config::PalaceConfig;
use rem_core::RemResult;

/// Build a material lookup from a Palace config.
pub fn build_domain_map(config: &PalaceConfig) -> RemResult<DomainMap> {
    DomainMap::from_config(config)
}
