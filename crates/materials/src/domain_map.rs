use crate::material::Material;
use rem_config::PalaceConfig;
use rem_core::RemResult;
use std::collections::HashMap;

/// Maps GMSH physical group tag → Material.
pub struct DomainMap {
    /// Indexed by material index (same order as `PalaceConfig::domains.materials`)
    pub materials: Vec<Material>,

    /// Physical tag → index into `materials`
    tag_to_mat: HashMap<u32, usize>,

    /// Default material (vacuum) used when a tag has no explicit assignment
    default_material: Material,
}

impl DomainMap {
    pub fn from_config(config: &PalaceConfig) -> RemResult<Self> {
        let materials: Vec<Material> = config
            .domains
            .materials
            .iter()
            .map(|spec| Material {
                permittivity: spec.permittivity,
                permeability: spec.permeability,
                conductivity: spec.conductivity,
                loss_tangent: spec.loss_tangent,
            })
            .collect();

        let mut tag_to_mat: HashMap<u32, usize> = HashMap::new();
        for (idx, spec) in config.domains.materials.iter().enumerate() {
            for &tag in &spec.attributes {
                if let Some(prev) = tag_to_mat.insert(tag, idx) {
                    log::warn!(
                        "Physical group {} is assigned to both material {} and {} — using the latter",
                        tag, prev, idx
                    );
                }
            }
        }

        Ok(DomainMap {
            materials,
            tag_to_mat,
            default_material: Material::default(),
        })
    }

    /// Look up material by physical group tag.
    /// Falls back to vacuum if the tag is not mapped.
    pub fn get(&self, tag: u32) -> &Material {
        self.tag_to_mat
            .get(&tag)
            .and_then(|&idx| self.materials.get(idx))
            .unwrap_or(&self.default_material)
    }

    /// Same as `get` but also returns the material index.
    pub fn get_indexed(&self, tag: u32) -> (usize, &Material) {
        if let Some(&idx) = self.tag_to_mat.get(&tag) {
            if let Some(m) = self.materials.get(idx) {
                return (idx, m);
            }
        }
        (usize::MAX, &self.default_material)
    }

    pub fn n_materials(&self) -> usize { self.materials.len() }

    /// Build a per-element epsilon coefficient array (for fem-rs assemblers).
    /// `element_tags` is a slice of physical group tags, one per volume element.
    pub fn epsilon_per_element(&self, element_tags: &[u32]) -> Vec<f64> {
        element_tags
            .iter()
            .map(|&tag| self.get(tag).epsilon_abs())
            .collect()
    }

    /// Build a per-element reluctivity coefficient array.
    pub fn reluctivity_per_element(&self, element_tags: &[u32]) -> Vec<f64> {
        element_tags
            .iter()
            .map(|&tag| self.get(tag).reluctivity())
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rem_config::{load_config_from_str, ConfigFormat};
    use rem_core::constants::EPS0;

    fn make_map(json: &str) -> DomainMap {
        let cfg = load_config_from_str(json, ConfigFormat::Json).unwrap();
        DomainMap::from_config(&cfg).unwrap()
    }

    #[test]
    fn single_material() {
        let map = make_map(r#"{
            "Problem": {"Type": "Electrostatic"},
            "Model":   {"Mesh": "x.msh"},
            "Domains": {
                "Materials": [{"Attributes": [1], "Permittivity": 4.5}]
            }
        }"#);
        assert_eq!(map.n_materials(), 1);
        let m = map.get(1);
        assert!((m.permittivity - 4.5).abs() < 1e-12);
        assert!((m.epsilon_abs() - 4.5 * EPS0).abs() < 1e-30);
    }

    #[test]
    fn unknown_tag_returns_vacuum() {
        let map = make_map(r#"{
            "Problem": {"Type": "Electrostatic"},
            "Model":   {"Mesh": "x.msh"},
            "Domains": {
                "Materials": [{"Attributes": [1], "Permittivity": 4.5}]
            }
        }"#);
        let m = map.get(99); // unknown tag
        assert!((m.permittivity - 1.0).abs() < 1e-15);
    }

    #[test]
    fn two_materials() {
        let map = make_map(r#"{
            "Problem": {"Type": "Magnetostatic"},
            "Model":   {"Mesh": "x.msh"},
            "Domains": {
                "Materials": [
                    {"Attributes": [1],   "Permittivity": 1.0, "Permeability": 1.0},
                    {"Attributes": [2,3], "Permittivity": 1.0, "Permeability": 1000.0}
                ]
            }
        }"#);
        assert!((map.get(2).permeability - 1000.0).abs() < 1e-9);
        assert!((map.get(3).permeability - 1000.0).abs() < 1e-9);
    }

    #[test]
    fn epsilon_per_element() {
        let map = make_map(r#"{
            "Problem": {"Type": "Electrostatic"},
            "Model":   {"Mesh": "x.msh"},
            "Domains": {
                "Materials": [{"Attributes": [1], "Permittivity": 2.0}]
            }
        }"#);
        let tags = vec![1u32, 1, 1];
        let eps = map.epsilon_per_element(&tags);
        assert_eq!(eps.len(), 3);
        for &e in &eps {
            assert!((e - 2.0 * EPS0).abs() < 1e-30);
        }
    }
}
