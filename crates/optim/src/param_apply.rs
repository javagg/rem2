//! Apply named parameter values to a cloned [`PalaceConfig`].
//!
//! The `apply` function modifies the relevant field in `MomSolverConfig` /
//! `SubstrateLayerConfig` / port impedance to a new scalar value.

use rem_config::{PalaceConfig, ParamTarget, SweepParam};
use rem_core::RemResult;

/// Apply a single `(param, value)` pair to a cloned config.
pub fn apply_param(config: &mut PalaceConfig, param: &SweepParam, value: f64) -> RemResult<()> {
    apply_target(config, &param.target, value)
}

/// Apply a `ParamTarget` → `value` to a mutable config.
pub fn apply_target(config: &mut PalaceConfig, target: &ParamTarget, value: f64) -> RemResult<()> {
    let mom_cfg = config.solver.mom.as_mut().ok_or_else(|| rem_core::RemError::Config(
        "Parametric: Solver.MoM section required".to_string(),
    ))?;

    match target {
        ParamTarget::SubstratePermittivity { layer } => {
            let layers = &mut mom_cfg
                .substrate.as_mut()
                .ok_or_else(|| rem_core::RemError::Config(
                    "Parametric: SubstratePermittivity requires Solver.MoM.Substrate".to_string()
                ))?
                .layers;
            if *layer >= layers.len() {
                return Err(rem_core::RemError::Config(format!(
                    "Parametric: SubstratePermittivity layer index {} out of range ({})",
                    layer, layers.len()
                )));
            }
            layers[*layer].permittivity = value;
        }

        ParamTarget::SubstrateThickness { layer } => {
            let layers = &mut mom_cfg
                .substrate.as_mut()
                .ok_or_else(|| rem_core::RemError::Config(
                    "Parametric: SubstrateThickness requires Solver.MoM.Substrate".to_string()
                ))?
                .layers;
            if *layer >= layers.len() {
                return Err(rem_core::RemError::Config(format!(
                    "Parametric: SubstrateThickness layer index {} out of range ({})",
                    layer, layers.len()
                )));
            }
            layers[*layer].thickness = value;
        }

        ParamTarget::SubstrateLossTangent { layer } => {
            let layers = &mut mom_cfg
                .substrate.as_mut()
                .ok_or_else(|| rem_core::RemError::Config(
                    "Parametric: SubstrateLossTangent requires Solver.MoM.Substrate".to_string()
                ))?
                .layers;
            if *layer >= layers.len() {
                return Err(rem_core::RemError::Config(format!(
                    "Parametric: SubstrateLossTangent layer index {} out of range ({})",
                    layer, layers.len()
                )));
            }
            layers[*layer].loss_tangent = value;
        }

        ParamTarget::PortZ0 { port } => {
            let p = mom_cfg.ports.iter_mut().find(|p| p.index as usize == *port)
                .ok_or_else(|| rem_core::RemError::Config(format!(
                    "Parametric: PortZ0 port index {} not found in Solver.MoM.Ports", port
                )))?;
            p.impedance = Some(value);
        }

        ParamTarget::FreqMin => { mom_cfg.freq_min = value; }
        ParamTarget::FreqMax => { mom_cfg.freq_max = value; }
    }
    Ok(())
}

/// Apply a full vector of (param, value) pairs to a cloned config.
pub fn apply_params(
    base_config: &PalaceConfig,
    params: &[&SweepParam],
    values: &[f64],
) -> RemResult<PalaceConfig> {
    let mut cfg = base_config.clone();
    for (param, &val) in params.iter().zip(values.iter()) {
        apply_param(&mut cfg, param, val)?;
    }
    Ok(cfg)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rem_config::{ParamTarget, SweepParam};

    fn make_param(target: ParamTarget) -> SweepParam {
        SweepParam {
            name: "p".to_string(),
            target,
            values: vec![],
            min: None, max: None, steps: None,
            initial: None, bounds: None,
        }
    }

    #[test]
    fn apply_freq_min() {
        use rem_config::{load_config_from_str, ConfigFormat};
        let json = r#"{
            "Problem": {"Type": "MoM"},
            "Model":   {"Mesh": "test.msh"},
            "Boundaries": {"PEC": {"Attributes": [1]}},
            "Solver": {"MoM": {
                "FreqMin": 1e9, "FreqMax": 2e9, "FreqStep": 1e9,
                "Ports": []
            }}
        }"#;
        let mut cfg = load_config_from_str(json, ConfigFormat::Json).expect("parse");
        let param = make_param(ParamTarget::FreqMin);
        apply_param(&mut cfg, &param, 3e9).expect("apply");
        assert_eq!(cfg.solver.mom.unwrap().freq_min, 3e9);
    }

    #[test]
    fn apply_substrate_permittivity() {
        use rem_config::{load_config_from_str, ConfigFormat};
        let json = r#"{
            "Problem": {"Type": "MoM"},
            "Model":   {"Mesh": "test.msh"},
            "Boundaries": {"PEC": {"Attributes": [1]}},
            "Solver": {"MoM": {
                "FreqMin": 1e9, "FreqMax": 2e9, "FreqStep": 1e9,
                "Ports": [],
                "Substrate": {
                    "Layers": [{"Permittivity": 4.0, "LossTangent": 0.02,
                                "Thickness": 1e-3}]
                }
            }}
        }"#;
        let mut cfg = load_config_from_str(json, ConfigFormat::Json).expect("parse");
        let param = make_param(ParamTarget::SubstratePermittivity { layer: 0 });
        apply_param(&mut cfg, &param, 6.0).expect("apply");
        let eps = cfg.solver.mom.unwrap().substrate.unwrap().layers[0].permittivity;
        assert!((eps - 6.0).abs() < 1e-9);
    }
}
