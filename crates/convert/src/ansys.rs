use anyhow::{bail, Context};
use serde_json::{json, Value};
use std::path::{Path, PathBuf};

pub fn convert_project_to_rem(
    project_path: &Path,
    out_config: &Path,
    out_msh: &Path,
) -> anyhow::Result<()> {
    let meta_path = resolve_metadata_path(project_path)?;
    let meta_text = std::fs::read_to_string(&meta_path)
        .with_context(|| format!("reading Ansys metadata: {}", meta_path.display()))?;
    let meta: Value = serde_json::from_str(&meta_text)
        .with_context(|| format!("parsing Ansys metadata JSON: {}", meta_path.display()))?;

    if let Some(parent) = out_msh.parent() {
        std::fs::create_dir_all(parent)
            .with_context(|| format!("creating output directory: {}", parent.display()))?;
    }
    if let Some(parent) = out_config.parent() {
        std::fs::create_dir_all(parent)
            .with_context(|| format!("creating output directory: {}", parent.display()))?;
    }

    ensure_mesh_available(&meta, &meta_path, out_msh)?;

    let freq_min = read_f64_with_default(&meta, &["solver", "freq_min"], 1.0e9);
    let freq_max = read_f64_with_default(&meta, &["solver", "freq_max"], freq_min);
    let freq_step = read_f64_with_default(&meta, &["solver", "freq_step"], 1.0e8);
    if !(freq_min > 0.0 && freq_max >= freq_min && freq_step > 0.0) {
        bail!(
            "invalid frequency sweep in Ansys metadata; require freq_min > 0, freq_max >= freq_min, freq_step > 0"
        );
    }

    let equation = read_str_with_default(&meta, &["solver", "equation"], "CFIE");
    let basis = read_str_with_default(&meta, &["solver", "basis"], "RWG");
    let alpha = read_f64_with_default(&meta, &["solver", "alpha"], 0.5);
    let fast_solver = read_str_with_default(&meta, &["solver", "fast_solver"], "ACA");
    let ref_impedance = read_f64_with_default(&meta, &["solver", "ref_impedance"], 50.0);

    let problem_output = read_str_with_default(&meta, &["problem", "output"], "./output/ansys");
    let verbose = read_u64_with_default(&meta, &["problem", "verbose"], 1);
    let l0 = read_f64_with_default(&meta, &["model", "l0"], 1.0);

    let ports = build_ports_json(&meta)?;
    let boundaries = build_boundaries_json(&meta);
    let domains = build_domains_json(&meta);
    let substrate = build_substrate_json(&meta);

    let mut mom = serde_json::Map::new();
    mom.insert("Equation".to_string(), json!(equation));
    mom.insert("Basis".to_string(), json!(basis));
    mom.insert("FreqMin".to_string(), json!(freq_min));
    mom.insert("FreqMax".to_string(), json!(freq_max));
    mom.insert("FreqStep".to_string(), json!(freq_step));
    mom.insert("Alpha".to_string(), json!(alpha));
    mom.insert("FastSolver".to_string(), json!(fast_solver));
    mom.insert("RefImpedance".to_string(), json!(ref_impedance));
    if !ports.is_empty() {
        mom.insert("Ports".to_string(), Value::Array(ports));
    }
    if let Some(substrate_v) = substrate {
        mom.insert("Substrate".to_string(), substrate_v);
    }

    let mut cfg = serde_json::Map::new();
    cfg.insert(
        "Problem".to_string(),
        json!({
            "Type": "MoM",
            "Verbose": verbose,
            "Output": problem_output
        }),
    );
    cfg.insert(
        "Model".to_string(),
        json!({
            "Mesh": out_msh.to_string_lossy(),
            "L0": l0
        }),
    );
    if let Some(domains_v) = domains {
        cfg.insert("Domains".to_string(), domains_v);
    }
    cfg.insert("Boundaries".to_string(), boundaries);
    cfg.insert(
        "Solver".to_string(),
        json!({
            "MoM": Value::Object(mom)
        }),
    );

    let cfg = Value::Object(cfg);

    let pretty = serde_json::to_string_pretty(&cfg)?;
    std::fs::write(out_config, pretty)
        .with_context(|| format!("writing config: {}", out_config.display()))?;

    log::info!(
        "Ansys project converted with external mesh: project={}, meta={}, config={}, mesh={}",
        project_path.display(),
        meta_path.display(),
        out_config.display(),
        out_msh.display()
    );

    Ok(())
}

fn resolve_metadata_path(project_path: &Path) -> anyhow::Result<PathBuf> {
    if project_path
        .extension()
        .and_then(|s| s.to_str())
        .is_some_and(|ext| ext.eq_ignore_ascii_case("json"))
        && project_path.exists()
    {
        return Ok(project_path.to_path_buf());
    }

    let mut candidates = Vec::new();
    candidates.push(project_path.with_extension("rem_ansys.json"));
    candidates.push(project_path.with_extension("json"));
    if let Some(parent) = project_path.parent() {
        candidates.push(parent.join("rem_ansys.json"));
    }

    for p in candidates {
        if p.exists() {
            return Ok(p);
        }
    }

    bail!(
        "Ansys metadata JSON not found for {}. Provide one of: <project>.rem_ansys.json, <project>.json, or sibling rem_ansys.json",
        project_path.display()
    )
}

fn ensure_mesh_available(meta: &Value, meta_path: &Path, out_msh: &Path) -> anyhow::Result<()> {
    if out_msh.exists() {
        return Ok(());
    }

    let mesh_path_str = read_str(meta, &["mesh"])
        .or_else(|| read_str(meta, &["model", "mesh"]))
        .ok_or_else(|| {
            anyhow::anyhow!(
                "output mesh {} does not exist and metadata has no 'mesh'/'model.mesh' path",
                out_msh.display()
            )
        })?;

    let source_path = resolve_relative_to(meta_path, Path::new(&mesh_path_str));
    if !source_path.exists() {
        bail!(
            "mesh source path from metadata does not exist: {}",
            source_path.display()
        );
    }

    std::fs::copy(&source_path, out_msh).with_context(|| {
        format!(
            "copying mesh from {} to {}",
            source_path.display(),
            out_msh.display()
        )
    })?;
    Ok(())
}

fn resolve_relative_to(base_file: &Path, path: &Path) -> PathBuf {
    if path.is_absolute() {
        return path.to_path_buf();
    }
    base_file
        .parent()
        .unwrap_or_else(|| Path::new("."))
        .join(path)
}

fn read_path<'a>(root: &'a Value, path: &[&str]) -> Option<&'a Value> {
    let mut cur = root;
    for key in path {
        cur = cur.get(*key)?;
    }
    Some(cur)
}

fn read_str(root: &Value, path: &[&str]) -> Option<String> {
    read_path(root, path)
        .and_then(Value::as_str)
        .map(ToOwned::to_owned)
}

fn read_str_with_default(root: &Value, path: &[&str], default: &str) -> String {
    read_str(root, path).unwrap_or_else(|| default.to_string())
}

fn read_f64_with_default(root: &Value, path: &[&str], default: f64) -> f64 {
    read_path(root, path)
        .and_then(Value::as_f64)
        .unwrap_or(default)
}

fn read_u64_with_default(root: &Value, path: &[&str], default: u64) -> u64 {
    read_path(root, path)
        .and_then(Value::as_u64)
        .unwrap_or(default)
}

fn read_u32_vec(root: &Value, path: &[&str]) -> Vec<u32> {
    read_path(root, path)
        .and_then(Value::as_array)
        .map(|arr| {
            arr.iter()
                .filter_map(|v| v.as_u64())
                .map(|v| v as u32)
                .collect::<Vec<u32>>()
        })
        .unwrap_or_default()
}

fn build_ports_json(meta: &Value) -> anyhow::Result<Vec<Value>> {
    let Some(ports) = read_path(meta, &["ports"]).and_then(Value::as_array) else {
        return Ok(Vec::new());
    };

    let mut out = Vec::with_capacity(ports.len());
    for (idx, p) in ports.iter().enumerate() {
        let index = p.get("index").and_then(Value::as_u64).unwrap_or((idx + 1) as u64) as u32;
        let attributes = p
            .get("attributes")
            .and_then(Value::as_array)
            .map(|arr| {
                arr.iter()
                    .filter_map(|v| v.as_u64())
                    .map(|v| v as u32)
                    .collect::<Vec<u32>>()
            })
            .unwrap_or_default();
        if attributes.is_empty() {
            bail!("port #{} must include non-empty 'attributes' array", index);
        }
        let direction = p
            .get("direction")
            .and_then(Value::as_str)
            .unwrap_or("x");
        let impedance = p
            .get("impedance")
            .and_then(Value::as_f64)
            .unwrap_or(50.0);

        out.push(json!({
            "Index": index,
            "Attributes": attributes,
            "Direction": direction,
            "Impedance": impedance
        }));
    }

    Ok(out)
}

fn build_boundaries_json(meta: &Value) -> Value {
    let pec = read_u32_vec(meta, &["boundaries", "pec"]);
    let pmc = read_u32_vec(meta, &["boundaries", "pmc"]);
    let radiation = read_u32_vec(meta, &["boundaries", "radiation"]);

    let mut map = serde_json::Map::new();
    if !pec.is_empty() {
        map.insert("PEC".to_string(), json!({ "Attributes": pec }));
    }
    if !pmc.is_empty() {
        map.insert("PMC".to_string(), json!({ "Attributes": pmc }));
    }
    if !radiation.is_empty() {
        map.insert("Radiation".to_string(), json!({ "Attributes": radiation }));
    }

    if map.is_empty() {
        json!({ "PEC": { "Attributes": [1] } })
    } else {
        Value::Object(map)
    }
}

fn build_domains_json(meta: &Value) -> Option<Value> {
    let Some(materials) = read_path(meta, &["materials"]).and_then(Value::as_array) else {
        return None;
    };

    let mut out = Vec::new();
    for m in materials {
        let attrs = m
            .get("attributes")
            .and_then(Value::as_array)
            .map(|arr| {
                arr.iter()
                    .filter_map(|v| v.as_u64())
                    .map(|v| v as u32)
                    .collect::<Vec<u32>>()
            })
            .unwrap_or_default();
        if attrs.is_empty() {
            continue;
        }

        let mut mm = serde_json::Map::new();
        mm.insert("Attributes".to_string(), json!(attrs));
        if let Some(v) = m.get("permittivity") {
            mm.insert("Permittivity".to_string(), v.clone());
        }
        if let Some(v) = m.get("permeability") {
            mm.insert("Permeability".to_string(), v.clone());
        }
        if let Some(v) = m.get("loss_tan") {
            mm.insert("LossTan".to_string(), v.clone());
        }
        out.push(Value::Object(mm));
    }

    if out.is_empty() {
        None
    } else {
        Some(json!({ "Materials": out }))
    }
}

fn build_substrate_json(meta: &Value) -> Option<Value> {
    let Some(layers) = read_path(meta, &["substrate", "layers"]).and_then(Value::as_array) else {
        return None;
    };

    let mut out_layers = Vec::new();
    for layer in layers {
        let permittivity = layer
            .get("permittivity")
            .and_then(Value::as_f64)
            .unwrap_or(1.0);
        let loss_tangent = layer
            .get("loss_tangent")
            .and_then(Value::as_f64)
            .unwrap_or(0.0);
        let permeability = layer
            .get("permeability")
            .and_then(Value::as_f64)
            .unwrap_or(1.0);
        let thickness = layer
            .get("thickness")
            .and_then(Value::as_f64)
            .unwrap_or(0.0);
        let name = layer
            .get("name")
            .and_then(Value::as_str)
            .unwrap_or("layer");

        out_layers.push(json!({
            "Permittivity": permittivity,
            "LossTangent": loss_tangent,
            "Permeability": permeability,
            "Thickness": thickness,
            "Name": name
        }));
    }

    if out_layers.is_empty() {
        None
    } else {
        Some(json!({
            "Layers": out_layers,
            "BottomPec": true
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn unique_temp_dir() -> PathBuf {
        let ts = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock should be monotonic")
            .as_nanos();
        let p = std::env::temp_dir().join(format!("rem_ansys_convert_test_{}", ts));
        std::fs::create_dir_all(&p).expect("temp dir should be created");
        p
    }

    #[test]
    fn converts_with_sidecar_metadata_and_external_mesh() {
        let dir = unique_temp_dir();
        let project = dir.join("demo.aedt");
        let sidecar = dir.join("demo.rem_ansys.json");
        let source_mesh = dir.join("rmsh_out.msh");
        let out_mesh = dir.join("out").join("demo.msh");
        let out_cfg = dir.join("out").join("demo.json");

        std::fs::write(&project, "dummy").expect("project stub should be writable");
        std::fs::write(&source_mesh, "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
            .expect("mesh should be writable");

        let meta = json!({
            "mesh": "rmsh_out.msh",
            "problem": { "output": "./output/ansys", "verbose": 1 },
            "solver": {
                "freq_min": 1.0e9,
                "freq_max": 1.5e9,
                "freq_step": 1.0e8,
                "equation": "CFIE",
                "basis": "RWG"
            },
            "boundaries": {
                "pec": [1, 2]
            },
            "materials": [
                {
                    "attributes": [10],
                    "permittivity": 2.2,
                    "permeability": 1.0,
                    "loss_tan": 0.001
                }
            ],
            "substrate": {
                "layers": [
                    {
                        "name": "FR4",
                        "permittivity": 4.2,
                        "permeability": 1.0,
                        "loss_tangent": 0.02,
                        "thickness": 1.6e-3
                    }
                ]
            },
            "ports": [
                {"index": 1, "attributes": [1001], "direction": "x", "impedance": 50.0}
            ]
        });
        std::fs::write(&sidecar, serde_json::to_string_pretty(&meta).unwrap())
            .expect("sidecar should be writable");

        convert_project_to_rem(&project, &out_cfg, &out_mesh).expect("conversion should succeed");

        assert!(out_mesh.exists(), "output mesh should exist");
        assert!(out_cfg.exists(), "output config should exist");

        let cfg_text = std::fs::read_to_string(&out_cfg).expect("config should be readable");
        let cfg: Value = serde_json::from_str(&cfg_text).expect("config should be valid json");
        assert_eq!(cfg["Problem"]["Type"], "MoM");
        assert_eq!(cfg["Solver"]["MoM"]["FreqMin"], 1.0e9);
        assert_eq!(cfg["Solver"]["MoM"]["Ports"][0]["Attributes"][0], 1001);
        assert_eq!(cfg["Domains"]["Materials"][0]["Attributes"][0], 10);
        assert_eq!(cfg["Solver"]["MoM"]["Substrate"]["Layers"][0]["Name"], "FR4");

        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn fails_on_invalid_sweep() {
        let dir = unique_temp_dir();
        let project = dir.join("bad_sweep.aedt");
        let sidecar = dir.join("bad_sweep.rem_ansys.json");
        let source_mesh = dir.join("rmsh_out.msh");
        let out_mesh = dir.join("out").join("demo.msh");
        let out_cfg = dir.join("out").join("demo.json");

        std::fs::write(&project, "dummy").expect("project stub should be writable");
        std::fs::write(&source_mesh, "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
            .expect("mesh should be writable");

        let meta = json!({
            "mesh": "rmsh_out.msh",
            "solver": {
                "freq_min": 2.0e9,
                "freq_max": 1.0e9,
                "freq_step": 1.0e8
            }
        });
        std::fs::write(&sidecar, serde_json::to_string_pretty(&meta).unwrap())
            .expect("sidecar should be writable");

        let err = convert_project_to_rem(&project, &out_cfg, &out_mesh)
            .expect_err("invalid sweep should fail");
        assert!(
            err.to_string().contains("invalid frequency sweep"),
            "unexpected error: {err}"
        );

        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn fails_on_port_without_attributes() {
        let dir = unique_temp_dir();
        let project = dir.join("bad_port.aedt");
        let sidecar = dir.join("bad_port.rem_ansys.json");
        let source_mesh = dir.join("rmsh_out.msh");
        let out_mesh = dir.join("out").join("demo.msh");
        let out_cfg = dir.join("out").join("demo.json");

        std::fs::write(&project, "dummy").expect("project stub should be writable");
        std::fs::write(&source_mesh, "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
            .expect("mesh should be writable");

        let meta = json!({
            "mesh": "rmsh_out.msh",
            "ports": [
                {"index": 1, "direction": "x", "impedance": 50.0}
            ]
        });
        std::fs::write(&sidecar, serde_json::to_string_pretty(&meta).unwrap())
            .expect("sidecar should be writable");

        let err = convert_project_to_rem(&project, &out_cfg, &out_mesh)
            .expect_err("missing port attributes should fail");
        assert!(
            err.to_string().contains("must include non-empty 'attributes'"),
            "unexpected error: {err}"
        );

        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn fails_when_no_mesh_source_is_available() {
        let dir = unique_temp_dir();
        let project = dir.join("missing_mesh.aedt");
        let sidecar = dir.join("missing_mesh.rem_ansys.json");
        let out_mesh = dir.join("out").join("demo.msh");
        let out_cfg = dir.join("out").join("demo.json");

        std::fs::write(&project, "dummy").expect("project stub should be writable");

        let meta = json!({
            "solver": {
                "freq_min": 1.0e9,
                "freq_max": 1.1e9,
                "freq_step": 1.0e8
            }
        });
        std::fs::write(&sidecar, serde_json::to_string_pretty(&meta).unwrap())
            .expect("sidecar should be writable");

        let err = convert_project_to_rem(&project, &out_cfg, &out_mesh)
            .expect_err("missing mesh source should fail");
        assert!(
            err.to_string().contains("metadata has no 'mesh'/'model.mesh' path"),
            "unexpected error: {err}"
        );

        let _ = std::fs::remove_dir_all(&dir);
    }
}
