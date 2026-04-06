use crate::examples;
use serde::Serialize;

#[derive(Clone, Debug, Serialize)]
pub struct SimResult {
    pub energy: f64,
    pub node_count: usize,
    pub max_e: Option<f64>,
    pub max_b: Option<f64>,
}

#[derive(Clone, Debug)]
pub struct OutputArtifact {
    pub file_name: String,
    pub content: String,
}

#[derive(Clone, Debug)]
pub struct SimRun {
    pub summary: SimResult,
    pub artifacts: Vec<OutputArtifact>,
}

fn phi_csv(phi: &[f64]) -> String {
    let mut out = String::from("index,phi\n");
    for (i, v) in phi.iter().enumerate() {
        out.push_str(&format!("{},{}\n", i, v));
    }
    out
}

fn vec3_csv(name: &str, field: &[[f64; 3]]) -> String {
    let mut out = format!("index,{}_x,{}_y,{}_z,{}_norm\n", name, name, name, name);
    for (i, v) in field.iter().enumerate() {
        let norm = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        out.push_str(&format!("{},{},{},{},{}\n", i, v[0], v[1], v[2], norm));
    }
    out
}

pub fn run_example(key: &str) -> Result<SimRun, String> {
    let example = examples::find_example(key)
        .ok_or_else(|| format!("Unknown example: {}", key))?;

    let mesh_bytes = examples::get_mesh_bytes(key);
    let value = rem_wasm::run_simulation(example.config_json, &mesh_bytes)
        .map_err(|e| format!("{} solve failed in WASM runtime: {:?}", example.problem_type, e))?;

    let result: rem_wasm::SimulationResult = serde_wasm_bindgen::from_value(value)
        .map_err(|e| format!("Failed to decode simulation result: {}", e))?;

    let max_e = result.e_field.as_ref()
        .map(|field| {
            field.iter()
                .map(|v| (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]).sqrt())
                .fold(0.0f64, f64::max)
        });
    let max_b = result.b_field.as_ref()
        .map(|field| {
            field.iter()
                .map(|v| (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]).sqrt())
                .fold(0.0f64, f64::max)
        });

    let mut artifacts = vec![OutputArtifact {
        file_name: "phi.csv".to_string(),
        content: phi_csv(&result.phi),
    }];

    if let Some(field) = &result.e_field {
        artifacts.push(OutputArtifact {
            file_name: "e_field.csv".to_string(),
            content: vec3_csv("e", field),
        });
    }

    if let Some(field) = &result.b_field {
        artifacts.push(OutputArtifact {
            file_name: "b_field.csv".to_string(),
            content: vec3_csv("b", field),
        });
    }

    Ok(SimRun {
        summary: SimResult {
            energy: result.energy,
            node_count: result.phi.len(),
            max_e,
            max_b,
        },
        artifacts,
    })
}
