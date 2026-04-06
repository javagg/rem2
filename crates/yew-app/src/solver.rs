use crate::examples;

#[derive(Clone, Debug)]
pub struct SimResult {
    pub energy: f64,
    pub node_count: usize,
    pub max_e: Option<f64>,
    pub max_b: Option<f64>,
}

pub fn run_example(key: &str) -> Result<SimResult, String> {
    let example = examples::find_example(key)
        .ok_or_else(|| format!("Unknown example: {}", key))?;

    let mesh_bytes = examples::get_mesh_bytes(key);
    let value = rem_wasm::run_simulation(example.config_json, &mesh_bytes)
        .map_err(|_| format!("{} solve failed in WASM runtime", example.problem_type))?;

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

    Ok(SimResult {
        energy: result.energy,
        node_count: result.phi.len(),
        max_e,
        max_b,
    })
}
