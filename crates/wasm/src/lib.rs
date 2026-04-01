use wasm_bindgen::prelude::*;
use rem_config::{load_config_from_str, ConfigFormat, ProblemType};
use rem_mesh::{load_mesh_from_bytes, gen::{annular_msh, rect_bimaterial_msh}};
use rem_materials::DomainMap;
use rem_parallel::WorldComm;
use rem_electrostatic::{solve_one as solve_es, postprocess as post_es};
use rem_magnetostatic::{solve_one as solve_ms};

extern crate console_error_panic_hook;

#[wasm_bindgen]
pub fn init_panic_hook() {
    console_error_panic_hook::set_once();
}

#[wasm_bindgen]
pub fn init_logger() {
    console_log::init_with_level(log::Level::Info).expect("error initializing logger");
}

#[derive(serde::Serialize)]
pub struct SimulationResult {
    pub phi: Vec<f64>,
    pub energy: f64,
    pub e_field: Option<Vec<[f64; 3]>>,
    pub b_field: Option<Vec<[f64; 3]>>,
}

#[wasm_bindgen]
pub fn run_simulation(config_json: &str, mesh_bytes: &[u8]) -> Result<JsValue, JsError> {
    let cfg = load_config_from_str(config_json, ConfigFormat::Json)
        .map_err(|e| JsError::new(&format!("Config error: {}", e)))?;

    let comm = WorldComm::new();

    let mesh = load_mesh_from_bytes(&cfg, mesh_bytes, &comm)
        .map_err(|e| JsError::new(&format!("Mesh error: {}", e)))?;

    let dm = DomainMap::from_config(&cfg)
        .map_err(|e| JsError::new(&format!("Domain error: {}", e)))?;

    match cfg.problem.problem_type {
        ProblemType::Electrostatic => {
            let phi = solve_es(&cfg, &mesh, &dm, Some(1), 1.0, &comm)
                .map_err(|e| JsError::new(&format!("Solve error: {}", e)))?;

            let eps_fn = |tag: u32| dm.get(tag).epsilon_abs();
            let energy = post_es::electrostatic_energy(&phi, &mesh, eps_fn);
            let e_field = Some(post_es::gradient_recovery(&phi, &mesh));

            let res = SimulationResult {
                phi,
                energy,
                e_field,
                b_field: None,
            };
            Ok(serde_wasm_bindgen::to_value(&res)?)
        }
        ProblemType::Magnetostatic => {
            let az = solve_ms(&cfg, &mesh, &dm, Some(1), &comm)
                .map_err(|e| JsError::new(&format!("Solve error: {}", e)))?;

            let nu_fn = |tag: u32| dm.get(tag).reluctivity();
            let energy = post_es::electrostatic_energy(&az, &mesh, nu_fn);
            let grad_az = post_es::gradient_recovery(&az, &mesh);
            let b_field: Vec<[f64; 3]> = grad_az.iter()
                .map(|g| [-g[1], g[0], 0.0])
                .collect();

            let res = SimulationResult {
                phi: az,
                energy,
                e_field: None,
                b_field: Some(b_field),
            };
            Ok(serde_wasm_bindgen::to_value(&res)?)
        }
        _ => Err(JsError::new("Unsupported problem type")),
    }
}

#[wasm_bindgen]
pub fn get_spheres_mesh() -> Vec<u8> {
    annular_msh(1.0, 4.0, 10, 32, 1, 2, 10).into_bytes()
}

#[wasm_bindgen]
pub fn get_rings_mesh() -> Vec<u8> {
    rect_bimaterial_msh(1.0, 1.0, 20, 20, 1, 2, 10, 20).into_bytes()
}
