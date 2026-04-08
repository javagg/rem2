use wasm_bindgen::prelude::*;
use rem_config::{load_config_from_str, ConfigFormat, ProblemType};
use rem_mesh::{load_mesh_from_bytes, gen::{annular_msh, rect_bimaterial_msh}};
use rem_materials::DomainMap;
use rem_parallel::WorldComm;
use rem_electrostatic::{solve_one as solve_es, postprocess as post_es};
use rem_magnetostatic::{solve_one as solve_ms};
use rem_eigenmode::solve as solve_eigen;

extern crate console_error_panic_hook;

#[wasm_bindgen]
pub fn init_panic_hook() {
    console_error_panic_hook::set_once();
}

#[wasm_bindgen]
pub fn init_logger() {
    console_log::init_with_level(log::Level::Info).expect("error initializing logger");
}

#[derive(serde::Serialize, serde::Deserialize, Clone, Debug)]
pub struct SParam {
    pub freq_hz: f64,
    pub s11_re:  f64,
    pub s11_im:  f64,
    /// |S11| in dB
    pub s11_db:  f64,
}

#[derive(serde::Serialize, serde::Deserialize, Clone, Debug)]
pub struct SimulationResult {
    pub phi: Vec<f64>,
    pub energy: f64,
    pub e_field: Option<Vec<[f64; 3]>>,
    pub b_field: Option<Vec<[f64; 3]>>,
    pub frequencies_hz: Option<Vec<f64>>,
    /// Driven: per-frequency S-parameter results
    pub s_params: Option<Vec<SParam>>,
    /// Transient: time points [s] and corresponding port voltages [V]
    pub time_points: Option<Vec<f64>>,
    pub port_voltages: Option<Vec<f64>>,
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
                frequencies_hz: None,
                s_params: None,
                time_points: None, port_voltages: None,
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
                frequencies_hz: None,
                s_params: None,
                time_points: None, port_voltages: None,
            };
            Ok(serde_wasm_bindgen::to_value(&res)?)
        }
        ProblemType::Driven => {
            let freq_results = rem_driven::run_with_mesh(&cfg, &mesh, &comm)
                .map_err(|e| JsError::new(&format!("Driven error: {}", e)))?;

            let s_params: Vec<SParam> = freq_results.iter().map(|r| {
                let mag = (r.s11_re * r.s11_re + r.s11_im * r.s11_im).sqrt();
                let s11_db = if mag > 1e-300 { 20.0 * mag.log10() } else { -300.0 };
                SParam { freq_hz: r.freq_hz, s11_re: r.s11_re, s11_im: r.s11_im, s11_db }
            }).collect();

            let res = SimulationResult {
                phi: vec![], energy: 0.0, e_field: None, b_field: None,
                frequencies_hz: None, s_params: Some(s_params),
                time_points: None, port_voltages: None,
            };
            Ok(serde_wasm_bindgen::to_value(&res)?)
        }
        ProblemType::MoM => {
            rem_mom::run_with_mesh(&cfg,
                cfg.solver.mom.as_ref()
                    .ok_or_else(|| JsError::new("MoM requires Solver.MoM section"))?,
                &mesh,
            ).map_err(|e| JsError::new(&format!("MoM error: {}", e)))?;
            let res = SimulationResult {
                phi: vec![], energy: 0.0, e_field: None, b_field: None,
                frequencies_hz: None, s_params: None,
                time_points: None, port_voltages: None,
            };
            Ok(serde_wasm_bindgen::to_value(&res)?)
        }
        ProblemType::SBR => {
            rem_sbr::run_with_mesh(&cfg,
                cfg.solver.sbr.as_ref()
                    .ok_or_else(|| JsError::new("SBR requires Solver.SBR section"))?,
                &mesh,
            ).map_err(|e| JsError::new(&format!("SBR error: {}", e)))?;
            let res = SimulationResult {
                phi: vec![], energy: 0.0, e_field: None, b_field: None,
                frequencies_hz: None, s_params: None,
                time_points: None, port_voltages: None,
            };
            Ok(serde_wasm_bindgen::to_value(&res)?)
        }
        ProblemType::Eigenmode => {
            let eigen = solve_eigen(&cfg, &mesh, &dm, &comm)
                .map_err(|e| JsError::new(&format!("Eigenmode error: {}", e)))?;

            let phi = eigen.eigenvectors.into_iter().next().unwrap_or_default();
            let res = SimulationResult {
                phi, energy: 0.0, e_field: None, b_field: None,
                frequencies_hz: Some(eigen.frequencies_hz), s_params: None,
                time_points: None, port_voltages: None,
            };
            Ok(serde_wasm_bindgen::to_value(&res)?)
        }
        ProblemType::Transient => {
            let (time_points, port_voltages) = rem_transient::run_with_mesh(&cfg, &mesh, &comm)
                .map_err(|e| JsError::new(&format!("Transient error: {}", e)))?;
            let res = SimulationResult {
                phi: vec![], energy: 0.0, e_field: None, b_field: None,
                frequencies_hz: None, s_params: None,
                time_points: Some(time_points), port_voltages: Some(port_voltages),
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

#[wasm_bindgen]
pub fn get_adapter_mesh() -> Vec<u8> {
    include_bytes!("../../../examples/adapter/mesh/adapter.msh").to_vec()
}

#[wasm_bindgen]
pub fn get_antenna_mesh() -> Vec<u8> {
    include_bytes!("../../../examples/antenna/mesh/antenna.msh").to_vec()
}

#[wasm_bindgen]
pub fn get_coaxial_mesh() -> Vec<u8> {
    include_bytes!("../../../examples/coaxial/mesh/coaxial.msh").to_vec()
}

#[wasm_bindgen]
pub fn get_cpw_mesh() -> Vec<u8> {
    include_bytes!("../../../examples/cpw/mesh/cpw_coax.msh").to_vec()
}

#[wasm_bindgen]
pub fn get_cylinder_mesh() -> Vec<u8> {
    include_bytes!("../../../examples/cylinder/mesh/cylinder_hex.msh").to_vec()
}

#[wasm_bindgen]
pub fn get_transmon_mesh() -> Vec<u8> {
    include_bytes!("../../../examples/transmon/mesh/transmon.msh2").to_vec()
}

#[wasm_bindgen]
pub fn get_spheres_mesh_v2() -> Vec<u8> {
    include_bytes!("../../../examples/spheres/mesh/spheres.msh").to_vec()
}
