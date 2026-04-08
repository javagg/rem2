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
pub struct RcsPoint {
    pub theta_deg: f64,
    pub phi_deg:   f64,
    pub rcs_m2:    f64,
    pub rcs_dbsm:  f64,
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
    /// Eigenmode: Q-factors from dielectric loss perturbation (None if lossless)
    pub q_factors: Option<Vec<f64>>,
    /// Eigenmode: all eigenvectors (one per mode); phi holds mode 0 for backwards compat
    pub eigenvectors: Option<Vec<Vec<f64>>>,
    /// MoM/SBR: RCS pattern data, one entry per frequency
    pub rcs_data: Option<Vec<(f64, Vec<RcsPoint>)>>,
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
                time_points: None, port_voltages: None, q_factors: None, rcs_data: None, eigenvectors: None,
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
                time_points: None, port_voltages: None, q_factors: None, rcs_data: None, eigenvectors: None,
            };
            Ok(serde_wasm_bindgen::to_value(&res)?)
        }
        ProblemType::Driven => {
            let driven = rem_driven::run_with_mesh(&cfg, &mesh, &comm)
                .map_err(|e| JsError::new(&format!("Driven error: {}", e)))?;

            let s_params: Vec<SParam> = driven.freq_results.iter().map(|r| {
                let mag = (r.s11_re * r.s11_re + r.s11_im * r.s11_im).sqrt();
                let s11_db = if mag > 1e-300 { 20.0 * mag.log10() } else { -300.0 };
                SParam { freq_hz: r.freq_hz, s11_re: r.s11_re, s11_im: r.s11_im, s11_db }
            }).collect();

            // E-field and energy at the peak-|S11| frequency
            let (e_field, energy) = if !driven.peak_phi.is_empty() {
                let e = post_es::gradient_recovery(&driven.peak_phi, &mesh);
                let eps_fn = |tag: u32| dm.get(tag).epsilon_abs();
                let u = post_es::electrostatic_energy(&driven.peak_phi, &mesh, eps_fn);
                (Some(e), u)
            } else {
                (None, 0.0)
            };

            let res = SimulationResult {
                phi: vec![], energy, e_field, b_field: None,
                frequencies_hz: None, s_params: Some(s_params),
                time_points: None, port_voltages: None, q_factors: None, rcs_data: None, eigenvectors: None,
            };
            Ok(serde_wasm_bindgen::to_value(&res)?)
        }
        ProblemType::MoM => {
            let mom_result = rem_mom::run_with_mesh(&cfg,
                cfg.solver.mom.as_ref()
                    .ok_or_else(|| JsError::new("MoM requires Solver.MoM section"))?,
                &mesh,
            ).map_err(|e| JsError::new(&format!("MoM error: {}", e)))?;
            let rcs_data: Vec<(f64, Vec<RcsPoint>)> = mom_result.rcs.into_iter().map(|(f, pts)| {
                let pts2 = pts.into_iter().map(|p| RcsPoint {
                    theta_deg: p.theta_deg, phi_deg: p.phi_deg,
                    rcs_m2: p.rcs_m2, rcs_dbsm: p.rcs_dbsm,
                }).collect();
                (f, pts2)
            }).collect();
            let res = SimulationResult {
                phi: vec![], energy: 0.0, e_field: None, b_field: None,
                frequencies_hz: None, s_params: None,
                time_points: None, port_voltages: None, q_factors: None,
                rcs_data: Some(rcs_data), eigenvectors: None,
            };
            Ok(serde_wasm_bindgen::to_value(&res)?)
        }
        ProblemType::SBR => {
            let sbr_result = rem_sbr::run_with_mesh(&cfg,
                cfg.solver.sbr.as_ref()
                    .ok_or_else(|| JsError::new("SBR requires Solver.SBR section"))?,
                &mesh,
            ).map_err(|e| JsError::new(&format!("SBR error: {}", e)))?;
            let rcs_data: Vec<(f64, Vec<RcsPoint>)> = sbr_result.rcs.into_iter().map(|(f, pts)| {
                let pts2 = pts.into_iter().map(|p| RcsPoint {
                    theta_deg: p.theta_deg, phi_deg: p.phi_deg,
                    rcs_m2: p.rcs_m2, rcs_dbsm: p.rcs_dbsm,
                }).collect();
                (f, pts2)
            }).collect();
            let res = SimulationResult {
                phi: vec![], energy: 0.0, e_field: None, b_field: None,
                frequencies_hz: None, s_params: None,
                time_points: None, port_voltages: None, q_factors: None,
                rcs_data: Some(rcs_data), eigenvectors: None,
            };
            Ok(serde_wasm_bindgen::to_value(&res)?)
        }
        ProblemType::Eigenmode => {
            let eigen = solve_eigen(&cfg, &mesh, &dm, &comm)
                .map_err(|e| JsError::new(&format!("Eigenmode error: {}", e)))?;

            let phi = eigen.eigenvectors.first().cloned().unwrap_or_default();
            let all_vecs = eigen.eigenvectors;
            let res = SimulationResult {
                phi, energy: 0.0, e_field: None, b_field: None,
                frequencies_hz: Some(eigen.frequencies_hz), s_params: None,
                time_points: None, port_voltages: None,
                q_factors: eigen.q_factors, rcs_data: None,
                eigenvectors: Some(all_vecs),
            };
            Ok(serde_wasm_bindgen::to_value(&res)?)
        }
        ProblemType::Transient => {
            let transient = rem_transient::run_with_mesh(&cfg, &mesh, &comm)
                .map_err(|e| JsError::new(&format!("Transient error: {}", e)))?;

            let (e_field, energy) = if !transient.peak_phi.is_empty() {
                let e = post_es::gradient_recovery(&transient.peak_phi, &mesh);
                let eps_fn = |tag: u32| dm.get(tag).epsilon_abs();
                let u = post_es::electrostatic_energy(&transient.peak_phi, &mesh, eps_fn);
                (Some(e), u)
            } else {
                (None, 0.0)
            };

            let res = SimulationResult {
                phi: transient.peak_phi, energy, e_field, b_field: None,
                frequencies_hz: None, s_params: None,
                time_points: Some(transient.time_points),
                port_voltages: Some(transient.port_voltages),
                q_factors: None, rcs_data: None, eigenvectors: None,
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
