use rem_config::{load_config_from_str, ConfigFormat, ProblemType};
use rem_mesh::load_mesh_from_bytes;
use rem_materials::DomainMap;
use rem_parallel::NoComm;
use rem_electrostatic::{solve_one as solve_es, postprocess as post_es};
use rem_magnetostatic::solve_one as solve_ms;

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

    let cfg = load_config_from_str(example.config_json, ConfigFormat::Json)
        .map_err(|e| format!("Config error: {}", e))?;

    let mesh_bytes = examples::get_mesh_bytes(key);
    let comm = NoComm;

    let mut mesh = load_mesh_from_bytes(&cfg, &mesh_bytes, &comm)
        .map_err(|e| format!("Mesh error: {}", e))?;
    mesh.partition(&comm);

    let dm = DomainMap::from_config(&cfg)
        .map_err(|e| format!("Domain error: {}", e))?;

    match cfg.problem.problem_type {
        ProblemType::Electrostatic => {
            let phi = solve_es(&cfg, &mesh, &dm, Some(1), 1.0, &comm)
                .map_err(|e| format!("Solve error: {}", e))?;

            let eps_fn = |tag: u32| dm.get(tag).epsilon_abs();
            let energy = post_es::electrostatic_energy(&phi, &mesh, eps_fn);
            let e_field = post_es::gradient_recovery(&phi, &mesh);
            let max_e = e_field.iter()
                .map(|v| (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]).sqrt())
                .fold(0.0f64, f64::max);

            Ok(SimResult {
                energy,
                node_count: phi.len(),
                max_e: Some(max_e),
                max_b: None,
            })
        }
        ProblemType::Magnetostatic => {
            let az = solve_ms(&cfg, &mesh, &dm, Some(1), &comm)
                .map_err(|e| format!("Solve error: {}", e))?;

            let nu_fn = |tag: u32| dm.get(tag).reluctivity();
            let energy = post_es::electrostatic_energy(&az, &mesh, nu_fn);
            let grad_az = post_es::gradient_recovery(&az, &mesh);
            let max_b = grad_az.iter()
                .map(|g| (g[0]*g[0] + g[1]*g[1]).sqrt())
                .fold(0.0f64, f64::max);

            Ok(SimResult {
                energy,
                node_count: az.len(),
                max_e: None,
                max_b: Some(max_b),
            })
        }
        _ => Err(format!("{} solver not yet implemented", example.problem_type)),
    }
}
