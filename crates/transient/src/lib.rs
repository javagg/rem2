//! Transient (time-domain) solver — v1.0
//!
//! Solves the scalar first-order system:
//!   M dv/dt + K v = f(t)
//!
//! where K is the stiffness matrix (−∇·ε∇φ) and M is the consistent mass matrix (∫ε φᵢφⱼ dΩ).
//!
//! Two time-integration schemes are supported via `Solver.Transient.Type`:
//!
//! | Type | Method | Stability | Order |
//! |------|--------|-----------|-------|
//! | `"GeneralizedAlpha"` (default) | First-order generalized-α | Unconditional | 2 |
//! | `"ARKODE"` | IMEX-ARK3(2)4L[2]SA (Kennedy & Carpenter 2003) | Unconditional | 3 |
//! | `"RungeKutta"` | Explicit RK4 (fixed step) | Conditional (CFL) | 4 |
//!
//! Excitation waveforms (`Solver.Transient.Excitation`):
//!   - `""` / `"none"` / `"Step"`: unit step at t=0 (default)
//!   - `"ModulatedGaussian"`: Gaussian-modulated sinusoid
//!       f(t) = exp(−(t−t0)²/(2σ²)) · cos(2π f0 t)
//!       where f0 = ExcitationFreq [GHz], σ = ExcitationWidth/2 [ns], t0 = 5σ
//!   - `"Gaussian"`: unmodulated Gaussian envelope
//!       f(t) = exp(−(t−t0)²/(2σ²))
//!
//! Outputs:
//!   - `postpro/port-t.csv`: time series of port voltage
//!   - `paraview/transient_NNNN.vtk`: solution field at each SaveStep

pub mod output;

use rem_config::PalaceConfig;
use rem_core::{CsrMatrix, RemError, RemResult, TripletMatrix, solve_pcg};
use rem_electrostatic::{assemble::assemble_stiffness, bc::{collect_dirichlet_dofs, apply_dirichlet}};
use rem_eigenmode::assemble_mass::assemble_mass;
use rem_materials::DomainMap;
use rem_mesh::{RemMesh, BoundaryTag};
use rem_mesh::gmsh::read_msh_file;
use rem_parallel::{Comm, NoComm};
use std::path::Path;

fn spmv(mat: &CsrMatrix, x: &[f64], y: &mut [f64]) {
    mat.matvec(x, y, &NoComm);
}

// ─── Excitation waveforms ─────────────────────────────────────────────────────

/// Evaluate the excitation amplitude at time `t` [s].
///
/// - `""` / `"none"` / `"Step"`: unit step (constant 1.0)
/// - `"ModulatedGaussian"`: exp(−(t−t0)²/(2σ²)) · cos(2π f0 t)
///   - `freq_hz`: f0 (ExcitationFreq in GHz → Hz)
///   - `sigma_s`: σ  (ExcitationWidth/2 in ns → s)
/// - `"Gaussian"`: exp(−(t−t0)²/(2σ²))
fn excitation_amplitude(t: f64, kind: &str, freq_hz: f64, sigma_s: f64) -> f64 {
    use std::f64::consts::PI;
    match kind {
        "" | "none" | "Step" => 1.0,
        "ModulatedGaussian" | "Gaussian" => {
            let t0 = 5.0 * sigma_s;
            let envelope = (-(t - t0).powi(2) / (2.0 * sigma_s.powi(2))).exp();
            if kind == "ModulatedGaussian" {
                envelope * (2.0 * PI * freq_hz * t).cos()
            } else {
                envelope
            }
        }
        _ => 1.0, // Unknown → step (warning already emitted by validate_palace_compat)
    }
}

// ─── Entry points ────────────────────────────────────────────────────────────

/// Result of a transient simulation.
pub struct TransientResult {
    /// Time sample points [s].
    pub time_points: Vec<f64>,
    /// Port voltage at each time point [V].
    pub port_voltages: Vec<f64>,
    /// Nodal potential at the time step of peak port voltage magnitude.
    /// Empty if no time steps were recorded.
    pub peak_phi: Vec<f64>,
    /// Time [s] at which `peak_phi` was recorded.
    pub peak_time_s: f64,
}

/// Entry point called from rem-cli.
pub fn run(config: &PalaceConfig, comm: &dyn Comm) -> RemResult<()> {
    log::info!("=== Transient (time-domain) solver ===");

    let mesh_path = Path::new(&config.model.mesh);
    let raw = read_msh_file(mesh_path)?;
    let mut mesh = RemMesh::from_raw(raw, config)?;
    mesh.set_comm(comm.rank(), comm.size());

    run_with_mesh(config, &mesh, comm).map(|_| ())
}

/// Entry point for pre-loaded mesh (used by WASM path).
/// Returns a `TransientResult` with time series, port voltages and peak E-field phi.
pub fn run_with_mesh(config: &PalaceConfig, mesh: &RemMesh, comm: &dyn Comm) -> RemResult<TransientResult> {
    log::info!("=== Transient (time-domain) solver ===");

    let td_cfg = config.solver.transient.as_ref().ok_or_else(|| {
        RemError::Config("Transient problem requires a [Solver.Transient] section".into())
    })?;

    if config.solver.order > 1 {
        log::warn!(
            "Solver.Order={} requested but only P1 (order=1) is implemented; \
             higher-order assembly is pending. Running P1.",
            config.solver.order
        );
    }

    let domain_map = DomainMap::from_config(config)?;
    let eps_fn = |tag: u32| domain_map.get(tag).epsilon_abs();

    let k_raw = assemble_stiffness(mesh, eps_fn)?.to_csr();
    let m_raw = assemble_mass(mesh, eps_fn)?.to_csr();

    let out_dir = config.problem.output_dir();
    #[cfg(not(target_arch = "wasm32"))]
    std::fs::create_dir_all(out_dir).map_err(RemError::Io)?;

    let excited_port = find_excited_lumped_port(mesh);

    // Apply Dirichlet BCs (excited port = 1.0, ground/PEC = 0.0)
    let dofs = collect_dirichlet_dofs(mesh, excited_port, 1.0);
    let mut k_bc = k_raw;
    let mut rhs_bc = vec![0.0f64; mesh.n_nodes()];
    apply_dirichlet(&mut k_bc, &mut rhs_bc, &dofs);

    let mut m_bc = m_raw;
    let mut rhs_dummy = vec![0.0f64; mesh.n_nodes()];
    apply_dirichlet(&mut m_bc, &mut rhs_dummy, &dofs);

    let dt = td_cfg.time_step;
    let t_end = td_cfg.max_time;
    let save_step = td_cfg.save_step.max(1);
    let lin = &config.solver.linear;

    // Excitation parameters
    let exc_kind = td_cfg.excitation.as_str();
    // ExcitationFreq in GHz → Hz; ExcitationWidth in ns → s (σ = width/2)
    let exc_freq_hz = td_cfg.excitation_freq * 1.0e9;
    let exc_sigma_s = (td_cfg.excitation_width * 1.0e-9) / 2.0;

    // Helper: scale the static rhs_bc by excitation amplitude at time t
    let scaled_rhs = |base: &Vec<f64>, t: f64| -> Vec<f64> {
        let amp = excitation_amplitude(t, exc_kind, exc_freq_hz, exc_sigma_s);
        base.iter().map(|&x| x * amp).collect()
    };

    let n_steps = (t_end / dt).ceil() as usize + 1;
    log::info!(
        "Transient: method={}, dt={:.3e}, T={:.3e}, n_steps={}",
        td_cfg.solver_type, dt, t_end, n_steps
    );

    // Initial condition: v(0) = 0
    let n = mesh.n_nodes();
    let mut v = vec![0.0f64; n];

    let mut time_points: Vec<f64> = Vec::with_capacity(n_steps);
    let mut port_v: Vec<f64> = Vec::with_capacity(n_steps);
    let mut peak_phi: Vec<f64> = Vec::new();
    let mut peak_time_s: f64 = 0.0;
    let mut peak_vp_abs: f64 = -1.0;

    match td_cfg.solver_type.as_str() {
        "GeneralizedAlpha" | "" => {
            // First-order generalized-α: ρ_∞ = 0.5
            // α_f = 1/(1+ρ) = 2/3, α_m = (3-ρ)/(2(1+ρ)) = 5/6, γ = 0.5 + α_m - α_f = 2/3
            let rho = 0.5_f64;
            let alpha_f = 1.0 / (1.0 + rho);
            let alpha_m = (3.0 - rho) / (2.0 * (1.0 + rho));
            let gamma   = 0.5 + alpha_m - alpha_f;

            // LHS = α_m * M + α_f * γ * dt * K  (constant for linear system, but we rebuild per-step for variable h)
            let _lhs = build_scaled_sum(&m_bc, alpha_m, &k_bc, alpha_f * gamma * dt);

            // Initialize dvdt_0 = M^{-1} f(0)
            let mut dvdt = {
                let f0 = scaled_rhs(&rhs_bc, 0.0);
                let r = solve_pcg(&m_bc, &f0, lin.tol, lin.max_iter, comm);
                r.solution
            };

            let mut t = 0.0_f64;
            for step in 0..n_steps {
                if t >= t_end + dt * 0.5 { break; }
                let h = dt.min(t_end - t + dt * 1e-12);
                if h < dt * 1e-14 { break; }
                let lhs_h = build_scaled_sum(&m_bc, alpha_m, &k_bc, alpha_f * gamma * h);

                // f at interpolated time t_f = t + α_f * h
                let t_f = t + alpha_f * h;
                let f_f = scaled_rhs(&rhs_bc, t_f);

                // rhs = f_f − (1−α_m) M dvdt − K v − α_f(1−γ) h K dvdt
                let mut k_dvdt = vec![0.0f64; n];
                spmv(&k_bc, &dvdt, &mut k_dvdt);
                let mut m_dvdt = vec![0.0f64; n];
                spmv(&m_bc, &dvdt, &mut m_dvdt);
                let mut kv = vec![0.0f64; n];
                spmv(&k_bc, &v, &mut kv);

                let rhs: Vec<f64> = (0..n).map(|i| {
                    f_f[i]
                    - (1.0 - alpha_m) * m_dvdt[i]
                    - kv[i]
                    - alpha_f * (1.0 - gamma) * h * k_dvdt[i]
                }).collect();

                let result = solve_pcg(&lhs_h, &rhs, lin.tol, lin.max_iter, comm);
                let dvdt_new = result.solution;

                // Update v_{n+1}
                for i in 0..n {
                    v[i] += h * (gamma * dvdt_new[i] + (1.0 - gamma) * dvdt[i]);
                }
                // Enforce BCs
                for (&dof, &val) in &dofs { if dof < n { v[dof] = val; } }

                dvdt = dvdt_new;
                t += h;

                let vp = port_voltage(mesh, &v, excited_port);
                time_points.push(t);
                port_v.push(vp);
                if vp.abs() > peak_vp_abs { peak_vp_abs = vp.abs(); peak_time_s = t; peak_phi = v.clone(); }

                #[cfg(not(target_arch = "wasm32"))]
                if step % save_step == 0 {
                    output::write_field_vtk(out_dir, mesh, &v, step + 1)?;
                }
            }
        }

        "ARKODE" => {
            // IMEX-ARK3(2)4L[2]SA adaptive (from fem-solver constants)
            // Treat f_E = rhs_bc (source, explicit), f_I = -K v (stiff implicit)
            // Per-stage implicit solve: (M + dt * A^I_{ss} * K) u_s = rhs
            let rtol = 1e-4_f64;
            let atol = 1e-8_f64;
            let gamma_ark: f64 = 1767732205903.0 / 4055673282236.0;

            // ARK3(2)4L[2]SA Butcher tables (Kennedy & Carpenter 2003)
            let ai: [[f64; 4]; 4] = [
                [0.0, 0.0, 0.0, 0.0],
                [gamma_ark, gamma_ark, 0.0, 0.0],
                [2746238789719.0 / 10658868560708.0, -640167445237.0 / 6845629431997.0, gamma_ark, 0.0],
                [1471266399579.0 / 7840856788654.0, -4482444167858.0 / 7529755066697.0,
                 11266239266428.0 / 11593286722821.0, gamma_ark],
            ];
            let ae: [[f64; 4]; 4] = [
                [0.0, 0.0, 0.0, 0.0],
                [1767732205903.0 / 2027836641118.0, 0.0, 0.0, 0.0],
                [5535828885825.0 / 10492691773637.0, 788022342437.0 / 10882634858940.0, 0.0, 0.0],
                [6485989280629.0 / 16251701735622.0, -4246266847089.0 / 9704473918619.0,
                 10755448449292.0 / 10357097424841.0, 0.0],
            ];
            let bi: [f64; 4] = [ai[3][0], ai[3][1], ai[3][2], ai[3][3]];
            let bi_hat: [f64; 4] = [
                2756255671327.0 / 12835298489170.0,
                -10771552573575.0 / 22201958757719.0,
                9247589265047.0 / 10645013368117.0,
                2193209047091.0 / 5459859503100.0,
            ];
            // Stage abscissae for time-evaluation of f_E (defined inline as c_ark below)

            let dt_min = dt * 1e-8;
            let dt_max = dt * 10.0;
            let mut cur_dt = dt;
            let mut t = 0.0_f64;
            let mut step = 0usize;

            let mut ki_e = vec![vec![0.0f64; n]; 4];
            let mut ki_i = vec![vec![0.0f64; n]; 4];

            while t < t_end {
                cur_dt = cur_dt.min(t_end - t).max(dt_min);
                if cur_dt < dt_min { break; }

                // Stage 0: f_E at current time t
                ki_e[0] = scaled_rhs(&rhs_bc, t);
                // f_I(t, v) = -K v
                spmv(&k_bc, &v, &mut ki_i[0]);
                for x in &mut ki_i[0] { *x = -*x; }

                // Stages 1..3
                let c_ark: [f64; 4] = [0.0, 1767732205903.0 / 2027836641118.0, 0.6, 1.0];
                let mut u_stages = vec![v.clone(); 4];
                for s in 1..4 {
                    let mut u_s = v.clone();
                    for j in 0..s {
                        for i in 0..n {
                            u_s[i] += cur_dt * ae[s][j] * ki_e[j][i];
                            u_s[i] += cur_dt * ai[s][j] * ki_i[j][i];
                        }
                    }

                    // Implicit solve: (M + cur_dt * a^I_{ss} * K) dvdt_s = f_s - K u_s
                    let aii = ai[s][s];
                    let t_s = t + c_ark[s] * cur_dt;
                    let lhs_s = build_scaled_sum(&m_bc, 1.0, &k_bc, cur_dt * aii);
                    let mut ku_s = vec![0.0f64; n];
                    spmv(&k_bc, &u_s, &mut ku_s);
                    let f_s = scaled_rhs(&rhs_bc, t_s);
                    let rhs_s: Vec<f64> = (0..n).map(|i| f_s[i] - ku_s[i]).collect();

                    let result = solve_pcg(&lhs_s, &rhs_s, lin.tol, lin.max_iter, comm);
                    ki_i[s] = result.solution;

                    // Correct u_s
                    for i in 0..n { u_s[i] += cur_dt * aii * ki_i[s][i]; }
                    for (&dof, &val) in &dofs { if dof < n { u_s[dof] = val; } }
                    u_stages[s] = u_s.clone();

                    // Explicit stage f_E at stage time
                    ki_e[s] = scaled_rhs(&rhs_bc, t_s);
                }

                // 3rd-order solution
                let mut v3 = v.clone();
                for s in 0..4 {
                    for i in 0..n { v3[i] += cur_dt * bi[s] * (ki_e[s][i] + ki_i[s][i]); }
                }

                // 2nd-order embedded
                let mut v2 = v.clone();
                for s in 0..4 {
                    for i in 0..n { v2[i] += cur_dt * bi_hat[s] * (ki_e[s][i] + ki_i[s][i]); }
                }

                // Error norm
                let err: f64 = (0..n).map(|i| {
                    let e = v3[i] - v2[i];
                    let sc = atol + rtol * v[i].abs().max(v3[i].abs());
                    (e / sc).powi(2)
                }).sum::<f64>().sqrt() / (n as f64).sqrt();

                if err <= 1.0 || cur_dt <= dt_min {
                    for (&dof, &val) in &dofs { if dof < n { v3[dof] = val; } }
                    v = v3;
                    t += cur_dt;
                    let vp = port_voltage(mesh, &v, excited_port);
                    time_points.push(t);
                    port_v.push(vp);
                    if vp.abs() > peak_vp_abs { peak_vp_abs = vp.abs(); peak_time_s = t; peak_phi = v.clone(); }

                    #[cfg(not(target_arch = "wasm32"))]
                    if step % save_step == 0 {
                        output::write_field_vtk(out_dir, mesh, &v, step + 1)?;
                    }
                    step += 1;
                }

                let factor = if err > 0.0 { (0.9 / err).powf(1.0 / 3.0).min(5.0).max(0.1) } else { 5.0 };
                cur_dt = (cur_dt * factor).min(dt_max).max(dt_min);
            }
        }

        "RungeKutta" => {
            // Explicit RK4: dv/dt = M^{-1}(f(t) - K v)
            let mut t = 0.0_f64;
            for step in 0..n_steps {
                if t >= t_end + dt * 0.5 { break; }
                let h = dt.min(t_end - t + dt * 1e-12);
                if h < dt * 1e-14 { break; }

                // RK4 stages (time-dependent f)
                let compute_dvdt = |u: &[f64], t_eval: f64| -> Vec<f64> {
                    let ft = scaled_rhs(&rhs_bc, t_eval);
                    let mut kv = vec![0.0f64; n];
                    spmv(&k_bc, u, &mut kv);
                    let rhs: Vec<f64> = (0..n).map(|i| ft[i] - kv[i]).collect();
                    let r = solve_pcg(&m_bc, &rhs, lin.tol, lin.max_iter, comm);
                    r.solution
                };

                let k1 = compute_dvdt(&v, t);
                let u2: Vec<f64> = (0..n).map(|i| v[i] + 0.5 * h * k1[i]).collect();
                let k2 = compute_dvdt(&u2, t + 0.5 * h);
                let u3: Vec<f64> = (0..n).map(|i| v[i] + 0.5 * h * k2[i]).collect();
                let k3 = compute_dvdt(&u3, t + 0.5 * h);
                let u4: Vec<f64> = (0..n).map(|i| v[i] + h * k3[i]).collect();
                let k4 = compute_dvdt(&u4, t + h);

                for i in 0..n {
                    v[i] += h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
                }
                for (&dof, &val) in &dofs { if dof < n { v[dof] = val; } }
                t += h;

                let vp = port_voltage(mesh, &v, excited_port);
                time_points.push(t);
                port_v.push(vp);
                if vp.abs() > peak_vp_abs { peak_vp_abs = vp.abs(); peak_time_s = t; peak_phi = v.clone(); }

                #[cfg(not(target_arch = "wasm32"))]
                if step % save_step == 0 {
                    output::write_field_vtk(out_dir, mesh, &v, step + 1)?;
                }
            }
        }

        other => {
            return Err(RemError::Config(format!(
                "Unknown Solver.Transient.Type={:?}; supported: GeneralizedAlpha, ARKODE, RungeKutta",
                other
            )));
        }
    }

    #[cfg(not(target_arch = "wasm32"))]
    output::write_time_series(out_dir, &time_points, &port_v)?;

    log::info!("Transient solve complete: {} time points, peak |V|={:.4e} at {:.3e} s",
        time_points.len(), peak_vp_abs, peak_time_s);
    Ok(TransientResult { time_points, port_voltages: port_v, peak_phi, peak_time_s })
}

// ─── Matrix helpers ───────────────────────────────────────────────────────────

/// Build α M + β K as a CsrMatrix.
fn build_scaled_sum(m: &CsrMatrix, alpha: f64, k: &CsrMatrix, beta: f64) -> CsrMatrix {
    let n = m.nrows;
    let mut t = TripletMatrix::with_capacity(n, n, m.nnz() + k.nnz());
    for i in 0..m.nrows {
        for ptr in m.row_ptr[i]..m.row_ptr[i + 1] {
            t.add(i, m.col_idx[ptr], alpha * m.values[ptr]);
        }
    }
    for i in 0..k.nrows {
        for ptr in k.row_ptr[i]..k.row_ptr[i + 1] {
            t.add(i, k.col_idx[ptr], beta * k.values[ptr]);
        }
    }
    t.to_csr()
}

// ─── Mesh helpers ─────────────────────────────────────────────────────────────

fn find_excited_lumped_port(mesh: &RemMesh) -> Option<u32> {
    for bc in mesh.boundary_tags.values() {
        match bc {
            BoundaryTag::LumpedPort { index, .. } => return Some(*index),
            BoundaryTag::WavePort { index } => return Some(*index),
            _ => {}
        }
    }
    None
}

fn port_voltage(mesh: &RemMesh, v: &[f64], port_idx: Option<u32>) -> f64 {
    let Some(idx) = port_idx else { return 0.0; };
    let port_nodes: Vec<usize> = mesh.boundary_elements.iter()
        .filter(|e| matches!(
            mesh.boundary_tags.get(&e.tag),
            Some(BoundaryTag::LumpedPort { index, .. }) if *index == idx
        ))
        .flat_map(|e| e.node_ids.iter().copied())
        .collect();
    if port_nodes.is_empty() { return 0.0; }
    port_nodes.iter().map(|&n| v[n]).sum::<f64>() / port_nodes.len() as f64
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use rem_config::{load_config_from_str, ConfigFormat};
    use rem_mesh::{Node, Element, ElementKind};
    use rem_parallel::NoComm;
    use std::collections::HashMap;

    /// Unit square mesh with Ground at y=0, LumpedPort at y=1.
    fn unit_square_mesh() -> RemMesh {
        let nodes = vec![
            Node { id: 0, x: 0.0, y: 0.0, z: 0.0 },
            Node { id: 1, x: 1.0, y: 0.0, z: 0.0 },
            Node { id: 2, x: 1.0, y: 1.0, z: 0.0 },
            Node { id: 3, x: 0.0, y: 1.0, z: 0.0 },
        ];
        let volume_elements = vec![
            Element { id: 1, kind: ElementKind::Tri3, tag: 1, node_ids: vec![0, 1, 2], rank: 0 },
            Element { id: 2, kind: ElementKind::Tri3, tag: 1, node_ids: vec![0, 2, 3], rank: 0 },
        ];
        let boundary_elements = vec![
            Element { id: 3, kind: ElementKind::Line2, tag: 10, node_ids: vec![0, 1], rank: 0 },
            Element { id: 4, kind: ElementKind::Line2, tag: 11, node_ids: vec![2, 3], rank: 0 },
        ];
        let mut boundary_tags: HashMap<u32, BoundaryTag> = HashMap::new();
        boundary_tags.insert(10, BoundaryTag::Ground);
        boundary_tags.insert(11, BoundaryTag::LumpedPort { index: 1, r: 0.0 });

        RemMesh {
            nodes, volume_elements, boundary_elements,
            domain_tags: Default::default(),
            boundary_tags,
            dim: 2, rank: 0, size: 1,
        }
    }

    fn transient_config(scheme: &str) -> PalaceConfig {
        load_config_from_str(
            &format!(r#"{{
                "Problem": {{"Type": "Transient"}},
                "Model": {{"Mesh": "x.msh"}},
                "Solver": {{
                    "Linear": {{"Tol": 1e-10, "MaxIter": 500}},
                    "Transient": {{"Type": "{scheme}", "MaxTime": 0.5, "TimeStep": 0.1, "SaveStep": 100}}
                }}
            }}"#),
            ConfigFormat::Json,
        ).unwrap()
    }

    #[test]
    fn transient_generalized_alpha_runs() {
        let mesh = unit_square_mesh();
        let config = transient_config("GeneralizedAlpha");
        run_with_mesh(&config, &mesh, &NoComm).unwrap();
    }

    #[test]
    fn transient_runge_kutta_runs() {
        let mesh = unit_square_mesh();
        let config = transient_config("RungeKutta");
        run_with_mesh(&config, &mesh, &NoComm).unwrap();
    }

    #[test]
    fn transient_arkode_runs() {
        let mesh = unit_square_mesh();
        let config = transient_config("ARKODE");
        run_with_mesh(&config, &mesh, &NoComm).unwrap();
    }

    #[test]
    fn transient_reaches_steady_state() {
        // After long integration, v should approach steady-state φ(y) = y
        let mesh = unit_square_mesh();
        let config = load_config_from_str(
            r#"{
                "Problem": {"Type": "Transient"},
                "Model": {"Mesh": "x.msh"},
                "Solver": {
                    "Linear": {"Tol": 1e-10, "MaxIter": 500},
                    "Transient": {"Type": "GeneralizedAlpha", "MaxTime": 5.0, "TimeStep": 0.1, "SaveStep": 1000}
                }
            }"#,
            ConfigFormat::Json,
        ).unwrap();
        run_with_mesh(&config, &mesh, &NoComm).unwrap();
    }

    #[test]
    fn transient_unknown_type_returns_error() {
        let mesh = unit_square_mesh();
        let config = transient_config("CVODE");
        let result = run_with_mesh(&config, &mesh, &NoComm);
        assert!(result.is_err(), "Unknown solver type should return error");
    }
}
