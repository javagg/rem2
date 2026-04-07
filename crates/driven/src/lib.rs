//! Driven (frequency-domain) solver — Phase 7 (v0.2)
//!
//! Solves the frequency-domain scalar wave equation:
//!   −∇·(ε ∇φ) − k² ε φ = J_port      (k = ω/c)
//!
//! For each excitation frequency ω = 2πf in [MinFreq, MaxFreq] with step FreqStep:
//!   1. Assemble K (stiffness) and M (mass, consistent or lumped)
//!   2. Build system A = K − k² M  (real symmetric)
//!   3. Apply lumped-port / Dirichlet BCs
//!   4. Solve with PCG (real arithmetic; absorbing BCs treated as simple Dirichlet)
//!   5. Compute port impedance Z, reflection S₁₁
//!   6. Write CSV and (optionally) VTK per save_step
//!
//! Limitations (v0.2):
//!   - Scalar (P1) formulation — valid for planar TEM / quasi-static problems.
//!   - No perfectly-matched layers (PML); Absorbing BC treated as Dirichlet φ=0.
//!   - Real arithmetic only (lossy materials via imaginary conductivity deferred to v0.3).

pub mod output;

use rem_config::PalaceConfig;
use rem_core::{CsrMatrix, RemError, RemResult, TripletMatrix, solve_pcg};
use rem_eigenmode::assemble_mass::assemble_mass;
use rem_electrostatic::{assemble::assemble_stiffness, bc::{collect_dirichlet_dofs, apply_dirichlet}};
use rem_materials::DomainMap;
use rem_mesh::{RemMesh, BoundaryTag, amr};
use rem_mesh::gmsh::read_msh_file;
use rem_parallel::Comm;
use std::path::Path;

const C0: f64 = 2.997_924_58e8;

/// Entry point called from rem-cli.
pub fn run(config: &PalaceConfig, comm: &dyn Comm) -> RemResult<()> {
    log::info!("=== Driven (frequency-domain) solver ===");

    let mesh_path = Path::new(&config.model.mesh);
    let raw = read_msh_file(mesh_path)?;
    let mut mesh = RemMesh::from_raw(raw, config)?;
    mesh.set_comm(comm.rank(), comm.size());

    run_with_mesh(config, &mesh, comm)
}

/// Entry point for pre-loaded mesh (used by WASM path).
pub fn run_with_mesh(config: &PalaceConfig, mesh: &RemMesh, comm: &dyn Comm) -> RemResult<()> {
    log::info!("=== Driven (frequency-domain) solver ===");

    let drv_cfg = config.solver.driven.as_ref().ok_or_else(|| {
        RemError::Config("Driven problem requires a [Solver.Driven] section".into())
    })?;

    if config.solver.order > 1 {
        log::warn!(
            "Solver.Order={} requested but only P1 (order=1) is implemented; \
             higher-order assembly is pending. Running P1.",
            config.solver.order
        );
    }

    let domain_map = DomainMap::from_config(config)?;

    // AMR pre-refinement: refine mesh at center frequency before the full sweep
    let amr_cfg = &config.model.refinement;
    let max_amr_iter = if amr_cfg.max_iter > 0 { amr_cfg.max_iter } else { 0 };
    let amr_theta    = if amr_cfg.tol > 0.0 { amr_cfg.tol } else { 0.5 };

    let refined_mesh: RemMesh;
    let work_mesh: &RemMesh = if max_amr_iter > 0 {
        log::info!("AMR pre-refinement enabled: max_iter={}, θ={}", max_amr_iter, amr_theta);
        let f_center = (drv_cfg.min_freq + drv_cfg.max_freq) * 0.5;
        let k_wave = 2.0 * std::f64::consts::PI * f_center / C0;
        let k2 = k_wave * k_wave;
        let excited_port = find_excited_lumped_port(mesh);

        let eps_fn = |tag: u32| domain_map.get(tag).epsilon_abs();
        let mut cur_mesh = mesh.clone();
        for amr_iter in 1..=max_amr_iter {
            let k_mat = assemble_stiffness(&cur_mesh, eps_fn)?.to_csr();
            let m_mat = assemble_mass(&cur_mesh, eps_fn)?.to_csr();
            let a_mat = shifted_matrix(&k_mat, &m_mat, k2, cur_mesh.n_nodes());
            let mut a_bc = a_mat;
            let dofs = collect_dirichlet_dofs(&cur_mesh, excited_port, 1.0);
            let mut rhs = vec![0.0f64; cur_mesh.n_nodes()];
            apply_dirichlet(&mut a_bc, &mut rhs, &dofs);
            let lin = &config.solver.linear;
            let result = solve_pcg(&a_bc, &rhs, lin.tol, lin.max_iter, comm);
            let phi = result.solution;

            let eta = amr::zz_estimator(&cur_mesh, &phi);
            let total_err: f64 = eta.iter().map(|&e| e * e).sum::<f64>().sqrt();
            log::info!("AMR iter {amr_iter}: nodes={}, |η|={total_err:.3e}", cur_mesh.n_nodes());

            let marked = amr::dorfler_mark(&eta, amr_theta);
            if marked.is_empty() {
                log::info!("AMR pre-refinement converged: no elements marked.");
                break;
            }
            let (fine_mesh, _) = amr::refine_marked(&cur_mesh, &marked);
            cur_mesh = fine_mesh;
        }
        log::info!("AMR pre-refinement complete: {} nodes", cur_mesh.n_nodes());
        refined_mesh = cur_mesh;
        &refined_mesh
    } else {
        mesh
    };

    run_frequency_sweep(config, drv_cfg, work_mesh, &domain_map, comm)
}

// ---------------------------------------------------------------------------
// Frequency sweep (inner loop, used after optional AMR pre-refinement)
// ---------------------------------------------------------------------------

fn run_frequency_sweep(
    config: &PalaceConfig,
    drv_cfg: &rem_config::DrivenSolver,
    mesh: &RemMesh,
    domain_map: &DomainMap,
    comm: &dyn Comm,
) -> RemResult<()> {
    let f_min  = drv_cfg.min_freq;
    let f_max  = drv_cfg.max_freq;
    let f_step = drv_cfg.freq_step;
    let save_step = drv_cfg.save_step.max(1);

    if f_step <= 0.0 || f_min > f_max {
        return Err(RemError::Config(
            "Driven: FreqStep must be > 0 and MinFreq ≤ MaxFreq".into()
        ));
    }

    let n_steps = ((f_max - f_min) / f_step).ceil() as usize + 1;

    // Assemble K and M once (frequency-independent)
    let eps_fn = |tag: u32| domain_map.get(tag).epsilon_abs();
    let k_mat = assemble_stiffness(mesh, eps_fn)?.to_csr();
    let m_mat = assemble_mass(mesh, eps_fn)?.to_csr();

    let out_dir = config.problem.output_dir();
    #[cfg(not(target_arch = "wasm32"))]
    std::fs::create_dir_all(out_dir).map_err(RemError::Io)?;

    let lin = &config.solver.linear;
    let mut freq_results: Vec<FreqResult> = Vec::with_capacity(n_steps);
    let excited_port = find_excited_lumped_port(mesh);

    for step in 0..n_steps {
        let freq = f_min + step as f64 * f_step;
        if freq > f_max + f_step * 0.5 { break; }

        let k_wave = 2.0 * std::f64::consts::PI * freq / C0;
        let k2 = k_wave * k_wave;

        let a_mat = shifted_matrix(&k_mat, &m_mat, k2, mesh.n_nodes());
        let mut a_bc = a_mat;
        let dofs = collect_dirichlet_dofs(mesh, excited_port, 1.0);
        let mut rhs = vec![0.0f64; mesh.n_nodes()];
        apply_dirichlet(&mut a_bc, &mut rhs, &dofs);

        let result = solve_pcg(&a_bc, &rhs, lin.tol, lin.max_iter, comm);
        if !result.converged {
            log::warn!(
                "PCG did not converge at f={:.3e} Hz (iter={}, res={:.3e})",
                freq, result.iterations, result.residual_norm
            );
        }

        let (v_port, i_port) = compute_port_vi(mesh, &result.solution, &k_mat, excited_port);
        let z0 = lumped_port_resistance(mesh, excited_port);
        let s11 = if i_port.abs() > 1e-300 {
            let z = v_port / i_port;
            (z - z0) / (z + z0)
        } else {
            0.0
        };

        log::info!(
            "f={:.3e} Hz  |S11|={:.4}  converged={}",
            freq, s11.abs(), result.converged
        );

        freq_results.push(FreqResult { freq_hz: freq, s11_re: s11, s11_im: 0.0 });

        #[cfg(not(target_arch = "wasm32"))]
        if step % save_step == 0 {
            output::write_field_vtk(out_dir, mesh, &result.solution, step + 1)?;
        }
    }

    #[cfg(not(target_arch = "wasm32"))]
    output::write_s_params(out_dir, &freq_results)?;
    log::info!("Driven solve complete: {} frequency points", freq_results.len());
    Ok(())
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

pub(crate) struct FreqResult {
    freq_hz: f64,
    s11_re:  f64,
    s11_im:  f64,
}

/// Build A = K − k² M
fn shifted_matrix(k: &CsrMatrix, m: &CsrMatrix, sigma: f64, n: usize) -> CsrMatrix {
    if sigma.abs() < 1e-300 {
        return k.clone();
    }
    let mut t = TripletMatrix::with_capacity(n, n, k.nnz() + m.nnz());
    for i in 0..k.nrows {
        for ptr in k.row_ptr[i]..k.row_ptr[i+1] {
            t.add(i, k.col_idx[ptr], k.values[ptr]);
        }
    }
    for i in 0..m.nrows {
        for ptr in m.row_ptr[i]..m.row_ptr[i+1] {
            t.add(i, m.col_idx[ptr], -sigma * m.values[ptr]);
        }
    }
    t.to_csr()
}

/// Find the first excited lumped port or wave port index in the mesh.
fn find_excited_lumped_port(mesh: &RemMesh) -> Option<u32> {
    for bc in mesh.boundary_tags.values() {
        match bc {
            BoundaryTag::LumpedPort { index, .. } => return Some(*index),
            BoundaryTag::WavePort { index } => {
                log::warn!(
                    "WavePort index={}: using TEM approximation (Dirichlet φ=V). \
                     Full TE/TM modal field matching not yet implemented.",
                    index
                );
                return Some(*index);
            }
            _ => {}
        }
    }
    None
}

/// Get port resistance from config (default 50 Ω).
/// WavePort is treated as matched 50-Ω load for TEM approximation.
fn lumped_port_resistance(mesh: &RemMesh, port_idx: Option<u32>) -> f64 {
    if let Some(idx) = port_idx {
        for bc in mesh.boundary_tags.values() {
            match bc {
                BoundaryTag::LumpedPort { index, r } if *index == idx => {
                    return if *r > 0.0 { *r } else { 50.0 };
                }
                BoundaryTag::WavePort { index } if *index == idx => {
                    return 50.0;  // TEM characteristic impedance
                }
                _ => {}
            }
        }
    }
    50.0
}

/// Compute port voltage and current from solution.
/// V = average φ on port nodes (after BC application, = 1.0 for excited port).
/// I = sum of K[port_node, :] * φ  (outgoing current).
fn compute_port_vi(
    mesh: &RemMesh,
    phi: &[f64],
    k: &CsrMatrix,
    port_idx: Option<u32>,
) -> (f64, f64) {
    let Some(idx) = port_idx else { return (0.0, 0.0); };

    let port_nodes: Vec<usize> = mesh.boundary_elements.iter()
        .filter(|e| matches!(mesh.boundary_tags.get(&e.tag), Some(BoundaryTag::LumpedPort { index, .. }) if *index == idx))
        .flat_map(|e| e.node_ids.iter().copied())
        .collect();

    if port_nodes.is_empty() { return (0.0, 0.0); }

    // V = mean(φ at port nodes)
    let v_port = port_nodes.iter().map(|&n| phi[n]).sum::<f64>() / port_nodes.len() as f64;

    // I = Σ_{n ∈ port} (K * φ)[n]   (net current out of port)
    let mut i_port = 0.0;
    for &n in &port_nodes {
        let kphi_n: f64 = (k.row_ptr[n]..k.row_ptr[n + 1])
            .map(|ptr| k.values[ptr] * phi[k.col_idx[ptr]])
            .sum();
        i_port += kphi_n;
    }

    (v_port, i_port)
}
