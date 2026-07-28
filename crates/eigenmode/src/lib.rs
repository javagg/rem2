//! Eigenmode solver — Phase 7 (v0.2)
//!
//! Solves the generalized eigenvalue problem:
//!   K x = λ M x
//!
//! where:
//!   K = stiffness matrix (−∇·ε∇φ = 0 Laplacian)
//!   M = mass matrix (∫ ε φ_i φ_j dΩ)
//!   λ = (2π f)² / c²   (eigenfrequency²  in units of [rad/m]²)
//!
//! Algorithm: Lanczos iteration with shift-invert.
//!   - Shift σ = (2π·target_freq)² · ε₀μ₀  (converts Hz → [rad/m]²)
//!   - Solve (K − σ M)^{-1} M v using PCG inner solver
//!   - Accumulate m-step Lanczos basis, build tridiagonal T
//!   - Solve T eigenvalue problem with nalgebra (small, dense)
//!   - Ritz eigenpairs → physical frequencies
//!
//! Output (CSV + VTK) mirrors the Palace port-post format.

pub mod assemble_mass;
mod hcurl;
pub mod output;

use rem_config::PalaceConfig;
use rem_core::{CsrMatrix, RemError, RemResult, TripletMatrix, report_peak_memory, solve_pcg, timing};
use rem_electrostatic::{assemble::assemble_stiffness, bc::{collect_dirichlet_dofs, apply_dirichlet, collect_periodic_node_pairs, apply_periodic}};
use rem_materials::DomainMap;
use rem_mesh::{RemMesh, ElementKind, amr, refine_marked_tri3};
use rem_mesh::gmsh::read_msh_file;
use rem_parallel::Comm;
use nalgebra::DMatrix;
use std::path::Path;

// Physical constants
const C0: f64 = 2.997_924_58e8;  // speed of light [m/s]

// AMR eigenfrequency convergence tolerance: stop if max relative freq change < this
const AMR_FREQ_TOL: f64 = 1e-4;
/// AMR absolute frequency convergence tolerance [Hz]: also stop if max absolute change < this
const AMR_FREQ_ABS_TOL: f64 = 1e6; // 1 MHz

/// Entry point called from rem-cli.
pub fn run(config: &PalaceConfig, comm: &dyn Comm) -> RemResult<()> {
    let _span = timing::span("eigenmode.total");
    log::info!("\n=== Eigenmode (frequency-domain) solver ===\n");
    if config.solver.uses_hcurl_for_eigenmode() {
        log::info!(
            "Eigenmode solver using HCurl/Nedelec edge-element discretization."
        );
    } else {
        log::warn!(
            "Eigenmode solver: scalar H1 nodal discretization selected (effective formulation=H1). \
             Non-physical spurious modes may appear in cavity/waveguide solutions. Use \"HCurl\" (the default)."
        );
    }

    let eig_cfg = config.solver.eigenmode.as_ref().ok_or_else(|| {
        RemError::Config("Eigenmode problem requires a [Solver.Eigenmode] section".into())
    })?;

    if config.solver.order > 2 {
        log::warn!(
            "Solver.Order={} requested; P1 and P2 (Tet10/Tri6) are implemented. \
             Order≥3 is not yet supported — running P2.",
            config.solver.order
        );
    } else if config.solver.order == 2 {
        log::info!("Solver.Order=2: using P2 quadratic assembly for Tet10/Tri6 elements.");
    }

    // Report solver configuration
    log::info!("Solver configuration:");
    log::info!("  Target frequency = {:.3e} Hz", eig_cfg.target);
    log::info!("  Number of modes  = {}", eig_cfg.n);
    log::info!("  Tolerance        = {:.3e}", eig_cfg.tol);
    log::info!("  Max iterations   = {}", eig_cfg.max_iter);
    log::info!("");

    let mesh_path = Path::new(&config.model.mesh);
    let raw = read_msh_file(mesh_path)?;
    let mut mesh = RemMesh::from_raw(raw, config.model.l0)?;
    mesh.set_comm(comm.rank(), comm.size());
    let domain_map = DomainMap::from_config(config)?;

    // AMR loop: refine mesh adaptively and re-solve
    let amr_cfg = &config.model.refinement;
    let max_amr_iter = if amr_cfg.max_iter > 0 { amr_cfg.max_iter } else { 0 };
    let amr_theta    = if amr_cfg.tol > 0.0 { amr_cfg.tol } else { 0.5 };
    if config.solver.uses_hcurl_for_eigenmode() && max_amr_iter > 0 {
        log::info!(
            "HCurl AMR enabled: using element-area indicator (frequency-convergence criterion). \
             Refinement proceeds until Δf/f < {:.1e} or max {} iterations.",
            AMR_FREQ_TOL, max_amr_iter
        );
    }

    let (final_mesh, result) = if max_amr_iter > 0 {
        log::info!("Adaptive mesh refinement (AMR):");
        log::info!("  Max iterations = {}", max_amr_iter);
        log::info!("  Dörfler marking = {:.1}%", amr_theta * 100.0);
        log::info!("");

        let mut cur_mesh = mesh;
        let mut result = solve(config, &cur_mesh, &domain_map, comm)?;
        let mut prev_freqs = result.frequencies_hz.clone();

        for amr_iter in 1..=max_amr_iter {
            // Choose error indicator based on discretisation:
            // · H1 path: ZZ gradient recovery on nodal solution.
            // · HCurl path: element-area-based indicator (no nodal field available;
            //   frequency-convergence criterion provides the stopping condition).
            let (eta, total_err) = if result.is_hcurl {
                let e = hcurl_element_area_indicator(&cur_mesh);
                let total = e.iter().map(|&v| v * v).sum::<f64>().sqrt();
                (e, total)
            } else if let Some(phi) = result.eigenvectors.first() {
                let e = amr::zz_estimator(&cur_mesh, phi);
                let total = e.iter().map(|&v| v * v).sum::<f64>().sqrt();
                (e, total)
            } else {
                break;
            };
            log::info!("AMR iteration {}: {} nodes, indicator = {:.3e}",
                amr_iter, cur_mesh.n_nodes(), total_err);

            let marked = amr::dorfler_mark(&eta, amr_theta);
            if marked.is_empty() {
                log::info!("  → No elements marked; stopping AMR.");
                break;
            }
            {
                let (fine_mesh, _midpoints) = refine_amr_mesh(&cur_mesh, &marked);
                result = solve(config, &fine_mesh, &domain_map, comm)?;
                cur_mesh = fine_mesh;
            }

            // Check eigenfrequency convergence between AMR iterations
            let max_rel_change = prev_freqs.iter()
                    .zip(result.frequencies_hz.iter())
                    .map(|(&f_old, &f_new)| {
                        if f_old.abs() > 1e-30 { ((f_new - f_old) / f_old).abs() } else { 0.0 }
                    })
                    .fold(0.0f64, f64::max);
            let max_abs_change_hz = prev_freqs.iter()
                    .zip(result.frequencies_hz.iter())
                    .map(|(&f_old, &f_new)| (f_new - f_old).abs())
                    .fold(0.0f64, f64::max);
            log::info!(
                    "  → Frequency change: {:.3e} (rel), {:.3e} Hz (abs)",
                    max_rel_change, max_abs_change_hz
                );
            if max_rel_change < AMR_FREQ_TOL || max_abs_change_hz < AMR_FREQ_ABS_TOL {
                log::info!(
                    "  → AMR converged (freq change < tolerance)"
                );
                break;
            }
            prev_freqs = result.frequencies_hz.clone();
        }
        (cur_mesh, result)
    } else {
        let result = solve(config, &mesh, &domain_map, comm)?;
        (mesh, result)
    };

    let out_dir = config.problem.output_dir();
    std::fs::create_dir_all(out_dir).map_err(RemError::Io)?;
    output::write_eigenfrequencies(out_dir, &result)?;

    let save_n = eig_cfg.save.min(result.frequencies_hz.len());
    for (mode_idx, phi) in result.eigenvectors.iter().enumerate().take(save_n) {
        if result.is_hcurl {
            let order = config.solver.eigenmode_hcurl_order().clamp(1, 2) as u8;
            output::write_mode_vector_vtk(out_dir, &final_mesh, phi, mode_idx + 1, order)?;
        } else {
            output::write_mode_vtk(out_dir, &final_mesh, phi, mode_idx + 1)?;
        }
    }

    // Field probes (Domains.Postprocessing.Probe) — one row per (mode, probe)
    if !result.is_hcurl {
        if let Some(dp) = &config.domains.postprocessing {
        if !dp.probe.is_empty() {
            let probes_input: Vec<(u32, [f64; 3])> = dp.probe.iter().map(|p| {
                let c = &p.center;
                let xyz = [c.first().copied().unwrap_or(0.0),
                           c.get(1).copied().unwrap_or(0.0),
                           c.get(2).copied().unwrap_or(0.0)];
                (p.index, xyz)
            }).collect();
            let n_modes = result.eigenvectors.len();
            let mut mode_probes: Vec<(usize, Vec<rem_electrostatic::postprocess::ProbeValue>)> =
                Vec::with_capacity(n_modes);
            for (mode_idx, phi) in result.eigenvectors.iter().enumerate() {
                let probe_vals = rem_electrostatic::postprocess::evaluate_probes(
                    phi, &final_mesh, &probes_input,
                );
                mode_probes.push((mode_idx + 1, probe_vals));
            }
            rem_electrostatic::postprocess::write_probe_modal_csv(
                std::path::Path::new(out_dir), &mode_probes,
            ).map_err(RemError::Io)?;
        }
        }
    } else if config.domains.postprocessing.as_ref().map(|p| !p.probe.is_empty()).unwrap_or(false) {
        log::warn!("Probe postprocessing is currently disabled for HCurl eigenvectors.");
    }

    log::info!("");
    log::info!("Eigenmode solve complete:");
    log::info!("  {} modes computed", result.frequencies_hz.len());
    log::info!("  {} modes saved to output/", save_n);
    report_peak_memory("Eigenmode solver");
    Ok(())
}

/// Element-area error indicator for HCurl (Nedelec) eigenmode AMR.
///
/// Returns the area (2-D) or volume (3-D) of each element as a proxy for
/// the local discretisation error.  Larger elements are marked first, driving
/// the mesh toward a uniform element size.  The physical stopping criterion
/// is the eigenfrequency convergence check in the outer AMR loop.
///
/// This replaces the ZZ gradient-recovery estimator which requires nodal H1
/// fields and cannot be applied to edge-DOF Nedelec solutions directly.
fn hcurl_element_area_indicator(mesh: &RemMesh) -> Vec<f64> {
    use rem_mesh::ElementKind;
    mesh.volume_elements.iter().map(|elem| {
        match elem.kind {
            ElementKind::Tri3 => {
                // |Tri| = ½ |det[p1-p0, p2-p0]|
                let nodes = &mesh.nodes;
                if elem.node_ids.len() < 3 { return 1.0; }
                let (i0, i1, i2) = (elem.node_ids[0], elem.node_ids[1], elem.node_ids[2]);
                let (n0, n1, n2) = (&nodes[i0], &nodes[i1], &nodes[i2]);
                let ax = n1.x - n0.x; let ay = n1.y - n0.y;
                let bx = n2.x - n0.x; let by = n2.y - n0.y;
                (ax * by - ay * bx).abs() * 0.5
            }
            ElementKind::Tet4 => {
                // |Tet| = ⅙ |det[p1-p0, p2-p0, p3-p0]|
                let nodes = &mesh.nodes;
                if elem.node_ids.len() < 4 { return 1.0; }
                let (i0, i1, i2, i3) = (elem.node_ids[0], elem.node_ids[1],
                                         elem.node_ids[2], elem.node_ids[3]);
                let (n0, n1, n2, n3) = (&nodes[i0], &nodes[i1], &nodes[i2], &nodes[i3]);
                let ax = n1.x-n0.x; let ay = n1.y-n0.y; let az = n1.z-n0.z;
                let bx = n2.x-n0.x; let by = n2.y-n0.y; let bz = n2.z-n0.z;
                let cx = n3.x-n0.x; let cy = n3.y-n0.y; let cz = n3.z-n0.z;
                let det = ax*(by*cz - bz*cy) - ay*(bx*cz - bz*cx) + az*(bx*cy - by*cx);
                det.abs() / 6.0
            }
            // For higher-order or other elements, return a unit weight (uniform refine).
            _ => 1.0,
        }
    }).collect()
}

fn refine_amr_mesh(
    mesh: &RemMesh,
    marked: &[usize],
) -> (RemMesh, std::collections::HashMap<(usize, usize), usize>) {
    if mesh.dim == 2
        && mesh.volume_elements.iter().all(|element| element.kind == ElementKind::Tri3)
        && mesh.boundary_elements.iter().all(|element| element.kind == ElementKind::Line2)
    {
        match refine_marked_tri3(mesh, marked) {
            Ok((fine_mesh, midpoint_map)) => {
                log::info!("AMR refine backend: fem-rs Tri3 bridge");
                return (fine_mesh, midpoint_map);
            }
            Err(err) => {
                log::warn!(
                    "fem-rs Tri3 bridge refinement failed ({}); falling back to legacy AMR",
                    err
                );
            }
        }
    }

    amr::refine_marked(mesh, marked)
}

/// Result of an eigenmode solve.
pub struct EigenResult {
    /// Eigenfrequencies in Hz.
    pub frequencies_hz: Vec<f64>,
    /// Corresponding eigenvectors (one per mode, length = n_nodes).
    pub eigenvectors: Vec<Vec<f64>>,
    /// Q-factors from dielectric loss perturbation (1/Q = tan_δ_eff).
    /// None if all materials are lossless.
    pub q_factors: Option<Vec<f64>>,
    /// True when eigenvectors are HCurl edge-DOF vectors (not nodal H1 fields).
    pub is_hcurl: bool,
}

/// Solve the generalized eigenvalue problem for `config` + pre-loaded mesh.
pub fn solve(
    config: &PalaceConfig,
    mesh: &RemMesh,
    domain_map: &DomainMap,
    comm: &dyn Comm,
) -> RemResult<EigenResult> {
    if config.solver.uses_hcurl_for_eigenmode() {
        return hcurl::solve_hcurl(config, mesh, domain_map, comm);
    }

    // P-refinement: if Order=2, promote the P1 mesh to P2 (Tri3→Tri6, Tet4→Tet10).
    // This wires the high-order FEM path without changing any assembly or BC code.
    let p2_mesh_owned: RemMesh;
    let work_mesh: &RemMesh = if config.solver.order >= 2 {
        let is_already_p2 = mesh.volume_elements.iter().all(|e| {
            matches!(e.kind, rem_mesh::ElementKind::Tri6 | rem_mesh::ElementKind::Tet10
                           | rem_mesh::ElementKind::Quad4 | rem_mesh::ElementKind::Hex8)
        });
        if is_already_p2 {
            mesh
        } else {
            log::info!(
                "Solver.Order=2: promoting P1 mesh ({} nodes) to P2 (adding edge midpoints).",
                mesh.n_nodes()
            );
            p2_mesh_owned = rem_mesh::p_refine_mesh(mesh);
            log::info!(
                "P2 mesh: {} nodes ({} added), {} elements.",
                p2_mesh_owned.n_nodes(),
                p2_mesh_owned.n_nodes() - mesh.n_nodes(),
                p2_mesh_owned.n_volume_elements()
            );
            &p2_mesh_owned
        }
    } else {
        mesh
    };

    let eig_cfg = config.solver.eigenmode.as_ref().ok_or_else(|| {
        RemError::Config("missing Eigenmode solver config".into())
    })?;

    let n = work_mesh.n_nodes();
    let n_modes = eig_cfg.n;
    let _tol    = eig_cfg.tol;  // reserved for iterative eigensolvers (not yet used)
    let target_hz = eig_cfg.target;

    // Angular frequency target → shift σ = (ω/c)²
    let sigma = if target_hz > 0.0 {
        let omega = 2.0 * std::f64::consts::PI * target_hz;
        (omega / C0) * (omega / C0)
    } else {
        0.0
    };

    // Assemble stiffness K and mass M — scalar or tensor path
    let mut k_triplet = if domain_map.any_anisotropic() {
        log::info!("Anisotropic material(s) detected — using tensor stiffness assembly.");
        use rem_electrostatic::assemble::assemble_stiffness_aniso;
        let tensor_fn = |tag: u32| domain_map.get(tag).epsilon_tensor;
        assemble_stiffness_aniso(work_mesh, tensor_fn)?
    } else {
        let eps_fn = |tag: u32| domain_map.get(tag).epsilon_abs();
        assemble_stiffness(work_mesh, eps_fn)?
    };
    let eps_fn = |tag: u32| domain_map.get(tag).epsilon_abs();
    let mut m_triplet = if domain_map.any_anisotropic() {
        let tensor_fn = |tag: u32| domain_map.get(tag).epsilon_tensor;
        assemble_mass::assemble_mass_aniso(work_mesh, tensor_fn)?
    } else {
        assemble_mass::assemble_mass(work_mesh, eps_fn)?
    };

    // Apply periodic node remapping before converting to CSR
    let periodic_pairs = collect_periodic_node_pairs(work_mesh, config);
    if !periodic_pairs.is_empty() {
        k_triplet.remap_periodic_nodes(&periodic_pairs);
        m_triplet.remap_periodic_nodes(&periodic_pairs);
    }

    let mut k_mat = k_triplet.to_csr();
    let m_mat = m_triplet.to_csr();

    // Collect Dirichlet DOFs (PEC / Ground → φ=0)
    let mut dofs = collect_dirichlet_dofs(work_mesh, None, 0.0);
    if !periodic_pairs.is_empty() {
        apply_periodic(&mut dofs, &periodic_pairs);
    }

    // Build shifted matrix A = K − σ M (before applying BCs to get scaling right)
    let a_mat = shifted_matrix(&k_mat, &m_mat, sigma, n);

    // Apply BCs to A and to K (for residual checks)
    let mut a_bc = a_mat;
    let mut rhs_dummy = vec![0.0f64; n];
    apply_dirichlet(&mut a_bc, &mut rhs_dummy, &dofs);
    apply_dirichlet(&mut k_mat, &mut rhs_dummy, &dofs);

    // Lanczos iteration: m steps (m = min(3*n_modes+10, n))
    let m = (3 * n_modes + 10).min(n);
    let lin = &config.solver.linear;

    let (t_alpha, t_beta, v_basis) =
        lanczos(&a_bc, &m_mat, &dofs, n, m, lin.tol, lin.max_iter, comm);

    let actual_m = t_alpha.len();
    if actual_m < n_modes {
        log::warn!(
            "Lanczos terminated early: {} steps completed, {} modes requested. \
             Invariant subspace found or PCG convergence issue.",
            actual_m, n_modes
        );
    }

    // Solve the m×m tridiagonal eigenvalue problem with nalgebra
    let (ritz_vals, ritz_vecs_small) = tridiag_eigen(&t_alpha, &t_beta);

    // Convert Ritz values μ back to λ = σ + 1/μ, then to Hz
    // Recover Ritz vectors: x_k = V * y_k  (V = Lanczos basis, y_k = k-th column of ritz_vecs_small)
    let mut eigenpairs: Vec<(f64, Vec<f64>)> = Vec::new();

    for (k, &mu) in ritz_vals.iter().enumerate().take(n_modes) {
        if mu.abs() < 1e-300 { continue; }
        let lambda = sigma + 1.0 / mu;
        if lambda <= 0.0 { continue; }
        let freq_hz = C0 * lambda.sqrt() / (2.0 * std::f64::consts::PI);

        // Ritz vector: x_k = V * y_k  where y_k is the k-th column of ritz_vecs_small
        let ritz_vec = if k < ritz_vecs_small.ncols() && !v_basis.is_empty() {
            let basis_m = v_basis.len();  // number of Lanczos vectors actually computed
            let coeff_m = ritz_vecs_small.nrows().min(basis_m);
            let mut x = vec![0.0f64; n];
            for j in 0..coeff_m {
                let y_jk = ritz_vecs_small[(j, k)];
                if y_jk.abs() < 1e-300 { continue; }
                let vj = &v_basis[j];
                for i in 0..n {
                    x[i] += y_jk * vj[i];
                }
            }
            // M-normalize the Ritz vector: ||x||_M = 1
            let mx: Vec<f64> = {
                let mut tmp = vec![0.0f64; n];
                m_mat.matvec(&x, &mut tmp, comm);
                tmp
            };
            let norm_sq: f64 = x.iter().zip(mx.iter()).map(|(a, b)| a * b).sum();
            if norm_sq > 1e-300 {
                let s = 1.0 / norm_sq.sqrt();
                x.iter_mut().for_each(|v| *v *= s);
            }
            x
        } else {
            vec![0.0f64; n]
        };

        eigenpairs.push((freq_hz, ritz_vec));
    }

    // Sort by frequency
    eigenpairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    // ── Dielectric loss perturbation → Q-factor ───────────────────────────────
    // For each mode φ:
    //   1/Q_dielectric = Σ_k (tanδ_k · φᵀ M_k φ) / (φᵀ M φ)
    // We assemble a loss-weighted mass matrix M_loss with ε·tanδ as the
    // weight function, then Q = (φᵀ M φ) / (φᵀ M_loss φ) for each mode.
    let any_lossy = work_mesh.domain_tags.keys()
        .any(|&tag| domain_map.get(tag).is_lossy());
    let q_factors: Option<Vec<f64>> = if any_lossy {
        let loss_fn = |tag: u32| {
            let mat = domain_map.get(tag);
            mat.epsilon_abs() * mat.loss_tangent   // ε₀ εᵣ tanδ
        };
        let m_loss_triplet = assemble_mass::assemble_mass(work_mesh, loss_fn)?;
        let m_loss = m_loss_triplet.to_csr();

        let qs: Vec<f64> = eigenpairs.iter().map(|(_, phi)| {
            // numerator: φᵀ M_loss φ
            let mut m_loss_phi = vec![0.0f64; n];
            for i in 0..m_loss.nrows {
                for ptr in m_loss.row_ptr[i]..m_loss.row_ptr[i + 1] {
                    m_loss_phi[i] += m_loss.values[ptr] * phi[m_loss.col_idx[ptr]];
                }
            }
            let numerator: f64 = phi.iter().zip(m_loss_phi.iter()).map(|(a, b)| a * b).sum();

            // denominator: φᵀ M φ
            let mut m_phi = vec![0.0f64; n];
            for i in 0..m_mat.nrows {
                for ptr in m_mat.row_ptr[i]..m_mat.row_ptr[i + 1] {
                    m_phi[i] += m_mat.values[ptr] * phi[m_mat.col_idx[ptr]];
                }
            }
            let denominator: f64 = phi.iter().zip(m_phi.iter()).map(|(a, b)| a * b).sum();

            if denominator.abs() > 1e-300 && numerator.abs() > 0.0 {
                denominator / numerator  // Q = 1 / tan_δ_eff
            } else {
                f64::INFINITY
            }
        }).collect();

        Some(qs)
    } else {
        None
    };

    // ── Conductor surface loss → Q_conductor ─────────────────────────────────
    // Uses perturbation: 1/Q_c = R_s · ∫_PEC |∇φ_t|² dS / (ω₀μ₀ · 2·∫_Ω ε|φ|² dΩ)
    // R_s = √(ω₀μ₀/(2σ_wall))  (surface resistance of conductor)
    //
    // Surface integral approximation: for each PEC boundary triangle (face),
    // compute area-weighted |∇φ|² using the element gradient of the adjacent volume element.
    // The denominator uses the eigenmode energy ∫ ε|∇φ|² dΩ (from stiffness Rayleigh quotient).
    let q_factors: Option<Vec<f64>> = {
        let sigma_wall = config.solver.eigenmode.as_ref()
            .map(|e| e.wall_conductivity)
            .unwrap_or(0.0);
        if sigma_wall > 0.0 {
            let freqs_for_q: Vec<f64> = eigenpairs.iter().map(|(f, _)| *f).collect();
            let qs_combined: Vec<f64> = eigenpairs.iter().zip(freqs_for_q.iter()).map(|((freq_hz, phi), _)| {
                // Rayleigh quotient denominator: φᵀ M φ (energy normalisation)
                let mut m_phi = vec![0.0f64; n];
                m_mat.matvec(phi, &mut m_phi, comm);
                let denom: f64 = phi.iter().zip(m_phi.iter()).map(|(a, b)| a * b).sum();

                // Surface integral: iterate over PEC boundary elements
                let omega = 2.0 * std::f64::consts::PI * freq_hz;
                let r_s = if omega > 0.0 {
                    (omega * 1.25663706212e-6 / (2.0 * sigma_wall)).sqrt()  // μ₀ = 4πe-7
                } else { 0.0 };

                let mut surf_integral = 0.0f64;
                for belem in &mesh.boundary_elements {
                    let bc = match mesh.boundary_tags.get(&belem.tag) {
                        Some(b) => b,
                        None => continue,
                    };
                    use rem_mesh::BoundaryTag;
                    if !matches!(bc, BoundaryTag::Pec | BoundaryTag::Ground) { continue; }
                    // For Line2 (2-D) and Tri3 (3-D surface) elements
                    let node_ids = &belem.node_ids;
                    let (grad_phi, area) = boundary_element_grad_and_area(mesh, node_ids, phi);
                    let grad_sq = grad_phi[0]*grad_phi[0] + grad_phi[1]*grad_phi[1] + grad_phi[2]*grad_phi[2];
                    surf_integral += grad_sq * area;
                }

                let q_diel = q_factors.as_ref().and_then(|qs| qs.get(eigenpairs.iter().position(|(f, _)| f == freq_hz).unwrap_or(0)).copied()).unwrap_or(f64::INFINITY);

                // Q_c from surface loss
                let q_c = if surf_integral > 1e-300 && denom > 1e-300 && r_s > 0.0 {
                    let mu0 = 1.25663706212e-6;
                    // 1/Q_c = (R_s * surf_integral) / (ω₀ μ₀ * denom)
                    let inv_qc = (r_s * surf_integral) / (omega * mu0 * denom);
                    if inv_qc > 1e-300 { 1.0 / inv_qc } else { f64::INFINITY }
                } else {
                    f64::INFINITY
                };

                // Combined: 1/Q_total = 1/Q_diel + 1/Q_c
                let inv_q_diel = if q_diel.is_infinite() { 0.0 } else { 1.0 / q_diel };
                let inv_qc = if q_c.is_infinite() { 0.0 } else { 1.0 / q_c };
                let inv_total = inv_q_diel + inv_qc;
                if inv_total > 1e-300 { 1.0 / inv_total } else { f64::INFINITY }
            }).collect();
            log::info!(
                "Q-factors: dielectric + conductor losses (σ_wall={:.3e} S/m, R_s={:.4} mΩ/□ at {:.3} GHz)",
                sigma_wall,
                {
                    let f0 = eigenpairs.first().map(|(f, _)| *f).unwrap_or(1e9);
                    let omega0 = 2.0 * std::f64::consts::PI * f0;
                    (omega0 * 1.25663706212e-6 / (2.0 * sigma_wall)).sqrt() * 1e3
                },
                eigenpairs.first().map(|(f, _)| *f / 1e9).unwrap_or(0.0)
            );
            Some(qs_combined)
        } else {
            if q_factors.is_some() {
                log::info!("Q-factors: dielectric loss only — set Solver.Eigenmode.WallConductivity for conductor losses.");
            }
            q_factors
        }
    };

    Ok(EigenResult {
        frequencies_hz: eigenpairs.iter().map(|(f, _)| *f).collect(),
        eigenvectors:   {
            let vecs: Vec<Vec<f64>> = eigenpairs.into_iter().map(|(_, mut v)| {
                if !periodic_pairs.is_empty() {
                    rem_electrostatic::bc::propagate_periodic(&mut v, &periodic_pairs);
                }
                v
            }).collect();
            vecs
        },
        q_factors,
        is_hcurl: false,
    })
}

// ---------------------------------------------------------------------------
// Shifted matrix A = K − σ M (CSR)
// ---------------------------------------------------------------------------

pub(crate) fn shifted_matrix(k: &CsrMatrix, m: &CsrMatrix, sigma: f64, n: usize) -> CsrMatrix {
    if sigma.abs() < 1e-300 {
        return k.clone();
    }
    // Compute K - σ M  in triplet form
    let mut t = TripletMatrix::with_capacity(n, n, k.nnz() + m.nnz());
    // Add all K entries
    for i in 0..k.nrows {
        for ptr in k.row_ptr[i]..k.row_ptr[i + 1] {
            t.add(i, k.col_idx[ptr], k.values[ptr]);
        }
    }
    // Subtract σ * M entries
    for i in 0..m.nrows {
        for ptr in m.row_ptr[i]..m.row_ptr[i + 1] {
            t.add(i, m.col_idx[ptr], -sigma * m.values[ptr]);
        }
    }
    t.to_csr()
}

/// Compute the average gradient of φ over a boundary element (Line2 or Tri3)
/// and the element's area (length for 1-D, area for 2-D surface).
///
/// Returns ([gx, gy, gz], measure) where the gradient is the arithmetic mean
/// of nodal gradients (first-order approximation for P1 elements).
fn boundary_element_grad_and_area(
    mesh: &RemMesh,
    node_ids: &[usize],
    phi: &[f64],
) -> ([f64; 3], f64) {
    let nodes: Vec<_> = node_ids.iter()
        .filter_map(|&id| mesh.nodes.get(id))
        .collect();
    let n = nodes.len();
    if n == 0 { return ([0.0; 3], 0.0); }

    match n {
        2 => {
            // Line2: length element (2-D mesh boundary)
            let dx = nodes[1].x - nodes[0].x;
            let dy = nodes[1].y - nodes[0].y;
            let len = (dx*dx + dy*dy).sqrt().max(1e-300);
            // Tangential direction: t = (dx, dy)/len
            // Normal: n = (-dy, dx)/len
            // Gradient of φ along tangent: (φ₁ - φ₀)/len
            let dphi = (phi[node_ids[1]] - phi[node_ids[0]]) / len;
            let gx = dphi * dx / len;
            let gy = dphi * dy / len;
            ([gx, gy, 0.0], len)
        }
        3 => {
            // Tri3: triangle surface element (3-D mesh PEC face)
            let v1 = [nodes[1].x - nodes[0].x, nodes[1].y - nodes[0].y, nodes[1].z - nodes[0].z];
            let v2 = [nodes[2].x - nodes[0].x, nodes[2].y - nodes[0].y, nodes[2].z - nodes[0].z];
            // Area = 0.5 |v1 × v2|
            let cx = v1[1]*v2[2] - v1[2]*v2[1];
            let cy = v1[2]*v2[0] - v1[0]*v2[2];
            let cz = v1[0]*v2[1] - v1[1]*v2[0];
            let area = 0.5 * (cx*cx + cy*cy + cz*cz).sqrt();
            // Gradient of P1 field on triangle:
            // ∇φ = (φ₁-φ₀)(∇λ₁) + (φ₂-φ₀)(∇λ₂)
            // For simplicity: average nodal gradient using finite differences
            let phi0 = phi.get(node_ids[0]).copied().unwrap_or(0.0);
            let phi1 = phi.get(node_ids[1]).copied().unwrap_or(0.0);
            let phi2 = phi.get(node_ids[2]).copied().unwrap_or(0.0);
            let a11 = v1[0]; let a12 = v1[1];
            let a21 = v2[0]; let a22 = v2[1];
            // Solve 2×2 system for (s, t) such that grad_φ ≈ (phi1-phi0)*grad_s + (phi2-phi0)*grad_t
            // Simplified: use barycentric gradient directly
            let det = a11*a22 - a12*a21;
            if det.abs() < 1e-300 {
                return ([0.0; 3], area);
            }
            let gx = ((phi1-phi0)*a22 - (phi2-phi0)*a21) / det;
            let gy = ((phi2-phi0)*a11 - (phi1-phi0)*a12) / det;
            ([gx, gy, 0.0], area)
        }
        _ => ([0.0; 3], 0.0),
    }
}
// ---------------------------------------------------------------------------
// Lanczos with shift-invert
// ---------------------------------------------------------------------------
//
// At each step solve: A^{-1} M v = w   (inner PCG solve)
// Builds the m-column Lanczos basis V and tridiagonal T = V^T M^{-1} A^{-1} V.

pub(crate) fn lanczos(
    a: &CsrMatrix,
    m: &CsrMatrix,
    dofs: &std::collections::HashMap<usize, f64>,
    n: usize,
    m_steps: usize,
    pcg_tol: f64,
    pcg_max_iter: usize,
    comm: &dyn Comm,
) -> (Vec<f64>, Vec<f64>, Vec<Vec<f64>>) {
    let mut alpha = Vec::with_capacity(m_steps);
    let mut beta: Vec<f64>  = Vec::with_capacity(m_steps.saturating_sub(1));
    let mut basis: Vec<Vec<f64>> = Vec::with_capacity(m_steps);

    // Initial vector: random-ish, zero on Dirichlet DOFs
    let mut v = initial_vector(n, dofs);
    m_normalize(&mut v, m, comm);

    basis.push(v.clone());

    let mut v_prev = vec![0.0f64; n];

    for j in 0..m_steps {
        // w = A^{-1} M v_j  (shift-invert apply)
        let mv = matvec_csr(m, &basis[j], comm);
        let pcg_result = solve_pcg(a, &mv, pcg_tol, pcg_max_iter, comm);
        let mut w = pcg_result.solution;

        // Zero Dirichlet DOFs in w
        for &d in dofs.keys() { if d < n { w[d] = 0.0; } }

        // Full reorthogonalization (double pass) against all previous basis vectors.
        // Prevents accumulation of rounding errors that cause spurious eigenvalues.
        for _pass in 0..2 {
            for vk in &basis {
                let mvk = matvec_csr(m, vk, comm);
                let c = comm.allreduce_f64(dot(&w, &mvk));
                axpy(-c, vk, &mut w);
                for &d in dofs.keys() { if d < n { w[d] = 0.0; } }
            }
        }

        // α_j = <v_j, w>_M  = v_j^T M w
        let mw = matvec_csr(m, &w, comm);
        let alpha_j = comm.allreduce_f64(dot(&basis[j], &mw));
        alpha.push(alpha_j);

        // w = w - α_j v_j - β_{j-1} v_{j-1}
        axpy(-alpha_j, &basis[j], &mut w);
        if j > 0 {
            axpy(-beta[j - 1], &v_prev, &mut w);
        }

        // β_j = ||w||_M
        let beta_j = m_norm(&w, m, comm);
        if j + 1 < m_steps {
            beta.push(beta_j);
            if beta_j < 1e-14 {
                // Invariant subspace found
                break;
            }
            let mut v_next = w;
            scale(1.0 / beta_j, &mut v_next);
            v_prev = basis[j].clone();
            basis.push(v_next);
        }
    }

    (alpha, beta, basis)
}

// ---------------------------------------------------------------------------
// Tridiagonal eigenvalue solve via nalgebra
// ---------------------------------------------------------------------------

/// Returns (eigenvalues, eigenvectors) of the m×m symmetric tridiagonal matrix.
/// Eigenvalues sorted largest-first (closest to shift target).
/// Eigenvectors are columns of the returned matrix (column-major, m×m).
pub(crate) fn tridiag_eigen(alpha: &[f64], beta: &[f64]) -> (Vec<f64>, DMatrix<f64>) {
    let m = alpha.len();
    if m == 0 {
        return (vec![], DMatrix::zeros(0, 0));
    }
    let mut mat = DMatrix::<f64>::zeros(m, m);
    for i in 0..m {
        mat[(i, i)] = alpha[i];
    }
    for i in 0..beta.len().min(m - 1) {
        mat[(i, i + 1)] = beta[i];
        mat[(i + 1, i)] = beta[i];
    }
    let sym = nalgebra::SymmetricEigen::new(mat);
    // Sort by eigenvalue descending (largest μ → smallest λ nearest target)
    let mut idx: Vec<usize> = (0..m).collect();
    let evals: Vec<f64> = sym.eigenvalues.iter().copied().collect();
    idx.sort_by(|&a, &b| evals[b].partial_cmp(&evals[a]).unwrap_or(std::cmp::Ordering::Equal));

    let sorted_vals: Vec<f64> = idx.iter().map(|&i| evals[i]).collect();
    // Build sorted eigenvector matrix: each column is the eigenvector for sorted_vals[k]
    let evecs = &sym.eigenvectors;
    let mut sorted_vecs = DMatrix::<f64>::zeros(m, m);
    for (new_col, &old_col) in idx.iter().enumerate() {
        for row in 0..m {
            sorted_vecs[(row, new_col)] = evecs[(row, old_col)];
        }
    }
    (sorted_vals, sorted_vecs)
}

// ---------------------------------------------------------------------------
// Small linear algebra helpers (WASM-safe, no external BLAS)
// ---------------------------------------------------------------------------

fn initial_vector(n: usize, dofs: &std::collections::HashMap<usize, f64>) -> Vec<f64> {
    let mut v = vec![1.0f64; n];
    // Pseudo-random via simple hash, so it's deterministic
    for i in 0..n {
        v[i] = ((i.wrapping_mul(2654435761)) as f64 / u32::MAX as f64) - 0.5;
    }
    for &d in dofs.keys() { if d < n { v[d] = 0.0; } }
    v
}

fn matvec_csr(m: &CsrMatrix, x: &[f64], comm: &dyn Comm) -> Vec<f64> {
    let mut y = vec![0.0f64; x.len()];
    m.matvec(x, &mut y, comm);
    y
}

fn m_inner(v: &[f64], mw: &[f64]) -> f64 {
    dot(v, mw)
}

fn m_norm(v: &[f64], m: &CsrMatrix, comm: &dyn Comm) -> f64 {
    let mv = matvec_csr(m, v, comm);
    comm.allreduce_f64(m_inner(v, &mv)).sqrt()
}

fn m_normalize(v: &mut [f64], m: &CsrMatrix, comm: &dyn Comm) {
    let n = m_norm(v, m, comm);
    if n > 1e-300 { scale(1.0 / n, v); }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(&ai, &bi)| ai * bi).sum()
}

fn axpy(alpha: f64, x: &[f64], y: &mut [f64]) {
    for (yi, &xi) in y.iter_mut().zip(x.iter()) {
        *yi += alpha * xi;
    }
}

fn scale(s: f64, v: &mut [f64]) {
    for vi in v.iter_mut() { *vi *= s; }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tridiag_eigenvalues_2x2() {
        // [[2, -1], [-1, 2]] → eigenvalues 1 and 3
        let alpha = vec![2.0, 2.0];
        let beta  = vec![-1.0];
        let (vals, vecs) = tridiag_eigen(&alpha, &beta);
        assert_eq!(vals.len(), 2);
        // Sorted largest first: [3.0, 1.0]
        assert!((vals[0] - 3.0).abs() < 1e-10, "val0={}", vals[0]);
        assert!((vals[1] - 1.0).abs() < 1e-10, "val1={}", vals[1]);
        // Eigenvectors should be orthonormal: y_0^T y_1 ≈ 0
        let dot = vecs[(0,0)]*vecs[(0,1)] + vecs[(1,0)]*vecs[(1,1)];
        assert!(dot.abs() < 1e-10, "eigenvecs not orthogonal: dot={}", dot);
    }
}
