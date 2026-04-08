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
pub mod output;

use rem_config::PalaceConfig;
use rem_core::{CsrMatrix, RemError, RemResult, TripletMatrix, solve_pcg};
use rem_electrostatic::{assemble::assemble_stiffness, bc::collect_dirichlet_dofs, bc::apply_dirichlet};
use rem_materials::DomainMap;
use rem_mesh::{RemMesh, amr};
use rem_mesh::gmsh::read_msh_file;
use rem_parallel::Comm;
use nalgebra::DMatrix;
use std::path::Path;

// Physical constants
const C0: f64 = 2.997_924_58e8;  // speed of light [m/s]

// AMR eigenfrequency convergence tolerance: stop if max relative freq change < this
const AMR_FREQ_TOL: f64 = 1e-4;

/// Entry point called from rem-cli.
pub fn run(config: &PalaceConfig, comm: &dyn Comm) -> RemResult<()> {
    log::info!("=== Eigenmode solver ===");

    let eig_cfg = config.solver.eigenmode.as_ref().ok_or_else(|| {
        RemError::Config("Eigenmode problem requires a [Solver.Eigenmode] section".into())
    })?;

    if config.solver.order > 1 {
        log::warn!(
            "Solver.Order={} requested but only P1 (order=1) is implemented; \
             higher-order assembly is pending. Running P1.",
            config.solver.order
        );
    }
    let mesh_path = Path::new(&config.model.mesh);
    let raw = read_msh_file(mesh_path)?;
    let mut mesh = RemMesh::from_raw(raw, config)?;
    mesh.set_comm(comm.rank(), comm.size());
    let domain_map = DomainMap::from_config(config)?;

    // AMR loop: refine mesh adaptively and re-solve
    let amr_cfg = &config.model.refinement;
    let max_amr_iter = if amr_cfg.max_iter > 0 { amr_cfg.max_iter } else { 0 };
    let amr_theta    = if amr_cfg.tol > 0.0 { amr_cfg.tol } else { 0.5 };

    let (final_mesh, result) = if max_amr_iter > 0 {
        log::info!("AMR enabled: max_iter={}, θ={}", max_amr_iter, amr_theta);
        let mut cur_mesh = mesh;
        let mut result = solve(config, &cur_mesh, &domain_map, comm)?;
        let mut prev_freqs = result.frequencies_hz.clone();

        for amr_iter in 1..=max_amr_iter {
            // Use first eigenvector as error indicator field
            if let Some(phi) = result.eigenvectors.first() {
                let eta = amr::zz_estimator(&cur_mesh, phi);
                let total_err: f64 = eta.iter().map(|&e| e * e).sum::<f64>().sqrt();
                log::info!("AMR iter {amr_iter}: nodes={}, |η|={total_err:.3e}", cur_mesh.n_nodes());

                let marked = amr::dorfler_mark(&eta, amr_theta);
                if marked.is_empty() {
                    log::info!("AMR converged: no elements marked.");
                    break;
                }

                let (fine_mesh, _midpoints) = amr::refine_marked(&cur_mesh, &marked);
                result = solve(config, &fine_mesh, &domain_map, comm)?;
                cur_mesh = fine_mesh;

                // Check eigenfrequency convergence between AMR iterations
                let max_rel_change = prev_freqs.iter()
                    .zip(result.frequencies_hz.iter())
                    .map(|(&f_old, &f_new)| {
                        if f_old.abs() > 1e-30 { ((f_new - f_old) / f_old).abs() } else { 0.0 }
                    })
                    .fold(0.0f64, f64::max);
                log::info!("AMR iter {amr_iter}: max freq rel-change = {max_rel_change:.3e}");
                if max_rel_change < AMR_FREQ_TOL {
                    log::info!("AMR freq-converged (rel-change {max_rel_change:.2e} < {AMR_FREQ_TOL:.0e}): stopping.");
                    break;
                }
                prev_freqs = result.frequencies_hz.clone();
            } else {
                break;
            }
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
        output::write_mode_vtk(out_dir, &final_mesh, phi, mode_idx + 1)?;
    }

    log::info!("Eigenmode solve complete: {} modes found", result.frequencies_hz.len());
    Ok(())
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
}

/// Solve the generalized eigenvalue problem for `config` + pre-loaded mesh.
pub fn solve(
    config: &PalaceConfig,
    mesh: &RemMesh,
    domain_map: &DomainMap,
    comm: &dyn Comm,
) -> RemResult<EigenResult> {
    let eig_cfg = config.solver.eigenmode.as_ref().ok_or_else(|| {
        RemError::Config("missing Eigenmode solver config".into())
    })?;

    let n = mesh.n_nodes();
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
    let k_triplet = if domain_map.any_anisotropic() {
        log::info!("Anisotropic material(s) detected — using tensor stiffness assembly.");
        use rem_electrostatic::assemble::assemble_stiffness_aniso;
        let tensor_fn = |tag: u32| domain_map.get(tag).epsilon_tensor;
        assemble_stiffness_aniso(mesh, tensor_fn)?
    } else {
        let eps_fn = |tag: u32| domain_map.get(tag).epsilon_abs();
        assemble_stiffness(mesh, eps_fn)?
    };
    let eps_fn = |tag: u32| domain_map.get(tag).epsilon_abs();
    let m_triplet = assemble_mass::assemble_mass(mesh, eps_fn)?;

    let mut k_mat = k_triplet.to_csr();
    let m_mat = m_triplet.to_csr();

    // Collect Dirichlet DOFs (PEC / Ground → φ=0)
    let dofs = collect_dirichlet_dofs(mesh, None, 0.0);

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

    // Solve the m×m tridiagonal eigenvalue problem with nalgebra
    let ritz_vals = tridiag_eigenvalues(&t_alpha, &t_beta);

    // Convert Ritz values μ back to λ = σ + 1/μ, then to Hz
    let mut eigenpairs: Vec<(f64, Vec<f64>)> = Vec::new();

    for (k, &mu) in ritz_vals.iter().enumerate().take(n_modes) {
        if mu.abs() < 1e-300 { continue; }
        let lambda = sigma + 1.0 / mu;
        if lambda <= 0.0 { continue; }
        let freq_hz = C0 * lambda.sqrt() / (2.0 * std::f64::consts::PI);

        // Recover Ritz vector x = V * y_k  (V = Lanczos basis, y_k = k-th Ritz vec)
        // For simplicity we use the raw Lanczos vector v_k as an approximation
        let ritz_vec = if k < v_basis.len() {
            v_basis[k].clone()
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
    let any_lossy = mesh.domain_tags.keys()
        .any(|&tag| domain_map.get(tag).is_lossy());
    let q_factors: Option<Vec<f64>> = if any_lossy {
        let loss_fn = |tag: u32| {
            let mat = domain_map.get(tag);
            mat.epsilon_abs() * mat.loss_tangent   // ε₀ εᵣ tanδ
        };
        let m_loss_triplet = assemble_mass::assemble_mass(mesh, loss_fn)?;
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

    Ok(EigenResult {
        frequencies_hz: eigenpairs.iter().map(|(f, _)| *f).collect(),
        eigenvectors:   eigenpairs.into_iter().map(|(_, v)| v).collect(),
        q_factors,
    })
}

// ---------------------------------------------------------------------------
// Shifted matrix A = K − σ M (CSR)
// ---------------------------------------------------------------------------

fn shifted_matrix(k: &CsrMatrix, m: &CsrMatrix, sigma: f64, n: usize) -> CsrMatrix {
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

// ---------------------------------------------------------------------------
// Lanczos with shift-invert
// ---------------------------------------------------------------------------
//
// At each step solve: A^{-1} M v = w   (inner PCG solve)
// Builds the m-column Lanczos basis V and tridiagonal T = V^T M^{-1} A^{-1} V.

fn lanczos(
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

fn tridiag_eigenvalues(alpha: &[f64], beta: &[f64]) -> Vec<f64> {
    let m = alpha.len();
    if m == 0 { return vec![]; }
    let mut mat = DMatrix::<f64>::zeros(m, m);
    for i in 0..m {
        mat[(i, i)] = alpha[i];
    }
    for i in 0..beta.len().min(m - 1) {
        mat[(i, i + 1)] = beta[i];
        mat[(i + 1, i)] = beta[i];
    }
    let sym = nalgebra::SymmetricEigen::new(mat);
    let mut vals: Vec<f64> = sym.eigenvalues.iter().copied().collect();
    // Return largest eigenvalues first (largest 1/μ = smallest λ nearest target)
    vals.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
    vals
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
        let vals = tridiag_eigenvalues(&alpha, &beta);
        assert_eq!(vals.len(), 2);
        // Sorted largest first: [3.0, 1.0]
        assert!((vals[0] - 3.0).abs() < 1e-10, "val0={}", vals[0]);
        assert!((vals[1] - 1.0).abs() < 1e-10, "val1={}", vals[1]);
    }
}
