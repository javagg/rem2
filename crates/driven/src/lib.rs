//! Driven (frequency-domain) solver — Phase 7 (v0.4)
//!
//! Solves the frequency-domain scalar wave equation:
//!   −∇·(ε ∇φ) − k² ε φ = J_port      (k = ω/c)
//!
//! For each excitation frequency ω = 2πf in [MinFreq, MaxFreq] with step FreqStep:
//!   1. Assemble K (stiffness) and M (mass, consistent or lumped)
//!   2. Build system A = K − k² M  (complex: real part K−k²M, zero imaginary part)
//!   3. Apply lumped-port / Dirichlet BCs (WavePort: modal φ=mode_shape × V)
//!   4. Solve with complex GMRES (nalgebra); correct at and above resonance
//!   5. Compute port impedance Z, reflection S₁₁ (WavePort: Z_TE = ωμ₀/k_z)
//!   6. Write CSV and (optionally) VTK per save_step
//!
//! WavePort modal analysis (v0.4):
//!   - Extract port cross-section 1-D mesh (Line2 elements on WavePort tag).
//!   - Solve 1-D Laplacian eigenvalue K_p x = λ M_p x with Dirichlet endpoints.
//!   - First eigenvalue λ₁ = k_c² (cutoff wavenumber²).
//!   - Use eigenvector as mode-shape Dirichlet profile on the port nodes.
//!   - Modal impedance Z_TE = ωμ₀/k_z, k_z = √(k²−k_c²).
//!   - Falls back to TEM (uniform φ=V, Z₀=50Ω) when below cutoff or on solve failure.
//!
//! v0.4 changes:
//!   - Replace real PCG with complex GMRES (nalgebra DMatrix<Complex64>).
//!   - S11 now carries both real and imaginary parts; |S11| and phase are correct.

pub mod output;
pub mod port_modal;

use nalgebra::{DMatrix, DVector};
use num_complex::Complex64;
use rem_config::{PalaceConfig, CurrentDipoleSpec};
use rem_core::{CsrMatrix, RemError, RemResult, TripletMatrix};
use rem_eigenmode::assemble_mass::assemble_mass;
use rem_electrostatic::{assemble::assemble_stiffness, bc::{collect_dirichlet_dofs, apply_dirichlet}};
use rem_materials::DomainMap;
use rem_mesh::{RemMesh, BoundaryTag, amr};
use rem_mesh::gmsh::read_msh_file;
use rem_parallel::Comm;
use port_modal::{PortMode, compute_wave_port_mode};
use std::collections::HashMap;
use std::path::Path;

const C0: f64 = 2.997_924_58e8;

/// Per-frequency S-parameter result.
pub struct FreqResult {
    pub freq_hz: f64,
    pub s11_re:  f64,
    pub s11_im:  f64,
}

/// Result of a driven frequency sweep.
pub struct DrivenResult {
    /// S-parameter results for each frequency point.
    pub freq_results: Vec<FreqResult>,
    /// Real part of the nodal potential at the frequency of peak |S11| response
    /// (i.e., worst reflection, which often corresponds to near-resonance).
    /// Empty if the sweep produced no frequencies.
    pub peak_phi: Vec<f64>,
    /// Frequency [Hz] at which `peak_phi` was recorded.
    pub peak_freq_hz: f64,
}

/// Entry point called from rem-cli.
pub fn run(config: &PalaceConfig, comm: &dyn Comm) -> RemResult<()> {
    log::info!("=== Driven (frequency-domain) solver ===");

    let mesh_path = Path::new(&config.model.mesh);
    let raw = read_msh_file(mesh_path)?;
    let mut mesh = RemMesh::from_raw(raw, config)?;
    mesh.set_comm(comm.rank(), comm.size());

    run_with_mesh(config, &mesh, comm).map(|_| ())
}

/// Entry point for pre-loaded mesh (used by WASM path).
/// Returns the driven frequency sweep result including S-params and peak E-field.
pub fn run_with_mesh(config: &PalaceConfig, mesh: &RemMesh, comm: &dyn Comm) -> RemResult<DrivenResult> {
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
        let (excited_port, _port_kind) = find_excited_port(mesh);

        let eps_fn = |tag: u32| domain_map.get(tag).epsilon_abs();
        let mut cur_mesh = mesh.clone();
        for amr_iter in 1..=max_amr_iter {
            let k_mat = assemble_stiffness(&cur_mesh, eps_fn)?.to_csr();
            let m_mat = assemble_mass(&cur_mesh, eps_fn)?.to_csr();
            // Use real system for AMR estimator (estimator only needs shape, not phase)
            let a_mat = shifted_matrix(&k_mat, &m_mat, k2, cur_mesh.n_nodes());
            let mut a_bc = a_mat;
            let dofs = collect_dirichlet_dofs(&cur_mesh, excited_port, 1.0);
            let mut rhs = vec![0.0f64; cur_mesh.n_nodes()];
            apply_dirichlet(&mut a_bc, &mut rhs, &dofs);

            // Use real solve for AMR only (estimator doesn't need complex accuracy)
            use rem_core::solve_pcg;
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

/// Build the list of frequencies to sweep from `DrivenSolver` config.
///
/// Precedence:
///   1. If `Samples` is non-empty, expand each sample spec (Linear/Log/Explicit).
///   2. Otherwise use MinFreq/MaxFreq/FreqStep.
fn build_freq_list(drv_cfg: &rem_config::DrivenSolver) -> RemResult<Vec<f64>> {
    if !drv_cfg.samples.is_empty() {
        let mut freqs: Vec<f64> = Vec::new();
        for spec in &drv_cfg.samples {
            match spec.sample_type.as_str() {
                "Linear" | "" => {
                    if spec.freq_step <= 0.0 || spec.min_freq > spec.max_freq {
                        return Err(RemError::Config(format!(
                            "Driven.Samples Linear: invalid MinFreq={}, MaxFreq={}, FreqStep={}",
                            spec.min_freq, spec.max_freq, spec.freq_step
                        )));
                    }
                    let n = ((spec.max_freq - spec.min_freq) / spec.freq_step).ceil() as usize + 1;
                    for i in 0..n {
                        let f = spec.min_freq + i as f64 * spec.freq_step;
                        if f <= spec.max_freq + spec.freq_step * 0.5 { freqs.push(f); }
                    }
                }
                "Log" => {
                    if spec.freq_step < 2.0 || spec.min_freq <= 0.0 || spec.max_freq <= 0.0 {
                        return Err(RemError::Config(
                            "Driven.Samples Log: FreqStep is points-per-decade; MinFreq/MaxFreq must be > 0".into()
                        ));
                    }
                    let n_decades = (spec.max_freq / spec.min_freq).log10();
                    let n_pts = (n_decades * spec.freq_step).ceil() as usize + 1;
                    for i in 0..n_pts {
                        let f = spec.min_freq * 10.0_f64.powf(i as f64 / spec.freq_step);
                        if f <= spec.max_freq * 1.001 { freqs.push(f); }
                    }
                }
                "Point" | "Explicit" => {
                    freqs.extend_from_slice(&spec.freq);
                }
                other => {
                    log::warn!("Driven.Samples: unknown Type={:?}; skipping this spec", other);
                }
            }
        }
        if freqs.is_empty() {
            return Err(RemError::Config("Driven.Samples produced no frequency points".into()));
        }
        freqs.sort_by(|a, b| a.partial_cmp(b).unwrap());
        freqs.dedup_by(|a, b| (*b - *a).abs() < 1.0); // remove duplicates closer than 1 Hz
        log::info!("Driven: {} frequency points from Samples array", freqs.len());
        Ok(freqs)
    } else {
        let f_min  = drv_cfg.min_freq;
        let f_max  = drv_cfg.max_freq;
        let f_step = drv_cfg.freq_step;
        if f_step <= 0.0 || f_min > f_max {
            return Err(RemError::Config(
                "Driven: FreqStep must be > 0 and MinFreq ≤ MaxFreq".into()
            ));
        }
        let n = ((f_max - f_min) / f_step).ceil() as usize + 1;
        let freqs: Vec<f64> = (0..n)
            .map(|i| f_min + i as f64 * f_step)
            .take_while(|&f| f <= f_max + f_step * 0.5)
            .collect();
        Ok(freqs)
    }
}

fn run_frequency_sweep(
    config: &PalaceConfig,
    drv_cfg: &rem_config::DrivenSolver,
    mesh: &RemMesh,
    domain_map: &DomainMap,
    comm: &dyn Comm,
) -> RemResult<DrivenResult> {
    let freqs = build_freq_list(drv_cfg)?;
    let save_step = drv_cfg.save_step.max(1);

    // Assemble K and M once (frequency-independent, real part)
    let eps_fn     = |tag: u32| domain_map.get(tag).epsilon_abs();
    let k_mat = assemble_stiffness(mesh, eps_fn)?.to_csr();
    let m_mat = assemble_mass(mesh, eps_fn)?.to_csr();

    // Pre-assemble loss matrices if any domain has dielectric loss or conductivity.
    // K_loss uses ε₀·εᵣ·tanδ;  M_loss is the same weight (for tanδ part).
    // K_cond/M_cond use σ [S/m]; conductivity adds −j·σ/ω to ε_eff (per-step).
    let any_lossy = mesh.domain_tags.keys()
        .any(|&tag| domain_map.get(tag).is_lossy());
    let (k_loss_dense, m_loss_dense, k_cond_dense, m_cond_dense) = if any_lossy {
        let tan_fn  = |tag: u32| { let m = domain_map.get(tag); m.epsilon_abs() * m.loss_tangent };
        let cond_fn = |tag: u32| { domain_map.get(tag).conductivity };   // σ [S/m]
        let n = mesh.n_nodes();
        let kl = assemble_stiffness(mesh, tan_fn)?.to_csr();
        let ml = assemble_mass(mesh,     tan_fn)?.to_csr();
        let kc = assemble_stiffness(mesh, cond_fn)?.to_csr();
        let mc = assemble_mass(mesh,     cond_fn)?.to_csr();
        (Some(csr_to_complex_dense(&kl, n)),
         Some(csr_to_complex_dense(&ml, n)),
         Some(csr_to_complex_dense(&kc, n)),
         Some(csr_to_complex_dense(&mc, n)))
    } else {
        (None, None, None, None)
    };
    if any_lossy {
        log::info!("Lossy materials detected: complex ε assembly enabled");
    }

    let out_dir = config.problem.output_dir();
    #[cfg(not(target_arch = "wasm32"))]
    std::fs::create_dir_all(out_dir).map_err(RemError::Io)?;

    let lin = &config.solver.linear;
    let mut freq_results: Vec<FreqResult> = Vec::with_capacity(freqs.len());
    // Track the phi at the frequency with maximum |S11| (peak reflection / near-resonance)
    let mut peak_phi: Vec<f64> = Vec::new();
    let mut peak_freq_hz: f64 = 0.0;
    let mut peak_s11_mag: f64 = -1.0;
    let (excited_port, port_kind) = find_excited_port(mesh);

    // Pre-compute wave-port mode shape (frequency-independent, geometry only)
    let wave_port_mode: Option<PortMode> = if let PortKind::Wave(idx) = port_kind {
        match compute_wave_port_mode(mesh, idx) {
            Some(m) => {
                log::info!(
                    "WavePort {idx}: TE modal excitation enabled (k_c={:.4e} rad/m)",
                    m.kc
                );
                Some(m)
            }
            None => {
                log::warn!(
                    "WavePort {idx}: modal solve failed; falling back to TEM (φ=V uniform)"
                );
                None
            }
        }
    } else {
        None
    };

    // Pre-convert real K and M to complex dense matrices for GMRES
    let n = mesh.n_nodes();
    let k_dense = csr_to_complex_dense(&k_mat, n);
    let m_dense = csr_to_complex_dense(&m_mat, n);

    for (step, &freq) in freqs.iter().enumerate() {
        let omega = 2.0 * std::f64::consts::PI * freq;
        let k_wave = omega / C0;
        let k2 = k_wave * k_wave;

        // Build complex system matrix:
        //   A = (K_re + j·K_loss) − k²·(M_re + j·M_loss) − j·(σ/ω)·M_cond
        // where the imaginary parts encode dielectric loss (tanδ) and conductivity.
        let mut a = k_dense.clone();
        for i in 0..n {
            for j in 0..n {
                a[(i, j)] -= Complex64::new(k2, 0.0) * m_dense[(i, j)];
            }
        }
        if let (Some(kl), Some(ml), Some(kc), Some(mc)) =
            (&k_loss_dense, &m_loss_dense, &k_cond_dense, &m_cond_dense)
        {
            for i in 0..n {
                for j in 0..n {
                    // tanδ part: A += j·(K_loss − k²·M_loss)
                    let tan_contrib = Complex64::new(0.0, 1.0)
                        * (kl[(i,j)] - Complex64::new(k2, 0.0) * ml[(i,j)]);
                    // conductivity part: A += (−j/ω)·(K_cond − k²·M_cond)
                    let cond_contrib = Complex64::new(0.0, -1.0 / omega)
                        * (kc[(i,j)] - Complex64::new(k2, 0.0) * mc[(i,j)]);
                    a[(i, j)] += tan_contrib + cond_contrib;
                }
            }
        }

        // Build Dirichlet DOF map
        let dofs: HashMap<usize, f64> = if let Some(mode) = &wave_port_mode {
            if mode.is_propagating(freq) {
                collect_dirichlet_dofs_modal(mesh, excited_port, mode)
            } else {
                log::warn!("f={freq:.3e} Hz is below WavePort cutoff (evanescent); using φ=0");
                collect_dirichlet_dofs(mesh, excited_port, 0.0)
            }
        } else {
            collect_dirichlet_dofs(mesh, excited_port, 1.0)
        };

        let mut rhs_c = vec![Complex64::ZERO; n];
        apply_dirichlet_complex(&mut a, &mut rhs_c, &dofs);

        // Add CurrentDipole source contributions to RHS
        //
        // For each Hertzian dipole at position r₀ with moment I·L [A·m] and direction d̂:
        //   f_i += jω μ₀ · Moment · |d̂|  at the nearest mesh node to r₀.
        //
        // In the scalar wave equation −∇·ε∇φ − k²εφ = S, the current source
        // maps to S_i = jω μ₀ J at DOF i.  We lump the point source onto the
        // single nearest node (zeroth-order approximation; adequate for far-field).
        if !config.domains.current_dipole.is_empty() {
            let jw_mu0 = Complex64::new(0.0, 2.0 * std::f64::consts::PI * freq)
                * Complex64::new(4.0 * std::f64::consts::PI * 1.0e-7, 0.0);
            for dipole in &config.domains.current_dipole {
                let node = nearest_node(mesh, dipole);
                if !dofs.contains_key(&node) {
                    let dir_mag = (dipole.direction.iter().map(|x| x * x).sum::<f64>()).sqrt();
                    let mag = if dir_mag > 1e-300 { dir_mag } else { 1.0 };
                    rhs_c[node] += jw_mu0 * Complex64::new(dipole.moment * mag, 0.0);
                }
            }
        }

        // Solve with complex GMRES
        let phi_c = gmres_complex(&a, &rhs_c, lin.tol, lin.max_iter)?;

        let phi_re: Vec<f64> = phi_c.iter().map(|x| x.re).collect();

        let (v_port, i_port) = compute_port_vi_complex(mesh, &phi_c, &k_dense, excited_port);

        // Reference impedance
        let z0 = if let Some(mode) = &wave_port_mode {
            let z = mode.te_impedance(freq);
            if z.is_finite() { z } else { 50.0 }
        } else {
            lumped_port_resistance(mesh, excited_port)
        };

        let s11 = if i_port.norm() > 1e-300 {
            let z = v_port / i_port;
            let z0c = Complex64::new(z0, 0.0);
            (z - z0c) / (z + z0c)
        } else {
            Complex64::ZERO
        };

        log::info!(
            "f={:.3e} Hz  |S11|={:.4}  ∠S11={:.2}°  Z0={:.1}Ω",
            freq, s11.norm(), s11.arg().to_degrees(), z0
        );

        freq_results.push(FreqResult { freq_hz: freq, s11_re: s11.re, s11_im: s11.im });

        // Track phi at the frequency with maximum |S11| (peak reflection)
        let s11_mag = s11.norm();
        if s11_mag > peak_s11_mag {
            peak_s11_mag = s11_mag;
            peak_freq_hz = freq;
            peak_phi = phi_re.clone();
        }

        #[cfg(not(target_arch = "wasm32"))]
        if step % save_step == 0 {
            output::write_field_vtk(out_dir, mesh, &phi_re, step + 1)?;
        }
        let _ = comm;
    }

    #[cfg(not(target_arch = "wasm32"))]
    output::write_s_params(out_dir, &freq_results)?;
    log::info!(
        "Driven solve complete: {} frequency points, peak |S11|={:.4} at {:.3e} Hz",
        freq_results.len(), peak_s11_mag, peak_freq_hz
    );
    Ok(DrivenResult { freq_results, peak_phi, peak_freq_hz })
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Whether the excited port is a LumpedPort or a WavePort.
#[derive(Debug, Clone, Copy)]
enum PortKind {
    Lumped,
    Wave(u32),
    None,
}

/// Build real A = K − k² M (used only for AMR estimator).
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

/// Convert a real CsrMatrix to a complex dense DMatrix<Complex64>.
fn csr_to_complex_dense(mat: &CsrMatrix, n: usize) -> DMatrix<Complex64> {
    let mut d = DMatrix::<Complex64>::zeros(n, n);
    for i in 0..mat.nrows {
        for ptr in mat.row_ptr[i]..mat.row_ptr[i + 1] {
            let j = mat.col_idx[ptr];
            d[(i, j)] += Complex64::new(mat.values[ptr], 0.0);
        }
    }
    d
}

/// Apply Dirichlet BCs to a complex dense matrix and RHS.
///
/// For each constrained DOF `(row, val)`:
///   - Zero the row and set diagonal to 1
///   - Subtract val * column from RHS, then zero the column
///   - Set rhs[row] = val
fn apply_dirichlet_complex(
    a: &mut DMatrix<Complex64>,
    rhs: &mut Vec<Complex64>,
    dofs: &HashMap<usize, f64>,
) {
    let n = a.nrows();
    // Subtract column contributions from RHS
    for (&row, &val) in dofs.iter() {
        let v = Complex64::new(val, 0.0);
        for i in 0..n {
            if !dofs.contains_key(&i) {
                rhs[i] -= a[(i, row)] * v;
            }
        }
    }
    // Zero rows and columns, set diagonal to 1
    for (&row, &val) in dofs.iter() {
        for j in 0..n { a[(row, j)] = Complex64::ZERO; }
        for i in 0..n { a[(i, row)] = Complex64::ZERO; }
        a[(row, row)] = Complex64::new(1.0, 0.0);
        rhs[row] = Complex64::new(val, 0.0);
    }
}

/// Complex GMRES solver using nalgebra (restart=30).
fn gmres_complex(
    a: &DMatrix<Complex64>,
    rhs: &[Complex64],
    tol: f64,
    max_iter: usize,
) -> RemResult<Vec<Complex64>> {
    let n = rhs.len();
    const RESTART: usize = 30;
    let max_outer = max_iter / RESTART + 1;

    let b = DVector::<Complex64>::from_iterator(n, rhs.iter().copied());
    let b_norm = b.norm();
    if b_norm < f64::EPSILON {
        return Ok(vec![Complex64::ZERO; n]);
    }

    let mut x = DVector::<Complex64>::zeros(n);

    for _outer in 0..max_outer {
        // r = b - A*x
        let r = &b - a * &x;
        let beta = r.norm();
        if beta / b_norm < tol {
            return Ok(x.iter().copied().collect());
        }

        let mut v: Vec<DVector<Complex64>> = Vec::with_capacity(RESTART + 1);
        v.push(r / Complex64::new(beta, 0.0));

        let mut h = vec![vec![Complex64::ZERO; RESTART]; RESTART + 1];
        let mut g = vec![Complex64::ZERO; RESTART + 1];
        let mut c = vec![0.0f64; RESTART];
        let mut s = vec![Complex64::ZERO; RESTART];
        g[0] = Complex64::new(beta, 0.0);

        let mut j_done = RESTART;
        for j in 0..RESTART {
            let w_full = a * &v[j];
            let mut w = w_full;
            for i in 0..=j {
                h[i][j] = v[i].dotc(&w);
                let hij = h[i][j];
                w -= &v[i] * hij;
            }
            let h_next = w.norm();
            h[j + 1][j] = Complex64::new(h_next, 0.0);

            if h_next > 1e-14 {
                v.push(w / Complex64::new(h_next, 0.0));
            }

            // Apply previous Givens rotations
            for i in 0..j {
                let tmp = c[i] * h[i][j] + s[i] * h[i + 1][j];
                h[i + 1][j] = -s[i].conj() * h[i][j] + c[i] * h[i + 1][j];
                h[i][j] = tmp;
            }
            // New Givens rotation
            let rr = (h[j][j].norm_sqr() + h[j + 1][j].norm_sqr()).sqrt();
            if rr > 1e-14 {
                c[j] = h[j][j].norm() / rr;
                s[j] = h[j + 1][j] * Complex64::new(h[j][j].norm() / rr / h[j][j].norm_sqr().max(1e-300), 0.0) * h[j][j].conj();
                h[j][j] = Complex64::new(rr, 0.0);
                h[j + 1][j] = Complex64::ZERO;
                let g_next = -s[j].conj() * g[j];
                g[j] = Complex64::new(c[j], 0.0) * g[j];
                g[j + 1] = g_next;
            }

            if g[j + 1].norm() / b_norm < tol {
                j_done = j + 1;
                break;
            }
        }

        // Back-substitution for y in H*y = g
        let m = j_done.min(RESTART);
        let mut y = vec![Complex64::ZERO; m];
        for i in (0..m).rev() {
            y[i] = g[i];
            for k in (i + 1)..m {
                let yk = y[k];
                y[i] -= h[i][k] * yk;
            }
            if h[i][i].norm() > 1e-300 {
                y[i] /= h[i][i];
            }
        }

        // Update x
        for j in 0..m {
            let yj = y[j];
            x += &v[j] * yj;
        }
    }

    // Return best solution even if not converged
    Ok(x.iter().copied().collect())
}

/// Find the first excited port and report its kind.
fn find_excited_port(mesh: &RemMesh) -> (Option<u32>, PortKind) {
    for bc in mesh.boundary_tags.values() {
        match bc {
            BoundaryTag::LumpedPort { index, .. } => return (Some(*index), PortKind::Lumped),
            BoundaryTag::WavePort { index } => return (Some(*index), PortKind::Wave(*index)),
            _ => {}
        }
    }
    (None, PortKind::None)
}

/// Get port resistance from config (default 50 Ω).
fn lumped_port_resistance(mesh: &RemMesh, port_idx: Option<u32>) -> f64 {
    if let Some(idx) = port_idx {
        for bc in mesh.boundary_tags.values() {
            match bc {
                BoundaryTag::LumpedPort { index, r } if *index == idx => {
                    return if *r > 0.0 { *r } else { 50.0 };
                }
                BoundaryTag::WavePort { index } if *index == idx => {
                    return 50.0;
                }
                _ => {}
            }
        }
    }
    50.0
}

/// Build Dirichlet DOF map using the modal shape as excitation profile.
fn collect_dirichlet_dofs_modal(
    mesh: &RemMesh,
    excited_index: Option<u32>,
    mode: &PortMode,
) -> HashMap<usize, f64> {
    let mut dofs: HashMap<usize, f64> = HashMap::new();

    for belem in &mesh.boundary_elements {
        if mesh.size > 1 && belem.rank != mesh.rank {
            continue;
        }
        let bc = match mesh.boundary_tags.get(&belem.tag) {
            Some(b) => b,
            None => continue,
        };

        match bc {
            BoundaryTag::Pec | BoundaryTag::Ground => {
                for &nid in &belem.node_ids {
                    dofs.entry(nid).or_insert(0.0);
                }
            }
            BoundaryTag::LumpedPort { index, .. } | BoundaryTag::Terminal { index } => {
                let val = if Some(*index) == excited_index { 1.0 } else { 0.0 };
                for &nid in &belem.node_ids {
                    dofs.entry(nid).or_insert(val);
                }
            }
            BoundaryTag::WavePort { index } => {
                if Some(*index) == excited_index {
                    for &nid in &belem.node_ids {
                        let val = mode.shape.get(&nid).copied().unwrap_or(0.0);
                        dofs.entry(nid).or_insert(val);
                    }
                } else {
                    for &nid in &belem.node_ids {
                        dofs.entry(nid).or_insert(0.0);
                    }
                }
            }
            _ => {}
        }
    }

    dofs
}

/// Compute complex port voltage and current from complex solution.
fn compute_port_vi_complex(
    mesh: &RemMesh,
    phi: &[Complex64],
    k: &DMatrix<Complex64>,
    port_idx: Option<u32>,
) -> (Complex64, Complex64) {
    let Some(idx) = port_idx else { return (Complex64::ZERO, Complex64::ZERO); };

    let port_nodes: Vec<usize> = mesh.boundary_elements.iter()
        .filter(|e| {
            match mesh.boundary_tags.get(&e.tag) {
                Some(BoundaryTag::LumpedPort { index, .. }) => *index == idx,
                Some(BoundaryTag::WavePort { index }) => *index == idx,
                _ => false,
            }
        })
        .flat_map(|e| e.node_ids.iter().copied())
        .collect();

    if port_nodes.is_empty() { return (Complex64::ZERO, Complex64::ZERO); }

    let v_port = port_nodes.iter().map(|&n| phi[n]).sum::<Complex64>()
        / Complex64::new(port_nodes.len() as f64, 0.0);

    let n = phi.len();
    let mut i_port = Complex64::ZERO;
    for &row in &port_nodes {
        let kphi_n: Complex64 = (0..n).map(|col| k[(row, col)] * phi[col]).sum();
        i_port += kphi_n;
    }

    (v_port, i_port)
}

/// Find the mesh node index nearest to the dipole center position.
fn nearest_node(mesh: &RemMesh, dipole: &CurrentDipoleSpec) -> usize {
    let cx = dipole.center.first().copied().unwrap_or(0.0);
    let cy = dipole.center.get(1).copied().unwrap_or(0.0);
    let cz = dipole.center.get(2).copied().unwrap_or(0.0);

    mesh.nodes.iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| {
            let da = (a.x - cx).powi(2) + (a.y - cy).powi(2) + (a.z - cz).powi(2);
            let db = (b.x - cx).powi(2) + (b.y - cy).powi(2) + (b.z - cz).powi(2);
            da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|(i, _)| i)
        .unwrap_or(0)
}
