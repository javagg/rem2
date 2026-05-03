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

pub mod far_field;
pub mod near_field;
pub mod output;
pub mod port_modal;
pub mod rom;
pub mod vf;

use nalgebra::{DMatrix, DVector};
use num_complex::Complex64;
use rem_config::{PalaceConfig, CurrentDipoleSpec};
use rem_core::{CsrMatrix, RemError, RemResult, TripletMatrix, report_peak_memory, solve_pcg_complex, CsrMatrixComplex};
use rem_eigenmode::assemble_mass::assemble_mass;
use rem_electrostatic::{assemble::{assemble_stiffness, assemble_stiffness_aniso}, bc::{collect_dirichlet_dofs, collect_dirichlet_dofs_open_circuit, apply_dirichlet}, postprocess};
use rem_materials::DomainMap;
use rem_mesh::{RemMesh, BoundaryTag, ElementKind, FemSubMesh2d, amr, extract_submesh_tri3};
use rem_mesh::gmsh::read_msh_file;
use rem_parallel::Comm;
use port_modal::{PortMode, PortSupportRegionSummary, collect_port_support_region, compute_wave_port_mode};
use std::collections::HashMap;
use std::path::Path;

const C0: f64 = 2.997_924_58e8;

/// Per-frequency S-parameter result.
pub struct FreqResult {
    pub freq_hz: f64,
    /// S11 kept for backward compatibility.
    pub s11_re:  f64,
    pub s11_im:  f64,
    /// Full N×N S-matrix (row i, col j → S[i][j]), indexed by port order in `port_list`.
    /// Empty when only one port is present (backward-compat path skips this).
    pub s_matrix: Vec<Vec<Complex64>>,
    /// Ordered port indices matching rows/cols of `s_matrix`.
    pub port_list: Vec<u32>,
}

/// Result of a driven frequency sweep.
pub struct DrivenResult {
    /// S-parameter results for each frequency point.
    pub freq_results: Vec<FreqResult>,
    /// Real part of the nodal potential at the frequency of peak |S11| response.
    pub peak_phi: Vec<f64>,
    /// Frequency [Hz] at which `peak_phi` was recorded.
    pub peak_freq_hz: f64,
    /// Far-field pattern at peak frequency (empty if FarField not configured).
    pub far_field_pattern: Vec<far_field::FarFieldPoint>,
    /// Vector Fitting circuit model (Some if CircuitSynthesis = true and fit succeeded).
    pub circuit_model: Option<vf::VfModel>,
}

impl DrivenResult {
    /// Generate circuit synthesis artifacts as strings.
    /// Returns `(touchstone_s1p, circuit_model_csv, spice_netlist)`.
    /// All three are `None` if `circuit_model` is `None`.
    pub fn circuit_artifacts(&self) -> (Option<String>, Option<String>, Option<String>) {
        match &self.circuit_model {
            None => (None, None, None),
            Some(m) => {
                let freqs: Vec<f64> = self.freq_results.iter().map(|r| r.freq_hz).collect();
                let s11: Vec<num_complex::Complex64> = self.freq_results.iter()
                    .map(|r| num_complex::Complex64::new(r.s11_re, r.s11_im))
                    .collect();
                (
                    Some(vf::write_touchstone_s1p(&freqs, &s11, 50.0)),
                    Some(vf::write_circuit_model_csv(m)),
                    Some(vf::write_spice_netlist(m, 50.0)),
                )
            }
        }
    }
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
    log::info!("\n=== Driven (frequency-domain) solver ===\n");

    let drv_cfg = config.solver.driven.as_ref().ok_or_else(|| {
        RemError::Config("Driven problem requires a [Solver.Driven] section".into())
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
    log::info!("  Frequency range    = {:.3e} to {:.3e} Hz", drv_cfg.min_freq, drv_cfg.max_freq);
    log::info!("  Frequency step     = {:.3e} Hz", drv_cfg.freq_step);
    log::info!("  Save step          = {}", drv_cfg.save_step);
    log::info!("");

    let domain_map = DomainMap::from_config(config)?;

    // AMR pre-refinement: refine mesh at center frequency before the full sweep
    let amr_cfg = &config.model.refinement;
    let max_amr_iter = if amr_cfg.max_iter > 0 { amr_cfg.max_iter } else { 0 };
    let amr_theta    = if amr_cfg.tol > 0.0 { amr_cfg.tol } else { 0.5 };

    let refined_mesh: RemMesh;
    let work_mesh: &RemMesh = if max_amr_iter > 0 {
        log::info!("Adaptive mesh refinement (AMR) — pre-sweep:");
        log::info!("  Max iterations = {}", max_amr_iter);
        log::info!("  Dörfler marking = {:.1}%", amr_theta * 100.0);
        log::info!("");

        let f_center = (drv_cfg.min_freq + drv_cfg.max_freq) * 0.5;
        let k_wave = 2.0 * std::f64::consts::PI * f_center / C0;
        let k2 = k_wave * k_wave;
        let (excited_port, _port_kind) = find_excited_port(mesh);

        let eps_fn = |tag: u32| domain_map.get(tag).epsilon_abs();
        let mut cur_mesh = mesh.clone();
        for amr_iter in 1..=max_amr_iter {
            let k_mat = if domain_map.any_anisotropic() {
                let tensor_fn = |tag: u32| domain_map.get(tag).epsilon_tensor;
                assemble_stiffness_aniso(&cur_mesh, tensor_fn)?.to_csr()
            } else {
                assemble_stiffness(&cur_mesh, eps_fn)?.to_csr()
            };
            let m_mat = assemble_mass(&cur_mesh, eps_fn)?.to_csr();
            // Use real system for AMR estimator (estimator only needs shape, not phase)
            let a_mat = shifted_matrix(&k_mat, &m_mat, k2, cur_mesh.n_nodes());
            let mut a_bc = a_mat;
            let dofs = collect_dirichlet_dofs(&cur_mesh, excited_port, 1.0);
            let mut rhs = vec![0.0f64; cur_mesh.n_nodes()];
            apply_dirichlet(&mut a_bc, &mut rhs, &dofs);

            // Use real solve for AMR only (estimator doesn't need complex accuracy)
            use rem_core::solve_spd;
            let lin = &config.solver.linear;
            let result = solve_spd(&a_bc, &rhs, lin.tol, lin.max_iter, comm);
            let phi = result.solution;

            let eta = amr::zz_estimator(&cur_mesh, &phi);
            let total_err: f64 = eta.iter().map(|&e| e * e).sum::<f64>().sqrt();
            log::info!("  Iteration {}: {} nodes, error = {:.3e}", amr_iter, cur_mesh.n_nodes(), total_err);

            let marked = amr::dorfler_mark(&eta, amr_theta);
            if marked.is_empty() {
                log::info!("  → Converged: no elements marked for refinement");
                break;
            }
            let (fine_mesh, _) = amr::refine_marked(&cur_mesh, &marked);
            cur_mesh = fine_mesh;
        }
        log::info!("  → Complete: {} nodes\n", cur_mesh.n_nodes());
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
    // Adaptive: if adaptive_tol > 0, densify near rapid S11 changes
    let adaptive_tol = drv_cfg.adaptive_tol;
    const ADAPTIVE_MAX_PASSES: usize = 3;
    const ADAPTIVE_MIN_STEP_RATIO: f64 = 0.1; // don't add points closer than 10% of original step

    // Assemble K and M once (frequency-independent, real part)
    let eps_fn = |tag: u32| domain_map.get(tag).epsilon_abs();
    let k_mat = if domain_map.any_anisotropic() {
        log::info!("Anisotropic material(s) detected — using tensor stiffness assembly.");
        let tensor_fn = |tag: u32| domain_map.get(tag).epsilon_tensor;
        assemble_stiffness_aniso(mesh, tensor_fn)?.to_csr()
    } else {
        assemble_stiffness(mesh, eps_fn)?.to_csr()
    };
    let m_mat = assemble_mass(mesh, eps_fn)?.to_csr();

    // Flag: does any material have Drude-Lorentz poles (frequency-dependent ε)?
    let any_freq_dep = domain_map.any_frequency_dependent();
    if any_freq_dep {
        log::info!("Drude-Lorentz material(s) detected: stiffness matrix will be rebuilt at each frequency step");
    }

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

    let n = mesh.n_nodes();
    let k_dense = csr_to_complex_dense(&k_mat, n);
    let m_dense = csr_to_complex_dense(&m_mat, n);

    // Collect all ports for multi-port S-matrix
    let all_ports = collect_all_ports(mesh);
    let n_ports = all_ports.len();
    // Pre-compute wave-port modes for all wave ports
    let mut wave_modes: HashMap<u32, PortMode> = HashMap::new();
    let mut wave_port_support_regions: Vec<PortSupportRegionSummary> = Vec::new();
    for (pidx, kind) in &all_ports {
        if let PortKind::Wave(idx) = kind {
            if let Some(summary) = collect_port_support_region(mesh, *idx) {
                wave_port_support_regions.push(summary);
            }
            if let Some(m) = compute_wave_port_mode(mesh, *idx) {
                log::info!("WavePort {idx}: k_c={:.4e} rad/m", m.kc);
                wave_modes.insert(*pidx, m);
            } else {
                log::warn!("WavePort {idx}: modal solve failed; falling back to TEM");
            }
        }
    }

    // Backward-compat single-port variables
    let (excited_port, port_kind) = if n_ports > 0 {
        let (idx, kind) = &all_ports[0];
        (Some(*idx), kind.clone())
    } else {
        find_excited_port(mesh)
    };
    let wave_port_mode: Option<PortMode> = match &port_kind {
        PortKind::Wave(idx) => wave_modes.get(idx).cloned(),
        _ => None,
    };

    if n_ports > 1 {
        log::info!("Multi-port S-matrix: {} ports — will run {} excitations per frequency", n_ports, n_ports);
    }

    // ── ROM basis construction (single-port only) ────────────────────────────
    // When rom_order > 0 and we have a single port, pre-compute full solutions at
    // `rom_order` expansion frequencies, build an orthonormal basis, and use the
    // reduced system for all non-expansion frequency points.
    let rom_order = drv_cfg.rom_order;
    let use_rom = rom_order >= 2 && n_ports <= 1 && !any_freq_dep;
    let rom_basis: Option<rom::RomBasis> = if use_rom {
        let f_min = freqs.first().copied().unwrap_or(0.0);
        let f_max = freqs.last().copied().unwrap_or(f_min);
        let exp_freqs = rom::choose_expansion_freqs(f_min, f_max, rom_order);
        log::info!("ROM: building basis from {} full solves at expansion frequencies", rom_order);

        let dofs_snap = collect_dirichlet_dofs(mesh, excited_port, 1.0);
        let mut snapshots: Vec<Vec<Complex64>> = Vec::with_capacity(rom_order);

        for &f_exp in &exp_freqs {
            let omega_e = 2.0 * std::f64::consts::PI * f_exp;
            let k_e = omega_e / C0;
            let k2_e = k_e * k_e;
            let mut a_e = k_dense.clone();
            for i in 0..n {
                for j in 0..n {
                    a_e[(i, j)] -= Complex64::new(k2_e, 0.0) * m_dense[(i, j)];
                }
            }
            if let (Some(kl), Some(ml), Some(kc), Some(mc)) =
                (&k_loss_dense, &m_loss_dense, &k_cond_dense, &m_cond_dense)
            {
                for i in 0..n {
                    for j in 0..n {
                        let tan = Complex64::new(0.0, 1.0)
                            * (kl[(i,j)] - Complex64::new(k2_e, 0.0) * ml[(i,j)]);
                        let cond = Complex64::new(0.0, -1.0 / omega_e)
                            * (kc[(i,j)] - Complex64::new(k2_e, 0.0) * mc[(i,j)]);
                        a_e[(i, j)] += tan + cond;
                    }
                }
            }
            let mut rhs_e = vec![Complex64::ZERO; n];
            apply_dirichlet_complex(&mut a_e, &mut rhs_e, &dofs_snap);
            let use_pcg_snap = std::env::var("REM_USE_PCG").is_ok();
            match if use_pcg_snap {
                solve_complex_helmholtz_adaptive(&a_e, &rhs_e, lin.tol, lin.max_iter, true)
            } else {
                gmres_complex(&a_e, &rhs_e, lin.tol, lin.max_iter)
            } {
                Ok(phi_c) => snapshots.push(phi_c),
                Err(e) => {
                    log::warn!("ROM: expansion solve at f={f_exp:.3e} Hz failed ({e}); disabling ROM");
                    return Err(e);
                }
            }
        }

        let basis = rom::RomBasis::from_snapshots(snapshots, 1e-12);
        log::info!("ROM basis: {} vectors (r={}) from {} snapshots", basis.r(), basis.r(), rom_order);
        Some(basis)
    } else {
        if rom_order > 0 {
            if n_ports > 1 {
                log::warn!("ROM: disabled for multi-port problems (not yet supported)");
            } else if any_freq_dep {
                log::warn!("ROM: disabled when Drude-Lorentz materials are present");
            } else {
                log::warn!("ROM: rom_order must be >= 2; disabling");
            }
        }
        None
    };

    // Set of expansion frequencies for quick lookup when ROM is active
    let rom_expansion_set: std::collections::HashSet<u64> = if let Some(_) = &rom_basis {
        let f_min = freqs.first().copied().unwrap_or(0.0);
        let f_max = freqs.last().copied().unwrap_or(f_min);
        rom::choose_expansion_freqs(f_min, f_max, rom_order)
            .iter()
            .map(|&f| f.to_bits())
            .collect()
    } else {
        std::collections::HashSet::new()
    };
    let rom_full_solves = if use_rom { rom_order } else { 0 };
    let rom_fast_solves = freqs.len().saturating_sub(rom_full_solves);
    if use_rom {
        log::info!("ROM: {rom_full_solves} full solves + {rom_fast_solves} reduced solves ({} total)", freqs.len());
    }

    for (step, &freq) in freqs.iter().enumerate() {
        let omega = 2.0 * std::f64::consts::PI * freq;
        let k_wave = omega / C0;
        let k2 = k_wave * k_wave;

        // Build base matrix A(ω) = K − k²M + losses (frequency-dependent, but port-independent)
        let a_base = {
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
                        let tan_contrib = Complex64::new(0.0, 1.0)
                            * (kl[(i,j)] - Complex64::new(k2, 0.0) * ml[(i,j)]);
                        let cond_contrib = Complex64::new(0.0, -1.0 / omega)
                            * (kc[(i,j)] - Complex64::new(k2, 0.0) * mc[(i,j)]);
                        a[(i, j)] += tan_contrib + cond_contrib;
                    }
                }
            }
            // Drude-Lorentz correction: Δε(ω) = ε_complex(ω) − ε_static
            // Assemble two extra stiffness matrices (re and im part of Δε) and add to A.
            if any_freq_dep {
                use rem_core::constants::EPS0 as EPS0_CONST;
                let dl_re_fn = |tag: u32| {
                    let m = domain_map.get(tag);
                    if m.has_drude_lorentz() {
                        let (re, _) = m.epsilon_complex(freq);
                        (re - m.permittivity) * EPS0_CONST
                    } else {
                        0.0
                    }
                };
                let dl_im_fn = |tag: u32| {
                    let m = domain_map.get(tag);
                    if m.has_drude_lorentz() {
                        let (_, im) = m.epsilon_complex(freq);
                        im * EPS0_CONST // already includes static loss tangent; subtract it to avoid double-counting
                        // NOTE: the static loss tangent is handled by k_loss_dense, so we subtract it here
                        - (-(m.permittivity * m.loss_tangent
                            + if omega > 0.0 { m.conductivity / (omega * EPS0_CONST) } else { 0.0 })) * EPS0_CONST
                    } else {
                        0.0
                    }
                };
                let k_dl_re = assemble_stiffness(mesh, dl_re_fn)?.to_csr();
                let k_dl_im = assemble_stiffness(mesh, dl_im_fn)?.to_csr();
                let k_dl_re_d = csr_to_complex_dense(&k_dl_re, n);
                let k_dl_im_d = csr_to_complex_dense(&k_dl_im, n);
                for i in 0..n {
                    for j in 0..n {
                        // Real Δε correction to stiffness
                        a[(i, j)] += k_dl_re_d[(i, j)];
                        // Imaginary Δε correction: jΔε_im · K_basis
                        a[(i, j)] += Complex64::new(0.0, 1.0) * k_dl_im_d[(i, j)];
                    }
                }
            }
            a
        };
        // Resistive thin-sheet surface contribution: A[i,j] += (jω/Rs) ∫_Γ φ_i φ_j dS
        let mut a_base = a_base;
        apply_resistive_sheet(&mut a_base, mesh, omega);
        // Surface impedance BC: A[i,j] += Ys(ω) ∫_Γ φ_i φ_j dS, Ys = 1/(Rs + jωLs + 1/(jωCs))
        apply_surface_impedance(&mut a_base, mesh, omega);
        // Silver-Müller ABC: A[i,j] += jk ∫_Γ φ_i φ_j dS  (1st-order absorbing BC)
        apply_absorbing_bc(&mut a_base, mesh, k_wave);

        // ── Multi-port S-matrix path ──────────────────────────────────────────
        let (s11, s_matrix, phi_re) = if n_ports > 1 {
            let omega = 2.0 * std::f64::consts::PI * freq;
            let z0_vec: Vec<Complex64> = all_ports.iter().map(|(pidx, kind)| {
                match kind {
                    PortKind::Wave(idx) => {
                        let z = wave_modes.get(idx).map(|m| m.impedance(freq)).unwrap_or(50.0);
                        Complex64::new(if z.is_finite() { z } else { 50.0 }, 0.0)
                    }
                    _ => lumped_port_impedance(mesh, Some(*pidx), omega),
                }
            }).collect();

            // N solves: excite each port in turn, short-circuit all others
            let mut z_cols: Vec<Vec<Complex64>> = Vec::with_capacity(n_ports);
            let mut first_phi_re = Vec::new();
            for (j, (exc_idx, exc_kind)) in all_ports.iter().enumerate() {
                let vols = solve_one_excitation(
                    &a_base, mesh, freq,
                    *exc_idx, exc_kind, &all_ports,
                    &wave_modes, config,
                )?;
                // Normalize: V_i when port j excited with V=1, others shorted
                // Z_ij = V_i (since I_j ≈ 1 for unit voltage, short-circuit other ports)
                z_cols.push(vols);
                if j == 0 {
                    // Extract phi_re for the first excitation (for peak field tracking)
                    let dofs_first = collect_dirichlet_dofs(mesh, Some(*exc_idx), 1.0);
                    let mut a_tmp = a_base.clone();
                    let mut rhs_tmp = vec![Complex64::ZERO; n];
                    apply_dirichlet_complex(&mut a_tmp, &mut rhs_tmp, &dofs_first);
                    let use_pcg_f = std::env::var("REM_USE_PCG").is_ok();
                    let phi_c_first = if use_pcg_f {
                        solve_complex_helmholtz_adaptive(&a_tmp, &rhs_tmp, lin.tol, lin.max_iter, true)?
                    } else {
                        gmres_complex(&a_tmp, &rhs_tmp, lin.tol, lin.max_iter)?
                    };
                    first_phi_re = phi_c_first.iter().map(|x| x.re).collect();
                }
            }

            // Debug: print Z-matrix before conversion
            eprintln!("\n[Driven] Z-matrix at f={:.3e} Hz (n_ports={}):", freq, n_ports);
            for i in 0..n_ports {
                let row_str = (0..n_ports).map(|j| format!("{:.4e}", z_cols[j][i]))
                    .collect::<Vec<_>>().join(", ");
                eprintln!("  Z[{}] = [{}]", i, row_str);
                println!("[Z_MAT] row={} {}", i, row_str);
            }

            let s_mat = z_to_s_matrix(&z_cols, &z0_vec);
            let s11_mp = s_mat[0][0];
            
            // Debug: print S-matrix after conversion
            eprintln!("[Driven] S-matrix at f={:.3e} Hz:", freq);
            for i in 0..n_ports {
                let row_str = (0..n_ports).map(|j| format!("{}", s_mat[i][j]))
                    .collect::<Vec<_>>().join(", ");
                eprintln!("  S[{}] = [{}]", i, row_str);
                println!("[S_MAT] row={} {}", i, row_str);
            }

            log::info!(
                "f={:.3e} Hz  |S11|={:.4}  (N={} ports)",
                freq, s11_mp.norm(), n_ports
            );
            (s11_mp, s_mat, first_phi_re)
        } else {
            // ── Single-port path (backward compat) ───────────────────────────
            let dofs: HashMap<usize, f64> = if let Some(mode) = &wave_port_mode {
                if mode.is_propagating(freq) {
                    collect_dirichlet_dofs_modal(mesh, excited_port, mode)
                } else {
                    log::warn!("f={freq:.3e} Hz is below WavePort cutoff (evanescent); using φ=0");
                    collect_dirichlet_dofs(mesh, excited_port, 0.0)
                }
            } else if let Some(nf_path) = &drv_cfg.near_field_source {
                // Near-field linked source: interpolate E from CSV onto port nodes
                let nf_path = Path::new(nf_path);
                let port_tags = {
                    let mut s = std::collections::HashSet::new();
                    for belem in &mesh.boundary_elements {
                        match mesh.boundary_tags.get(&belem.tag) {
                            Some(BoundaryTag::LumpedPort { index, .. }) => {
                                if Some(*index) == excited_port { s.insert(belem.tag); }
                            }
                            Some(BoundaryTag::WavePort { index }) => {
                                if Some(*index) == excited_port { s.insert(belem.tag); }
                            }
                            _ => {}
                        }
                    }
                    s
                };
                near_field::build_near_field_dirichlet(mesh, nf_path, &port_tags)?
            } else {
                collect_dirichlet_dofs(mesh, excited_port, 1.0)
            };
            let mut a = a_base.clone();
            let mut rhs_c = vec![Complex64::ZERO; n];
            apply_dirichlet_complex(&mut a, &mut rhs_c, &dofs);
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

            // ── ROM fast path ─────────────────────────────────────────────────
            // If ROM basis is available and this is not an expansion frequency,
            // solve the cheap reduced system instead of the full GMRES.
            let phi_c = if let Some(basis) = &rom_basis {
                let is_expansion = rom_expansion_set.contains(&freq.to_bits());
                if is_expansion {
                    // Full solve — result is already in the snapshots used for basis
                    // construction, but we re-solve here for correct a_base(ω) with BCs.
                    let use_pcg_r = std::env::var("REM_USE_PCG").is_ok();
                    if use_pcg_r {
                        solve_complex_helmholtz_adaptive(&a, &rhs_c, lin.tol, lin.max_iter, true)?
                    } else {
                        gmres_complex(&a, &rhs_c, lin.tol, lin.max_iter)?
                    }
                } else {
                    // ROM solve: project A(ω) and b down to r×r, solve, expand back.
                    let b_r = basis.project_rhs(&rhs_c);
                    let a_r = basis.project_matrix_mv(|v| rom::dense_matvec(&a, v));
                    match rom::solve_reduced(a_r, b_r) {
                        Some(x_r) => basis.expand(&x_r),
                        None => {
                            log::warn!("ROM: reduced system singular at f={freq:.3e} Hz; falling back to full solve");
                            let use_pcg_rs = std::env::var("REM_USE_PCG").is_ok();
                            if use_pcg_rs {
                                solve_complex_helmholtz_adaptive(&a, &rhs_c, lin.tol, lin.max_iter, true)?
                            } else {
                                gmres_complex(&a, &rhs_c, lin.tol, lin.max_iter)?
                            }
                        }
                    }
                }
            } else {
                let use_pcg = std::env::var("REM_USE_PCG").is_ok();
                if use_pcg {
                    log::info!("Phase2: single-port main-loop PCG solve at f={freq:.3e} Hz");
                    solve_complex_helmholtz_adaptive(&a, &rhs_c, lin.tol, lin.max_iter, true)?
                } else {
                    gmres_complex(&a, &rhs_c, lin.tol, lin.max_iter)?
                }
            };

            let phi_re: Vec<f64> = phi_c.iter().map(|x| x.re).collect();
            let (v_port, i_kphi) = compute_port_vi_complex(mesh, &phi_c, &k_dense, excited_port);
            let omega = 2.0 * std::f64::consts::PI * freq;
            let (z0, i_port) = if let Some(mode) = &wave_port_mode {
                let z = mode.impedance(freq);
                let z_use = if z.is_finite() { z } else { 50.0 };
                let z0c = Complex64::new(z_use, 0.0);
                let i = if z_use.abs() > 1e-300 { v_port / z0c } else { Complex64::ZERO };
                (z0c, i)
            } else {
                let z0c = lumped_port_impedance(mesh, excited_port, omega);
                (z0c, i_kphi)
            };
            let s11 = if i_port.norm() > 1e-300 {
                let z = v_port / i_port;
                let z0c = z0;
                (z - z0c) / (z + z0c)
            } else {
                Complex64::ZERO
            };
            log::info!(
                "f={:.3e} Hz  |S11|={:.4}  ∠S11={:.2}°  Z0={:.1}Ω",
                freq, s11.norm(), s11.arg().to_degrees(), z0
            );
            (s11, vec![], phi_re)
        };

        let port_list: Vec<u32> = all_ports.iter().map(|(i, _)| *i).collect();
        freq_results.push(FreqResult {
            freq_hz: freq,
            s11_re: s11.re,
            s11_im: s11.im,
            s_matrix,
            port_list,
        });

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

    // ── Adaptive frequency densification ──────────────────────────────────────
    // If adaptive_tol > 0: find intervals where |ΔS11| > tol * max|ΔS11| and
    // bisect them (up to ADAPTIVE_MAX_PASSES times, min spacing enforced).
    if adaptive_tol > 0.0 && freq_results.len() >= 2 {
        let orig_step = (freq_results.last().unwrap().freq_hz
                       - freq_results.first().unwrap().freq_hz)
                       / (freq_results.len() as f64 - 1.0).max(1.0);
        let min_spacing = orig_step * ADAPTIVE_MIN_STEP_RATIO;

        for pass in 0..ADAPTIVE_MAX_PASSES {
            // Compute |S11| array
            let mags: Vec<f64> = freq_results.iter()
                .map(|r| (r.s11_re*r.s11_re + r.s11_im*r.s11_im).sqrt())
                .collect();
            let max_delta = mags.windows(2)
                .map(|w| (w[1] - w[0]).abs())
                .fold(0.0f64, f64::max);
            if max_delta < 1e-15 { break; }
            let threshold = adaptive_tol * max_delta;

            let mut extra_freqs: Vec<f64> = Vec::new();
            for i in 0..mags.len().saturating_sub(1) {
                let delta = (mags[i+1] - mags[i]).abs();
                if delta > threshold {
                    let f_mid = 0.5 * (freq_results[i].freq_hz + freq_results[i+1].freq_hz);
                    let gap = (freq_results[i+1].freq_hz - freq_results[i].freq_hz).abs();
                    if gap > min_spacing * 2.0 { extra_freqs.push(f_mid); }
                }
            }
            if extra_freqs.is_empty() { break; }
            log::info!("Adaptive pass {}: inserting {} extra frequency points", pass + 1, extra_freqs.len());

            for &freq in &extra_freqs {
                let omega_a = 2.0 * std::f64::consts::PI * freq;
                let k_wave_a = omega_a / C0;
                let k2_a = k_wave_a * k_wave_a;
                // Build A_base for this adaptive frequency
                let a_base_a = {
                    let mut a = k_dense.clone();
                    for i in 0..n { for j in 0..n {
                        a[(i,j)] -= Complex64::new(k2_a, 0.0) * m_dense[(i,j)];
                    }}
                    if let (Some(kl), Some(ml), Some(kc), Some(mc)) =
                        (&k_loss_dense, &m_loss_dense, &k_cond_dense, &m_cond_dense)
                    {
                        for i in 0..n { for j in 0..n {
                            let tc = Complex64::new(0.0, 1.0) * (kl[(i,j)] - Complex64::new(k2_a, 0.0) * ml[(i,j)]);
                            let cc = Complex64::new(0.0, -1.0 / omega_a) * (kc[(i,j)] - Complex64::new(k2_a, 0.0) * mc[(i,j)]);
                            a[(i,j)] += tc + cc;
                        }}
                    }
                    a
                };
                let (s11, s_mat, phi_re) = if n_ports > 1 {
                    let omega = 2.0 * std::f64::consts::PI * freq;
                    let z0_vec: Vec<Complex64> = all_ports.iter().map(|(pidx, kind)| match kind {
                        PortKind::Wave(idx) => { let z = wave_modes.get(idx).map(|m| m.impedance(freq)).unwrap_or(50.0); Complex64::new(if z.is_finite() { z } else { 50.0 }, 0.0) }
                        _ => lumped_port_impedance(mesh, Some(*pidx), omega),
                    }).collect();
                    let mut z_cols = Vec::with_capacity(n_ports);
                    let mut first_phi_re = Vec::new();
                    for (j, (exc_idx, exc_kind)) in all_ports.iter().enumerate() {
                        let vols = solve_one_excitation(&a_base_a, mesh, freq, *exc_idx, exc_kind, &all_ports, &wave_modes, config)?;
                        z_cols.push(vols);
                        if j == 0 {
                            let dofs_f = collect_dirichlet_dofs(mesh, Some(*exc_idx), 1.0);
                            let mut at = a_base_a.clone();
                            let mut rt = vec![Complex64::ZERO; n];
                            apply_dirichlet_complex(&mut at, &mut rt, &dofs_f);
                            let use_pcg_af = std::env::var("REM_USE_PCG").is_ok();
                            let pc = if use_pcg_af {
                                solve_complex_helmholtz_adaptive(&at, &rt, lin.tol, lin.max_iter, true)?
                            } else {
                                gmres_complex(&at, &rt, lin.tol, lin.max_iter)?
                            };
                            first_phi_re = pc.iter().map(|x| x.re).collect();
                        }
                    }
                    let s_mat = z_to_s_matrix(&z_cols, &z0_vec);
                    let s11 = s_mat[0][0];
                    (s11, s_mat, first_phi_re)
                } else {
                    let dofs = if let Some(mode) = &wave_port_mode {
                        if mode.is_propagating(freq) { collect_dirichlet_dofs_modal(mesh, excited_port, mode) }
                        else { collect_dirichlet_dofs(mesh, excited_port, 0.0) }
                    } else { collect_dirichlet_dofs(mesh, excited_port, 1.0) };
                    let mut a = a_base_a;
                    let mut rhs_c = vec![Complex64::ZERO; n];
                    apply_dirichlet_complex(&mut a, &mut rhs_c, &dofs);
                    if !config.domains.current_dipole.is_empty() {
                        let jw_mu0 = Complex64::new(0.0, omega_a) * Complex64::new(4.0e-7 * std::f64::consts::PI, 0.0);
                        for dipole in &config.domains.current_dipole {
                            let node = nearest_node(mesh, dipole);
                            if !dofs.contains_key(&node) {
                                let dir_mag = (dipole.direction.iter().map(|x| x * x).sum::<f64>()).sqrt();
                                rhs_c[node] += jw_mu0 * Complex64::new(dipole.moment * dir_mag.max(1.0), 0.0);
                            }
                        }
                    }
                    let use_pcg = std::env::var("REM_USE_PCG").is_ok();
                    let phi_c = if use_pcg {
                        log::info!("Phase2: single-port Helmholtz attempting PCG solve");
                        solve_complex_helmholtz_adaptive(&a, &rhs_c, lin.tol, lin.max_iter, true)?
                    } else {
                        gmres_complex(&a, &rhs_c, lin.tol, lin.max_iter)?
                    };
                    let phi_re: Vec<f64> = phi_c.iter().map(|x| x.re).collect();
                    let (v_port, i_kphi) = compute_port_vi_complex(mesh, &phi_c, &k_dense, excited_port);
                    let omega = 2.0 * std::f64::consts::PI * freq;
                    let (z0, i_port) = if let Some(mode) = &wave_port_mode {
                        let z = mode.impedance(freq);
                        let z_use = if z.is_finite() { z } else { 50.0 };
                        let z0c = Complex64::new(z_use, 0.0);
                        (z0c, if z_use.abs() > 1e-300 { v_port / z0c } else { Complex64::ZERO })
                    } else {
                        let z0c = lumped_port_impedance(mesh, excited_port, omega);
                        (z0c, i_kphi)
                    };
                    let s11 = if i_port.norm() > 1e-300 {
                        let z = v_port / i_port; let z0c = z0;
                        (z - z0c) / (z + z0c)
                    } else { Complex64::ZERO };
                    (s11, vec![], phi_re)
                };
                log::info!("Adaptive f={:.3e} Hz  |S11|={:.4}", freq, s11.norm());
                let port_list: Vec<u32> = all_ports.iter().map(|(i, _)| *i).collect();
                freq_results.push(FreqResult { freq_hz: freq, s11_re: s11.re, s11_im: s11.im, s_matrix: s_mat, port_list });
                let s11_mag = s11.norm();
                if s11_mag > peak_s11_mag {
                    peak_s11_mag = s11_mag;
                    peak_freq_hz = freq;
                    peak_phi = phi_re;
                }
            }
            // Re-sort by frequency after inserting adaptive points
            freq_results.sort_by(|a, b| a.freq_hz.partial_cmp(&b.freq_hz).unwrap());
        }
    }

    #[cfg(not(target_arch = "wasm32"))]
    output::write_s_params(out_dir, &freq_results)?;
    #[cfg(not(target_arch = "wasm32"))]
    if !wave_port_support_regions.is_empty() {
        output::write_wave_port_support_regions(out_dir, &wave_port_support_regions)?;
    }
    #[cfg(not(target_arch = "wasm32"))]
    if !peak_phi.is_empty() {
        let peak_energy = postprocess::electrostatic_energy(&peak_phi, mesh, |tag| domain_map.get(tag).epsilon_abs());
        let domain_energies = peak_domain_energy_records(mesh, domain_map, &peak_phi, peak_energy);
        if !domain_energies.is_empty() {
            output::write_peak_domain_energy(out_dir, peak_freq_hz, &domain_energies)?;
        }
    }

    // Near-to-far-field transform (if configured)
    let mut far_field_pattern: Vec<far_field::FarFieldPoint> = Vec::new();
    #[cfg(not(target_arch = "wasm32"))]
    if let Some(ff_cfg) = &config.solver.far_field {
        if !peak_phi.is_empty() {
            log::info!("[REM] Computing far-field pattern at peak frequency {:.3e} Hz", peak_freq_hz);
            far_field_pattern = far_field::compute_far_field(mesh, &peak_phi, peak_freq_hz, ff_cfg);
            if !far_field_pattern.is_empty() {
                far_field::write_far_field_csv(out_dir, &far_field_pattern, peak_freq_hz)?;
            }
        }
    }
    #[cfg(target_arch = "wasm32")]
    if let Some(ff_cfg) = &config.solver.far_field {
        if !peak_phi.is_empty() {
            far_field_pattern = far_field::compute_far_field(mesh, &peak_phi, peak_freq_hz, ff_cfg);
        }
    }

    // Near-field export (if configured)
    #[cfg(not(target_arch = "wasm32"))]
    if let Some(nf_cfg) = &config.postprocessing.near_field {
        if !peak_phi.is_empty() {
            log::info!("[REM] Exporting near-field data on {} boundary attributes", nf_cfg.attributes.len());
            let nf_points = near_field::export_near_field(mesh, &peak_phi, nf_cfg)?;
            if !nf_points.is_empty() {
                near_field::write_near_field(Path::new(out_dir), &nf_points, nf_cfg)?;
            }
        }
    }

    log::info!(
        "Driven solve complete: {} frequency points, peak |S11|={:.4} at {:.3e} Hz",
        freq_results.len(), peak_s11_mag, peak_freq_hz
    );
    report_peak_memory("Driven solver");

    // ── Vector Fitting circuit synthesis ────────────────────────────────────
    let circuit_model: Option<vf::VfModel> = if drv_cfg.circuit_synthesis {
        if freq_results.len() < 4 {
            log::warn!("CircuitSynthesis: need ≥ 4 frequency points ({} available); skipping",
                freq_results.len());
            None
        } else {
            let freqs_hz: Vec<f64> = freq_results.iter().map(|r| r.freq_hz).collect();
            let s11_data: Vec<Complex64> = freq_results.iter()
                .map(|r| Complex64::new(r.s11_re, r.s11_im))
                .collect();
            let n_poles = if drv_cfg.rom_order >= 2 {
                drv_cfg.rom_order.min(32)
            } else {
                (freq_results.len() / 4).clamp(4, 16)
            };
            log::info!("CircuitSynthesis: VF with {} poles over {} points", n_poles, freqs_hz.len());
            match vf::vector_fit(&freqs_hz, &s11_data, n_poles, 10, 1e-6) {
                Some(m) => {
                    log::info!("CircuitSynthesis: VF converged, RMS = {:.4e}", m.rms_error);

                    // Write files on native (non-WASM) path
                    #[cfg(not(target_arch = "wasm32"))]
                    {
                        use std::io::Write as _;
                        let ts = vf::write_touchstone_s1p(&freqs_hz, &s11_data, 50.0);
                        let ts_path = Path::new(out_dir).join("s_params.s1p");
                        if let Ok(mut f) = std::fs::File::create(&ts_path) {
                            let _ = f.write_all(ts.as_bytes());
                            log::info!("Wrote {}", ts_path.display());
                        }
                        let csv = vf::write_circuit_model_csv(&m);
                        let csv_path = Path::new(out_dir).join("circuit_model.csv");
                        if let Ok(mut f) = std::fs::File::create(&csv_path) {
                            let _ = f.write_all(csv.as_bytes());
                            log::info!("Wrote {}", csv_path.display());
                        }
                        let spice = vf::write_spice_netlist(&m, 50.0);
                        let spice_path = Path::new(out_dir).join("equivalent_circuit.cir");
                        if let Ok(mut f) = std::fs::File::create(&spice_path) {
                            let _ = f.write_all(spice.as_bytes());
                            log::info!("Wrote {}", spice_path.display());
                        }
                    }

                    Some(m)
                }
                None => {
                    log::warn!("CircuitSynthesis: VF solve failed; no circuit model produced");
                    None
                }
            }
        }
    } else {
        None
    };

    Ok(DrivenResult { freq_results, peak_phi, peak_freq_hz, far_field_pattern, circuit_model })
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

/// Try to solve with PCG, fallback to GMRES if conversion or convergence fails.
/// Enable PCG attempt via use_pcg=true parameter.
fn solve_complex_helmholtz_adaptive(
    a: &DMatrix<Complex64>,
    rhs: &[Complex64],
    tol: f64,
    max_iter: usize,
    use_pcg: bool,
) -> RemResult<Vec<Complex64>> {
    if !use_pcg {
        return gmres_complex(a, rhs, tol, max_iter);
    }

    // Convert dense DMatrix to CSR format
    let mat_csr = CsrMatrixComplex::from_dense(a);
    
    // Attempt PCG solve
    let result = solve_pcg_complex(&mat_csr, rhs, tol, max_iter);
    
    if result.converged {
        log::info!("PCG: converged in {} iterations (residual {:.3e})", 
                    result.iterations, result.residual_norm);
        return Ok(result.solution);
    }

    // PCG diverged or hit max iterations; fallback to GMRES
    log::info!("PCG: no convergence after {} iterations (residual {:.3e}); using GMRES",
                result.iterations, result.residual_norm);
    gmres_complex(a, rhs, tol, max_iter)
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

/// Collect all distinct LumpedPort/WavePort indices in stable order.
fn collect_all_ports(mesh: &RemMesh) -> Vec<(u32, PortKind)> {
    let mut seen: Vec<(u32, PortKind)> = Vec::new();
    for bc in mesh.boundary_tags.values() {
        match bc {
            BoundaryTag::LumpedPort { index, .. } => {
                if !seen.iter().any(|(i, _)| *i == *index) {
                    seen.push((*index, PortKind::Lumped));
                }
            }
            BoundaryTag::WavePort { index } => {
                if !seen.iter().any(|(i, _)| *i == *index) {
                    seen.push((*index, PortKind::Wave(*index)));
                }
            }
            _ => {}
        }
    }
    seen.sort_by_key(|(i, _)| *i);
    seen
}

/// Add Silver-Müller first-order absorbing BC:
/// A[i,j] += jk · ∫_Γ φ_i φ_j dS
///
/// For order=2, adds a second-order correction term (not yet implemented; falls back to order=1).
fn apply_absorbing_bc(a: &mut DMatrix<Complex64>, mesh: &RemMesh, k: f64) {
    let jk = Complex64::new(0.0, k);
    for belem in &mesh.boundary_elements {
        let _order = match mesh.boundary_tags.get(&belem.tag) {
            Some(BoundaryTag::Absorbing { order }) => *order,
            _ => continue,
        };
        let nids = &belem.node_ids;
        match nids.len() {
            2 => {
                let (p0, p1) = (&mesh.nodes[nids[0]], &mesh.nodes[nids[1]]);
                let l = ((p1.x-p0.x).powi(2) + (p1.y-p0.y).powi(2) + (p1.z-p0.z).powi(2)).sqrt();
                a[(nids[0], nids[0])] += jk * l / 3.0;
                a[(nids[1], nids[1])] += jk * l / 3.0;
                a[(nids[0], nids[1])] += jk * l / 6.0;
                a[(nids[1], nids[0])] += jk * l / 6.0;
            }
            3 => {
                let (p0, p1, p2) = (&mesh.nodes[nids[0]], &mesh.nodes[nids[1]], &mesh.nodes[nids[2]]);
                let v1 = [p1.x-p0.x, p1.y-p0.y, p1.z-p0.z];
                let v2 = [p2.x-p0.x, p2.y-p0.y, p2.z-p0.z];
                let cx = v1[1]*v2[2] - v1[2]*v2[1];
                let cy = v1[2]*v2[0] - v1[0]*v2[2];
                let cz = v1[0]*v2[1] - v1[1]*v2[0];
                let area = 0.5 * (cx*cx + cy*cy + cz*cz).sqrt();
                for &ni in nids { a[(ni, ni)] += jk * area / 6.0; }
                for ii in 0..3 {
                    for jj in 0..3 {
                        if ii != jj { a[(nids[ii], nids[jj])] += jk * area / 12.0; }
                    }
                }
            }
            _ => {}
        }
    }
}

/// Add resistive thin-sheet surface contribution to system matrix:
/// A[i,j] += (jω/Rs) · ∫_Γ φ_i φ_j dS
fn apply_resistive_sheet(a: &mut DMatrix<Complex64>, mesh: &RemMesh, omega: f64) {
    let j = Complex64::new(0.0, 1.0);
    for belem in &mesh.boundary_elements {
        let rs = match mesh.boundary_tags.get(&belem.tag) {
            Some(BoundaryTag::ResistiveSheet { rs }) if *rs > 0.0 => *rs,
            _ => continue,
        };
        let scale = j * omega / rs;
        let nids = &belem.node_ids;
        match nids.len() {
            2 => {
                let (p0, p1) = (&mesh.nodes[nids[0]], &mesh.nodes[nids[1]]);
                let l = ((p1.x-p0.x).powi(2) + (p1.y-p0.y).powi(2) + (p1.z-p0.z).powi(2)).sqrt();
                a[(nids[0], nids[0])] += scale * l / 3.0;
                a[(nids[1], nids[1])] += scale * l / 3.0;
                a[(nids[0], nids[1])] += scale * l / 6.0;
                a[(nids[1], nids[0])] += scale * l / 6.0;
            }
            3 => {
                let (p0, p1, p2) = (&mesh.nodes[nids[0]], &mesh.nodes[nids[1]], &mesh.nodes[nids[2]]);
                let v1 = [p1.x-p0.x, p1.y-p0.y, p1.z-p0.z];
                let v2 = [p2.x-p0.x, p2.y-p0.y, p2.z-p0.z];
                let cx = v1[1]*v2[2] - v1[2]*v2[1];
                let cy = v1[2]*v2[0] - v1[0]*v2[2];
                let cz = v1[0]*v2[1] - v1[1]*v2[0];
                let area = 0.5 * (cx*cx + cy*cy + cz*cz).sqrt();
                for &ni in nids { a[(ni, ni)] += scale * area / 6.0; }
                for ii in 0..3 {
                    for jj in 0..3 {
                        if ii != jj { a[(nids[ii], nids[jj])] += scale * area / 12.0; }
                    }
                }
            }
            _ => {}
        }
    }
}

/// Add surface impedance BC contribution to system matrix:
/// A[i,j] += Ys(ω) · ∫_Γ φ_i φ_j dS
/// where Ys(ω) = 1 / (Rs + jωLs + 1/(jωCs))
fn apply_surface_impedance(a: &mut DMatrix<Complex64>, mesh: &RemMesh, omega: f64) {
    let j = Complex64::new(0.0, 1.0);
    for belem in &mesh.boundary_elements {
        let (rs, ls, cs) = match mesh.boundary_tags.get(&belem.tag) {
            Some(BoundaryTag::Impedance { rs, ls, cs }) => (*rs, *ls, *cs),
            _ => continue,
        };
        // Zs(ω) = Rs + jωLs + 1/(jωCs)
        let mut zs = Complex64::new(rs.max(0.0), 0.0);
        if ls > 0.0 { zs += j * omega * ls; }
        if cs > 0.0 { zs += Complex64::new(1.0, 0.0) / (j * omega * cs); }
        if zs.norm() < 1e-30 { continue; }
        let ys = Complex64::new(1.0, 0.0) / zs;
        let nids = &belem.node_ids;
        match nids.len() {
            2 => {
                let (p0, p1) = (&mesh.nodes[nids[0]], &mesh.nodes[nids[1]]);
                let l = ((p1.x-p0.x).powi(2) + (p1.y-p0.y).powi(2) + (p1.z-p0.z).powi(2)).sqrt();
                a[(nids[0], nids[0])] += ys * l / 3.0;
                a[(nids[1], nids[1])] += ys * l / 3.0;
                a[(nids[0], nids[1])] += ys * l / 6.0;
                a[(nids[1], nids[0])] += ys * l / 6.0;
            }
            3 => {
                let (p0, p1, p2) = (&mesh.nodes[nids[0]], &mesh.nodes[nids[1]], &mesh.nodes[nids[2]]);
                let v1 = [p1.x-p0.x, p1.y-p0.y, p1.z-p0.z];
                let v2 = [p2.x-p0.x, p2.y-p0.y, p2.z-p0.z];
                let cx = v1[1]*v2[2] - v1[2]*v2[1];
                let cy = v1[2]*v2[0] - v1[0]*v2[2];
                let cz = v1[0]*v2[1] - v1[1]*v2[0];
                let area = 0.5 * (cx*cx + cy*cy + cz*cz).sqrt();
                for &ni in nids { a[(ni, ni)] += ys * area / 6.0; }
                for ii in 0..3 {
                    for jj in 0..3 {
                        if ii != jj { a[(nids[ii], nids[jj])] += ys * area / 12.0; }
                    }
                }
            }
            _ => {}
        }
    }
}

/// Compute frequency-dependent lumped port impedance Z(ω) = R + jωL + 1/(jωC).
fn lumped_port_impedance(mesh: &RemMesh, port_idx: Option<u32>, omega: f64) -> Complex64 {
    let j = Complex64::new(0.0, 1.0);
    if let Some(idx) = port_idx {
        for bc in mesh.boundary_tags.values() {
            match bc {
                BoundaryTag::LumpedPort { index, r, l, c } if *index == idx => {
                    let r_val = if *r > 0.0 { *r } else { 50.0 };
                    let mut z = Complex64::new(r_val, 0.0);
                    if *l > 0.0 { z += j * omega * l; }
                    if *c > 0.0 { z += Complex64::new(1.0, 0.0) / (j * omega * c); }
                    return z;
                }
                BoundaryTag::WavePort { index } if *index == idx => {
                    return Complex64::new(50.0, 0.0);
                }
                _ => {}
            }
        }
    }
    Complex64::new(50.0, 0.0)
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

/// Compute Z-matrix row for port `exc_idx` excited at V=1, all other ports terminated
/// in Z0 (matched load — implemented via a Robin term on unexcited ports).
///
/// Returns the open-circuit port voltages Vec<Complex64> indexed by `all_ports`.
fn solve_one_excitation(
    a_base: &DMatrix<Complex64>,
    mesh: &RemMesh,
    freq: f64,
    exc_idx: u32,
    exc_kind: &PortKind,
    all_ports: &[(u32, PortKind)],
    wave_modes: &HashMap<u32, PortMode>,
    config: &PalaceConfig,
) -> RemResult<Vec<Complex64>> {
    let n = a_base.nrows();
    let lin = &config.solver.linear;
    let mut a = a_base.clone();

    // Build Dirichlet map for multi-port Z-parameter extraction:
    // Note: We're using the SHORT-CIRCUIT method (non-excited ports → 0.0).
    // This directly gives Z-parameters. An alternative would be OPEN-CIRCUIT method
    // (non-excited ports → natural BC) which gives Y-parameters requiring matrix inversion.
    // Current approach: Short-circuit is more standard for FEM solvers.
    let excited_mode = if let PortKind::Wave(idx) = exc_kind {
        wave_modes.get(idx)
    } else {
        None
    };

    let dofs: HashMap<usize, f64> = if let Some(mode) = excited_mode {
        if mode.is_propagating(freq) {
            collect_dirichlet_dofs_modal(mesh, Some(exc_idx), mode)
        } else {
            // For evanescent wave ports, also short-circuit non-excited ports
            let mut dofs = collect_dirichlet_dofs_open_circuit(mesh, Some(exc_idx), 1.0);
            // Short-circuit all other ports
            for belem in &mesh.boundary_elements {
                let bc = match mesh.boundary_tags.get(&belem.tag) {
                    Some(b) => b,
                    None => continue,
                };
                match bc {
                    BoundaryTag::LumpedPort { index, .. }
                    | BoundaryTag::WavePort { index } => {
                        if *index != exc_idx {
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
    } else {
        // Lumped ports: SHORT-CIRCUIT METHOD
        // Excited port → 1.0, all others → 0.0 (short-circuit), PEC/Ground → 0.0
        let mut dofs: HashMap<usize, f64> = HashMap::new();
        for belem in &mesh.boundary_elements {
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
                BoundaryTag::LumpedPort { index, .. }
                | BoundaryTag::Terminal { index }
                | BoundaryTag::WavePort { index } => {
                    let val = if *index == exc_idx { 1.0 } else { 0.0 };
                    for &nid in &belem.node_ids {
                        dofs.entry(nid).or_insert(val);
                    }
                }
                _ => {}
            }
        }
        dofs
    };

    let mut rhs_c = vec![Complex64::ZERO; n];
    apply_dirichlet_complex(&mut a, &mut rhs_c, &dofs);

    // CurrentDipole sources
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

    // Use PCG if enabled, fallback to GMRES on divergence
    let use_pcg = std::env::var("REM_USE_PCG").is_ok();
    let phi_c = if use_pcg {
        solve_complex_helmholtz_adaptive(&a, &rhs_c, lin.tol, lin.max_iter, true)?
    } else {
        gmres_complex(&a, &rhs_c, lin.tol, lin.max_iter)?
    };

    // Helper: collect unique node IDs for a given port index.
    let port_node_ids = |pidx: u32| -> Vec<usize> {
        mesh.boundary_elements.iter()
            .filter(|e| match mesh.boundary_tags.get(&e.tag) {
                Some(BoundaryTag::LumpedPort { index, .. }) => *index == pidx,
                Some(BoundaryTag::WavePort { index }) => *index == pidx,
                _ => false,
            })
            .flat_map(|e| e.node_ids.iter().copied())
            .collect::<std::collections::HashSet<_>>()
            .into_iter()
            .collect()
    };

    // Compute port current I_j for the excited port using the *original* (pre-BC) matrix:
    //   I_j = Σ_{i ∈ port_j} (A_base · φ_c)[i]
    // This is the short-circuit reaction current needed to form Z_ij = V_i / I_j.
    let exc_nodes = port_node_ids(exc_idx);
    let i_exc = if exc_nodes.is_empty() {
        Complex64::new(1.0, 0.0) // degenerate: no normalization
    } else {
        let phi_vec: nalgebra::DVector<Complex64> =
            nalgebra::DVector::from_iterator(n, phi_c.iter().copied());
        let a_phi = a_base * phi_vec;
        let current: Complex64 = exc_nodes.iter().map(|&row| a_phi[row]).sum();
        if current.norm() < 1e-300 { Complex64::new(1.0, 0.0) } else { current }
    };

    log::debug!("[Driven] Excited port {} with {} nodes, I_exc = {:.6e}", exc_idx, exc_nodes.len(), i_exc);

    // Extract Z-matrix column: Z_ij = V_i / I_j
    let mut voltages = Vec::with_capacity(all_ports.len());
    for (pidx, _) in all_ports {
        let port_nodes = port_node_ids(*pidx);
        let v = if port_nodes.is_empty() {
            Complex64::ZERO
        } else {
            port_nodes.iter().map(|&nd| phi_c[nd]).sum::<Complex64>()
                / Complex64::new(port_nodes.len() as f64, 0.0)
        };
        // Normalize by drive current → Z_ij [Ω]
        let z_ij = v / i_exc;
        log::debug!("  Port {}: {} nodes, V_i = {:.6e}, Z_ij = {:.6e}", pidx, port_nodes.len(), v, z_ij);
        println!("[Z] exc={} read={}: {} nodes, V={:.6e}, Z={:.6e}", exc_idx, pidx, port_nodes.len(), v, z_ij);
        voltages.push(z_ij);
    }
    Ok(voltages)
}

/// Build full S-matrix from Z-matrix using S = (Z−Z₀)(Z+Z₀)⁻¹ for the diagonal Z₀ case.
/// With Z0_vec[i] = reference impedance of port i, for N ports:
///   S_ij = 2·√(Z0_i·Z0_j) · [Y · Z0]_ij − δ_ij
/// where Y = (Z + diag(Z0))⁻¹ ... simplified for diagonal Z0:
///   S = diag(1/√Z0) · (Z − diag(Z0)) · (Z + diag(Z0))⁻¹ · diag(√Z0)
/// We use the straightforward element-wise formula for the 2-port case and full
/// matrix inversion for N > 2.
fn z_to_s_matrix(
    z_cols: &[Vec<Complex64>],   // z_cols[j] = j-th column of Z (voltage at each port when port j excited)
    z0: &[Complex64],            // reference impedance per port (complex for RLC ports)
) -> Vec<Vec<Complex64>> {
    let n = z0.len();
    // Build Z matrix
    let mut z = vec![vec![Complex64::ZERO; n]; n];
    for j in 0..n {
        for i in 0..n {
            z[i][j] = z_cols[j][i];
        }
    }
    // S = (Z − Z0)(Z + Z0)⁻¹ using nalgebra for inversion
    let z0_diag: Vec<Complex64> = z0.to_vec();
    let mut zp = DMatrix::<Complex64>::zeros(n, n);
    let mut zm = DMatrix::<Complex64>::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            zp[(i, j)] = z[i][j];
            zm[(i, j)] = z[i][j];
        }
        zp[(i, i)] += z0_diag[i];
        zm[(i, i)] -= z0_diag[i];
    }
    // Invert (Z + Z0)
    let zp_inv = match zp.try_inverse() {
        Some(inv) => inv,
        None => return vec![vec![Complex64::ZERO; n]; n],
    };
    let s_mat = zm * zp_inv;
    let mut s = vec![vec![Complex64::ZERO; n]; n];
    for i in 0..n {
        for j in 0..n {
            s[i][j] = s_mat[(i, j)];
        }
    }
    s
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
        .collect::<std::collections::HashSet<_>>()
        .into_iter()
        .collect();

    if port_nodes.is_empty() { return (Complex64::ZERO, Complex64::ZERO); }

    // Port voltage: arithmetic mean of φ over port nodes (uniform TEM or lumped-port)
    let v_port = port_nodes.iter().map(|&n| phi[n]).sum::<Complex64>()
        / Complex64::new(port_nodes.len() as f64, 0.0);

    // Port current: computed via K·φ residual over port rows.
    // For a LumpedPort this gives the net conduction current into the port.
    // For WavePort the caller overrides Z0 with Z_TE and computes I = V/Z0 instead
    // (see S11 computation in the sweep loop).
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

fn peak_domain_energy_records(
    mesh: &RemMesh,
    domain_map: &DomainMap,
    phi: &[f64],
    total_energy: f64,
) -> Vec<output::DomainEnergyRecord> {
    let mut domain_tags: Vec<u32> = mesh.domain_tags.keys().copied().collect();
    domain_tags.sort_unstable();

    if mesh.dim == 2
        && mesh.volume_elements.iter().all(|element| element.kind == ElementKind::Tri3)
        && mesh.boundary_elements.iter().all(|element| element.kind == ElementKind::Line2)
    {
        return domain_tags
            .into_iter()
            .map(|tag| {
                let (material_index, _) = domain_map.get_indexed(tag);
                let energy = extract_domain_submesh(mesh, tag)
                    .map(|submesh| {
                        let sub_phi = submesh.transfer_from_parent(phi);
                        postprocess::electrostatic_energy(
                            &sub_phi,
                            &submesh.mesh,
                            |sub_tag| domain_map.get(sub_tag).epsilon_abs(),
                        )
                    })
                    .unwrap_or(0.0);
                output::DomainEnergyRecord {
                    domain_tag: tag,
                    material_index: (material_index != usize::MAX).then_some(material_index),
                    energy,
                    fraction: if total_energy.abs() > 1e-300 { energy / total_energy } else { 0.0 },
                }
            })
            .collect();
    }

    domain_tags
        .into_iter()
        .map(|tag| {
            let (material_index, _) = domain_map.get_indexed(tag);
            let energy = postprocess::electrostatic_energy(phi, mesh, |elem_tag| {
                if elem_tag == tag {
                    domain_map.get(elem_tag).epsilon_abs()
                } else {
                    0.0
                }
            });
            output::DomainEnergyRecord {
                domain_tag: tag,
                material_index: (material_index != usize::MAX).then_some(material_index),
                energy,
                fraction: if total_energy.abs() > 1e-300 { energy / total_energy } else { 0.0 },
            }
        })
        .collect()
}

fn extract_domain_submesh(mesh: &RemMesh, domain_tag: u32) -> Option<FemSubMesh2d> {
    match extract_submesh_tri3(mesh, &[domain_tag]) {
        Ok(submesh) if !submesh.mesh.volume_elements.is_empty() => Some(submesh),
        Ok(_) => None,
        Err(err) => {
            log::warn!(
                "fem-rs Tri3 bridge submesh extraction failed for driven domain tag {} ({})",
                domain_tag,
                err
            );
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rem_materials::Material;
    use rem_mesh::{Element, Node};

    fn unit_square_mesh() -> RemMesh {
        let nodes = vec![
            Node { id: 0, x: 0.0, y: 0.0, z: 0.0 },
            Node { id: 1, x: 1.0, y: 0.0, z: 0.0 },
            Node { id: 2, x: 1.0, y: 1.0, z: 0.0 },
            Node { id: 3, x: 0.0, y: 1.0, z: 0.0 },
        ];
        let volume_elements = vec![
            Element { id: 1, kind: ElementKind::Tri3, tag: 1, node_ids: vec![0, 1, 2], rank: 0 },
            Element { id: 2, kind: ElementKind::Tri3, tag: 2, node_ids: vec![0, 2, 3], rank: 0 },
        ];
        let boundary_elements = vec![
            Element { id: 3, kind: ElementKind::Line2, tag: 10, node_ids: vec![0, 1], rank: 0 },
            Element { id: 4, kind: ElementKind::Line2, tag: 11, node_ids: vec![1, 2], rank: 0 },
            Element { id: 5, kind: ElementKind::Line2, tag: 12, node_ids: vec![2, 3], rank: 0 },
            Element { id: 6, kind: ElementKind::Line2, tag: 13, node_ids: vec![3, 0], rank: 0 },
        ];
        RemMesh {
            nodes,
            volume_elements,
            boundary_elements,
            domain_tags: [(1u32, 0usize), (2u32, 0usize)].into_iter().collect(),
            boundary_tags: [
                (10u32, BoundaryTag::Ground),
                (11u32, BoundaryTag::Ground),
                (12u32, BoundaryTag::Ground),
                (13u32, BoundaryTag::Ground),
            ].into_iter().collect(),
            dim: 2,
            rank: 0,
            size: 1,
        }
    }

    #[test]
    fn peak_domain_energy_breakdown_matches_total() {
        use rem_core::constants::EPS0;

        let mesh = unit_square_mesh();
        let domain_map = DomainMap::from_materials(vec![Material::default()], [(1u32, 0usize)]);
        let phi: Vec<f64> = mesh.nodes.iter().map(|node| node.y).collect();
        let total = postprocess::electrostatic_energy(&phi, &mesh, |_| EPS0);
        let parts = peak_domain_energy_records(&mesh, &domain_map, &phi, total);

        assert_eq!(parts.len(), 2);
        let summed: f64 = parts.iter().map(|record| record.energy).sum();
        assert!((summed - total).abs() < 1e-30, "summed={summed:.6e}, total={total:.6e}");
        assert_eq!(parts[0].material_index, Some(0));
        assert_eq!(parts[1].material_index, None);
    }
}
