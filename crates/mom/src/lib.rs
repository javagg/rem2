//! rem-mom �?Method of Moments solver for rem2
//!
//! Implements full-wave electromagnetic scattering and radiation using the
//! Method of Moments (MoM) with RWG basis functions and CFIE formulation.
//!
//! # Architecture
//! ```text
//! RemMesh  �? SurfaceMesh  �? RwgBases
//!                  �?              �?
//!             quadrature      assemble Z (dense matrix)
//!                  �?              �?
//!             singular.rs    LU solve �?currents
//!                                  �?
//!                            postprocess �?RCS CSV
//! ```

pub mod surface_mesh;
pub mod quadrature;
pub mod green;
pub mod green_wrapper;
pub mod singular;
pub mod assemble;
pub mod excitation;
pub mod postprocess;
pub mod basis;
pub mod mie;
pub mod gpu;
pub mod aca;
pub mod pmchwt;
pub mod loop_star;
pub mod port;
pub mod sparams;
pub mod sibc;
pub mod fft_accel;
pub mod fmm;
pub mod mlfma;
pub mod amr;
pub mod rom;

// Public re-exports for cross-crate solver integration
pub use assemble::{gmres_solve, gmres_solve_generic, gmres_solve_op, aca_gmres_solve, gmres_generic_with_aca};

use rem_config::{PalaceConfig, MomSolverConfig};

/// Default multipole order for MLFMA (matches `mlfma::DEFAULT_P`).
const DEFAULT_MLFMA_P: usize = 6;
use rem_core::RemResult;
use rem_mesh::RemMesh;
use rem_parallel::NoComm;
use rem_layered_green::{GreenFunction, FreeSpaceGreen, LayeredGreen, DielectricLayer};

/// Build a boxed [`GreenFunction`] from the MoM solver config at a given frequency.
///
/// If `mom_cfg.substrate` is set, returns a [`LayeredGreen`]; otherwise free-space.
fn build_green(mom_cfg: &MomSolverConfig, freq: f64) -> Box<dyn GreenFunction> {
    use std::f64::consts::PI;
    let k0 = 2.0 * PI * freq / rem_core::C0;
    if let Some(sub) = &mom_cfg.substrate {
        let layers: Vec<DielectricLayer> = sub.layers.iter().map(|l| DielectricLayer {
            eps_r: l.permittivity,
            loss_tan: l.loss_tangent,
            mu_r: l.permeability,
            thickness_m: l.thickness,
            // Frequency-dependent or anisotropic lateral ε: evaluated at this frequency
            eps_r_complex_override: Some(l.eps_r_complex(freq)),
            // Anisotropic vertical ε_zz (None = isotropic)
            eps_r_z: l.eps_r_z_complex(freq),
        }).collect();
        log::info!(
            "MoM: using layered Green function ({} layer(s), bottom_pec={}, f={:.3e} Hz)",
            layers.len(), sub.bottom_pec, freq
        );
        Box::new(LayeredGreen::new(layers, k0, 0.0, false))
    } else {
        Box::new(FreeSpaceGreen::new(k0))
    }
}

/// One observation angle's RCS result.
#[derive(Debug, Clone)]
pub struct RcsPoint {
    pub theta_deg: f64,
    pub phi_deg:   f64,
    pub rcs_m2:    f64,
    /// RCS in dBsm = 10·log10(rcs_m2); -300 if rcs_m2 �?0
    pub rcs_dbsm:  f64,
}

/// Result returned by `run_with_mesh`.
#[derive(Debug, Clone)]
pub struct MomResult {
    /// Per-frequency RCS pattern data
    pub rcs: Vec<(f64, Vec<RcsPoint>)>,   // (freq_hz, points)
}

/// Entry point called from the CLI for `Problem.Type = "MoM"`.
pub fn run(config: &PalaceConfig) -> RemResult<()> {
    let mom_cfg = config.solver.mom.as_ref()
        .ok_or_else(|| rem_core::RemError::Config(
            "Problem.Type = \"MoM\" requires a Solver.MoM section".to_string()
        ))?;

    let mesh = rem_mesh::load_mesh(config, &NoComm)?;
    run_with_mesh(config, mom_cfg, &mesh).map(|_| ())
}

/// Run MoM solve on an already-loaded mesh (also used from WASM / tests).
pub fn run_with_mesh(
    config: &PalaceConfig,
    mom_cfg: &MomSolverConfig,
    mesh: &RemMesh,
) -> RemResult<MomResult> {
    use std::f64::consts::PI;

    log::info!("\n=== Method of Moments (MoM) solver ===\n");

    // Collect PEC surface attribute IDs
    let pec_attrs: Vec<u32> = config.boundaries.pec
        .as_ref()
        .map(|p| p.attributes.clone())
        .unwrap_or_default();

    if pec_attrs.is_empty() {
        return Err(rem_core::RemError::Config(
            "MoM solver requires at least one PEC boundary (Boundaries.PEC.Attributes)".to_string()
        ));
    }

    // Build surface mesh
    let surf = surface_mesh::SurfaceMesh::extract(mesh, &pec_attrs)?;
    log::info!("Surface mesh (RWG basis):");
    log::info!("  {} triangular faces", surf.faces.len());
    log::info!("  {} interior edges (basis functions)", surf.edges.len());
    log::info!("");

    // ── Port path (S-parameter sweep) ─────────────────────────────────────
    if !mom_cfg.ports.is_empty() {
        return run_s_param_sweep(config, mom_cfg, &surf);
    }

    // ── RCS path (original plane-wave path) ───────────────────────────────

    // Frequency sweep
    let freq_min  = mom_cfg.freq_min;
    let freq_max  = mom_cfg.freq_max;
    let freq_step = mom_cfg.freq_step;
    let output_dir = std::path::Path::new(config.problem.output_dir());
    #[cfg(not(target_arch = "wasm32"))]
    std::fs::create_dir_all(output_dir.join("postpro"))?;

    // RCS angles
    let (theta_deg, phi_deg) = if let Some(rcs) = &config.postprocessing.rcs {
        (rcs.theta_deg.clone(), rcs.phi_deg.clone())
    } else {
        // Default: bistatic cut in E-plane
        let theta: Vec<f64> = (0..=180).step_by(5).map(|i| i as f64).collect();
        (theta, vec![0.0])
    };

    // Gaussian quadrature (7-point, degree 5)
    let quad = quadrature::TriQuad::new(5);

    let mut freq = freq_min;
    let mut all_rcs: Vec<(f64, Vec<RcsPoint>)> = Vec::new();

    while freq <= freq_max + 1e-3 * freq_step {
        log::info!("MoM solve at f = {:.3e} Hz", freq);

        let k = 2.0 * PI * freq / rem_core::C0;

        // Incident plane wave OR near-field source
        let rhs = if let Some(nf_path) = &mom_cfg.near_field_source {
            let nf_path = std::path::Path::new(nf_path);
            log::info!("MoM: loading near-field source from {}", nf_path.display());
            let nf_points = rem_core::read_near_field_csv(nf_path)?;
            log::info!("MoM: loaded {} near-field points", nf_points.len());
            excitation::near_field_rhs(&surf, k, &nf_points, &mom_cfg.basis)
        } else {
            let wave = excitation::PlaneWave {
                theta_inc: mom_cfg.theta_inc_deg.to_radians(),
                phi_inc:   mom_cfg.phi_inc_deg.to_radians(),
                pol:       mom_cfg.polarization.clone(),
            };
            excitation::plane_wave_rhs_general(&surf, k, &wave, &mom_cfg.basis)
        };

        // PMCHWT path: dielectric target (J + M unknowns, 2N×2N system)
        let currents = if mom_cfg.equation.to_uppercase() == "PMCHWT" {
            let (eps_r, mu_r) = config.domains.materials.first()
                .map(|m| (m.permittivity, m.permeability))
                .unwrap_or((2.0, 1.0));
            let mat = pmchwt::DielectricMaterial::new(eps_r, mu_r);
            log::info!("MoM PMCHWT: ε_r={eps_r:.2}, μ_r={mu_r:.2}, f={freq:.3e} Hz");
            // PMCHWT path: for now use plane wave excitation only
            let wave = excitation::PlaneWave {
                theta_inc: mom_cfg.theta_inc_deg.to_radians(),
                phi_inc:   mom_cfg.phi_inc_deg.to_radians(),
                pol:       mom_cfg.polarization.clone(),
            };
            if mom_cfg.near_field_source.is_some() {
                log::warn!("MoM: NearFieldSource is ignored for PMCHWT equation; using plane wave");
            }
            let (j_coeffs, _m_coeffs) = pmchwt::solve_pmchwt(
                &surf, mat, freq, &wave, &quad, &mom_cfg.fast_solver,
            )?;
            j_coeffs
        } else {
            // Assemble impedance matrix Z  (PEC EFIE/MFIE/CFIE path)
            // Early-exit for FMM/MLFMA: builds a matrix-free operator, skips full Z assembly.
            let fast = mom_cfg.fast_solver.to_uppercase();
            let is_fmm   = fast == "FMM"   && !mom_cfg.basis.eq_ignore_ascii_case("Pulse");
            let is_mlfma = fast == "MLFMA" && !mom_cfg.basis.eq_ignore_ascii_case("Pulse");
            if is_mlfma {
                let bases = basis::rwg::generate_rwg_bases(&surf);
                log::info!(
                    "MoM MLFMA: building multilevel FMM operator (N={}, P={})",
                    bases.len(), DEFAULT_MLFMA_P
                );
                let green_ml = build_green(mom_cfg, freq);
                let quad_ml  = quadrature::TriQuad::new(3);
                let rhs_dv = nalgebra::DVector::from_vec(rhs.clone());
                mlfma::mlfma_solve(
                    &surf, &bases, green_ml.as_ref(),
                    freq, mom_cfg.alpha, &quad_ml, &rhs_dv, DEFAULT_MLFMA_P,
                )?.as_slice().to_vec()
            } else if is_fmm {
                let bases = basis::rwg::generate_rwg_bases(&surf);
                log::info!("MoM FMM: building 3-D FFT monopole FMM (N={})", bases.len());
                let green_fmm = build_green(mom_cfg, freq);
                let quad_fmm  = quadrature::TriQuad::new(3);
                let fmm_op = fmm::FmmMomSolver::build(
                    &surf, &bases, green_fmm.as_ref(),
                    freq, mom_cfg.alpha, &quad_fmm,
                )?;
                let rhs_dv = nalgebra::DVector::from_vec(rhs.clone());
                assemble::gmres_solve_op(&fmm_op, &rhs_dv)?.as_slice().to_vec()
            } else {
            let z_mat = match mom_cfg.basis.as_str() {
                "Pulse" | "pulse" => {
                    let mut z = assemble::assemble_efie_pulse(&surf, freq, &quad, mom_cfg.singular_tol)?;
                    if mom_cfg.wall_conductivity > 0.0 {
                        sibc::apply_sibc_pulse(&mut z, &surf, freq, mom_cfg.wall_conductivity);
                        log::info!(
                            "MoM SIBC (Pulse): sigma_wall={:.3e} S/m",
                            mom_cfg.wall_conductivity
                        );
                    }
                    Ok::<nalgebra::DMatrix<num_complex::Complex64>, rem_core::RemError>(z)
                }
                _ => {
                    let bases = basis::rwg::generate_rwg_bases(&surf);
                    let green = build_green(mom_cfg, freq);
                    let mut z = assemble::assemble_cfie_rwg_green(
                        &surf, &bases, green.as_ref(), freq, mom_cfg.alpha, &quad, mom_cfg.singular_tol,
                    )?;
                    if mom_cfg.wall_conductivity > 0.0 {
                        sibc::apply_sibc_rwg(&mut z, &surf, &bases, freq, mom_cfg.wall_conductivity, &quad);
                        log::info!(
                            "MoM SIBC enabled: sigma_wall={:.3e} S/m",
                            mom_cfg.wall_conductivity
                        );
                    }
                    Ok(z)
                }
            }?;
            match mom_cfg.fast_solver.to_uppercase().as_str() {
                "GMRES" => assemble::gmres_solve(&z_mat, &rhs)?,
                "ACA" => {
                    log::info!("MoM: using ACA+GMRES (tol_aca=1e-4, tol_gmres=1e-8)");
                    assemble::aca_gmres_solve(&z_mat, &rhs, 1e-4, 1e-8)?
                }
                "FFT" => {
                    if fft_accel::FftMomSolver::is_applicable(&surf.nodes) {
                        log::info!("MoM FFT: planar mesh detected, building FFT operator (N={})", surf.nodes.len());
                        let fft_op = fft_accel::FftMomSolver::build(&surf.nodes, k)?;
                        let rhs_dv = nalgebra::DVector::from_vec(rhs.clone());
                        assemble::gmres_solve_op(&fft_op, &rhs_dv)?.as_slice().to_vec()
                    } else {
                        log::warn!("MoM FFT: mesh is not planar �?falling back to GMRES with dense matrix");
                        assemble::gmres_solve(&z_mat, &rhs)?
                    }
                }
                _ => assemble::lu_solve(&z_mat, &rhs)?,
            }
            } // end else (non-FMM)
        };

        // Compute RCS pattern (always, not just for file output)
        let rcs_grid = postprocess::rcs_pattern(&currents, &surf, k, &theta_deg, &phi_deg);
        let mut pts = Vec::new();
        for (ti, &th) in theta_deg.iter().enumerate() {
            for (pi, &ph) in phi_deg.iter().enumerate() {
                let rcs_m2 = rcs_grid[ti][pi];
                let rcs_dbsm = if rcs_m2 > 1e-300 { 10.0 * rcs_m2.log10() } else { -300.0 };
                pts.push(RcsPoint { theta_deg: th, phi_deg: ph, rcs_m2, rcs_dbsm });
            }
        }
        log::info!("  RCS computed: {} angle pairs", pts.len());
        all_rcs.push((freq, pts));

        #[cfg(not(target_arch = "wasm32"))]
        {
            postprocess::write_rcs(output_dir, freq, &currents, &surf, k, &theta_deg, &phi_deg)?;
            let vtk_path = output_dir
                .join("postpro")
                .join(format!("surface_current_{:.3e}Hz.vtk", freq));
            postprocess::write_surface_vtk(&vtk_path, &currents, &surf)?;

            // Near-field export (if configured)
            if let Some(nf_cfg) = &config.postprocessing.near_field {
                let nf_points = postprocess::compute_near_field(&currents, &surf, k);
                postprocess::write_near_field_csv(output_dir, &nf_points, nf_cfg.output_file.as_deref())?;
            }
        }

        freq += freq_step;
    }

    #[cfg(not(target_arch = "wasm32"))]
    log::info!("MoM solve complete. Results in {}", output_dir.display());
    #[cfg(target_arch = "wasm32")]
    log::info!("MoM solve complete.");
    Ok(MomResult { rcs: all_rcs })
}

/// Run the S-parameter sweep (port-excited MoM path).
fn run_s_param_sweep(
    config: &PalaceConfig,
    mom_cfg: &MomSolverConfig,
    surf: &surface_mesh::SurfaceMesh,
) -> RemResult<MomResult> {
    use port::MomLumpedPort;

    let output_dir = std::path::Path::new(config.problem.output_dir());
    #[cfg(not(target_arch = "wasm32"))]
    std::fs::create_dir_all(output_dir.join("postpro"))?;

    let bases = basis::rwg::generate_rwg_bases(surf);
    let quad  = quadrature::TriQuad::new(5);

    // Build physical ports (Lumped/WavePort currently share the same MoM
    // surface-excitation kernel; WavePort mode metadata is preserved for
    // output/de-embedding bookkeeping).
    let lumped_ports: Vec<MomLumpedPort> = mom_cfg.ports.iter().map(|p| {
        let z0 = p.impedance.unwrap_or(mom_cfg.ref_impedance);
        MomLumpedPort::from_surface(
            surf,
            &bases,
            &p.attributes,
            p.index,
            &p.direction,
            &p.port_type,
            p.mode,
            z0,
        )
    }).collect::<RemResult<_>>()?;

    for (cfg_p, port_p) in mom_cfg.ports.iter().zip(lumped_ports.iter()) {
        if cfg_p.port_type.eq_ignore_ascii_case("waveport") {
            if port_p.modal_profile.is_some() {
                log::info!(
                    "MoM WavePort {} mode {}: modal profile enabled",
                    cfg_p.index,
                    cfg_p.mode
                );
            } else {
                log::warn!(
                    "MoM WavePort {} mode {}: modal solve unavailable on current port patch, fallback to uniform excitation",
                    cfg_p.index,
                    cfg_p.mode
                );
            }
        }
    }

    // Frequency sweep �?use ROM if RomOrder > 0
    let freq = mom_cfg.freq_min;
    let freq_max  = mom_cfg.freq_max;
    let freq_step = mom_cfg.freq_step;
    let mut all_matrices: Vec<sparams::SMatrix> = Vec::new();

    // Collect all sweep frequencies up front (needed for ROM)
    let mut freq_list: Vec<f64> = Vec::new();
    let mut f = freq;
    while f <= freq_max + 1e-3 * freq_step {
        freq_list.push(f);
        f += freq_step;
    }

    if mom_cfg.rom_order > 0 && freq_list.len() > mom_cfg.rom_order {
        // ── ROM-accelerated sweep ─────────────────────────────────────────
        log::info!(
            "MoM ROM sweep: {} anchor points over {} frequencies",
            mom_cfg.rom_order, freq_list.len()
        );
        let alpha   = mom_cfg.alpha;
        let sing_tol = mom_cfg.singular_tol;
        let sigma   = mom_cfg.wall_conductivity;
        let build_z = |fq: f64| -> RemResult<nalgebra::DMatrix<num_complex::Complex64>> {
            let green = build_green(mom_cfg, fq);
            let mut z = assemble::assemble_cfie_rwg_green(
                surf, &bases, green.as_ref(), fq, alpha, &quad, sing_tol,
            )?;
            if sigma > 0.0 {
                sibc::apply_sibc_rwg(&mut z, surf, &bases, fq, sigma, &quad);
            }
            Ok(z)
        };
        all_matrices = rom::mom_rom_sweep(
            surf, &bases, &lumped_ports,
            &freq_list, mom_cfg.rom_order, 1e-10,
            &build_z,
        )?;
    } else {
        // ── Direct sweep ──────────────────────────────────────────────────
        for &freq in &freq_list {
            log::info!("MoM S-param solve at f = {:.3e} Hz", freq);

            let green = build_green(mom_cfg, freq);

            let sm = if mom_cfg.fast_solver.eq_ignore_ascii_case("FFT")
                && fft_accel::FftMomSolver::is_applicable(&surf.nodes)
            {
                // FFT-accelerated path: avoid building full Z matrix
                let k = 2.0 * std::f64::consts::PI * freq / rem_core::C0;
                log::info!("MoM FFT S-param: planar mesh N={}, building FFT operator", surf.nodes.len());
                let fft_op = fft_accel::FftMomSolver::build(&surf.nodes, k)?;
                sparams::compute_s_matrix_op(surf, &bases, &lumped_ports, &fft_op, freq)?
            } else if mom_cfg.fast_solver.eq_ignore_ascii_case("FMM") {
                let k = 2.0 * std::f64::consts::PI * freq / rem_core::C0;
                let quad_fmm = quadrature::TriQuad::new(3);
                log::info!("MoM FMM S-param: building 3-D FFT monopole FMM (N={})", bases.len());
                let fmm_op = fmm::FmmMomSolver::build(
                    surf, &bases, green.as_ref(), freq, mom_cfg.alpha, &quad_fmm,
                )?;
                let _ = k;
                sparams::compute_s_matrix_op(surf, &bases, &lumped_ports, &fmm_op, freq)?
            } else {
                let z_mat = assemble::assemble_cfie_rwg_green(
                    surf, &bases, green.as_ref(), freq, mom_cfg.alpha, &quad, mom_cfg.singular_tol,
                )?;
                let mut z_mat = z_mat;
                if mom_cfg.wall_conductivity > 0.0 {
                    sibc::apply_sibc_rwg(&mut z_mat, surf, &bases, freq, mom_cfg.wall_conductivity, &quad);
                    log::info!(
                        "MoM SIBC enabled: sigma_wall={:.3e} S/m",
                        mom_cfg.wall_conductivity
                    );
                }
                sparams::compute_s_matrix(surf, &bases, &lumped_ports, &z_mat, freq)?
            };

            log::info!("  S-matrix computed: {}×{}", sm.n_ports, sm.n_ports);
            all_matrices.push(sm);
        }
    }

    // Optional per-port de-embedding before writing outputs.
    let deembed_lengths: Vec<f64> = mom_cfg.ports.iter().map(|p| p.deembed_length).collect();
    if deembed_lengths.iter().any(|&l| l.abs() > 0.0) {
        // Determine whether any port is a WavePort with a modal eigenvalue.
        // When available, use precise modal de-embedding; otherwise fall back
        // to scalar phase-delay approximation.
        let has_modal = lumped_ports.iter().any(|p| p.modal_eigenvalue.is_some());
        if has_modal {
            let modal_data: Vec<sparams::ModalPortData> = lumped_ports.iter()
                .map(|p| p.compute_modal_data(
                    surf,
                    freq, // use first sweep frequency for de-embedding γ
                    mom_cfg.deembed_eps_eff,
                    mom_cfg.deembed_alpha_np_per_m,
                ))
                .collect();
            all_matrices = all_matrices.iter()
                .map(|s| sparams::apply_modal_deembed(s, &deembed_lengths, &modal_data))
                .collect::<RemResult<_>>()?;
            log::info!(
                "MoM WavePort modal de-embedding: {} ports with eigenvalue-derived γ",
                lumped_ports.iter().filter(|p| p.modal_eigenvalue.is_some()).count()
            );
        } else {
            all_matrices = all_matrices.iter()
                .map(|s| sparams::apply_reference_plane_deembed(
                    s,
                    &deembed_lengths,
                    mom_cfg.deembed_eps_eff,
                    mom_cfg.deembed_alpha_np_per_m,
                ))
                .collect::<RemResult<_>>()?;
            log::info!(
                "MoM reference-plane de-embedding enabled: eps_eff={:.4}, alpha={:.4e} Np/m",
                mom_cfg.deembed_eps_eff,
                mom_cfg.deembed_alpha_np_per_m,
            );
        }
    }

    // Build differential-pair index list (0-based positions in `lumped_ports`).
    let mut differential_pairs: Vec<(usize, usize)> = Vec::new();
    if !mom_cfg.ports.is_empty() {
        let mut by_index = std::collections::HashMap::<u32, usize>::new();
        for (i, p) in mom_cfg.ports.iter().enumerate() {
            by_index.insert(p.index, i);
        }
        let mut paired_member = vec![false; mom_cfg.ports.len()];
        for (i, p) in mom_cfg.ports.iter().enumerate() {
            if let Some(pair_idx) = p.pair_with {
                if let Some(&j) = by_index.get(&pair_idx) {
                    let reciprocal = mom_cfg.ports[j].pair_with == Some(p.index);
                    if i < j && reciprocal {
                        differential_pairs.push((i, j));
                        paired_member[i] = true;
                        paired_member[j] = true;
                    }
                }
            }
        }
        for (i, p) in mom_cfg.ports.iter().enumerate() {
            if p.pair_with.is_some() && !paired_member[i] {
                log::warn!(
                    "MoM differential pair ignored for port {}: PairWith requires reciprocal reference",
                    p.index
                );
            }
        }
    }

    // Write outputs
    #[cfg(not(target_arch = "wasm32"))]
    {
        let ts_ext = format!("s{}p", lumped_ports.len());
        let ts_path = output_dir.join("postpro").join(format!("s_params.{}", ts_ext));
        sparams::write_touchstone(&all_matrices, &ts_path, mom_cfg.ref_impedance)?;
        log::info!("MoM S-param output: {}", ts_path.display());

        let csv_path = output_dir.join("postpro").join("port-S.csv");
        sparams::append_palace_csv(&all_matrices, &csv_path)?;

        // Mixed-mode differential output.
        if !differential_pairs.is_empty() {
            // Full mixed-mode matrix when all ports are paired.
            if differential_pairs.len() * 2 == lumped_ports.len() {
                let mixed_matrices: Vec<sparams::SMatrix> = all_matrices.iter()
                    .map(|s| sparams::single_ended_to_mixed_mode(s, &differential_pairs))
                    .collect::<RemResult<_>>()?;
                let mm_ts_ext = format!("s{}p", lumped_ports.len());
                let mm_ts_path = output_dir.join("postpro").join(format!("s_params_mixed.{}", mm_ts_ext));
                sparams::write_touchstone(&mixed_matrices, &mm_ts_path, mom_cfg.ref_impedance)?;
                sparams::append_palace_csv(&mixed_matrices, &output_dir.join("postpro").join("port-S-mixed.csv"))?;
                log::info!(
                    "MoM full mixed-mode S output written for {} differential pairs",
                    differential_pairs.len()
                );
            } else {
                log::info!(
                    "MoM partial mixed-mode: {} paired ports over {} total ports",
                    differential_pairs.len() * 2,
                    lumped_ports.len()
                );
            }

            // Per-pair 2x2 mixed-mode blocks are always emitted for each valid pair.
            for &(pi, ni) in &differential_pairs {
                let pair_idx_p = mom_cfg.ports[pi].index;
                let pair_idx_n = mom_cfg.ports[ni].index;
                let pair_matrices: Vec<sparams::SMatrix> = all_matrices.iter()
                    .map(|s| sparams::pair_mixed_mode_block(s, (pi, ni)))
                    .collect::<RemResult<_>>()?;

                let pair_ts_path = output_dir.join("postpro").join(format!(
                    "s_params_mixed_p{}_p{}.s2p",
                    pair_idx_p, pair_idx_n,
                ));
                sparams::write_touchstone(&pair_matrices, &pair_ts_path, mom_cfg.ref_impedance)?;

                let pair_csv_path = output_dir.join("postpro").join(format!(
                    "port-S-mixed-p{}_p{}.csv",
                    pair_idx_p, pair_idx_n,
                ));
                sparams::append_palace_csv(&pair_matrices, &pair_csv_path)?;
            }
        }

        // ── Z/Y matrices ──────────────────────────────────────────────────
        let z_matrices: Vec<sparams::SMatrix> = all_matrices.iter()
            .map(|s| sparams::s_to_z(s, mom_cfg.ref_impedance))
            .collect::<RemResult<_>>()?;
        sparams::write_param_csv(&z_matrices,
            &output_dir.join("postpro").join("port-Z.csv"), "Z")?;

        let y_matrices: Vec<sparams::SMatrix> = z_matrices.iter()
            .map(sparams::z_to_y)
            .collect::<RemResult<_>>()?;
        sparams::write_param_csv(&y_matrices,
            &output_dir.join("postpro").join("port-Y.csv"), "Y")?;
        log::info!("MoM Z/Y-param output written");

        // ── Transmission-line RLGC (2-port, if tline_length > 0) ─────────
        if mom_cfg.tline_length > 0.0 && lumped_ports.len() == 2 {
            let tline = sparams::extract_tline_rlgc(
                &all_matrices, mom_cfg.ref_impedance, mom_cfg.tline_length,
            );
            sparams::write_tline_csv(&tline,
                &output_dir.join("postpro").join("tline_params.csv"))?;
            log::info!("MoM RLGC output: {} freq points (�?{:.4e} m)",
                tline.len(), mom_cfg.tline_length);
        }

        // ── Near-field probes (all ports) ─────────────────────────────────
        if !mom_cfg.near_field_probes.is_empty() {
            let probes: Vec<[f64; 3]> = mom_cfg.near_field_probes.iter()
                .map(|p| [p.x, p.y, p.z])
                .collect();
            let bases = basis::rwg::generate_rwg_bases(surf);
            let quad  = quadrature::TriQuad::new(5);

            // Re-solve for each port's excitation at every frequency.
            // Output: probe_e_field_portN.csv  (N = 1-based port index).
            // For single-port (or when only port 1 needed), also write the
            // legacy probe_e_field.csv for backward compatibility.
            for (port_idx, port) in lumped_ports.iter().enumerate() {
                let port_label = port.index;
                let mut freq_e: Vec<(f64, Vec<[num_complex::Complex64; 3]>)> = Vec::new();
                for &fq in &freq_list {
                    let green = build_green(mom_cfg, fq);
                    let mut z = assemble::assemble_cfie_rwg_green(
                        surf, &bases, green.as_ref(), fq, mom_cfg.alpha, &quad, mom_cfg.singular_tol,
                    )?;
                    if mom_cfg.wall_conductivity > 0.0 {
                        sibc::apply_sibc_rwg(&mut z, surf, &bases, fq, mom_cfg.wall_conductivity, &quad);
                    }
                    let n_rwg = bases.len();
                    let v0 = num_complex::Complex64::new(1.0, 0.0);
                    let rhs = port.excitation_rhs(surf, &bases, n_rwg, v0);
                    let currents = assemble::lu_solve(&z, &rhs)?;
                    let k = 2.0 * std::f64::consts::PI * fq / rem_core::C0;
                    let e_vals = postprocess::compute_e_at_probes(surf, &bases, &currents, &probes, k);
                    freq_e.push((fq, e_vals));
                }
                let csv_name = if lumped_ports.len() == 1 {
                    "probe_e_field.csv".to_string()
                } else {
                    format!("probe_e_field_port{}.csv", port_label)
                };
                postprocess::write_probe_e_field_csv(
                    &output_dir.join("postpro").join(&csv_name),
                    &probes, &freq_e,
                )?;
                // Backward-compatibility alias for port 1
                if port_idx == 0 && lumped_ports.len() > 1 {
                    postprocess::write_probe_e_field_csv(
                        &output_dir.join("postpro").join("probe_e_field.csv"),
                        &probes, &freq_e,
                    )?;
                }
            }
            log::info!("MoM near-field probe output: {} probes × {} frequencies × {} ports",
                probes.len(), freq_list.len(), lumped_ports.len());
        }

        // ── Far-field radiation pattern (per port, if FarField config present) ─
        if let Some(ff_cfg) = &config.solver.far_field {
            let bases = basis::rwg::generate_rwg_bases(surf);
            let quad  = quadrature::TriQuad::new(5);

            let n_theta = ff_cfg.n_theta.max(2);
            let n_phi   = ff_cfg.n_phi.max(2);
            let theta_list: Vec<f64> = (0..n_theta)
                .map(|i| i as f64 * 180.0 / (n_theta - 1) as f64)
                .collect();
            let phi_list: Vec<f64> = (0..n_phi)
                .map(|i| i as f64 * 360.0 / n_phi as f64)
                .collect();

            log::info!(
                "MoM far-field pattern: {} θ × {} φ points per port",
                n_theta, n_phi
            );

            for (port_idx, port) in lumped_ports.iter().enumerate() {
                let port_label = port.index;
                // Use last-frequency solve for the far-field pattern (or
                // sweep over all frequencies if needed).  Here we compute at
                // every sweep frequency and append to CSV.
                let csv_name = if lumped_ports.len() == 1 {
                    "far_field.csv".to_string()
                } else {
                    format!("far_field_port{}.csv", port_label)
                };
                let ff_csv = output_dir.join("postpro").join(&csv_name);
                // Remove stale file so write_radiation_pattern_csv writes a
                // fresh header on the first frequency.
                let _ = std::fs::remove_file(&ff_csv);

                for &fq in &freq_list {
                    let green = build_green(mom_cfg, fq);
                    let mut z = assemble::assemble_cfie_rwg_green(
                        surf, &bases, green.as_ref(), fq, mom_cfg.alpha, &quad, mom_cfg.singular_tol,
                    )?;
                    if mom_cfg.wall_conductivity > 0.0 {
                        sibc::apply_sibc_rwg(&mut z, surf, &bases, fq, mom_cfg.wall_conductivity, &quad);
                    }
                    let n_rwg = bases.len();
                    let v0 = num_complex::Complex64::new(1.0, 0.0);
                    let rhs = port.excitation_rhs(surf, &bases, n_rwg, v0);
                    let currents = assemble::lu_solve(&z, &rhs)?;
                    let k = 2.0 * std::f64::consts::PI * fq / rem_core::C0;

                    let pattern = postprocess::compute_radiation_pattern_rwg(
                        &currents, surf, &bases, k, &theta_list, &phi_list,
                    );
                    postprocess::write_radiation_pattern_csv(&ff_csv, &pattern, fq)?;

                    // Per-port surface current VTK at each frequency
                    let vtk_name = if lumped_ports.len() == 1 {
                        format!("surface_current_{:.3e}Hz.vtk", fq)
                    } else {
                        format!("surface_current_port{}_{:.3e}Hz.vtk", port_label, fq)
                    };
                    let vtk_path = output_dir.join("postpro").join(vtk_name);
                    postprocess::write_surface_current_vtk_rwg(&vtk_path, surf, &bases, &currents)?;
                }
                log::info!(
                    "MoM far-field pattern port {}: {} frequencies �?{}",
                    port_label, freq_list.len(), ff_csv.display()
                );
                // backward-compat alias for port 1
                if port_idx == 0 && lumped_ports.len() > 1 {
                    let _ = std::fs::remove_file(&output_dir.join("postpro").join("far_field.csv"));
                    std::fs::copy(&ff_csv, output_dir.join("postpro").join("far_field.csv")).ok();
                }
            }
        }
    }

    log::info!("MoM S-param sweep complete. {} frequency points.", all_matrices.len());
    Ok(MomResult { rcs: vec![] })  // no RCS in S-param mode
}

/// Compute S-parameter sweep without writing any output files.
///
/// Used by `rem-optim` for parametric sweeps and gradient optimization.
/// Returns one [`sparams::SMatrix`] per frequency point in ascending order.
///
/// The caller must supply:
/// - a [`PalaceConfig`] with `Solver.MoM.Ports` populated (non-empty), and
/// - the already-loaded [`RemMesh`] for the geometry.
///
/// De-embedding (if configured via `DeembedLength`) is applied on the returned matrices.
pub fn compute_s_param_sweep_for_optim(
    config: &PalaceConfig,
    mom_cfg: &MomSolverConfig,
    mesh: &RemMesh,
) -> RemResult<Vec<sparams::SMatrix>> {
    use port::MomLumpedPort;

    let pec_attrs: Vec<u32> = config.boundaries.pec
        .as_ref()
        .map(|p| p.attributes.clone())
        .unwrap_or_default();
    let surf = surface_mesh::SurfaceMesh::extract(mesh, &pec_attrs)?;
    let bases = basis::rwg::generate_rwg_bases(&surf);
    let quad  = quadrature::TriQuad::new(5);

    let lumped_ports: Vec<MomLumpedPort> = mom_cfg.ports.iter().map(|p| {
        let z0 = p.impedance.unwrap_or(mom_cfg.ref_impedance);
        MomLumpedPort::from_surface(
            &surf, &bases, &p.attributes, p.index,
            &p.direction, &p.port_type, p.mode, z0,
        )
    }).collect::<RemResult<_>>()?;

    let mut freq_list: Vec<f64> = Vec::new();
    let mut f = mom_cfg.freq_min;
    while f <= mom_cfg.freq_max + 1e-3 * mom_cfg.freq_step {
        freq_list.push(f);
        f += mom_cfg.freq_step;
    }

    let mut all_matrices: Vec<sparams::SMatrix> = Vec::new();

    if mom_cfg.rom_order > 0 && freq_list.len() > mom_cfg.rom_order {
        let alpha    = mom_cfg.alpha;
        let sing_tol = mom_cfg.singular_tol;
        let sigma    = mom_cfg.wall_conductivity;
        let build_z = |fq: f64| -> RemResult<nalgebra::DMatrix<num_complex::Complex64>> {
            let green = build_green(mom_cfg, fq);
            let mut z = assemble::assemble_cfie_rwg_green(
                &surf, &bases, green.as_ref(), fq, alpha, &quad, sing_tol,
            )?;
            if sigma > 0.0 {
                sibc::apply_sibc_rwg(&mut z, &surf, &bases, fq, sigma, &quad);
            }
            Ok(z)
        };
        all_matrices = rom::mom_rom_sweep(
            &surf, &bases, &lumped_ports,
            &freq_list, mom_cfg.rom_order, 1e-10, &build_z,
        )?;
    } else {
        for &freq in &freq_list {
            let green = build_green(mom_cfg, freq);
            let sm = if mom_cfg.fast_solver.eq_ignore_ascii_case("FFT")
                && fft_accel::FftMomSolver::is_applicable(&surf.nodes)
            {
                let k = 2.0 * std::f64::consts::PI * freq / rem_core::C0;
                let fft_op = fft_accel::FftMomSolver::build(&surf.nodes, k)?;
                sparams::compute_s_matrix_op(&surf, &bases, &lumped_ports, &fft_op, freq)?
            } else if mom_cfg.fast_solver.eq_ignore_ascii_case("FMM") {
                let quad_fmm = quadrature::TriQuad::new(3);
                log::info!("MoM FMM S-param (ROM path): building FMM (N={})", bases.len());
                let fmm_op = fmm::FmmMomSolver::build(
                    &surf, &bases, green.as_ref(), freq, mom_cfg.alpha, &quad_fmm,
                )?;
                sparams::compute_s_matrix_op(&surf, &bases, &lumped_ports, &fmm_op, freq)?
            } else {
                let mut z_mat = assemble::assemble_cfie_rwg_green(
                    &surf, &bases, green.as_ref(), freq, mom_cfg.alpha, &quad, mom_cfg.singular_tol,
                )?;
                if mom_cfg.wall_conductivity > 0.0 {
                    sibc::apply_sibc_rwg(&mut z_mat, &surf, &bases, freq, mom_cfg.wall_conductivity, &quad);
                }
                sparams::compute_s_matrix(&surf, &bases, &lumped_ports, &z_mat, freq)?
            };
            all_matrices.push(sm);
        }
    }

    // De-embedding (mirrors run_s_param_sweep logic)
    let deembed_lengths: Vec<f64> = mom_cfg.ports.iter().map(|p| p.deembed_length).collect();
    if deembed_lengths.iter().any(|&l| l.abs() > 0.0) {
        let has_modal = lumped_ports.iter().any(|p| p.modal_eigenvalue.is_some());
        if has_modal {
            let first_freq = freq_list.first().copied().unwrap_or(mom_cfg.freq_min);
            let modal_data: Vec<sparams::ModalPortData> = lumped_ports.iter()
                .map(|p| p.compute_modal_data(&surf, first_freq,
                    mom_cfg.deembed_eps_eff, mom_cfg.deembed_alpha_np_per_m))
                .collect();
            all_matrices = all_matrices.iter()
                .map(|s| sparams::apply_modal_deembed(s, &deembed_lengths, &modal_data))
                .collect::<RemResult<_>>()?;
        } else {
            all_matrices = all_matrices.iter()
                .map(|s| sparams::apply_reference_plane_deembed(
                    s, &deembed_lengths,
                    mom_cfg.deembed_eps_eff, mom_cfg.deembed_alpha_np_per_m,
                ))
                .collect::<RemResult<_>>()?;
        }
    }

    Ok(all_matrices)
}

