//! rem-mom — Method of Moments solver for rem2
//!
//! Implements full-wave electromagnetic scattering and radiation using the
//! Method of Moments (MoM) with RWG basis functions and CFIE formulation.
//!
//! # Architecture
//! ```text
//! RemMesh  →  SurfaceMesh  →  RwgBases
//!                  ↓               ↓
//!             quadrature      assemble Z (dense matrix)
//!                  ↓               ↓
//!             singular.rs    LU solve → currents
//!                                  ↓
//!                            postprocess → RCS CSV
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
pub mod port;
pub mod sparams;
pub mod sibc;
pub mod fft_accel;

// Public re-exports for cross-crate solver integration
pub use assemble::{gmres_solve, gmres_solve_generic, gmres_solve_op, aca_gmres_solve, gmres_generic_with_aca};

use rem_config::{PalaceConfig, MomSolverConfig};
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
        }).collect();
        log::info!(
            "MoM: using layered Green function ({} layer(s), bottom_pec={})",
            layers.len(), sub.bottom_pec
        );
        Box::new(LayeredGreen::new(layers, k0))
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
    /// RCS in dBsm = 10·log10(rcs_m2); -300 if rcs_m2 ≈ 0
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
            let z_mat = match mom_cfg.basis.as_str() {
                "Pulse" | "pulse" => {
                    if mom_cfg.wall_conductivity > 0.0 {
                        log::warn!(
                            "MoM SIBC: WallConductivity is currently supported only for RWG basis; ignoring for Pulse basis"
                        );
                    }
                    assemble::assemble_efie_pulse(&surf, freq, &quad, mom_cfg.singular_tol)
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
                        log::warn!("MoM FFT: mesh is not planar — falling back to GMRES with dense matrix");
                        assemble::gmres_solve(&z_mat, &rhs)?
                    }
                }
                "FMM" => {
                    return Err(rem_core::RemError::Config(
                        "FastSolver \"FMM\" is not yet implemented; use \"Direct\", \"GMRES\", \"ACA\", or \"FFT\"".to_string()
                    ));
                }
                _ => assemble::lu_solve(&z_mat, &rhs)?,
            }
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

    // Build lumped ports
    let lumped_ports: Vec<MomLumpedPort> = mom_cfg.ports.iter().map(|p| {
        let z0 = p.impedance.unwrap_or(mom_cfg.ref_impedance);
        MomLumpedPort::from_surface(surf, &bases, &p.attributes, p.index, &p.direction, z0)
    }).collect::<RemResult<_>>()?;

    // Frequency sweep
    let mut freq = mom_cfg.freq_min;
    let freq_max  = mom_cfg.freq_max;
    let freq_step = mom_cfg.freq_step;
    let mut all_matrices: Vec<sparams::SMatrix> = Vec::new();

    while freq <= freq_max + 1e-3 * freq_step {
        log::info!("MoM S-param solve at f = {:.3e} Hz", freq);

        let green = build_green(mom_cfg, freq);
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

        let sm = sparams::compute_s_matrix(surf, &bases, &lumped_ports, &z_mat, freq)?;
        log::info!("  S-matrix computed: {}×{}", sm.n_ports, sm.n_ports);
        all_matrices.push(sm);

        freq += freq_step;
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
    }

    log::info!("MoM S-param sweep complete. {} frequency points.", all_matrices.len());
    Ok(MomResult { rcs: vec![] })  // no RCS in S-param mode
}
