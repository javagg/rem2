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
pub mod singular;
pub mod assemble;
pub mod excitation;
pub mod postprocess;
pub mod basis;
pub mod mie;
pub mod aca;
pub mod pmchwt;

use rem_config::{PalaceConfig, MomSolverConfig};
use rem_core::RemResult;
use rem_mesh::RemMesh;
use rem_parallel::NoComm;

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
    log::info!("MoM surface mesh: {} faces, {} interior edges (RWG bases)",
        surf.faces.len(), surf.edges.len());

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

        // Incident plane wave
        let wave = excitation::PlaneWave {
            theta_inc: mom_cfg.theta_inc_deg.to_radians(),
            phi_inc:   mom_cfg.phi_inc_deg.to_radians(),
            pol:       mom_cfg.polarization.clone(),
        };

        // PMCHWT path: dielectric target (J + M unknowns, 2N×2N system)
        let currents = if mom_cfg.equation.to_uppercase() == "PMCHWT" {
            let (eps_r, mu_r) = config.domains.materials.first()
                .map(|m| (m.permittivity, m.permeability))
                .unwrap_or((2.0, 1.0));
            let mat = pmchwt::DielectricMaterial::new(eps_r, mu_r);
            log::info!("MoM PMCHWT: ε_r={eps_r:.2}, μ_r={mu_r:.2}, f={freq:.3e} Hz");
            let (j_coeffs, _m_coeffs) = pmchwt::solve_pmchwt(
                &surf, mat, freq, &wave, &quad, &mom_cfg.fast_solver,
            )?;
            j_coeffs
        } else {
            // Assemble impedance matrix Z  (PEC EFIE/MFIE/CFIE path)
            let z_mat = match mom_cfg.basis.as_str() {
                "Pulse" | "pulse" => {
                    assemble::assemble_efie_pulse(&surf, freq, &quad, mom_cfg.singular_tol)
                }
                _ => {
                    let bases = basis::rwg::generate_rwg_bases(&surf);
                    assemble::assemble_cfie_rwg(&surf, &bases, freq, mom_cfg.alpha, &quad, mom_cfg.singular_tol)
                }
            }?;
            let rhs = excitation::plane_wave_rhs_general(&surf, k, &wave, &mom_cfg.basis);
            match mom_cfg.fast_solver.to_uppercase().as_str() {
                "GMRES" => assemble::gmres_solve(&z_mat, &rhs)?,
                "ACA" => {
                    log::info!("MoM: using ACA+GMRES (tol_aca=1e-4, tol_gmres=1e-8)");
                    assemble::aca_gmres_solve(&z_mat, &rhs, 1e-4, 1e-8)?
                }
                "FMM" => {
                    return Err(rem_core::RemError::Config(
                        "FastSolver \"FMM\" is not yet implemented; use \"Direct\", \"GMRES\", or \"ACA\"".to_string()
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
        }

        freq += freq_step;
    }

    #[cfg(not(target_arch = "wasm32"))]
    log::info!("MoM solve complete. Results in {}", output_dir.display());
    #[cfg(target_arch = "wasm32")]
    log::info!("MoM solve complete.");
    Ok(MomResult { rcs: all_rcs })
}
