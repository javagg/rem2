//! rem-mom — Method of Moments solver for rem2
//!
//! Implements full-wave electromagnetic scattering and radiation using the
//! Method of Moments (MoM) with RWG basis functions and CFIE formulation.
//!
//! # Architecture
//! ```text
//! RemMesh  →  SurfaceMesh  →  RwgBases
//!                  ↓               ↓
//!             quadrature      assemble Z (faer dense)
//!                  ↓               ↓
//!             singular.rs    faer LU solve → currents
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

use rem_config::{PalaceConfig, MomSolverConfig};
use rem_core::RemResult;
use rem_mesh::RemMesh;
use rem_parallel::NoComm;

/// Entry point called from the CLI for `Problem.Type = "MoM"`.
pub fn run(config: &PalaceConfig) -> RemResult<()> {
    let mom_cfg = config.solver.mom.as_ref()
        .ok_or_else(|| rem_core::RemError::Config(
            "Problem.Type = \"MoM\" requires a Solver.MoM section".to_string()
        ))?;

    let mesh = rem_mesh::load_mesh(config, &NoComm)?;
    run_with_mesh(config, mom_cfg, &mesh)
}

/// Run MoM solve on an already-loaded mesh (also used from WASM / tests).
pub fn run_with_mesh(
    config: &PalaceConfig,
    mom_cfg: &MomSolverConfig,
    mesh: &RemMesh,
) -> RemResult<()> {
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
    while freq <= freq_max + 1e-3 * freq_step {
        log::info!("MoM solve at f = {:.3e} Hz", freq);

        let k = 2.0 * PI * freq / rem_core::C0;

        // Assemble impedance matrix Z
        let z_mat = match mom_cfg.basis.as_str() {
            "Pulse" | "pulse" => {
                assemble::assemble_efie_pulse(&surf, freq, &quad, mom_cfg.singular_tol)
            }
            _ => {
                // RWG default
                let bases = basis::rwg::generate_rwg_bases(&surf);
                assemble::assemble_cfie_rwg(&surf, &bases, freq, mom_cfg.alpha, &quad, mom_cfg.singular_tol)
            }
        }?;

        // Incident plane wave excitation (+z direction, x-polarized)
        let rhs = excitation::plane_wave_rhs(&surf, k, &mom_cfg.basis);

        // Solve Z·I = V
        let currents = assemble::lu_solve(&z_mat, &rhs)?;

        // Post-process: RCS
        postprocess::write_rcs(
            output_dir,
            freq,
            &currents,
            &surf,
            k,
            &theta_deg,
            &phi_deg,
        )?;

        freq += freq_step;
    }

    log::info!("MoM solve complete. Results in {}", output_dir.display());
    Ok(())
}
