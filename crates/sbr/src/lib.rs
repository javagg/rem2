//! rem-sbr — Shooting and Bouncing Rays + Physical Optics solver.
//!
//! Entry point: `run(config)` dispatched from `rem-cli` for
//! `Problem.Type = "SBR"`.
//!
//! # Algorithm
//! ```text
//! 1. Load mesh → extract PEC surface → build BVH
//! 2. For each frequency:
//!    a. Launch aperture rays (uniform grid ⊥ k̂_inc)
//!    b. For each ray, trace bounces:
//!         intersect BVH → compute J_PO = 2 n̂ × H_inc → update currents
//!         reflect ray (Fresnel / PEC mirror) → next bounce
//!    c. Far-field PO integral → RCS CSV + surface VTK
//! ```

pub mod ray;
pub mod bvh;
pub mod excitation;
pub mod fresnel;
pub mod po_integral;
pub mod output;

use std::f64::consts::PI;
use std::sync::Arc;

use rem_config::{PalaceConfig, SbrSolverConfig};
use rem_core::{C0, RemError, RemResult};
use rem_mom::surface_mesh::SurfaceMesh;
use rem_parallel::NoComm;

use bvh::Bvh;
use excitation::{PlaneWave, incident_fields, launch_aperture_rays};
use fresnel::{Interface, reflect_field, po_current_pec};
use po_integral::{zero_currents, CurrentMap};
use output::{write_rcs, write_surface_vtk};

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

/// Run SBR+ solver. Called from CLI for `Problem.Type = "SBR"`.
pub fn run(config: &PalaceConfig) -> RemResult<()> {
    let sbr_cfg = config.solver.sbr.as_ref()
        .ok_or_else(|| RemError::Config(
            "Problem.Type = \"SBR\" requires a Solver.SBR section".to_string()
        ))?;

    let mesh = rem_mesh::load_mesh(config, &NoComm)?;
    run_with_mesh(config, sbr_cfg, &mesh)
}

/// Inner function (also callable from tests / WASM).
pub fn run_with_mesh(
    config: &PalaceConfig,
    sbr_cfg: &SbrSolverConfig,
    mesh: &rem_mesh::RemMesh,
) -> RemResult<()> {
    // ── 1. Surface mesh + BVH ─────────────────────────────────────────────
    let pec_attrs: Vec<u32> = config.boundaries.pec
        .as_ref()
        .map(|p| p.attributes.clone())
        .unwrap_or_default();

    if pec_attrs.is_empty() {
        return Err(RemError::Config(
            "SBR solver requires at least one PEC boundary (Boundaries.PEC.Attributes)".to_string()
        ));
    }

    let surf = Arc::new(SurfaceMesh::extract(mesh, &pec_attrs)?);
    log::info!("SBR surface mesh: {} faces", surf.faces.len());

    let bvh = Bvh::build(Arc::clone(&surf));
    log::info!("BVH built with {} nodes", surf.faces.len());

    // ── 2. Output directory ───────────────────────────────────────────────
    let output_dir = std::path::Path::new(config.problem.output_dir());
    std::fs::create_dir_all(output_dir.join("postpro"))?;

    // ── 3. RCS observation angles ─────────────────────────────────────────
    let (theta_deg, phi_deg) = if let Some(rcs) = &config.postprocessing.rcs {
        (rcs.theta_deg.clone(), rcs.phi_deg.clone())
    } else {
        let theta: Vec<f64> = (0..=180).step_by(5).map(|i| i as f64).collect();
        (theta, vec![0.0])
    };

    // ── 4. Frequency sweep ────────────────────────────────────────────────
    let freq_step = if sbr_cfg.freq_step > 0.0 { sbr_cfg.freq_step } else { sbr_cfg.freq_max - sbr_cfg.freq_min + 1.0 };
    let mut freq = sbr_cfg.freq_min;
    while freq <= sbr_cfg.freq_max + 1e-3 * freq_step {
        log::info!("SBR+ solve at f = {:.3e} Hz", freq);

        let k = 2.0 * PI * freq / C0;

        // Incident plane wave
        let wave = PlaneWave {
            theta_inc: sbr_cfg.theta_inc_deg.to_radians(),
            phi_inc:   sbr_cfg.phi_inc_deg.to_radians(),
            pol:       sbr_cfg.polarization.clone(),
        };

        // Launch aperture rays
        let init_rays = launch_aperture_rays(&wave, &surf, k, sbr_cfg.ray_density, freq);
        log::info!("  {} aperture rays launched", init_rays.len());

        // Trace all rays and accumulate PO currents
        let currents = trace_all_rays(
            init_rays, &bvh, &surf, &wave, k, sbr_cfg,
        );

        // Far-field RCS → CSV
        write_rcs(output_dir, freq, &currents, &surf, k, &theta_deg, &phi_deg)?;

        // Surface current VTK
        let vtk_path = output_dir
            .join("postpro")
            .join(format!("sbr_{:.3e}Hz.vtk", freq));
        write_surface_vtk(&vtk_path, &currents, &surf)?;

        freq += freq_step;
    }

    log::info!("SBR+ complete. Results in {}", output_dir.display());
    Ok(())
}

// ---------------------------------------------------------------------------
// Ray tracing kernel
// ---------------------------------------------------------------------------

fn trace_all_rays(
    init_rays: Vec<ray::Ray>,
    bvh: &Bvh,
    surf: &SurfaceMesh,
    wave: &PlaneWave,
    k: f64,
    cfg: &SbrSolverConfig,
) -> CurrentMap {
    let mut currents = zero_currents(surf);

    let trace_one = |mut r: ray::Ray| -> Vec<(usize, [num_complex::Complex64; 3])> {
        let mut contributions = Vec::new();

        loop {
            if !r.weight.is_finite() || r.weight < cfg.weight_thresh {
                break;
            }
            if r.bounce >= cfg.max_bounces {
                break;
            }

            // Intersect BVH
            let hit = match bvh.intersect(&r.origin, &r.dir) {
                Some(h) => h,
                None    => break,
            };

            // Offset hit point slightly along normal to avoid self-intersection
            let eps_offset = 1e-9 * surf.faces.iter()
                .map(|f| f.area.sqrt())
                .fold(f64::INFINITY, f64::min)
                .max(1e-12);
            let p_offset = ray::add3(&hit.point, &ray::scale3(&hit.normal, eps_offset));

            // --- PO current accumulation (PEC) ---
            // Compute incident H at hit point
            let (_, h_hit) = incident_fields(wave, k, &hit.point);
            let j_po = po_current_pec(&h_hit, &hit.normal, &r.dir);
            contributions.push((hit.face_idx, j_po));

            // --- Reflection ---
            let iface = Interface::pec();
            let (e_refl, new_dir) = reflect_field(&r.e_field, &r.h_field, &r.dir, &hit.normal, &iface);

            // Compute reflected H from reflected E via H = (k̂ × E)/η₀
            use rem_core::ETA0;
            let mut h_refl = [num_complex::Complex64::ZERO; 3];
            let nd = new_dir;
            h_refl[0] = (nd[1]*e_refl[2] - nd[2]*e_refl[1]) / ETA0;
            h_refl[1] = (nd[2]*e_refl[0] - nd[0]*e_refl[2]) / ETA0;
            h_refl[2] = (nd[0]*e_refl[1] - nd[1]*e_refl[0]) / ETA0;

            // PEC: total reflection, weight reduced by 1 (stays at 1 ideally)
            // For coated/dielectric targets, multiply by |Γ|
            r.origin  = p_offset;
            r.dir     = new_dir;
            r.e_field = e_refl;
            r.h_field = h_refl;
            r.bounce += 1;
            // PEC reflection coefficient magnitude = 1; weight unchanged
        }

        contributions
    };

    // Parallel dispatch on non-WASM targets
    #[cfg(not(target_arch = "wasm32"))]
    {
        use rayon::prelude::*;
        let all: Vec<_> = init_rays.into_par_iter().map(trace_one).collect();
        for ray_hits in all {
            for (fi, j) in ray_hits {
                for i in 0..3 { currents[fi].j[i] += j[i]; }
            }
        }
    }

    #[cfg(target_arch = "wasm32")]
    {
        for r in init_rays {
            for (fi, j) in trace_one(r) {
                for i in 0..3 { currents[fi].j[i] += j[i]; }
            }
        }
    }

    currents
}
