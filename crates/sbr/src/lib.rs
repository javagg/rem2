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
    use num_complex::Complex64;
    use rem_core::ETA0;

    let mut currents = zero_currents(surf);

    // Pre-compute self-intersection offset: 1e-5 × √(min face area)
    // Moved outside the closure so it isn't recomputed per ray per bounce.
    let eps_offset = 1e-5 * surf.faces.iter()
        .map(|f| f.area.sqrt())
        .fold(f64::INFINITY, f64::min)
        .max(1e-10);

    // Each ray returns a list of (face_index, J contribution)
    let trace_one = |mut r: ray::Ray| -> Vec<(usize, [Complex64; 3])> {
        let mut contributions = Vec::new();

        loop {
            // Termination checks
            if !r.weight.is_finite() || r.weight < cfg.weight_thresh {
                break;
            }
            if r.bounce >= cfg.max_bounces {
                break;
            }

            // ── Nearest intersection ──────────────────────────────────────
            let hit = match bvh.intersect(&r.origin, &r.dir) {
                Some(h) => h,
                None    => break,
            };

            // Ensure face normal faces the incoming ray
            let face_normal = {
                let fn_ = surf.faces[hit.face_idx].normal;
                if ray::dot3(&r.dir, &fn_) < 0.0 { fn_ } else {
                    [-fn_[0], -fn_[1], -fn_[2]]
                }
            };

            // ── PO current: J = 2 n̂ × H at hit point ─────────────────────
            // For bounce 0 (incident ray): H from the incident plane wave.
            // For bounce n>0: H is stored in the ray (from previous reflection).
            let h_at_hit: [Complex64; 3] = if r.bounce == 0 {
                let (_, h) = incident_fields(wave, k, &hit.point);
                h
            } else {
                // Propagate ray H by the phase delay from r.origin to hit.point
                let dist = hit.t; // parametric t = physical distance (unit dir)
                let phase_delay = Complex64::new(0.0, k * dist).exp();
                [r.h_field[0] * phase_delay,
                 r.h_field[1] * phase_delay,
                 r.h_field[2] * phase_delay]
            };

            let j_po = po_current_pec(&h_at_hit, &face_normal, &r.dir);
            contributions.push((hit.face_idx, j_po));

            // ── Reflection ────────────────────────────────────────────────
            // Propagate E to hit point
            let dist = hit.t;
            let phase_delay = Complex64::new(0.0, k * dist).exp();
            let e_at_hit = [r.e_field[0] * phase_delay,
                            r.e_field[1] * phase_delay,
                            r.e_field[2] * phase_delay];
            let h_at_hit_e = [h_at_hit[0], h_at_hit[1], h_at_hit[2]];

            let iface = Interface::pec();
            let (e_refl, new_dir) = reflect_field(
                &e_at_hit, &h_at_hit_e, &r.dir, &face_normal, &iface,
            );

            // Derive reflected H from reflected E: H = (k̂_refl × E_refl) / η₀
            let h_refl = [
                (new_dir[1] * e_refl[2] - new_dir[2] * e_refl[1]) / ETA0,
                (new_dir[2] * e_refl[0] - new_dir[0] * e_refl[2]) / ETA0,
                (new_dir[0] * e_refl[1] - new_dir[1] * e_refl[0]) / ETA0,
            ];

            // Offset origin to avoid re-hitting the same face
            let p_offset = ray::add3(&hit.point, &ray::scale3(&face_normal, eps_offset));

            r.origin  = p_offset;
            r.dir     = new_dir;
            r.e_field = e_refl;
            r.h_field = h_refl;
            r.bounce += 1;
            // PEC: |Γ| = 1, weight unchanged; for dielectrics multiply by |Γ|
        }

        contributions
    };

    // Parallel dispatch (non-WASM) / serial (WASM)
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
