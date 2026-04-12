//! rem-sbr — Shooting and Bouncing Rays + Physical Optics solver.
//!
//! Entry point: `run(config)` dispatched from `rem-cli` for
//! `Problem.Type = "SBR"`.
//!
//! # Algorithm
//! ```text
//! 1. Load mesh → extract PEC surface → build BVH
//! 2. For each frequency:
//!    a. First-bounce PO (per-face): J = 2 n̂ × H_inc, shadow-tested via BVH
//!    b. Multi-bounce (per-ray): trace reflected rays (bounce ≥ 1),
//!       J += (A_ray / A_face) × 2 n̂ × H_ray  (flux-to-density conversion)
//!    c. Far-field PO integral → RCS CSV + surface VTK
//! ```
//!
//! ## Why two stages?
//!
//! First-bounce PO should be computed once per illuminated face (not once per
//! ray that happens to hit the face), because the ray density just controls
//! angular sampling resolution.  Accumulating first-bounce J from rays leads to
//! J ∝ ray_density, which then over-counts in the far-field integral.
//!
//! For bounces ≥ 1 the ray carries a flux per unit cross-sectional area
//! (A_ray = spacing²); to convert it to a surface current *density* [A/m] we
//! scale by A_ray / A_face.

pub mod ray;
pub mod bvh;
pub mod excitation;
pub mod fresnel;
pub mod po_integral;
pub mod output;
pub mod ptd;

use std::f64::consts::PI;
use std::sync::Arc;

use num_complex::Complex64;
use rem_config::{PalaceConfig, SbrSolverConfig};
use rem_core::{C0, ETA0, RemError, RemResult};
use rem_mom::surface_mesh::SurfaceMesh;
use rem_parallel::NoComm;

use bvh::Bvh;
use excitation::{PlaneWave, incident_fields, launch_aperture_rays};
use fresnel::{Interface, reflect_field, po_current_pec};
use po_integral::{zero_currents, CurrentMap, rcs_pattern_with_ptd};
use ptd::{extract_boundary_edges};
use output::{write_rcs_with_ptd, write_surface_vtk};

// ---------------------------------------------------------------------------
// Result types
// ---------------------------------------------------------------------------

/// One RCS observation point.
#[derive(Debug, Clone)]
pub struct RcsPoint {
    pub theta_deg: f64,
    pub phi_deg:   f64,
    pub rcs_m2:    f64,
    pub rcs_dbsm:  f64,
}

/// Result returned by `run_with_mesh`.
#[derive(Debug, Clone)]
pub struct SbrResult {
    pub rcs: Vec<(f64, Vec<RcsPoint>)>,   // (freq_hz, points)
}

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
    run_with_mesh(config, sbr_cfg, &mesh).map(|_| ())
}

/// Inner function (also callable from tests / WASM).
pub fn run_with_mesh(
    config: &PalaceConfig,
    sbr_cfg: &SbrSolverConfig,
    mesh: &rem_mesh::RemMesh,
) -> RemResult<SbrResult> {
    log::info!("\n=== Shooting and Bouncing Rays (SBR) solver ===\n");

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
    log::info!("Surface mesh (PEC scatterer):");
    log::info!("  {} triangular faces", surf.faces.len());

    // Pre-compute boundary edges for PTD correction
    let ptd_edges = extract_boundary_edges(&surf);
    log::info!("PTD boundary diffraction:");
    log::info!("  {} boundary edges", ptd_edges.len());

    let bvh = Bvh::build(Arc::clone(&surf));
    log::info!("BVH acceleration structure built");
    log::info!("");

    // ── 2. Output directory ───────────────────────────────────────────────
    let output_dir = std::path::Path::new(config.problem.output_dir());
    #[cfg(not(target_arch = "wasm32"))]
    std::fs::create_dir_all(output_dir.join("postpro"))?;

    // ── 3. RCS observation angles ─────────────────────────────────────────
    let (theta_deg, phi_deg) = if let Some(rcs) = &config.postprocessing.rcs {
        (rcs.theta_deg.clone(), rcs.phi_deg.clone())
    } else {
        let theta: Vec<f64> = (0..=180).step_by(5).map(|i| i as f64).collect();
        (theta, vec![0.0])
    };

    // ── 4. Frequency sweep ────────────────────────────────────────────────
    let freq_step = if sbr_cfg.freq_step > 0.0 {
        sbr_cfg.freq_step
    } else {
        sbr_cfg.freq_max - sbr_cfg.freq_min + 1.0
    };
    let mut freq = sbr_cfg.freq_min;
    let mut all_rcs: Vec<(f64, Vec<RcsPoint>)> = Vec::new();
    while freq <= sbr_cfg.freq_max + 1e-3 * freq_step {
        log::info!("SBR+ solve at f = {:.3e} Hz", freq);

        let k = 2.0 * PI * freq / C0;

        // Incident plane wave
        let wave = PlaneWave {
            theta_inc: sbr_cfg.theta_inc_deg.to_radians(),
            phi_inc:   sbr_cfg.phi_inc_deg.to_radians(),
            pol:       sbr_cfg.polarization.clone(),
        };

        // ── Stage 1: first-bounce PO (per-face, ray-density independent) ──
        let mut currents = first_bounce_po(&surf, &bvh, &wave, k);

        // ── Stage 2: multi-bounce rays (bounce ≥ 1) ───────────────────────
        if sbr_cfg.max_bounces > 1 {
            // Ray area = spacing² (aperture grid cell area [m²])
            let spacing = 1.0 / sbr_cfg.ray_density.sqrt();
            let a_ray = spacing * spacing;

            let init_rays = launch_aperture_rays(&wave, &surf, k, sbr_cfg.ray_density, freq);
            log::info!("  {} aperture rays launched for multi-bounce", init_rays.len());

            let mb = multibounce_rays(init_rays, &bvh, &surf, &wave, k, sbr_cfg, a_ray);
            for (fi, fc) in mb.iter().enumerate() {
                for i in 0..3 {
                    currents[fi].j[i] += fc.j[i];
                }
            }
        }

        // ── Compute RCS pattern (always, not just non-WASM) ──────────────────
        let wave_ptd = wave.clone();
        let k_ptd = k;
        let e_fn = move |r: &[f64; 3]| -> [Complex64; 3] {
            let (e, _h) = incident_fields(&wave_ptd, k_ptd, r);
            e
        };
        let pattern = rcs_pattern_with_ptd(
            &currents, &surf, k, &theta_deg, &phi_deg,
            &wave, &ptd_edges, &e_fn,
        );
        let mut pts = Vec::new();
        for (ti, &th) in theta_deg.iter().enumerate() {
            for (pi, &ph) in phi_deg.iter().enumerate() {
                let rcs_m2 = pattern[ti][pi];
                let rcs_dbsm = if rcs_m2 > 1e-40 { 10.0 * rcs_m2.log10() } else { -300.0 };
                pts.push(RcsPoint { theta_deg: th, phi_deg: ph, rcs_m2, rcs_dbsm });
            }
        }
        all_rcs.push((freq, pts));

        #[cfg(not(target_arch = "wasm32"))]
        {
            // Build incident field closure for PTD (captures k, wave)
            let wave_ptd2 = wave.clone();
            let k_ptd2 = k;
            let e_fn2 = move |r: &[f64; 3]| -> [Complex64; 3] {
                let (e, _h) = incident_fields(&wave_ptd2, k_ptd2, r);
                e
            };

            write_rcs_with_ptd(
                output_dir, freq, &currents, &surf, k,
                &theta_deg, &phi_deg, &wave, &ptd_edges, &e_fn2,
            )?;

            let vtk_path = output_dir
                .join("postpro")
                .join(format!("sbr_{:.3e}Hz.vtk", freq));
            write_surface_vtk(&vtk_path, &currents, &surf)?;
        }

        freq += freq_step;
    }

    #[cfg(not(target_arch = "wasm32"))]
    log::info!("SBR+ complete. Results in {}", output_dir.display());
    #[cfg(target_arch = "wasm32")]
    log::info!("SBR+ complete.");
    Ok(SbrResult { rcs: all_rcs })
}

// ---------------------------------------------------------------------------
// Stage 1 – First-bounce PO  (per-face, shadow-tested)
// ---------------------------------------------------------------------------

/// Compute PO surface currents for the direct (first-bounce) illumination.
///
/// For each face in `surf`:
///   1. Geometric visibility: dot(n̂_face, −k̂_inc) > 0 (face points toward source)
///   2. Shadow test: cast a ray from the face centroid toward the source;
///      if it hits another face the centroid is in shadow → skip.
///   3. If illuminated: `J = 2 n̂ × H_inc(centroid)`
///
/// This is *independent of ray density* — every illuminated face gets exactly
/// one contribution from the incident wave.
fn first_bounce_po(
    surf: &SurfaceMesh,
    bvh: &Bvh,
    wave: &PlaneWave,
    k: f64,
) -> CurrentMap {
    let kh = wave.k_hat(); // unit incident direction
    let neg_kh = [-kh[0], -kh[1], -kh[2]]; // toward source

    // Self-intersection offset (move centroid slightly toward source)
    let eps = 1e-6 * surf.faces.iter()
        .map(|f| f.area.sqrt())
        .fold(f64::INFINITY, f64::min)
        .max(1e-10);

    let mut currents = zero_currents(surf);

    for (fi, face) in surf.faces.iter().enumerate() {
        // 1. Geometric test: face must point toward incoming wave
        let cos_inc = ray::dot3(&face.normal, &neg_kh);
        if cos_inc <= 0.0 {
            continue; // back-face or edge-on
        }

        // 2. Shadow test: is centroid visible from source direction?
        //
        // IMPORTANT: the centroid of a curved triangular face lies slightly
        // INSIDE the surface (vertices are on the sphere, centroid is not).
        // We must offset along the outward face normal (not along neg_kh) to
        // place the shadow-ray origin outside the surface; otherwise the shadow
        // ray starts inside the closed mesh and immediately hits the back face.
        let offset = (0.01 * face.area.sqrt()).max(eps);
        let shadow_origin = ray::add3(&face.centroid, &ray::scale3(&face.normal, offset));
        if bvh.any_hit(&shadow_origin, &neg_kh, f64::INFINITY) {
            continue; // shadowed by another part of the surface
        }

        // 3. PO current: J = 2 n̂ × H_inc
        let (_e_inc, h_inc) = incident_fields(wave, k, &face.centroid);
        let j = po_current_pec(&h_inc, &face.normal, &kh);
        currents[fi].j = j;
    }

    currents
}

// ---------------------------------------------------------------------------
// Stage 2 – Multi-bounce rays  (bounce ≥ 1)
// ---------------------------------------------------------------------------

/// Trace reflected rays and accumulate PO contributions from bounces ≥ 1.
///
/// Each ray carries a flux tube of cross-sectional area `a_ray` [m²].
/// When it hits a face of area `A_face`, the surface current density is:
///
/// ```text
/// ΔJ = (a_ray / A_face) × 2 n̂ × H_ray
/// ```
///
/// This ensures that the total current on a face scales with the illuminated
/// power per unit area, independent of how many rays happen to hit it.
fn multibounce_rays(
    init_rays: Vec<ray::Ray>,
    bvh: &Bvh,
    surf: &SurfaceMesh,
    _wave: &PlaneWave,
    k: f64,
    cfg: &SbrSolverConfig,
    a_ray: f64,
) -> CurrentMap {
    // Self-intersection offset
    let eps_offset = 1e-6 * surf.faces.iter()
        .map(|f| f.area.sqrt())
        .fold(f64::INFINITY, f64::min)
        .max(1e-10);

    let trace_one = |mut r: ray::Ray| -> Vec<(usize, [Complex64; 3])> {
        let mut contributions = Vec::new();

        // First intersection = first bounce (bounce 0 is handled by first_bounce_po)
        // We trace the ray and skip accumulating J at bounce 0, only bounces ≥ 1.
        // However, we still need to *reflect* at bounce 0 to get the right direction.
        // Strategy: trace normally, but only accumulate for bounce > 0.

        loop {
            if !r.weight.is_finite() || r.weight < cfg.weight_thresh {
                break;
            }
            if r.bounce >= cfg.max_bounces {
                break;
            }

            let hit = match bvh.intersect(&r.origin, &r.dir) {
                Some(h) => h,
                None    => break,
            };

            let face = &surf.faces[hit.face_idx];

            // Ensure face normal faces incoming ray
            let face_normal: [f64; 3] = {
                let fn_ = face.normal;
                if ray::dot3(&r.dir, &fn_) < 0.0 { fn_ } else {
                    [-fn_[0], -fn_[1], -fn_[2]]
                }
            };

            // Propagate E and H to hit point via phase delay
            let dist = hit.t;
            let phase_delay = Complex64::new(0.0, k * dist).exp();
            let e_at_hit = [
                r.e_field[0] * phase_delay,
                r.e_field[1] * phase_delay,
                r.e_field[2] * phase_delay,
            ];
            let h_at_hit = [
                r.h_field[0] * phase_delay,
                r.h_field[1] * phase_delay,
                r.h_field[2] * phase_delay,
            ];

            // Accumulate PO current only for bounce ≥ 1
            if r.bounce >= 1 {
                let j_po = po_current_pec(&h_at_hit, &face_normal, &r.dir);
                // Scale by a_ray / A_face to convert flux to surface current density
                let scale = a_ray / face.area.max(1e-30);
                let j_scaled = [
                    j_po[0] * scale,
                    j_po[1] * scale,
                    j_po[2] * scale,
                ];
                contributions.push((hit.face_idx, j_scaled));
            }

            // Reflect
            let iface = Interface::pec();
            let (e_refl, new_dir) = reflect_field(
                &e_at_hit, &h_at_hit, &r.dir, &face_normal, &iface,
            );

            // Derive reflected H from reflected E: H = (k̂_refl × E_refl) / η₀
            let h_refl = [
                (new_dir[1] * e_refl[2] - new_dir[2] * e_refl[1]) / ETA0,
                (new_dir[2] * e_refl[0] - new_dir[0] * e_refl[2]) / ETA0,
                (new_dir[0] * e_refl[1] - new_dir[1] * e_refl[0]) / ETA0,
            ];

            let p_offset = ray::add3(&hit.point, &ray::scale3(&face_normal, eps_offset));

            r.origin  = p_offset;
            r.dir     = new_dir;
            r.e_field = e_refl;
            r.h_field = h_refl;
            r.bounce += 1;
        }

        contributions
    };

    let mut currents = zero_currents(surf);

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
