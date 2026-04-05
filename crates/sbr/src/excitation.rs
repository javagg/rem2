//! Plane-wave excitation and aperture ray launch for SBR+.
//!
//! Reuses the `PlaneWave` definition from `rem_mom::excitation` and adds the
//! aperture-based ray grid needed by the SBR+ algorithm.

use num_complex::Complex64;
use rem_core::{C0, ETA0};
use rem_mom::surface_mesh::SurfaceMesh;
use crate::ray::{Ray, dot3, cross3, sub3, add3, scale3, normalize3};

// ---------------------------------------------------------------------------
// Re-export PlaneWave from MoM (same definition, same convention)
// ---------------------------------------------------------------------------

pub use rem_mom::excitation::PlaneWave;

// ---------------------------------------------------------------------------
// Incident field evaluation
// ---------------------------------------------------------------------------

/// Compute E_inc and H_inc at position `r` for a plane wave at wave-number `k`.
///
/// ```text
/// E_inc(r) = ê * exp(-jk k̂·r)
/// H_inc(r) = (k̂ × ê) / η₀ * exp(-jk k̂·r)
/// ```
pub fn incident_fields(
    wave: &PlaneWave,
    k: f64,
    r: &[f64; 3],
) -> ([Complex64; 3], [Complex64; 3]) {
    let kh = wave.k_hat();
    let eh = wave.e_hat();
    let hh = cross3(&kh, &eh); // H polarization unit vector

    let phase = k * dot3(&kh, r);
    let phasor = Complex64::new(0.0, -phase).exp(); // exp(-jk k̂·r)

    let e = [
        Complex64::new(eh[0], 0.0) * phasor,
        Complex64::new(eh[1], 0.0) * phasor,
        Complex64::new(eh[2], 0.0) * phasor,
    ];
    let h = [
        Complex64::new(hh[0] / ETA0, 0.0) * phasor,
        Complex64::new(hh[1] / ETA0, 0.0) * phasor,
        Complex64::new(hh[2] / ETA0, 0.0) * phasor,
    ];
    (e, h)
}

// ---------------------------------------------------------------------------
// Aperture ray launch
// ---------------------------------------------------------------------------

/// Launch a uniform grid of rays from a virtual aperture plane perpendicular
/// to the incident direction.
///
/// The aperture is sized to enclose the target bounding box (× 1.2 margin)
/// projected onto the plane normal to `k̂`.  Rays are placed on a square
/// grid with spacing `Δ = 1 / sqrt(ray_density)` [m].
pub fn launch_aperture_rays(
    wave: &PlaneWave,
    surf: &SurfaceMesh,
    k: f64,
    ray_density: f64,
    freq: f64,
) -> Vec<Ray> {
    let kh = wave.k_hat();

    // ── 1. Target bounding box ────────────────────────────────────────────
    let (mut bbox_min, mut bbox_max) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    for &[x, y, z] in &surf.nodes {
        bbox_min[0] = bbox_min[0].min(x); bbox_max[0] = bbox_max[0].max(x);
        bbox_min[1] = bbox_min[1].min(y); bbox_max[1] = bbox_max[1].max(y);
        bbox_min[2] = bbox_min[2].min(z); bbox_max[2] = bbox_max[2].max(z);
    }
    let center = scale3(&add3(&bbox_min, &bbox_max), 0.5);
    let extent = sub3(&bbox_max, &bbox_min);
    let radius = (extent[0]*extent[0] + extent[1]*extent[1] + extent[2]*extent[2]).sqrt() * 0.6;

    // ── 2. Build local coordinate system in aperture plane ───────────────
    // Pick an "up" vector not parallel to k̂
    let up_candidate = if kh[0].abs() < 0.9 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] };
    let u_hat = normalize3(cross3(&kh, &up_candidate));
    let v_hat = cross3(&kh, &u_hat);

    // ── 3. Grid spacing from ray density ─────────────────────────────────
    let lambda = C0 / freq;
    let spacing = 1.0 / ray_density.sqrt(); // [m]
    let _ = lambda; // available for future wavelength-adaptive spacing

    let n_side = (2.0 * radius / spacing).ceil() as i32 + 1;
    let aperture_origin = sub3(&center, &scale3(&kh, radius * 2.0)); // far upstream

    // ── 4. Emit one ray per grid cell ─────────────────────────────────────
    let mut rays = Vec::with_capacity((n_side * n_side) as usize);

    for iu in (-n_side / 2)..=(n_side / 2) {
        for iv in (-n_side / 2)..=(n_side / 2) {
            let pu = iu as f64 * spacing;
            let pv = iv as f64 * spacing;

            // Discard grid points outside the circular aperture
            if pu*pu + pv*pv > radius*radius { continue; }

            let origin = add3(
                &aperture_origin,
                &add3(&scale3(&u_hat, pu), &scale3(&v_hat, pv)),
            );

            let (e_field, h_field) = incident_fields(wave, k, &origin);
            rays.push(Ray::new(origin, kh, e_field, h_field));
        }
    }
    rays
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use rem_mom::excitation::PlaneWave;

    #[test]
    fn h_field_orthogonal_to_e() {
        let wave = PlaneWave { theta_inc: 0.0, phi_inc: 0.0, pol: "theta".to_string() };
        let (e, h) = incident_fields(&wave, 10.0, &[0.0, 0.0, 0.0]);
        // E · H* should be zero (orthogonal polarisations)
        let dot: Complex64 = (0..3).map(|i| e[i] * h[i].conj()).sum();
        assert!(dot.norm() < 1e-12, "E·H = {:?}", dot);
    }

    #[test]
    fn aperture_rays_nonempty() {
        use rem_mom::surface_mesh::{SurfaceMesh, TriFace};
        use rem_mom::surface_mesh::tri_geometry;

        let nodes = vec![[0.0f64,0.0,0.0],[0.1,0.0,0.0],[0.0,0.1,0.0]];
        let (c,n,a) = tri_geometry(&nodes[0],&nodes[1],&nodes[2]);
        let surf = SurfaceMesh {
            nodes, faces: vec![TriFace{nodes:[0,1,2],centroid:c,normal:n,area:a}],
            edges: vec![], boundary_edges: vec![],
        };
        let wave = PlaneWave { theta_inc:0.0, phi_inc:0.0, pol:"theta".to_string() };
        let rays = launch_aperture_rays(&wave, &surf, 20.0, 1e4, 1e9);
        assert!(!rays.is_empty(), "aperture should generate rays");
    }
}
