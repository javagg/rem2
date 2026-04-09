//! Near-to-far-field transform for driven (frequency-domain) solver.
//!
//! Uses the Kirchhoff/Stratton-Chu surface integral approximation:
//!
//!   F(r̂) = ∫_S **E**(r') e^{jk r̂·r'} dS'
//!
//! where **E** = −∇φ is the recovered electric field from the scalar FEM solution,
//! r̂ = (sin θ cos φ, sin θ sin φ, cos θ) is the observation direction,
//! k = ω/c is the wavenumber at the peak frequency.
//!
//! For a real-valued near field (taking Re(phi)), the phase e^{jk r̂·r'} makes F complex.
//! The amplitude |F(r̂)|² is proportional to the radiated power density in direction r̂.
//!
//! Gain is normalized to isotropic radiator (dBi):
//!   D(r̂) = 4π |F(r̂)|² / ∫ |F(r̂)|² dΩ
//!   G_dBi = 10 log10(D(r̂))

use rem_config::FarFieldConfig;
use rem_core::constants::C0;
use rem_mesh::{RemMesh, ElementKind};
use rem_electrostatic::postprocess::gradient_recovery;
use std::f64::consts::PI;

/// Far-field pattern point.
#[derive(Debug, Clone)]
pub struct FarFieldPoint {
    /// Elevation angle θ [degrees] (0 = +z axis, 180 = −z axis)
    pub theta_deg: f64,
    /// Azimuth angle φ [degrees] (0 = +x axis)
    pub phi_deg: f64,
    /// |F|² (unnormalized, proportional to power density in this direction)
    pub power_linear: f64,
    /// Directivity relative to isotropic (dBi)
    pub gain_dbi: f64,
}

/// Compute the far-field radiation pattern from a driven solver solution.
///
/// - `phi`: real-part nodal potential at the peak response frequency
/// - `freq_hz`: the frequency of the solution
/// - `cfg`: far-field configuration (angles, attributes)
pub fn compute_far_field(
    mesh: &RemMesh,
    phi: &[f64],
    freq_hz: f64,
    cfg: &FarFieldConfig,
) -> Vec<FarFieldPoint> {
    if phi.is_empty() {
        return Vec::new();
    }

    let omega = 2.0 * PI * freq_hz;
    let k = omega / C0;

    // Recover E-field (−∇φ) at all nodes
    let e_field = gradient_recovery(phi, mesh);
    // Negate: E = −∇φ
    let e_field: Vec<[f64; 3]> = e_field.iter().map(|e| [-e[0], -e[1], -e[2]]).collect();

    // Collect boundary element contributions
    // For each boundary element, compute the centroid, area-weighted normal E contribution
    let mut belem_data: Vec<([f64; 3], [f64; 3], f64)> = Vec::new(); // (centroid, E_avg, area)

    for belem in &mesh.boundary_elements {
        // Filter by attributes if specified
        if !cfg.attributes.is_empty() && !cfg.attributes.contains(&belem.tag) {
            continue;
        }

        match belem.kind {
            ElementKind::Line2 => {
                if belem.node_ids.len() < 2 { continue; }
                let n0 = &mesh.nodes[belem.node_ids[0]];
                let n1 = &mesh.nodes[belem.node_ids[1]];
                let cx = (n0.x + n1.x) * 0.5;
                let cy = (n0.y + n1.y) * 0.5;
                let cz = (n0.z + n1.z) * 0.5;
                let dx = n1.x - n0.x;
                let dy = n1.y - n0.y;
                let len = (dx * dx + dy * dy).sqrt();
                // Average E at centroid
                let e0 = e_field.get(belem.node_ids[0]).copied().unwrap_or([0.0; 3]);
                let e1 = e_field.get(belem.node_ids[1]).copied().unwrap_or([0.0; 3]);
                let e_avg = [(e0[0]+e1[0])*0.5, (e0[1]+e1[1])*0.5, (e0[2]+e1[2])*0.5];
                belem_data.push(([cx, cy, cz], e_avg, len));
            }
            ElementKind::Tri3 => {
                if belem.node_ids.len() < 3 { continue; }
                let n0 = &mesh.nodes[belem.node_ids[0]];
                let n1 = &mesh.nodes[belem.node_ids[1]];
                let n2 = &mesh.nodes[belem.node_ids[2]];
                let cx = (n0.x + n1.x + n2.x) / 3.0;
                let cy = (n0.y + n1.y + n2.y) / 3.0;
                let cz = (n0.z + n1.z + n2.z) / 3.0;
                // Area via cross product
                let ax = n1.x - n0.x; let ay = n1.y - n0.y; let az = n1.z - n0.z;
                let bx = n2.x - n0.x; let by = n2.y - n0.y; let bz = n2.z - n0.z;
                let cx2 = ay*bz - az*by;
                let cy2 = az*bx - ax*bz;
                let cz2 = ax*by - ay*bx;
                let area = 0.5 * (cx2*cx2 + cy2*cy2 + cz2*cz2).sqrt();
                let e0 = e_field.get(belem.node_ids[0]).copied().unwrap_or([0.0; 3]);
                let e1 = e_field.get(belem.node_ids[1]).copied().unwrap_or([0.0; 3]);
                let e2 = e_field.get(belem.node_ids[2]).copied().unwrap_or([0.0; 3]);
                let e_avg = [
                    (e0[0]+e1[0]+e2[0])/3.0,
                    (e0[1]+e1[1]+e2[1])/3.0,
                    (e0[2]+e1[2]+e2[2])/3.0,
                ];
                belem_data.push(([cx, cy, cz], e_avg, area));
            }
            _ => {} // Skip Quad4 etc. (not common as boundary elements)
        }
    }

    if belem_data.is_empty() {
        log::warn!("[REM] FarField: no boundary elements found for integration. Check FarField.Attributes.");
        return Vec::new();
    }

    log::info!("[REM] FarField: integrating over {} boundary element(s), k={:.4e} rad/m", belem_data.len(), k);

    // Compute F over the grid
    let mut ff_values: Vec<(f64, f64, f64)> = Vec::new(); // (theta_deg, phi_deg, |F|^2)
    let mut total_solid_angle = 0.0;
    let mut weighted_power = 0.0;

    for it in 0..cfg.n_theta {
        let theta = it as f64 * PI / (cfg.n_theta - 1).max(1) as f64;
        let theta_deg = theta.to_degrees();
        let sin_t = theta.sin();
        let cos_t = theta.cos();

        for ip in 0..cfg.n_phi {
            let phi_ang = ip as f64 * 2.0 * PI / cfg.n_phi as f64;
            let phi_deg = phi_ang.to_degrees();

            // Observation unit vector r̂
            let rx = sin_t * phi_ang.cos();
            let ry = sin_t * phi_ang.sin();
            let rz = cos_t;

            // F(r̂) = ∫ E(r') e^{jk r̂·r'} dS'
            // = Σ_e E_e · area_e · e^{jk (r̂·c_e)}
            // where c_e is centroid of element e.
            // F is a 3-vector; compute |F|² = |Fx|² + |Fy|² + |Fz|²
            let mut fx_re = 0.0; let mut fx_im = 0.0;
            let mut fy_re = 0.0; let mut fy_im = 0.0;
            let mut fz_re = 0.0; let mut fz_im = 0.0;

            for &(centroid, e_avg, area) in &belem_data {
                let phase = k * (rx * centroid[0] + ry * centroid[1] + rz * centroid[2]);
                let cos_p = phase.cos();
                let sin_p = phase.sin();
                let w = area;
                // F += E_e * area * e^{jk r̂·r'} = E_e * area * (cos_p + j sin_p)
                fx_re += e_avg[0] * w * cos_p;
                fx_im += e_avg[0] * w * sin_p;
                fy_re += e_avg[1] * w * cos_p;
                fy_im += e_avg[1] * w * sin_p;
                fz_re += e_avg[2] * w * cos_p;
                fz_im += e_avg[2] * w * sin_p;
            }

            let f2 = fx_re*fx_re + fx_im*fx_im
                   + fy_re*fy_re + fy_im*fy_im
                   + fz_re*fz_re + fz_im*fz_im;

            ff_values.push((theta_deg, phi_deg, f2));

            // Solid angle weight for spherical average: dΩ = sin(θ) dθ dφ
            let d_theta = PI / (cfg.n_theta - 1).max(1) as f64;
            let d_phi = 2.0 * PI / cfg.n_phi as f64;
            let d_omega = sin_t * d_theta * d_phi;
            total_solid_angle += d_omega;
            weighted_power += f2 * d_omega;
        }
    }

    // Normalize to dBi
    let avg_power = if total_solid_angle > 1e-300 {
        weighted_power / total_solid_angle
    } else {
        1.0
    };

    ff_values.into_iter().map(|(theta_deg, phi_deg, f2)| {
        let directivity = if avg_power > 1e-300 { f2 / avg_power } else { 0.0 };
        let gain_dbi = if directivity > 1e-300 { 10.0 * directivity.log10() } else { -100.0 };
        FarFieldPoint { theta_deg, phi_deg, power_linear: f2, gain_dbi }
    }).collect()
}

/// Write far-field pattern to CSV.
pub fn write_far_field_csv(out_dir: &str, pattern: &[FarFieldPoint], freq_hz: f64) -> rem_core::RemResult<()> {
    use std::io::Write;
    let path = std::path::Path::new(out_dir).join("far_field.csv");
    let mut f = std::fs::File::create(&path).map_err(rem_core::RemError::Io)?;
    writeln!(f, "freq_hz,theta_deg,phi_deg,power_linear,gain_dbi").map_err(rem_core::RemError::Io)?;
    for pt in pattern {
        writeln!(f, "{:.6e},{:.2},{:.2},{:.6e},{:.4}",
            freq_hz, pt.theta_deg, pt.phi_deg, pt.power_linear, pt.gain_dbi)
            .map_err(rem_core::RemError::Io)?;
    }
    log::info!("[REM] Far-field pattern written to {}", path.display());
    Ok(())
}
