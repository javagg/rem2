//! Post-processing: RCS computation and CSV output.

use crate::surface_mesh::SurfaceMesh;
use num_complex::Complex64;
use rem_core::RemResult;
use std::f64::consts::PI;
use std::path::Path;

/// Compute bistatic RCS and write `{output_dir}/postpro/rcs.csv`.
///
/// For pulse basis, `currents[m]` is the surface current density on face m [A/m].
/// For RWG basis, `currents[n]` is the coefficient of the n-th RWG basis function.
///
/// Far-field formula (pulse approximation):
/// F(r̂) = ∫_S J_s(r') exp(jk r̂·r') dS'
/// σ_bi(r̂) = 4π|F|² / |E_inc|²   [m²]
pub fn write_rcs(
    output_dir: &Path,
    freq: f64,
    currents: &[Complex64],
    surf: &SurfaceMesh,
    k: f64,
    theta_deg: &[f64],
    phi_deg: &[f64],
) -> RemResult<()> {
    let path = output_dir.join("postpro").join("rcs.csv");
    let write_header = !path.exists();

    let mut file = std::fs::OpenOptions::new()
        .create(true).append(true).open(&path)?;

    use std::io::Write;
    if write_header {
        writeln!(file, "Freq (GHz),Theta (deg),Phi (deg),RCS (dBsm)")?;
    }

    let freq_ghz = freq / 1.0e9;

    for &phi_d in phi_deg {
        for &theta_d in theta_deg {
            let theta = theta_d.to_radians();
            let phi   = phi_d.to_radians();

            // Observation unit vector r̂
            let rx = theta.sin() * phi.cos();
            let ry = theta.sin() * phi.sin();
            let rz = theta.cos();

            // Far-field integral F = Σ_m J_m * area_m * exp(jk r̂·r_m)
            let fx: Complex64 = currents.iter().zip(surf.faces.iter())
                .map(|(&jm, face)| {
                    let phase = k * (rx*face.centroid[0] + ry*face.centroid[1] + rz*face.centroid[2]);
                    jm * Complex64::new(0.0, phase).exp() * face.area
                })
                .sum();

            // |F|² = fx² (x-polarised current assumed)
            let f_sq = fx.norm_sqr();
            let rcs_m2 = 4.0 * PI * f_sq;
            let rcs_dbsm = if rcs_m2 > 1e-40 {
                10.0 * rcs_m2.log10()
            } else {
                -999.9
            };

            writeln!(file, "{:.6e},{:.1},{:.1},{:.4}",
                freq_ghz, theta_d, phi_d, rcs_dbsm)?;
        }
    }

    Ok(())
}
