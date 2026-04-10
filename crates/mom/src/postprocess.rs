//! Post-processing: RCS computation and CSV output.

use crate::surface_mesh::SurfaceMesh;
use num_complex::Complex64;
use rem_core::{RemResult, ETA0};
use std::f64::consts::PI;
use std::path::Path;

/// Compute bistatic RCS [m²] at given (theta, phi) angles.
/// Returns a 2-D array: `result[i_theta][i_phi]`.
///
/// For x-polarized pulse-basis currents (scalar J_m on each face):
///   N_x(r̂) = Σ_m J_m · A_m · exp(jk r̂·r_m)
///   σ(r̂) = k²η₀²/(4π) · |N_x|²   [m²]
pub fn rcs_pattern(
    currents: &[Complex64],
    surf: &SurfaceMesh,
    k: f64,
    theta_deg: &[f64],
    phi_deg: &[f64],
) -> Vec<Vec<f64>> {
    let prefactor = k * k * ETA0 * ETA0 / (4.0 * PI);
    theta_deg.iter().map(|&theta_d| {
        let theta = theta_d.to_radians();
        phi_deg.iter().map(|&phi_d| {
            let phi = phi_d.to_radians();
            let rx = theta.sin() * phi.cos();
            let ry = theta.sin() * phi.sin();
            let rz = theta.cos();
            // Radiation vector (x-component)
            let nx: Complex64 = currents.iter().zip(surf.faces.iter())
                .map(|(&jm, face)| {
                    let phase = k * (rx*face.centroid[0] + ry*face.centroid[1] + rz*face.centroid[2]);
                    jm * Complex64::new(0.0, phase).exp() * face.area
                })
                .sum();
            prefactor * nx.norm_sqr()
        }).collect()
    }).collect()
}


/// Write surface current distribution to VTK Legacy ASCII format.
///
/// Generates a `.vtk` file with the triangular surface mesh and
/// the following cell data (one value per triangle):
/// - `J_mag` : |J_m| magnitude [A/m]
/// - `J_real`: Re(J_m) component
/// - `J_imag`: Im(J_m) component
///
/// Compatible with ParaView and VisIt.
pub fn write_surface_vtk(
    path: &Path,
    currents: &[Complex64],
    surf: &SurfaceMesh,
) -> RemResult<()> {
    use std::io::Write;
    let mut f = std::fs::File::create(path)?;

    let n_pts = surf.nodes.len();
    let n_cells = surf.faces.len();

    // Header
    writeln!(f, "# vtk DataFile Version 3.0")?;
    writeln!(f, "MoM surface current")?;
    writeln!(f, "ASCII")?;
    writeln!(f, "DATASET UNSTRUCTURED_GRID")?;
    writeln!(f)?;

    // Points
    writeln!(f, "POINTS {} float", n_pts)?;
    for &[x,y,z] in &surf.nodes {
        writeln!(f, "{:.8e} {:.8e} {:.8e}", x, y, z)?;
    }
    writeln!(f)?;

    // Cells: each row = "3 i0 i1 i2" for a triangle
    writeln!(f, "CELLS {} {}", n_cells, n_cells * 4)?;
    for face in &surf.faces {
        writeln!(f, "3 {} {} {}", face.nodes[0], face.nodes[1], face.nodes[2])?;
    }
    writeln!(f)?;

    // Cell types: 5 = VTK_TRIANGLE
    writeln!(f, "CELL_TYPES {}", n_cells)?;
    for _ in 0..n_cells {
        writeln!(f, "5")?;
    }
    writeln!(f)?;

    // Cell data
    writeln!(f, "CELL_DATA {}", n_cells)?;

    // |J| magnitude
    writeln!(f, "SCALARS J_mag float 1")?;
    writeln!(f, "LOOKUP_TABLE default")?;
    for &j in currents {
        writeln!(f, "{:.8e}", j.norm())?;
    }

    // Re(J)
    writeln!(f, "SCALARS J_real float 1")?;
    writeln!(f, "LOOKUP_TABLE default")?;
    for &j in currents {
        writeln!(f, "{:.8e}", j.re)?;
    }

    // Im(J)
    writeln!(f, "SCALARS J_imag float 1")?;
    writeln!(f, "LOOKUP_TABLE default")?;
    for &j in currents {
        writeln!(f, "{:.8e}", j.im)?;
    }

    Ok(())
}


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

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::surface_mesh::{tri_geometry, TriFace};

    fn mini_surf() -> (SurfaceMesh, Vec<Complex64>) {
        let nodes = vec![[0.0,0.0,0.0],[1.0,0.0,0.0],[0.0,1.0,0.0]];
        let (c,n,a) = tri_geometry(&nodes[0],&nodes[1],&nodes[2]);
        let surf = SurfaceMesh {
            nodes,
            faces: vec![TriFace{nodes:[0,1,2],centroid:c,normal:n,area:a}],
            edges: vec![],
            boundary_edges: vec![],
            face_attrs: vec![0],
        };
        let currents = vec![Complex64::new(1.0, 0.5)];
        (surf, currents)
    }

    #[test]
    fn vtk_output_creates_file() {
        let (surf, currents) = mini_surf();
        let tmp = std::env::temp_dir().join("test_surface_current.vtk");
        write_surface_vtk(&tmp, &currents, &surf).expect("VTK write failed");
        assert!(tmp.exists(), "VTK file not created");
        let content = std::fs::read_to_string(&tmp).unwrap();
        assert!(content.contains("vtk DataFile"), "Missing VTK header");
        assert!(content.contains("J_mag"), "Missing J_mag scalar");
        assert!(content.contains("CELL_DATA"), "Missing cell data section");
        let _ = std::fs::remove_file(&tmp);
    }
}
