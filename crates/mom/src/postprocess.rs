//! Post-processing: RCS computation and CSV output.

use crate::surface_mesh::SurfaceMesh;
use crate::basis::rwg::{generate_rwg_bases, RwgBasis};
use crate::green::green3d;
use num_complex::Complex64;
use rem_core::{NearFieldPoint, RemResult, ETA0, MU0};
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
// Near-field export
// ---------------------------------------------------------------------------

/// Compute near-field (E and H) at the centroids of all surface faces.
///
/// Uses the electric field integral representation for PEC:
///   E(r) = -jωμ A(r) - ∇(∇·A(r)) / (jωε₀)
/// where A(r) = Σ_n I_n ∫ G(r,r') f_n(r') dS'
///
/// We evaluate A at face centroids using centroid quadrature and approximate
/// E ≈ -jωμ A (the scalar-potential term is neglected at the centroid level).
/// H is obtained from the surface current via n̂ × J_s / η₀.
pub fn compute_near_field(
    currents: &[Complex64],
    surf: &SurfaceMesh,
    k: f64,
) -> Vec<NearFieldPoint> {
    let omega = k * rem_core::C0;
    let jw_mu = Complex64::new(0.0, omega * MU0);

    let bases = generate_rwg_bases(surf);
    let n_faces = surf.faces.len();
    let mut points: Vec<NearFieldPoint> = Vec::with_capacity(n_faces);

    for face_idx in 0..n_faces {
        let face = &surf.faces[face_idx];
        let r = &face.centroid;

        // Vector potential A = Σ_n I_n ∫ G(r,r') f_n(r') dS'
        let mut ax = Complex64::ZERO;
        let mut ay = Complex64::ZERO;
        let mut az = Complex64::ZERO;

        for (n, base) in bases.iter().enumerate() {
            let i_n = currents[n];
            for &(fi, in_plus) in &[(base.plus_face, true), (base.minus_face, false)] {
                let src_face = &surf.faces[fi];
                let g = green3d(r, &src_face.centroid, k);
                if g.norm() < 1e-300 {
                    continue;
                }
                let f_n = base.eval(&src_face.centroid, surf, in_plus);
                let contrib = i_n * g * src_face.area;
                ax += f_n[0] * contrib;
                ay += f_n[1] * contrib;
                az += f_n[2] * contrib;
            }
        }

        // E ≈ -jωμ A
        let ex = -jw_mu * ax;
        let ey = -jw_mu * ay;
        let ez = -jw_mu * az;

        // Surface current at this face: J_s = Σ_{n: n touches face} I_n · f_n(r)
        let mut jx = Complex64::ZERO;
        let mut jy = Complex64::ZERO;
        let mut jz = Complex64::ZERO;
        for (n, base) in bases.iter().enumerate() {
            let i_n = currents[n];
            for &(fi, in_plus) in &[(base.plus_face, true), (base.minus_face, false)] {
                if fi == face_idx {
                    let f_n = base.eval(r, surf, in_plus);
                    jx += i_n * f_n[0];
                    jy += i_n * f_n[1];
                    jz += i_n * f_n[2];
                }
            }
        }
        // H = n̂ × Re(J_s) / η₀  (real-valued for PEC surface)
        let nn = &face.normal;
        let hx = (nn[1]*jz.re - nn[2]*jy.re) / ETA0;
        let hy = (nn[2]*jx.re - nn[0]*jz.re) / ETA0;
        let hz = (nn[0]*jy.re - nn[1]*jx.re) / ETA0;

        points.push(NearFieldPoint::from_complex(
            r[0], r[1], r[2],
            ex, ey, ez,
            Complex64::new(hx, 0.0),
            Complex64::new(hy, 0.0),
            Complex64::new(hz, 0.0),
        ));
    }

    points
}

/// Write near-field data to CSV.
pub fn write_near_field_csv(
    output_dir: &Path,
    points: &[NearFieldPoint],
    output_file: Option<&str>,
) -> RemResult<()> {
    let filename = output_file.unwrap_or("postpro/near_field.csv");
    let path = output_dir.join(filename);
    rem_core::write_near_field_csv(&path, points)
}

// ---------------------------------------------------------------------------
// Near-field at arbitrary probe points (S-param / RWG current mode)
// ---------------------------------------------------------------------------

/// Compute the E-field vector at a list of arbitrary probe points.
///
/// Uses the electric-field integral representation (centroid quadrature):
///   E(r) ≈ −jωμ₀ Σₙ Iₙ Σ_{face ∈ support(n)} G(r, r'_c) f_n(r'_c) A_face
///
/// `currents` are the RWG basis-function coefficients (length = `bases.len()`).
/// Returns one `[Complex64; 3]` per probe point (Ex, Ey, Ez).
pub fn compute_e_at_probes(
    surf: &SurfaceMesh,
    bases: &[RwgBasis],
    currents: &[Complex64],
    probes: &[[f64; 3]],
    k: f64,
) -> Vec<[Complex64; 3]> {
    let omega   = k * rem_core::C0;
    let jw_mu   = Complex64::new(0.0, omega * MU0);

    probes.iter().map(|r| {
        let mut ax = Complex64::ZERO;
        let mut ay = Complex64::ZERO;
        let mut az = Complex64::ZERO;

        for (n, base) in bases.iter().enumerate() {
            let i_n = currents[n];
            for &(fi, in_plus) in &[(base.plus_face, true), (base.minus_face, false)] {
                let src_face = &surf.faces[fi];
                let g = green3d(r, &src_face.centroid, k);
                let f_v = base.eval(&src_face.centroid, surf, in_plus);
                let contrib = i_n * g * src_face.area;
                ax += f_v[0] * contrib;
                ay += f_v[1] * contrib;
                az += f_v[2] * contrib;
            }
        }
        [-jw_mu * ax, -jw_mu * ay, -jw_mu * az]
    }).collect()
}

/// Write near-field probe results to CSV.
///
/// Columns: `freq_hz, x, y, z, Ex_re, Ex_im, Ey_re, Ey_im, Ez_re, Ez_im, |E| (dBV/m)`
pub fn write_probe_e_field_csv(
    path: &Path,
    probe_xyz:  &[[f64; 3]],         // probe coordinates
    freq_e:     &[(f64, Vec<[Complex64; 3]>)], // (freq_hz, e_per_probe)
) -> RemResult<()> {
    use std::io::Write;
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() { std::fs::create_dir_all(parent)?; }
    }
    let mut f = std::fs::File::create(path)?;
    writeln!(f, "Freq (GHz),x,y,z,Ex_re,Ex_im,Ey_re,Ey_im,Ez_re,Ez_im,|E| (dBV/m)")?;
    for &(freq_hz, ref e_vals) in freq_e {
        for (idx, e) in e_vals.iter().enumerate() {
            let [x, y, z] = probe_xyz[idx];
            let [ex, ey, ez] = *e;
            let e_mag = (ex.norm_sqr() + ey.norm_sqr() + ez.norm_sqr()).sqrt();
            let e_db = if e_mag > 1e-300 { 20.0 * e_mag.log10() } else { -999.0 };
            writeln!(f,
                "{:.9e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.4}",
                freq_hz / 1e9, x, y, z,
                ex.re, ex.im, ey.re, ey.im, ez.re, ez.im, e_db,
            )?;
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
            global_node_ids: vec![],
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
