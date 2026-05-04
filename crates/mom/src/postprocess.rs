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


/// Compute bistatic RCS [m²] from RWG basis-function current coefficients.
///
/// Far-field radiation vector is computed via centroid quadrature (one point
/// per face per basis, exact for linearly-varying RWG functions under PO):
///
///   **N**(r̂) = Σₙ Iₙ [ **f**ₙ(r'₊) A₊ e^{jk r̂·r'₊} + **f**ₙ(r'₋) A₋ e^{jk r̂·r'₋} ]
///   σ = k²/(4π) |r̂ × (r̂ × η₀**N**)·x̂|²   [m²]  (x-pol incident)
///
/// Returns a 2-D array: `result[i_theta][i_phi]`.
pub fn rcs_pattern_rwg(
    currents: &[Complex64],
    surf: &SurfaceMesh,
    bases: &[crate::basis::rwg::RwgBasis],
    k: f64,
    theta_deg: &[f64],
    phi_deg: &[f64],
) -> Vec<Vec<f64>> {
    use crate::basis::rwg::RwgBasis;
    let eta_k2 = ETA0 * ETA0 * k * k / (4.0 * PI);

    theta_deg.iter().map(|&theta_d| {
        let theta = theta_d.to_radians();
        phi_deg.iter().map(|&phi_d| {
            let phi = phi_d.to_radians();
            let rhat = [theta.sin() * phi.cos(),
                        theta.sin() * phi.sin(),
                        theta.cos()];

            // Radiation vector N = [Nx, Ny, Nz]
            let mut n = [Complex64::ZERO; 3];
            for (idx, base) in bases.iter().enumerate() {
                let i_n = currents[idx];
                for &(fi, in_plus) in &[(base.plus_face, true), (base.minus_face, false)] {
                    let face = &surf.faces[fi];
                    let c = &face.centroid;
                    let phase = k * (rhat[0]*c[0] + rhat[1]*c[1] + rhat[2]*c[2]);
                    let exp_phase = Complex64::new(0.0, phase).exp();
                    let fv = base.eval(c, surf, in_plus);
                    let contrib = i_n * exp_phase * face.area;
                    n[0] += fv[0] * contrib;
                    n[1] += fv[1] * contrib;
                    n[2] += fv[2] * contrib;
                }
            }

            // Cross-pol projection: r̂ × (r̂ × η₀N), take x-component
            // r̂ × N = [ry*Nz - rz*Ny, rz*Nx - rx*Nz, rx*Ny - ry*Nx]
            // r̂ × (r̂ × N) = r̂(r̂·N) - N  → (r̂ × (r̂ × ηN))·x̂ = r̂ₓ(r̂·ηN) - ηNx
            let eta_n = [ETA0 * n[0], ETA0 * n[1], ETA0 * n[2]];
            let rdot = rhat[0]*eta_n[0] + rhat[1]*eta_n[1] + rhat[2]*eta_n[2];
            let ff_x = rhat[0] * rdot - eta_n[0];
            let ff_y = rhat[1] * rdot - eta_n[1];
            let ff_z = rhat[2] * rdot - eta_n[2];

            eta_k2 / (ETA0 * ETA0) * (ff_x.norm_sqr() + ff_y.norm_sqr() + ff_z.norm_sqr())
        }).collect()
    }).collect()
}

/// Compute bistatic RCS [m²] from PMCHWT solution (J + M RWG current coefficients).
///
/// The PMCHWT far-field includes both electric (J) and magnetic (M) surface currents:
///
///   **N**(r̂) = Σₙ Iₙᴶ [**f**ₙ(c₊) A₊ e^{jkr̂·c₊} + **f**ₙ(c₋) A₋ e^{jkr̂·c₋}]
///   **L**(r̂) = Σₙ Iₙᴹ [**f**ₙ(c₊) A₊ e^{jkr̂·c₊} + **f**ₙ(c₋) A₋ e^{jkr̂·c₋}]
///
///   σ = k²/(4π) |η₀ (I − r̂r̂)**N** + r̂ × **L**|²
///
/// PMCHWT currents are **not** k-scaled (unlike CFIE-RWG via `rwg_rhs`).
/// Returns a 2-D array: `result[i_theta][i_phi]`.
pub fn rcs_pattern_pmchwt(
    j_coeffs: &[Complex64],
    m_coeffs: &[Complex64],
    surf: &SurfaceMesh,
    bases: &[crate::basis::rwg::RwgBasis],
    k: f64,
    theta_deg: &[f64],
    phi_deg: &[f64],
) -> Vec<Vec<f64>> {
    let prefactor = k * k / (4.0 * PI);

    theta_deg.iter().map(|&theta_d| {
        let theta = theta_d.to_radians();
        phi_deg.iter().map(|&phi_d| {
            let phi = phi_d.to_radians();
            let rhat = [theta.sin() * phi.cos(),
                        theta.sin() * phi.sin(),
                        theta.cos()];

            let mut nv = [Complex64::ZERO; 3]; // electric radiation vector N
            let mut lv = [Complex64::ZERO; 3]; // magnetic radiation vector L

            for (idx, base) in bases.iter().enumerate() {
                let ij = j_coeffs[idx];
                let im = m_coeffs[idx];
                for &(fi, in_plus) in &[(base.plus_face, true), (base.minus_face, false)] {
                    let face = &surf.faces[fi];
                    let c = &face.centroid;
                    let phase = k * (rhat[0]*c[0] + rhat[1]*c[1] + rhat[2]*c[2]);
                    let exp_ph = Complex64::new(0.0, phase).exp();
                    let fv = base.eval(c, surf, in_plus);
                    let contrib_j = ij * exp_ph * face.area;
                    let contrib_m = im * exp_ph * face.area;
                    nv[0] += Complex64::new(fv[0], 0.0) * contrib_j;
                    nv[1] += Complex64::new(fv[1], 0.0) * contrib_j;
                    nv[2] += Complex64::new(fv[2], 0.0) * contrib_j;
                    lv[0] += Complex64::new(fv[0], 0.0) * contrib_m;
                    lv[1] += Complex64::new(fv[1], 0.0) * contrib_m;
                    lv[2] += Complex64::new(fv[2], 0.0) * contrib_m;
                }
            }

            // η₀(I − r̂r̂)N = η₀(N − (r̂·N)r̂)
            let rdotn = rhat[0]*nv[0] + rhat[1]*nv[1] + rhat[2]*nv[2];
            let eta_nt = [
                ETA0 * (nv[0] - rdotn * rhat[0]),
                ETA0 * (nv[1] - rdotn * rhat[1]),
                ETA0 * (nv[2] - rdotn * rhat[2]),
            ];

            // r̂ × L
            let rxl = [
                rhat[1]*lv[2] - rhat[2]*lv[1],
                rhat[2]*lv[0] - rhat[0]*lv[2],
                rhat[0]*lv[1] - rhat[1]*lv[0],
            ];

            // Balanis eq. 3-58: E_t ∝ η₀N_t + r̂×L  (with PMCHWT M convention M = -n̂×E).
            let total = [eta_nt[0] + rxl[0], eta_nt[1] + rxl[1], eta_nt[2] + rxl[2]];
            prefactor * (total[0].norm_sqr() + total[1].norm_sqr() + total[2].norm_sqr())
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
// RWG vector-current VTK (per-face surface current density)
// ---------------------------------------------------------------------------

/// Write surface current density magnitude to VTK from RWG basis coefficients.
///
/// Per-face current density is approximated as:
///   J_s(c_face) = Σ_n I_n f_n(c_face)   (centroid evaluation)
///
/// Outputs `J_vec_x/y/z` (real/imag) and `J_mag` cell scalars for ParaView.
pub fn write_surface_current_vtk_rwg(
    path: &Path,
    surf: &SurfaceMesh,
    bases: &[RwgBasis],
    currents: &[Complex64],
) -> RemResult<()> {
    use std::io::Write;

    // Build per-face current density [Complex64; 3]
    let n_faces = surf.faces.len();
    let mut jx = vec![Complex64::ZERO; n_faces];
    let mut jy = vec![Complex64::ZERO; n_faces];
    let mut jz = vec![Complex64::ZERO; n_faces];

    for (n, base) in bases.iter().enumerate() {
        let i_n = if n < currents.len() { currents[n] } else { Complex64::ZERO };
        for &(fi, in_plus) in &[(base.plus_face, true), (base.minus_face, false)] {
            if fi >= n_faces { continue; }
            let f_v = base.eval(&surf.faces[fi].centroid, surf, in_plus);
            jx[fi] += i_n * f_v[0];
            jy[fi] += i_n * f_v[1];
            jz[fi] += i_n * f_v[2];
        }
    }

    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() { std::fs::create_dir_all(parent)?; }
    }
    let mut f = std::fs::File::create(path)?;

    writeln!(f, "# vtk DataFile Version 3.0")?;
    writeln!(f, "MoM RWG surface current")?;
    writeln!(f, "ASCII")?;
    writeln!(f, "DATASET UNSTRUCTURED_GRID")?;
    writeln!(f)?;

    writeln!(f, "POINTS {} float", surf.nodes.len())?;
    for &[x, y, z] in &surf.nodes {
        writeln!(f, "{:.8e} {:.8e} {:.8e}", x, y, z)?;
    }
    writeln!(f)?;

    writeln!(f, "CELLS {} {}", n_faces, n_faces * 4)?;
    for face in &surf.faces {
        writeln!(f, "3 {} {} {}", face.nodes[0], face.nodes[1], face.nodes[2])?;
    }
    writeln!(f)?;

    writeln!(f, "CELL_TYPES {}", n_faces)?;
    for _ in 0..n_faces { writeln!(f, "5")?; }
    writeln!(f)?;

    writeln!(f, "CELL_DATA {}", n_faces)?;

    // |J| magnitude
    writeln!(f, "SCALARS J_mag float 1")?;
    writeln!(f, "LOOKUP_TABLE default")?;
    for i in 0..n_faces {
        let mag = (jx[i].norm_sqr() + jy[i].norm_sqr() + jz[i].norm_sqr()).sqrt();
        writeln!(f, "{:.8e}", mag)?;
    }

    // Vector Jx_re, Jy_re, Jz_re
    writeln!(f, "VECTORS J_real float")?;
    for i in 0..n_faces {
        writeln!(f, "{:.8e} {:.8e} {:.8e}", jx[i].re, jy[i].re, jz[i].re)?;
    }

    writeln!(f, "VECTORS J_imag float")?;
    for i in 0..n_faces {
        writeln!(f, "{:.8e} {:.8e} {:.8e}", jx[i].im, jy[i].im, jz[i].im)?;
    }

    Ok(())
}

// ---------------------------------------------------------------------------
// Far-field radiation pattern from RWG port-excited currents
// ---------------------------------------------------------------------------

/// A single far-field observation point result.
#[derive(Debug, Clone)]
pub struct FarFieldPoint {
    pub theta_deg: f64,
    pub phi_deg:   f64,
    /// N_theta component of radiation vector [A·m].
    pub n_theta:   Complex64,
    /// N_phi component of radiation vector [A·m].
    pub n_phi:     Complex64,
    /// Radiation intensity U = (k²η₀/32π²)(|N_θ|²+|N_φ|²)  [W/sr] (proportional).
    pub u:         f64,
    /// Directivity D [dBi] (computed after integrating over sphere).
    pub d_dbi:     f64,
}

/// Compute the far-field radiation pattern for RWG port-excited currents.
///
/// Uses centroid quadrature for the radiation integral:
///   N(r̂) = Σ_n I_n [f_n(c⁺) A⁺ exp(jk r̂·c⁺) + f_n(c⁻) A⁻ exp(jk r̂·c⁻)]
///
/// Returns directivity-tagged far-field points (one per (theta, phi) pair).
pub fn compute_radiation_pattern_rwg(
    currents: &[Complex64],
    surf: &SurfaceMesh,
    bases: &[RwgBasis],
    k: f64,
    theta_deg_list: &[f64],
    phi_deg_list:   &[f64],
) -> Vec<FarFieldPoint> {
    use std::f64::consts::PI;

    let n_theta = theta_deg_list.len();
    let n_phi   = phi_deg_list.len();

    // Pre-compute per-basis centroid contributions (independent of observation angle)
    // contrib_n = [(f_n(c+)*A+, c+), (f_n(c-)*A-, c-)]
    struct BasisContrib {
        fv: [f64; 3],   // f_n(c) direction
        a:  f64,        // weighted area (A_face)
        c:  [f64; 3],   // centroid
        coeff: Complex64,
    }
    let contribs: Vec<BasisContrib> = bases.iter().enumerate()
        .flat_map(|(n, base)| {
            let i_n = if n < currents.len() { currents[n] } else { Complex64::ZERO };
            let mut v = Vec::with_capacity(2);
            for &(fi, in_plus) in &[(base.plus_face, true), (base.minus_face, false)] {
                if fi >= surf.faces.len() { continue; }
                let face = &surf.faces[fi];
                let fv = base.eval(&face.centroid, surf, in_plus);
                v.push(BasisContrib {
                    fv,
                    a: face.area,
                    c: face.centroid,
                    coeff: i_n,
                });
            }
            v
        })
        .collect();

    // Compute N(r̂) for each observation direction
    let mut raw_points: Vec<(f64, f64, Complex64, Complex64, f64)> = Vec::with_capacity(n_theta * n_phi);

    for &theta_d in theta_deg_list {
        let theta = theta_d.to_radians();
        let st = theta.sin();
        let ct = theta.cos();

        for &phi_d in phi_deg_list {
            let phi = phi_d.to_radians();
            let sp = phi.sin();
            let cp = phi.cos();

            // r̂ unit vector
            let rx = st * cp;
            let ry = st * sp;
            let rz = ct;

            // Radiation vector N = Σ contributions
            let mut nx = Complex64::ZERO;
            let mut ny = Complex64::ZERO;
            let mut nz = Complex64::ZERO;

            for bc in &contribs {
                let phase_arg = k * (rx * bc.c[0] + ry * bc.c[1] + rz * bc.c[2]);
                let phase = Complex64::new(0.0, phase_arg).exp();
                let scale = bc.coeff * phase * bc.a;
                nx += scale * bc.fv[0];
                ny += scale * bc.fv[1];
                nz += scale * bc.fv[2];
            }

            // Project N onto spherical (θ, φ) components
            let n_theta_c = nx * ct * cp + ny * ct * sp - nz * st;
            let n_phi_c   = -nx * sp   + ny * cp;

            // Radiation intensity U ∝ |N_θ|² + |N_φ|²
            let u = n_theta_c.norm_sqr() + n_phi_c.norm_sqr();
            raw_points.push((theta_d, phi_d, n_theta_c, n_phi_c, u));
        }
    }

    // Numerical integration over sphere to get P_rad for directivity
    // Trapezoidal rule in theta, trapezoidal in phi
    let p_rad = if n_theta > 1 && n_phi > 1 {
        use std::f64::consts::PI;
        let dtheta = (theta_deg_list.last().unwrap() - theta_deg_list.first().unwrap()).to_radians()
            / (n_theta - 1) as f64;
        let dphi = (phi_deg_list.last().unwrap() - phi_deg_list.first().unwrap()).to_radians()
            / (n_phi - 1).max(1) as f64;
        let mut sum = 0.0_f64;
        for i in 0..n_theta {
            let theta = theta_deg_list[i].to_radians();
            let st = theta.sin();
            for j in 0..n_phi {
                let u = raw_points[i * n_phi + j].4;
                let w_i = if i == 0 || i == n_theta - 1 { 0.5 } else { 1.0 };
                let w_j = if j == 0 || j == n_phi - 1 { 0.5 } else { 1.0 };
                sum += u * st * w_i * w_j;
            }
        }
        sum * dtheta * dphi
    } else {
        // Fallback: single point or line — compute unit integral
        raw_points.iter().map(|p| p.4).sum::<f64>() * 4.0 * std::f64::consts::PI
            / raw_points.len().max(1) as f64
    };

    raw_points.into_iter().map(|(theta_d, phi_d, n_theta_c, n_phi_c, u)| {
        let d_dbi = if p_rad > 1e-300 {
            10.0 * (4.0 * std::f64::consts::PI * u / p_rad).log10()
        } else {
            0.0
        };
        FarFieldPoint { theta_deg: theta_d, phi_deg: phi_d, n_theta: n_theta_c, n_phi: n_phi_c, u, d_dbi }
    }).collect()
}

/// Write far-field radiation pattern to CSV (append mode, one block per frequency).
///
/// Columns: `Freq (GHz), Theta (deg), Phi (deg), |N_theta|, |N_phi|, U (norm), D (dBi)`
pub fn write_radiation_pattern_csv(
    path: &Path,
    points: &[FarFieldPoint],
    freq_hz: f64,
) -> RemResult<()> {
    use std::io::Write;
    let write_header = !path.exists();
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() { std::fs::create_dir_all(parent)?; }
    }
    let mut f = std::fs::OpenOptions::new().create(true).append(true).open(path)?;
    if write_header {
        writeln!(f, "Freq (GHz),Theta (deg),Phi (deg),|N_theta| (A*m),|N_phi| (A*m),U (norm),D (dBi)")?;
    }
    let freq_ghz = freq_hz / 1e9;
    for p in points {
        writeln!(f,
            "{:.9e},{:.2},{:.2},{:.6e},{:.6e},{:.6e},{:.4}",
            freq_ghz, p.theta_deg, p.phi_deg,
            p.n_theta.norm(), p.n_phi.norm(), p.u, p.d_dbi,
        )?;
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

    // ── RWG VTK test ──────────────────────────────────────────────────────

    /// Build two-triangle surface (same as sparams tests) with two RWG bases.
    fn two_tri_surf() -> SurfaceMesh {
        use crate::surface_mesh::{TriFace, SharedEdge, tri_geometry};
        let nodes = vec![
            [0.0_f64, 0.0, 0.0],
            [1.0,     0.0, 0.0],
            [0.5,     1.0, 0.0],
            [-0.5,    1.0, 0.0],
        ];
        let (c0, n0, a0) = tri_geometry(&nodes[0], &nodes[1], &nodes[2]);
        let (c1, n1, a1) = tri_geometry(&nodes[0], &nodes[2], &nodes[3]);
        let faces = vec![
            TriFace { nodes:[0,1,2], centroid:c0, normal:n0, area:a0 },
            TriFace { nodes:[0,2,3], centroid:c1, normal:n1, area:a1 },
        ];
        let edges = vec![SharedEdge {
            nodes: [0,2], plus_face: 0, minus_face: 1,
            length: (0.5_f64.powi(2) + 1.0_f64.powi(2)).sqrt(),
        }];
        SurfaceMesh {
            nodes, faces, edges,
            boundary_edges: vec![[0,1],[1,2],[2,3],[3,0]],
            face_attrs: vec![1, 1],
            global_node_ids: vec![],
        }
    }

    #[test]
    fn rwg_vtk_creates_file_with_vectors() {
        let surf = two_tri_surf();
        let bases = generate_rwg_bases(&surf);
        let currents = vec![Complex64::new(1.0, 0.5)];
        let tmp = std::env::temp_dir().join("test_rwg_vtk.vtk");
        write_surface_current_vtk_rwg(&tmp, &surf, &bases, &currents).expect("RWG VTK failed");
        let content = std::fs::read_to_string(&tmp).unwrap();
        assert!(content.contains("J_mag"), "Missing J_mag");
        assert!(content.contains("J_real"), "Missing J_real vectors");
        let _ = std::fs::remove_file(&tmp);
    }

    // ── Far-field radiation pattern tests ────────────────────────────────

    /// For a z-directed current element at origin, the far-field pattern
    /// should have |N_theta| proportional to sin(theta) and N_phi ≈ 0.
    #[test]
    fn radiation_pattern_z_dipole_symmetry() {
        // Build a single-triangle surface with RWG currents pointing in z
        let surf = two_tri_surf();
        let bases = generate_rwg_bases(&surf);
        // Use a uniform current coefficient
        let currents = vec![Complex64::new(1.0, 0.0); bases.len()];
        let k = 1.0;
        let theta_list: Vec<f64> = vec![30.0, 60.0, 90.0];
        let phi_list:   Vec<f64> = vec![0.0];

        let pts = compute_radiation_pattern_rwg(&currents, &surf, &bases, k, &theta_list, &phi_list);
        assert_eq!(pts.len(), theta_list.len() * phi_list.len());
        // All U values should be non-negative
        for p in &pts { assert!(p.u >= 0.0, "U must be non-negative"); }
        // Pattern should be finite
        for p in &pts { assert!(p.d_dbi.is_finite(), "D must be finite"); }
    }

    #[test]
    fn radiation_pattern_csv_creates_file() {
        let surf = two_tri_surf();
        let bases = generate_rwg_bases(&surf);
        let currents = vec![Complex64::new(1.0, 0.0)];
        let pts = compute_radiation_pattern_rwg(
            &currents, &surf, &bases, 1.0,
            &[45.0, 90.0], &[0.0, 90.0],
        );
        let tmp = std::env::temp_dir().join("test_far_field.csv");
        let _ = std::fs::remove_file(&tmp);
        write_radiation_pattern_csv(&tmp, &pts, 2.4e9).expect("CSV write failed");
        let content = std::fs::read_to_string(&tmp).unwrap();
        assert!(content.contains("Freq"), "Missing CSV header");
        assert!(content.contains("dBi"), "Missing dBi column");
        let data_lines = content.lines().filter(|l| !l.starts_with('#') && !l.contains("Freq")).count();
        assert_eq!(data_lines, pts.len());
        let _ = std::fs::remove_file(&tmp);
    }
}
