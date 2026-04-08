//! PMCHWT formulation for dielectric targets in MoM.
//!
//! Poggio-Miller-Chang-Harrington-Wu-Tsai (PMCHWT) integral equation for a
//! homogeneous dielectric body with relative permittivity ε_r and permeability μ_r.
//!
//! # Physical model
//!
//! A closed dielectric surface Γ separates two regions:
//! - Region 1 (exterior): free space (ε₀, μ₀)
//! - Region 2 (interior): homogeneous dielectric (ε₀ ε_r, μ₀ μ_r)
//!
//! Unknown surface quantities: electric current **J** and magnetic current **M**.
//!
//! # PMCHWT system (2N × 2N)
//!
//! With N RWG basis functions, the unknowns are expanded as J = Σ aₙ fₙ and
//! M = Σ bₙ fₙ.  Galerkin testing with the same RWG set gives:
//!
//! ```text
//! ┌                    ┐ ┌   ┐   ┌           ┐
//! │ T₁+T₂    K₁+K₂    │ │ a │   │ -⟨f, E_i⟩ │
//! │                    │ │   │ = │           │
//! │ -(K₁+K₂) T₁/η₁+T₂/η₂ │ │ b │   │ -⟨f, H_i⟩ │
//! └                    ┘ └   ┘   └           ┘
//! ```
//!
//! where Tₚ is the EFIE operator for medium p and Kₚ is the MFIE operator,
//! η₁ = η₀ (free space), η₂ = η₀ √(μ_r / ε_r).
//!
//! # Reference
//! Rao, Wilton, Glisson (1982); Peterson et al. "Computational Methods for
//! Electromagnetics" §8.5.

use nalgebra::DMatrix;
use num_complex::Complex64;
use rem_core::{RemResult, EPS0, MU0, C0};
use std::f64::consts::PI;

use crate::surface_mesh::SurfaceMesh;
use crate::basis::rwg::{RwgBasis, generate_rwg_bases};
use crate::quadrature::TriQuad;
use crate::assemble::{lu_solve, gmres_solve};
use crate::excitation::PlaneWave;

/// Material parameters for the interior dielectric region.
#[derive(Debug, Clone, Copy)]
pub struct DielectricMaterial {
    /// Relative permittivity ε_r  (≥ 1)
    pub eps_r: f64,
    /// Relative permeability μ_r  (≥ 1)
    pub mu_r: f64,
}

impl DielectricMaterial {
    pub fn new(eps_r: f64, mu_r: f64) -> Self {
        Self { eps_r, mu_r }
    }

    /// Wave number in medium: k = ω √(ε μ)
    pub fn wave_number(&self, omega: f64) -> f64 {
        omega * (EPS0 * self.eps_r * MU0 * self.mu_r).sqrt()
    }

    /// Wave impedance η = √(μ / ε)
    pub fn impedance(&self) -> f64 {
        (MU0 * self.mu_r / (EPS0 * self.eps_r)).sqrt()
    }
}

// ---------------------------------------------------------------------------
// PMCHWT system assembly
// ---------------------------------------------------------------------------

/// Assemble the 2N×2N PMCHWT impedance matrix and 2N RHS for a dielectric body.
///
/// Returns `(z_pmchwt, rhs_pmchwt)` where `z_pmchwt` is 2N×2N and `rhs_pmchwt`
/// has length 2N.  The first N unknowns are J coefficients; the last N are M.
///
/// `mat` — interior dielectric material
/// `freq` — frequency [Hz]
/// `wave` — incident plane wave
/// `quad` — quadrature rule
/// `fast_solver` — "Direct" (LU) or "GMRES"
pub fn assemble_pmchwt(
    surf: &SurfaceMesh,
    mat: DielectricMaterial,
    freq: f64,
    wave: &PlaneWave,
    quad: &TriQuad,
) -> RemResult<(DMatrix<Complex64>, Vec<Complex64>)> {
    let bases = generate_rwg_bases(surf);
    let n = bases.len();
    if n == 0 {
        return Err(rem_core::RemError::Mesh(
            "PMCHWT: No RWG bases found — check surface mesh".to_string()
        ));
    }

    let omega  = 2.0 * PI * freq;
    let k1     = omega / C0;                    // exterior (free space)
    let k2     = mat.wave_number(omega);        // interior
    let eta1   = (MU0 / EPS0).sqrt();           // η₁ = η₀
    let eta2   = mat.impedance();               // η₂ = η₀ √(μ_r/ε_r)

    // Assemble T (EFIE-like) and K (MFIE-like) for both media
    let t1 = assemble_t_matrix(surf, &bases, k1, omega, MU0, EPS0, quad)?;
    let t2 = assemble_t_matrix(surf, &bases, k2, omega, MU0 * mat.mu_r, EPS0 * mat.eps_r, quad)?;
    let k1m = assemble_k_matrix(surf, &bases, k1, quad)?;
    let k2m = assemble_k_matrix(surf, &bases, k2, quad)?;

    // Build 2N×2N block matrix
    let two_n = 2 * n;
    let mut z = DMatrix::<Complex64>::zeros(two_n, two_n);

    // Block (0,0): T₁ + T₂
    // Block (0,1): K₁ + K₂
    // Block (1,0): -(K₁ + K₂)
    // Block (1,1): T₁/η₁² + T₂/η₂²  (normalization for H equation)
    for i in 0..n {
        for j in 0..n {
            let t_sum  = t1[(i,j)] + t2[(i,j)];
            let k_sum  = k1m[(i,j)] + k2m[(i,j)];

            // Row block 0 (E equation, test fn f_i)
            z[(i, j)]     = t_sum;          // [0,0] block
            z[(i, j+n)]   = k_sum;          // [0,1] block

            // Row block 1 (H equation, test fn f_i)
            let t_h_sum = t1[(i,j)] / Complex64::new(eta1 * eta1, 0.0)
                        + t2[(i,j)] / Complex64::new(eta2 * eta2, 0.0);
            z[(i+n, j)]   = -k_sum;         // [1,0] block
            z[(i+n, j+n)] = t_h_sum;        // [1,1] block
        }
    }

    // Build 2N RHS: [-⟨f_m, E_inc⟩, -⟨f_m, H_inc⟩]
    let rhs = build_pmchwt_rhs(surf, &bases, k1, wave, quad);

    Ok((z, rhs))
}

/// Solve the PMCHWT system and return (J_coeffs, M_coeffs), each of length N.
pub fn solve_pmchwt(
    surf: &SurfaceMesh,
    mat: DielectricMaterial,
    freq: f64,
    wave: &PlaneWave,
    quad: &TriQuad,
    fast_solver: &str,
) -> RemResult<(Vec<Complex64>, Vec<Complex64>)> {
    let (z, rhs) = assemble_pmchwt(surf, mat, freq, wave, quad)?;
    let n = surf.edges.len(); // = N (RWG count)
    let two_n = 2 * n;

    let x = match fast_solver.to_uppercase().as_str() {
        "GMRES" => gmres_solve(&z, &rhs)?,
        _       => lu_solve(&z, &rhs)?,
    };

    let j_coeffs = x[..two_n/2].to_vec();
    let m_coeffs = x[two_n/2..].to_vec();

    Ok((j_coeffs, m_coeffs))
}

// ---------------------------------------------------------------------------
// T matrix (EFIE operator, medium with wave number k, permittivity eps, permeability mu)
// ---------------------------------------------------------------------------

/// Assemble EFIE-type T matrix for medium with parameters (k, mu, eps).
///
/// T[m,n] = jωμ ⟨f_m, L(f_n)⟩ where L is the EFIE integral operator.
fn assemble_t_matrix(
    surf: &SurfaceMesh,
    bases: &[RwgBasis],
    k: f64,
    omega: f64,
    mu: f64,
    eps: f64,
    quad: &TriQuad,
) -> RemResult<DMatrix<Complex64>> {
    use crate::green::green3d;

    let n = bases.len();
    let omega_mu      = omega * mu;
    let inv_omega_eps = 1.0 / (omega * eps);
    let jomega_mu     = Complex64::new(0.0, omega_mu);

    let mut t = DMatrix::<Complex64>::zeros(n, n);

    for ni in 0..n {
        let bn = &bases[ni];
        for mi in 0..n {
            let bm = &bases[mi];
            let mut val = Complex64::ZERO;

            for &(m_face, m_plus) in &[(bm.plus_face, true), (bm.minus_face, false)] {
                for &(n_face, n_plus) in &[(bn.plus_face, true), (bn.minus_face, false)] {
                    let face_m = &surf.faces[m_face];
                    let face_n = &surf.faces[n_face];
                    let div_n = bn.divergence(surf, n_plus);

                    for (bm_pt, &wm) in quad.bary.iter().zip(quad.weights.iter()) {
                        let rm = crate::quadrature::TriQuad::global_point(bm_pt, face_m, &surf.nodes);
                        let fm = bm.eval(&rm, surf, m_plus);

                        for (bn_pt, &wn) in quad.bary.iter().zip(quad.weights.iter()) {
                            let rn = crate::quadrature::TriQuad::global_point(bn_pt, face_n, &surf.nodes);
                            let fn_ = bn.eval(&rn, surf, n_plus);
                            let g = green3d(&rm, &rn, k);

                            let dot_ff = fm[0]*fn_[0] + fm[1]*fn_[1] + fm[2]*fn_[2];
                            let div_m = bm.divergence(surf, m_plus);

                            let integrand = g * (dot_ff - inv_omega_eps / omega_mu * div_m * div_n);
                            val += integrand * (wm * wn * 4.0 * face_m.area * face_n.area);
                        }
                    }
                }
            }

            t[(mi, ni)] = jomega_mu * val;
        }
    }

    Ok(t)
}

// ---------------------------------------------------------------------------
// K matrix (MFIE operator)
// ---------------------------------------------------------------------------

/// Assemble MFIE-type K matrix for medium with wave number k.
///
/// K[m,n] = ⟨f_m, (1/2 δ_{mn} + K_op) f_n⟩ where K_op uses the
/// curl of the Green's function.
fn assemble_k_matrix(
    surf: &SurfaceMesh,
    bases: &[RwgBasis],
    k: f64,
    quad: &TriQuad,
) -> RemResult<DMatrix<Complex64>> {
    let n = bases.len();
    let mut km = DMatrix::<Complex64>::zeros(n, n);

    for ni in 0..n {
        let bn = &bases[ni];
        for mi in 0..n {
            let bm = &bases[mi];

            // Identity term: δ_{mn}/2 * ⟨f_m, f_n⟩ = δ_{mn}/2 * overlap
            let identity_term = if bm.edge_idx == bn.edge_idx {
                let mut overlap = 0.0f64;
                for &(face_idx, in_plus) in &[(bm.plus_face, true), (bm.minus_face, false)] {
                    let face = &surf.faces[face_idx];
                    for (b_pt, &w) in quad.bary.iter().zip(quad.weights.iter()) {
                        let r = crate::quadrature::TriQuad::global_point(b_pt, face, &surf.nodes);
                        let f = bm.eval(&r, surf, in_plus);
                        overlap += (f[0]*f[0] + f[1]*f[1] + f[2]*f[2]) * (w * 2.0 * face.area);
                    }
                }
                Complex64::new(0.5 * overlap, 0.0)
            } else {
                Complex64::ZERO
            };

            // Curl-Green integral
            let mut curl_term = Complex64::ZERO;
            for &(m_face, m_plus) in &[(bm.plus_face, true), (bm.minus_face, false)] {
                for &(n_face, n_plus) in &[(bn.plus_face, true), (bn.minus_face, false)] {
                    let face_m = &surf.faces[m_face];
                    let face_n = &surf.faces[n_face];
                    let nm = face_m.normal;

                    for (bm_pt, &wm) in quad.bary.iter().zip(quad.weights.iter()) {
                        let rm = crate::quadrature::TriQuad::global_point(bm_pt, face_m, &surf.nodes);
                        let fm = bm.eval(&rm, surf, m_plus);

                        for (bn_pt, &wn) in quad.bary.iter().zip(quad.weights.iter()) {
                            let rn = crate::quadrature::TriQuad::global_point(bn_pt, face_n, &surf.nodes);
                            let fn_ = bn.eval(&rn, surf, n_plus);
                            let grad_g = green_gradient_k(&rm, &rn, k);

                            let fn_c = [
                                Complex64::new(fn_[0], 0.0),
                                Complex64::new(fn_[1], 0.0),
                                Complex64::new(fn_[2], 0.0),
                            ];
                            let curl_gfn = cross_c(&grad_g, &fn_c);

                            let nm_c = [
                                Complex64::new(nm[0], 0.0),
                                Complex64::new(nm[1], 0.0),
                                Complex64::new(nm[2], 0.0),
                            ];
                            let n_x_curl = cross_c(&nm_c, &curl_gfn);
                            let dot_val = Complex64::new(
                                fm[0]*n_x_curl[0].re + fm[1]*n_x_curl[1].re + fm[2]*n_x_curl[2].re,
                                fm[0]*n_x_curl[0].im + fm[1]*n_x_curl[1].im + fm[2]*n_x_curl[2].im,
                            );

                            curl_term += dot_val
                                       * (wm * wn * 4.0 * face_m.area * face_n.area);
                        }
                    }
                }
            }

            km[(mi, ni)] = identity_term + curl_term;
        }
    }

    Ok(km)
}

// ---------------------------------------------------------------------------
// RHS: PMCHWT excitation
// ---------------------------------------------------------------------------

/// Build the 2N PMCHWT RHS vector.
///
/// The first N entries are  -⟨f_m, E_tan_inc⟩  (E-field equation).
/// The last  N entries are  -⟨f_m, H_tan_inc⟩  (H-field equation).
fn build_pmchwt_rhs(
    surf: &SurfaceMesh,
    bases: &[RwgBasis],
    k: f64,
    wave: &PlaneWave,
    quad: &TriQuad,
) -> Vec<Complex64> {
    let n = bases.len();
    let mut rhs = vec![Complex64::ZERO; 2 * n];

    let kh = wave.k_hat();
    let eh = wave.e_hat();

    // H_inc direction: H_hat = k̂ × ê / η₀
    let eta0 = (MU0 / EPS0).sqrt();
    let hh = [
        (kh[1]*eh[2] - kh[2]*eh[1]) / eta0,
        (kh[2]*eh[0] - kh[0]*eh[2]) / eta0,
        (kh[0]*eh[1] - kh[1]*eh[0]) / eta0,
    ];

    for (mi, bm) in bases.iter().enumerate() {
        let mut ve = Complex64::ZERO; // E-equation entry
        let mut vh = Complex64::ZERO; // H-equation entry

        for &(face_idx, in_plus) in &[(bm.plus_face, true), (bm.minus_face, false)] {
            let face = &surf.faces[face_idx];
            for (b_pt, &w) in quad.bary.iter().zip(quad.weights.iter()) {
                let r = crate::quadrature::TriQuad::global_point(b_pt, face, &surf.nodes);
                let fm = bm.eval(&r, surf, in_plus);

                // Phase factor for incident plane wave
                let phase = k * (kh[0]*r[0] + kh[1]*r[1] + kh[2]*r[2]);
                let phasor = Complex64::new(0.0, -phase).exp();

                // E_inc(r) = ê * exp(-jk k̂·r)
                let e_inc = [
                    Complex64::new(eh[0], 0.0) * phasor,
                    Complex64::new(eh[1], 0.0) * phasor,
                    Complex64::new(eh[2], 0.0) * phasor,
                ];
                // H_inc(r) = ĥ * exp(-jk k̂·r)
                let h_inc = [
                    Complex64::new(hh[0], 0.0) * phasor,
                    Complex64::new(hh[1], 0.0) * phasor,
                    Complex64::new(hh[2], 0.0) * phasor,
                ];

                let dot_fe: Complex64 = fm[0]*e_inc[0] + fm[1]*e_inc[1] + fm[2]*e_inc[2];
                let dot_fh: Complex64 = fm[0]*h_inc[0] + fm[1]*h_inc[1] + fm[2]*h_inc[2];

                ve += dot_fe * (w * 2.0 * face.area);
                vh += dot_fh * (w * 2.0 * face.area);
            }
        }

        rhs[mi]     = -ve;
        rhs[mi + n] = -vh;
    }

    rhs
}

// ---------------------------------------------------------------------------
// Helpers (duplicated locally to avoid coupling with assemble.rs internals)
// ---------------------------------------------------------------------------

fn green_gradient_k(r: &[f64; 3], rp: &[f64; 3], k: f64) -> [Complex64; 3] {
    use crate::green::green3d;
    let dr = [r[0]-rp[0], r[1]-rp[1], r[2]-rp[2]];
    let dist = (dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]).sqrt();
    if dist < 1e-14 { return [Complex64::ZERO; 3]; }
    let g = green3d(r, rp, k);
    let factor = g * Complex64::new(1.0/dist, k) / dist;
    [factor * dr[0], factor * dr[1], factor * dr[2]]
}

fn cross_c(a: &[Complex64; 3], b: &[Complex64; 3]) -> [Complex64; 3] {
    [
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0],
    ]
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::surface_mesh::{SurfaceMesh, TriFace};

    /// Build a minimal two-triangle surface mesh (shared edge) for testing.
    fn two_tri_surf() -> SurfaceMesh {
        use crate::surface_mesh::SharedEdge;
        let nodes = vec![
            [0.0_f64, 0.0, 0.0],
            [1.0,     0.0, 0.0],
            [0.5,     1.0, 0.0],
            [-0.5,    1.0, 0.0],
        ];
        let faces = vec![
            TriFace {
                nodes: [0, 1, 2],
                centroid: [0.5, 1.0/3.0, 0.0],
                normal: [0.0, 0.0, 1.0],
                area: 0.5,
            },
            TriFace {
                nodes: [0, 2, 3],
                centroid: [-0.5/3.0 + 0.5/3.0, 2.0/3.0, 0.0],
                normal: [0.0, 0.0, 1.0],
                area: 0.5,
            },
        ];
        // Interior edge shared by both faces: [0, 2]
        let edges = vec![SharedEdge {
            nodes: [0, 2],
            plus_face: 0,
            minus_face: 1,
            length: 1.118,  // approx sqrt(1² + 1²) / sqrt(2)
        }];
        let boundary_edges = vec![[0usize, 1], [1, 2], [2, 3], [3, 0]];
        SurfaceMesh { nodes, faces, edges, boundary_edges }
    }

    #[test]
    fn pmchwt_matrix_assembles_finite() {
        let surf = two_tri_surf();
        let mat  = DielectricMaterial::new(4.0, 1.0); // glass-like
        let wave = PlaneWave { theta_inc: 0.0, phi_inc: 0.0, pol: "theta".to_string() };
        let quad = crate::quadrature::TriQuad::new(5);

        let freq = 1e9;
        let result = assemble_pmchwt(&surf, mat, freq, &wave, &quad);
        assert!(result.is_ok(), "PMCHWT assembly failed: {:?}", result.err());

        let (z, rhs) = result.unwrap();
        // System should be 2N × 2N where N = #RWG = 1 edge → 2×2
        assert_eq!(z.nrows(), 2);
        assert_eq!(z.ncols(), 2);
        assert_eq!(rhs.len(), 2);

        // All entries should be finite
        for i in 0..z.nrows() {
            for j in 0..z.ncols() {
                assert!(z[(i,j)].re.is_finite() && z[(i,j)].im.is_finite(),
                    "z[{i},{j}] not finite");
            }
        }
        for (i, &r) in rhs.iter().enumerate() {
            assert!(r.re.is_finite() && r.im.is_finite(), "rhs[{i}] not finite");
        }
    }

    #[test]
    fn pmchwt_rhs_nonzero() {
        let surf = two_tri_surf();
        let mat  = DielectricMaterial::new(2.25, 1.0);
        let wave = PlaneWave { theta_inc: std::f64::consts::PI/4.0, phi_inc: 0.0, pol: "theta".to_string() };
        let quad = crate::quadrature::TriQuad::new(5);

        let (_, rhs) = assemble_pmchwt(&surf, mat, 3e9, &wave, &quad).unwrap();
        let rhs_norm: f64 = rhs.iter().map(|x| x.norm_sqr()).sum::<f64>().sqrt();
        assert!(rhs_norm > 1e-20, "RHS should be non-zero for oblique incidence");
    }

    #[test]
    fn dielectric_material_wave_number() {
        let mat = DielectricMaterial::new(4.0, 1.0);
        let omega = 2.0 * std::f64::consts::PI * 1e9;
        let k = mat.wave_number(omega);
        let k_free = omega / C0;
        // k_medium = sqrt(eps_r) * k_free for non-magnetic
        let expected = (4.0_f64).sqrt() * k_free;
        assert!((k - expected).abs() / expected < 1e-10,
            "k_medium = {k:.4e}, expected {expected:.4e}");
    }

    #[test]
    fn pmchwt_solve_small() {
        let surf = two_tri_surf();
        let mat  = DielectricMaterial::new(2.0, 1.0);
        let wave = PlaneWave { theta_inc: 0.0, phi_inc: 0.0, pol: "theta".to_string() };
        let quad = crate::quadrature::TriQuad::new(5);

        let result = solve_pmchwt(&surf, mat, 1e9, &wave, &quad, "Direct");
        assert!(result.is_ok(), "solve_pmchwt failed: {:?}", result.err());
        let (j, m) = result.unwrap();
        assert_eq!(j.len(), 1);
        assert_eq!(m.len(), 1);
        assert!(j[0].re.is_finite() && m[0].re.is_finite());
    }
}
