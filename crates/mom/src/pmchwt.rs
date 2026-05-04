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
    use crate::assemble::{assemble_efie_rwg_medium, assemble_mfie_rwg_k};

    let bases = generate_rwg_bases(surf);
    let n = bases.len();
    if n == 0 {
        return Err(rem_core::RemError::Mesh(
            "PMCHWT: No RWG bases found — check surface mesh".to_string()
        ));
    }

    let omega  = 2.0 * PI * freq;
    let k1     = omega / C0;
    let k2     = mat.wave_number(omega);
    let eta1   = (MU0 / EPS0).sqrt();
    let eta2   = mat.impedance();

    // T matrices via validated assemble_efie_rwg_medium (proper singularity treatment,
    // correct -jωμ sign convention matching zmn_efie_rwg)
    let t1 = assemble_efie_rwg_medium(surf, &bases, k1, omega, MU0, EPS0, quad)?;
    let t2 = assemble_efie_rwg_medium(surf, &bases, k2, omega, MU0 * mat.mu_r, EPS0 * mat.eps_r, quad)?;
    // K matrices via validated assemble_mfie_rwg_k
    let k1m = assemble_mfie_rwg_k(surf, &bases, k1, quad)?;
    let k2m = assemble_mfie_rwg_k(surf, &bases, k2, quad)?;

    // Build 2N×2N PMCHWT block matrix.
    // Peterson "Computational Methods for EM" (2001) §8.5 eqs (8.46)-(8.47):
    //   E equation: (T₁+T₂)J + (K₁+K₂)M = -E_i^tan
    //   H equation: (K₁+K₂)J − (T₁/η₁²+T₂/η₂²)M = -H_i^tan
    let two_n = 2 * n;
    let mut z = DMatrix::<Complex64>::zeros(two_n, two_n);

    for i in 0..n {
        for j in 0..n {
            let t_sum   = t1[(i,j)] + t2[(i,j)];
            let k_sum   = k1m[(i,j)] + k2m[(i,j)];
            let t_h_sum = t1[(i,j)] / Complex64::new(eta1 * eta1, 0.0)
                        + t2[(i,j)] / Complex64::new(eta2 * eta2, 0.0);

            z[(i,   j  )] =  t_sum;    // [0,0]: T₁+T₂
            z[(i,   j+n)] =  k_sum;    // [0,1]: K₁+K₂
            z[(i+n, j  )] =  k_sum;    // [1,0]: K₁+K₂
            z[(i+n, j+n)] = -t_h_sum;  // [1,1]: -(T₁/η₁²+T₂/η₂²)
        }
    }

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
        SurfaceMesh { nodes, faces, edges, boundary_edges, face_attrs: vec![0, 0], global_node_ids: vec![] }
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
