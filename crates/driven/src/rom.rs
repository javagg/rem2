//! Snapshot-based Reduced-Order Model (ROM) for the frequency-domain driven solver.
//!
//! # Algorithm
//!
//! Given the parametric system  A(ω) x(ω) = b,  where A(ω) = K − k(ω)² M + losses(ω),
//! build a reduced basis V ∈ ℂ^{n×r} from `r` full solves at expansion frequencies,
//! then for every other frequency evaluate the cheap reduced system:
//!
//!   A_r(ω) = V† A(ω) V  ∈ ℂ^{r×r}
//!   b_r    = V† b        ∈ ℂ^r
//!   x_r(ω) = A_r(ω)⁻¹ b_r
//!   x(ω)  ≈ V x_r(ω)
//!
//! Orthonormalization uses modified Gram-Schmidt so the reduced system is
//! well-conditioned even when snapshot frequencies are close together.
//!
//! # Accuracy
//!
//! The approximation is exact at expansion frequencies and interpolates between them.
//! For a smooth S-parameter sweep (no sharp resonances within the band), `r` = 4–8
//! expansion points typically give < 1 % error in |S11|.
//! Near sharp resonances, increase `r` or fall back to the full sweep.

use nalgebra::{DMatrix, DVector};
use num_complex::Complex64;

/// Orthonormal column basis V ∈ ℂ^{n×r} built from snapshots.
pub struct RomBasis {
    /// Columns of V (each column is a unit-norm orthogonal vector of length n).
    pub cols: Vec<Vec<Complex64>>,
    /// Number of full-system DOFs.
    pub n: usize,
}

impl RomBasis {
    /// Build an orthonormal basis from a set of snapshot solution vectors using
    /// modified Gram-Schmidt.  Linearly dependent vectors (norm < tol) are skipped.
    pub fn from_snapshots(snapshots: Vec<Vec<Complex64>>, tol: f64) -> Self {
        let n = snapshots.first().map(|v| v.len()).unwrap_or(0);
        let mut cols: Vec<Vec<Complex64>> = Vec::new();

        for mut v in snapshots {
            // Orthogonalize against existing basis vectors
            for q in &cols {
                let dot: Complex64 = q.iter().zip(v.iter()).map(|(a, b)| a.conj() * b).sum();
                for (vi, qi) in v.iter_mut().zip(q.iter()) {
                    *vi -= dot * qi;
                }
            }
            // Normalize
            let norm: f64 = v.iter().map(|x| x.norm_sqr()).sum::<f64>().sqrt();
            if norm < tol {
                continue; // linearly dependent — skip
            }
            let inv_norm = Complex64::new(1.0 / norm, 0.0);
            for vi in &mut v {
                *vi *= inv_norm;
            }
            cols.push(v);
        }

        RomBasis { cols, n }
    }

    /// Number of basis vectors (reduced dimension r).
    pub fn r(&self) -> usize {
        self.cols.len()
    }

    /// Project a full vector: b_r = V† b  (length r).
    pub fn project_rhs(&self, b: &[Complex64]) -> Vec<Complex64> {
        self.cols.iter().map(|q| {
            q.iter().zip(b.iter()).map(|(a, x)| a.conj() * x).sum()
        }).collect()
    }

    /// Form the reduced matrix: A_r = V† A V  (r×r dense).
    ///
    /// Rather than forming A explicitly, we accept a closure that applies A to a vector.
    pub fn project_matrix_mv<F>(&self, matvec: F) -> DMatrix<Complex64>
    where
        F: Fn(&[Complex64]) -> Vec<Complex64>,
    {
        let r = self.r();
        let mut ar = DMatrix::<Complex64>::zeros(r, r);
        for j in 0..r {
            let av = matvec(&self.cols[j]);
            for i in 0..r {
                let dot: Complex64 = self.cols[i].iter().zip(av.iter()).map(|(a, b)| a.conj() * b).sum();
                ar[(i, j)] = dot;
            }
        }
        ar
    }

    /// Expand a reduced solution x_r (length r) back to full DOFs.
    pub fn expand(&self, x_r: &[Complex64]) -> Vec<Complex64> {
        let mut x = vec![Complex64::ZERO; self.n];
        for (j, &xj) in x_r.iter().enumerate() {
            for (xi, &vji) in x.iter_mut().zip(self.cols[j].iter()) {
                *xi += xj * vji;
            }
        }
        x
    }
}

/// Solve the reduced system A_r x_r = b_r using LU decomposition (nalgebra).
/// Returns None if the system is singular.
pub fn solve_reduced(ar: DMatrix<Complex64>, br: Vec<Complex64>) -> Option<Vec<Complex64>> {
    let b_dvec = DVector::from_vec(br);
    let lu = ar.lu();
    let x_dvec = lu.solve(&b_dvec)?;
    Some(x_dvec.iter().copied().collect())
}

/// Apply a dense complex matrix to a vector: y = A x.
pub fn dense_matvec(a: &DMatrix<Complex64>, x: &[Complex64]) -> Vec<Complex64> {
    let n = x.len();
    let mut y = vec![Complex64::ZERO; n];
    for i in 0..n {
        for j in 0..n {
            y[i] += a[(i, j)] * x[j];
        }
    }
    y
}

/// Choose `r` expansion frequencies uniformly spread over [f_min, f_max] (log scale if span > 10×).
pub fn choose_expansion_freqs(f_min: f64, f_max: f64, r: usize) -> Vec<f64> {
    if r == 0 { return Vec::new(); }
    if r == 1 { return vec![(f_min + f_max) * 0.5]; }
    let log_scale = f_max / f_min > 10.0 && f_min > 0.0;
    (0..r).map(|i| {
        let t = i as f64 / (r - 1) as f64;
        if log_scale {
            f_min * (f_max / f_min).powf(t)
        } else {
            f_min + t * (f_max - f_min)
        }
    }).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Solve a small 3×3 system exactly with the full solver,
    /// then verify that a ROM basis built from 3 snapshots (= full rank) reproduces it exactly.
    #[test]
    fn test_rom_exact_for_full_rank_basis() {
        // System: A x = b  with  A = diag(1+j, 2, 3+2j),  b = [1, 1, 1]
        let n = 3usize;
        let a_diag = [
            Complex64::new(1.0, 1.0),
            Complex64::new(2.0, 0.0),
            Complex64::new(3.0, 2.0),
        ];
        // Exact solution: x[i] = b[i] / a_diag[i]
        let b = vec![Complex64::new(1.0, 0.0); n];
        let x_exact: Vec<Complex64> = b.iter().zip(a_diag.iter()).map(|(bi, ai)| bi / ai).collect();

        // Build 3 "snapshots" = standard basis vectors (full rank) × exact solution
        let snap0 = vec![x_exact[0], Complex64::ZERO, Complex64::ZERO];
        let snap1 = vec![Complex64::ZERO, x_exact[1], Complex64::ZERO];
        let snap2 = vec![Complex64::ZERO, Complex64::ZERO, x_exact[2]];
        let basis = RomBasis::from_snapshots(vec![snap0, snap1, snap2], 1e-14);
        assert_eq!(basis.r(), 3);

        // Build full dense matrix A
        let mut a_mat = DMatrix::<Complex64>::zeros(n, n);
        for i in 0..n { a_mat[(i, i)] = a_diag[i]; }

        // Project and solve
        let b_r = basis.project_rhs(&b);
        let a_r = basis.project_matrix_mv(|v| dense_matvec(&a_mat, v));
        let x_r = solve_reduced(a_r, b_r).expect("ROM solve should succeed");
        let x_approx = basis.expand(&x_r);

        for i in 0..n {
            let err = (x_approx[i] - x_exact[i]).norm();
            assert!(err < 1e-10, "Component {i}: error={err:.2e}");
        }
    }

    /// Verify that Gram-Schmidt produces an orthonormal set.
    #[test]
    fn test_gram_schmidt_orthonormality() {
        let snaps = vec![
            vec![Complex64::new(1.0, 0.0), Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0)],
            vec![Complex64::new(1.0, 0.0), Complex64::new(-1.0, 0.0), Complex64::new(0.0, 0.0)],
            vec![Complex64::new(0.0, 0.0), Complex64::new(0.0, 0.0), Complex64::new(1.0, 1.0)],
        ];
        let basis = RomBasis::from_snapshots(snaps, 1e-12);
        let r = basis.r();
        assert_eq!(r, 3);
        // Check orthonormality: <q_i, q_j> = δ_{ij}
        for i in 0..r {
            for j in 0..r {
                let dot: Complex64 = basis.cols[i].iter().zip(basis.cols[j].iter())
                    .map(|(a, b)| a.conj() * b).sum();
                if i == j {
                    assert!((dot.re - 1.0).abs() < 1e-12, "diagonal {i}: {dot}");
                    assert!(dot.im.abs() < 1e-12, "diagonal {i}: {dot}");
                } else {
                    assert!(dot.norm() < 1e-12, "off-diag ({i},{j}): {dot}");
                }
            }
        }
    }

    /// Verify expansion frequency selection covers endpoints.
    #[test]
    fn test_expansion_freqs_endpoints() {
        let freqs = choose_expansion_freqs(1e9, 3e9, 4);
        assert_eq!(freqs.len(), 4);
        assert!((freqs[0] - 1e9).abs() < 1.0);
        assert!((freqs[3] - 3e9).abs() < 1.0);
    }
}
