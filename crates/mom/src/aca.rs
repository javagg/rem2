//! Adaptive Cross Approximation (ACA) for MoM impedance matrices.
//!
//! ACA builds a low-rank approximation of the off-diagonal blocks of the MoM
//! Z matrix:  Z_block ≈ U · V^T  where U ∈ ℂ^{m×r}, V ∈ ℂ^{n×r}, r << m,n.
//!
//! Note: The MoM impedance matrix Z is complex *symmetric* (Z = Z^T), not
//! Hermitian.  Therefore we use the non-conjugated transpose throughout.
//!
//! This reduces:
//!  - Assembly from O(N²) full evaluations → O(N r) evaluations  (r ~ log N)
//!  - Storage from O(N²) → O(N r)
//!  - Matrix-vector products from O(N²) → O(N r)  (used in GMRES)
//!
//! # Algorithm (partially-pivoted ACA, Bebendorf 2000)
//!
//! For a block Z[I, J] (|I|=m, |J|=n):
//! 1. Choose pivot row i* = first row of I.
//! 2. Compute the full i*-th row: r_i = Z[i*, J].
//! 3. Find column pivot j* = argmax |r_i(j) − accumulated|.
//! 4. Compute the full j*-th column: c_j = Z[I, j*].
//! 5. Scale: u = c_j / r_i(j*),   v = r_i.
//! 6. Update residual: R -= u v^T.
//! 7. Stopping criterion: ‖u‖ · ‖v‖ < ε · ‖Z_approx‖_F.

use nalgebra::{DMatrix, DVector};
use num_complex::Complex64;

/// Low-rank approximation of a complex matrix: Z ≈ U · V^T  (no conjugation).
///
/// `U` has shape (nrows × rank), `V` has shape (ncols × rank).
#[derive(Debug, Clone)]
pub struct AcaMatrix {
    /// Left factor U  [nrows × rank]
    pub u: DMatrix<Complex64>,
    /// Right factor V [ncols × rank]
    pub v: DMatrix<Complex64>,
    /// Original matrix dimensions
    pub nrows: usize,
    pub ncols: usize,
    /// Achieved approximation rank
    pub rank: usize,
    /// Final Frobenius error estimate
    pub rel_error: f64,
}

impl AcaMatrix {
    /// Build ACA from a fully assembled dense matrix.
    pub fn compress(z: &DMatrix<Complex64>, tol: f64, max_rank: usize) -> Self {
        let m = z.nrows();
        let n = z.ncols();
        let row_fn = |i: usize| -> Vec<Complex64> {
            (0..n).map(|j| z[(i, j)]).collect()
        };
        let col_fn = |j: usize| -> Vec<Complex64> {
            (0..m).map(|i| z[(i, j)]).collect()
        };
        Self::assemble(m, n, &row_fn, &col_fn, tol, max_rank)
    }

    /// Build ACA by evaluating individual rows and columns on demand.
    ///
    /// `row_fn(i)` returns the full i-th row (length n).
    /// `col_fn(j)` returns the full j-th column (length m).
    pub fn assemble(
        m: usize,
        n: usize,
        row_fn: &dyn Fn(usize) -> Vec<Complex64>,
        col_fn: &dyn Fn(usize) -> Vec<Complex64>,
        tol: f64,
        max_rank: usize,
    ) -> Self {
        let max_rank = max_rank.min(m).min(n);

        let mut u_cols: Vec<Vec<Complex64>> = Vec::with_capacity(max_rank);
        let mut v_rows: Vec<Vec<Complex64>> = Vec::with_capacity(max_rank);

        // Residual tracked implicitly via U·V^T
        let mut used_rows = vec![false; m];
        let mut used_cols = vec![false; n];

        let mut frob_sq: f64 = 0.0;

        for _iter in 0..max_rank {
            // Pick next unused pivot row
            let i_star = match used_rows.iter().position(|&u| !u) {
                Some(i) => i,
                None => break,
            };
            used_rows[i_star] = true;

            // Residual row: r̃_i = row_fn(i*) − Σ_k u_k[i*] · v_k  (no conj)
            let mut r_row = row_fn(i_star);
            for k in 0..u_cols.len() {
                let u_ik = u_cols[k][i_star];
                if u_ik == Complex64::ZERO { continue; }
                for j in 0..n {
                    r_row[j] -= u_ik * v_rows[k][j];  // V^T, no conjugation
                }
            }

            // Column pivot: argmax |r̃_i(j)| over unused columns
            let j_star = r_row.iter().enumerate()
                .filter(|(j, _)| !used_cols[*j])
                .max_by(|(_, a), (_, b)| a.norm_sqr().partial_cmp(&b.norm_sqr()).unwrap())
                .map(|(j, _)| j);

            let j_star = match j_star {
                Some(j) => j,
                None => break,
            };
            used_cols[j_star] = true;

            let pivot = r_row[j_star];
            if pivot.norm() < 1e-300 { break; }

            // v_new = r̃_i  (row)
            let v_new = r_row;

            // Residual column: c̃_j = col_fn(j*) − Σ_k u_k · v_k[j*]
            let mut c_col = col_fn(j_star);
            for k in 0..u_cols.len() {
                let v_kj = v_rows[k][j_star];  // no conjugation
                if v_kj == Complex64::ZERO { continue; }
                for i in 0..m {
                    c_col[i] -= u_cols[k][i] * v_kj;
                }
            }

            // u_new = c̃_j / pivot
            let inv_pivot = Complex64::new(1.0, 0.0) / pivot;
            let u_new: Vec<Complex64> = c_col.iter().map(|&x| x * inv_pivot).collect();

            // Stopping criterion
            let u_norm2: f64 = u_new.iter().map(|x| x.norm_sqr()).sum();
            let v_norm2: f64 = v_new.iter().map(|x| x.norm_sqr()).sum();
            let delta_frob_sq = u_norm2 * v_norm2;

            u_cols.push(u_new);
            v_rows.push(v_new);

            frob_sq += delta_frob_sq;

            let rel_err = if frob_sq > 0.0 {
                (delta_frob_sq / frob_sq).sqrt()
            } else {
                1.0
            };

            if rel_err < tol { break; }
        }

        let rank = u_cols.len();

        let mut u = DMatrix::<Complex64>::zeros(m, rank);
        let mut v = DMatrix::<Complex64>::zeros(n, rank);
        for k in 0..rank {
            for i in 0..m { u[(i, k)] = u_cols[k][i]; }
            for j in 0..n { v[(j, k)] = v_rows[k][j]; }
        }

        let rel_error = if frob_sq > 0.0 && rank > 0 {
            let u_n: f64 = (0..m).map(|i| u[(i, rank-1)].norm_sqr()).sum::<f64>().sqrt();
            let v_n: f64 = (0..n).map(|j| v[(j, rank-1)].norm_sqr()).sum::<f64>().sqrt();
            u_n * v_n / frob_sq.sqrt()
        } else {
            0.0
        };

        Self { u, v, nrows: m, ncols: n, rank, rel_error }
    }

    /// Matrix-vector product: y = Z x ≈ U (V^T x).  O(N · rank).
    ///
    /// For complex symmetric matrices use standard transpose (no conjugation).
    pub fn matvec(&self, x: &DVector<Complex64>) -> DVector<Complex64> {
        // tmp = V^T x  (standard transpose, no conjugation)
        let tmp = self.v.transpose() * x;
        &self.u * tmp
    }

    /// Reconstruct full dense matrix Z ≈ U V^T  (no conjugation).
    pub fn to_dense(&self) -> DMatrix<Complex64> {
        &self.u * self.v.transpose()
    }
}

/// Partition the N×N Z matrix into blocks and apply ACA to each far-field block.
///
/// Near-field blocks (|block_i − block_j| ≤ near_thresh) are assembled exactly.
pub fn aca_partition(
    n: usize,
    block_size: usize,
    near_thresh: usize,
    tol: f64,
    max_rank: usize,
    entry_fn: &dyn Fn(usize, usize) -> Complex64,
) -> (Vec<(usize, usize, Complex64)>, Vec<(usize, usize, AcaMatrix)>) {
    let n_blocks = (n + block_size - 1) / block_size;
    let mut near = Vec::new();
    let mut far  = Vec::new();

    for bi in 0..n_blocks {
        let i0 = bi * block_size;
        let i1 = n.min(i0 + block_size);
        let m_block = i1 - i0;

        for bj in 0..n_blocks {
            let j0 = bj * block_size;
            let j1 = n.min(j0 + block_size);
            let n_block = j1 - j0;

            let block_dist = (bi as isize - bj as isize).unsigned_abs();

            if block_dist <= near_thresh {
                for i in i0..i1 {
                    for j in j0..j1 {
                        near.push((i, j, entry_fn(i, j)));
                    }
                }
            } else {
                let row_fn = |li: usize| -> Vec<Complex64> {
                    (j0..j1).map(|j| entry_fn(i0 + li, j)).collect()
                };
                let col_fn = |lj: usize| -> Vec<Complex64> {
                    (i0..i1).map(|i| entry_fn(i, j0 + lj)).collect()
                };
                let aca = AcaMatrix::assemble(m_block, n_block, &row_fn, &col_fn, tol, max_rank);
                far.push((i0, j0, aca));
            }
        }
    }

    (near, far)
}

/// Apply the ACA-compressed matrix to a vector x → y = Z x.
pub fn aca_matvec(
    n: usize,
    near: &[(usize, usize, Complex64)],
    far: &[(usize, usize, AcaMatrix)],
    x: &[Complex64],
) -> Vec<Complex64> {
    let mut y = vec![Complex64::ZERO; n];

    for &(i, j, z_ij) in near {
        y[i] += z_ij * x[j];
    }

    for (i0, j0, aca) in far {
        let n_block = aca.ncols;
        let m_block = aca.nrows;
        let x_block = DVector::from_iterator(n_block, (0..n_block).map(|k| x[j0 + k]));
        let y_block = aca.matvec(&x_block);
        for i in 0..m_block {
            y[i0 + i] += y_block[i];
        }
    }

    y
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DMatrix;

    /// Rank-2 complex symmetric matrix: Z = u1 v1^T + u2 v2^T
    fn rank2_matrix(m: usize, n: usize) -> DMatrix<Complex64> {
        let u1: Vec<Complex64> = (0..m).map(|i| Complex64::new(i as f64 + 1.0, 0.0)).collect();
        let v1: Vec<Complex64> = (0..n).map(|j| Complex64::new(1.0, j as f64 * 0.1)).collect();
        let u2: Vec<Complex64> = (0..m).map(|i| Complex64::new(0.0, i as f64 * 0.5 + 0.1)).collect();
        let v2: Vec<Complex64> = (0..n).map(|j| Complex64::new((j as f64 + 1.0).recip(), 0.0)).collect();
        let mut z = DMatrix::<Complex64>::zeros(m, n);
        for i in 0..m {
            for j in 0..n {
                z[(i, j)] = u1[i] * v1[j] + u2[i] * v2[j];  // V^T (no conj)
            }
        }
        z
    }

    #[test]
    fn aca_compress_rank2_exact() {
        let z = rank2_matrix(10, 8);
        let aca = AcaMatrix::compress(&z, 1e-10, 20);
        let z_approx = aca.to_dense();

        let mut err = 0.0_f64;
        let mut z_frob = 0.0_f64;
        for i in 0..z.nrows() {
            for j in 0..z.ncols() {
                err    += (z[(i,j)] - z_approx[(i,j)]).norm_sqr();
                z_frob += z[(i,j)].norm_sqr();
            }
        }
        let err_rel = (err / z_frob).sqrt();
        assert!(err_rel < 1e-6, "relative error = {err_rel:.3e}, rank = {}", aca.rank);
        assert!(aca.rank <= 4, "rank should be ≤ 4 for rank-2 matrix, got {}", aca.rank);
    }

    #[test]
    fn aca_matvec_matches_dense() {
        let z = rank2_matrix(12, 10);
        let aca = AcaMatrix::compress(&z, 1e-10, 20);

        let x: Vec<Complex64> = (0..10).map(|i| Complex64::new(i as f64 * 0.1, 1.0)).collect();
        let x_dv = DVector::from_vec(x.clone());

        let y_dense = &z * &x_dv;
        let y_aca   = aca.matvec(&x_dv);

        for i in 0..12 {
            let diff = (y_dense[i] - y_aca[i]).norm();
            assert!(diff < 1e-8, "y[{i}]: dense={:.4e}, aca={:.4e}, diff={:.4e}",
                y_dense[i].norm(), y_aca[i].norm(), diff);
        }
    }

    #[test]
    fn aca_partition_matvec() {
        let n = 20_usize;
        let z = DMatrix::<Complex64>::from_fn(n, n, |i, j| {
            if i == j {
                Complex64::new(10.0, 0.0)
            } else {
                Complex64::new(1.0 / ((i as f64 - j as f64).abs() + 1.0), 0.0)
            }
        });

        let entry_fn = |i: usize, j: usize| z[(i, j)];
        let (near, far) = aca_partition(n, 5, 0, 1e-6, 10, &entry_fn);

        let x: Vec<Complex64> = (0..n).map(|i| Complex64::new(1.0 + i as f64 * 0.1, 0.0)).collect();
        let x_dv = DVector::from_vec(x.clone());

        let y_dense = &z * &x_dv;
        let y_aca   = aca_matvec(n, &near, &far, &x);

        let mut max_err = 0.0_f64;
        for i in 0..n {
            let diff = (y_dense[i] - y_aca[i]).norm();
            if diff > max_err { max_err = diff; }
        }
        assert!(max_err < 1e-4, "max matvec error = {max_err:.3e}");
    }

    #[test]
    fn aca_zero_matrix_rank_zero() {
        let z = DMatrix::<Complex64>::zeros(8, 8);
        let aca = AcaMatrix::compress(&z, 1e-10, 10);
        assert!(aca.rank <= 2, "rank should be near 0 for zero matrix, got {}", aca.rank);
    }
}
