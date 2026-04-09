/// Minimal COO → CSR sparse matrix and PCG solver for FEM.
///
/// Design goals:
/// - No external dependencies (WASM-compatible)
/// - Simple enough to audit easily
/// - Correct for symmetric positive-definite systems arising from FEM

// ---------------------------------------------------------------------------
// Triplet (COO) format — used during assembly
// ---------------------------------------------------------------------------

/// Sparse matrix in triplet (COO) format.
/// Duplicate entries are allowed and will be summed during `to_csr()`.
#[derive(Debug, Clone)]
pub struct TripletMatrix {
    pub nrows: usize,
    pub ncols: usize,
    rows: Vec<usize>,
    cols: Vec<usize>,
    vals: Vec<f64>,
}

impl TripletMatrix {
    pub fn new(nrows: usize, ncols: usize) -> Self {
        TripletMatrix { nrows, ncols, rows: Vec::new(), cols: Vec::new(), vals: Vec::new() }
    }

    pub fn with_capacity(nrows: usize, ncols: usize, cap: usize) -> Self {
        TripletMatrix {
            nrows,
            ncols,
            rows: Vec::with_capacity(cap),
            cols: Vec::with_capacity(cap),
            vals: Vec::with_capacity(cap),
        }
    }

    #[inline]
    pub fn add(&mut self, row: usize, col: usize, val: f64) {
        debug_assert!(row < self.nrows, "row {} >= nrows {}", row, self.nrows);
        debug_assert!(col < self.ncols, "col {} >= ncols {}", col, self.ncols);
        self.rows.push(row);
        self.cols.push(col);
        self.vals.push(val);
    }

    pub fn nnz(&self) -> usize { self.rows.len() }

    /// Remap node indices for periodic boundary conditions.
    ///
    /// For each (donor, receiver) pair, replaces all occurrences of `receiver`
    /// with `donor` in the row and column index arrays.  When `to_csr()` is
    /// subsequently called, the receiver contributions accumulate into the donor
    /// DOF, implementing the Γ-point periodic constraint φ[recv] = φ[donor].
    ///
    /// The `nrows` / `ncols` are NOT reduced — the receiver DOFs still exist
    /// in the assembled matrix but with only a diagonal = 1 placeholder added
    /// by `apply_dirichlet` when `recv` is put into the Dirichlet map with value 0.
    pub fn remap_periodic_nodes(&mut self, pairs: &[(usize, usize)]) {
        if pairs.is_empty() {
            return;
        }
        // Build a substitution map: recv → donor (follow chains)
        let n = self.nrows;
        let mut subst: Vec<usize> = (0..n).collect();
        for &(donor, recv) in pairs {
            subst[recv] = donor;
        }
        // Apply substitution (one pass is enough since donor is never a receiver in well-formed input)
        for r in &mut self.rows {
            *r = subst[*r];
        }
        for c in &mut self.cols {
            *c = subst[*c];
        }
    }

    /// Convert to CSR format, summing duplicate (row, col) entries.
    pub fn to_csr(self) -> CsrMatrix {
        let n = self.rows.len();
        if n == 0 {
            return CsrMatrix {
                nrows: self.nrows,
                ncols: self.ncols,
                row_ptr: vec![0; self.nrows + 1],
                col_idx: Vec::new(),
                values: Vec::new(),
            };
        }

        // Sort indices by (row, col)
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_unstable_by(|&a, &b| {
            self.rows[a]
                .cmp(&self.rows[b])
                .then(self.cols[a].cmp(&self.cols[b]))
        });

        let mut col_idx: Vec<usize> = Vec::with_capacity(n);
        let mut values: Vec<f64> = Vec::with_capacity(n);
        let mut row_ptr = vec![0usize; self.nrows + 1];

        let first = order[0];
        let mut cur_row = self.rows[first];
        let mut cur_col = self.cols[first];
        let mut cur_val = self.vals[first];

        for k in 1..n {
            let i = order[k];
            let r = self.rows[i];
            let c = self.cols[i];
            let v = self.vals[i];

            if r == cur_row && c == cur_col {
                cur_val += v; // merge duplicate
            } else {
                // commit current entry
                col_idx.push(cur_col);
                values.push(cur_val);

                if r != cur_row {
                    // mark end of cur_row
                    row_ptr[cur_row + 1] = col_idx.len();
                    cur_row = r;
                }
                cur_col = c;
                cur_val = v;
            }
        }
        // commit final entry
        col_idx.push(cur_col);
        values.push(cur_val);
        row_ptr[cur_row + 1] = col_idx.len();

        // forward-fill row_ptr for empty rows
        for i in 1..=self.nrows {
            if row_ptr[i] == 0 {
                row_ptr[i] = row_ptr[i - 1];
            }
        }

        CsrMatrix { nrows: self.nrows, ncols: self.ncols, row_ptr, col_idx, values }
    }
}

// ---------------------------------------------------------------------------
// CSR sparse matrix
// ---------------------------------------------------------------------------

/// Sparse matrix in Compressed Sparse Row (CSR) format.
#[derive(Debug, Clone)]
pub struct CsrMatrix {
    pub nrows: usize,
    pub ncols: usize,
    /// Row pointers: row i spans `col_idx[row_ptr[i]..row_ptr[i+1]]`.
    pub row_ptr: Vec<usize>,
    /// Column indices (sorted within each row).
    pub col_idx: Vec<usize>,
    /// Non-zero values.
    pub values: Vec<f64>,
}

impl CsrMatrix {
    pub fn nnz(&self) -> usize { self.values.len() }

    /// Sparse matrix–vector product: y = A * x
    pub fn matvec(&self, x: &[f64], y: &mut [f64], comm: &dyn rem_parallel::Comm) {
        assert_eq!(x.len(), self.ncols);
        assert_eq!(y.len(), self.nrows);
        for i in 0..self.nrows {
            let mut s = 0.0;
            for k in self.row_ptr[i]..self.row_ptr[i + 1] {
                s += self.values[k] * x[self.col_idx[k]];
            }
            y[i] = s;
        }
        if comm.size() > 1 {
            comm.allreduce_f64_vec(y);
        }
    }

    /// Extract the main diagonal.
    pub fn diagonal(&self) -> Vec<f64> {
        let mut d = vec![0.0; self.nrows.min(self.ncols)];
        for i in 0..d.len() {
            for k in self.row_ptr[i]..self.row_ptr[i + 1] {
                if self.col_idx[k] == i {
                    d[i] = self.values[k];
                    break;
                }
            }
        }
        d
    }

    /// Set a single entry (row, col). Panics if the entry does not exist in the sparsity pattern.
    pub fn set(&mut self, row: usize, col: usize, val: f64) {
        for k in self.row_ptr[row]..self.row_ptr[row + 1] {
            if self.col_idx[k] == col {
                self.values[k] = val;
                return;
            }
        }
        panic!("entry ({},{}) not found in sparsity pattern", row, col);
    }

    /// Zero out an entire row and set the diagonal to `diag_val`.
    pub fn zero_row_set_diag(&mut self, row: usize, diag_val: f64) {
        for k in self.row_ptr[row]..self.row_ptr[row + 1] {
            if self.col_idx[k] == row {
                self.values[k] = diag_val;
            } else {
                self.values[k] = 0.0;
            }
        }
        // Note: if the row has no diagonal entry (isolated DOF not in any element),
        // we simply zero the row. PCG will see a 0-diagonal and skip the Jacobi update.
    }

    /// Read the diagonal entry K[row,row] without modifying the matrix.
    pub fn diagonal_entry(&self, row: usize) -> f64 {
        for k in self.row_ptr[row]..self.row_ptr[row + 1] {
            if self.col_idx[k] == row {
                return self.values[k];
            }
        }
        0.0
    }

    /// For each row `i`, zero the entry at column `dof` and return the original value.
    /// Used during symmetric Dirichlet elimination to modify the RHS.
    pub fn zero_col_entry(&mut self, row: usize, col: usize) -> f64 {
        for k in self.row_ptr[row]..self.row_ptr[row + 1] {
            if self.col_idx[k] == col {
                let v = self.values[k];
                self.values[k] = 0.0;
                return v;
            }
        }
        0.0
    }
}

// ---------------------------------------------------------------------------
// Preconditioned Conjugate Gradient solver
// ---------------------------------------------------------------------------

/// Result of a PCG solve.
#[derive(Debug)]
pub struct SolveResult {
    pub solution: Vec<f64>,
    pub iterations: usize,
    pub residual_norm: f64,
    pub converged: bool,
}

/// Solve A x = b using Preconditioned Conjugate Gradient with SSOR
/// (Symmetric Successive Over-Relaxation, ω = 1.5) preconditioner.
///
/// SSOR is significantly more effective than Jacobi for FEM stiffness matrices,
/// typically reducing iteration counts by 3–5×.
///
/// - `tol`: convergence tolerance (relative: ||r|| / ||b|| < tol)
/// - `max_iter`: maximum number of CG iterations
use rem_parallel::Comm;

pub fn solve_pcg(mat: &CsrMatrix, b: &[f64], tol: f64, max_iter: usize, comm: &dyn Comm) -> SolveResult {
    let n = b.len();
    assert_eq!(mat.nrows, n);
    assert_eq!(mat.ncols, n);

    let diag = mat.diagonal();

    let b_norm = comm.allreduce_f64(dot(b, b)).sqrt();
    if b_norm < 1e-300 {
        return SolveResult {
            solution: vec![0.0; n],
            iterations: 0,
            residual_norm: 0.0,
            converged: true,
        };
    }

    let mut x = vec![0.0f64; n];
    let mut r = b.to_vec(); // r = b - A*x = b (x=0)
    let mut z = apply_ssor(mat, &diag, &r, 1.5);
    let mut p = z.clone();
    let mut rz = comm.allreduce_f64(dot(&r, &z));

    let mut ap = vec![0.0f64; n];

    for iter in 0..max_iter {
        mat.matvec(&p, &mut ap, comm);
        let pap = comm.allreduce_f64(dot(&p, &ap));
        if pap.abs() < 1e-300 {
            break;
        }
        let alpha = rz / pap;

        // x = x + alpha * p
        axpy(alpha, &p, &mut x);
        // r = r - alpha * A*p
        axpy(-alpha, &ap, &mut r);

        let r_norm = comm.allreduce_f64(dot(&r, &r)).sqrt();
        if r_norm < tol * b_norm {
            return SolveResult {
                solution: x,
                iterations: iter + 1,
                residual_norm: r_norm,
                converged: true,
            };
        }

        z = apply_ssor(mat, &diag, &r, 1.5);
        let rz_new = comm.allreduce_f64(dot(&r, &z));
        let beta = rz_new / rz;
        // p = z + beta * p
        for i in 0..n {
            p[i] = z[i] + beta * p[i];
        }
        rz = rz_new;
    }

    let r_norm = comm.allreduce_f64(dot(&r, &r)).sqrt();
    SolveResult {
        solution: x,
        iterations: max_iter,
        residual_norm: r_norm,
        converged: r_norm < tol * b_norm,
    }
}

// ---------------------------------------------------------------------------
// BLAS-like helpers (no allocation)
// ---------------------------------------------------------------------------

#[inline]
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(&x, &y)| x * y).sum()
}

/// y += alpha * x
#[inline]
fn axpy(alpha: f64, x: &[f64], y: &mut [f64]) {
    for (yi, &xi) in y.iter_mut().zip(x.iter()) {
        *yi += alpha * xi;
    }
}

/// SSOR preconditioner: M^{-1} r ≈ (D/ω + L)^{-1} · D/(2-ω) · (D/ω + U)^{-1} r
/// where L, D, U are the strict lower, diagonal, strict upper parts of A.
/// Implemented as: forward SOR sweep then backward SOR sweep.
fn apply_ssor(mat: &CsrMatrix, diag: &[f64], r: &[f64], omega: f64) -> Vec<f64> {
    let n = mat.nrows;
    let mut z = vec![0.0f64; n];

    // Forward sweep: (D/ω + L) z = r
    for i in 0..n {
        let mut s = r[i];
        for k in mat.row_ptr[i]..mat.row_ptr[i + 1] {
            let j = mat.col_idx[k];
            if j < i {
                s -= mat.values[k] * z[j];
            }
        }
        let d = diag[i];
        z[i] = if d.abs() > 1e-300 { s * omega / d } else { 0.0 };
    }

    // Scale by D * (2-ω)/ω
    for i in 0..n {
        let d = diag[i];
        z[i] *= d * (2.0 - omega) / omega;
    }

    // Backward sweep: (D/ω + U) z = (scaled z from above)
    for i in (0..n).rev() {
        let mut s = z[i];
        for k in mat.row_ptr[i]..mat.row_ptr[i + 1] {
            let j = mat.col_idx[k];
            if j > i {
                s -= mat.values[k] * z[j];
            }
        }
        let d = diag[i];
        z[i] = if d.abs() > 1e-300 { s * omega / d } else { 0.0 };
    }

    z
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn build_laplacian_1d(n: usize) -> CsrMatrix {
        // 1-D Laplacian: tridiagonal [−1, 2, −1]
        let mut t = TripletMatrix::new(n, n);
        for i in 0..n {
            t.add(i, i, 2.0);
            if i > 0     { t.add(i, i - 1, -1.0); }
            if i + 1 < n { t.add(i, i + 1, -1.0); }
        }
        t.to_csr()
    }

    #[test]
    fn triplet_to_csr_basic() {
        let mut t = TripletMatrix::new(3, 3);
        t.add(0, 0, 1.0);
        t.add(1, 1, 2.0);
        t.add(2, 2, 3.0);
        t.add(0, 1, 0.5);
        t.add(1, 0, 0.5);
        let csr = t.to_csr();
        assert_eq!(csr.nnz(), 5);
        // diagonal
        let d = csr.diagonal();
        assert!((d[0] - 1.0).abs() < 1e-14);
        assert!((d[1] - 2.0).abs() < 1e-14);
        assert!((d[2] - 3.0).abs() < 1e-14);
    }

    #[test]
    fn triplet_sums_duplicates() {
        let mut t = TripletMatrix::new(2, 2);
        t.add(0, 0, 1.0);
        t.add(0, 0, 2.0); // duplicate
        t.add(1, 1, 4.0);
        let csr = t.to_csr();
        let d = csr.diagonal();
        assert!((d[0] - 3.0).abs() < 1e-14);
        assert!((d[1] - 4.0).abs() < 1e-14);
    }

    #[test]
    fn matvec_laplacian() {
        let mat = build_laplacian_1d(4);
        let x = vec![1.0, 2.0, 3.0, 4.0];
        let mut y = vec![0.0; 4];
        mat.matvec(&x, &mut y, &rem_parallel::NoComm);
        // K = [2,-1,0,0; -1,2,-1,0; 0,-1,2,-1; 0,0,-1,2]
        // y[0] = 2*1 - 1*2 = 0
        // y[1] = -1*1 + 2*2 - 1*3 = 0
        // y[2] = -1*2 + 2*3 - 1*4 = 0
        // y[3] = -1*3 + 2*4 = 5
        assert!((y[0] - 0.0).abs() < 1e-14);
        assert!((y[1] - 0.0).abs() < 1e-14);
        assert!((y[2] - 0.0).abs() < 1e-14);
        assert!((y[3] - 5.0).abs() < 1e-14);
    }

    #[test]
    fn pcg_solves_laplacian() {
        // Solve L x = b on 5-node 1D Laplacian with Dirichlet at both ends.
        // K = [[2,-1,0],[-1,2,-1],[0,-1,2]], b = [0, 1, 0] → x = [0.5, 1.0, 0.5]
        let mut t = TripletMatrix::new(3, 3);
        t.add(0, 0, 2.0); t.add(0, 1, -1.0);
        t.add(1, 0, -1.0); t.add(1, 1, 2.0); t.add(1, 2, -1.0);
        t.add(2, 1, -1.0); t.add(2, 2, 2.0);
        let mat = t.to_csr();
        let b = vec![0.0, 1.0, 0.0];
        let res = solve_pcg(&mat, &b, 1e-12, 100, &rem_parallel::NoComm);
        assert!(res.converged, "PCG did not converge");
        assert!((res.solution[0] - 0.5).abs() < 1e-10);
        assert!((res.solution[1] - 1.0).abs() < 1e-10);
        assert!((res.solution[2] - 0.5).abs() < 1e-10);
    }

    #[test]
    fn empty_triplet() {
        let t = TripletMatrix::new(5, 5);
        let csr = t.to_csr();
        assert_eq!(csr.nnz(), 0);
        assert_eq!(csr.row_ptr, vec![0; 6]);
    }

    #[test]
    fn empty_rows_in_csr() {
        // Only row 0 and row 3 have entries
        let mut t = TripletMatrix::new(4, 4);
        t.add(0, 0, 1.0);
        t.add(3, 3, 4.0);
        let csr = t.to_csr();
        assert_eq!(csr.row_ptr, vec![0, 1, 1, 1, 2]);
        assert_eq!(csr.col_idx, vec![0, 3]);
    }
}
