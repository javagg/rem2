/// COO → CSR sparse matrix and PCG solver for FEM.
///
/// Design goals:
/// - WASM-compatible (no unsafe, no FFI)
/// - Correct for symmetric positive-definite systems arising from FEM
/// - Interoperable with `fem_linalg` via `to_fem_csr()` / `from_fem_csr()`
///   so physics crates can pass assembled matrices to fem-rs solvers.

// ---------------------------------------------------------------------------
// fem-linalg interop
// ---------------------------------------------------------------------------

use fem_linalg::{CooMatrix as FemCoo, CsrMatrix as FemCsr};


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

    /// Convert to a `fem_linalg::CooMatrix<f64>` for use with fem-rs assemblers.
    ///
    /// The returned matrix owns fresh storage; `self` is unchanged.
    pub fn to_fem_coo(&self) -> FemCoo<f64> {
        let mut out = FemCoo::new(self.nrows, self.ncols);
        out.reserve(self.rows.len());
        for i in 0..self.rows.len() {
            out.add(self.rows[i], self.cols[i], self.vals[i]);
        }
        out
    }

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

    // -----------------------------------------------------------------------
    // fem-linalg interop
    // -----------------------------------------------------------------------

    /// Convert to `fem_linalg::CsrMatrix<f64>` for use with fem-rs solvers.
    ///
    /// `col_idx` is widened from `usize` to `u32`.  Panics in debug builds if
    /// any column index overflows `u32` (> 4 billion columns is not realistic).
    pub fn to_fem_csr(&self) -> FemCsr<f64> {
        let col_idx_u32: Vec<u32> = self.col_idx.iter().map(|&c| {
            debug_assert!(c <= u32::MAX as usize, "column index {} overflows u32", c);
            c as u32
        }).collect();
        FemCsr {
            nrows:   self.nrows,
            ncols:   self.ncols,
            row_ptr: self.row_ptr.clone(),
            col_idx: col_idx_u32,
            values:  self.values.clone(),
        }
    }

    /// Construct from a `fem_linalg::CsrMatrix<f64>`, narrowing `u32` → `usize`.
    pub fn from_fem_csr(src: FemCsr<f64>) -> Self {
        let col_idx: Vec<usize> = src.col_idx.iter().map(|&c| c as usize).collect();
        CsrMatrix {
            nrows:   src.nrows,
            ncols:   src.ncols,
            row_ptr: src.row_ptr,
            col_idx,
            values:  src.values,
        }
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

/// Result from complex linear system solver.
///
/// Stores solution vector (Complex64), iteration count, residual norm, and convergence flag.
#[derive(Debug, Clone)]
pub struct ComplexSolveResult {
    pub solution: Vec<Complex64>,
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
// fem-solver backends
// ---------------------------------------------------------------------------

/// Solve A x = b using fem-rs ILU(0)-preconditioned CG.
///
/// ILU(0) is a stronger preconditioner than SSOR for large, sparse FEM
/// stiffness matrices.  Falls back to returning `Err` on failure so the
/// caller can retry with the built-in SSOR-PCG.
///
/// - `tol`      — relative residual tolerance
/// - `max_iter` — maximum iterations
///
/// Not available on `wasm32` targets (fem-solver links against linger which
/// requires native threading).
#[cfg(not(target_arch = "wasm32"))]
pub fn solve_pcg_ilu0(
    mat: &CsrMatrix,
    b: &[f64],
    tol: f64,
    max_iter: usize,
) -> Result<SolveResult, String> {
    let fem_mat = mat.to_fem_csr();
    let cfg = fem_solver::SolverConfig {
        rtol: tol,
        max_iter,
        ..fem_solver::SolverConfig::default()
    };
    let mut x = vec![0.0f64; b.len()];
    match fem_solver::solve_pcg_ilu0(&fem_mat, b, &mut x, &cfg) {
        Ok(r) => Ok(SolveResult {
            solution: x,
            iterations: r.iterations,
            residual_norm: r.final_residual,
            converged: r.converged,
        }),
        Err(e) => Err(e.to_string()),
    }
}

/// Solve A x = b using AMG-preconditioned CG (fem-rs AMG backend).
///
/// Algebraic Multigrid typically gives 5-30× speedup over Jacobi/SSOR
/// for large SPD systems from FEM discretisations.
#[cfg(not(target_arch = "wasm32"))]
pub fn solve_pcg_amg(
    mat: &CsrMatrix,
    b: &[f64],
    tol: f64,
    max_iter: usize,
) -> Result<SolveResult, String> {
    let fem_mat = mat.to_fem_csr();
    let solver_cfg = fem_solver::SolverConfig {
        rtol: tol,
        max_iter,
        ..fem_solver::SolverConfig::default()
    };
    let amg_cfg = fem_amg::AmgConfig::default();
    let mut x = vec![0.0f64; b.len()];
    match fem_amg::solve_amg_cg(&fem_mat, b, &mut x, &amg_cfg, &solver_cfg) {
        Ok(r) => Ok(SolveResult {
            solution: x,
            iterations: r.iterations,
            residual_norm: r.final_residual,
            converged: r.converged,
        }),
        Err(e) => Err(e.to_string()),
    }
}

/// Solve A x = b using a sparse direct Cholesky factorisation (fem-rs backend).
///
/// Suitable for small to medium symmetric positive-definite systems where
/// iterative methods converge slowly (ill-conditioned problems, eigenmode
/// shift-invert, etc.).
///
/// Not available on `wasm32` targets.
#[cfg(not(target_arch = "wasm32"))]
pub fn solve_cholesky(mat: &CsrMatrix, b: &[f64]) -> Result<SolveResult, String> {
    let fem_mat = mat.to_fem_csr();
    match fem_solver::solve_sparse_cholesky(&fem_mat, b) {
        Ok(x) => Ok(SolveResult {
            solution: x,
            iterations: 0,
            residual_norm: 0.0,
            converged: true,
        }),
        Err(e) => Err(e.to_string()),
    }
}

/// Solve A x = b using fem-rs Conjugate Gradient with a matrix-free operator.
///
/// The `apply` callback computes `y <- A * x`.
#[cfg(not(target_arch = "wasm32"))]
pub fn solve_cg_operator<F>(
    nrows: usize,
    ncols: usize,
    apply: F,
    b: &[f64],
    tol: f64,
    max_iter: usize,
) -> Result<SolveResult, String>
where
    F: Fn(&[f64], &mut [f64]),
{
    let cfg = fem_solver::SolverConfig {
        rtol: tol,
        max_iter,
        ..fem_solver::SolverConfig::default()
    };
    let mut x = vec![0.0f64; ncols];
    match fem_solver::solve_cg_operator(nrows, ncols, apply, b, &mut x, &cfg) {
        Ok(r) => Ok(SolveResult {
            solution: x,
            iterations: r.iterations,
            residual_norm: r.final_residual,
            converged: r.converged,
        }),
        Err(e) => Err(e.to_string()),
    }
}

/// Solve A x = b using fem-rs restarted GMRES with a matrix-free operator.
///
/// The `apply` callback computes `y <- A * x`.
#[cfg(not(target_arch = "wasm32"))]
pub fn solve_gmres_operator<F>(
    nrows: usize,
    ncols: usize,
    apply: F,
    b: &[f64],
    restart: usize,
    tol: f64,
    max_iter: usize,
) -> Result<SolveResult, String>
where
    F: Fn(&[f64], &mut [f64]),
{
    let cfg = fem_solver::SolverConfig {
        rtol: tol,
        max_iter,
        ..fem_solver::SolverConfig::default()
    };
    let mut x = vec![0.0f64; ncols];
    match fem_solver::solve_gmres_operator(nrows, ncols, apply, b, &mut x, restart, &cfg) {
        Ok(r) => Ok(SolveResult {
            solution: x,
            iterations: r.iterations,
            residual_norm: r.final_residual,
            converged: r.converged,
        }),
        Err(e) => Err(e.to_string()),
    }
}

/// Solve A x = b using fem-rs BiCGSTAB with a matrix-free operator.
///
/// The `apply` callback computes `y <- A * x`.
#[cfg(not(target_arch = "wasm32"))]
pub fn solve_bicgstab_operator<F>(
    nrows: usize,
    ncols: usize,
    apply: F,
    b: &[f64],
    tol: f64,
    max_iter: usize,
) -> Result<SolveResult, String>
where
    F: Fn(&[f64], &mut [f64]),
{
    let cfg = fem_solver::SolverConfig {
        rtol: tol,
        max_iter,
        ..fem_solver::SolverConfig::default()
    };
    let mut x = vec![0.0f64; ncols];
    match fem_solver::solve_bicgstab_operator(nrows, ncols, apply, b, &mut x, &cfg) {
        Ok(r) => Ok(SolveResult {
            solution: x,
            iterations: r.iterations,
            residual_norm: r.final_residual,
            converged: r.converged,
        }),
        Err(e) => Err(e.to_string()),
    }
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

    #[cfg(not(target_arch = "wasm32"))]
    fn laplacian_apply(x: &[f64], y: &mut [f64]) {
        assert_eq!(x.len(), 3);
        assert_eq!(y.len(), 3);
        y[0] = 2.0 * x[0] - x[1];
        y[1] = -x[0] + 2.0 * x[1] - x[2];
        y[2] = -x[1] + 2.0 * x[2];
    }

    #[cfg(not(target_arch = "wasm32"))]
    #[test]
    fn cg_operator_solves_laplacian() {
        let b = vec![0.0, 1.0, 0.0];
        let res = solve_cg_operator(3, 3, laplacian_apply, &b, 1e-12, 100)
            .expect("CG operator solve should succeed");
        assert!(res.converged);
        assert!((res.solution[0] - 0.5).abs() < 1e-10);
        assert!((res.solution[1] - 1.0).abs() < 1e-10);
        assert!((res.solution[2] - 0.5).abs() < 1e-10);
    }

    #[cfg(not(target_arch = "wasm32"))]
    #[test]
    fn gmres_operator_solves_laplacian() {
        let b = vec![0.0, 1.0, 0.0];
        let res = solve_gmres_operator(3, 3, laplacian_apply, &b, 8, 1e-12, 100)
            .expect("GMRES operator solve should succeed");
        assert!(res.converged);
        assert!((res.solution[0] - 0.5).abs() < 1e-10);
        assert!((res.solution[1] - 1.0).abs() < 1e-10);
        assert!((res.solution[2] - 0.5).abs() < 1e-10);
    }

    #[cfg(not(target_arch = "wasm32"))]
    #[test]
    fn bicgstab_operator_solves_laplacian() {
        let b = vec![0.0, 1.0, 0.0];
        let res = solve_bicgstab_operator(3, 3, laplacian_apply, &b, 1e-12, 100)
            .expect("BiCGSTAB operator solve should succeed");
        assert!(res.converged);
        assert!((res.solution[0] - 0.5).abs() < 1e-10);
        assert!((res.solution[1] - 1.0).abs() < 1e-10);
        assert!((res.solution[2] - 0.5).abs() < 1e-10);
    }
}

// ---------------------------------------------------------------------------
// CsrMatrixComplex: Complex sparse matrix with LinearOperator support
// ---------------------------------------------------------------------------

use num_complex::Complex64;

/// Sparse matrix in CSR format with Complex64 values.
/// Supports LinearOperator interface for use with generic GMRES and other iterative solvers.
#[derive(Debug, Clone)]
pub struct CsrMatrixComplex {
    pub nrows: usize,
    pub ncols: usize,
    pub row_ptr: Vec<usize>,
    pub col_idx: Vec<usize>,
    pub values: Vec<Complex64>,
}

impl CsrMatrixComplex {
    /// Construct a zero matrix.
    pub fn new(nrows: usize, ncols: usize) -> Self {
        CsrMatrixComplex {
            nrows,
            ncols,
            row_ptr: vec![0; nrows + 1],
            col_idx: Vec::new(),
            values: Vec::new(),
        }
    }

    /// Number of non-zero entries.
    pub fn nnz(&self) -> usize {
        self.values.len()
    }

    /// Convert a dense DMatrix<Complex64> to CSR format (skipping zero entries).
    pub fn from_dense(mat: &nalgebra::DMatrix<Complex64>) -> Self {
        let (nrows, ncols) = mat.shape();
        let mut row_ptr = Vec::with_capacity(nrows + 1);
        let mut col_idx = Vec::new();
        let mut values = Vec::new();

        row_ptr.push(0);
        for i in 0..nrows {
            for j in 0..ncols {
                if mat[(i, j)].norm() > 1e-16 {
                    col_idx.push(j);
                    values.push(mat[(i, j)]);
                }
            }
            row_ptr.push(values.len());
        }

        CsrMatrixComplex {
            nrows,
            ncols,
            row_ptr,
            col_idx,
            values,
        }
    }

    /// Sparse matrix–vector product: y = A * x
    pub fn matvec(&self, x: &nalgebra::DVector<Complex64>, y: &mut nalgebra::DVector<Complex64>) -> Result<(), String> {
        if x.len() != self.ncols || y.len() != self.nrows {
            return Err(format!(
                "CsrMatrixComplex::matvec dimension mismatch: matrix {}×{}, x.len()={}, y.len()={}",
                self.nrows, self.ncols, x.len(), y.len()
            ));
        }
        
        for i in 0..self.nrows {
            let mut s = Complex64::new(0.0, 0.0);
            for k in self.row_ptr[i]..self.row_ptr[i + 1] {
                let j = self.col_idx[k];
                s += self.values[k] * x[j];
            }
            y[i] = s;
        }
        Ok(())
    }

    /// Extract the main diagonal.
    pub fn diagonal(&self) -> nalgebra::DVector<Complex64> {
        let len = self.nrows.min(self.ncols);
        let mut d = nalgebra::DVector::zeros(len);
        for i in 0..len {
            for k in self.row_ptr[i]..self.row_ptr[i + 1] {
                if self.col_idx[k] == i {
                    d[i] = self.values[k];
                    break;
                }
            }
        }
        d
    }
}

/// Implement LinearOperator trait for CsrMatrixComplex.
impl crate::operator::LinearOperator<Complex64> for CsrMatrixComplex {
    fn size(&self) -> (usize, usize) {
        (self.nrows, self.ncols)
    }

    fn matvec(&self, x: &nalgebra::DVector<Complex64>, y: &mut nalgebra::DVector<Complex64>) -> Result<(), String> {
        CsrMatrixComplex::matvec(self, x, y)
    }

    fn matvec_adjoint(&self, x: &nalgebra::DVector<Complex64>, y: &mut nalgebra::DVector<Complex64>) -> Result<(), String> {
        if x.len() != self.nrows || y.len() != self.ncols {
            return Err(format!(
                "CsrMatrixComplex::matvec_adjoint dimension mismatch: matrix {}×{}, x.len()={}, y.len()={}",
                self.nrows, self.ncols, x.len(), y.len()
            ));
        }
        
        // Clear y first
        for i in 0..y.len() {
            y[i] = Complex64::new(0.0, 0.0);
        }
        
        // Accumulate: y[j] += conj(A[i,j]) * x[i]
        for i in 0..self.nrows {
            for k in self.row_ptr[i]..self.row_ptr[i + 1] {
                let j = self.col_idx[k];
                y[j] += self.values[k].conj() * x[i];
            }
        }
        Ok(())
    }

    fn diagonal(&self) -> Option<nalgebra::DVector<Complex64>> {
        Some(CsrMatrixComplex::diagonal(self))
    }

    fn density(&self) -> f64 {
        let total = self.nrows * self.ncols;
        if total == 0 {
            0.0
        } else {
            self.nnz() as f64 / total as f64
        }
    }
}

// ---------------------------------------------------------------------------
// Complex iterative solver for non-Hermitian Helmholtz systems
// ---------------------------------------------------------------------------

/// Solve complex linear system using right-preconditioned BiCGSTAB.
///
/// Despite the historical name (`solve_pcg_complex`), this routine is intended
/// for non-Hermitian systems arising from frequency-domain Helmholtz FEM.
///
/// Preconditioner: Jacobi (inverse diagonal, with safe fallback for near-zero diagonal).
pub fn solve_pcg_complex(
    mat: &CsrMatrixComplex,
    b: &[Complex64],
    tol: f64,
    max_iter: usize,
) -> ComplexSolveResult {
    let n = b.len();
    if n == 0 {
        return ComplexSolveResult {
            solution: vec![],
            iterations: 0,
            residual_norm: 0.0,
            converged: true,
        };
    }

    if mat.nrows != n || mat.ncols != n {
        return ComplexSolveResult {
            solution: vec![],
            iterations: 0,
            residual_norm: f64::NAN,
            converged: false,
        };
    }

    let b_vec = nalgebra::DVector::from_row_slice(b);
    let b_norm = b_vec.norm();
    
    if b_norm < f64::EPSILON {
        return ComplexSolveResult {
            solution: vec![Complex64::ZERO; n],
            iterations: 0,
            residual_norm: 0.0,
            converged: true,
        };
    }

    // Jacobi preconditioner M^{-1} = diag(A)^{-1}
    let diag = mat.diagonal();
    let mut minv = nalgebra::DVector::<Complex64>::zeros(n);
    for i in 0..n {
        if diag[i].norm() > 1e-30 {
            minv[i] = Complex64::new(1.0, 0.0) / diag[i];
        } else {
            minv[i] = Complex64::new(1.0, 0.0);
        }
    }

    let mut x = nalgebra::DVector::<Complex64>::zeros(n);
    let mut r = b_vec.clone(); // x = 0 => r = b
    let r_hat = r.clone();

    let mut p = nalgebra::DVector::<Complex64>::zeros(n);
    let mut v = nalgebra::DVector::<Complex64>::zeros(n);
    let mut s = nalgebra::DVector::<Complex64>::zeros(n);
    let mut t = nalgebra::DVector::<Complex64>::zeros(n);
    let mut p_hat = nalgebra::DVector::<Complex64>::zeros(n);
    let mut s_hat = nalgebra::DVector::<Complex64>::zeros(n);

    let mut rho_prev = Complex64::new(1.0, 0.0);
    let mut alpha = Complex64::new(1.0, 0.0);
    let mut omega = Complex64::new(1.0, 0.0);

    let mut iterations = 0;
    let mut converged = false;
    let mut residual_norm = r.norm() / b_norm;
    if residual_norm <= tol {
        return ComplexSolveResult {
            solution: x.iter().copied().collect(),
            iterations,
            residual_norm,
            converged: true,
        };
    }

    for iter in 0..max_iter {
        iterations = iter + 1;

        let rho = r_hat.dotc(&r);
        if rho.norm() < 1e-30 {
            break;
        }

        let beta = if iter == 0 {
            Complex64::new(0.0, 0.0)
        } else {
            (rho / rho_prev) * (alpha / omega)
        };

        for i in 0..n {
            p[i] = r[i] + beta * (p[i] - omega * v[i]);
            p_hat[i] = minv[i] * p[i];
        }

        if mat.matvec(&p_hat, &mut v).is_err() {
            break;
        }
        let denom = r_hat.dotc(&v);
        if denom.norm() < 1e-30 {
            break;
        }
        alpha = rho / denom;

        for i in 0..n {
            s[i] = r[i] - alpha * v[i];
        }

        let s_rel = s.norm() / b_norm;
        if s_rel <= tol {
            for i in 0..n {
                x[i] += alpha * p_hat[i];
            }
            residual_norm = s_rel;
            converged = true;
            break;
        }

        for i in 0..n {
            s_hat[i] = minv[i] * s[i];
        }
        if mat.matvec(&s_hat, &mut t).is_err() {
            break;
        }

        let tt = t.dotc(&t);
        if tt.norm() < 1e-30 {
            break;
        }
        omega = t.dotc(&s) / tt;
        if omega.norm() < 1e-30 {
            break;
        }

        for i in 0..n {
            x[i] += alpha * p_hat[i] + omega * s_hat[i];
            r[i] = s[i] - omega * t[i];
        }

        residual_norm = r.norm() / b_norm;
        if residual_norm <= tol {
            converged = true;
            break;
        }

        rho_prev = rho;
    }

    ComplexSolveResult {
        solution: x.iter().copied().collect(),
        iterations,
        residual_norm,
        converged,
    }
}

#[cfg(test)]
mod tests_csr_complex {
    use super::*;

    #[test]
    fn csr_complex_matvec_basic() {
        let mut mat = CsrMatrixComplex::new(2, 2);
        mat.row_ptr = vec![0, 1, 2];
        mat.col_idx = vec![0, 1];
        mat.values = vec![Complex64::new(1.0, 1.0), Complex64::new(2.0, 0.0)];
        
        let x = nalgebra::DVector::from_vec(vec![
            Complex64::new(1.0, 0.0),
            Complex64::new(0.0, 1.0),
        ]);
        let mut y = nalgebra::DVector::zeros(2);
        
        mat.matvec(&x, &mut y).unwrap();
        
        // y[0] = (1+i) * (1+0i) = 1+i
        assert!((y[0].re - 1.0).abs() < 1e-14);
        assert!((y[0].im - 1.0).abs() < 1e-14);
        
        // y[1] = 2 * (0+i) = 0+2i
        assert!((y[1].re - 0.0).abs() < 1e-14);
        assert!((y[1].im - 2.0).abs() < 1e-14);
    }

    #[test]
    fn csr_complex_adjoint_basic() {
        use crate::operator::LinearOperator;

        let mut mat = CsrMatrixComplex::new(2, 2);
        mat.row_ptr = vec![0, 1, 2];
        mat.col_idx = vec![0, 1];
        mat.values = vec![Complex64::new(1.0, 1.0), Complex64::new(2.0, 0.0)];
        
        let x = nalgebra::DVector::from_vec(vec![
            Complex64::new(1.0, 0.0),
            Complex64::new(0.0, 1.0),
        ]);
        let mut y = nalgebra::DVector::zeros(2);
        
        mat.matvec_adjoint(&x, &mut y).unwrap();
        
        // y[0] = conj(1+i) * (1+0i) = (1-i) * 1 = 1-i
        assert!((y[0].re - 1.0).abs() < 1e-14);
        assert!((y[0].im + 1.0).abs() < 1e-14);
        
        // y[1] = conj(2) * (0+i) = 2 * (0+i) = 0+2i
        assert!((y[1].re - 0.0).abs() < 1e-14);
        assert!((y[1].im - 2.0).abs() < 1e-14);
    }

    #[test]
    fn csr_complex_implements_linear_operator() {

        let mut mat = CsrMatrixComplex::new(3, 3);
        mat.row_ptr = vec![0, 1, 2, 3];
        mat.col_idx = vec![0, 1, 2];
        mat.values = vec![
            Complex64::new(1.0, 0.0),
            Complex64::new(2.0, 0.0),
            Complex64::new(3.0, 0.0),
        ];
        
        let x = nalgebra::DVector::from_vec(vec![
            Complex64::new(1.0, 0.0),
            Complex64::new(1.0, 0.0),
            Complex64::new(1.0, 0.0),
        ]);
        let mut y = nalgebra::DVector::zeros(3);
        
        // Via LinearOperator trait
        mat.matvec(&x, &mut y).unwrap();
        
        assert!((y[0].re - 1.0).abs() < 1e-14);
        assert!((y[1].re - 2.0).abs() < 1e-14);
        assert!((y[2].re - 3.0).abs() < 1e-14);
    }

    #[test]
    fn solve_pcg_complex_diagonal() {
        // Test on a simple diagonal system: diag(1, 2, 3) * x = (1, 2, 3)
        // Solution: x = (1, 1, 1)
        let mut mat = CsrMatrixComplex::new(3, 3);
        mat.row_ptr = vec![0, 1, 2, 3];
        mat.col_idx = vec![0, 1, 2];
        mat.values = vec![
            Complex64::new(1.0, 0.0),
            Complex64::new(2.0, 0.0),
            Complex64::new(3.0, 0.0),
        ];
        
        let b = vec![
            Complex64::new(1.0, 0.0),
            Complex64::new(2.0, 0.0),
            Complex64::new(3.0, 0.0),
        ];
        
        let result = solve_pcg_complex(&mat, &b, 1e-10, 100);
        
        assert!(result.converged, "Solver should converge");
        assert_eq!(result.solution.len(), 3);
        assert!((result.solution[0].re - 1.0).abs() < 1e-8);
        assert!((result.solution[1].re - 1.0).abs() < 1e-8);
        assert!((result.solution[2].re - 1.0).abs() < 1e-8);
    }

    #[test]
    fn solve_pcg_complex_helmholtz_2x2() {
        // Minimal Helmholtz-like system: -Δu - k²u = f on unit interval
        // Discretized 2x2 system with complex wavenumber
        let mut mat = CsrMatrixComplex::new(2, 2);
        // Row 0: [2+i, -1]
        // Row 1: [-1, 2+i]
        mat.row_ptr = vec![0, 2, 4];
        mat.col_idx = vec![0, 1, 0, 1];
        mat.values = vec![
            Complex64::new(2.0, 1.0),
            Complex64::new(-1.0, 0.0),
            Complex64::new(-1.0, 0.0),
            Complex64::new(2.0, 1.0),
        ];
        
        let b = vec![
            Complex64::new(1.0, 0.0),
            Complex64::new(1.0, 0.0),
        ];
        
        let result = solve_pcg_complex(&mat, &b, 1e-8, 100);
        
        assert!(result.converged, "Solver should converge for Helmholtz system");
        assert_eq!(result.solution.len(), 2);
        assert!(result.residual_norm < 1e-8);
    }
}
