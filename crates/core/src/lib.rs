pub mod constants;
pub mod error;
pub mod memory;
pub mod near_field;
pub mod sparse;
pub mod operator;
pub mod terminal_ui;
pub mod timing;

pub use constants::{EPS0, MU0, C0, ETA0, NU0};
pub use error::{RemError, RemResult};
pub use memory::report_peak_memory;
pub use near_field::{NearFieldPoint, write_near_field_csv, read_near_field_csv, interpolate_e_at, interpolate_e_vec_at};
pub use sparse::{CsrMatrix, CsrMatrixComplex, TripletMatrix, SolveResult, ComplexSolveResult, solve_pcg, solve_pcg_complex, combine_csr_complex};
pub use operator::{LinearOperator, LinearSolver};

#[cfg(not(target_arch = "wasm32"))]
pub use sparse::{solve_bicgstab_operator, solve_cg_operator, solve_gmres_operator};

// ---------------------------------------------------------------------------
// Unified solver entry point
// ---------------------------------------------------------------------------

/// Solve the symmetric positive-definite system A x = b.
///
/// On native targets, attempts ILU(0)-preconditioned CG from fem-rs (typically
/// 3–10× fewer iterations than SSOR-PCG for FEM stiffness matrices).
/// If ILU(0) fails or the target is `wasm32`, falls back to the built-in
/// SSOR-PCG automatically.
pub fn solve_spd(
    mat: &CsrMatrix,
    b: &[f64],
    tol: f64,
    max_iter: usize,
    comm: &dyn rem_parallel::Comm,
) -> SolveResult {
    let _span = timing::span("solve_spd");
    // Native single-rank path: try AMG-PCG first (fastest), fall back to ILU(0)-PCG
    #[cfg(not(target_arch = "wasm32"))]
    if comm.size() == 1 {
        // 1. AMG-PCG (best for large systems)
        match sparse::solve_pcg_amg(mat, b, tol, max_iter) {
            Ok(r) if r.converged => return r,
            Ok(r) => {
                log::debug!(
                    "AMG-PCG did not converge ({} iters, res={:.3e}), retrying with ILU(0)-PCG",
                    r.iterations, r.residual_norm
                );
            }
            Err(e) => {
                log::debug!("AMG-PCG error: {e}, falling back to ILU(0)-PCG");
            }
        }

        // 2. ILU(0)-PCG (reliable for medium systems)
        match sparse::solve_pcg_ilu0(mat, b, tol, max_iter) {
            Ok(r) if r.converged => return r,
            Ok(r) => {
                log::debug!(
                    "ILU(0)-PCG did not converge ({} iters, res={:.3e}), retrying with SSOR-PCG",
                    r.iterations, r.residual_norm
                );
            }
            Err(e) => {
                log::debug!("ILU(0)-PCG error: {e}, falling back to SSOR-PCG");
            }
        }
    }

    // Fallback (WASM, MPI, or all preconditioned solvers failed)
    solve_pcg(mat, b, tol, max_iter, comm)
}

/// Solve the symmetric positive-definite system A x = b using a matrix-free
/// operator callback.
///
/// This is the operator counterpart to [`solve_spd`]. It currently requires a
/// native, single-rank execution context because the backend is routed through
/// fem-rs operator solvers.
#[cfg(not(target_arch = "wasm32"))]
pub fn solve_spd_operator<F>(
    nrows: usize,
    ncols: usize,
    apply: F,
    b: &[f64],
    tol: f64,
    max_iter: usize,
    comm: &dyn rem_parallel::Comm,
) -> Result<SolveResult, String>
where
    F: Fn(&[f64], &mut [f64]),
{
    if comm.size() != 1 {
        return Err("solve_spd_operator currently requires a single-rank communicator".to_string());
    }
    sparse::solve_cg_operator(nrows, ncols, apply, b, tol, max_iter)
}
