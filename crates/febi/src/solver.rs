//! FE-BI 线性系统求解器

use nalgebra::{DMatrix, DVector};
use num_complex::Complex64;
use rem_core::{RemError, RemResult};

/// FE-BI 系统求解：LU 分解（MVP 实现）
pub fn solve_febi(
    mat: &DMatrix<Complex64>,
    rhs: &DVector<Complex64>,
) -> RemResult<DVector<Complex64>> {
    let n = mat.nrows();
    if n == 0 {
        return Ok(DVector::zeros(0));
    }

    let lu = mat.clone().lu();
    let sol = lu.solve(rhs)
        .ok_or_else(|| RemError::Config("FE-BI system is singular".to_string()))?;

    log::info!("FE-BI solve (LU): {} DOFs", n);
    Ok(sol)
}

/// FE-BI 系统求解：GMRES 迭代求解（基于 LinearOperator trait）
///
/// 用于大规模系统，避免 O(N³) LU 分解成本。
pub fn solve_febi_gmres(
    mat: &DMatrix<Complex64>,
    rhs: &DVector<Complex64>,
    tol: f64,
    max_iters: usize,
) -> RemResult<DVector<Complex64>> {
    let n = mat.nrows();
    if n == 0 {
        return Ok(DVector::zeros(0));
    }

    log::info!("FE-BI solve (GMRES): {} DOFs, tol={:.2e}, max_iters={}", n, tol, max_iters);
    
    // Use gmres_solve_generic from rem-mom, which works with LinearOperator
    rem_mom::gmres_solve_generic(mat, rhs, 30, tol, max_iters)
}
