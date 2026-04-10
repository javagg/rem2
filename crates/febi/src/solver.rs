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

    log::info!("FE-BI solve: {} DOFs", n);
    Ok(sol)
}
