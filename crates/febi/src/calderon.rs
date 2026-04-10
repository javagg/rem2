//! Calderón 算子组装 — FE-BI 中的边界积分（BI）部分
//!
//! 在辐射截断面 Γ 上组装 EFIE 阻抗矩阵 Z_bi，
//! 复用 rem-mom 的格林函数和奇异积分处理。

use nalgebra::DMatrix;
use num_complex::Complex64;
use rem_core::RemResult;
use rem_mom::surface_mesh::SurfaceMesh;

/// 在辐射截断面 Γ 上组装 Calderón（EFIE）矩阵。
///
/// 返回 N×N 稠密复数矩阵，N = surf.n_rwg()（RWG DOF 数）。
/// 当 aca_tol > 0 时仅返回全矩阵（ACA 压缩留作后续扩展）。
pub fn assemble_calderon(
    surf: &SurfaceMesh,
    freq: f64,
    _aca_tol: f64,
) -> RemResult<DMatrix<Complex64>> {
    use rem_mom::assemble::assemble_efie_pulse;
    use rem_mom::quadrature::TriQuad;

    // 脉冲基函数 EFIE（后续可升级为 RWG）
    let quad = TriQuad::new(4);
    let z = assemble_efie_pulse(surf, freq, &quad, 1e-6)?;
    Ok(z)
}
