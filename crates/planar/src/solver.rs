//! 平面 MoM 求解器：组装阻抗矩阵并求解 [Z]{I} = {V}

use crate::grid::PlanarGrid;
use crate::impedance::{fill_impedance_naive, ImpedanceParams};
use nalgebra::{DMatrix, DVector};
use num_complex::Complex64;

/// 平面 MoM 求解结果
pub struct MomSolution {
    /// 各基函数的展开系数（电流）
    pub coefficients: DVector<Complex64>,
    /// 频率 (Hz)
    pub frequency: f64,
}

/// 求解器配置
pub struct SolverConfig {
    /// 使用 FFT 加速卷积（仅适用于均匀网格）
    pub use_fft: bool,
    /// 迭代求解最大迭代次数（0 = 直接求解）
    pub max_iter: usize,
    /// 迭代收敛容差
    pub tol: f64,
}

impl Default for SolverConfig {
    fn default() -> Self {
        Self {
            use_fft: true,
            max_iter: 0,
            tol: 1e-6,
        }
    }
}

/// 平面 MoM 求解器
pub struct PlanarMomSolver {
    pub grid: PlanarGrid,
    pub config: SolverConfig,
}

impl PlanarMomSolver {
    pub fn new(grid: PlanarGrid, config: SolverConfig) -> Self {
        Self { grid, config }
    }

    /// 求解：给定激励向量 V（长度 = n_basis），返回电流系数
    pub fn solve(&self, frequency: f64, excitation: &DVector<Complex64>) -> MomSolution {
        let params = ImpedanceParams::new(frequency);
        let zmat = fill_impedance_naive(&self.grid, &params);
        let n = zmat.n;

        // 转换为 nalgebra DMatrix
        let z_na = DMatrix::from_fn(n, n, |i, j| zmat.get(i, j));

        let coefficients = if self.config.max_iter == 0 {
            solve_direct(&z_na, excitation)
        } else {
            solve_steepest_descent(&z_na, excitation, self.config.max_iter, self.config.tol)
        };

        MomSolution {
            coefficients,
            frequency,
        }
    }
}

/// 直接 LU 分解求解（小规模，nalgebra 内置）
fn solve_direct(z: &DMatrix<Complex64>, v: &DVector<Complex64>) -> DVector<Complex64> {
    // nalgebra LU 分解
    match z.clone().lu().solve(v) {
        Some(x) => x,
        None => DVector::zeros(v.len()),
    }
}

/// 最速下降迭代（大规模近似，实际项目中应用 GMRES）
fn solve_steepest_descent(
    z: &DMatrix<Complex64>,
    v: &DVector<Complex64>,
    max_iter: usize,
    tol: f64,
) -> DVector<Complex64> {
    let n = v.len();
    let mut x = DVector::<Complex64>::zeros(n);

    let b_norm = v.norm();
    if b_norm < 1e-30 {
        return x;
    }

    for _iter in 0..max_iter {
        let r = v - z * &x;
        let r_norm = r.norm();
        if r_norm / b_norm < tol {
            break;
        }
        let ar = z * &r;
        let rr: Complex64 = r.iter().zip(r.iter()).map(|(a, b)| a.conj() * b).sum();
        let rar: Complex64 = r.iter().zip(ar.iter()).map(|(a, b)| a.conj() * b).sum();
        if rar.norm() < 1e-30 {
            break;
        }
        let alpha = rr / rar;
        x += &r * alpha;
    }

    x
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grid::PlanarGrid;

    #[test]
    fn test_solver_small() {
        let grid = PlanarGrid::new(1.0, 1.0, 2, 2);
        let solver = PlanarMomSolver::new(grid, SolverConfig::default());
        let n = solver.grid.edges.len();
        let excitation = DVector::from_element(n, Complex64::new(1.0, 0.0));
        let sol = solver.solve(1e9, &excitation);
        assert_eq!(sol.coefficients.len(), n);
    }
}
