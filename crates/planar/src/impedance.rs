//! 阻抗矩阵填充：平面 RWG 基函数的 Z 矩阵元素计算
//!
//! 对于平面结构（z=0），利用 2D FFT 卷积加速相同网格平移不变部分，
//! 近场修正使用精确数值积分。

use num_complex::Complex64;
use crate::grid::PlanarGrid;

/// 阻抗矩阵填充参数
pub struct ImpedanceParams {
    /// 工作频率 (Hz)
    pub freq: f64,
    /// 近场修正半径（单元格数）
    pub near_field_radius: usize,
}

impl ImpedanceParams {
    pub fn new(freq: f64) -> Self {
        Self {
            freq,
            near_field_radius: 3,
        }
    }

    /// 角频率
    pub fn omega(&self) -> f64 {
        2.0 * std::f64::consts::PI * self.freq
    }

    /// 真空中波数 k = omega * sqrt(eps0 * mu0)
    pub fn k0(&self) -> f64 {
        const C0: f64 = 2.997_924_58e8;
        self.omega() / C0
    }
}

/// 格林函数 G(r) = exp(-jkr) / (4 pi r)
pub fn greens_function(r: f64, k: f64) -> Complex64 {
    if r < 1e-12 {
        return Complex64::new(0.0, 0.0);
    }
    let phase = Complex64::new(0.0, -k * r).exp();
    phase / (4.0 * std::f64::consts::PI * r)
}

/// 计算两点之间的距离
pub fn dist(p1: [f64; 3], p2: [f64; 3]) -> f64 {
    let dx = p1[0] - p2[0];
    let dy = p1[1] - p2[1];
    let dz = p1[2] - p2[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// 阻抗矩阵（稠密，用于小规模验证）
pub struct ImpedanceMatrix {
    pub data: Vec<Complex64>,
    pub n: usize,
}

impl ImpedanceMatrix {
    pub fn new(n: usize) -> Self {
        Self {
            data: vec![Complex64::new(0.0, 0.0); n * n],
            n,
        }
    }

    pub fn get(&self, i: usize, j: usize) -> Complex64 {
        self.data[i * self.n + j]
    }

    pub fn set(&mut self, i: usize, j: usize, val: Complex64) {
        self.data[i * self.n + j] = val;
    }

    /// 使用 Gauss-Seidel 迭代求解 Z * x = b（简单实现，生产中用 GMRES）
    pub fn solve_iterative(&self, b: &[Complex64], max_iter: usize) -> Vec<Complex64> {
        let n = self.n;
        let mut x = vec![Complex64::new(0.0, 0.0); n];
        for _iter in 0..max_iter {
            for i in 0..n {
                let mut sum = b[i];
                for j in 0..n {
                    if j != i {
                        sum -= self.get(i, j) * x[j];
                    }
                }
                let diag = self.get(i, i);
                if diag.norm() > 1e-15 {
                    x[i] = sum / diag;
                }
            }
        }
        x
    }
}

/// 填充平面网格的阻抗矩阵（朴素 O(N^2) 版本，用于参考）
pub fn fill_impedance_naive(
    grid: &PlanarGrid,
    params: &ImpedanceParams,
) -> ImpedanceMatrix {
    let n = grid.edges.len();
    let mut zmat = ImpedanceMatrix::new(n);
    let k = params.k0();
    let mu0 = 4.0 * std::f64::consts::PI * 1e-7_f64;
    let eps0 = 8.854_187_817e-12_f64;
    let omega = params.omega();
    let jomegamu = Complex64::new(0.0, omega * mu0);
    let inv_jomegaeps = Complex64::new(0.0, -1.0 / (omega * eps0));

    for i in 0..n {
        let ei = &grid.edges[i];
        for j in 0..n {
            let ej = &grid.edges[j];
            // 使用中心点近似（实际应用需高斯积分）
            let ci = ei.center;
            let cj = ej.center;
            let r = dist(ci, cj);
            let g = greens_function(r, k);

            // EFIE 矩阵元素（简化版）:
            // Z_ij = jω μ₀ ∫∫ f_i · f_j G dr dr' - (1/jω ε₀) ∫∫ (∇·f_i)(∇·f_j) G dr dr'
            let dot_ff = ei.tangent[0] * ej.tangent[0]
                + ei.tangent[1] * ej.tangent[1];
            let div_product = ei.div_rho * ej.div_rho;

            let z_ij = jomegamu * dot_ff * g * ei.length * ej.length
                + inv_jomegaeps * div_product * g * ei.length * ej.length;
            zmat.set(i, j, z_ij);
        }
    }
    zmat
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grid::{PlanarGrid, PlanarEdge};

    fn make_simple_grid() -> PlanarGrid {
        // 两条简单边
        let edges = vec![
            PlanarEdge {
                center: [0.0, 0.0, 0.0],
                tangent: [1.0, 0.0, 0.0],
                length: 0.1,
                div_rho: 1.0,
            },
            PlanarEdge {
                center: [0.5, 0.0, 0.0],
                tangent: [1.0, 0.0, 0.0],
                length: 0.1,
                div_rho: 1.0,
            },
        ];
        PlanarGrid { edges, nx: 1, ny: 1, dx: 0.5, dy: 0.5 }
    }

    #[test]
    fn test_impedance_shape() {
        let grid = make_simple_grid();
        let params = ImpedanceParams::new(1e9);
        let zmat = fill_impedance_naive(&grid, &params);
        assert_eq!(zmat.n, 2);
    }

    #[test]
    fn test_greens_function_far() {
        let g = greens_function(1.0, 1.0);
        // |G| ≈ 1/(4π) ≈ 0.0796
        assert!((g.norm() - 1.0 / (4.0 * std::f64::consts::PI)).abs() < 1e-3);
    }
}
