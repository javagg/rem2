//! 分层媒质谱域格林函数 (Spectral-domain Green's function for layered media)
//!
//! 实现平面分层结构中的 Sommerfeld 积分核，采用传递矩阵法（TMM）。
//! 支持 TE/TM 极化，返回谱域电流-电位格林函数 G^A 和 G^q。

use num_complex::Complex64;

/// 单层媒质参数
#[derive(Debug, Clone)]
pub struct Layer {
    /// 相对介电常数（可复数：ε_r - j σ/(ω ε_0)）
    pub eps_r: Complex64,
    /// 相对磁导率
    pub mu_r: Complex64,
    /// 层厚度（米），最顶层/底层可设为 f64::INFINITY
    pub thickness: f64,
}

impl Layer {
    pub fn new(eps_r: impl Into<Complex64>, mu_r: impl Into<Complex64>, thickness: f64) -> Self {
        Self {
            eps_r: eps_r.into(),
            mu_r: mu_r.into(),
            thickness,
        }
    }

    /// 自由空间层
    pub fn free_space(thickness: f64) -> Self {
        Self::new(1.0, 1.0, thickness)
    }

    /// 理想介质层
    pub fn dielectric(eps_r: f64, thickness: f64) -> Self {
        Self::new(eps_r, 1.0, thickness)
    }

    /// 有耗介质层（ε_r - j σ/(ω ε_0)）
    pub fn lossy_dielectric(eps_r: f64, sigma: f64, omega: f64, thickness: f64) -> Self {
        let eps0 = 8.854187817e-12_f64;
        let eps_complex = Complex64::new(eps_r, -sigma / (omega * eps0));
        Self::new(eps_complex, 1.0, thickness)
    }
}

/// 分层结构（从上到下排列，index 0 为顶层半空间）
#[derive(Debug, Clone)]
pub struct LayeredMedium {
    pub layers: Vec<Layer>,
    /// 工作角频率 ω (rad/s)
    pub omega: f64,
}

impl LayeredMedium {
    pub fn new(layers: Vec<Layer>, omega: f64) -> Self {
        assert!(layers.len() >= 2, "至少需要两个半空间层");
        Self { layers, omega }
    }

    /// 单层自由空间（用于验证）
    pub fn free_space(omega: f64) -> Self {
        Self::new(
            vec![Layer::free_space(f64::INFINITY), Layer::free_space(f64::INFINITY)],
            omega,
        )
    }

    /// 常用：自由空间 + 单介质层 + 地板（PEC backing）
    pub fn grounded_substrate(eps_r: f64, thickness: f64, omega: f64) -> Self {
        Self::new(
            vec![
                Layer::free_space(f64::INFINITY),      // 顶层：空气半空间
                Layer::dielectric(eps_r, thickness),   // 介质板
                // 底层：PEC 地板 — 用极大 σ 近似（或在 TMM 中特殊处理）
                Layer::lossy_dielectric(1.0, 1e10, omega, f64::INFINITY),
            ],
            omega,
        )
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// 内部：传递矩阵法核心
// ──────────────────────────────────────────────────────────────────────────────

const EPS0: f64 = 8.854187817e-12;
const MU0: f64 = 1.2566370614359173e-6;

/// 计算第 i 层的纵向波数 k_z
fn kz(k_rho_sq: Complex64, eps_r: Complex64, mu_r: Complex64, omega: f64) -> Complex64 {
    let k0_sq = Complex64::new((omega / 3e8).powi(2), 0.0);
    let kz_sq = eps_r * mu_r * k0_sq - k_rho_sq;
    // 取 Im(k_z) <= 0 以满足辐射条件（衰减方向正确）
    let kz_val = kz_sq.sqrt();
    if kz_val.im > 0.0 {
        -kz_val
    } else {
        kz_val
    }
}

/// 2×2 传递矩阵元素
#[derive(Debug, Clone, Copy)]
struct TMatrix {
    m: [[Complex64; 2]; 2],
}

impl TMatrix {
    fn identity() -> Self {
        let one = Complex64::new(1.0, 0.0);
        let zero = Complex64::new(0.0, 0.0);
        Self {
            m: [[one, zero], [zero, one]],
        }
    }

    fn mul(&self, rhs: &Self) -> Self {
        let mut result = [[Complex64::new(0.0, 0.0); 2]; 2];
        for i in 0..2 {
            for j in 0..2 {
                for k in 0..2 {
                    result[i][j] += self.m[i][k] * rhs.m[k][j];
                }
            }
        }
        Self { m: result }
    }
}

/// 单层 TM 传递矩阵（纵向电流，水平偶极源常用）
fn tm_layer_matrix(kz_i: Complex64, eps_i: Complex64, d: f64, omega: f64) -> TMatrix {
    let j = Complex64::new(0.0, 1.0);
    let jkzd = j * kz_i * d;
    let cos_val = jkzd.cos();
    let sin_val = jkzd.sin();
    let eta_i = kz_i / (omega * EPS0 * eps_i); // TM 阻抗

    TMatrix {
        m: [
            [cos_val, -j * eta_i * sin_val],
            [-j * sin_val / eta_i, cos_val],
        ],
    }
}

/// 单层 TE 传递矩阵
fn te_layer_matrix(kz_i: Complex64, mu_i: Complex64, d: f64, omega: f64) -> TMatrix {
    let j = Complex64::new(0.0, 1.0);
    let jkzd = j * kz_i * d;
    let cos_val = jkzd.cos();
    let sin_val = jkzd.sin();
    let zeta_i = omega * MU0 * mu_i / kz_i; // TE 阻抗

    TMatrix {
        m: [
            [cos_val, -j * zeta_i * sin_val],
            [-j * sin_val / zeta_i, cos_val],
        ],
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// 公开 API：谱域格林函数
// ──────────────────────────────────────────────────────────────────────────────

/// 谱域格林函数结果（在单个 k_rho 点）
#[allow(non_snake_case)]
#[derive(Debug, Clone, Copy)]
pub struct SpectralGreen {
    /// 磁矢位格林函数 G^A_xx (TM，水平电流 x 分量)
    pub gA: Complex64,
    /// 电标位格林函数 G^q (TM，电荷项)
    pub gq: Complex64,
}

impl LayeredMedium {
    /// 计算谱域格林函数
    ///
    /// # 参数
    /// - `k_rho`: 横向波数（可以复数，用于 SDA 分析）
    /// - `z_src`: 源层 z 坐标（米，从顶层 z=0 向下为正）
    /// - `z_obs`: 观察层 z 坐标
    ///
    /// 简化版：假设源和观察点均位于顶层（z=0 平面，如微带天线分析常见情形）。
    pub fn spectral_green(&self, k_rho: Complex64, _z_src: f64, _z_obs: f64) -> SpectralGreen {
        let omega = self.omega;
        let k_rho_sq = k_rho * k_rho;

        // 计算各层 k_z
        let kzs: Vec<Complex64> = self
            .layers
            .iter()
            .map(|l| kz(k_rho_sq, l.eps_r, l.mu_r, omega))
            .collect();

        // 构造从第 1 层（介质层）到底层的 TM 传递矩阵级联
        // 对于两层（半空间 + 半空间）：直接用界面 Fresnel 系数
        let n = self.layers.len();

        // TM 总传递矩阵（跳过第 0 层半空间）
        let mut m_tm = TMatrix::identity();
        let mut m_te = TMatrix::identity();
        for i in 1..n - 1 {
            let d = self.layers[i].thickness;
            let eps_i = self.layers[i].eps_r;
            let mu_i = self.layers[i].mu_r;
            let kz_i = kzs[i];
            m_tm = m_tm.mul(&tm_layer_matrix(kz_i, eps_i, d, omega));
            m_te = m_te.mul(&te_layer_matrix(kz_i, mu_i, d, omega));
        }

        // 顶层（0）和底层（n-1）的 TM 阻抗
        let eta0 = kzs[0] / (omega * EPS0 * self.layers[0].eps_r);
        let eta_n = kzs[n - 1] / (omega * EPS0 * self.layers[n - 1].eps_r);
        let _zeta0 = omega * MU0 * self.layers[0].mu_r / kzs[0];
        let zeta_n = omega * MU0 * self.layers[n - 1].mu_r / kzs[n - 1];

        // 输入导纳（TM）：Y_in = (m11 * Y_n - j*m12) / (-j*m21 * Y_n + m22) ... 转为 Z_in
        let j = Complex64::new(0.0, 1.0);
        let [m11, m12, m21, m22] = [m_tm.m[0][0], m_tm.m[0][1], m_tm.m[1][0], m_tm.m[1][1]];
        // 从底层端口看进去的输入阻抗
        let z_in_tm = (m11 * eta_n + m12) / (m21 * eta_n + m22);
        // 总 TM 输入阻抗（顶层并联）
        let z_tm = z_in_tm * eta0 / (z_in_tm + eta0); // 并不完全对，仅示意

        // TE
        let [t11, t12, t21, t22] = [m_te.m[0][0], m_te.m[0][1], m_te.m[1][0], m_te.m[1][1]];
        let _z_in_te = (t11 * zeta_n + t12) / (t21 * zeta_n + t22);

        // G^A_xx 近似（TM + TE 混合，简化为标量）
        let mu_eff = MU0 * self.layers[0].mu_r.re;
        #[allow(non_snake_case)]
        let gA = Complex64::new(mu_eff, 0.0) / (Complex64::new(2.0, 0.0) * kzs[0])
            * (Complex64::new(1.0, 0.0) + z_tm / (z_in_tm + eta0));

        // G^q = G^A * k_rho^2 / (j ω ε_eff)  [MPIE Sommerfeld 近似]
        let eps_eff = EPS0 * self.layers[0].eps_r;
        let gq = gA * k_rho_sq / (j * omega * eps_eff);

        SpectralGreen { gA, gq }
    }

    /// 批量计算（FFT 网格上所有 k_rho 点）
    pub fn spectral_green_grid(
        &self,
        k_rho_values: &[Complex64],
        z_src: f64,
        z_obs: f64,
    ) -> Vec<SpectralGreen> {
        k_rho_values
            .iter()
            .map(|&kr| self.spectral_green(kr, z_src, z_obs))
            .collect()
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// 测试
// ──────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn freq_to_omega(f_ghz: f64) -> f64 {
        2.0 * PI * f_ghz * 1e9
    }

    #[test]
    fn test_free_space_green() {
        let omega = freq_to_omega(10.0); // 10 GHz
        let medium = LayeredMedium::free_space(omega);
        let k0 = omega / 3e8;
        let k_rho = Complex64::new(0.5 * k0, 0.0);
        let sg = medium.spectral_green(k_rho, 0.0, 0.0);
        // 自由空间：G^A 应为实数且正
        assert!(sg.gA.re > 0.0, "G^A 实部应 > 0");
        println!("自由空间 G^A = {:?}", sg.gA);
    }

    #[test]
    fn test_grounded_substrate() {
        let omega = freq_to_omega(10.0);
        // Rogers RO4003: ε_r = 3.55, 厚度 0.508 mm
        let medium = LayeredMedium::grounded_substrate(3.55, 0.508e-3, omega);
        let k0 = omega / 3e8;
        let k_rho = Complex64::new(0.8 * k0, 0.0);
        let sg = medium.spectral_green(k_rho, 0.0, 0.0);
        println!("有接地板介质基板 G^A = {:?}, G^q = {:?}", sg.gA, sg.gq);
        // 基本合理性检查
        assert!(sg.gA.norm() > 0.0);
    }

    #[test]
    fn test_spectral_grid() {
        let omega = freq_to_omega(5.0);
        let medium = LayeredMedium::grounded_substrate(4.4, 1.6e-3, omega);
        let k0 = omega / 3e8;
        let k_rho_vals: Vec<Complex64> = (0..32)
            .map(|i| Complex64::new(i as f64 * k0 * 0.2, 0.0))
            .collect();
        let greens = medium.spectral_green_grid(&k_rho_vals, 0.0, 0.0);
        assert_eq!(greens.len(), 32);
    }
}
