//! Anderson 加速（Anderson Mixing）用于加速 Schwarz DDM 迭代收敛。
//!
//! # 算法
//!
//! 给定迭代映射 G：xₙ₊₁ = G(xₙ)，Anderson(m) 维护一个深度为 m 的历史窗口：
//!
//!   ΔX = [xₙ−xₙ₋₁ | ... | xₙ₋ₘ₊₁−xₙ₋ₘ]   (n×m)
//!   ΔF = [fₙ−fₙ₋₁ | ... ]                  f = G(x)−x（残差）
//!
//! 求解最小二乘问题：
//!   θ* = argmin ||fₙ − ΔF·θ||₂
//!
//! 新迭代：
//!   x̃ₙ₊₁ = G(xₙ) − (ΔX + ΔF)·θ*
//!          = xₙ₊₁ − (ΔX + ΔF)·θ*
//!
//! # 参考文献
//! - Walker & Ni (2011) "Anderson Acceleration for Fixed-Point Iterations", SIAM J. Numer. Anal.
//! - Fang & Saad (2009) "Two Classes of Multisecant Methods for Nonlinear Acceleration"
//!
//! # 使用
//! ```rust,ignore
//! let mut aa = AndersonAccelerator::new(depth, n_dofs);
//! // 在 Schwarz 迭代循环中：
//! let x_new = solver(x_old);
//! let x_acc = aa.apply(&x_old, &x_new); // 混合更新
//! x_current = x_acc;
//! ```

use num_complex::Complex64;

/// Anderson 加速器状态
pub struct AndersonAccelerator {
    /// 最大历史深度 m
    depth: usize,
    /// 历史 x 迭代序列（ring buffer，newest-first）
    x_hist: Vec<Vec<Complex64>>,
    /// 历史残差 f = G(x) - x
    f_hist: Vec<Vec<Complex64>>,
    /// 当前窗口大小（≤ depth）
    window: usize,
}

impl AndersonAccelerator {
    /// 创建新的加速器。`depth=0` 退化为无加速的简单迭代。
    pub fn new(depth: usize) -> Self {
        Self {
            depth,
            x_hist: Vec::with_capacity(depth + 1),
            f_hist: Vec::with_capacity(depth + 1),
            window: 0,
        }
    }

    /// 将历史归零（用于重启）。
    pub fn reset(&mut self) {
        self.x_hist.clear();
        self.f_hist.clear();
        self.window = 0;
    }

    /// 接受当前迭代 `x_cur` 和固定点映射输出 `g_cur = G(x_cur)`，
    /// 返回 Anderson 混合后的下一迭代点。
    ///
    /// 当 `depth == 0` 或历史不足时，直接返回 `g_cur`（无加速）。
    pub fn apply(&mut self, x_cur: &[Complex64], g_cur: &[Complex64]) -> Vec<Complex64> {
        if self.depth == 0 {
            return g_cur.to_vec();
        }

        let n = x_cur.len();
        let f_cur: Vec<Complex64> = g_cur.iter().zip(x_cur.iter())
            .map(|(&g, &x)| g - x)
            .collect();

        // Store current (x, f) into history ring buffer (prepend = newest first).
        self.x_hist.insert(0, x_cur.to_vec());
        self.f_hist.insert(0, f_cur.clone());

        // Keep at most depth+1 entries (we need pairs of consecutive entries).
        if self.x_hist.len() > self.depth + 1 {
            self.x_hist.pop();
            self.f_hist.pop();
        }
        self.window = self.x_hist.len().saturating_sub(1);

        if self.window == 0 {
            // Only one entry — no differences yet; return plain G(x).
            return g_cur.to_vec();
        }

        let m = self.window;

        // Build ΔF (n×m) and ΔX (n×m) column-major.
        // ΔF[:,k] = f[k] - f[k+1]   (newest at k=0)
        // ΔX[:,k] = x[k] - x[k+1]
        let mut df: Vec<Vec<Complex64>> = Vec::with_capacity(m);
        let mut dx: Vec<Vec<Complex64>> = Vec::with_capacity(m);
        for k in 0..m {
            let col_df: Vec<Complex64> = (0..n)
                .map(|i| self.f_hist[k][i] - self.f_hist[k + 1][i])
                .collect();
            let col_dx: Vec<Complex64> = (0..n)
                .map(|i| self.x_hist[k][i] - self.x_hist[k + 1][i])
                .collect();
            df.push(col_df);
            dx.push(col_dx);
        }

        // Solve: ΔF · θ ≈ f_cur  in least-squares sense.
        // Use normal equations: (ΔFᴴ ΔF) θ = ΔFᴴ f_cur.
        // For small m (≤ 10), direct Gram-matrix solve is fine.
        let theta = match lstsq_gram(&df, &f_cur, m) {
            Some(t) => t,
            None    => return g_cur.to_vec(), // fallback
        };

        // x̃ₙ₊₁ = G(xₙ) - (ΔX + ΔF)·θ
        let mut x_mix: Vec<Complex64> = g_cur.to_vec();
        for k in 0..m {
            for i in 0..n {
                x_mix[i] -= (dx[k][i] + df[k][i]) * theta[k];
            }
        }
        x_mix
    }

    /// Returns current window depth (0 until enough history is accumulated).
    pub fn current_depth(&self) -> usize { self.window }
}

// ── Least-squares solver via Gram matrix (normal equations) ─────────────────

/// Solve the m×m normal-equation system (ΔFᴴ ΔF) θ = ΔFᴴ f  by Cholesky-like
/// LDLᴴ factorisation.  For small m (≤ 15) this is fast and numerically stable.
///
/// Returns `None` if the Gram matrix is singular (rank-deficient history).
fn lstsq_gram(df: &[Vec<Complex64>], f: &[Complex64], m: usize) -> Option<Vec<Complex64>> {
    let n = f.len();

    // Gram matrix G[i,j] = Σₖ conj(ΔF[k,i]) · ΔF[k,j]
    let mut gram: Vec<Complex64> = vec![Complex64::ZERO; m * m];
    let mut rhs:  Vec<Complex64> = vec![Complex64::ZERO; m];

    for i in 0..m {
        for j in 0..m {
            let mut s = Complex64::ZERO;
            for k in 0..n {
                s += df[i][k].conj() * df[j][k];
            }
            gram[i * m + j] = s;
        }
        // rhs[i] = Σₖ conj(ΔF[k,i]) · f[k]
        let mut r = Complex64::ZERO;
        for k in 0..n {
            r += df[i][k].conj() * f[k];
        }
        rhs[i] = r;
    }

    // Gaussian elimination with partial pivoting on the m×m system.
    gauss_elim_complex(&mut gram, &mut rhs, m)
}

/// Gaussian elimination with partial (column) pivoting for an m×m complex system.
/// Solves in-place; returns the solution vector or `None` if singular.
fn gauss_elim_complex(a: &mut [Complex64], b: &mut [Complex64], m: usize) -> Option<Vec<Complex64>> {
    let mut perm: Vec<usize> = (0..m).collect();

    for col in 0..m {
        // Find pivot row (largest absolute value in column `col` from row `col` onward)
        let pivot_row = (col..m).max_by(|&r1, &r2| {
            a[r1 * m + col].norm().partial_cmp(&a[r2 * m + col].norm()).unwrap()
        })?;

        // Swap rows
        if pivot_row != col {
            for j in 0..m { a.swap(col * m + j, pivot_row * m + j); }
            b.swap(col, pivot_row);
            perm.swap(col, pivot_row);
        }

        let pivot = a[col * m + col];
        if pivot.norm() < 1e-30 {
            return None; // singular
        }

        for row in (col + 1)..m {
            let factor = a[row * m + col] / pivot;
            for j in col..m {
                let sub = factor * a[col * m + j];
                a[row * m + j] -= sub;
            }
            b[row] -= factor * b[col];
        }
    }

    // Back-substitution
    let mut x = vec![Complex64::ZERO; m];
    for i in (0..m).rev() {
        let mut s = b[i];
        for j in (i + 1)..m {
            s -= a[i * m + j] * x[j];
        }
        x[i] = s / a[i * m + i];
    }
    Some(x)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn c(re: f64, im: f64) -> Complex64 { Complex64::new(re, im) }

    /// Anderson(0) is identity (no acceleration).
    #[test]
    fn depth_zero_is_passthrough() {
        let mut aa = AndersonAccelerator::new(0);
        let x = vec![c(1.0, 0.0), c(2.0, 0.0)];
        let g = vec![c(1.5, 0.0), c(2.5, 0.0)];
        let out = aa.apply(&x, &g);
        assert_eq!(out, g);
    }

    /// With a single history entry (window=0 after first call), should also pass through.
    #[test]
    fn first_call_with_depth_gt0_is_passthrough() {
        let mut aa = AndersonAccelerator::new(3);
        let x = vec![c(1.0, 0.0)];
        let g = vec![c(1.1, 0.0)];
        let out = aa.apply(&x, &g);
        assert_eq!(out, g);
    }

    /// For a scalar fixed-point iteration x → 0.5x (converges to 0),
    /// Anderson should produce a mixed iterate closer to 0 than plain G(x).
    #[test]
    fn anderson_accelerates_scalar_contraction() {
        let mut aa = AndersonAccelerator::new(2);
        let mut x = vec![c(1.0, 0.0)];

        // Run 5 steps of G(x) = 0.5x and check Anderson output converges faster.
        let mut plain_x = vec![c(1.0, 0.0)];
        for _ in 0..5 {
            let g: Vec<Complex64> = plain_x.iter().map(|&v| v * 0.5).collect();
            plain_x = g;
        }

        for _ in 0..5 {
            let g: Vec<Complex64> = x.iter().map(|&v| v * c(0.5, 0.0)).collect();
            x = aa.apply(&x, &g);
        }

        // Anderson-accelerated x[0] should be closer to 0
        assert!(x[0].norm() <= plain_x[0].norm() + 1e-10,
            "Anderson x={:.4e}, plain={:.4e}", x[0].norm(), plain_x[0].norm());
    }

    /// Gram/Gauss solver: 2×2 identity system.
    #[test]
    fn gauss_elim_identity_2x2() {
        let mut a = vec![c(1.0,0.0), c(0.0,0.0),
                         c(0.0,0.0), c(1.0,0.0)];
        let mut b = vec![c(3.0,0.0), c(4.0,0.0)];
        let x = gauss_elim_complex(&mut a, &mut b, 2).unwrap();
        assert!((x[0] - c(3.0,0.0)).norm() < 1e-12);
        assert!((x[1] - c(4.0,0.0)).norm() < 1e-12);
    }

    /// Gauss returns None for singular matrix.
    #[test]
    fn gauss_elim_singular_returns_none() {
        let mut a = vec![c(0.0,0.0); 4];
        let mut b = vec![c(1.0,0.0), c(0.0,0.0)];
        assert!(gauss_elim_complex(&mut a, &mut b, 2).is_none());
    }
}
