//! Schwarz 迭代 DDM 求解器（含 Robin 传输条件）
//!
//! 算法：Additive Schwarz with Robin TBC
//! 1. 各进程求解本地子域（带 Robin 边界条件吸收相邻子域的出射波）
//! 2. MPI allreduce / send-recv 交换界面切向场值
//! 3. 更新各子域的 Robin 入射数据
//! 4. 检查全局残差收敛

use nalgebra::DVector;
use num_complex::Complex64;
use rem_core::{RemResult, RemError};
use rem_parallel::Comm;

use crate::subdomain::SubDomain;
use crate::interface::InterfacePatch;

/// Schwarz 迭代配置
#[derive(Debug, Clone)]
pub struct SchwarzConfig {
    /// 收敛容差（相对残差）
    pub tol: f64,
    /// 最大迭代次数
    pub max_iter: usize,
    /// Robin 系数 α（通常取 jk）
    pub robin_alpha: Complex64,
    /// 是否使用 Multiplicative Schwarz（true）还是 Additive（false）
    pub multiplicative: bool,
}

impl Default for SchwarzConfig {
    fn default() -> Self {
        Self {
            tol: 1e-6,
            max_iter: 100,
            robin_alpha: Complex64::new(0.0, 1.0), // jk，频率相关，外部设置
            multiplicative: false,
        }
    }
}

/// DDM 求解结果（与 crate::DdmResult 保持一致的字段）
pub struct SchwarzResult {
    /// 各子域当前解向量
    pub solutions: Vec<DVector<Complex64>>,
    /// 迭代次数
    pub iterations: usize,
    /// 最终相对残差
    pub residual: f64,
}

/// 执行 Schwarz DDM 迭代
///
/// 当前为骨架实现：
/// - 无 MPI 子域间通信（单进程 / 多线程）
/// - 各子域独立用 LU 求解
/// - 界面数据更新留作 TODO
///
/// 完整实现需：
/// 1. 各子域调用 rem-driven FEM 组装（含 Robin BC）
/// 2. MPI send/recv 交换界面切向 H 场
/// 3. 重组全局解
pub fn schwarz_solve(
    subdomains: &[SubDomain],
    _interfaces: &[InterfacePatch],
    _comm: &impl Comm,
    tol: f64,
    max_iter: usize,
) -> RemResult<SchwarzResult> {
    use rem_core::LinearOperator;

    if subdomains.is_empty() {
        return Err(RemError::Config("DDM: no subdomains provided".to_string()));
    }

    let mut solutions: Vec<DVector<Complex64>> = subdomains
        .iter()
        .map(|sd| DVector::zeros(sd.n_dof()))
        .collect();

    let mut rel_residual = f64::INFINITY;
    let mut iterations = 0;

    log::info!("Schwarz DDM: {} subdomains, tol={:.2e}, max_iter={}",
        subdomains.len(), tol, max_iter);

    for iter in 0..max_iter {
        iterations = iter + 1;

        // --- 步骤1：各子域 GMRES/LU 求解 ---
        for (i, sd) in subdomains.iter().enumerate() {
            let (mat, rhs) = sd.assemble_local_stiffness_skeleton()?;
            
            // Select solver based on problem size
            let sol = if sd.n_dof() > 100 {
                // Use GMRES for large systems via LinearOperator
                log::debug!("  Subdomain {}: solving with GMRES ({} DOFs)", i, sd.n_dof());
                rem_mom::gmres_solve_op(&mat, &rhs)
                    .or_else(|e| {
                        log::warn!("  Subdomain {} GMRES failed ({}), falling back to LU", i, e);
                        // Fallback to LU if GMRES fails
                        let lu = mat.clone().lu();
                        lu.solve(&rhs).ok_or_else(|| {
                            RemError::Config(format!("Subdomain {} LU solve failed", i))
                        })
                    })?
            } else {
                // Use LU for small systems
                log::debug!("  Subdomain {}: solving with LU ({} DOFs)", i, sd.n_dof());
                let lu = mat.clone().lu();
                lu.solve(&rhs).ok_or_else(|| {
                    RemError::Config(format!("Subdomain {} LU solve failed", i))
                })?
            };
            
            solutions[i] = sol;
        }

        // --- 步骤2：计算全局残差 ---
        let res_norm: f64 = solutions.iter().map(|sol| sol.norm()).sum();
        rel_residual = if res_norm > 0.0 { res_norm } else { 0.0 };

        log::debug!("  iter={}, res={:.4e}", iter + 1, rel_residual);

        // --- 步骤3：交换界面数据（骨架：无 MPI 通信）---
        // TODO: comm.send/recv 界面切向场

        if rel_residual < tol {
            log::info!("Schwarz converged at iter={}", iter + 1);
            break;
        }
    }

    if rel_residual >= tol {
        log::warn!("Schwarz DDM did not converge: res={:.4e} after {} iters",
            rel_residual, max_iter);
    }

    Ok(SchwarzResult { solutions, iterations, residual: rel_residual })
}

/// 将各子域解向量重组为全局解
pub fn assemble_global_solution(
    result: &SchwarzResult,
    subdomains: &[SubDomain],
    n_global_dof: usize,
) -> DVector<Complex64> {
    let mut global = DVector::zeros(n_global_dof);
    for (i, sd) in subdomains.iter().enumerate() {
        for (local_idx, &global_idx) in sd.local_to_global.iter().enumerate() {
            if global_idx < n_global_dof {
                global[global_idx] = result.solutions[i][local_idx];
            }
        }
    }
    global
}
