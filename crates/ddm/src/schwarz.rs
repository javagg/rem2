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
use crate::interface::{apply_robin_to_diagonal, InterfaceExchange, InterfacePatch};

fn apply_robin_diagonal_terms(
    mat: &mut nalgebra::DMatrix<Complex64>,
    patches: &[&InterfacePatch],
) {
    let n = mat.nrows().min(mat.ncols());
    if n == 0 || patches.is_empty() {
        return;
    }

    let mut diag: Vec<Complex64> = (0..n).map(|i| mat[(i, i)]).collect();
    for patch in patches {
        let unit_areas = vec![1.0_f64; patch.n_dofs()];
        apply_robin_to_diagonal(&mut diag, patch, &unit_areas);
    }
    for i in 0..n {
        mat[(i, i)] = diag[i];
    }
}

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
    interfaces: &[InterfacePatch],
    comm: &impl Comm,
    tol: f64,
    max_iter: usize,
) -> RemResult<SchwarzResult> {
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
    log::info!("Schwarz DDM: {} directed interface patches", interfaces.len());

    for iter in 0..max_iter {
        iterations = iter + 1;
        let prev_solutions = solutions.clone();

        // --- 步骤1：各子域 GMRES/LU 求解 ---
        for (i, sd) in subdomains.iter().enumerate() {
            let (mut mat, mut rhs) = sd.assemble_local_stiffness_skeleton()?;
            let owner_patches: Vec<&InterfacePatch> = interfaces
                .iter()
                .filter(|p| p.owner_rank == i as i32)
                .collect();

            // Robin diagonal term on local operator (unit-area placeholder in skeleton mode).
            apply_robin_diagonal_terms(&mut mat, &owner_patches);

            // Robin interface update: use previous-iteration neighbor trace as incoming field.
            for patch in owner_patches {
                let neighbor = patch.neighbor_rank as usize;
                if neighbor >= subdomains.len() {
                    continue;
                }

                let incoming_e: Vec<Complex64> = patch.global_node_ids.iter()
                    .map(|gid| {
                        subdomains[neighbor]
                            .global_to_local
                            .get(gid)
                            .and_then(|&lid| prev_solutions.get(neighbor).and_then(|v| v.get(lid)).copied())
                            .unwrap_or(Complex64::ZERO)
                    })
                    .collect();

                let exch = InterfaceExchange {
                    incoming_e,
                    incoming_h: vec![Complex64::ZERO; patch.n_dofs()],
                };
                let contrib = exch.robin_rhs_contribution(patch.robin_alpha);

                for (&ldof, &val) in patch.local_dofs.iter().zip(contrib.iter()) {
                    if ldof < rhs.len() {
                        rhs[ldof] += val;
                    }
                }
            }
            
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

        // --- 步骤2：计算界面更新残差 ---
        let mut delta_sq = 0.0_f64;
        let mut ref_sq = 0.0_f64;
        for patch in interfaces {
            let owner = patch.owner_rank as usize;
            if owner >= subdomains.len() {
                continue;
            }
            for &ldof in &patch.local_dofs {
                if ldof < solutions[owner].len() {
                    let cur = solutions[owner][ldof];
                    let prev = prev_solutions[owner][ldof];
                    delta_sq += (cur - prev).norm_sqr();
                    ref_sq += cur.norm_sqr();
                }
            }
        }
        let delta_sq_glob = comm.allreduce_f64(delta_sq);
        let ref_sq_glob = comm.allreduce_f64(ref_sq);
        rel_residual = (delta_sq_glob.sqrt()) / (ref_sq_glob.sqrt().max(1e-30));

        log::debug!("  iter={}, res={:.4e}", iter + 1, rel_residual);

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

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    fn mock_subdomain(id: usize, gid: usize) -> SubDomain {
        let mut global_to_local = HashMap::new();
        global_to_local.insert(gid, 0);
        SubDomain {
            id,
            volume_elements: vec![],
            boundary_elements: vec![],
            global_to_local,
            local_to_global: vec![gid],
            interface_nodes: vec![],
            interface_neighbor: vec![],
        }
    }

    #[test]
    fn robin_diagonal_terms_are_applied_on_owner_dofs() {
        let mut mat = nalgebra::DMatrix::<Complex64>::identity(3, 3);
        let patch = InterfacePatch::new(
            0,
            vec![0, 2],
            vec![10, 20],
            1,
            Complex64::new(0.0, 2.0),
        );
        apply_robin_diagonal_terms(&mut mat, &[&patch]);

        assert_eq!(mat[(0, 0)], Complex64::new(1.0, 2.0));
        assert_eq!(mat[(1, 1)], Complex64::new(1.0, 0.0));
        assert_eq!(mat[(2, 2)], Complex64::new(1.0, 2.0));
    }

    #[test]
    fn schwarz_runs_with_bidirectional_interface_patches() {
        let subdomains = vec![mock_subdomain(0, 0), mock_subdomain(1, 0)];
        let interfaces = vec![
            InterfacePatch::new(0, vec![0], vec![0], 1, Complex64::new(0.0, 1.0)),
            InterfacePatch::new(1, vec![0], vec![0], 0, Complex64::new(0.0, 1.0)),
        ];

        let out = schwarz_solve(&subdomains, &interfaces, &rem_parallel::NoComm, 1e-9, 5)
            .expect("schwarz solve should complete");

        assert_eq!(out.solutions.len(), 2);
        assert!(out.iterations >= 1);
        assert!(out.residual.is_finite());
    }
}
