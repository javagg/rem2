//! Schwarz 迭代 DDM 求解器（含 Robin 传输条件）
//!
//! # 实现特性
//!
//! | 特性 | 状态 |
//! |------|------|
//! | Additive Schwarz + Robin TBC | ✅ |
//! | **Multiplicative Schwarz** | ✅ (新) |
//! | **Anderson 加速** | ✅ (新) |
//! | **P1 Helmholtz FEM 本地组装** | ✅ (新，via SubDomain::assemble_local_p1_helmholtz) |
//! | MPI allreduce 残差聚合 | ✅ |
//!
//! ## Multiplicative Schwarz
//! 每当子域 i 完成求解后，立即将最新解传递给后续子域 i+1 使用（而不是等到本轮迭代结束）。
//! 收敛速度通常比 Additive 快 2–5 倍，代价是子域求解串行化。
//!
//! ## Anderson(m) 加速
//! 维护深度为 m 的历史窗口，用小型最小二乘问题在界面 DOF 向量上做混合，
//! 典型地减少 50–70% 的 Schwarz 外迭代次数。

use nalgebra::DVector;
use num_complex::Complex64;
use rem_core::{RemResult, RemError};
use rem_parallel::Comm;

use crate::anderson::AndersonAccelerator;
use crate::interface::{apply_robin_to_diagonal, InterfaceExchange, InterfacePatch};
use crate::subdomain::SubDomain;

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

/// Schwarz 求解配置
#[derive(Debug, Clone)]
pub struct SchwarzConfig {
    /// 收敛容差（相对残差）
    pub tol: f64,
    /// 最大迭代次数
    pub max_iter: usize,
    /// Robin 系数 α（通常取 jk₀）
    pub robin_alpha: Complex64,
    /// Multiplicative Schwarz（true）或 Additive（false）
    pub multiplicative: bool,
    /// Anderson 加速深度（0 = 关闭）
    pub anderson_depth: usize,
    /// 频率 [Hz]，用于真实 Helmholtz FEM 组装
    pub freq_hz: f64,
    /// 相对介电常数
    pub eps_r: f64,
    /// 相对磁导率
    pub mu_r: f64,
}

impl Default for SchwarzConfig {
    fn default() -> Self {
        Self {
            tol: 1e-6,
            max_iter: 100,
            robin_alpha: Complex64::new(0.0, 1.0),
            multiplicative: false,
            anderson_depth: 5,
            freq_hz: 1e9,
            eps_r: 1.0,
            mu_r: 1.0,
        }
    }
}

/// DDM 求解结果
pub struct SchwarzResult {
    pub solutions: Vec<DVector<Complex64>>,
    pub iterations: usize,
    pub residual: f64,
}

/// 执行 Schwarz DDM 迭代
///
/// 支持：
/// - Additive / Multiplicative Schwarz（由 `cfg.multiplicative` 控制）
/// - Anderson(m) 加速（由 `cfg.anderson_depth` 控制，0 = 关闭）
/// - 真实 P1 Helmholtz FEM 本地组装（通过 `mesh` + `cfg.freq_hz`）
///
/// 当 `mesh` 为 `None` 时退化为骨架（单位矩阵）求解，用于测试。
pub fn schwarz_solve(
    subdomains: &[SubDomain],
    interfaces: &[InterfacePatch],
    comm: &impl Comm,
    tol: f64,
    max_iter: usize,
) -> RemResult<SchwarzResult> {
    let cfg = SchwarzConfig { tol, max_iter, ..Default::default() };
    schwarz_solve_cfg(subdomains, interfaces, comm, &cfg, None)
}

/// 完整配置版入口（提供网格和配置）
pub fn schwarz_solve_full<'m>(
    subdomains: &[SubDomain],
    interfaces: &[InterfacePatch],
    comm: &impl Comm,
    cfg: &SchwarzConfig,
    mesh: Option<&'m rem_mesh::RemMesh>,
) -> RemResult<SchwarzResult> {
    schwarz_solve_cfg(subdomains, interfaces, comm, cfg, mesh)
}

fn schwarz_solve_cfg(
    subdomains: &[SubDomain],
    interfaces: &[InterfacePatch],
    comm: &impl Comm,
    cfg: &SchwarzConfig,
    mesh: Option<&rem_mesh::RemMesh>,
) -> RemResult<SchwarzResult> {
    if subdomains.is_empty() {
        return Err(RemError::Config("DDM: no subdomains provided".to_string()));
    }

    let n_sub = subdomains.len();
    let mut solutions: Vec<DVector<Complex64>> = subdomains
        .iter()
        .map(|sd| DVector::zeros(sd.n_dof()))
        .collect();

    // Anderson accelerator operates on the flattened interface-DOF vector.
    let total_iface_dofs: usize = interfaces.iter().map(|p| p.n_dofs()).sum();
    let mut anderson = AndersonAccelerator::new(
        if total_iface_dofs > 0 { cfg.anderson_depth } else { 0 }
    );

    let mut rel_residual = f64::INFINITY;
    let mut iterations = 0;

    log::info!(
        "Schwarz DDM: {} subdomains | {} | Anderson({}), tol={:.2e}, max_iter={}",
        n_sub,
        if cfg.multiplicative { "Multiplicative" } else { "Additive" },
        cfg.anderson_depth,
        cfg.tol,
        cfg.max_iter,
    );

    for iter in 0..cfg.max_iter {
        iterations = iter + 1;
        let prev_solutions = solutions.clone();

        // ── Per-subdomain solve ──────────────────────────────────────────────
        for i in 0..n_sub {
            let sd = &subdomains[i];

            // Assemble local system: real FEM if mesh available, else skeleton.
            let (mut mat, mut rhs) = if let Some(m) = mesh {
                sd.assemble_local_p1_helmholtz(m, cfg.freq_hz, cfg.eps_r, cfg.mu_r)?
            } else {
                sd.assemble_local_stiffness_skeleton()?
            };

            let owner_patches: Vec<&InterfacePatch> = interfaces
                .iter()
                .filter(|p| p.owner_rank == i as i32)
                .collect();

            // Robin diagonal correction on interface DOFs.
            apply_robin_diagonal_terms(&mut mat, &owner_patches);

            // Robin interface update from neighbor traces.
            // Multiplicative: use the already-updated solutions[neighbor] (latest).
            // Additive: use prev_solutions[neighbor] (frozen at start of this iter).
            for patch in &owner_patches {
                let neighbor = patch.neighbor_rank as usize;
                if neighbor >= n_sub { continue; }

                let src = if cfg.multiplicative { &solutions } else { &prev_solutions };
                let incoming_e: Vec<Complex64> = patch.global_node_ids.iter()
                    .map(|gid| {
                        subdomains[neighbor]
                            .global_to_local
                            .get(gid)
                            .and_then(|&lid| src.get(neighbor).and_then(|v| v.get(lid)).copied())
                            .unwrap_or(Complex64::ZERO)
                    })
                    .collect();

                let exch = InterfaceExchange {
                    incoming_e,
                    incoming_h: vec![Complex64::ZERO; patch.n_dofs()],
                };
                let contrib = exch.robin_rhs_contribution(patch.robin_alpha);
                for (&ldof, &val) in patch.local_dofs.iter().zip(contrib.iter()) {
                    if ldof < rhs.len() { rhs[ldof] += val; }
                }
            }

            // Solve subdomain system (LU → GMRES fallback).
            let sol = if sd.n_dof() > 100 {
                log::debug!("  Subdomain {i}: GMRES ({} DOFs)", sd.n_dof());
                rem_mom::gmres_solve_op(&mat, &rhs)
                    .or_else(|e| {
                        log::warn!("  Subdomain {i} GMRES failed ({e}), falling back to LU");
                        mat.clone().lu().solve(&rhs).ok_or_else(|| {
                            RemError::Config(format!("Subdomain {i} LU solve failed"))
                        })
                    })?
            } else {
                log::debug!("  Subdomain {i}: LU ({} DOFs)", sd.n_dof());
                mat.clone().lu().solve(&rhs).ok_or_else(|| {
                    RemError::Config(format!("Subdomain {i} LU solve failed"))
                })?
            };
            solutions[i] = sol;
        }

        // ── Anderson acceleration on interface DOF vector ────────────────────
        if cfg.anderson_depth > 0 && total_iface_dofs > 0 {
            // Extract old and new interface DOF values.
            let x_iface = collect_iface_dofs(interfaces, subdomains, &prev_solutions);
            let g_iface = collect_iface_dofs(interfaces, subdomains, &solutions);

            let x_mixed = anderson.apply(&x_iface, &g_iface);

            // Write mixed interface DOF values back into solutions.
            scatter_iface_dofs(&x_mixed, interfaces, subdomains, &mut solutions);
        }

        // ── Convergence check ────────────────────────────────────────────────
        let mut delta_sq = 0.0_f64;
        let mut ref_sq   = 0.0_f64;
        for patch in interfaces {
            let owner = patch.owner_rank as usize;
            if owner >= n_sub { continue; }
            for &ldof in &patch.local_dofs {
                if ldof < solutions[owner].len() {
                    let cur  = solutions[owner][ldof];
                    let prev = prev_solutions[owner][ldof];
                    delta_sq += (cur - prev).norm_sqr();
                    ref_sq   += cur.norm_sqr();
                }
            }
        }
        let delta_sq_glob = comm.allreduce_f64(delta_sq);
        let ref_sq_glob   = comm.allreduce_f64(ref_sq);
        rel_residual = delta_sq_glob.sqrt() / ref_sq_glob.sqrt().max(1e-30);

        log::debug!("  iter={}, res={:.4e}, Anderson_depth={}",
            iter + 1, rel_residual, anderson.current_depth());

        if rel_residual < cfg.tol {
            log::info!("Schwarz converged at iter={} (res={:.4e})", iter + 1, rel_residual);
            break;
        }
    }

    if rel_residual >= cfg.tol {
        log::warn!("Schwarz DDM did not converge: res={:.4e} after {} iters",
            rel_residual, iterations);
    }

    Ok(SchwarzResult { solutions, iterations, residual: rel_residual })
}

// ── Anderson interface-DOF helpers ───────────────────────────────────────────

/// Flatten all owner interface DOF values into a single vector.
fn collect_iface_dofs(
    interfaces: &[InterfacePatch],
    subdomains: &[SubDomain],
    solutions: &[DVector<Complex64>],
) -> Vec<Complex64> {
    let mut out = Vec::new();
    for patch in interfaces {
        let owner = patch.owner_rank as usize;
        if owner >= subdomains.len() { continue; }
        for &ldof in &patch.local_dofs {
            let v = solutions.get(owner)
                .and_then(|s| s.get(ldof).copied())
                .unwrap_or(Complex64::ZERO);
            out.push(v);
        }
    }
    out
}

/// Scatter mixed interface DOF values back into subdomain solution vectors.
fn scatter_iface_dofs(
    iface_vals: &[Complex64],
    interfaces: &[InterfacePatch],
    subdomains: &[SubDomain],
    solutions: &mut [DVector<Complex64>],
) {
    let mut idx = 0;
    for patch in interfaces {
        let owner = patch.owner_rank as usize;
        if owner >= subdomains.len() { continue; }
        for &ldof in &patch.local_dofs {
            if idx >= iface_vals.len() { break; }
            if let Some(s) = solutions.get_mut(owner) {
                if ldof < s.len() {
                    s[ldof] = iface_vals[idx];
                }
            }
            idx += 1;
        }
    }
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
    use crate::interface::InterfacePatch;
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
            0, vec![0, 2], vec![10, 20], 1, Complex64::new(0.0, 2.0),
        );
        apply_robin_diagonal_terms(&mut mat, &[&patch]);
        assert_eq!(mat[(0, 0)], Complex64::new(1.0, 2.0));
        assert_eq!(mat[(1, 1)], Complex64::new(1.0, 0.0));
        assert_eq!(mat[(2, 2)], Complex64::new(1.0, 2.0));
    }

    #[test]
    fn additive_schwarz_runs_bidirectional() {
        let subs = vec![mock_subdomain(0, 0), mock_subdomain(1, 0)];
        let ifaces = vec![
            InterfacePatch::new(0, vec![0], vec![0], 1, Complex64::new(0.0, 1.0)),
            InterfacePatch::new(1, vec![0], vec![0], 0, Complex64::new(0.0, 1.0)),
        ];
        let out = schwarz_solve(&subs, &ifaces, &rem_parallel::NoComm, 1e-9, 5)
            .expect("additive schwarz should complete");
        assert_eq!(out.solutions.len(), 2);
        assert!(out.residual.is_finite());
    }

    #[test]
    fn multiplicative_schwarz_runs_bidirectional() {
        let subs = vec![mock_subdomain(0, 0), mock_subdomain(1, 0)];
        let ifaces = vec![
            InterfacePatch::new(0, vec![0], vec![0], 1, Complex64::new(0.0, 1.0)),
            InterfacePatch::new(1, vec![0], vec![0], 0, Complex64::new(0.0, 1.0)),
        ];
        let cfg = SchwarzConfig {
            multiplicative: true,
            tol: 1e-9, max_iter: 5,
            anderson_depth: 0,
            ..Default::default()
        };
        let out = schwarz_solve_full(&subs, &ifaces, &rem_parallel::NoComm, &cfg, None)
            .expect("multiplicative schwarz should complete");
        assert_eq!(out.solutions.len(), 2);
        assert!(out.residual.is_finite());
    }

    #[test]
    fn anderson_accelerated_schwarz_runs() {
        let subs = vec![mock_subdomain(0, 0), mock_subdomain(1, 0)];
        let ifaces = vec![
            InterfacePatch::new(0, vec![0], vec![0], 1, Complex64::new(0.0, 1.0)),
            InterfacePatch::new(1, vec![0], vec![0], 0, Complex64::new(0.0, 1.0)),
        ];
        let cfg = SchwarzConfig {
            anderson_depth: 3, tol: 1e-9, max_iter: 10,
            ..Default::default()
        };
        let out = schwarz_solve_full(&subs, &ifaces, &rem_parallel::NoComm, &cfg, None)
            .expect("anderson schwarz should complete");
        assert!(out.residual.is_finite());
    }
}
