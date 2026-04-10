//! rem-ddm — Domain Decomposition Method (DDM) for large-scale EM simulation
//!
//! 将大计算域分割为多个子域并行/迭代求解，突破单机内存限制。
//!
//! # 方法
//! - 基于 Robin 传输条件（TBC）的 Schwarz 迭代
//! - METIS 子域划分（复用 vendor/rmetis）
//! - MPI 通信（复用 rem-parallel 的 Comm trait）
//!
//! # 配置示例
//! ```json
//! {
//!   "Problem": { "Type": "Driven" },
//!   "Solver": {
//!     "DDM": {
//!       "NumSubdomains": 4,
//!       "Method": "Schwarz",
//!       "Tolerance": 1e-6,
//!       "MaxIter": 100,
//!       "RobinOrder": 1
//!     }
//!   }
//! }
//! ```

pub mod partition;
pub mod subdomain;
pub mod interface;
pub mod schwarz;
pub mod postprocess;

use num_complex::Complex64;
use rem_config::{PalaceConfig, DdmSolverConfig};
use rem_core::RemResult;
use rem_mesh::RemMesh;
use rem_parallel::{NoComm, Comm};

use subdomain::SubDomain;
use interface::InterfacePatch;

/// DDM 求解结果
#[derive(Debug, Clone)]
pub struct DdmResult {
    /// 每个子域的解向量（体 DOF）
    pub subdomain_solutions: Vec<Vec<Complex64>>,
    /// 迭代次数
    pub iterations: usize,
    /// 最终残差
    pub residual: f64,
}

/// CLI 入口：DDM 作为求解器加速器（Problem.Type = "Driven" + Solver.DDM）
pub fn run(config: &PalaceConfig) -> RemResult<()> {
    let ddm_cfg = config.solver.ddm.as_ref()
        .ok_or_else(|| rem_core::RemError::Config(
            "DDM solver requires a Solver.DDM section".to_string()
        ))?;

    let mesh = rem_mesh::load_mesh(config, &NoComm)?;
    run_with_mesh(config, ddm_cfg, &mesh, &NoComm).map(|_| ())
}

/// 在已加载网格上运行 DDM 求解（供测试/WASM 调用）
pub fn run_with_mesh(
    _config: &PalaceConfig,
    ddm_cfg: &DdmSolverConfig,
    mesh: &RemMesh,
    comm: &impl Comm,
) -> RemResult<DdmResult> {
    let n_sub = ddm_cfg.num_subdomains.max(1);
    log::info!("DDM solver start — {} subdomains, method = {}",
        n_sub, ddm_cfg.method);

    // 1. METIS 子域划分
    let part = partition::partition_mesh(mesh, n_sub)?;
    log::info!("Mesh partitioned: {} elements → {} subdomains",
        mesh.volume_elements.len(), n_sub);

    // 2. 构建子域数据结构
    let subdomains: Vec<SubDomain> = (0..n_sub)
        .map(|id| SubDomain::build(id, mesh, &part))
        .collect();
    log::info!("Subdomains built: avg {} elements each",
        mesh.volume_elements.len() / n_sub.max(1));

    // 3. 识别子域界面 DOF（骨架：空列表）
    // TODO: 实际界面检测需要子域间共享节点识别
    let interfaces: Vec<InterfacePatch> = Vec::new();
    log::info!("Interfaces: {} interface pairs", interfaces.len());

    // 4. Schwarz 迭代求解
    let schwarz_result = schwarz::schwarz_solve(
        &subdomains,
        &interfaces,
        comm,
        ddm_cfg.tolerance,
        ddm_cfg.max_iter,
    )?;

    log::info!("DDM converged in {} iterations, residual = {:.3e}",
        schwarz_result.iterations, schwarz_result.residual);

    Ok(DdmResult {
        subdomain_solutions: schwarz_result.solutions.into_iter()
            .map(|v| v.iter().copied().collect())
            .collect(),
        iterations: schwarz_result.iterations,
        residual: schwarz_result.residual,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ddm_config_defaults() {
        let cfg = DdmSolverConfig {
            num_subdomains: 4,
            method: "Schwarz".to_string(),
            robin_order: 1,
            tolerance: 1e-6,
            max_iter: 100,
            partition_type: "Dual".to_string(),
        };
        assert_eq!(cfg.num_subdomains, 4);
        assert_eq!(cfg.method, "Schwarz");
        assert!(cfg.tolerance < 1e-4);
    }
}
