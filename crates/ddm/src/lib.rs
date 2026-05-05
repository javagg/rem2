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
pub mod anderson;

use num_complex::Complex64;
use rem_config::{PalaceConfig, DdmSolverConfig};
use rem_core::RemResult;
use rem_mesh::RemMesh;
use rem_parallel::{NoComm, Comm};
use std::collections::{BTreeMap, BTreeSet};

use subdomain::SubDomain;
use interface::InterfacePatch;

/// Build interface patches from shared volume nodes across subdomains.
///
/// For each shared node, creates directed interface records (owner -> neighbor)
/// so Schwarz updates can apply per-neighbor Robin coupling later.
fn build_interfaces(
    subdomains: &mut [SubDomain],
    mesh: &RemMesh,
    partition: &[i32],
    robin_alpha: Complex64,
) -> Vec<InterfacePatch> {
    let n_sub = subdomains.len();
    let mut node_subdomains: Vec<BTreeSet<usize>> = vec![BTreeSet::new(); mesh.nodes.len()];

    for (ei, elem) in mesh.volume_elements.iter().enumerate() {
        let sid = partition.get(ei).copied().unwrap_or(0).max(0) as usize;
        if sid >= n_sub {
            continue;
        }
        for &nid in &elem.node_ids {
            if nid < node_subdomains.len() {
                node_subdomains[nid].insert(sid);
            }
        }
    }

    let mut node_pairs: BTreeMap<(usize, usize), Vec<usize>> = BTreeMap::new();
    for (nid, owners) in node_subdomains.iter().enumerate() {
        if owners.len() < 2 {
            continue;
        }
        let owners_vec: Vec<usize> = owners.iter().copied().collect();
        for &owner in &owners_vec {
            for &neighbor in &owners_vec {
                if owner != neighbor {
                    node_pairs.entry((owner, neighbor)).or_default().push(nid);
                }
            }
        }
    }

    let mut interfaces = Vec::new();

    // Deduplicate and construct InterfacePatch.
    for ((owner, neighbor), mut global_nodes) in node_pairs {
        if owner >= n_sub || neighbor >= n_sub {
            continue;
        }

        global_nodes.sort_unstable();
        global_nodes.dedup();

        let mut local_dofs = Vec::with_capacity(global_nodes.len());
        let mut owner_iface_pairs: BTreeSet<(usize, usize)> = BTreeSet::new();

        for &gid in &global_nodes {
            if let Some(&ldof) = subdomains[owner].global_to_local.get(&gid) {
                local_dofs.push(ldof);
                owner_iface_pairs.insert((ldof, neighbor));
            }
        }

        if local_dofs.is_empty() {
            continue;
        }

        // Keep SubDomain interface metadata in sync with generated patches.
        for (ldof, nbr) in owner_iface_pairs {
            subdomains[owner].interface_nodes.push(ldof);
            subdomains[owner].interface_neighbor.push(nbr);
        }

        interfaces.push(InterfacePatch::new(
            owner as i32,
            local_dofs,
            global_nodes,
            neighbor as i32,
            robin_alpha,
        ));
    }

    interfaces
}

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
    let mut subdomains: Vec<SubDomain> = (0..n_sub)
        .map(|id| SubDomain::build(id, mesh, &part))
        .collect();
    log::info!("Subdomains built: avg {} elements each",
        mesh.volume_elements.len() / n_sub.max(1));

    // 3. 识别子域界面 DOF（基于共享体节点）
    let omega = 2.0 * std::f64::consts::PI * ddm_cfg.freq_hz;
    let k0 = omega / rem_core::C0;
    let robin_alpha = num_complex::Complex64::new(0.0, k0 * ddm_cfg.robin_order.max(1) as f64);
    let interfaces = build_interfaces(&mut subdomains, mesh, &part, robin_alpha);
    log::info!("Interfaces: {} interface pairs, α = {:.3e}j (k₀={:.3e}, f={:.3e} Hz)",
        interfaces.len(), robin_alpha.im, k0, ddm_cfg.freq_hz);

    // 4. Schwarz 迭代求解（使用完整配置）
    let schwarz_cfg = schwarz::SchwarzConfig {
        tol: ddm_cfg.tolerance,
        max_iter: ddm_cfg.max_iter,
        robin_alpha,
        multiplicative: ddm_cfg.multiplicative,
        anderson_depth: ddm_cfg.anderson_depth,
        freq_hz: ddm_cfg.freq_hz,
        eps_r: ddm_cfg.eps_r,
        mu_r: ddm_cfg.mu_r,
    };
    let schwarz_result = schwarz::schwarz_solve_full(
        &subdomains,
        &interfaces,
        comm,
        &schwarz_cfg,
        Some(mesh),
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
    use rem_mesh::{Element, ElementKind, Node, RemMesh};
    use std::collections::HashMap;

    #[test]
    fn test_ddm_config_defaults() {
        let cfg = DdmSolverConfig {
            num_subdomains: 4,
            method: "Schwarz".to_string(),
            robin_order: 1,
            tolerance: 1e-6,
            max_iter: 100,
            partition_type: "Dual".to_string(),
            ..Default::default()
        };
        assert_eq!(cfg.num_subdomains, 4);
        assert_eq!(cfg.method, "Schwarz");
        assert!(cfg.tolerance < 1e-4);
    }

    #[test]
    fn test_build_interfaces_minimal_shared_node() {
        let mesh = RemMesh {
            nodes: vec![
                Node { id: 1, x: 0.0, y: 0.0, z: 0.0 },
                Node { id: 2, x: 1.0, y: 0.0, z: 0.0 },
                Node { id: 3, x: 0.0, y: 1.0, z: 0.0 },
                Node { id: 4, x: 0.0, y: 0.0, z: 1.0 },
                Node { id: 5, x: 2.0, y: 0.0, z: 0.0 },
                Node { id: 6, x: 0.0, y: 2.0, z: 0.0 },
                Node { id: 7, x: 0.0, y: 0.0, z: 2.0 },
            ],
            volume_elements: vec![
                Element {
                    id: 1,
                    kind: ElementKind::Tet4,
                    tag: 1,
                    node_ids: vec![0, 1, 2, 3],
                    rank: 0,
                },
                Element {
                    id: 2,
                    kind: ElementKind::Tet4,
                    tag: 1,
                    node_ids: vec![0, 4, 5, 6],
                    rank: 0,
                },
            ],
            boundary_elements: vec![],
            domain_tags: HashMap::new(),
            boundary_tags: HashMap::new(),
            dim: 3,
            rank: 0,
            size: 1,
        };

        let partition = vec![0_i32, 1_i32];
        let mut subdomains: Vec<SubDomain> = (0..2)
            .map(|id| SubDomain::build(id, &mesh, &partition))
            .collect();
        let robin_alpha = Complex64::new(0.0, 1.0);

        let interfaces = build_interfaces(&mut subdomains, &mesh, &partition, robin_alpha);

        assert_eq!(interfaces.len(), 2, "expected two directed interface patches");
        assert_eq!(interfaces.iter().filter(|p| p.owner_rank == 0 && p.neighbor_rank == 1).count(), 1);
        assert_eq!(interfaces.iter().filter(|p| p.owner_rank == 1 && p.neighbor_rank == 0).count(), 1);

        for patch in &interfaces {
            assert_eq!(patch.global_node_ids, vec![0]);
            assert_eq!(patch.local_dofs.len(), 1);
            assert_eq!(patch.robin_alpha, robin_alpha);
        }

        assert_eq!(subdomains[0].interface_nodes.len(), 1);
        assert_eq!(subdomains[1].interface_nodes.len(), 1);
        assert_eq!(subdomains[0].interface_neighbor, vec![1]);
        assert_eq!(subdomains[1].interface_neighbor, vec![0]);
    }
}
