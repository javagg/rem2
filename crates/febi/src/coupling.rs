//! FEM-BI 耦合矩阵组装
//!
//! 组装混合系统（Schur 补近似）：
//!   A_eff = A_fem + B^T · Z_bi · B
//!
//! 其中：
//!   A_fem = K_fem − k²M_fem  （FEM 体域刚度矩阵，n_vol × n_vol）
//!   Z_bi                      （BI 阻抗矩阵，n_rwg × n_rwg）
//!   B                         （耦合矩阵，n_rwg × n_vol）
//!
//! 耦合矩阵 B 的构造：
//!   对每条 RWG 边 e（正面 f+，负面 f−）：
//!     B[e, global_node] = +area_f+ / 3  （f+ 的三个节点）
//!     B[e, global_node] = −area_f− / 3  （f− 的三个节点）

use rem_config::FeBiSolverConfig;
use rem_core::{RemResult, CsrMatrix};
use rem_mesh::{RemMesh, BoundaryTag};
use rem_mom::surface_mesh::SurfaceMesh;
use rem_electrostatic::assemble::assemble_stiffness;
use rem_eigenmode::assemble_mass::assemble_mass;
use nalgebra::{DMatrix, DVector};
use num_complex::Complex64;

/// 混合系统表示
pub struct FebiSystem {
    /// 有效系统矩阵 A_eff = A_fem + B^T·Z_bi·B（n_vol × n_vol）
    pub z_bi: DMatrix<Complex64>,
    /// 右端向量（n_vol）
    pub rhs: DVector<Complex64>,
    /// DOF 数量（= n_vol）
    pub n_dof: usize,
}

/// 组装 FE-BI 混合系统
pub fn assemble_febi_system(
    febi_cfg: &FeBiSolverConfig,
    mesh: &RemMesh,
    surf: &SurfaceMesh,
    z_bi: &DMatrix<Complex64>,
    freq: f64,
) -> RemResult<FebiSystem> {
    let omega = 2.0 * std::f64::consts::PI * freq;
    let c0 = 299_792_458.0;
    let k0 = omega / c0;
    let eps_r = febi_cfg.exterior_eps_r;
    let mu_r = febi_cfg.exterior_mu_r;
    let k = k0 * (eps_r * mu_r).sqrt();
    let k2 = k * k;
    let n_vol = mesh.n_nodes();
    let n_rwg = z_bi.nrows();

    // 1. FEM 系统矩阵 A_fem = K − k²M
    let eps_fn = |_tag: u32| eps_r;
    let k_csr = assemble_stiffness(mesh, eps_fn)?.to_csr();
    let m_csr = assemble_mass(mesh, eps_fn)?.to_csr();
    let mut a_fem = csr_to_complex_dense(&k_csr, n_vol);
    let m_dense = csr_to_complex_dense(&m_csr, n_vol);
    for i in 0..n_vol {
        for j in 0..n_vol {
            a_fem[(i, j)] -= Complex64::new(k2, 0.0) * m_dense[(i, j)];
        }
    }

    // 2. 耦合矩阵 B（n_rwg × n_vol）
    let b = build_coupling_matrix(surf, n_vol, n_rwg);

    // 3. A_eff = A_fem + B^T · Z_bi · B
    if n_rwg > 0 {
        let bt_zbi = b.transpose() * z_bi;   // n_vol × n_rwg
        let correction = bt_zbi * &b;         // n_vol × n_vol
        for i in 0..n_vol {
            for j in 0..n_vol {
                a_fem[(i, j)] += correction[(i, j)];
            }
        }
    }

    // 4. RHS：端口激励（在集总端口节点上设 φ = 1）
    let rhs = build_febi_rhs(febi_cfg, mesh, n_vol);

    log::info!("FE-BI system: n_vol={}, n_rwg={}", n_vol, n_rwg);

    Ok(FebiSystem { z_bi: a_fem, rhs, n_dof: n_vol })
}

// ── 辅助函数 ──────────────────────────────────────────────────────────────────

fn csr_to_complex_dense(csr: &CsrMatrix, n: usize) -> DMatrix<Complex64> {
    let mut m = DMatrix::<Complex64>::zeros(n, n);
    for i in 0..n {
        for k in csr.row_ptr[i]..csr.row_ptr[i + 1] {
            let j = csr.col_idx[k];
            m[(i, j)] = Complex64::new(csr.values[k], 0.0);
        }
    }
    m
}

fn build_coupling_matrix(surf: &SurfaceMesh, n_vol: usize, n_rwg: usize) -> DMatrix<Complex64> {
    let mut b = DMatrix::<Complex64>::zeros(n_rwg, n_vol);
    for (e_idx, edge) in surf.edges.iter().enumerate() {
        let f_plus  = &surf.faces[edge.plus_face];
        let f_minus = &surf.faces[edge.minus_face];
        for &ln in &f_plus.nodes {
            let gn = surf.global_node_ids[ln];
            if gn < n_vol { b[(e_idx, gn)] += Complex64::new(f_plus.area / 3.0, 0.0); }
        }
        for &ln in &f_minus.nodes {
            let gn = surf.global_node_ids[ln];
            if gn < n_vol { b[(e_idx, gn)] -= Complex64::new(f_minus.area / 3.0, 0.0); }
        }
    }
    b
}

fn build_febi_rhs(febi_cfg: &FeBiSolverConfig, mesh: &RemMesh, n_dof: usize) -> DVector<Complex64> {
    let mut rhs = DVector::zeros(n_dof);
    for port in &febi_cfg.ports {
        for belem in &mesh.boundary_elements {
            if let Some(BoundaryTag::LumpedPort { index, .. }) = mesh.boundary_tags.get(&belem.tag) {
                if *index == port.index {
                    for &nid in &belem.node_ids {
                        if nid < n_dof { rhs[nid] = Complex64::new(1.0, 0.0); }
                    }
                }
            }
        }
    }
    rhs
}
