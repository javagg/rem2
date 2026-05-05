//! 子域数据结构与本地 FEM 组装

use nalgebra::{DMatrix, DVector};
use num_complex::Complex64;
use rem_core::RemResult;
use rem_mesh::RemMesh;
use std::collections::HashMap;

/// 单个子域的数据
pub struct SubDomain {
    /// 子域编号（0-based）
    pub id: usize,
    /// 本子域内的体单元索引（在原始 RemMesh 中的下标）
    pub volume_elements: Vec<usize>,
    /// 本子域内的边界单元索引（在原始 RemMesh 中的下标）
    pub boundary_elements: Vec<usize>,
    /// 本子域节点（全局节点编号 → 本地编号 映射）
    pub global_to_local: HashMap<usize, usize>,
    /// 本地编号 → 全局编号
    pub local_to_global: Vec<usize>,
    /// 与其他子域共享的界面节点（本地编号）
    pub interface_nodes: Vec<usize>,
    /// 每个界面节点对应的相邻子域编号
    pub interface_neighbor: Vec<usize>,
}

impl SubDomain {
    /// 从分区结果构建子域
    pub fn build(
        id: usize,
        mesh: &RemMesh,
        partition: &[i32],
    ) -> Self {
        // 收集本子域内的体单元
        let volume_elements: Vec<usize> = mesh.volume_elements.iter()
            .enumerate()
            .filter(|(i, _)| partition[*i] as usize == id)
            .map(|(i, _)| i)
            .collect();

        // 收集本子域内涉及的所有节点
        let mut node_set = std::collections::BTreeSet::new();
        for &ei in &volume_elements {
            for &nid in &mesh.volume_elements[ei].node_ids {
                node_set.insert(nid);
            }
        }

        let local_to_global: Vec<usize> = node_set.into_iter().collect();
        let global_to_local: HashMap<usize, usize> = local_to_global.iter()
            .enumerate()
            .map(|(li, &gi)| (gi, li))
            .collect();

        // 收集本子域内的边界单元
        let boundary_elements: Vec<usize> = mesh.boundary_elements.iter()
            .enumerate()
            .filter(|(_, elem)| {
                elem.node_ids.iter().all(|nid| global_to_local.contains_key(nid))
            })
            .map(|(i, _)| i)
            .collect();

        SubDomain {
            id,
            volume_elements,
            boundary_elements,
            global_to_local,
            local_to_global,
            interface_nodes: Vec::new(),
            interface_neighbor: Vec::new(),
        }
    }

    /// 本子域 DOF 数量
    pub fn n_dof(&self) -> usize {
        self.local_to_global.len()
    }

    /// 组装本地刚度矩阵（P1 FEM 占位）
    ///
    /// 当前为骨架：返回单位矩阵。
    /// TODO：接入 rem-driven 的真实 FEM 装配。
    pub fn assemble_local_stiffness_skeleton(
        &self,
    ) -> RemResult<(DMatrix<Complex64>, DVector<Complex64>)> {
        let n = self.n_dof();
        let mat = DMatrix::identity(n, n).map(|x: f64| Complex64::new(x, 0.0));
        let rhs = DVector::zeros(n);
        Ok((mat, rhs))
    }

    /// 组装本地刚度矩阵（P1 FEM 占位）
    ///
    /// 当前为骨架：返回单位矩阵。
    /// TODO：接入 rem-driven 的真实 FEM 装配。
    pub fn assemble_local_stiffness(
        &self,
        _mesh: &RemMesh,
        _freq: f64,
    ) -> RemResult<(DMatrix<Complex64>, DVector<Complex64>)> {
        self.assemble_local_stiffness_skeleton()
    }
}
