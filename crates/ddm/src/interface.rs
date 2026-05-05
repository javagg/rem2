//! 界面自由度管理与 Robin 传输边界条件
//!
//! 在子域 Ωᵢ 和 Ωⱼ 的共享界面 Γᵢⱼ 上施加 Robin TBC：
//!   (nᵢ × H + α·E)|_Ωᵢ = −(nⱼ × H + α·E)|_Ωⱼ
//!
//! 其中 α = jk（一阶最优 OSRC 条件）

use num_complex::Complex64;

/// 一个子域与相邻子域的界面描述
#[derive(Debug, Clone)]
pub struct InterfacePatch {
    /// 本界面补丁所属子域编号（本地子域 ID）
    pub owner_rank: i32,
    /// 本子域中属于该界面的自由度（局部索引）
    pub local_dofs: Vec<usize>,
    /// 对应的全局节点 ID（用于与邻域匹配）
    pub global_node_ids: Vec<usize>,
    /// 相邻子域的 MPI rank
    pub neighbor_rank: i32,
    /// Robin 条件系数 α = jk₀（波数相关）
    pub robin_alpha: Complex64,
}

impl InterfacePatch {
    pub fn new(
        owner_rank: i32,
        local_dofs: Vec<usize>,
        global_node_ids: Vec<usize>,
        neighbor_rank: i32,
        robin_alpha: Complex64,
    ) -> Self {
        Self { owner_rank, local_dofs, global_node_ids, neighbor_rank, robin_alpha }
    }

    pub fn n_dofs(&self) -> usize { self.local_dofs.len() }
}

/// Robin 条件贡献：将界面项添加到子域矩阵对角线
///
/// Z[i,i] += α · area_i（简化标量形式，完整版需积分 RWG 基函数）
pub fn apply_robin_to_diagonal(
    diag: &mut [Complex64],
    patch: &InterfacePatch,
    face_areas: &[f64],
) {
    for (&ldof, &area) in patch.local_dofs.iter().zip(face_areas.iter()) {
        if ldof < diag.len() {
            diag[ldof] += patch.robin_alpha * Complex64::new(area, 0.0);
        }
    }
}

/// 从相邻子域接收到的界面场值（用于更新 RHS）
#[derive(Debug, Clone)]
pub struct InterfaceExchange {
    /// 来自相邻子域的界面 E 切向分量
    pub incoming_e: Vec<Complex64>,
    /// 来自相邻子域的界面 H 切向分量
    pub incoming_h: Vec<Complex64>,
}

impl InterfaceExchange {
    pub fn zeros(n: usize) -> Self {
        Self {
            incoming_e: vec![Complex64::ZERO; n],
            incoming_h: vec![Complex64::ZERO; n],
        }
    }

    /// 计算 Robin RHS 贡献：f_robin = 2α·E_in + 2·(n×H_in)
    pub fn robin_rhs_contribution(&self, alpha: Complex64) -> Vec<Complex64> {
        self.incoming_e.iter().zip(self.incoming_h.iter())
            .map(|(&e, &h)| {
                let two = Complex64::new(2.0, 0.0);
                two * alpha * e + two * h
            })
            .collect()
    }
}
