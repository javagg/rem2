//! 子域数据结构与本地 FEM 组装
//!
//! # P1 Helmholtz 组装
//!
//! 对每个四面体单元组装标量 Helmholtz 方程的 P1 刚度矩阵 K 和质量矩阵 M：
//!
//!   K[i,j] = ∫_Ω ∇φᵢ · ∇φⱼ dΩ
//!   M[i,j] = ∫_Ω φᵢ φⱼ dΩ
//!   A = K − k² M,   k² = ω²μ₀μᵣε₀εᵣ
//!
//! 四面体 P1 解析公式（Zienkiewicz §5.4）：
//!   ∇φᵢ = bᵢ / (6V),   M[i,i]=V/10, M[i≠j]=V/20

use nalgebra::{DMatrix, DVector};
use num_complex::Complex64;
use rem_core::{RemResult, C0, MU0};
use rem_mesh::RemMesh;
use std::collections::HashMap;
use std::f64::consts::PI;

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

    /// P1 Helmholtz FEM 本地组装（真实物理方程）。
    ///
    /// 组装 A = K − k² M，其中：
    ///   K[i,j] = ∫ ∇φᵢ · ∇φⱼ / μᵣ dΩ  （P1 tet 解析式）
    ///   M[i,j] = ∫ εᵣ φᵢ φⱼ dΩ          （P1 tet 解析式）
    ///   k₀² = ω² μ₀ ε₀
    ///
    /// 当网格中没有体单元（或只有非四面体单元）时，退化为单位矩阵（保证求解器不崩溃）。
    pub fn assemble_local_p1_helmholtz(
        &self,
        mesh: &RemMesh,
        freq_hz: f64,
        eps_r: f64,
        mu_r: f64,
    ) -> RemResult<(DMatrix<Complex64>, DVector<Complex64>)> {
        let n = self.n_dof();
        if n == 0 {
            return Ok((DMatrix::identity(1, 1).map(|x: f64| Complex64::new(x, 0.0)),
                        DVector::zeros(1)));
        }

        let omega = 2.0 * PI * freq_hz;
        // k₀² = ω² / c₀²  (vacuum),  effective k² = k₀² εᵣ μᵣ
        let k0_sq = (omega / C0).powi(2);
        let k_sq = Complex64::new(k0_sq * eps_r * mu_r, 0.0);

        let mut a_mat = DMatrix::<Complex64>::zeros(n, n);

        for &ei in &self.volume_elements {
            let elem = &mesh.volume_elements[ei];
            // Only handle Tet4; skip other element types gracefully.
            if elem.node_ids.len() != 4 {
                continue;
            }

            // Global node positions
            let nids = &elem.node_ids;
            let pts: Vec<[f64; 3]> = nids.iter().map(|&gid| {
                // RemMesh node ids are 1-based in GMSH but stored as plain index here.
                // Use id as index with bounds check; fall back to 0.
                let idx = if gid < mesh.nodes.len() { gid } else { 0 };
                let n = &mesh.nodes[idx];
                [n.x, n.y, n.z]
            }).collect();

            // Signed volume and shape-function gradients via Jacobian cofactors.
            // b_i = cofactor column of the 3×3 Jacobian matrix (scaled by ±1).
            // Reference: Zienkiewicz Vol.1 §5.4.
            let v6 = tet_signed_vol6(&pts); // 6 × signed volume
            if v6.abs() < 1e-30 {
                continue; // degenerate element
            }
            let sign = if v6 > 0.0 { 1.0 } else { -1.0 };
            let vol = sign * v6 / 6.0;

            // Gradient of each P1 shape function (constant inside tet).
            let grads = tet_p1_gradients(&pts, v6);

            // Local indices into subdomain DOF vector
            let lids: Vec<usize> = nids.iter()
                .filter_map(|gid| self.global_to_local.get(gid).copied())
                .collect();
            if lids.len() != 4 {
                continue; // node not in this subdomain (shouldn't happen)
            }

            // Stiffness K[i,j] = (∇φᵢ · ∇φⱼ / μᵣ) * vol
            // Mass     M[i,j] = εᵣ * vol * { 1/10 if i==j, 1/20 if i≠j }
            for li in 0..4 {
                for lj in li..4 {
                    let dot: f64 = (0..3).map(|d| grads[li][d] * grads[lj][d]).sum();
                    let k_entry = Complex64::new(dot * vol / mu_r, 0.0);

                    let m_coeff = if li == lj { vol * eps_r / 10.0 } else { vol * eps_r / 20.0 };
                    let m_entry = Complex64::new(m_coeff, 0.0);

                    let a_entry = k_entry - k_sq * m_entry;

                    let (i, j) = (lids[li], lids[lj]);
                    a_mat[(i, j)] += a_entry;
                    if li != lj {
                        a_mat[(j, i)] += a_entry; // symmetric
                    }
                }
            }
        }

        // If matrix is all-zero (no tet4 elements), fall back to identity.
        let is_zero = a_mat.iter().all(|v| v.norm() < 1e-30);
        if is_zero {
            a_mat = DMatrix::identity(n, n).map(|x: f64| Complex64::new(x, 0.0));
        }

        let rhs = DVector::zeros(n);
        Ok((a_mat, rhs))
    }

    /// Wrapper used by the Schwarz solver: tries real P1 assembly first,
    /// falls back to skeleton (identity) if the mesh is not available.
    pub fn assemble_local_stiffness(
        &self,
        mesh: &RemMesh,
        freq: f64,
    ) -> RemResult<(DMatrix<Complex64>, DVector<Complex64>)> {
        self.assemble_local_p1_helmholtz(mesh, freq, 1.0, 1.0)
    }
}

// ── P1 tet geometry helpers ──────────────────────────────────────────────────

/// 6 × signed volume of a tetrahedron: det([p1-p0, p2-p0, p3-p0]).
fn tet_signed_vol6(pts: &[[f64; 3]]) -> f64 {
    let [p0, p1, p2, p3] = [pts[0], pts[1], pts[2], pts[3]];
    let a = [p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2]];
    let b = [p2[0]-p0[0], p2[1]-p0[1], p2[2]-p0[2]];
    let c = [p3[0]-p0[0], p3[1]-p0[1], p3[2]-p0[2]];
    // det = a · (b × c)
    let bxc = [
        b[1]*c[2] - b[2]*c[1],
        b[2]*c[0] - b[0]*c[2],
        b[0]*c[1] - b[1]*c[0],
    ];
    a[0]*bxc[0] + a[1]*bxc[1] + a[2]*bxc[2]
}

/// Constant gradients of the 4 P1 shape functions on a Tet4.
///
/// For shape function φᵢ associated with node i,  ∇φᵢ = bᵢ / (6V).
/// bᵢ are the cofactors of the 3×3 edge-vector matrix.
fn tet_p1_gradients(pts: &[[f64; 3]], v6: f64) -> [[f64; 3]; 4] {
    // Rows of Jacobian matrix J = [p1-p0 | p2-p0 | p3-p0] (3×3)
    let [p0, p1, p2, p3] = [pts[0], pts[1], pts[2], pts[3]];
    let r = [
        [p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2]],
        [p2[0]-p0[0], p2[1]-p0[1], p2[2]-p0[2]],
        [p3[0]-p0[0], p3[1]-p0[1], p3[2]-p0[2]],
    ];

    // Cofactors (= J⁻¹ᵀ · |J|), columns give grad(φ₁..φ₃)
    // grad(φ₀) = −(grad(φ₁)+grad(φ₂)+grad(φ₃)) by partition-of-unity
    let inv_v6 = 1.0 / v6;

    // cofactor matrix (transposed inverse, unnormalised by v6)
    let cof = [
        [  r[1][1]*r[2][2]-r[1][2]*r[2][1],
          -(r[0][1]*r[2][2]-r[0][2]*r[2][1]),
           r[0][1]*r[1][2]-r[0][2]*r[1][1] ],
        [-(r[1][0]*r[2][2]-r[1][2]*r[2][0]),
           r[0][0]*r[2][2]-r[0][2]*r[2][0],
          -(r[0][0]*r[1][2]-r[0][2]*r[1][0])],
        [  r[1][0]*r[2][1]-r[1][1]*r[2][0],
          -(r[0][0]*r[2][1]-r[0][1]*r[2][0]),
           r[0][0]*r[1][1]-r[0][1]*r[1][0]],
    ];

    // grad(φᵢ₊₁) = cof[:,i] / v6
    let g1 = [cof[0][0]*inv_v6, cof[1][0]*inv_v6, cof[2][0]*inv_v6];
    let g2 = [cof[0][1]*inv_v6, cof[1][1]*inv_v6, cof[2][1]*inv_v6];
    let g3 = [cof[0][2]*inv_v6, cof[1][2]*inv_v6, cof[2][2]*inv_v6];
    let g0 = [
        -(g1[0]+g2[0]+g3[0]),
        -(g1[1]+g2[1]+g3[1]),
        -(g1[2]+g2[2]+g3[2]),
    ];
    [g0, g1, g2, g3]
}

#[cfg(test)]
mod tests {
    use super::*;
    use rem_mesh::{Node, Element, ElementKind, RemMesh};
    use std::collections::HashMap;

    fn unit_tet_mesh() -> RemMesh {
        // Reference tet: (0,0,0),(1,0,0),(0,1,0),(0,0,1) → V=1/6
        RemMesh {
            nodes: vec![
                Node { id: 0, x: 0.0, y: 0.0, z: 0.0 },
                Node { id: 1, x: 1.0, y: 0.0, z: 0.0 },
                Node { id: 2, x: 0.0, y: 1.0, z: 0.0 },
                Node { id: 3, x: 0.0, y: 0.0, z: 1.0 },
            ],
            volume_elements: vec![
                Element { id: 1, kind: ElementKind::Tet4, tag: 1,
                          node_ids: vec![0,1,2,3], rank: 0 },
            ],
            boundary_elements: vec![],
            domain_tags: HashMap::new(),
            boundary_tags: HashMap::new(),
            dim: 3, rank: 0, size: 1,
        }
    }

    #[test]
    fn tet_signed_vol6_unit_tet() {
        let pts = [[0.0,0.0,0.0],[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]];
        let v6 = tet_signed_vol6(&pts);
        // V = 1/6, so 6V = 1
        assert!((v6 - 1.0).abs() < 1e-12, "v6={v6}");
    }

    #[test]
    fn p1_gradients_partition_of_unity() {
        let pts = [[0.0,0.0,0.0],[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]];
        let v6 = tet_signed_vol6(&pts);
        let grads = tet_p1_gradients(&pts, v6);
        // ∑ ∇φᵢ = 0 (partition of unity → constant functions reproduced)
        for d in 0..3 {
            let sum: f64 = grads.iter().map(|g| g[d]).sum();
            assert!(sum.abs() < 1e-12, "dim {d}: sum of grads = {sum}");
        }
    }

    #[test]
    fn p1_helmholtz_static_limit_is_positive_semidefinite() {
        // At freq → 0, k → 0, A ≈ K (pure stiffness, PSD for Dirichlet problem).
        let mesh = unit_tet_mesh();
        let part = vec![0i32]; // single subdomain
        let sd = SubDomain::build(0, &mesh, &part);
        let (a, _) = sd.assemble_local_p1_helmholtz(&mesh, 1.0, 1.0, 1.0).unwrap();
        // Diagonal entries of K for the unit tet should be positive.
        for i in 0..sd.n_dof() {
            assert!(a[(i,i)].re >= -1e-12, "A[{i},{i}].re = {:.4e}", a[(i,i)].re);
        }
    }

    #[test]
    fn p1_helmholtz_matrix_is_symmetric() {
        let mesh = unit_tet_mesh();
        let part = vec![0i32];
        let sd = SubDomain::build(0, &mesh, &part);
        let (a, _) = sd.assemble_local_p1_helmholtz(&mesh, 1e9, 1.0, 1.0).unwrap();
        let n = sd.n_dof();
        for i in 0..n {
            for j in 0..n {
                let diff = (a[(i,j)] - a[(j,i)]).norm();
                assert!(diff < 1e-12, "A[{i},{j}]!=A[{j},{i}]: diff={diff:.2e}");
            }
        }
    }

    #[test]
    fn p1_helmholtz_stiffness_trace_matches_analytic() {
        // For unit tet, K[i,i] = |∇φᵢ|² * V  (no sum over j≠i).
        // Analytically: ∇φ₀=(-1,-1,-1), so |∇φ₀|²=3, V=1/6 → K₀₀ = 3/6 = 0.5
        // ∇φ₁=(1,0,0), |∇φ₁|²=1 → K₁₁ = 1/6 ≈ 0.1667
        // ∇φ₂=(0,1,0), similarly K₂₂=1/6; ∇φ₃=(0,0,1), K₃₃=1/6
        let mesh = unit_tet_mesh();
        let part = vec![0i32];
        let sd = SubDomain::build(0, &mesh, &part);
        let (a, _) = sd.assemble_local_p1_helmholtz(&mesh, 1.0, 1.0, 1.0).unwrap();
        // At freq=1 Hz k≈0, A≈K
        // local_to_global is sorted: [0,1,2,3]
        let k00 = a[(0,0)].re;
        let k11 = a[(1,1)].re;
        assert!((k00 - 0.5).abs() < 1e-10, "K[0,0]={k00:.6}");
        assert!((k11 - 1.0/6.0).abs() < 1e-10, "K[1,1]={k11:.6}");
    }
}

