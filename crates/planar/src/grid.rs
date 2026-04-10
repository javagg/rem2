/// 平面矩形网格（RWG 基函数用）
///
/// 网格覆盖 x ∈ [0, Lx]，y ∈ [0, Ly]，
/// 按 nx×ny 均匀分割为矩形单元，每对相邻单元共享一条内边，
/// 对应一个 RWG 基函数。

use nalgebra::Vector3;

pub type Point3 = Vector3<f64>;

/// 单个三角形面元（RWG 用半-基函数对应的三角形）
#[derive(Debug, Clone)]
pub struct Triangle {
    pub v: [Point3; 3], // 顶点
}

/// RWG 基函数：由两个共享公共边的三角形组成
#[derive(Debug, Clone)]
pub struct RwgBasis {
    /// 正三角形（+侧）
    pub tri_p: Triangle,
    /// 负三角形（-侧）
    pub tri_m: Triangle,
    /// 公共边长度
    pub edge_len: f64,
    /// 公共边中点
    pub edge_mid: Point3,
    /// 自由顶点 p（+侧，不在公共边上）
    pub free_p: Point3,
    /// 自由顶点 m（-侧，不在公共边上）
    pub free_m: Point3,
}

/// 在 z=0 平面生成均匀矩形网格，返回 RWG 基函数列表
///
/// # 参数
/// - `lx`: x 方向长度（单位：m）
/// - `ly`: y 方向长度（单位：m）
/// - `nx`: x 方向分段数
/// - `ny`: y 方向分段数
///
/// 矩形单元按对角线分成两个三角形（右下对角）。
/// 内部共享边产生 RWG 基函数；边界边不产生基函数（PEC 边界）。
pub fn build_planar_rwg(lx: f64, ly: f64, nx: usize, ny: usize) -> Vec<RwgBasis> {
    let dx = lx / nx as f64;
    let dy = ly / ny as f64;

    // 顶点网格坐标 (nx+1)×(ny+1)
    let vert = |i: usize, j: usize| -> Point3 {
        Vector3::new(i as f64 * dx, j as f64 * dy, 0.0)
    };

    // 每个矩形单元 (i,j) 被对角线分成两个三角形：
    //   T0: (i,j), (i+1,j), (i,j+1)
    //   T1: (i+1,j), (i+1,j+1), (i,j+1)
    // 公共斜边：(i+1,j) -- (i,j+1)
    //
    // 水平内边：(i,j+1)--(i+1,j+1)  属于单元 (i,j) 的 T1 和 单元 (i,j+1) 的 T0
    // 垂直内边：(i+1,j)--(i+1,j+1)  属于单元 (i,j) 的 T1 和 单元 (i+1,j) 的 T0

    let mut bases = Vec::new();

    // 斜向内边（每个矩形内部）：不跨单元，为矩形内部公共边
    // 每个矩形 (i,j) 的斜边 (i+1,j)--(i,j+1) 是内边
    for j in 0..ny {
        for i in 0..nx {
            // T0 自由顶点 = (i,j)
            // T1 自由顶点 = (i+1,j+1)
            let v00 = vert(i, j);
            let v10 = vert(i + 1, j);
            let v01 = vert(i, j + 1);
            let v11 = vert(i + 1, j + 1);

            let edge_len = (v10 - v01).norm();
            let edge_mid = (v10 + v01) * 0.5;

            bases.push(RwgBasis {
                tri_p: Triangle { v: [v10, v01, v00] }, // +侧，自由顶点 v00
                tri_m: Triangle { v: [v01, v10, v11] }, // -侧，自由顶点 v11
                edge_len,
                edge_mid,
                free_p: v00,
                free_m: v11,
            });
        }
    }

    // 水平内边：(i,j+1)--(i+1,j+1)，跨 (i,j) 的 T1 和 (i,j+1) 的 T0
    for j in 0..ny - 1 {
        for i in 0..nx {
            let v01 = vert(i, j + 1);
            let v11 = vert(i + 1, j + 1);
            let v10 = vert(i + 1, j);
            let v02 = vert(i, j + 2);

            let edge_len = (v11 - v01).norm();
            let edge_mid = (v01 + v11) * 0.5;

            // 下方三角形 T1 of (i,j): (v10,v11,v01) 自由顶点 v10
            // 上方三角形 T0 of (i,j+1): (v01,v11,v02)... 需重新确定方向
            // T0 of (i,j+1): (v01,v11,v02) → 自由顶点 v02
            bases.push(RwgBasis {
                tri_p: Triangle { v: [v01, v11, v10] }, // +侧，自由顶点 v10
                tri_m: Triangle { v: [v11, v01, v02] }, // -侧，自由顶点 v02
                edge_len,
                edge_mid,
                free_p: v10,
                free_m: v02,
            });
        }
    }

    // 垂直内边：(i+1,j)--(i+1,j+1)，跨 (i,j) 的 T1 和 (i+1,j) 的 T0
    for j in 0..ny {
        for i in 0..nx - 1 {
            let v10 = vert(i + 1, j);
            let v11 = vert(i + 1, j + 1);
            let v01 = vert(i, j + 1);
            let v20 = vert(i + 2, j);

            let edge_len = (v11 - v10).norm();
            let edge_mid = (v10 + v11) * 0.5;

            // 左侧 T1 of (i,j): (v10,v11,v01) 自由顶点 v01
            // 右侧 T0 of (i+1,j): (v20,v10,v11)... 自由顶点 = (i+2,j+1)? 不对
            // T0 of (i+1,j): vertices (i+1,j),(i+2,j),(i+1,j+1) 自由顶点 (i+2,j)=v20
            bases.push(RwgBasis {
                tri_p: Triangle { v: [v10, v11, v01] }, // +侧，自由顶点 v01
                tri_m: Triangle { v: [v11, v10, v20] }, // -侧，自由顶点 v20（实为(i+2,j)）
                edge_len,
                edge_mid,
                free_p: v01,
                free_m: v20,
            });
        }
    }

    bases
}

// ──────────────────────────────────────────────────────────────────────────────
// PlanarGrid：impedance.rs / solver.rs 所需的平面网格抽象
// ──────────────────────────────────────────────────────────────────────────────

/// 平面网格中的一条 "边"（对应一个标量基函数，简化 RWG 为矩形像素边）
#[derive(Debug, Clone)]
pub struct PlanarEdge {
    /// 边中点坐标 [x, y, z]
    pub center: [f64; 3],
    /// 切向方向（单位向量）[tx, ty, 0]
    pub tangent: [f64; 3],
    /// 边长度（m）
    pub length: f64,
    /// 散度系数 ρ = +1 / -1（RWG ±侧）
    pub div_rho: f64,
}

/// 均匀矩形平面网格（整合 RWG 基函数边表 + 网格几何）
#[derive(Debug, Clone)]
pub struct PlanarGrid {
    /// 所有基函数边列表
    pub edges: Vec<PlanarEdge>,
    /// x 方向分段数
    pub nx: usize,
    /// y 方向分段数
    pub ny: usize,
    /// x 方向单元尺寸（m）
    pub dx: f64,
    /// y 方向单元尺寸（m）
    pub dy: f64,
}

impl PlanarGrid {
    /// 从几何参数构造均匀平面网格
    ///
    /// 每个内部水平/垂直边产生两个 PlanarEdge（+/- 两侧）
    pub fn new(lx: f64, ly: f64, nx: usize, ny: usize) -> Self {
        let dx = lx / nx as f64;
        let dy = ly / ny as f64;
        let mut edges = Vec::new();

        // 水平边（切向 x）
        for j in 0..=ny {
            for i in 0..nx {
                let cx = (i as f64 + 0.5) * dx;
                let cy = j as f64 * dy;
                let is_boundary = j == 0 || j == ny;
                // 内部水平边产生两个像素基函数（±）
                if !is_boundary {
                    edges.push(PlanarEdge {
                        center: [cx, cy, 0.0],
                        tangent: [1.0, 0.0, 0.0],
                        length: dx,
                        div_rho: 1.0,
                    });
                    edges.push(PlanarEdge {
                        center: [cx, cy, 0.0],
                        tangent: [1.0, 0.0, 0.0],
                        length: dx,
                        div_rho: -1.0,
                    });
                }
            }
        }

        // 垂直边（切向 y）
        for j in 0..ny {
            for i in 0..=nx {
                let cx = i as f64 * dx;
                let cy = (j as f64 + 0.5) * dy;
                let is_boundary = i == 0 || i == nx;
                if !is_boundary {
                    edges.push(PlanarEdge {
                        center: [cx, cy, 0.0],
                        tangent: [0.0, 1.0, 0.0],
                        length: dy,
                        div_rho: 1.0,
                    });
                    edges.push(PlanarEdge {
                        center: [cx, cy, 0.0],
                        tangent: [0.0, 1.0, 0.0],
                        length: dy,
                        div_rho: -1.0,
                    });
                }
            }
        }

        Self { edges, nx, ny, dx, dy }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_1x1_grid() {
        // 1×1 矩形只有 1 条斜内边
        let bases = build_planar_rwg(1.0, 1.0, 1, 1);
        assert_eq!(bases.len(), 1);
        assert!((bases[0].edge_len - std::f64::consts::SQRT_2).abs() < 1e-10);
    }

    #[test]
    fn test_2x2_grid() {
        // 2×2：斜边 4，水平内边 2×1=2，垂直内边 1×2=2，共 8
        let bases = build_planar_rwg(1.0, 1.0, 2, 2);
        assert_eq!(bases.len(), 4 + 2 + 2);
    }
}
