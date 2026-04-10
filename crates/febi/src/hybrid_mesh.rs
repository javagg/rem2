//! 辐射边界提取 — 从 RemMesh 中识别 FE-BI 截断面 Γ
//!
//! 复用 rem-mom 的 SurfaceMesh，只需按辐射边界属性过滤。

use rem_core::RemResult;
use rem_mesh::RemMesh;
use rem_mom::surface_mesh::SurfaceMesh;

/// 辐射边界网格（BI 截断面 Γ）
/// 本质上是 rem-mom 的 SurfaceMesh，但按辐射边界属性过滤。
pub type RadiationMesh = SurfaceMesh;

/// 从体网格中提取辐射边界面（FE-BI 截断面 Γ）。
///
/// `rad_attrs`：在配置文件 `Solver.FEBI.RadiationBoundary.Attributes` 中指定的
/// 物理组属性 ID 列表，对应网格中的辐射截断外表面。
pub fn extract_radiation_boundary(
    mesh: &RemMesh,
    rad_attrs: &[u32],
) -> RemResult<RadiationMesh> {
    // 复用 MoM 的 SurfaceMesh::extract，按辐射边界属性提取表面
    SurfaceMesh::extract(mesh, rad_attrs)
}
