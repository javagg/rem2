//! 子域划分 — 调用 rmetis 将体网格分成 n_parts 个子域。

use rem_core::{RemError, RemResult};
use rem_mesh::RemMesh;

/// 将网格划分为 n_parts 个子域。
/// 返回每个体单元的子域编号（0-based），长度 = mesh.volume_elements.len()。
pub fn partition_mesh(mesh: &RemMesh, n_parts: usize) -> RemResult<Vec<i32>> {
    let n_elems = mesh.volume_elements.len();
    if n_elems == 0 {
        return Err(RemError::Config("DDM: mesh has no volume elements".to_string()));
    }
    if n_parts <= 1 {
        return Ok(vec![0i32; n_elems]);
    }

    // 构建单元-节点邻接表（用于 METIS 双偶图划分）
    // 简化实现：直接按单元索引均匀分配
    // TODO: 调用 rmetis::partition_mesh_dual() 做真正的几何平衡划分
    let partition: Vec<i32> = (0..n_elems)
        .map(|i| (i * n_parts / n_elems) as i32)
        .collect();

    log::info!("DDM partition: {} elements → {} subdomains", n_elems, n_parts);
    Ok(partition)
}

/// 统计每个子域的单元数量。
pub fn partition_stats(partition: &[i32], n_parts: usize) -> Vec<usize> {
    let mut counts = vec![0usize; n_parts];
    for &p in partition {
        if (p as usize) < n_parts {
            counts[p as usize] += 1;
        }
    }
    counts
}
