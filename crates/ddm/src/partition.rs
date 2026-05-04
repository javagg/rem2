//! 子域划分 — 调用 rmetis 将体网格分成 n_parts 个子域。

use rem_core::{RemError, RemResult};
use rem_mesh::RemMesh;

/// Partition mesh elements into `n_parts` subdomains using METIS dual-graph
/// partitioning. Falls back to strip partitioning if METIS is unavailable.
///
/// Returns a per-element subdomain index (0-based), length = n_volume_elements.
pub fn partition_mesh(mesh: &RemMesh, n_parts: usize) -> RemResult<Vec<i32>> {
    let n_elems = mesh.volume_elements.len();
    if n_elems == 0 {
        return Err(RemError::Config("DDM: mesh has no volume elements".to_string()));
    }
    if n_parts <= 1 {
        return Ok(vec![0i32; n_elems]);
    }

    // Try METIS dual-graph partitioning
    match partition_with_rmetis(mesh, n_elems, n_parts) {
        Ok(partition) => {
            log::info!(
                "DDM rmetis partition: {} elements → {} subdomains",
                n_elems, n_parts
            );
            let stats = partition_stats(&partition, n_parts);
            log::info!("  Subdomain sizes: {:?}", stats);
            return Ok(partition);
        }
        Err(e) => {
            log::warn!("rmetis partitioning failed ({}), using strip fallback", e);
        }
    }

    // Fallback: simple strip partition
    let partition: Vec<i32> = (0..n_elems)
        .map(|i| (i * n_parts / n_elems) as i32)
        .collect();
    log::info!(
        "DDM strip partition (fallback): {} elements → {} subdomains",
        n_elems, n_parts
    );
    Ok(partition)
}

/// Build a dual-graph from mesh elements and partition via rmetis.
fn partition_with_rmetis(
    mesh: &RemMesh,
    n_elems: usize,
    n_parts: usize,
) -> RemResult<Vec<i32>> {
    // Build element adjacency: elements sharing a face are neighbors
    let mut adj_map: Vec<Vec<usize>> = vec![Vec::new(); n_elems];

    // For tetrahedral meshes, elements sharing 3 nodes share a face
    use std::collections::HashMap;
    let mut face_to_elems: HashMap<[usize; 3], Vec<usize>> = HashMap::new();

    for (ei, elem) in mesh.volume_elements.iter().enumerate() {
        let nodes = &elem.node_ids;
        let n = nodes.len();
        if n < 3 {
            continue;
        }
        // Extract triangular faces
        let faces: Vec<[usize; 3]> = match n {
            4 => {
                // Tetrahedron faces: (0,1,2), (0,1,3), (0,2,3), (1,2,3)
                vec![
                    [nodes[0], nodes[1], nodes[2]],
                    [nodes[0], nodes[1], nodes[3]],
                    [nodes[0], nodes[2], nodes[3]],
                    [nodes[1], nodes[2], nodes[3]],
                ]
            }
            _ => continue,
        };
        for mut face in faces {
            face.sort_unstable();
            face_to_elems.entry(face).or_default().push(ei);
        }
    }

    for (_face, elems) in face_to_elems {
        if elems.len() == 2 {
            adj_map[elems[0]].push(elems[1]);
            adj_map[elems[1]].push(elems[0]);
        }
    }

    // Build METIS graph in CSR format
    let mut xadj = vec![0i32];
    let mut adjncy = Vec::new();
    for adj in &adj_map {
        let mut sorted: Vec<i32> = adj.iter().map(|&v| v as i32).collect();
        sorted.sort_unstable();
        sorted.dedup();
        adjncy.extend(&sorted);
        xadj.push(adjncy.len() as i32);
    }

    // Call rmetis
    let graph = rmetis::graph::Graph::new(n_elems as i32, &xadj, &adjncy)
        .map_err(|e| RemError::Other(format!("rmetis graph build: {e}")))?;

    let options = rmetis::types::Options::default();
    let result = rmetis::part_graph_kway(&graph, n_parts, None, None, &options)
        .map_err(|e| RemError::Other(format!("rmetis partition: {e}")))?;

    Ok(result.part)
}

/// Count elements in each subdomain.
pub fn partition_stats(partition: &[i32], n_parts: usize) -> Vec<usize> {
    let mut counts = vec![0usize; n_parts];
    for &p in partition {
        if (p as usize) < n_parts {
            counts[p as usize] += 1;
        }
    }
    counts
}

#[cfg(test)]
mod tests {
    use super::*;
    use rem_mesh::{Node, Element, ElementKind, RemMesh};
    use std::collections::HashMap;

    fn two_tet_mesh() -> RemMesh {
        RemMesh {
            nodes: vec![
                Node { id: 0, x: 0.0, y: 0.0, z: 0.0 },
                Node { id: 1, x: 1.0, y: 0.0, z: 0.0 },
                Node { id: 2, x: 0.0, y: 1.0, z: 0.0 },
                Node { id: 3, x: 0.0, y: 0.0, z: 1.0 },
                Node { id: 4, x: 1.0, y: 0.0, z: 1.0 },
            ],
            volume_elements: vec![
                Element { id: 1, kind: ElementKind::Tet4, tag: 1, node_ids: vec![0, 1, 2, 3], rank: 0 },
                Element { id: 2, kind: ElementKind::Tet4, tag: 1, node_ids: vec![1, 3, 2, 4], rank: 0 },
            ],
            boundary_elements: vec![],
            domain_tags: HashMap::new(),
            boundary_tags: HashMap::new(),
            dim: 3,
            rank: 0,
            size: 1,
        }
    }

    #[test]
    fn single_part_returns_all_zero() {
        let mesh = two_tet_mesh();
        let part = partition_mesh(&mesh, 1).unwrap();
        assert_eq!(part, vec![0, 0]);
    }

    #[test]
    fn two_parts_produces_valid_partition() {
        let mesh = two_tet_mesh();
        let part = partition_mesh(&mesh, 2).unwrap();
        assert_eq!(part.len(), 2);
        let stats = partition_stats(&part, 2);
        // Both subdomains should have elements
        assert!(stats.iter().all(|&c| c > 0));
    }
}
