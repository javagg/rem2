use crate::{Element, ElementKind, Node, RemMesh};
use std::collections::HashMap;

const MIDPOINT_TOL: f64 = 1.0e-12;

/// 2-D submesh extracted through fem-rs `SubMesh`, converted back into a `RemMesh`.
#[derive(Clone)]
pub struct FemSubMesh2d {
    pub mesh: RemMesh,
    pub parent_elem_ids: Vec<usize>,
    pub parent_node_ids: Vec<usize>,
}

impl FemSubMesh2d {
    /// Pull nodal values from the parent mesh onto the extracted submesh.
    pub fn transfer_from_parent(&self, parent_values: &[f64]) -> Vec<f64> {
        self.parent_node_ids
            .iter()
            .map(|&parent_node| parent_values[parent_node])
            .collect()
    }

    /// Scatter nodal values from the submesh back to the parent mesh by averaging.
    pub fn transfer_to_parent(&self, sub_values: &[f64], parent_n_nodes: usize) -> Vec<f64> {
        assert_eq!(sub_values.len(), self.parent_node_ids.len());
        let mut out = vec![0.0; parent_n_nodes];
        let mut counts = vec![0usize; parent_n_nodes];

        for (sub_node, &parent_node) in self.parent_node_ids.iter().enumerate() {
            out[parent_node] += sub_values[sub_node];
            counts[parent_node] += 1;
        }

        for parent_node in 0..parent_n_nodes {
            if counts[parent_node] > 0 {
                out[parent_node] /= counts[parent_node] as f64;
            }
        }

        out
    }
}

/// Refine a 2-D `Tri3` `RemMesh` through fem-rs' conforming refinement path.
///
/// Returns the refined mesh and a midpoint map compatible with rem2's
/// `prolongate_p1` helper.
pub fn refine_marked_tri3(
    mesh: &RemMesh,
    marked: &[usize],
) -> Result<(RemMesh, HashMap<(usize, usize), usize>), String> {
    ensure_tri3_mesh(mesh)?;

    let simplex = mesh.to_simplex_mesh_2d();
    let marked_fem: Vec<u32> = marked.iter().map(|&element_id| element_id as u32).collect();
    let refined = fem_mesh::refine_marked(&simplex, &marked_fem);
    let midpoint_map = derive_midpoint_map(mesh, &refined)?;
    let rem_mesh = remesh_from_simplex_2d(&refined, mesh)?;

    Ok((rem_mesh, midpoint_map))
}

/// Extract a 2-D `Tri3` submesh by volume tag through fem-rs `extract_submesh`.
pub fn extract_submesh_tri3(mesh: &RemMesh, element_tags: &[u32]) -> Result<FemSubMesh2d, String> {
    ensure_tri3_mesh(mesh)?;

    let simplex = mesh.to_simplex_mesh_2d();
    let fem_tags = element_tags
        .iter()
        .map(|&tag| i32::try_from(tag).map_err(|_| format!("element tag {tag} does not fit into i32")))
        .collect::<Result<Vec<_>, _>>()?;
    let sub = fem_mesh::extract_submesh(&simplex, &fem_tags);
    let rem_mesh = remesh_from_simplex_2d(&sub.mesh, mesh)?;

    Ok(FemSubMesh2d {
        mesh: rem_mesh,
        parent_elem_ids: sub.parent_elem_ids.iter().map(|&id| id as usize).collect(),
        parent_node_ids: sub.parent_node_of_sub.iter().map(|&id| id as usize).collect(),
    })
}

/// Extract a 2-D `Tri3` submesh by explicit parent volume element indices.
pub fn extract_submesh_by_element_ids_tri3(
    mesh: &RemMesh,
    element_ids: &[usize],
) -> Result<FemSubMesh2d, String> {
    ensure_tri3_mesh(mesh)?;

    if element_ids.is_empty() {
        return Err("at least one element id is required for submesh extraction".to_string());
    }

    let simplex = mesh.to_simplex_mesh_2d();
    let mut tagged = simplex.clone();
    tagged.elem_tags.fill(0);

    for &element_id in element_ids {
        let tag = tagged
            .elem_tags
            .get_mut(element_id)
            .ok_or_else(|| format!("element id {} is out of range", element_id))?;
        *tag = 1;
    }

    let sub = fem_mesh::extract_submesh(&tagged, &[1]);
    let rem_mesh = remesh_from_simplex_2d(&sub.mesh, mesh)?;

    Ok(FemSubMesh2d {
        mesh: rem_mesh,
        parent_elem_ids: sub.parent_elem_ids.iter().map(|&id| id as usize).collect(),
        parent_node_ids: sub.parent_node_of_sub.iter().map(|&id| id as usize).collect(),
    })
}

fn ensure_tri3_mesh(mesh: &RemMesh) -> Result<(), String> {
    if mesh.dim != 2 {
        return Err("fem_bridge currently supports only 2-D meshes".to_string());
    }
    if mesh.volume_elements.iter().any(|element| element.kind != ElementKind::Tri3) {
        return Err("fem_bridge currently supports only Tri3 volume meshes".to_string());
    }
    if mesh.boundary_elements.iter().any(|element| element.kind != ElementKind::Line2) {
        return Err("fem_bridge currently supports only Line2 boundary meshes".to_string());
    }
    Ok(())
}

fn remesh_from_simplex_2d(simplex: &fem_mesh::SimplexMesh<2>, template: &RemMesh) -> Result<RemMesh, String> {
    if simplex.elem_type != fem_mesh::ElementType::Tri3 {
        return Err(format!(
            "expected Tri3 simplex mesh, got {:?}",
            simplex.elem_type
        ));
    }
    if simplex.face_type != fem_mesh::ElementType::Line2 {
        return Err(format!(
            "expected Line2 simplex boundary mesh, got {:?}",
            simplex.face_type
        ));
    }

    let nodes = simplex
        .coords
        .chunks_exact(2)
        .enumerate()
        .map(|(idx, xy)| Node {
            id: idx,
            x: xy[0],
            y: xy[1],
            z: 0.0,
        })
        .collect::<Vec<_>>();

    let mut volume_elements = Vec::with_capacity(simplex.n_elems());
    for (idx, conn) in simplex.conn.chunks_exact(3).enumerate() {
        let tag = u32::try_from(simplex.elem_tags[idx])
            .map_err(|_| format!("negative element tag {} is unsupported in RemMesh", simplex.elem_tags[idx]))?;
        volume_elements.push(Element {
            id: idx + 1,
            kind: ElementKind::Tri3,
            tag,
            node_ids: conn.iter().map(|&node| node as usize).collect(),
            rank: template.rank,
        });
    }

    let mut boundary_elements = Vec::with_capacity(simplex.n_faces());
    for (idx, conn) in simplex.face_conn.chunks_exact(2).enumerate() {
        let tag = u32::try_from(simplex.face_tags[idx])
            .map_err(|_| format!("negative boundary tag {} is unsupported in RemMesh", simplex.face_tags[idx]))?;
        boundary_elements.push(Element {
            id: idx + 1,
            kind: ElementKind::Line2,
            tag,
            node_ids: conn.iter().map(|&node| node as usize).collect(),
            rank: template.rank,
        });
    }

    Ok(RemMesh {
        nodes,
        volume_elements,
        boundary_elements,
        domain_tags: template.domain_tags.clone(),
        boundary_tags: template.boundary_tags.clone(),
        dim: 2,
        rank: template.rank,
        size: template.size,
    })
}

fn derive_midpoint_map(
    coarse: &RemMesh,
    refined: &fem_mesh::SimplexMesh<2>,
) -> Result<HashMap<(usize, usize), usize>, String> {
    let coarse_n_nodes = coarse.n_nodes();
    if refined.n_nodes() < coarse_n_nodes {
        return Err("refined mesh has fewer nodes than the coarse mesh".to_string());
    }

    let mut midpoint_lookup: HashMap<(i64, i64), (usize, usize)> = HashMap::new();
    for element in &coarse.volume_elements {
        let ids = &element.node_ids;
        for &(a_local, b_local) in &[(0usize, 1usize), (1, 2), (0, 2)] {
            let a = ids[a_local];
            let b = ids[b_local];
            let key = midpoint_coord_key(&coarse.nodes[a], &coarse.nodes[b]);
            midpoint_lookup.entry(key).or_insert(edge_key(a, b));
        }
    }

    let mut midpoint_map = HashMap::new();
    for fine_node in coarse_n_nodes..refined.n_nodes() {
        let xy = refined.coords_of(fine_node as u32);
        let key = coord_key(xy[0], xy[1]);
        let edge = midpoint_lookup.get(&key).copied().ok_or_else(|| {
            format!(
                "could not recover coarse edge for refined midpoint node {} at ({:.6e}, {:.6e})",
                fine_node, xy[0], xy[1]
            )
        })?;
        midpoint_map.insert(edge, fine_node);
    }

    Ok(midpoint_map)
}

fn midpoint_coord_key(a: &Node, b: &Node) -> (i64, i64) {
    coord_key(0.5 * (a.x + b.x), 0.5 * (a.y + b.y))
}

fn coord_key(x: f64, y: f64) -> (i64, i64) {
    (
        (x / MIDPOINT_TOL).round() as i64,
        (y / MIDPOINT_TOL).round() as i64,
    )
}

fn edge_key(a: usize, b: usize) -> (usize, usize) {
    if a <= b { (a, b) } else { (b, a) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::BoundaryTag;
    use std::collections::HashMap;

    fn unit_square_mesh() -> RemMesh {
        let nodes = vec![
            Node { id: 0, x: 0.0, y: 0.0, z: 0.0 },
            Node { id: 1, x: 1.0, y: 0.0, z: 0.0 },
            Node { id: 2, x: 1.0, y: 1.0, z: 0.0 },
            Node { id: 3, x: 0.0, y: 1.0, z: 0.0 },
        ];
        let volume_elements = vec![
            Element { id: 1, kind: ElementKind::Tri3, tag: 1, node_ids: vec![0, 1, 2], rank: 0 },
            Element { id: 2, kind: ElementKind::Tri3, tag: 2, node_ids: vec![0, 2, 3], rank: 0 },
        ];
        let boundary_elements = vec![
            Element { id: 3, kind: ElementKind::Line2, tag: 10, node_ids: vec![0, 1], rank: 0 },
            Element { id: 4, kind: ElementKind::Line2, tag: 11, node_ids: vec![1, 2], rank: 0 },
            Element { id: 5, kind: ElementKind::Line2, tag: 12, node_ids: vec![2, 3], rank: 0 },
            Element { id: 6, kind: ElementKind::Line2, tag: 13, node_ids: vec![3, 0], rank: 0 },
        ];
        let mut boundary_tags = HashMap::new();
        boundary_tags.insert(10, BoundaryTag::Ground);
        boundary_tags.insert(11, BoundaryTag::Ground);
        boundary_tags.insert(12, BoundaryTag::Terminal { index: 1 });
        boundary_tags.insert(13, BoundaryTag::Terminal { index: 2 });

        RemMesh {
            nodes,
            volume_elements,
            boundary_elements,
            domain_tags: HashMap::from([(1, 0usize), (2, 1usize)]),
            boundary_tags,
            dim: 2,
            rank: 0,
            size: 1,
        }
    }

    #[test]
    fn refine_tri3_returns_midpoints_and_more_elements() {
        let mesh = unit_square_mesh();
        let (fine, midpoint_map) = refine_marked_tri3(&mesh, &[0, 1]).expect("fem bridge refinement should succeed");
        assert!(fine.n_nodes() > mesh.n_nodes());
        assert!(fine.n_volume_elements() > mesh.n_volume_elements());
        assert!(!midpoint_map.is_empty());
    }

    #[test]
    fn extract_submesh_preserves_parent_mapping() {
        let mesh = unit_square_mesh();
        let sub = extract_submesh_tri3(&mesh, &[2]).expect("fem bridge submesh extraction should succeed");
        assert_eq!(sub.parent_elem_ids, vec![1]);
        assert!(sub.mesh.volume_elements.iter().all(|element| element.tag == 2));
        let parent_values: Vec<f64> = (0..mesh.n_nodes()).map(|node| node as f64).collect();
        let sub_values = sub.transfer_from_parent(&parent_values);
        let scattered = sub.transfer_to_parent(&sub_values, mesh.n_nodes());
        for &parent_node in &sub.parent_node_ids {
            assert!((scattered[parent_node] - parent_values[parent_node]).abs() < 1.0e-12);
        }
    }

    #[test]
    fn extract_submesh_by_element_ids_preserves_local_selection() {
        let mesh = unit_square_mesh();
        let sub = extract_submesh_by_element_ids_tri3(&mesh, &[0])
            .expect("element-id submesh extraction should succeed");

        assert_eq!(sub.parent_elem_ids, vec![0]);
        assert_eq!(sub.mesh.n_volume_elements(), 1);
        assert_eq!(sub.mesh.volume_elements[0].tag, 1);
    }
}