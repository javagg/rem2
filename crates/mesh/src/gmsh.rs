/// Gmsh reader adapter built on rmsh-io.
///
/// Keeps rem's `RawMesh` interface stable while delegating MSH parsing to
/// `rmsh-io` (including v2.2 ASCII/binary and v4.1 ASCII support).
use rem_core::{RemError, RemResult};
use std::path::Path;

#[derive(Debug)]
pub struct RawElement {
    pub id: usize,
    pub elem_type: u32,
    pub phys_tag: u32,
    pub node_ids: Vec<usize>, // 1-based GMSH node indices
}

#[derive(Debug)]
pub struct RawMesh {
    /// (id, x, y, z) in GMSH units (before L0 scaling)
    pub nodes: Vec<(usize, f64, f64, f64)>,
    pub elements: Vec<RawElement>,
}

pub fn read_msh_file(path: &Path) -> RemResult<RawMesh> {
    let mesh = rmsh_io::load_msh_from_path(path)
        .map_err(|e| RemError::Mesh(format!("rmsh-io load_msh_from_path failed: {}", e)))?;
    Ok(to_raw_mesh(mesh))
}

pub fn read_msh_str(text: &str) -> RemResult<RawMesh> {
    read_msh_bytes(text.as_bytes())
}

pub fn read_msh_bytes(bytes: &[u8]) -> RemResult<RawMesh> {
    let mesh = rmsh_io::load_msh_from_bytes(bytes)
        .map_err(|e| RemError::Mesh(format!("rmsh-io load_msh_from_bytes failed: {}", e)))?;
    Ok(to_raw_mesh(mesh))
}

fn to_raw_mesh(mesh: rmsh_model::Mesh) -> RawMesh {
    let mut nodes: Vec<(usize, f64, f64, f64)> = mesh
        .nodes
        .values()
        .map(|n| (n.id as usize, n.position.x, n.position.y, n.position.z))
        .collect();
    // Preserve legacy expectation that node index roughly follows node tag ordering.
    nodes.sort_by_key(|(id, _, _, _)| *id);

    let mut elements: Vec<RawElement> = mesh
        .elements
        .into_iter()
        .map(|e| RawElement {
            id: e.id as usize,
            elem_type: gmsh_type_for_element(e.etype),
            phys_tag: e.physical_tag.unwrap_or(0).max(0) as u32,
            node_ids: e.node_ids.into_iter().map(|id| id as usize).collect(),
        })
        .collect();
    elements.sort_by_key(|e| e.id);

    RawMesh { nodes, elements }
}

fn gmsh_type_for_element(etype: rmsh_model::ElementType) -> u32 {
    match etype {
        rmsh_model::ElementType::Line2 => 1,
        rmsh_model::ElementType::Triangle3 => 2,
        rmsh_model::ElementType::Quad4 => 3,
        rmsh_model::ElementType::Tetrahedron4 => 4,
        rmsh_model::ElementType::Hexahedron8 => 5,
        rmsh_model::ElementType::Prism6 => 6,
        rmsh_model::ElementType::Pyramid5 => 7,
        rmsh_model::ElementType::Point1 => 15,
        rmsh_model::ElementType::Unknown(code) => code.max(0) as u32,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const SIMPLE_MSH: &str = r#"$MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 4 1 4
2 1 0 4
1
2
3
4
0.0 0.0 0.0
1.0 0.0 0.0
1.0 1.0 0.0
0.0 1.0 0.0
$EndNodes
$Elements
1 2 1 2
2 7 2 2
1 1 2 3
2 2 3 4
$EndElements
"#;

    #[test]
    fn parse_simple_mesh_via_rmsh() {
        let raw = read_msh_str(SIMPLE_MSH).expect("simple mesh should parse");
        assert_eq!(raw.nodes.len(), 4);
        assert_eq!(raw.elements.len(), 2);
        assert_eq!(raw.elements[0].elem_type, 2);
        assert_eq!(raw.elements[0].phys_tag, 7);
        assert_eq!(raw.elements[0].node_ids, vec![1, 2, 3]);
    }

    #[test]
    fn parse_v2_ascii_mesh_via_rmsh() {
        let v2 = r#"$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
3
1 0.0 0.0 0.0
2 1.0 0.0 0.0
3 1.0 1.0 0.0
$EndNodes
$Elements
1
1 2 2 10 1 1 2 3
$EndElements
"#;
        let raw = read_msh_str(v2).unwrap();
        assert_eq!(raw.nodes.len(), 3);
        assert_eq!(raw.elements.len(), 1);
        assert_eq!(raw.elements[0].phys_tag, 10);
        assert_eq!(raw.elements[0].node_ids, vec![1, 2, 3]);
    }
}
