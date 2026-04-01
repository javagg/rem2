/// Minimal GMSH v4.1 ASCII reader.
///
/// Supports `$MeshFormat 4.1`, `$Nodes`, `$Elements`, `$PhysicalNames`.
/// Does NOT support binary format or MSH 2.x.
use rem_core::{RemError, RemResult};
use std::path::Path;

// ---------------------------------------------------------------------------
// Raw data structures (before unit scaling / tag resolution)
// ---------------------------------------------------------------------------

#[derive(Debug)]
pub struct RawElement {
    pub id:        usize,
    pub elem_type: u32,
    pub phys_tag:  u32,
    pub node_ids:  Vec<usize>, // 1-based GMSH node indices
}

#[derive(Debug)]
pub struct RawMesh {
    /// (id, x, y, z) in GMSH units (before L0 scaling)
    pub nodes:    Vec<(usize, f64, f64, f64)>,
    pub elements: Vec<RawElement>,
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

pub fn read_msh_file(path: &Path) -> RemResult<RawMesh> {
    let content = std::fs::read_to_string(path).map_err(RemError::Io)?;
    read_msh_str(&content)
}

pub fn read_msh_str(text: &str) -> RemResult<RawMesh> {
    let mut lines = text.lines().peekable();
    let mut nodes: Vec<(usize, f64, f64, f64)> = Vec::new();
    let mut elements: Vec<RawElement> = Vec::new();
    // physical group tag → entity tag → set of element tags (for tag resolution)
    // We store a simple tag→phys_tag map built from $Entities or $PhysicalNames.
    // In practice we read the entity tag embedded in $Elements block.

    while let Some(line) = lines.next() {
        let line = line.trim();
        match line {
            "$MeshFormat" => {
                let fmt_line = lines.next().ok_or_else(|| bad("missing $MeshFormat data"))?;
                let parts: Vec<&str> = fmt_line.split_whitespace().collect();
                if parts.len() < 3 {
                    return Err(bad("malformed $MeshFormat"));
                }
                let version = parts[0];
                if !version.starts_with("4") {
                    return Err(RemError::Mesh(format!(
                        "Only GMSH v4.x is supported, got version {}",
                        version
                    )));
                }
                // skip to $EndMeshFormat
                skip_to_end(&mut lines, "EndMeshFormat")?;
            }
            "$Nodes" => {
                nodes = parse_nodes_v4(&mut lines)?;
            }
            "$Elements" => {
                elements = parse_elements_v4(&mut lines)?;
            }
            "$PhysicalNames" | "$Entities" | "$Partitioned" | "$Periodic"
            | "$GhostElements" | "$Parametrizations" | "$NodeData"
            | "$ElementData" | "$ElementNodeData" | "$InterpolationScheme" => {
                // skip these sections
                let end_tag = format!("End{}", &line[1..]);
                skip_to_end(&mut lines, &end_tag)?;
            }
            _ => {
                // ignore unknown lines
            }
        }
    }

    if nodes.is_empty() {
        return Err(bad("no nodes found in mesh file"));
    }
    if elements.is_empty() {
        return Err(bad("no elements found in mesh file"));
    }

    Ok(RawMesh { nodes, elements })
}

// ---------------------------------------------------------------------------
// Section parsers
// ---------------------------------------------------------------------------

fn parse_nodes_v4(
    lines: &mut std::iter::Peekable<std::str::Lines>,
) -> RemResult<Vec<(usize, f64, f64, f64)>> {
    // Header: numEntityBlocks numNodes minNodeTag maxNodeTag
    let header = lines.next().ok_or_else(|| bad("missing $Nodes header"))?;
    let hparts: Vec<&str> = header.split_whitespace().collect();
    if hparts.len() < 2 {
        return Err(bad("malformed $Nodes header"));
    }
    let _n_blocks: usize = parse_usize(hparts[0])?;
    let n_nodes: usize   = parse_usize(hparts[1])?;

    let mut nodes = Vec::with_capacity(n_nodes);

    loop {
        // Entity block header: entityDim entityTag parametric numNodesInBlock
        let block_header = lines.next().ok_or_else(|| bad("unexpected end in $Nodes"))?;
        let block_header = block_header.trim();
        if block_header.starts_with("$EndNodes") {
            break;
        }
        let bp: Vec<&str> = block_header.split_whitespace().collect();
        if bp.len() < 4 {
            return Err(bad("malformed node block header"));
        }
        let _parametric: u8  = bp[2].parse().unwrap_or(0);
        let n_in_block: usize = parse_usize(bp[3])?;

        // Read node tags (one per line)
        let mut tags = Vec::with_capacity(n_in_block);
        for _ in 0..n_in_block {
            let tl = lines.next().ok_or_else(|| bad("unexpected end reading node tags"))?;
            tags.push(parse_usize(tl.trim())?);
        }
        // Read coordinates (one per line)
        for tag in tags {
            let cl = lines.next().ok_or_else(|| bad("unexpected end reading node coords"))?;
            let cp: Vec<&str> = cl.split_whitespace().collect();
            if cp.len() < 3 {
                return Err(bad("malformed node coordinate line"));
            }
            let x: f64 = parse_f64(cp[0])?;
            let y: f64 = parse_f64(cp[1])?;
            let z: f64 = parse_f64(cp[2])?;
            nodes.push((tag, x, y, z));
        }
    }

    Ok(nodes)
}

fn parse_elements_v4(
    lines: &mut std::iter::Peekable<std::str::Lines>,
) -> RemResult<Vec<RawElement>> {
    // Header: numEntityBlocks numElements minElementTag maxElementTag
    let header = lines.next().ok_or_else(|| bad("missing $Elements header"))?;
    let hparts: Vec<&str> = header.split_whitespace().collect();
    if hparts.len() < 2 {
        return Err(bad("malformed $Elements header"));
    }
    let _n_blocks: usize  = parse_usize(hparts[0])?;
    let n_elements: usize = parse_usize(hparts[1])?;

    let mut elements = Vec::with_capacity(n_elements);

    loop {
        let block_header = lines.next().ok_or_else(|| bad("unexpected end in $Elements"))?;
        let block_header = block_header.trim();
        if block_header.starts_with("$EndElements") {
            break;
        }
        let bp: Vec<&str> = block_header.split_whitespace().collect();
        if bp.len() < 4 {
            return Err(bad("malformed element block header"));
        }
        // entityDim  entityTag  elementType  numElementsInBlock
        let phys_tag: u32    = parse_u32(bp[1])?;  // use entity tag as phys tag
        let elem_type: u32   = parse_u32(bp[2])?;
        let n_in_block: usize = parse_usize(bp[3])?;

        let n_nodes_per = n_nodes_for_type(elem_type).unwrap_or(0);

        for _ in 0..n_in_block {
            let el = lines.next().ok_or_else(|| bad("unexpected end reading elements"))?;
            let ep: Vec<&str> = el.split_whitespace().collect();
            if ep.is_empty() {
                return Err(bad("empty element line"));
            }
            let id: usize = parse_usize(ep[0])?;

            let node_ids: Vec<usize> = if n_nodes_per > 0 {
                (1..=n_nodes_per)
                    .map(|i| {
                        ep.get(i)
                            .ok_or_else(|| bad("too few node IDs in element"))
                            .and_then(|s| parse_usize(s))
                    })
                    .collect::<RemResult<_>>()?
            } else {
                // unknown type — read remaining tokens
                ep[1..].iter().map(|s| parse_usize(s)).collect::<RemResult<_>>()?
            };

            elements.push(RawElement { id, elem_type, phys_tag, node_ids });
        }
    }

    Ok(elements)
}

// ---------------------------------------------------------------------------
// Utilities
// ---------------------------------------------------------------------------

fn skip_to_end(
    lines: &mut std::iter::Peekable<std::str::Lines>,
    end_tag: &str,
) -> RemResult<()> {
    for line in lines.by_ref() {
        if line.trim().starts_with('$') && line.trim()[1..] == *end_tag {
            return Ok(());
        }
    }
    Err(bad(&format!("missing ${}", end_tag)))
}

fn bad(msg: &str) -> RemError {
    RemError::Mesh(msg.to_string())
}

fn parse_usize(s: &str) -> RemResult<usize> {
    s.trim().parse().map_err(|_| bad(&format!("expected integer, got '{}'", s)))
}
fn parse_u32(s: &str) -> RemResult<u32> {
    s.trim().parse().map_err(|_| bad(&format!("expected u32, got '{}'", s)))
}
fn parse_f64(s: &str) -> RemResult<f64> {
    s.trim().parse().map_err(|_| bad(&format!("expected float, got '{}'", s)))
}

/// Returns the node count for known GMSH element types.
fn n_nodes_for_type(t: u32) -> Option<usize> {
    match t {
        1  => Some(2),   // Line 2
        2  => Some(3),   // Tri 3
        3  => Some(4),   // Quad 4
        4  => Some(4),   // Tet 4
        5  => Some(8),   // Hex 8
        6  => Some(6),   // Prism 6
        7  => Some(5),   // Pyramid 5
        8  => Some(3),   // Line 3
        9  => Some(6),   // Tri 6
        10 => Some(9),   // Quad 9
        11 => Some(10),  // Tet 10
        15 => Some(1),   // Point
        _  => None,
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Minimal GMSH v4.1 ASCII mesh: 4 nodes, 2 triangles
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
2 1 2 2
1 1 2 3
2 2 3 4
$EndElements
"#;

    #[test]
    fn parse_simple_mesh() {
        let raw = read_msh_str(SIMPLE_MSH).expect("simple mesh should parse");
        assert_eq!(raw.nodes.len(), 4);
        assert_eq!(raw.elements.len(), 2);
        assert_eq!(raw.elements[0].elem_type, 2); // Tri3
        assert_eq!(raw.elements[0].node_ids, vec![1, 2, 3]);
        assert_eq!(raw.elements[1].node_ids, vec![2, 3, 4]);
    }

    #[test]
    fn node_coords_correct() {
        let raw = read_msh_str(SIMPLE_MSH).unwrap();
        let (_, x, y, z) = raw.nodes[0];
        assert!((x - 0.0).abs() < 1e-15);
        assert!((y - 0.0).abs() < 1e-15);
        assert!((z - 0.0).abs() < 1e-15);
        let (_, x, _, _) = raw.nodes[1];
        assert!((x - 1.0).abs() < 1e-15);
    }

    #[test]
    fn reject_v2_mesh() {
        let v2 = "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n";
        assert!(read_msh_str(v2).is_err());
    }
}
