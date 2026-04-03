/// Minimal GMSH v2.2 and v4.1 ASCII/binary reader.
///
/// Supports `$MeshFormat 2.2` (ASCII and binary) and `4.1` (ASCII only).
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
    let bytes = std::fs::read(path).map_err(RemError::Io)?;
    read_msh_bytes(&bytes)
}

pub fn read_msh_str(text: &str) -> RemResult<RawMesh> {
    read_msh_bytes(text.as_bytes())
}

pub fn read_msh_bytes(bytes: &[u8]) -> RemResult<RawMesh> {
    // Detect binary v2: look for "2.2 1 " or "2.0 1 " anywhere in the first ~100 bytes
    // (the format line is always near the top, ASCII-readable)
    let header_region = &bytes[..bytes.len().min(128)];
    let is_binary_v2 = contains_bytes(header_region, b"$MeshFormat")
        && (contains_bytes(header_region, b"2.2 1 ")
            || contains_bytes(header_region, b"2.0 1 "));

    if is_binary_v2 {
        read_msh_v2_binary(bytes)
    } else {
        let text = std::str::from_utf8(bytes)
            .map_err(|_| RemError::Mesh("mesh file is not valid UTF-8 and not binary v2.2".into()))?;
        read_msh_ascii(text)
    }
}

fn contains_bytes(haystack: &[u8], needle: &[u8]) -> bool {
    haystack.windows(needle.len()).any(|w| w == needle)
}

/// Scan forward in `bytes` starting at `*pos` until the `marker` byte sequence is found.
/// Advances `*pos` past the marker and any trailing newline.
fn skip_to_marker(pos: &mut usize, bytes: &[u8], marker: &[u8]) {
    while *pos + marker.len() <= bytes.len() {
        if bytes[*pos..].starts_with(marker) {
            *pos += marker.len();
            if *pos < bytes.len() && bytes[*pos] == b'\n' { *pos += 1; }
            return;
        }
        *pos += 1;
    }
}

fn read_msh_ascii(text: &str) -> RemResult<RawMesh> {
    let mut lines = text.lines().peekable();
    let mut nodes: Vec<(usize, f64, f64, f64)> = Vec::new();
    let mut elements: Vec<RawElement> = Vec::new();
    let mut is_v2 = false;

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
                if version.starts_with('2') {
                    is_v2 = true;
                } else if !version.starts_with('4') {
                    return Err(RemError::Mesh(format!(
                        "Only GMSH v2.x or v4.x are supported, got version {}",
                        version
                    )));
                }
                // skip to $EndMeshFormat
                skip_to_end(&mut lines, "EndMeshFormat")?;
            }
            "$Nodes" => {
                nodes = if is_v2 {
                    parse_nodes_v2(&mut lines)?
                } else {
                    parse_nodes_v4(&mut lines)?
                };
            }
            "$Elements" => {
                elements = if is_v2 {
                    parse_elements_v2(&mut lines)?
                } else {
                    parse_elements_v4(&mut lines)?
                };
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

fn parse_nodes_v2(
    lines: &mut std::iter::Peekable<std::str::Lines>,
) -> RemResult<Vec<(usize, f64, f64, f64)>> {
    let header = lines.next().ok_or_else(|| bad("missing $Nodes header"))?;
    let n_nodes: usize = parse_usize(header.trim())?;
    let mut nodes = Vec::with_capacity(n_nodes);
    for _ in 0..n_nodes {
        let line = lines.next().ok_or_else(|| bad("unexpected end in $Nodes"))?;
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 4 {
            return Err(bad("malformed node line"));
        }
        let id = parse_usize(parts[0])?;
        let x  = parse_f64(parts[1])?;
        let y  = parse_f64(parts[2])?;
        let z  = parse_f64(parts[3])?;
        nodes.push((id, x, y, z));
    }
    skip_to_end(lines, "EndNodes")?;
    Ok(nodes)
}

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

fn parse_elements_v2(
    lines: &mut std::iter::Peekable<std::str::Lines>,
) -> RemResult<Vec<RawElement>> {
    let header = lines.next().ok_or_else(|| bad("missing $Elements header"))?;
    let n_elements: usize = parse_usize(header.trim())?;
    let mut elements = Vec::with_capacity(n_elements);
    for _ in 0..n_elements {
        let line = lines.next().ok_or_else(|| bad("unexpected end in $Elements"))?;
        let ep: Vec<&str> = line.split_whitespace().collect();
        if ep.len() < 5 {
            return Err(bad("malformed element line"));
        }
        let id: usize        = parse_usize(ep[0])?;
        let elem_type: u32   = parse_u32(ep[1])?;
        let n_tags: usize    = parse_usize(ep[2])?;
        if ep.len() < 3 + n_tags {
            return Err(bad("element line too short for tags"));
        }
        // GMSH v2 tags: first tag is usually physical group, second is elementary entity
        let phys_tag = if n_tags > 0 { parse_u32(ep[3])? } else { 0 };
        let n_nodes_per = n_nodes_for_type(elem_type).unwrap_or(0);
        if ep.len() < 3 + n_tags + n_nodes_per {
            return Err(bad("element line too short for nodes"));
        }
        let mut node_ids = Vec::with_capacity(n_nodes_per);
        for i in 0..n_nodes_per {
            node_ids.push(parse_usize(ep[3 + n_tags + i])?);
        }
        elements.push(RawElement { id, elem_type, phys_tag, node_ids });
    }
    skip_to_end(lines, "EndElements")?;
    Ok(elements)
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
        12 => Some(20),  // Hex 20
        13 => Some(15),  // Prism 15
        14 => Some(13),  // Pyramid 13
        15 => Some(1),   // Point
        16 => Some(8),   // Quad 8
        17 => Some(20),  // Hex 20 (serendipity)
        18 => Some(15),  // Prism 15
        19 => Some(13),  // Pyramid 13
        20 => Some(9),   // Tri 9
        21 => Some(10),  // Tri 10
        22 => Some(12),  // Tri 12 (incomplete)
        23 => Some(15),  // Tri 15 (serendipity)
        26 => Some(4),   // Line 4
        29 => Some(20),  // Tet 20
        36 => Some(16),  // Quad 16
        37 => Some(25),  // Hex 25
        92 => Some(27),  // Hex 27
        93 => Some(18),  // Prism 18
        _  => None,
    }
}

// ---------------------------------------------------------------------------
// GMSH v2.2 binary reader
// ---------------------------------------------------------------------------
//
// Format: text header sections with "$SectionName" / "$EndSectionName" markers,
// but $Nodes and $Elements bodies are encoded as packed binary (little-endian
// 32-bit ints + 64-bit doubles).
//
// $Nodes body layout:
//   int32 n_nodes
//   for each node:
//     int32 id,  double x, double y, double z   (total 28 bytes)
//
// $Elements body layout:
//   int32 n_element_types
//   for each element-type block:
//     int32 elem_type, int32 n_elems, int32 n_tags
//     for each element in block:
//       int32 id, int32[n_tags] tags, int32[n_nodes_per] node_ids

fn read_msh_v2_binary(bytes: &[u8]) -> RemResult<RawMesh> {
    let mut pos = 0usize;

    // Helper closures to read from byte slice
    let read_line = |pos: &mut usize| -> Option<&str> {
        let start = *pos;
        let slice = &bytes[start..];
        if let Some(nl) = slice.iter().position(|&b| b == b'\n') {
            let line = std::str::from_utf8(&slice[..nl]).ok()?.trim();
            *pos = start + nl + 1;
            Some(line)
        } else if !slice.is_empty() {
            let line = std::str::from_utf8(slice).ok()?.trim();
            *pos = bytes.len();
            Some(line)
        } else {
            None
        }
    };

    let read_i32_le = |pos: &mut usize, bytes: &[u8]| -> RemResult<i32> {
        if *pos + 4 > bytes.len() {
            return Err(bad("binary: unexpected end reading int32"));
        }
        let v = i32::from_le_bytes([bytes[*pos], bytes[*pos+1], bytes[*pos+2], bytes[*pos+3]]);
        *pos += 4;
        Ok(v)
    };

    let read_f64_le = |pos: &mut usize, bytes: &[u8]| -> RemResult<f64> {
        if *pos + 8 > bytes.len() {
            return Err(bad("binary: unexpected end reading float64"));
        }
        let v = f64::from_le_bytes([
            bytes[*pos],   bytes[*pos+1], bytes[*pos+2], bytes[*pos+3],
            bytes[*pos+4], bytes[*pos+5], bytes[*pos+6], bytes[*pos+7],
        ]);
        *pos += 8;
        Ok(v)
    };

    let mut nodes: Vec<(usize, f64, f64, f64)> = Vec::new();
    let mut elements: Vec<RawElement> = Vec::new();
    let mut phys_names: std::collections::HashMap<String, u32> = std::collections::HashMap::new();

    loop {
        // Skip blank lines / find section tag
        let Some(line) = read_line(&mut pos) else { break };
        let line = line.trim();
        if line.is_empty() { continue; }

        match line {
            "$MeshFormat" => {
                // Read the format line (ASCII)
                let Some(_fmt_line) = read_line(&mut pos) else { break };
                // After the ascii header line there is a binary int32 (=1) then newline
                // Advance past it using byte-level marker search
                skip_to_marker(&mut pos, bytes, b"$EndMeshFormat");
            }
            "$PhysicalNames" => {
                let Some(count_line) = read_line(&mut pos) else { break };
                let n_names: usize = count_line.trim().parse().unwrap_or(0);
                for _ in 0..n_names {
                    let Some(_name_line) = read_line(&mut pos) else { break };
                    // We don't need to map physical names in binary mode
                    // (element tags are numeric from the data)
                }
                skip_to_marker(&mut pos, bytes, b"$EndPhysicalNames");
            }
            "$Nodes" => {
                // ASCII count line
                let Some(count_line) = read_line(&mut pos) else { break };
                let n_nodes: usize = count_line.trim().parse()
                    .map_err(|_| bad("binary: bad node count"))?;
                nodes.reserve(n_nodes);
                for _ in 0..n_nodes {
                    let id = read_i32_le(&mut pos, bytes)? as usize;
                    let x  = read_f64_le(&mut pos, bytes)?;
                    let y  = read_f64_le(&mut pos, bytes)?;
                    let z  = read_f64_le(&mut pos, bytes)?;
                    nodes.push((id, x, y, z));
                }
                // Skip optional newline + mandatory $EndNodes\n (do NOT use marker search
                // because the literal bytes could appear inside binary data)
                if pos < bytes.len() && bytes[pos] == b'\n' { pos += 1; }
                if pos + 9 <= bytes.len() && &bytes[pos..pos+9] == b"$EndNodes" {
                    pos += 9;
                    if pos < bytes.len() && bytes[pos] == b'\n' { pos += 1; }
                }
            }
            "$Elements" => {
                // ASCII count line: total number of elements
                let Some(count_line) = read_line(&mut pos) else { break };
                let _n_total: usize = count_line.trim().parse().unwrap_or(0);

                // Read binary element-type blocks until $EndElements
                // Strategy: read binary blocks until we hit the $EndElements marker.
                // We use byte-level detection of the marker.
                loop {
                    if pos >= bytes.len() { break; }
                    // Check if we're at the $EndElements marker
                    if bytes[pos..].starts_with(b"$EndElements") {
                        pos += b"$EndElements".len();
                        // skip trailing newline
                        if pos < bytes.len() && bytes[pos] == b'\n' { pos += 1; }
                        break;
                    }
                    // Skip any leading newlines
                    if bytes[pos] == b'\n' { pos += 1; continue; }

                    // block header: int32 elem_type, int32 n_elems, int32 n_tags
                    let elem_type = read_i32_le(&mut pos, bytes)? as u32;
                    let n_elems   = read_i32_le(&mut pos, bytes)? as usize;
                    let n_tags    = read_i32_le(&mut pos, bytes)? as usize;
                    let n_nodes_per_opt = n_nodes_for_type(elem_type);

                    if n_elems == 0 { continue; }

                    match n_nodes_per_opt {
                        None => {
                            // Unknown element type: can't determine record size.
                            // Scan forward to $EndElements to avoid corrupting the stream.
                            log::debug!(
                                "binary v2.2: unknown element type {}, scanning to $EndElements",
                                elem_type
                            );
                            skip_to_marker(&mut pos, bytes, b"$EndElements");
                            break;
                        }
                        Some(n_nodes_per) => {
                            for _ in 0..n_elems {
                                let id = read_i32_le(&mut pos, bytes)? as usize;
                                let mut phys_tag: u32 = 0;
                                for ti in 0..n_tags {
                                    let tag = read_i32_le(&mut pos, bytes)? as u32;
                                    if ti == 0 { phys_tag = tag; }
                                }
                                let mut node_ids = Vec::with_capacity(n_nodes_per);
                                for _ in 0..n_nodes_per {
                                    node_ids.push(read_i32_le(&mut pos, bytes)? as usize);
                                }
                                if n_nodes_per > 0 {
                                    elements.push(RawElement { id, elem_type, phys_tag, node_ids });
                                }
                            }
                        }
                    }
                }
            }
            _ if line.starts_with('$') => {
                // Skip unknown section using byte-level marker search
                let end_marker = format!("$End{}", &line[1..]);
                skip_to_marker(&mut pos, bytes, end_marker.as_bytes());
            }
            _ => { /* ignore */ }
        }
    }

    if nodes.is_empty() {
        return Err(bad("binary v2.2: no nodes found"));
    }
    if elements.is_empty() {
        return Err(bad("binary v2.2: no elements found"));
    }
    Ok(RawMesh { nodes, elements })
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
    fn parse_v2_mesh() {
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
