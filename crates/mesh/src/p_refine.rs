//! P-refinement: promote a P1 mesh to P2 by adding edge midpoint nodes.
//!
//! Converts each element kind:
//!   - `Tri3`  → `Tri6`   (3 corners + 3 edge midpoints)
//!   - `Tet4`  → `Tet10`  (4 corners + 6 edge midpoints)
//!   - `Line2` → `Line3`  (2 endpoints + 1 midpoint — for boundary elements)
//!
//! Other element types (`Quad4`, `Hex8`, existing `Tri6`/`Tet10`/`Line3`) are
//! returned unchanged.
//!
//! # DOF compatibility
//!
//! Original node indices 0..n_old are **preserved exactly**.  New midpoint
//! nodes are appended starting at index `n_old`.  This means:
//!
//! - `collect_dirichlet_dofs` — iterates `boundary_elements[*].node_ids`; the
//!   Line3 midpoints are included automatically.
//! - `assemble_stiffness` / `assemble_mass` in `rem_electrostatic::assemble`
//!   dispatch on `ElementKind::Tri6` / `Tet10` and use all 6/10 node ids.
//!
//! After `p_refine_mesh`, calling `assemble_stiffness(p2_mesh, ...)` correctly
//! assembles the P2 stiffness with the proper quadratic shape functions.

use std::collections::HashMap;

use crate::mesh_data::{Element, ElementKind, Node, RemMesh};

/// Unique sorted edge key.
#[inline]
fn edge_key(a: usize, b: usize) -> (usize, usize) {
    if a < b { (a, b) } else { (b, a) }
}

/// Promote a P1 mesh to P2 by adding edge midpoint DOF nodes.
///
/// Elements that are already P2 (`Tri6`, `Tet10`, `Line3`) or unsupported
/// (`Quad4`, `Hex8`) are passed through without modification.
///
/// Returns the new P2 `RemMesh`.  The original mesh is unmodified.
pub fn p_refine_mesh(mesh: &RemMesh) -> RemMesh {
    // --- Phase 1: collect all edges that need midpoints ---
    let mut midpoint_map: HashMap<(usize, usize), usize> = HashMap::new();
    let mut new_nodes: Vec<Node> = mesh.nodes.clone();
    let mut next_id = new_nodes.len();

    // Helper: ensure a midpoint for edge (a, b) exists, return its index.
    let mut get_or_insert = |a: usize, b: usize, new_nodes: &mut Vec<Node>, next_id: &mut usize| -> usize {
        let key = edge_key(a, b);
        if let Some(&mid_idx) = midpoint_map.get(&key) {
            return mid_idx;
        }
        let na = &mesh.nodes[key.0];
        let nb = &mesh.nodes[key.1];
        let mid = Node {
            id: *next_id,
            x: 0.5 * (na.x + nb.x),
            y: 0.5 * (na.y + nb.y),
            z: 0.5 * (na.z + nb.z),
        };
        new_nodes.push(mid);
        midpoint_map.insert(key, *next_id);
        *next_id += 1;
        *next_id - 1
    };

    // Scan volume elements first to decide which edges need midpoints
    for elem in &mesh.volume_elements {
        match elem.kind {
            ElementKind::Tri3 => {
                let [n0, n1, n2] = [elem.node_ids[0], elem.node_ids[1], elem.node_ids[2]];
                get_or_insert(n0, n1, &mut new_nodes, &mut next_id);
                get_or_insert(n1, n2, &mut new_nodes, &mut next_id);
                get_or_insert(n0, n2, &mut new_nodes, &mut next_id);
            }
            ElementKind::Tet4 => {
                let [n0, n1, n2, n3] = [
                    elem.node_ids[0],
                    elem.node_ids[1],
                    elem.node_ids[2],
                    elem.node_ids[3],
                ];
                get_or_insert(n0, n1, &mut new_nodes, &mut next_id);
                get_or_insert(n0, n2, &mut new_nodes, &mut next_id);
                get_or_insert(n0, n3, &mut new_nodes, &mut next_id);
                get_or_insert(n1, n2, &mut new_nodes, &mut next_id);
                get_or_insert(n1, n3, &mut new_nodes, &mut next_id);
                get_or_insert(n2, n3, &mut new_nodes, &mut next_id);
            }
            // Already P2 or unsupported → no midpoints needed from this path
            _ => {}
        }
    }

    // --- Phase 2: rebuild volume elements with new kinds / node_ids ---
    let volume_elements: Vec<Element> = mesh
        .volume_elements
        .iter()
        .map(|elem| match elem.kind {
            ElementKind::Tri3 => {
                let [n0, n1, n2] = [elem.node_ids[0], elem.node_ids[1], elem.node_ids[2]];
                let m01 = midpoint_map[&edge_key(n0, n1)];
                let m12 = midpoint_map[&edge_key(n1, n2)];
                let m02 = midpoint_map[&edge_key(n0, n2)];
                // GMSH / FEM standard Tri6 ordering: [v0,v1,v2, m01,m12,m02]
                Element {
                    id: elem.id,
                    kind: ElementKind::Tri6,
                    tag: elem.tag,
                    node_ids: vec![n0, n1, n2, m01, m12, m02],
                    rank: elem.rank,
                }
            }
            ElementKind::Tet4 => {
                let [n0, n1, n2, n3] = [
                    elem.node_ids[0],
                    elem.node_ids[1],
                    elem.node_ids[2],
                    elem.node_ids[3],
                ];
                let m01 = midpoint_map[&edge_key(n0, n1)];
                let m02 = midpoint_map[&edge_key(n0, n2)];
                let m03 = midpoint_map[&edge_key(n0, n3)];
                let m12 = midpoint_map[&edge_key(n1, n2)];
                let m13 = midpoint_map[&edge_key(n1, n3)];
                let m23 = midpoint_map[&edge_key(n2, n3)];
                // GMSH Tet10 ordering: [v0,v1,v2,v3, m01,m12,m02,m03,m13,m23]
                Element {
                    id: elem.id,
                    kind: ElementKind::Tet10,
                    tag: elem.tag,
                    node_ids: vec![n0, n1, n2, n3, m01, m12, m02, m03, m13, m23],
                    rank: elem.rank,
                }
            }
            // All other kinds pass through unchanged
            _ => elem.clone(),
        })
        .collect();

    // --- Phase 3: rebuild boundary elements (Line2 → Line3) ---
    let boundary_elements: Vec<Element> = mesh
        .boundary_elements
        .iter()
        .map(|belem| match belem.kind {
            ElementKind::Line2 => {
                let [n0, n1] = [belem.node_ids[0], belem.node_ids[1]];
                let key = edge_key(n0, n1);
                if let Some(&mid) = midpoint_map.get(&key) {
                    // Promote to Line3: [n0, n1, midpoint]
                    Element {
                        id: belem.id,
                        kind: ElementKind::Line3,
                        tag: belem.tag,
                        node_ids: vec![n0, n1, mid],
                        rank: belem.rank,
                    }
                } else {
                    // Edge not found in volume (isolated boundary?); keep Line2
                    belem.clone()
                }
            }
            // Tri3 boundary faces (3-D surface BCs) → Tri6
            ElementKind::Tri3 => {
                let [n0, n1, n2] = [belem.node_ids[0], belem.node_ids[1], belem.node_ids[2]];
                let m01 = midpoint_map.get(&edge_key(n0, n1)).copied();
                let m12 = midpoint_map.get(&edge_key(n1, n2)).copied();
                let m02 = midpoint_map.get(&edge_key(n0, n2)).copied();
                if let (Some(m01), Some(m12), Some(m02)) = (m01, m12, m02) {
                    Element {
                        id: belem.id,
                        kind: ElementKind::Tri6,
                        tag: belem.tag,
                        node_ids: vec![n0, n1, n2, m01, m12, m02],
                        rank: belem.rank,
                    }
                } else {
                    belem.clone()
                }
            }
            _ => belem.clone(),
        })
        .collect();

    RemMesh {
        nodes: new_nodes,
        volume_elements,
        boundary_elements,
        domain_tags: mesh.domain_tags.clone(),
        boundary_tags: mesh.boundary_tags.clone(),
        dim: mesh.dim,
        rank: mesh.rank,
        size: mesh.size,
    }
}

// ---------------------------------------------------------------------------
// P3 refinement: Tri3 → Tri10 (cubic triangle)
// ---------------------------------------------------------------------------

/// Sorted edge key with direction: returns (min, max, flipped) where flipped=true when a>b.
#[allow(dead_code)]
#[inline]
fn edge_key_dir(a: usize, b: usize) -> (usize, usize, bool) {
    if a < b { (a, b, false) } else { (b, a, true) }
}

/// Promote a P1 mesh (Tri3) to P3 (Tri10) by adding cubic DOF nodes.
///
/// For each triangular element the following new nodes are created:
///   - 2 nodes per edge (at 1/3 and 2/3 of the edge, shared between adjacent elements)
///   - 1 interior bubble node per element (at the centroid)
///
/// Other element types are passed through unchanged.
///
/// GMSH/FEM Tri10 node ordering (barycentric coords):
///   0:(1,0,0)  1:(0,1,0)  2:(0,0,1)   ← corners
///   3:(2/3,1/3,0)  4:(1/3,2/3,0)      ← edge 0-1
///   5:(0,2/3,1/3)  6:(0,1/3,2/3)      ← edge 1-2
///   7:(1/3,0,2/3)  8:(2/3,0,1/3)      ← edge 0-2 (reversed: near v2, near v0)
///   9:(1/3,1/3,1/3)                    ← interior
pub fn p3_refine_mesh(mesh: &RemMesh) -> RemMesh {
    // Map from sorted edge (min,max) to [near_min_node, near_max_node]
    let mut edge_nodes: HashMap<(usize, usize), [usize; 2]> = HashMap::new();
    let mut new_nodes: Vec<Node> = mesh.nodes.clone();
    let mut next_id = new_nodes.len();

    // Helper: insert a point at (t · node_a + (1-t) · node_b), referencing original mesh.
    let insert_node = |a: usize, b: usize, t: f64, nodes: &mut Vec<Node>, nid: &mut usize| -> usize {
        let na = &mesh.nodes[a];
        let nb = &mesh.nodes[b];
        let nd = Node {
            id: *nid,
            x: t * na.x + (1.0 - t) * nb.x,
            y: t * na.y + (1.0 - t) * nb.y,
            z: t * na.z + (1.0 - t) * nb.z,
        };
        nodes.push(nd);
        *nid += 1;
        *nid - 1
    };

    // Pass 1: create 2 nodes per unique edge.
    for elem in &mesh.volume_elements {
        if elem.kind != ElementKind::Tri3 { continue; }
        let [n0, n1, n2] = [elem.node_ids[0], elem.node_ids[1], elem.node_ids[2]];
        for (a, b) in [(n0, n1), (n1, n2), (n0, n2)] {
            let key = edge_key(a, b);
            if edge_nodes.contains_key(&key) { continue; }
            // node near a (t=2/3 from a = 1/3 from b)  → at 2/3·a + 1/3·b
            let id_near_a = insert_node(a, b, 2.0/3.0, &mut new_nodes, &mut next_id);
            // node near b (t=1/3 from a = 2/3 from b) → at 1/3·a + 2/3·b
            let id_near_b = insert_node(a, b, 1.0/3.0, &mut new_nodes, &mut next_id);
            // store as [near_min, near_max]
            // edge_key() gives sorted order (min,max); track which endpoint is which
            let (near_min, near_max) = if a < b {
                (id_near_a, id_near_b)
            } else {
                (id_near_b, id_near_a)
            };
            edge_nodes.insert(key, [near_min, near_max]);
        }
    }

    // Pass 2: build volume elements.
    let volume_elements: Vec<Element> = mesh.volume_elements.iter().map(|elem| {
        if elem.kind != ElementKind::Tri3 {
            return elem.clone();
        }
        let [n0, n1, n2] = [elem.node_ids[0], elem.node_ids[1], elem.node_ids[2]];

        // Edge 0-1: [near n0, near n1]
        let e01 = edge_nodes[&edge_key(n0, n1)];
        let (e01_near0, e01_near1) = if n0 < n1 { (e01[0], e01[1]) } else { (e01[1], e01[0]) };

        // Edge 1-2: [near n1, near n2]
        let e12 = edge_nodes[&edge_key(n1, n2)];
        let (e12_near1, e12_near2) = if n1 < n2 { (e12[0], e12[1]) } else { (e12[1], e12[0]) };

        // Edge 0-2: [near n0, near n2]
        let e02 = edge_nodes[&edge_key(n0, n2)];
        let (e02_near0, e02_near2) = if n0 < n2 { (e02[0], e02[1]) } else { (e02[1], e02[0]) };

        // Interior bubble node at centroid (1/3, 1/3, 1/3)
        let nb = {
            let a = &mesh.nodes[n0];
            let b = &mesh.nodes[n1];
            let c = &mesh.nodes[n2];
            let nd = Node {
                id: next_id,
                x: (a.x + b.x + c.x) / 3.0,
                y: (a.y + b.y + c.y) / 3.0,
                z: (a.z + b.z + c.z) / 3.0,
            };
            new_nodes.push(nd);
            next_id += 1;
            next_id - 1
        };

        // Tri10 ordering (GMSH convention):
        //  [v0,v1,v2, e01_near0,e01_near1, e12_near1,e12_near2, e02_near2,e02_near0, interior]
        //
        // In GMSH barycentric notation (λ1=l1, λ2=l2, λ3=l3):
        //   node 3: (2/3,1/3,0) → near v0 on edge v0-v1
        //   node 4: (1/3,2/3,0) → near v1 on edge v0-v1
        //   node 5: (0,2/3,1/3) → near v1 on edge v1-v2
        //   node 6: (0,1/3,2/3) → near v2 on edge v1-v2
        //   node 7: (1/3,0,2/3) → near v2 on edge v0-v2 (note: reversed order for GMSH)
        //   node 8: (2/3,0,1/3) → near v0 on edge v0-v2
        //   node 9: (1/3,1/3,1/3) → interior
        Element {
            id: elem.id,
            kind: ElementKind::Tri10,
            tag: elem.tag,
            node_ids: vec![
                n0, n1, n2,
                e01_near0, e01_near1,
                e12_near1, e12_near2,
                e02_near2, e02_near0,
                nb,
            ],
            rank: elem.rank,
        }
    }).collect();

    // Pass 3: boundary elements — keep as-is (Tri10 FEM doesn't need special boundary promotion).
    let boundary_elements = mesh.boundary_elements.clone();

    RemMesh {
        nodes: new_nodes,
        volume_elements,
        boundary_elements,
        domain_tags: mesh.domain_tags.clone(),
        boundary_tags: mesh.boundary_tags.clone(),
        dim: mesh.dim,
        rank: mesh.rank,
        size: mesh.size,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesh_data::{BoundaryTag, ElementKind};
    use std::collections::HashMap;

    fn make_single_tri3_mesh() -> RemMesh {
        // One Tri3, one Line2 boundary edge
        RemMesh {
            nodes: vec![
                Node { id: 0, x: 0.0, y: 0.0, z: 0.0 },
                Node { id: 1, x: 1.0, y: 0.0, z: 0.0 },
                Node { id: 2, x: 0.0, y: 1.0, z: 0.0 },
            ],
            volume_elements: vec![Element {
                id: 0, kind: ElementKind::Tri3, tag: 1,
                node_ids: vec![0, 1, 2], rank: 0,
            }],
            boundary_elements: vec![Element {
                id: 1, kind: ElementKind::Line2, tag: 2,
                node_ids: vec![0, 1], rank: 0,
            }],
            domain_tags: HashMap::from([(1u32, 0usize)]),
            boundary_tags: HashMap::from([(2u32, BoundaryTag::Pec)]),
            dim: 2, rank: 0, size: 1,
        }
    }

    fn make_single_tet4_mesh() -> RemMesh {
        RemMesh {
            nodes: vec![
                Node { id: 0, x: 0.0, y: 0.0, z: 0.0 },
                Node { id: 1, x: 1.0, y: 0.0, z: 0.0 },
                Node { id: 2, x: 0.0, y: 1.0, z: 0.0 },
                Node { id: 3, x: 0.0, y: 0.0, z: 1.0 },
            ],
            volume_elements: vec![Element {
                id: 0, kind: ElementKind::Tet4, tag: 1,
                node_ids: vec![0, 1, 2, 3], rank: 0,
            }],
            boundary_elements: vec![Element {
                id: 1, kind: ElementKind::Tri3, tag: 2,
                node_ids: vec![0, 1, 2], rank: 0,
            }],
            domain_tags: HashMap::from([(1u32, 0usize)]),
            boundary_tags: HashMap::from([(2u32, BoundaryTag::Pec)]),
            dim: 3, rank: 0, size: 1,
        }
    }

    #[test]
    fn p_refine_tri3_to_tri6() {
        let p1 = make_single_tri3_mesh();
        let p2 = p_refine_mesh(&p1);

        // Tri3 has 3 edges → 3 midpoints added
        assert_eq!(p2.nodes.len(), 6, "P2 Tri3→Tri6 should have 6 nodes");
        assert_eq!(p2.volume_elements.len(), 1);
        assert_eq!(p2.volume_elements[0].kind, ElementKind::Tri6);
        assert_eq!(p2.volume_elements[0].node_ids.len(), 6);
        // Original corners preserved at same indices
        assert_eq!(p2.nodes[0].x, 0.0);
        assert_eq!(p2.nodes[1].x, 1.0);
        assert_eq!(p2.nodes[2].y, 1.0);
        // Boundary edge → Line3
        assert_eq!(p2.boundary_elements[0].kind, ElementKind::Line3);
        assert_eq!(p2.boundary_elements[0].node_ids.len(), 3);
    }

    #[test]
    fn p_refine_tet4_to_tet10() {
        let p1 = make_single_tet4_mesh();
        let p2 = p_refine_mesh(&p1);

        // Tet4 has 6 edges → 6 midpoints added
        assert_eq!(p2.nodes.len(), 10, "P2 Tet4→Tet10 should have 10 nodes");
        assert_eq!(p2.volume_elements[0].kind, ElementKind::Tet10);
        assert_eq!(p2.volume_elements[0].node_ids.len(), 10);
        // Boundary face → Tri6
        assert_eq!(p2.boundary_elements[0].kind, ElementKind::Tri6);
        assert_eq!(p2.boundary_elements[0].node_ids.len(), 6);
    }

    #[test]
    fn p_refine_midpoint_coordinates_correct() {
        let p1 = make_single_tri3_mesh();
        let p2 = p_refine_mesh(&p1);

        // Edge (0,1): midpoint should be at (0.5, 0.0, 0.0)
        // Edge (1,2): midpoint at (0.5, 0.5, 0.0)
        // Edge (0,2): midpoint at (0.0, 0.5, 0.0)
        let new_nodes = &p2.nodes[3..]; // midpoints start at index 3
        for n in new_nodes {
            // All midpoints should be at 0.5 coordinates (in {0, 0.5, 1})
            assert!(
                n.x == 0.0 || n.x == 0.5 || n.x == 1.0,
                "unexpected x = {}", n.x
            );
            assert!(
                n.y == 0.0 || n.y == 0.5 || n.y == 1.0,
                "unexpected y = {}", n.y
            );
        }
    }

    #[test]
    fn p_refine_idempotent_for_tri6() {
        // Applying p_refine on an already-P2 mesh should leave volume elements unchanged
        let p1 = make_single_tri3_mesh();
        let p2 = p_refine_mesh(&p1);
        let p2_again = p_refine_mesh(&p2);
        // Volume element should still be Tri6 (not converted again)
        assert_eq!(p2_again.volume_elements[0].kind, ElementKind::Tri6);
    }
}
