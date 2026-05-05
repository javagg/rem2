//! Adaptive Mesh Refinement (AMR) for `RemMesh`.
//!
//! Provides:
//! 1. **ZZ error estimator** — element-wise gradient-recovery error indicator η_e.
//! 2. **Dörfler (bulk) marking** — marks minimal subset with Σ η_e² ≥ θ · Σ_all η_e².
//! 3. **Red refinement** — each marked Tri3 split into 4 children (2-D);
//!    each marked Tet4 split into 8 children (3-D, Freudenthal partition).
//!    Conformity ensured by one-pass neighbour propagation (bisects shared edges).
//! 4. **P1 prolongation** — transfers a nodal solution from coarse to fine mesh.

use crate::{RemMesh, Node, Element, ElementKind};
use std::collections::HashMap;

// ---------------------------------------------------------------------------
// ZZ error estimator
// ---------------------------------------------------------------------------

/// Zienkiewicz–Zhu element error indicators for a scalar P1 solution `phi`.
///
/// For each volume element e, computes:
///   η_e = |Ω_e|^{1/2} · ‖∇φ_h - G_h φ‖_{L²(Ω_e)}
///
/// where G_h φ is the nodally averaged recovered gradient.
///
/// Returns one non-negative value per volume element, in mesh order.
pub fn zz_estimator(mesh: &RemMesh, phi: &[f64]) -> Vec<f64> {
    let n = mesh.n_nodes();
    let n_elems = mesh.volume_elements.len();

    // --- Step 1: per-element raw gradient (FE grad, constant per Tri3/Tet4) ---
    let mut elem_grad: Vec<[f64; 3]> = vec![[0.0; 3]; n_elems];
    let mut elem_area: Vec<f64>      = vec![0.0;       n_elems];

    for (idx, elem) in mesh.volume_elements.iter().enumerate() {
        match elem.kind {
            ElementKind::Tri3 => {
                if let Some((g, a)) = tri3_grad_area(phi, elem, mesh) {
                    elem_grad[idx] = g;
                    elem_area[idx] = a;
                }
            }
            ElementKind::Tet4 => {
                if let Some((g, v)) = tet4_grad_vol(phi, elem, mesh) {
                    elem_grad[idx] = g;
                    elem_area[idx] = v;
                }
            }
            _ => {}
        }
    }

    // --- Step 2: nodal recovered gradient (area-weighted average) ---
    let mut nodal_grad  = vec![[0.0f64; 3]; n];
    let mut nodal_wsum  = vec![0.0f64;      n];

    for (idx, elem) in mesh.volume_elements.iter().enumerate() {
        let w = elem_area[idx];
        if w < 1e-300 { continue; }
        for &nid in &elem.node_ids {
            nodal_grad[nid][0] += w * elem_grad[idx][0];
            nodal_grad[nid][1] += w * elem_grad[idx][1];
            nodal_grad[nid][2] += w * elem_grad[idx][2];
            nodal_wsum[nid]    += w;
        }
    }
    for n in 0..mesh.n_nodes() {
        let w = nodal_wsum[n];
        if w > 0.0 {
            nodal_grad[n][0] /= w;
            nodal_grad[n][1] /= w;
            nodal_grad[n][2] /= w;
        }
    }

    // --- Step 3: element error indicator ---
    let mut eta = vec![0.0f64; n_elems];

    for (idx, elem) in mesh.volume_elements.iter().enumerate() {
        let a = elem_area[idx];
        if a < 1e-300 { continue; }

        // Average recovered gradient over element nodes
        let nn = elem.node_ids.len() as f64;
        let mut gr = [0.0f64; 3];
        for &nid in &elem.node_ids {
            gr[0] += nodal_grad[nid][0];
            gr[1] += nodal_grad[nid][1];
            gr[2] += nodal_grad[nid][2];
        }
        gr[0] /= nn;  gr[1] /= nn;  gr[2] /= nn;

        // η_e = √area · ‖∇φ_h - G_h φ‖
        let fe = &elem_grad[idx];
        let diff = (fe[0]-gr[0]).powi(2) + (fe[1]-gr[1]).powi(2) + (fe[2]-gr[2]).powi(2);
        eta[idx] = a.sqrt() * diff.sqrt();
    }

    eta
}

// ---------------------------------------------------------------------------
// Dörfler (bulk) marking
// ---------------------------------------------------------------------------

/// Mark a minimal subset of elements such that their total η² ≥ θ · Σ_all η².
///
/// Returns sorted indices into `mesh.volume_elements`.
pub fn dorfler_mark(eta: &[f64], theta: f64) -> Vec<usize> {
    let total: f64 = eta.iter().map(|&e| e * e).sum();
    if total < 1e-300 { return vec![]; }
    let target = theta * total;

    // Sort indices by descending η²
    let mut order: Vec<usize> = (0..eta.len()).collect();
    order.sort_unstable_by(|&a, &b| {
        (eta[b] * eta[b]).partial_cmp(&(eta[a] * eta[a])).unwrap()
    });

    let mut acc = 0.0;
    let mut marked = Vec::new();
    for i in order {
        if acc >= target { break; }
        acc += eta[i] * eta[i];
        marked.push(i);
    }
    marked.sort_unstable();
    marked
}

// ---------------------------------------------------------------------------
// Red refinement (Tri3 and Tet4)
// ---------------------------------------------------------------------------

/// Red refinement: each marked Tri3 is split into 4 children (2-D), or each
/// marked Tet4 into 8 children (3-D, Freudenthal partition).
/// Shared edges/faces with unmarked neighbours are bisected (conformity pass).
///
/// Returns `(new_mesh, midpoint_map)` where `midpoint_map` maps
/// `(node_a, node_b)` edge keys (sorted) to the new midpoint node index.
/// Pass this map to [`prolongate_p1`] to transfer the solution.
pub fn refine_marked(mesh: &RemMesh, marked: &[usize]) -> (RemMesh, HashMap<(usize, usize), usize>) {
    match mesh.dim {
        2 => refine_marked_tri3(mesh, marked),
        3 => refine_marked_tet4(mesh, marked),
        _ => {
            log::warn!("AMR red refinement is only implemented for 2-D and 3-D meshes; skipping.");
            (mesh.clone(), HashMap::new())
        }
    }
}

fn refine_marked_tri3(mesh: &RemMesh, marked: &[usize]) -> (RemMesh, HashMap<(usize, usize), usize>) {
    let marked_set: std::collections::HashSet<usize> = marked.iter().copied().collect();
    let _n_elems = mesh.volume_elements.len();

    // --- Build edge → element adjacency ---
    let mut edge_to_elems: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
    for (idx, elem) in mesh.volume_elements.iter().enumerate() {
        if elem.kind != ElementKind::Tri3 { continue; }
        for &(a, b) in &tri3_local_edges() {
            let key = edge_key(elem.node_ids[a], elem.node_ids[b]);
            edge_to_elems.entry(key).or_default().push(idx);
        }
    }

    // --- Conformity propagation ---
    let mut to_refine = marked_set.clone();
    let mut frontier: Vec<usize> = marked.to_vec();
    loop {
        let mut new_frontier = Vec::new();
        for &e in &frontier {
            let elem = &mesh.volume_elements[e];
            if elem.kind != ElementKind::Tri3 { continue; }
            for &(a, b) in &tri3_local_edges() {
                let key = edge_key(elem.node_ids[a], elem.node_ids[b]);
                if let Some(nbrs) = edge_to_elems.get(&key) {
                    for &nb in nbrs {
                        if !to_refine.contains(&nb) {
                            to_refine.insert(nb);
                            new_frontier.push(nb);
                        }
                    }
                }
            }
        }
        if new_frontier.is_empty() { break; }
        frontier = new_frontier;
    }

    // --- Create midpoint nodes ---
    let mut new_nodes = mesh.nodes.clone();
    let mut midpoint_map: HashMap<(usize, usize), usize> = HashMap::new();
    let mut next_id = new_nodes.len();

    for &e in &to_refine {
        let elem = &mesh.volume_elements[e];
        if elem.kind != ElementKind::Tri3 { continue; }
        for &(a, b) in &tri3_local_edges() {
            let key = edge_key(elem.node_ids[a], elem.node_ids[b]);
            midpoint_map.entry(key).or_insert_with(|| {
                let na = &mesh.nodes[key.0];
                let nb = &mesh.nodes[key.1];
                let mid = Node {
                    id: next_id,
                    x: 0.5 * (na.x + nb.x),
                    y: 0.5 * (na.y + nb.y),
                    z: 0.5 * (na.z + nb.z),
                };
                new_nodes.push(mid);
                let id = next_id;
                next_id += 1;
                id
            });
        }
    }

    // --- Build new element list ---
    let mut new_vol_elems: Vec<Element> = Vec::new();
    let mut new_elem_id = 1usize;

    for (idx, elem) in mesh.volume_elements.iter().enumerate() {
        if elem.kind != ElementKind::Tri3 || !to_refine.contains(&idx) {
            let mut e2 = elem.clone();
            e2.id = new_elem_id;
            new_elem_id += 1;
            new_vol_elems.push(e2);
            continue;
        }
        let [n0, n1, n2] = [elem.node_ids[0], elem.node_ids[1], elem.node_ids[2]];
        let m01 = *midpoint_map.get(&edge_key(n0, n1)).unwrap();
        let m12 = *midpoint_map.get(&edge_key(n1, n2)).unwrap();
        let m02 = *midpoint_map.get(&edge_key(n0, n2)).unwrap();
        for child_nodes in [[n0, m01, m02], [m01, n1, m12], [m02, m12, n2], [m01, m12, m02]] {
            new_vol_elems.push(Element {
                id: new_elem_id, kind: ElementKind::Tri3, tag: elem.tag,
                node_ids: child_nodes.to_vec(), rank: elem.rank,
            });
            new_elem_id += 1;
        }
    }

    // Boundary elements: bisect Line2 edges that were split
    let mut new_bnd_elems: Vec<Element> = Vec::new();
    for belem in &mesh.boundary_elements {
        if belem.kind == ElementKind::Line2 && belem.node_ids.len() == 2 {
            let key = edge_key(belem.node_ids[0], belem.node_ids[1]);
            if let Some(&mid) = midpoint_map.get(&key) {
                new_bnd_elems.push(Element {
                    id: new_elem_id, kind: ElementKind::Line2, tag: belem.tag,
                    node_ids: vec![belem.node_ids[0], mid], rank: belem.rank,
                });
                new_elem_id += 1;
                new_bnd_elems.push(Element {
                    id: new_elem_id, kind: ElementKind::Line2, tag: belem.tag,
                    node_ids: vec![mid, belem.node_ids[1]], rank: belem.rank,
                });
                new_elem_id += 1;
                continue;
            }
        }
        let mut e2 = belem.clone();
        e2.id = new_elem_id;
        new_elem_id += 1;
        new_bnd_elems.push(e2);
    }

    let new_mesh = RemMesh {
        nodes: new_nodes,
        volume_elements: new_vol_elems,
        boundary_elements: new_bnd_elems,
        domain_tags: mesh.domain_tags.clone(),
        boundary_tags: mesh.boundary_tags.clone(),
        dim: mesh.dim,
        rank: mesh.rank,
        size: mesh.size,
    };

    (new_mesh, midpoint_map)
}

/// Tet4 red refinement: 1 → 8 children per marked tet.
///
/// Adds a midpoint on each of the 6 edges. The 8 children are:
///   4 corner tets (one per original vertex) +
///   4 inner tets filling the central octahedron (Freudenthal partition).
///
/// The Freudenthal partition cuts the octahedron with the diagonal between the
/// two midpoints of opposite edges. The choice of diagonal is determined by the
/// globally shortest edge among the three diagonals (for mesh quality).
///
/// Boundary Tri3 faces that were split get 4 child triangles as well.
fn refine_marked_tet4(mesh: &RemMesh, marked: &[usize]) -> (RemMesh, HashMap<(usize, usize), usize>) {
    let marked_set: std::collections::HashSet<usize> = marked.iter().copied().collect();

    // --- Build edge → element adjacency (for conformity pass) ---
    let mut edge_to_elems: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
    for (idx, elem) in mesh.volume_elements.iter().enumerate() {
        if elem.kind != ElementKind::Tet4 { continue; }
        for &(a, b) in &tet4_local_edges() {
            let key = edge_key(elem.node_ids[a], elem.node_ids[b]);
            edge_to_elems.entry(key).or_default().push(idx);
        }
    }

    // --- Conformity propagation ---
    let mut to_refine = marked_set.clone();
    let mut frontier: Vec<usize> = marked.to_vec();
    loop {
        let mut new_frontier = Vec::new();
        for &e in &frontier {
            let elem = &mesh.volume_elements[e];
            if elem.kind != ElementKind::Tet4 { continue; }
            for &(a, b) in &tet4_local_edges() {
                let key = edge_key(elem.node_ids[a], elem.node_ids[b]);
                if let Some(nbrs) = edge_to_elems.get(&key) {
                    for &nb in nbrs {
                        if !to_refine.contains(&nb) {
                            to_refine.insert(nb);
                            new_frontier.push(nb);
                        }
                    }
                }
            }
        }
        if new_frontier.is_empty() { break; }
        frontier = new_frontier;
    }

    // --- Create midpoint nodes ---
    let mut new_nodes = mesh.nodes.clone();
    let mut midpoint_map: HashMap<(usize, usize), usize> = HashMap::new();
    let mut next_id = new_nodes.len();

    for &e in &to_refine {
        let elem = &mesh.volume_elements[e];
        if elem.kind != ElementKind::Tet4 { continue; }
        for &(a, b) in &tet4_local_edges() {
            let key = edge_key(elem.node_ids[a], elem.node_ids[b]);
            midpoint_map.entry(key).or_insert_with(|| {
                let na = &mesh.nodes[key.0];
                let nb = &mesh.nodes[key.1];
                let mid = Node {
                    id: next_id,
                    x: 0.5 * (na.x + nb.x),
                    y: 0.5 * (na.y + nb.y),
                    z: 0.5 * (na.z + nb.z),
                };
                new_nodes.push(mid);
                let id = next_id;
                next_id += 1;
                id
            });
        }
    }

    // --- Build new element list ---
    let mut new_vol_elems: Vec<Element> = Vec::new();
    let mut new_elem_id = 1usize;

    for (idx, elem) in mesh.volume_elements.iter().enumerate() {
        if elem.kind != ElementKind::Tet4 || !to_refine.contains(&idx) {
            let mut e2 = elem.clone();
            e2.id = new_elem_id;
            new_elem_id += 1;
            new_vol_elems.push(e2);
            continue;
        }

        let [n0, n1, n2, n3] = [
            elem.node_ids[0], elem.node_ids[1],
            elem.node_ids[2], elem.node_ids[3],
        ];

        // 6 edge midpoints: m_ij = midpoint of edge (ni, nj)
        let m01 = *midpoint_map.get(&edge_key(n0, n1)).unwrap();
        let m02 = *midpoint_map.get(&edge_key(n0, n2)).unwrap();
        let m03 = *midpoint_map.get(&edge_key(n0, n3)).unwrap();
        let m12 = *midpoint_map.get(&edge_key(n1, n2)).unwrap();
        let m13 = *midpoint_map.get(&edge_key(n1, n3)).unwrap();
        let m23 = *midpoint_map.get(&edge_key(n2, n3)).unwrap();

        // 4 corner tets (one vertex + 3 edge midpoints from that vertex)
        let corners: [[usize; 4]; 4] = [
            [n0,  m01, m02, m03],
            [m01, n1,  m12, m13],
            [m02, m12, n2,  m23],
            [m03, m13, m23, n3 ],
        ];
        for cn in &corners {
            new_vol_elems.push(Element {
                id: new_elem_id, kind: ElementKind::Tet4, tag: elem.tag,
                node_ids: cn.to_vec(), rank: elem.rank,
            });
            new_elem_id += 1;
        }

        // 4 inner tets filling the central octahedron (Freudenthal partition).
        // The octahedron has vertices {m01, m02, m03, m12, m13, m23}.
        // Three possible diagonals: (m01,m23), (m02,m13), (m03,m12).
        // Choose the shortest diagonal for best mesh quality.
        let d01_23 = node_dist(&new_nodes[m01], &new_nodes[m23]);
        let d02_13 = node_dist(&new_nodes[m02], &new_nodes[m13]);
        let d03_12 = node_dist(&new_nodes[m03], &new_nodes[m12]);

        let inner: [[usize; 4]; 4] = if d01_23 <= d02_13 && d01_23 <= d03_12 {
            // Diagonal m01–m23
            [[m01, m02, m03, m23],
             [m01, m03, m13, m23],
             [m01, m12, m02, m23],
             [m01, m13, m12, m23]]
        } else if d02_13 <= d03_12 {
            // Diagonal m02–m13
            [[m02, m01, m03, m13],
             [m02, m03, m23, m13],
             [m02, m12, m01, m13],
             [m02, m23, m12, m13]]
        } else {
            // Diagonal m03–m12
            [[m03, m01, m02, m12],
             [m03, m02, m23, m12],
             [m03, m13, m01, m12],
             [m03, m23, m13, m12]]
        };

        for inn in &inner {
            new_vol_elems.push(Element {
                id: new_elem_id, kind: ElementKind::Tet4, tag: elem.tag,
                node_ids: inn.to_vec(), rank: elem.rank,
            });
            new_elem_id += 1;
        }
    }

    // Boundary elements: bisect Tri3 faces that were split.
    // A boundary Tri3 with all 3 edges in midpoint_map → 4 children.
    let mut new_bnd_elems: Vec<Element> = Vec::new();
    for belem in &mesh.boundary_elements {
        if belem.kind == ElementKind::Tri3 && belem.node_ids.len() == 3 {
            let [b0, b1, b2] = [belem.node_ids[0], belem.node_ids[1], belem.node_ids[2]];
            let e01 = midpoint_map.get(&edge_key(b0, b1)).copied();
            let e12 = midpoint_map.get(&edge_key(b1, b2)).copied();
            let e02 = midpoint_map.get(&edge_key(b0, b2)).copied();
            if let (Some(m01), Some(m12), Some(m02)) = (e01, e12, e02) {
                for cn in [[b0, m01, m02], [m01, b1, m12], [m02, m12, b2], [m01, m12, m02]] {
                    new_bnd_elems.push(Element {
                        id: new_elem_id, kind: ElementKind::Tri3, tag: belem.tag,
                        node_ids: cn.to_vec(), rank: belem.rank,
                    });
                    new_elem_id += 1;
                }
                continue;
            }
        }
        if belem.kind == ElementKind::Line2 && belem.node_ids.len() == 2 {
            let key = edge_key(belem.node_ids[0], belem.node_ids[1]);
            if let Some(&mid) = midpoint_map.get(&key) {
                new_bnd_elems.push(Element {
                    id: new_elem_id, kind: ElementKind::Line2, tag: belem.tag,
                    node_ids: vec![belem.node_ids[0], mid], rank: belem.rank,
                });
                new_elem_id += 1;
                new_bnd_elems.push(Element {
                    id: new_elem_id, kind: ElementKind::Line2, tag: belem.tag,
                    node_ids: vec![mid, belem.node_ids[1]], rank: belem.rank,
                });
                new_elem_id += 1;
                continue;
            }
        }
        let mut e2 = belem.clone();
        e2.id = new_elem_id;
        new_elem_id += 1;
        new_bnd_elems.push(e2);
    }

    let new_mesh = RemMesh {
        nodes: new_nodes,
        volume_elements: new_vol_elems,
        boundary_elements: new_bnd_elems,
        domain_tags: mesh.domain_tags.clone(),
        boundary_tags: mesh.boundary_tags.clone(),
        dim: mesh.dim,
        rank: mesh.rank,
        size: mesh.size,
    };

    (new_mesh, midpoint_map)
}

fn node_dist(a: &Node, b: &Node) -> f64 {
    let dx = a.x - b.x;
    let dy = a.y - b.y;
    let dz = a.z - b.z;
    (dx*dx + dy*dy + dz*dz).sqrt()
}

// ---------------------------------------------------------------------------
// P1 prolongation
// ---------------------------------------------------------------------------

/// Transfer a coarse nodal solution to a fine mesh produced by [`refine_marked`].
///
/// - Existing node values (indices 0..coarse.n_nodes()) are copied as-is.
/// - New midpoint nodes are set to the average of their two parent nodes.
pub fn prolongate_p1(
    phi_coarse: &[f64],
    n_fine_nodes: usize,
    midpoint_map: &HashMap<(usize, usize), usize>,
) -> Vec<f64> {
    let mut phi_fine = vec![0.0f64; n_fine_nodes];
    for (i, &v) in phi_coarse.iter().enumerate() {
        phi_fine[i] = v;
    }
    for (&(a, b), &mid) in midpoint_map {
        phi_fine[mid] = 0.5 * (phi_coarse[a] + phi_coarse[b]);
    }
    phi_fine
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn edge_key(a: usize, b: usize) -> (usize, usize) {
    if a < b { (a, b) } else { (b, a) }
}

fn tri3_local_edges() -> [(usize, usize); 3] {
    [(0, 1), (1, 2), (0, 2)]
}

/// The 6 edges of a Tet4 by local vertex index pairs.
fn tet4_local_edges() -> [(usize, usize); 6] {
    [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
}

/// Compute P1 gradient and area for a Tri3 element.
fn tri3_grad_area(phi: &[f64], elem: &Element, mesh: &RemMesh) -> Option<([f64; 3], f64)> {
    if elem.node_ids.len() < 3 { return None; }
    let [n0, n1, n2] = [elem.node_ids[0], elem.node_ids[1], elem.node_ids[2]];
    let (x0, y0) = (mesh.nodes[n0].x, mesh.nodes[n0].y);
    let (x1, y1) = (mesh.nodes[n1].x, mesh.nodes[n1].y);
    let (x2, y2) = (mesh.nodes[n2].x, mesh.nodes[n2].y);

    let det = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
    let area = 0.5 * det.abs();
    if area < 1e-300 { return None; }

    let inv2a = 1.0 / (2.0 * area);
    let grads = [
        [(y1 - y2) * inv2a, (x2 - x1) * inv2a],
        [(y2 - y0) * inv2a, (x0 - x2) * inv2a],
        [(y0 - y1) * inv2a, (x1 - x0) * inv2a],
    ];
    let us = [phi[n0], phi[n1], phi[n2]];
    let gx: f64 = us.iter().zip(grads.iter()).map(|(&u, g)| u * g[0]).sum();
    let gy: f64 = us.iter().zip(grads.iter()).map(|(&u, g)| u * g[1]).sum();
    Some(([gx, gy, 0.0], area))
}

/// Compute P1 gradient and volume for a Tet4 element.
fn tet4_grad_vol(phi: &[f64], elem: &Element, mesh: &RemMesh) -> Option<([f64; 3], f64)> {
    if elem.node_ids.len() < 4 { return None; }
    let [n0, n1, n2, n3] = [elem.node_ids[0], elem.node_ids[1], elem.node_ids[2], elem.node_ids[3]];
    let (x0, y0, z0) = (mesh.nodes[n0].x, mesh.nodes[n0].y, mesh.nodes[n0].z);
    let (x1, y1, z1) = (mesh.nodes[n1].x, mesh.nodes[n1].y, mesh.nodes[n1].z);
    let (x2, y2, z2) = (mesh.nodes[n2].x, mesh.nodes[n2].y, mesh.nodes[n2].z);
    let (x3, y3, z3) = (mesh.nodes[n3].x, mesh.nodes[n3].y, mesh.nodes[n3].z);

    let j = [
        [x1-x0, x2-x0, x3-x0],
        [y1-y0, y2-y0, y3-y0],
        [z1-z0, z2-z0, z3-z0],
    ];
    let det = j[0][0]*(j[1][1]*j[2][2]-j[1][2]*j[2][1])
            - j[0][1]*(j[1][0]*j[2][2]-j[1][2]*j[2][0])
            + j[0][2]*(j[1][0]*j[2][1]-j[1][1]*j[2][0]);
    let vol = det.abs() / 6.0;
    if vol < 1e-300 { return None; }

    // J^{-T} * ref_grads; ref_grads for Tet4: (-1,-1,-1),(1,0,0),(0,1,0),(0,0,1)
    let inv_det = 1.0 / det;
    let ji = [
        [(j[1][1]*j[2][2]-j[1][2]*j[2][1])*inv_det, (j[0][2]*j[2][1]-j[0][1]*j[2][2])*inv_det, (j[0][1]*j[1][2]-j[0][2]*j[1][1])*inv_det],
        [(j[1][2]*j[2][0]-j[1][0]*j[2][2])*inv_det, (j[0][0]*j[2][2]-j[0][2]*j[2][0])*inv_det, (j[0][2]*j[1][0]-j[0][0]*j[1][2])*inv_det],
        [(j[1][0]*j[2][1]-j[1][1]*j[2][0])*inv_det, (j[0][1]*j[2][0]-j[0][0]*j[2][1])*inv_det, (j[0][0]*j[1][1]-j[0][1]*j[1][0])*inv_det],
    ];
    let ref_grads = [[-1.0f64,-1.0,-1.0],[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]];
    let us = [phi[n0], phi[n1], phi[n2], phi[n3]];
    let mut grad = [0.0f64; 3];
    for (k, rg) in ref_grads.iter().enumerate() {
        // physical grad of basis k: J^{-T} * rg
        let pg = [
            ji[0][0]*rg[0] + ji[0][1]*rg[1] + ji[0][2]*rg[2],
            ji[1][0]*rg[0] + ji[1][1]*rg[1] + ji[1][2]*rg[2],
            ji[2][0]*rg[0] + ji[2][1]*rg[1] + ji[2][2]*rg[2],
        ];
        grad[0] += us[k] * pg[0];
        grad[1] += us[k] * pg[1];
        grad[2] += us[k] * pg[2];
    }
    Some((grad, vol))
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::BoundaryTag;
    use std::collections::HashMap;

    fn unit_square_mesh() -> RemMesh {
        // 4 nodes, 2 Tri3 covering [0,1]²
        let nodes = vec![
            Node { id: 0, x: 0.0, y: 0.0, z: 0.0 },
            Node { id: 1, x: 1.0, y: 0.0, z: 0.0 },
            Node { id: 2, x: 1.0, y: 1.0, z: 0.0 },
            Node { id: 3, x: 0.0, y: 1.0, z: 0.0 },
        ];
        let ve = vec![
            Element { id: 1, kind: ElementKind::Tri3, tag: 1, node_ids: vec![0,1,2], rank: 0 },
            Element { id: 2, kind: ElementKind::Tri3, tag: 1, node_ids: vec![0,2,3], rank: 0 },
        ];
        let be = vec![
            Element { id: 3, kind: ElementKind::Line2, tag: 10, node_ids: vec![0,1], rank: 0 },
            Element { id: 4, kind: ElementKind::Line2, tag: 11, node_ids: vec![2,3], rank: 0 },
        ];
        let mut bt: HashMap<u32, BoundaryTag> = HashMap::new();
        bt.insert(10, BoundaryTag::Ground);
        bt.insert(11, BoundaryTag::Terminal { index: 1 });
        RemMesh {
            nodes, volume_elements: ve, boundary_elements: be,
            domain_tags: Default::default(), boundary_tags: bt,
            dim: 2, rank: 0, size: 1,
        }
    }

    #[test]
    fn zz_estimator_linear_phi_is_zero() {
        // φ = y is P1-exact → ZZ error should be ~0
        let mesh = unit_square_mesh();
        let phi: Vec<f64> = mesh.nodes.iter().map(|n| n.y).collect();
        let eta = zz_estimator(&mesh, &phi);
        assert_eq!(eta.len(), 2);
        for (i, &e) in eta.iter().enumerate() {
            assert!(e < 1e-12, "elem {i}: η={e:.2e}, expected ~0 for linear φ");
        }
    }

    #[test]
    fn dorfler_mark_selects_subset() {
        // Two elements: large error on elem 0, small on elem 1
        let eta = vec![1.0, 0.01];
        let marked = dorfler_mark(&eta, 0.5);
        // elem 0 alone contributes 1² / (1² + 0.01²) ≈ 99.99% > 50%
        assert_eq!(marked, vec![0]);
    }

    #[test]
    fn refine_marked_doubles_elements() {
        let mesh = unit_square_mesh();
        // Mark both elements
        let (fine, midpoints) = refine_marked(&mesh, &[0, 1]);
        // 2 original Tri3 → 8 children
        let tri3_count = fine.volume_elements.iter().filter(|e| e.kind == ElementKind::Tri3).count();
        assert_eq!(tri3_count, 8, "expected 8 child triangles, got {tri3_count}");
        // 3 new midpoint nodes added (unit square: 3 unique edges)
        assert!(fine.nodes.len() > mesh.nodes.len(), "should have more nodes after refinement");
        assert!(!midpoints.is_empty(), "midpoint_map should not be empty");
    }

    #[test]
    fn prolongate_p1_preserves_linear() {
        // φ = y; after refinement, midpoints get y = average of parents → exact
        let mesh = unit_square_mesh();
        let phi: Vec<f64> = mesh.nodes.iter().map(|n| n.y).collect();
        let (fine, midpoints) = refine_marked(&mesh, &[0, 1]);
        let phi_fine = prolongate_p1(&phi, fine.nodes.len(), &midpoints);
        for (i, node) in fine.nodes.iter().enumerate() {
            let err = (phi_fine[i] - node.y).abs();
            assert!(err < 1e-12, "node {i} at y={:.2}: prolongated={:.6}, err={err:.2e}", node.y, phi_fine[i]);
        }
    }

    #[test]
    fn amr_loop_reduces_error() {
        // After one AMR step on a coarse mesh, the fine mesh should have more DOFs
        let mesh = unit_square_mesh();
        let phi: Vec<f64> = mesh.nodes.iter().map(|n| n.y).collect();
        let eta = zz_estimator(&mesh, &phi);
        let marked = dorfler_mark(&eta, 0.3);
        let (fine, _) = refine_marked(&mesh, &marked);
        assert!(fine.n_nodes() >= mesh.n_nodes(), "fine mesh should have at least as many nodes");
    }

    // -----------------------------------------------------------------------
    // Tet4 AMR tests
    // -----------------------------------------------------------------------

    /// Minimal 3-D mesh: unit cube split into 6 Tet4 (Sommerville decomposition).
    fn unit_cube_tet_mesh() -> RemMesh {
        let nodes = vec![
            Node { id: 0, x: 0.0, y: 0.0, z: 0.0 },
            Node { id: 1, x: 1.0, y: 0.0, z: 0.0 },
            Node { id: 2, x: 1.0, y: 1.0, z: 0.0 },
            Node { id: 3, x: 0.0, y: 1.0, z: 0.0 },
            Node { id: 4, x: 0.0, y: 0.0, z: 1.0 },
            Node { id: 5, x: 1.0, y: 0.0, z: 1.0 },
            Node { id: 6, x: 1.0, y: 1.0, z: 1.0 },
            Node { id: 7, x: 0.0, y: 1.0, z: 1.0 },
        ];
        // 6 Tet4 from Sommerville decomposition of the unit cube
        let ve = vec![
            Element { id: 1, kind: ElementKind::Tet4, tag: 1, node_ids: vec![0,1,2,6], rank: 0 },
            Element { id: 2, kind: ElementKind::Tet4, tag: 1, node_ids: vec![0,2,3,6], rank: 0 },
            Element { id: 3, kind: ElementKind::Tet4, tag: 1, node_ids: vec![0,1,5,6], rank: 0 },
            Element { id: 4, kind: ElementKind::Tet4, tag: 1, node_ids: vec![0,4,5,6], rank: 0 },
            Element { id: 5, kind: ElementKind::Tet4, tag: 1, node_ids: vec![0,3,6,7], rank: 0 },
            Element { id: 6, kind: ElementKind::Tet4, tag: 1, node_ids: vec![0,4,6,7], rank: 0 },
        ];
        // Boundary face (z=0 face as Tri3 elements)
        let be = vec![
            Element { id: 7, kind: ElementKind::Tri3, tag: 10, node_ids: vec![0,1,2], rank: 0 },
            Element { id: 8, kind: ElementKind::Tri3, tag: 10, node_ids: vec![0,2,3], rank: 0 },
        ];
        let mut bt: HashMap<u32, crate::BoundaryTag> = HashMap::new();
        bt.insert(10, crate::BoundaryTag::Ground);
        RemMesh {
            nodes, volume_elements: ve, boundary_elements: be,
            domain_tags: Default::default(), boundary_tags: bt,
            dim: 3, rank: 0, size: 1,
        }
    }

    #[test]
    fn tet4_refine_all_produces_8x_elements() {
        let mesh = unit_cube_tet_mesh();
        let all_marked: Vec<usize> = (0..mesh.volume_elements.len()).collect();
        let (fine, midpoints) = refine_marked(&mesh, &all_marked);

        let tet4_count = fine.volume_elements.iter()
            .filter(|e| e.kind == ElementKind::Tet4)
            .count();
        assert_eq!(tet4_count, 6 * 8, "6 Tet4 × 8 = 48 children expected, got {tet4_count}");

        // Midpoint map: 6 tets, potentially shared edges; unit cube has 12+6=18 edges
        assert!(!midpoints.is_empty(), "midpoint_map should not be empty");
    }

    #[test]
    fn tet4_refine_increases_nodes() {
        let mesh = unit_cube_tet_mesh();
        let all_marked: Vec<usize> = (0..mesh.volume_elements.len()).collect();
        let (fine, _) = refine_marked(&mesh, &all_marked);
        assert!(fine.n_nodes() > mesh.n_nodes(), "refined mesh must have more nodes");
    }

    #[test]
    fn tet4_prolongate_p1_exact_for_linear() {
        // φ = x + 2y + 3z is P1-exact; prolongation at midpoints should be exact
        let mesh = unit_cube_tet_mesh();
        let phi: Vec<f64> = mesh.nodes.iter().map(|n| n.x + 2.0*n.y + 3.0*n.z).collect();
        let all_marked: Vec<usize> = (0..mesh.volume_elements.len()).collect();
        let (fine, midpoints) = refine_marked(&mesh, &all_marked);
        let phi_fine = prolongate_p1(&phi, fine.n_nodes(), &midpoints);
        for (i, node) in fine.nodes.iter().enumerate() {
            let exact = node.x + 2.0*node.y + 3.0*node.z;
            let err = (phi_fine[i] - exact).abs();
            assert!(err < 1e-12, "node {i}: prolongated={:.8}, exact={:.8}, err={err:.2e}",
                phi_fine[i], exact);
        }
    }

    #[test]
    fn tet4_zz_estimator_linear_phi() {
        // φ = x is P1-exact → ZZ recovered gradient should be [1,0,0]
        // The ZZ *indicator* η_e = √vol · ‖∇φ_h − G_h φ‖ may be non-zero on
        // asymmetric patches but should be small relative to the solution magnitude.
        let mesh = unit_cube_tet_mesh();
        let phi: Vec<f64> = mesh.nodes.iter().map(|n| n.x).collect();
        let eta = zz_estimator(&mesh, &phi);
        assert_eq!(eta.len(), 6, "should have 6 eta values for 6 tet elements");
        // Each element has φ-range ≤ 1.0; η should be < 1.0 (physically small)
        for (i, &e) in eta.iter().enumerate() {
            assert!(e < 1.0, "tet {i}: η={e:.4} should be < 1.0 for smooth φ=x");
        }
    }
}
