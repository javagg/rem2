/// Dirichlet boundary condition application.
///
/// Strategy: symmetric elimination.
///  For DOF `d` with prescribed value `v`:
///    1. For every free row `i`: rhs[i] -= K[i,d] * v, then K[i,d] = 0
///    2. Zero row `d`: K[d,j] = 0 for j ≠ d, K[d,d] = 1
///    3. rhs[d] = v

use rem_config::PalaceConfig;
use rem_core::CsrMatrix;
use rem_mesh::{RemMesh, BoundaryTag};
use std::collections::HashMap;

/// Identifies all Dirichlet DOFs from the mesh boundary conditions.
///
/// Returns a `HashMap<node_index → prescribed_value>` where:
/// - PEC and Ground nodes → 0.0
/// - Terminal with `index == excited_index` → `excitation_val`, others → 0.0
/// - LumpedPort with `index == excited_index` → `excitation_val`, others → 0.0
///
/// `excited_index`: the Terminal/LumpedPort **index** (not physical tag) to excite.
pub fn collect_dirichlet_dofs(
    mesh: &RemMesh,
    excited_index: Option<u32>,
    excitation_val: f64,
) -> HashMap<usize, f64> {
    let mut dofs: HashMap<usize, f64> = HashMap::new();

    for belem in &mesh.boundary_elements {
        if mesh.size > 1 && belem.rank != mesh.rank {
            continue;
        }
        let bc = match mesh.boundary_tags.get(&belem.tag) {
            Some(b) => b,
            None => continue,
        };

        let val = match bc {
            BoundaryTag::Pec | BoundaryTag::Ground => 0.0,
            BoundaryTag::Terminal { index } => {
                if Some(*index) == excited_index { excitation_val } else { 0.0 }
            }
            BoundaryTag::LumpedPort { index, .. } => {
                if Some(*index) == excited_index { excitation_val } else { 0.0 }
            }
            // WavePort (TEM approximation): treat as LumpedPort — apply Dirichlet φ=V on port face.
            // This is exact for TEM modes and a reasonable approximation for quasi-TEM modes.
            // Full TE/TM modal field matching is deferred to a future release.
            BoundaryTag::WavePort { index } => {
                if Some(*index) == excited_index { excitation_val } else { 0.0 }
            }
            _ => continue,
        };

        for &nid in &belem.node_ids {
            dofs.entry(nid).or_insert(val);
        }
    }

    dofs
}

/// Apply Dirichlet BCs to the assembled stiffness matrix and RHS vector.
///
/// Uses symmetric elimination so the resulting system remains SPD.
/// `dofs`: map of node-index → prescribed value.
/// Identifies Dirichlet DOFs for **open-circuit Z-matrix extraction**.
///
/// Unlike `collect_dirichlet_dofs`, this function does NOT set non-excited
/// LumpedPort / WavePort / Terminal nodes to φ=0.  Unexcited ports are left
/// with the natural Neumann BC (∂φ/∂n = 0 → no current, i.e. open circuit),
/// so the solution φ carries the open-circuit port voltage needed for Z_ij.
///
/// Only PEC/Ground and the single excited port receive Dirichlet values.
pub fn collect_dirichlet_dofs_open_circuit(
    mesh: &RemMesh,
    excited_index: Option<u32>,
    excitation_val: f64,
) -> HashMap<usize, f64> {
    let mut dofs: HashMap<usize, f64> = HashMap::new();

    for belem in &mesh.boundary_elements {
        if mesh.size > 1 && belem.rank != mesh.rank {
            continue;
        }
        let bc = match mesh.boundary_tags.get(&belem.tag) {
            Some(b) => b,
            None => continue,
        };

        match bc {
            BoundaryTag::Pec | BoundaryTag::Ground => {
                for &nid in &belem.node_ids {
                    dofs.entry(nid).or_insert(0.0);
                }
            }
            BoundaryTag::Terminal { index }
            | BoundaryTag::LumpedPort { index, .. }
            | BoundaryTag::WavePort { index } => {
                if Some(*index) == excited_index {
                    for &nid in &belem.node_ids {
                        dofs.entry(nid).or_insert(excitation_val);
                    }
                }
                // Non-excited ports: no entry → natural BC (open circuit)
            }
            _ => {}
        }
    }
    dofs
}

pub fn apply_dirichlet(
    mat: &mut CsrMatrix,
    rhs: &mut Vec<f64>,
    dofs: &HashMap<usize, f64>,
) {
    // Step 1: modify RHS for non-Dirichlet rows (subtract K[i,d]*v)
    // We iterate over the Dirichlet DOFs and zero out their column entries
    // in all other rows.
    for (&d, &val) in dofs {
        for row in 0..mat.nrows {
            if dofs.contains_key(&row) {
                continue; // will handle in step 2
            }
            let k_id = mat.zero_col_entry(row, d);
            if k_id.abs() > 0.0 {
                rhs[row] -= k_id * val;
            }
        }
    }

    // Step 2: for Dirichlet DOF d, preserve its original diagonal so all equations stay
    // on the same scale (O(ε·area·grad²)).  Setting K[d,d]=1.0 would create a 1e12
    // mismatch vs free rows (K~1e-14), causing PCG to falsely converge in 1 iteration.
    for (&d, &val) in dofs {
        let k_diag = mat.diagonal_entry(d);
        let effective = if k_diag.abs() > 1e-300 { k_diag } else { 1.0 };
        mat.zero_row_set_diag(d, effective);
        rhs[d] = effective * val;
    }
}

// ---------------------------------------------------------------------------
// Periodic (Floquet) boundary conditions
// ---------------------------------------------------------------------------

/// Collect (donor_node, receiver_node) pairs for periodic boundaries.
///
/// For each `Boundaries.Periodic` spec with a `Translation` vector,
/// this function iterates over donor boundary faces and matches each donor
/// node to the receiver node closest to `donor_pos + translation`.
///
/// **Γ-point only**: when `FloquetWaveVector` is zero (or absent) the BC
/// reduces to a standard periodic (mirror) constraint  φ[recv] = φ[donor],
/// handled by folding `recv` into the Dirichlet map.
///
/// If any `FloquetWaveVector` entry is non-zero the function logs a warning
/// and returns an empty vec — complex phase-shift BCs require a complex solver
/// and are not yet supported.
pub fn collect_periodic_node_pairs(
    mesh: &RemMesh,
    config: &PalaceConfig,
) -> Vec<(usize, usize)> {
    let mut pairs: Vec<(usize, usize)> = Vec::new();

    for periodic in &config.boundaries.periodic {
        // Check for complex Floquet wave vector
        let k_nonzero = periodic.floquet_wave_vector.iter().any(|&v| v.abs() > 1e-14);
        if k_nonzero {
            log::warn!(
                "[REM] Floquet periodic BC: non-zero FloquetWaveVector — complex phase-shift BCs \
                 are not yet supported. Ignoring periodic constraint for this pair. \
                 Results for photonic crystals / metamaterials will be approximate."
            );
            continue;
        }

        for pair in &periodic.boundary_pairs {
            // Translation vector (default zero = mirror symmetry)
            let tx = pair.translation.first().copied().unwrap_or(0.0) * config.model.l0;
            let ty = pair.translation.get(1).copied().unwrap_or(0.0) * config.model.l0;
            let tz = pair.translation.get(2).copied().unwrap_or(0.0) * config.model.l0;

            // Collect receiver nodes (sorted into a spatial lookup)
            let mut recv_nodes: Vec<usize> = Vec::new();
            for belem in &mesh.boundary_elements {
                if pair.receiver_attributes.contains(&belem.tag) {
                    for &nid in &belem.node_ids {
                        if !recv_nodes.contains(&nid) {
                            recv_nodes.push(nid);
                        }
                    }
                }
            }

            // For each donor node, find the receiver node nearest to donor + translation
            for belem in &mesh.boundary_elements {
                if !pair.donor_attributes.contains(&belem.tag) {
                    continue;
                }
                for &dnid in &belem.node_ids {
                    let n = &mesh.nodes[dnid];
                    let tx_x = n.x + tx;
                    let tx_y = n.y + ty;
                    let tx_z = n.z + tz;

                    let closest = recv_nodes.iter().copied().min_by(|&a, &b| {
                        let na = &mesh.nodes[a];
                        let nb = &mesh.nodes[b];
                        let da = (na.x - tx_x).powi(2) + (na.y - tx_y).powi(2) + (na.z - tx_z).powi(2);
                        let db = (nb.x - tx_x).powi(2) + (nb.y - tx_y).powi(2) + (nb.z - tx_z).powi(2);
                        da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
                    });

                    if let Some(rnid) = closest {
                        // Only add if the receiver node is genuinely close (within 1% of cell size)
                        let rn = &mesh.nodes[rnid];
                        let dist2 = (rn.x - tx_x).powi(2) + (rn.y - tx_y).powi(2) + (rn.z - tx_z).powi(2);
                        // Use a generous tolerance: 1 nm or 0.1% of translation magnitude
                        let tol_sq = {
                            let tmag = (tx * tx + ty * ty + tz * tz).sqrt();
                            let tol = if tmag > 1e-12 { tmag * 1e-3 } else { 1e-9 };
                            tol * tol
                        };
                        if dist2 < tol_sq && dnid != rnid {
                            if !pairs.iter().any(|&(d, r)| d == dnid && r == rnid) {
                                pairs.push((dnid, rnid));
                            }
                        }
                    }
                }
            }
        }
    }

    log::info!("[REM] Periodic BC: {} node pairs matched (Γ-point, real constraint)", pairs.len());
    pairs
}

/// Apply periodic (Γ-point) constraints φ[recv] = φ[donor].
///
/// **Strategy**: before this function is called, the stiffness / mass
/// triplet matrices should have been remapped with
/// `TripletMatrix::remap_periodic_nodes(pairs)`, which folds all receiver
/// contributions into the donor DOF.  This function then:
///   1. Inserts each receiver into `dofs` with value 0.0, so that
///      `apply_dirichlet` subsequently zeros the recv row and sets K[recv,recv]=1.
///   2. After the solve, the caller must copy  φ[recv] = φ[donor].
///
/// If the donor is already a Dirichlet DOF (e.g. PEC φ=0), the receiver
/// inherits the same prescribed value (both are constrained; merge is still
/// correct because the receiver row carries no free contributions after remapping).
pub fn apply_periodic(
    dofs: &mut HashMap<usize, f64>,
    pairs: &[(usize, usize)],
) {
    for &(donor, recv) in pairs {
        if dofs.contains_key(&recv) {
            continue; // already constrained — Dirichlet wins
        }
        let donor_val = dofs.get(&donor).copied().unwrap_or(0.0);
        // Mark receiver as Dirichlet = 0 (all coupling already remapped to donor).
        // The actual solution value will be restored to φ[donor] post-solve.
        let _ = donor_val;
        dofs.insert(recv, 0.0);
    }
}

/// After solving, propagate periodic DOF values from donor to receiver.
///
/// Call this once the solver has produced `phi[donor]`.
pub fn propagate_periodic(phi: &mut Vec<f64>, pairs: &[(usize, usize)]) {
    for &(donor, recv) in pairs {
        if recv < phi.len() && donor < phi.len() {
            phi[recv] = phi[donor];
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use rem_core::{TripletMatrix, solve_pcg};
    use rem_parallel::NoComm;

    /// Build 3×3 system [[2,-1,0],[-1,2,-1],[0,-1,2]] and apply φ[0]=0, φ[2]=1.
    /// Exact solution: φ = [0, 0.5, 1].
    #[test]
    fn dirichlet_1d_bar() {
        let mut t = TripletMatrix::new(3, 3);
        t.add(0, 0, 2.0); t.add(0, 1, -1.0);
        t.add(1, 0, -1.0); t.add(1, 1, 2.0); t.add(1, 2, -1.0);
        t.add(2, 1, -1.0); t.add(2, 2, 2.0);
        let mut mat = t.to_csr();
        let mut rhs = vec![0.0f64; 3];

        let dofs: HashMap<usize, f64> = [(0, 0.0), (2, 1.0)].iter().copied().collect();
        apply_dirichlet(&mut mat, &mut rhs, &dofs);

        let result = solve_pcg(&mat, &rhs, 1e-12, 100, &NoComm);
        assert!(result.converged);
        assert!((result.solution[0] - 0.0).abs() < 1e-10, "φ[0]={}", result.solution[0]);
        assert!((result.solution[1] - 0.5).abs() < 1e-10, "φ[1]={}", result.solution[1]);
        assert!((result.solution[2] - 1.0).abs() < 1e-10, "φ[2]={}", result.solution[2]);
    }

    /// WavePort BC: identical DOF values to LumpedPort with same index.
    #[test]
    fn waveport_treated_as_dirichlet() {
        use rem_mesh::{Node, Element, ElementKind, RemMesh};
        use std::collections::HashMap;

        let nodes = vec![
            Node { id: 0, x: 0.0, y: 0.0, z: 0.0 },
            Node { id: 1, x: 1.0, y: 0.0, z: 0.0 },
            Node { id: 2, x: 1.0, y: 1.0, z: 0.0 },
            Node { id: 3, x: 0.0, y: 1.0, z: 0.0 },
        ];
        let volume_elements = vec![
            Element { id: 1, kind: ElementKind::Tri3, tag: 1, node_ids: vec![0, 1, 2], rank: 0 },
            Element { id: 2, kind: ElementKind::Tri3, tag: 1, node_ids: vec![0, 2, 3], rank: 0 },
        ];
        let boundary_elements = vec![
            Element { id: 3, kind: ElementKind::Line2, tag: 10, node_ids: vec![0, 1], rank: 0 },
            Element { id: 4, kind: ElementKind::Line2, tag: 11, node_ids: vec![2, 3], rank: 0 },
        ];
        let mut boundary_tags: HashMap<u32, rem_mesh::BoundaryTag> = HashMap::new();
        boundary_tags.insert(10, rem_mesh::BoundaryTag::Ground);
        boundary_tags.insert(11, rem_mesh::BoundaryTag::WavePort { index: 1 });

        let mesh = RemMesh {
            nodes, volume_elements, boundary_elements,
            domain_tags: Default::default(), boundary_tags,
            dim: 2, rank: 0, size: 1,
        };

        let dofs = collect_dirichlet_dofs(&mesh, Some(1), 1.0);

        // Bottom (Ground) nodes 0,1 → 0.0; top (WavePort) nodes 2,3 → 1.0
        assert_eq!(dofs.get(&0), Some(&0.0), "node 0 should be Ground=0");
        assert_eq!(dofs.get(&1), Some(&0.0), "node 1 should be Ground=0");
        assert_eq!(dofs.get(&2), Some(&1.0), "node 2 should be WavePort=1");
        assert_eq!(dofs.get(&3), Some(&1.0), "node 3 should be WavePort=1");
    }
}
