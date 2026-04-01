/// Dirichlet boundary condition application.
///
/// Strategy: symmetric elimination.
///  For DOF `d` with prescribed value `v`:
///    1. For every free row `i`: rhs[i] -= K[i,d] * v, then K[i,d] = 0
///    2. Zero row `d`: K[d,j] = 0 for j ≠ d, K[d,d] = 1
///    3. rhs[d] = v

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
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use rem_core::{TripletMatrix, solve_pcg};

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

        let result = solve_pcg(&mat, &rhs, 1e-12, 100);
        assert!(result.converged);
        assert!((result.solution[0] - 0.0).abs() < 1e-10, "φ[0]={}", result.solution[0]);
        assert!((result.solution[1] - 0.5).abs() < 1e-10, "φ[1]={}", result.solution[1]);
        assert!((result.solution[2] - 1.0).abs() < 1e-10, "φ[2]={}", result.solution[2]);
    }
}
