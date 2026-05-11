//! HCurl (Nedelec edge-element) magnetostatic solver for 3-D problems.
//!
//! Solves the vector magnetic potential equation:
//!   curl(ν curl A) = J_s    in Ω
//!   n × A = 0               on Γ_PEC/Ground (tangential A vanishes)
//!
//! Discretisation: Nedelec ND1/ND2 edge elements (H(curl) conforming).
//!
//! ## Gauge regularisation
//!
//! The curl-curl operator K has a non-trivial null space (gradient fields:
//! K ∇φ = 0 for any φ).  To make the system uniquely solvable we add a
//! small diagonal shift: K_reg = K + ε_gauge · I, with ε_gauge = 1e-10.
//! For physical problems where div J = 0 the regularisation does not affect
//! the physical solution (B = curl A is unchanged by a gradient correction
//! of order ε_gauge).
//!
//! ## Limitations (first implementation)
//!
//! - Surface-current source is modelled as a Neumann contribution on the
//!   excited boundary: f_i = |edge_i| for each edge on the SurfaceCurrent tag.
//! - B = curl A post-processing is returned as zero-length (not yet wired
//!   through the scalar output path); the raw edge-DOF vector A_h is stored
//!   in the result.
//! - 2-D problems still use the scalar H1 path (A_z) which is physically
//!   correct; HCurl is only meaningful for 3-D magnetic vector potential.

use fem_assembly::coefficient::CtxFnCoeff;
use fem_assembly::standard::CurlCurlIntegrator;
use fem_assembly::{
    VectorAssembler, VectorBoundaryAssembler, VectorBoundaryLinearIntegrator, VectorBdQpData,
};
use fem_space::{boundary_dofs_hcurl, FESpace, HCurlSpace};
use rem_config::PalaceConfig;
use rem_core::{RemResult, solve_pcg};
use rem_materials::DomainMap;
use rem_mesh::{BoundaryTag, RemMesh};
use rem_parallel::Comm;

use crate::MagnetostaticResult;

/// Entry point for 3-D HCurl magnetostatic solve.
///
/// Returns a `MagnetostaticResult` whose `a_vec` contains the Nedelec
/// edge-DOF vector A_h (length = n_edges).
///
/// ## Solver strategy
///
/// The curl-curl operator is symmetric positive semi-definite with a
/// non-trivial nullspace (gradient fields).  We add a small gauge
/// regularisation ε·I (ε = 1e-10) to make the system positive-definite, then
/// solve with SSOR-preconditioned CG (real arithmetic).  SSOR-CG converges
/// reliably for the near-singular curl-curl system whereas the Jacobi-BiCGSTAB
/// solver historically used by `solve_pcg_complex` often diverges.
pub fn run_hcurl_3d(
    config: &PalaceConfig,
    mesh: &RemMesh,
    domain_map: &DomainMap,
    comm: &dyn Comm,
) -> RemResult<MagnetostaticResult> {
    let order = config.solver.order.clamp(1, 2) as u8;
    if config.solver.order > 2 {
        log::warn!(
            "HCurl magnetostatic supports order 1/2 (ND1/ND2); Solver.Order={} → using order 2.",
            config.solver.order
        );
    }
    log::info!("HCurl magnetostatic (3-D), ND{}", order);

    // ── Build HCurl space ────────────────────────────────────────────────────
    let simplex = mesh.to_simplex_mesh();
    let space   = HCurlSpace::new(simplex, order);
    let n       = space.n_dofs();
    log::info!("  Edge DOFs: {}", n);

    // ── Assemble curl-curl stiffness K ───────────────────────────────────────
    let inv_mu = CtxFnCoeff(|ctx: &fem_assembly::coefficient::CoeffCtx<'_>| {
        domain_map.get(ctx.elem_tag as u32).reluctivity()
    });
    let curl_curl = CurlCurlIntegrator { mu: inv_mu };
    let k_fem = VectorAssembler::assemble_bilinear(&space, &[&curl_curl], 4);
    let mut k_csr = rem_core::CsrMatrix::from_fem_csr(k_fem);

    // ── Gauge regularisation: K_reg = K + ε_gauge · I ────────────────────────
    const EPS_GAUGE: f64 = 1e-10;
    for i in 0..n {
        add_to_diagonal(&mut k_csr, i, EPS_GAUGE);
    }

    // ── RHS: surface-current source ──────────────────────────────────────────
    //   f_i = ∫_Γ J_s · φ_i dS   over SurfaceCurrent boundaries
    let mut rhs = vec![0.0f64; n];
    build_surface_current_rhs(config, mesh, &space, n, &mut rhs);

    // ── Dirichlet BCs: PEC / Ground → zero tangential A ─────────────────────
    let pec_tags: Vec<i32> = mesh
        .boundary_tags
        .iter()
        .filter_map(|(tag, bc)| {
            if matches!(bc, BoundaryTag::Pec | BoundaryTag::Ground) {
                Some(*tag as i32)
            } else {
                None
            }
        })
        .collect();
    let pec_dofs: Vec<usize> = boundary_dofs_hcurl(space.mesh(), &space, &pec_tags)
        .into_iter()
        .map(|d| d as usize)
        .collect();
    log::info!("  PEC-constrained edge DOFs: {}", pec_dofs.len());

    apply_zero_dirichlet(&mut k_csr, &mut rhs, &pec_dofs);

    // ── Solve (real SSOR-PCG) ────────────────────────────────────────────────
    let lin = &config.solver.linear;
    let result = solve_pcg(&k_csr, &rhs, lin.tol, lin.max_iter, comm);
    if result.converged {
        log::info!("  SSOR-PCG converged in {} iterations (|r|={:.2e})", result.iterations, result.residual_norm);
    } else {
        log::warn!("  SSOR-PCG did not converge after {} iterations (|r|={:.2e})", result.iterations, result.residual_norm);
    }
    let a_h = if result.solution.is_empty() {
        vec![0.0f64; n]
    } else {
        result.solution
    };

    // ── Energy = ½ x^T K_curl x ──────────────────────────────────────────────
    let energy = compute_magnetic_energy(&k_csr, &a_h, n, comm);
    log::info!("  Magnetic energy W_m = {:.4e} J", energy);

    // ── B-field (curl A) at element centroids ─────────────────────────────────
    let b_raw = fem_assembly::postprocess::compute_element_curl(&space, &a_h);
    let b_field: Vec<[f64; 3]> = b_raw
        .iter()
        .map(|v| {
            let mut arr = [0.0f64; 3];
            for (d, &val) in v.iter().enumerate().take(3) {
                arr[d] = val;
            }
            arr
        })
        .collect();

    Ok(MagnetostaticResult {
        a_vec: a_h,
        b_field,
        energy,
        is_hcurl: true,
    })
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Add `val` to the diagonal entry of row `i` in a real CSR matrix.
fn add_to_diagonal(k: &mut rem_core::CsrMatrix, i: usize, val: f64) {
    for p in k.row_ptr[i]..k.row_ptr[i + 1] {
        if k.col_idx[p] == i {
            k.values[p] += val;
            return;
        }
    }
}

/// Zero out rows and columns of Dirichlet DOFs in a real CSR matrix,
/// set diagonal = 1, rhs = 0.
fn apply_zero_dirichlet(
    k: &mut rem_core::CsrMatrix,
    rhs: &mut Vec<f64>,
    dofs: &[usize],
) {
    let dof_set: std::collections::HashSet<usize> = dofs.iter().copied().collect();
    let n = k.nrows;

    // Zero rows of constrained DOFs, set diagonal = 1
    for &d in dofs {
        if d >= n { continue; }
        for p in k.row_ptr[d]..k.row_ptr[d + 1] {
            k.values[p] = if k.col_idx[p] == d { 1.0 } else { 0.0 };
        }
        rhs[d] = 0.0;
    }

    // Zero columns of constrained DOFs in un-constrained rows
    for row in 0..n {
        if dof_set.contains(&row) { continue; }
        for p in k.row_ptr[row]..k.row_ptr[row + 1] {
            if dof_set.contains(&k.col_idx[p]) {
                k.values[p] = 0.0;
            }
        }
    }
}

// ─── SurfaceCurrent integrator ────────────────────────────────────────────────

/// Boundary linear integrator: f_i = ∫_Γ J_s · φ_i dS
///
/// The surface current density J_s is treated as spatially constant over
/// each boundary face (the value at the face centroid is used).
struct SurfaceCurrentIntegrator {
    js: [f64; 3],
}

impl VectorBoundaryLinearIntegrator for SurfaceCurrentIntegrator {
    fn add_to_face_vector(&self, qp: &VectorBdQpData<'_>, f_face: &mut [f64]) {
        let n = qp.n_dofs;
        let dim = qp.dim;
        for i in 0..n {
            let mut s = 0.0;
            for c in 0..dim {
                s += self.js[c] * qp.phi_vec[i * dim + c];
            }
            f_face[i] += qp.weight * s;
        }
    }
}

/// Build the RHS vector from SurfaceCurrent boundaries using the
/// H(curl) boundary linear form: f_i = ∫_Γ J_s · φ_i dS.
fn build_surface_current_rhs<S: FESpace>(
    config: &PalaceConfig,
    mesh: &RemMesh,
    space: &S,
    _n: usize,
    rhs: &mut Vec<f64>,
) {
    // Collect SurfaceCurrent boundary tags grouped by config index.
    let mut sc_by_index: std::collections::HashMap<u32, Vec<i32>> =
        std::collections::HashMap::new();
    for (&tag, bc) in &mesh.boundary_tags {
        if let BoundaryTag::SurfaceCurrent { index } = bc {
            sc_by_index.entry(*index).or_default().push(tag as i32);
        }
    }

    if sc_by_index.is_empty() {
        log::warn!("HCurl magnetostatic: no SurfaceCurrent boundary found; RHS = 0.");
        return;
    }

    // Build index → direction map from config.
    let mut dir_map: std::collections::HashMap<u32, [f64; 3]> =
        std::collections::HashMap::new();
    for sc in &config.boundaries.surface_current {
        dir_map.insert(sc.index, parse_direction(&sc.direction));
    }

    // Assemble boundary linear form for each SurfaceCurrent index.
    for (index, tags) in &sc_by_index {
        let js = dir_map.get(index).copied().unwrap_or_else(|| {
            log::warn!(
                "HCurl magnetostatic: SurfaceCurrent index {} not found in config, using +Z",
                index
            );
            [0.0, 0.0, 1.0]
        });
        let integrator = SurfaceCurrentIntegrator { js };
        let f = VectorBoundaryAssembler::assemble_boundary_linear(space, &[&integrator], tags, 4);
        for (i, val) in f.iter().enumerate() {
            rhs[i] += val;
        }
    }
}

/// Parse a direction string (from `SurfaceCurrentSpec.direction`) into a unit vector.
///
/// Accepted formats: `"+X"`, `"-X"`, `"+Y"`, `"-Y"`, `"+Z"`, `"-Z"` (case-insensitive).
fn parse_direction(s: &str) -> [f64; 3] {
    match s.trim().to_uppercase().as_str() {
        "+X" | "X" => [1.0, 0.0, 0.0],
        "-X" => [-1.0, 0.0, 0.0],
        "+Y" | "Y" => [0.0, 1.0, 0.0],
        "-Y" => [0.0, -1.0, 0.0],
        "+Z" | "Z" => [0.0, 0.0, 1.0],
        "-Z" => [0.0, 0.0, -1.0],
        _ => {
            log::warn!(
                "HCurl magnetostatic: unrecognised SurfaceCurrent direction '{s}', defaulting to +Z"
            );
            [0.0, 0.0, 1.0]
        }
    }
}

/// Compute magnetic energy W_m = ½ A^T K_curl A.
fn compute_magnetic_energy(
    k: &rem_core::CsrMatrix,
    a: &[f64],
    n: usize,
    _comm: &dyn Comm,
) -> f64 {
    let mut ka = vec![0.0f64; n];
    k.matvec(a, &mut ka, &rem_parallel::NoComm);
    0.5 * a.iter().zip(ka.iter()).map(|(ai, kai)| ai * kai).sum::<f64>()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use rem_config::{load_config_from_str, ConfigFormat};
    use rem_mesh::{Element, ElementKind, Node, RemMesh};
    use rem_parallel::NoComm;
    use std::collections::HashMap;

    /// Build a 6-tet unit cube mesh for HCurl magnetostatic testing.
    ///
    /// Boundary conditions:
    ///   tag 10 (face z=0): Ground
    ///   tag 11 (face z=1): SurfaceCurrent { index: 1 }
    fn unit_cube_mesh() -> RemMesh {
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
        let tets = [
            [0usize, 1, 3, 4],
            [1, 3, 4, 5],
            [3, 4, 5, 7],
            [1, 2, 3, 5],
            [2, 3, 5, 6],
            [3, 5, 6, 7],
        ];
        let volume_elements: Vec<Element> = tets
            .iter()
            .enumerate()
            .map(|(i, ns)| Element {
                id: i + 1,
                kind: ElementKind::Tet4,
                tag: 1,
                node_ids: ns.to_vec(),
                rank: 0,
            })
            .collect();

        // NOTE: boundary triangulation must match volume element faces.
        // The 6-tet Sommerville decomposition splits:
        //   z=0 quad 0-1-2-3 into tri {0,1,3} (Tet 1) + {1,2,3} (Tet 4)
        //   z=1 quad 4-5-6-7 into tri {4,5,7} (Tet 3) + {5,6,7} (Tet 6)
        let boundary_elements = vec![
            Element { id: 100, kind: ElementKind::Tri3, tag: 10, node_ids: vec![0, 1, 3], rank: 0 },
            Element { id: 101, kind: ElementKind::Tri3, tag: 10, node_ids: vec![1, 2, 3], rank: 0 },
            Element { id: 102, kind: ElementKind::Tri3, tag: 11, node_ids: vec![4, 5, 7], rank: 0 },
            Element { id: 103, kind: ElementKind::Tri3, tag: 11, node_ids: vec![5, 6, 7], rank: 0 },
        ];
        let mut boundary_tags: HashMap<u32, BoundaryTag> = HashMap::new();
        boundary_tags.insert(10, BoundaryTag::Ground);
        boundary_tags.insert(11, BoundaryTag::SurfaceCurrent { index: 1 });

        RemMesh {
            nodes,
            volume_elements,
            boundary_elements,
            domain_tags: Default::default(),
            boundary_tags,
            dim: 3,
            rank: 0,
            size: 1,
        }
    }

    #[test]
    fn hcurl_magnetostatic_3d_b_field_and_energy() {
        let mesh = unit_cube_mesh();
        let json = format!(
            r#"{{
                "Problem": {{"Type": "Magnetostatic"}},
                "Model":   {{"Mesh": "cube.msh", "L0": 1.0}},
                "Domains": {{
                    "Materials": [{{"Attributes": [1], "Permeability": 1.0}}]
                }},
                "Boundaries": {{
                    "Ground":         {{"Attributes": [10]}},
                    "SurfaceCurrent": [{{"Index": 1, "Attributes": [11], "Direction": "+X"}}]
                }},
                "Solver": {{
                    "Order": 1,
                    "Linear": {{"Tol": 1e-12, "MaxIter": 500}}
                }}
            }}"#
        );
        let config =
            load_config_from_str(&json, ConfigFormat::Json).expect("config should parse");
        let domain_map = DomainMap::from_config(&config).expect("DomainMap::from_config failed");

        let result = run_hcurl_3d(&config, &mesh, &domain_map, &NoComm)
            .expect("HCurl magnetostatic solve failed");

        // Magnetic energy must be positive for a non-trivial source
        assert!(
            result.energy > 0.0,
            "Magnetic energy should be positive, got {:.4e}",
            result.energy
        );

        // B-field should be non-zero somewhere
        let max_b: f64 = result
            .b_field
            .iter()
            .map(|b| b[0] * b[0] + b[1] * b[1] + b[2] * b[2])
            .fold(0.0_f64, f64::max)
            .sqrt();
        assert!(
            max_b > 1e-15,
            "B-field should be non-zero, max |B| = {:.2e}",
            max_b
        );

        // The HCurl path flag must be set
        assert!(result.is_hcurl, "solver should report HCurl path");
    }
}
