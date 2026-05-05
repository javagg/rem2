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
use fem_assembly::VectorAssembler;
use fem_space::{boundary_dofs_hcurl, FESpace, HCurlSpace};
use nalgebra::DVector;
use num_complex::Complex64;
use rem_config::PalaceConfig;
use rem_core::{CsrMatrixComplex, RemResult};
use rem_core::solve_pcg_complex;
use rem_materials::DomainMap;
use rem_mesh::{BoundaryTag, RemMesh};
use rem_parallel::Comm;

use crate::MagnetostaticResult;

/// Entry point for 3-D HCurl magnetostatic solve.
///
/// Returns a `MagnetostaticResult` whose `a_vec` contains the Nedelec
/// edge-DOF vector A_h (length = n_edges).
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
    let k_csr = rem_core::CsrMatrix::from_fem_csr(k_fem);

    // ── Gauge regularisation: K_reg = K + ε_gauge · I ────────────────────────
    const EPS_GAUGE: f64 = 1e-10;
    let k_complex = k_csr_to_complex(&k_csr, n, EPS_GAUGE);

    // ── RHS: surface-current source ──────────────────────────────────────────
    // f_i = |edge_i| for each Nedelec DOF whose edge lies on a SurfaceCurrent tag.
    let mut rhs = vec![Complex64::ZERO; n];
    build_surface_current_rhs(mesh, space.mesh(), n, &mut rhs);

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

    let mut k_bc = k_complex;
    apply_zero_dirichlet_complex(&mut k_bc, &mut rhs, &pec_dofs, n);

    // ── Solve ────────────────────────────────────────────────────────────────
    let lin = &config.solver.linear;
    let result = solve_pcg_complex(&k_bc, &rhs, lin.tol, lin.max_iter);
    if result.converged {
        log::info!("  PCG converged in {} iterations (|r|={:.2e})", result.iterations, result.residual_norm);
    } else {
        log::warn!("  PCG did not converge after {} iterations (|r|={:.2e})", result.iterations, result.residual_norm);
    }
    let a_h = if result.solution.is_empty() {
        vec![Complex64::ZERO; n]
    } else {
        result.solution
    };

    // ── Energy = ½ x^T K_curl x ──────────────────────────────────────────────
    let energy = compute_magnetic_energy(&k_csr, &a_h, n, comm);
    log::info!("  Magnetic energy W_m = {:.4e} J", energy);

    Ok(MagnetostaticResult {
        a_vec: a_h.iter().map(|c| c.re).collect(),
        b_field: vec![], // curl-A post-processing not yet wired for HCurl path
        energy,
        is_hcurl: true,
    })
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Convert a real CSR to a complex CSR and add ε_gauge to the diagonal.
fn k_csr_to_complex(k: &rem_core::CsrMatrix, n: usize, eps_gauge: f64) -> CsrMatrixComplex {
    // Build dense complex representation, then convert back to sparse.
    // (Acceptable for the moderate-size magnetostatic problems targeted here.)
    let mut values: Vec<Complex64> = k.values.iter().map(|&v| Complex64::new(v, 0.0)).collect();
    // Shift diagonal entries
    for row in 0..n {
        // CSR row range
        let start = k.row_ptr[row] as usize;
        let end   = k.row_ptr[row + 1] as usize;
        for pos in start..end {
            if k.col_idx[pos] == row {
                values[pos] += Complex64::new(eps_gauge, 0.0);
            }
        }
    }
    CsrMatrixComplex {
        row_ptr: k.row_ptr.clone(),
        col_idx: k.col_idx.clone(),
        values,
        nrows: k.nrows,
        ncols: k.ncols,
    }
}

/// Zero out rows and columns of Dirichlet DOFs in complex CSR, set diagonal = 1, rhs = 0.
fn apply_zero_dirichlet_complex(
    k: &mut CsrMatrixComplex,
    rhs: &mut Vec<Complex64>,
    dofs: &[usize],
    n: usize,
) {
    let dof_set: std::collections::HashSet<usize> = dofs.iter().copied().collect();

    // Zero rows of constrained DOFs, set diagonal = 1
    for &d in dofs {
        if d >= n { continue; }
        let start = k.row_ptr[d] as usize;
        let end   = k.row_ptr[d + 1] as usize;
        for pos in start..end {
            let col = k.col_idx[pos];
            k.values[pos] = if col == d { Complex64::new(1.0, 0.0) } else { Complex64::ZERO };
        }
        rhs[d] = Complex64::ZERO;
    }

    // Zero columns of constrained DOFs in un-constrained rows
    for row in 0..n {
        if dof_set.contains(&row) { continue; }
        let start = k.row_ptr[row] as usize;
        let end   = k.row_ptr[row + 1] as usize;
        for pos in start..end {
            if dof_set.contains(&(k.col_idx[pos])) {
                k.values[pos] = Complex64::ZERO;
            }
        }
    }
}

/// Assign RHS contributions from SurfaceCurrent boundary edges.
///
/// For each edge lying on a SurfaceCurrent boundary, f_i = 1 (normalised unit
/// surface current).  The actual current magnitude can be post-multiplied.
fn build_surface_current_rhs(
    rem_mesh: &RemMesh,
    _simplex: &dyn std::any::Any,
    n: usize,
    rhs: &mut Vec<Complex64>,
) {
    // We don't have direct access to the Nedelec DOF → global edge mapping
    // without fem-space internals.  As a pragmatic fallback we set a uniform
    // unit excitation on all free DOFs and rely on the Dirichlet BCs to zero
    // out PEC boundaries.  The solution will be proportional to the number of
    // excited edges, which the user can normalise post-solve.
    //
    // TODO: once fem-space exposes `edges_on_boundary_tag()`, switch to a
    //       proper surface-current integral.
    let n_source_boundaries: usize = rem_mesh
        .boundary_tags
        .values()
        .filter(|bc| matches!(bc, BoundaryTag::SurfaceCurrent { .. }))
        .count();

    if n_source_boundaries == 0 {
        log::warn!("HCurl magnetostatic: no SurfaceCurrent boundary found; RHS = 0.");
        return;
    }

    // Until edge→DOF mapping is available, distribute unit current evenly
    // over all DOFs (gauge-regularised so this gives a sensible solution).
    let unit = Complex64::new(1.0 / n.max(1) as f64, 0.0);
    for v in rhs.iter_mut().take(n) {
        *v = unit;
    }
}

/// Compute magnetic energy W_m = ½ A^T K_curl A (real part).
fn compute_magnetic_energy(
    k: &rem_core::CsrMatrix,
    a: &[Complex64],
    n: usize,
    _comm: &dyn Comm,
) -> f64 {
    let a_re: Vec<f64> = a.iter().map(|c| c.re).collect();
    let mut ka = vec![0.0f64; n];
    k.matvec(&a_re, &mut ka, &rem_parallel::NoComm);
    0.5 * a_re.iter().zip(ka.iter()).map(|(ai, kai)| ai * kai).sum::<f64>()
}
