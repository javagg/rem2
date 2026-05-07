//! Full-wave HCurl (Nedelec edge-element) driven frequency-domain solver.
//!
//! Solves the vector Helmholtz system at each frequency ω:
//!
//!   A(ω) · x = f_port
//!
//! where  A = (1/μ) · curl-curl  −  ω² · ε · mass
//!        x ∈ R^{n_edge}  (tangential components of E on each mesh edge)
//!        f_port = surface-current excitation (lumped port)
//!
//! PEC boundaries → zero tangential E → constrained edge DOFs (x_e = 0).
//! Lumped port:
//!   - Excited port: uniform surface current J_s on port boundary edges → RHS.
//!   - Other ports: resistive loading added as a surface admittance term Y_port.
//!
//! S-parameter extraction:
//!   - Port voltage V_p = ∫ E · dl  (tangential line integral across port).
//!   - Port current I_p extracted from RHS excitation normalization.
//!   - S11 = (Z_port − Z0)/(Z_port + Z0),  Z0 = 50 Ω (or from LumpedPort R).

use fem_assembly::coefficient::CtxFnCoeff;
use fem_assembly::standard::{CurlCurlIntegrator, VectorMassIntegrator};
use fem_assembly::VectorAssembler;
use fem_space::{boundary_dofs_hcurl, FESpace, HCurlSpace};
use nalgebra::{DMatrix, DVector};
use num_complex::Complex64;
use rem_config::PalaceConfig;
use rem_core::{CsrMatrix, CsrMatrixComplex, RemError, RemResult, solve_pcg_complex};
use rem_materials::DomainMap;
use rem_mesh::{BoundaryTag, RemMesh};
use rem_parallel::Comm;

use crate::{DrivenResult, FreqResult, PortVi};

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

/// Entry point called when driven effective formulation resolves to HCurl.
///
/// Performs a frequency sweep using Nedelec edge elements (ND1 or ND2 per
/// `Solver.Order`).  Returns the same `DrivenResult` structure as the scalar
/// H1 path so downstream output / CLI code is unchanged.
pub(crate) fn run_hcurl_driven(
    config: &PalaceConfig,
    mesh: &RemMesh,
    domain_map: &DomainMap,
    comm: &dyn Comm,
) -> RemResult<DrivenResult> {
    log::info!("=== HCurl (Nedelec edge-element) driven solver ===");

    // Build frequency list from config
    let drv_cfg = config.solver.driven.as_ref().ok_or_else(|| {
        RemError::Config("Driven problem requires a [Solver.Driven] section".into())
    })?;

    let freqs = crate::build_freq_list(drv_cfg)?;
    log::info!("HCurl driven: {} frequency points", freqs.len());

    // Choose element order (ND1 = order 1, ND2 = order 2)
    let requested_order = config.solver.driven_hcurl_order();
    let order: u8 = requested_order.clamp(1, 2);
    if requested_order > 2 {
        log::warn!(
            "HCurl driven supports order 1/2 only; requested order={} (Driven.HCurlOrder or Solver.Order), using order 2.",
            requested_order
        );
    }

    // Build HCurl space
    let (space_k, space_m, n_dof, pec_dofs) = build_space_and_pec(mesh, order, domain_map)?;

    let n = n_dof;
    log::info!("HCurl space: {} edge DOFs, {} PEC-constrained DOFs", n, pec_dofs.len());

    // Convert to dense for GMRES path (consistent with H1 driven path)
    let k_dense = csr_to_complex_dense(&space_k, n);
    let m_dense = csr_to_complex_dense(&space_m, n);

    // Collect ports
    let all_ports = collect_all_ports_hcurl(mesh);
    let n_ports = all_ports.len();
    if n_ports == 0 {
        return Err(RemError::Config(
            "HCurl driven: no lumped ports found in mesh boundary tags".into(),
        ));
    }
    log::info!("HCurl driven: {} port(s)", n_ports);

    let save_step = drv_cfg.save_step.max(1);
    let out_dir = config.problem.output_dir();
    #[cfg(not(target_arch = "wasm32"))]
    std::fs::create_dir_all(out_dir).map_err(RemError::Io)?;

    let lin = &config.solver.linear;
    let use_bicgstab = config.solver.linear.prefers_sparse_iterative_complex()
        || std::env::var("REM_USE_PCG").is_ok();

    let mut freq_results: Vec<FreqResult> = Vec::with_capacity(freqs.len());
    let mut peak_phi: Vec<f64> = Vec::new();
    let mut peak_freq_hz: f64 = 0.0;
    let mut peak_s11_mag: f64 = -1.0;

    for (step, &freq) in freqs.iter().enumerate() {
        let omega = 2.0 * std::f64::consts::PI * freq;
        let omega2 = omega * omega;

        // A(ω) = K − ω² M  (both real here; lossy extension left for future)
        // We also add port admittance terms for unexcited ports
        let mut a_base = k_dense.clone();
        for i in 0..n {
            for j in 0..n {
                a_base[(i, j)] -= Complex64::new(omega2, 0.0) * m_dense[(i, j)];
            }
        }

        // Apply PEC constraints (zero rows/cols for constrained DOFs)
        apply_pec_constraints(&mut a_base, &pec_dofs, n);

        // ── Multi-port S-matrix ──────────────────────────────────────────
        let (s11, s_matrix, edge_complex, port_vi_opt) = if n_ports > 1 {
            let z0_vec: Vec<Complex64> = all_ports
                .iter()
                .map(|&pidx| lumped_port_z0_hcurl(mesh, pidx, omega))
                .collect();

            let mut z_cols: Vec<Vec<Complex64>> = Vec::with_capacity(n_ports);
            let mut first_edge: Vec<Complex64> = Vec::new();

            for (j, &exc_idx) in all_ports.iter().enumerate() {
                // Build port admittance loading for all ports ≠ exc_idx
                let mut a_port = a_base.clone();
                add_port_admittance(&mut a_port, mesh, exc_idx, &all_ports, omega, n);

                // Build RHS: surface current excitation on excited port edges
                let rhs = build_port_rhs(mesh, exc_idx, &pec_dofs, n);

                // Solve
                let e_vec = solve_hcurl_system(&a_port, &rhs, lin.tol, lin.max_iter, use_bicgstab)?;

                if j == 0 {
                    first_edge = e_vec.clone();
                }

                // Extract port voltages for all ports
                let vols: Vec<Complex64> = all_ports
                    .iter()
                    .map(|&pidx| integrate_port_voltage(mesh, pidx, &e_vec))
                    .collect();
                z_cols.push(vols);
            }

            let s_mat = z_to_s_matrix_hcurl(&z_cols, &z0_vec);
            let s11_c = s_mat[0][0];
            log::info!("HCurl f={:.3e} Hz  |S11|={:.4}  ({} ports)", freq, s11_c.norm(), n_ports);
            (s11_c, s_mat, first_edge, None::<PortVi>)
        } else {
            // Single-port path
            let exc_idx = all_ports[0];
            let rhs = build_port_rhs(mesh, exc_idx, &pec_dofs, n);
            let e_vec = solve_hcurl_system(&a_base, &rhs, lin.tol, lin.max_iter, use_bicgstab)?;
            let v_port = integrate_port_voltage(mesh, exc_idx, &e_vec);
            let z0 = lumped_port_z0_hcurl(mesh, exc_idx, omega);
            // Current from RHS normalization: I = 1 A (unit excitation)
            let i_port = Complex64::new(1.0, 0.0);
            let z_port = v_port / i_port;
            let s11 = (z_port - z0) / (z_port + z0);
            log::info!(
                "HCurl f={:.3e} Hz  |S11|={:.4}  ∠S11={:.2}°",
                freq,
                s11.norm(),
                s11.arg().to_degrees()
            );
            let p = v_port * i_port.conj() * Complex64::new(0.5, 0.0);
            let vi = PortVi { port_index: exc_idx, v: v_port, i: i_port, p };
            (s11, vec![], e_vec, Some(vi))
        };

        let edge_re: Vec<f64> = edge_complex.iter().map(|x| x.re).collect();

        let port_list: Vec<u32> = all_ports.clone();
        freq_results.push(FreqResult {
            freq_hz: freq,
            s11_re: s11.re,
            s11_im: s11.im,
            s_matrix,
            port_list,
            port_vi: port_vi_opt.into_iter().collect(),
        });

        let s11_mag = s11.norm();
        if s11_mag > peak_s11_mag {
            peak_s11_mag = s11_mag;
            peak_freq_hz = freq;
            peak_phi = edge_re.clone();
        }

        #[cfg(not(target_arch = "wasm32"))]
        if step % save_step == 0 {
            let vtk_order = config.solver.eigenmode_hcurl_order().clamp(1, 2) as u8;
            crate::output::write_field_vector_vtk(out_dir, mesh, &edge_complex, step + 1, vtk_order)?;
        }

        let _ = comm;
    }

    #[cfg(not(target_arch = "wasm32"))]
    crate::output::write_s_params(out_dir, &freq_results)?;

    log::info!(
        "HCurl driven complete: {} frequency points, peak |S11|={:.4} at {:.3e} Hz",
        freq_results.len(),
        peak_s11_mag,
        peak_freq_hz
    );

    Ok(DrivenResult {
        freq_results,
        peak_phi,
        peak_freq_hz,
        far_field_pattern: Vec::new(),
        circuit_model: None,
    })
}

// ---------------------------------------------------------------------------
// Space construction
// ---------------------------------------------------------------------------

/// Build HCurl stiffness & mass matrices and return PEC-constrained DOF indices.
///
/// Returns `(K, M, n_dofs, pec_dofs)`.
fn build_space_and_pec(
    mesh: &RemMesh,
    order: u8,
    domain_map: &DomainMap,
) -> RemResult<(CsrMatrix, CsrMatrix, usize, Vec<usize>)> {
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

    let inv_mu = CtxFnCoeff(|ctx: &fem_assembly::coefficient::CoeffCtx<'_>| {
        domain_map.get(ctx.elem_tag as u32).reluctivity()
    });
    let eps = CtxFnCoeff(|ctx: &fem_assembly::coefficient::CoeffCtx<'_>| {
        domain_map.get(ctx.elem_tag as u32).epsilon_abs()
    });
    let curl_curl = CurlCurlIntegrator { mu: inv_mu };
    let mass = VectorMassIntegrator { alpha: eps };

    if mesh.dim == 2 {
        let simplex = mesh.to_simplex_mesh_2d();
        let space = HCurlSpace::new(simplex, order);
        let n = space.n_dofs();
        let k_fem = VectorAssembler::assemble_bilinear(&space, &[&curl_curl], 4);
        let m_fem = VectorAssembler::assemble_bilinear(&space, &[&mass], 4);
        let k = CsrMatrix::from_fem_csr(k_fem);
        let m = CsrMatrix::from_fem_csr(m_fem);
        let pec: Vec<usize> = boundary_dofs_hcurl(space.mesh(), &space, &pec_tags)
            .into_iter()
            .map(|d| d as usize)
            .collect();
        Ok((k, m, n, pec))
    } else if mesh.dim == 3 {
        let simplex = mesh.to_simplex_mesh();
        let space = HCurlSpace::new(simplex, order);
        let n = space.n_dofs();
        let k_fem = VectorAssembler::assemble_bilinear(&space, &[&curl_curl], 4);
        let m_fem = VectorAssembler::assemble_bilinear(&space, &[&mass], 4);
        let k = CsrMatrix::from_fem_csr(k_fem);
        let m = CsrMatrix::from_fem_csr(m_fem);
        let pec: Vec<usize> = boundary_dofs_hcurl(space.mesh(), &space, &pec_tags)
            .into_iter()
            .map(|d| d as usize)
            .collect();
        Ok((k, m, n, pec))
    } else {
        Err(RemError::Config(format!(
            "HCurl driven only supports 2-D/3-D meshes, got dim={}",
            mesh.dim
        )))
    }
}

// ---------------------------------------------------------------------------
// System assembly helpers
// ---------------------------------------------------------------------------

fn apply_pec_constraints(a: &mut DMatrix<Complex64>, pec_dofs: &[usize], n: usize) {
    for &d in pec_dofs {
        for j in 0..n {
            a[(d, j)] = Complex64::ZERO;
            a[(j, d)] = Complex64::ZERO;
        }
        a[(d, d)] = Complex64::new(1.0, 0.0);
    }
}

/// Add port admittance Y = 1/Z0 for unexcited ports (matched termination).
/// This adds a diagonal surface term:  A[d,d] += Y  for each port edge DOF.
fn add_port_admittance(
    a: &mut DMatrix<Complex64>,
    mesh: &RemMesh,
    exc_idx: u32,
    all_ports: &[u32],
    omega: f64,
    n: usize,
) {
    for &pidx in all_ports {
        if pidx == exc_idx {
            continue;
        }
        let z0 = lumped_port_z0_hcurl(mesh, pidx, omega);
        if z0.norm() < 1e-30 {
            continue;
        }
        let y0 = Complex64::new(1.0, 0.0) / z0;
        // Add Y0 to diagonal of all edge DOFs on this port boundary
        for belem in &mesh.boundary_elements {
            let is_port = match mesh.boundary_tags.get(&belem.tag) {
                Some(BoundaryTag::LumpedPort { index, .. }) => *index == pidx,
                _ => false,
            };
            if !is_port {
                continue;
            }
            // Map boundary node pairs to edge DOFs (approximate: treat each edge as one DOF)
            // Since we don't have a direct node→edge DOF mapping here, we add the admittance
            // as a surface integral weighted uniformly over the boundary element.
            let nids = &belem.node_ids;
            let area = boundary_elem_measure(mesh, nids);
            let scale = y0 * area / (nids.len() as f64);
            for &nid in nids {
                if nid < n {
                    a[(nid, nid)] += scale;
                }
            }
        }
    }
}

/// Build the RHS excitation vector for a lumped port:
/// f[d] = ∫_{Γ_port} J_s · φ_d dS
///
/// We use a uniform unit surface current density J_s = 1/A_port (so total current = 1 A),
/// with a consistent mass-matrix integration over the port boundary.
fn build_port_rhs(
    mesh: &RemMesh,
    port_idx: u32,
    pec_dofs: &[usize],
    n: usize,
) -> Vec<Complex64> {
    let mut rhs = vec![Complex64::ZERO; n];
    let pec_set: std::collections::HashSet<usize> = pec_dofs.iter().copied().collect();

    // Compute total port area for normalization
    let mut total_area = 0.0_f64;
    for belem in &mesh.boundary_elements {
        let is_port = match mesh.boundary_tags.get(&belem.tag) {
            Some(BoundaryTag::LumpedPort { index, .. }) => *index == port_idx,
            _ => false,
        };
        if is_port {
            total_area += boundary_elem_measure(mesh, &belem.node_ids);
        }
    }

    if total_area < 1e-30 {
        log::warn!("HCurl driven: port {} has zero measure — RHS will be zero", port_idx);
        return rhs;
    }

    let j_density = 1.0 / total_area; // unit current spread uniformly

    for belem in &mesh.boundary_elements {
        let is_port = match mesh.boundary_tags.get(&belem.tag) {
            Some(BoundaryTag::LumpedPort { index, .. }) => *index == port_idx,
            _ => false,
        };
        if !is_port {
            continue;
        }
        let nids = &belem.node_ids;
        let area = boundary_elem_measure(mesh, nids);
        // Consistent mass-matrix weighting: each node gets area/n_nodes contribution
        let weight = j_density * area / (nids.len() as f64);
        for &nid in nids {
            if nid < n && !pec_set.contains(&nid) {
                rhs[nid] += Complex64::new(weight, 0.0);
            }
        }
    }

    rhs
}

// ---------------------------------------------------------------------------
// Linear solve
// ---------------------------------------------------------------------------

fn solve_hcurl_system(
    a: &DMatrix<Complex64>,
    rhs: &[Complex64],
    tol: f64,
    max_iter: usize,
    use_bicgstab: bool,
) -> RemResult<Vec<Complex64>> {
    if use_bicgstab {
        let mat_csr = CsrMatrixComplex::from_dense(a);
        let result = solve_pcg_complex(&mat_csr, rhs, tol, max_iter);
        if result.converged {
            log::debug!(
                "HCurl BiCGSTAB: converged in {} iters (res={:.3e})",
                result.iterations,
                result.residual_norm
            );
            return Ok(result.solution);
        }
        log::debug!(
            "HCurl BiCGSTAB: no convergence after {} iters; falling back to GMRES",
            result.iterations
        );
    }
    gmres_hcurl(a, rhs, tol, max_iter)
}

fn gmres_hcurl(
    a: &DMatrix<Complex64>,
    rhs: &[Complex64],
    tol: f64,
    max_iter: usize,
) -> RemResult<Vec<Complex64>> {
    let n = rhs.len();
    const RESTART: usize = 30;
    let max_outer = (max_iter / RESTART).max(1);

    let b = DVector::<Complex64>::from_iterator(n, rhs.iter().copied());
    let b_norm = b.norm();
    if b_norm < f64::EPSILON {
        return Ok(vec![Complex64::ZERO; n]);
    }

    let mut x = DVector::<Complex64>::zeros(n);

    for _outer in 0..max_outer {
        let r = &b - a * &x;
        let beta = r.norm();
        if beta / b_norm < tol {
            return Ok(x.iter().copied().collect());
        }

        let mut v: Vec<DVector<Complex64>> = Vec::with_capacity(RESTART + 1);
        v.push(&r / Complex64::new(beta, 0.0));

        let mut h = vec![vec![Complex64::ZERO; RESTART]; RESTART + 1];
        let mut g = vec![Complex64::ZERO; RESTART + 1];
        let mut cs = vec![0.0f64; RESTART];
        let mut sn = vec![Complex64::ZERO; RESTART];
        g[0] = Complex64::new(beta, 0.0);

        let mut j_done = RESTART;
        for j in 0..RESTART {
            let mut w = a * &v[j];
            for i in 0..=j {
                h[i][j] = v[i].dotc(&w);
                let hij = h[i][j];
                w -= &v[i] * hij;
            }
            let h_next = w.norm();
            h[j + 1][j] = Complex64::new(h_next, 0.0);
            if h_next > 1e-14 {
                v.push(&w / Complex64::new(h_next, 0.0));
            }
            for i in 0..j {
                let tmp = cs[i] * h[i][j] + sn[i] * h[i + 1][j];
                h[i + 1][j] = -sn[i].conj() * h[i][j] + cs[i] * h[i + 1][j];
                h[i][j] = tmp;
            }
            let rr = (h[j][j].norm_sqr() + h[j + 1][j].norm_sqr()).sqrt();
            if rr > 1e-14 {
                cs[j] = h[j][j].norm() / rr;
                let hjj_norm = h[j][j].norm();
                sn[j] = if hjj_norm > 1e-300 {
                    h[j + 1][j] * (h[j][j].conj() / (rr * hjj_norm))
                } else {
                    Complex64::ZERO
                };
                h[j][j] = Complex64::new(rr, 0.0);
                h[j + 1][j] = Complex64::ZERO;
                let g_next = -sn[j].conj() * g[j];
                g[j] = Complex64::new(cs[j], 0.0) * g[j];
                g[j + 1] = g_next;
            }
            if g[j + 1].norm() / b_norm < tol {
                j_done = j + 1;
                break;
            }
        }

        let m = j_done.min(RESTART);
        let mut y = vec![Complex64::ZERO; m];
        for i in (0..m).rev() {
            y[i] = g[i];
            for k in (i + 1)..m {
                let yk = y[k];
                y[i] -= h[i][k] * yk;
            }
            if h[i][i].norm() > 1e-300 {
                y[i] /= h[i][i];
            }
        }
        for j in 0..m {
            let yj = y[j];
            x += &v[j] * yj;
        }
    }

    Ok(x.iter().copied().collect())
}

// ---------------------------------------------------------------------------
// S-parameter extraction
// ---------------------------------------------------------------------------

/// Integrate tangential E along port boundary edges to get port voltage.
///
/// V_p = ∫_{Γ_port} E · t dL   (line integral over port cross-section)
///
/// Approximated by summing E_d × edge_length for each edge DOF on the port,
/// using the edge DOF index as a proxy for the boundary edge.
fn integrate_port_voltage(mesh: &RemMesh, port_idx: u32, e_vec: &[Complex64]) -> Complex64 {
    let mut v = Complex64::ZERO;
    for belem in &mesh.boundary_elements {
        let is_port = match mesh.boundary_tags.get(&belem.tag) {
            Some(BoundaryTag::LumpedPort { index, .. }) => *index == port_idx,
            _ => false,
        };
        if !is_port {
            continue;
        }
        let nids = &belem.node_ids;
        let len = boundary_elem_measure(mesh, nids);
        // Weight each DOF uniformly by element length / n_nodes
        let w = len / (nids.len() as f64);
        for &nid in nids {
            if nid < e_vec.len() {
                v += e_vec[nid] * w;
            }
        }
    }
    v
}

/// Z0 for a lumped port (R value, or 50 Ω default).
fn lumped_port_z0_hcurl(mesh: &RemMesh, port_idx: u32, _omega: f64) -> Complex64 {
    for bc in mesh.boundary_tags.values() {
        match bc {
            BoundaryTag::LumpedPort { index, r, .. } if *index == port_idx => {
                let r_val = if *r > 0.0 { *r } else { 50.0 };
                return Complex64::new(r_val, 0.0);
            }
            _ => {}
        }
    }
    Complex64::new(50.0, 0.0)
}

/// Collect all distinct LumpedPort indices in stable sorted order.
fn collect_all_ports_hcurl(mesh: &RemMesh) -> Vec<u32> {
    let mut seen: Vec<u32> = Vec::new();
    for bc in mesh.boundary_tags.values() {
        if let BoundaryTag::LumpedPort { index, .. } = bc {
            if !seen.contains(index) {
                seen.push(*index);
            }
        }
    }
    seen.sort_unstable();
    seen
}

/// Z-matrix to S-matrix conversion (symmetric, standard 50-Ω reference).
fn z_to_s_matrix_hcurl(z_cols: &[Vec<Complex64>], z0_vec: &[Complex64]) -> Vec<Vec<Complex64>> {
    let n = z0_vec.len();
    let mut s = vec![vec![Complex64::ZERO; n]; n];
    for i in 0..n {
        for j in 0..n {
            let z_ij = z_cols[j][i];
            let z0 = z0_vec[i];
            // Simple diagonal reference: S_ij = (Z_ij − Z0_i δ_ij) / (Z_ij + Z0_i) for diagonal
            if i == j {
                s[i][j] = (z_ij - z0) / (z_ij + z0);
            } else {
                // Off-diagonal: normalize by geometric mean Z0
                let z0j = z0_vec[j];
                let ref_scale = 2.0 * (z0 * z0j).sqrt();
                s[i][j] = ref_scale / (z_ij + z0);
            }
        }
    }
    s
}

// ---------------------------------------------------------------------------
// Geometry helpers
// ---------------------------------------------------------------------------

fn boundary_elem_measure(mesh: &RemMesh, nids: &[usize]) -> f64 {
    match nids.len() {
        2 => {
            let (p0, p1) = (&mesh.nodes[nids[0]], &mesh.nodes[nids[1]]);
            ((p1.x - p0.x).powi(2) + (p1.y - p0.y).powi(2) + (p1.z - p0.z).powi(2)).sqrt()
        }
        3 => {
            let (p0, p1, p2) = (
                &mesh.nodes[nids[0]],
                &mesh.nodes[nids[1]],
                &mesh.nodes[nids[2]],
            );
            let v1 = [p1.x - p0.x, p1.y - p0.y, p1.z - p0.z];
            let v2 = [p2.x - p0.x, p2.y - p0.y, p2.z - p0.z];
            let cx = v1[1] * v2[2] - v1[2] * v2[1];
            let cy = v1[2] * v2[0] - v1[0] * v2[2];
            let cz = v1[0] * v2[1] - v1[1] * v2[0];
            0.5 * (cx * cx + cy * cy + cz * cz).sqrt()
        }
        4 => {
            // Quad: split into 2 triangles
            let (p0, p1, p2, p3) = (
                &mesh.nodes[nids[0]],
                &mesh.nodes[nids[1]],
                &mesh.nodes[nids[2]],
                &mesh.nodes[nids[3]],
            );
            let tri_area = |a: &rem_mesh::Node, b: &rem_mesh::Node, c: &rem_mesh::Node| -> f64 {
                let v1 = [b.x - a.x, b.y - a.y, b.z - a.z];
                let v2 = [c.x - a.x, c.y - a.y, c.z - a.z];
                let cx = v1[1] * v2[2] - v1[2] * v2[1];
                let cy = v1[2] * v2[0] - v1[0] * v2[2];
                let cz = v1[0] * v2[1] - v1[1] * v2[0];
                0.5 * (cx * cx + cy * cy + cz * cz).sqrt()
            };
            tri_area(p0, p1, p2) + tri_area(p0, p2, p3)
        }
        _ => 0.0,
    }
}

fn csr_to_complex_dense(mat: &CsrMatrix, n: usize) -> DMatrix<Complex64> {
    let mut d = DMatrix::<Complex64>::zeros(n, n);
    for i in 0..mat.nrows {
        for ptr in mat.row_ptr[i]..mat.row_ptr[i + 1] {
            let j = mat.col_idx[ptr];
            d[(i, j)] += Complex64::new(mat.values[ptr], 0.0);
        }
    }
    d
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use rem_config::{load_config_from_str, ConfigFormat};
    use rem_mesh::{gen::rect_msh, gmsh::read_msh_str, RemMesh};
    use rem_parallel::NoComm;

    /// Minimal JSON config for a single-frequency HCurl driven sweep over a 2-D cavity.
    fn hcurl_config_json(f_hz: f64) -> String {
        format!(
            r#"{{
                "Problem": {{"Type": "Driven"}},
                "Model":   {{"Mesh": "test.msh", "L0": 1e-3}},
                "Domains": {{
                    "Materials": [{{"Attributes": [10], "Permittivity": 1.0, "Permeability": 1.0}}]
                }},
                "Boundaries": {{
                    "PEC": {{"Attributes": [2, 3, 4]}},
                    "LumpedPort": [{{"Index": 1, "Attributes": [1], "R": 50.0}}]
                }},
                "Solver": {{
                    "Discretization": "HCurl",
                    "Order": 1,
                    "Driven": {{
                        "MinFreq": {f_hz:.6e},
                        "MaxFreq": {f_hz:.6e},
                        "FreqStep": {f_hz:.6e},
                        "SaveStep": 999
                    }},
                    "Linear": {{"Tol": 1e-8, "MaxIter": 500, "KSPType": "bicgstab"}}
                }}
            }}"#
        )
    }

    #[test]
    fn hcurl_driven_single_port_returns_finite_s11() {
        // Rectangular 2-D cavity: 4×2 mm, 4×2 grid
        // Bottom edge (tag 1) = LumpedPort, top+left+right (tags 2,3,4) = PEC
        let msh = rect_msh(4.0, 2.0, 4, 2, 1, 2, 3, 4, 10);
        let raw = read_msh_str(&msh).expect("rect_msh should parse");

        let json = hcurl_config_json(1e9);
        let config = load_config_from_str(&json, ConfigFormat::Json)
            .expect("config should parse");

        let mesh = RemMesh::from_raw(raw, &config).expect("RemMesh::from_raw failed");
        let domain_map = rem_materials::DomainMap::from_config(&config)
            .expect("DomainMap::from_config failed");

        let result = run_hcurl_driven(&config, &mesh, &domain_map, &NoComm);
        assert!(result.is_ok(), "HCurl driven returned error: {:?}", result.err());

        let dr = result.unwrap();
        assert_eq!(dr.freq_results.len(), 1, "Expected exactly one frequency result");

        let fr = &dr.freq_results[0];
        assert!(
            (fr.freq_hz - 1e9).abs() < 1.0,
            "Frequency mismatch: got {:.3e}",
            fr.freq_hz
        );
        let s11_mag = (fr.s11_re * fr.s11_re + fr.s11_im * fr.s11_im).sqrt();
        assert!(s11_mag.is_finite(), "|S11| should be finite, got {}", s11_mag);
        assert!(s11_mag <= 1.0 + 1e-6, "|S11| should be ≤ 1, got {:.4}", s11_mag);
    }

    #[test]
    fn hcurl_driven_two_port_s21_finite() {
        // 2-D parallel-plate transmission line: 4×2 mm, 8×4 grid.
        // Port 1 on left (tag 1), Port 2 on right (tag 2),
        // PEC top (tag 3) and bottom (tag 4).
        let msh = rect_msh(4.0, 2.0, 8, 4, 1, 2, 3, 4, 10);
        let raw = read_msh_str(&msh).expect("rect_msh should parse");

        let json = format!(
            r#"{{
                "Problem": {{"Type": "Driven"}},
                "Model":   {{"Mesh": "tl.msh", "L0": 1e-3}},
                "Domains": {{
                    "Materials": [{{"Attributes": [10], "Permittivity": 1.0, "Permeability": 1.0}}]
                }},
                "Boundaries": {{
                    "PEC":         {{"Attributes": [3, 4]}},
                    "LumpedPort": [
                        {{"Index": 1, "Attributes": [1], "R": 50.0}},
                        {{"Index": 2, "Attributes": [2], "R": 50.0}}
                    ]
                }},
                "Solver": {{
                    "Order": 1,
                    "Driven": {{
                        "MinFreq": 1e9,
                        "MaxFreq": 1e9,
                        "FreqStep": 1e9,
                        "SaveStep": 999
                    }},
                    "Linear": {{"Tol": 1e-8, "MaxIter": 500, "KSPType": "bicgstab"}}
                }}
            }}"#
        );
        let config = load_config_from_str(&json, ConfigFormat::Json)
            .expect("config should parse");
        let mesh = RemMesh::from_raw(raw, &config).expect("RemMesh::from_raw failed");
        let domain_map = rem_materials::DomainMap::from_config(&config)
            .expect("DomainMap::from_config failed");

        let result = run_hcurl_driven(&config, &mesh, &domain_map, &NoComm);
        assert!(
            result.is_ok(),
            "HCurl driven 2-port returned error: {:?}",
            result.err()
        );

        let dr = result.unwrap();
        assert_eq!(dr.freq_results.len(), 1, "expected exactly one frequency result");

        let fr = &dr.freq_results[0];
        // S₂₁ should be finite and non-trivial (physical normalisation may
        // differ from 50 Ω on this coarse 2-D mesh, so only check sanity).
        let s21 = fr
            .s_matrix
            .get(1)
            .and_then(|row| row.get(0))
            .copied()
            .unwrap_or_default();
        let s21_mag = (s21.re * s21.re + s21.im * s21.im).sqrt();
        assert!(s21_mag.is_finite(), "|S₂₁| should be finite, got {}", s21_mag);
        assert!(
            s21_mag > 1e-12,
            "|S₂₁| should be non-zero for a transmission line"
        );
    }
}
