use std::collections::HashMap;

use fem_assembly::coefficient::CtxFnCoeff;
use fem_assembly::standard::{CurlCurlIntegrator, VectorMassIntegrator};
use fem_assembly::VectorAssembler;
use fem_space::{boundary_dofs_hcurl, FESpace, HCurlSpace};
use rem_config::PalaceConfig;
use rem_core::{CsrMatrix, RemError, RemResult};
use rem_materials::DomainMap;
use rem_mesh::{BoundaryTag, RemMesh};
use rem_parallel::Comm;

use crate::{lanczos, shifted_matrix, tridiag_eigen, C0, EigenResult};

pub(crate) fn solve_hcurl(
    config: &PalaceConfig,
    mesh: &RemMesh,
    domain_map: &DomainMap,
    comm: &dyn Comm,
) -> RemResult<EigenResult> {
    if !config.boundaries.periodic.is_empty() {
        return Err(RemError::Config(
            "HCurl eigenmode path does not yet support periodic boundary constraints".into(),
        ));
    }

    let eig_cfg = config.solver.eigenmode.as_ref().ok_or_else(|| {
        RemError::Config("missing Eigenmode solver config".into())
    })?;

    let requested_order = config.solver.eigenmode_hcurl_order();
    let order = requested_order.clamp(1, 2);
    if requested_order > 2 {
        log::warn!(
            "HCurl currently supports order 1/2 only; requested order={} (Eigenmode.HCurlOrder or Solver.Order), using order 2.",
            requested_order
        );
    }

    let target_hz = eig_cfg.target;
    let sigma = if target_hz > 0.0 {
        let omega = 2.0 * std::f64::consts::PI * target_hz;
        (omega / C0) * (omega / C0)
    } else {
        0.0
    };

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

    let (k_mat, m_mat, constrained_dofs, n) = if mesh.dim == 2 {
        let simplex = mesh.to_simplex_mesh_2d();
        let space = HCurlSpace::new(simplex, order);
        let (k, m) = assemble_hcurl_system(&space, domain_map);
        let constrained = boundary_dofs_hcurl(space.mesh(), &space, &pec_tags)
            .into_iter()
            .map(|d| d as usize)
            .collect::<Vec<_>>();
        (k, m, constrained, space.n_dofs())
    } else if mesh.dim == 3 {
        let simplex = mesh.to_simplex_mesh();
        let space = HCurlSpace::new(simplex, order);
        let (k, m) = assemble_hcurl_system(&space, domain_map);
        let constrained = boundary_dofs_hcurl(space.mesh(), &space, &pec_tags)
            .into_iter()
            .map(|d| d as usize)
            .collect::<Vec<_>>();
        (k, m, constrained, space.n_dofs())
    } else {
        return Err(RemError::Config(format!(
            "HCurl eigenmode only supports 2-D/3-D meshes, got dim={}",
            mesh.dim
        )));
    };

    let mut dofs = HashMap::new();
    for d in constrained_dofs {
        dofs.insert(d, 0.0f64);
    }

    let mut k_bc = k_mat.clone();
    let a_mat = shifted_matrix(&k_mat, &m_mat, sigma, n);
    let mut a_bc = a_mat;
    let mut rhs_dummy = vec![0.0f64; n];
    rem_electrostatic::bc::apply_dirichlet(&mut a_bc, &mut rhs_dummy, &dofs);
    rem_electrostatic::bc::apply_dirichlet(&mut k_bc, &mut rhs_dummy, &dofs);

    let n_modes = eig_cfg.n;
    let m_steps = (3 * n_modes + 10).min(n);
    let lin = &config.solver.linear;

    let (t_alpha, t_beta, v_basis) =
        lanczos(&a_bc, &m_mat, &dofs, n, m_steps, lin.tol, lin.max_iter, comm);

    let (ritz_vals, ritz_vecs_small) = tridiag_eigen(&t_alpha, &t_beta);

    let mut eigenpairs: Vec<(f64, Vec<f64>)> = Vec::new();
    for (k, &mu) in ritz_vals.iter().enumerate().take(n_modes) {
        if mu.abs() < 1e-300 {
            continue;
        }
        let lambda = sigma + 1.0 / mu;
        if lambda <= 0.0 {
            continue;
        }
        let freq_hz = C0 * lambda.sqrt() / (2.0 * std::f64::consts::PI);

        let mut x = vec![0.0f64; n];
        if k < ritz_vecs_small.ncols() && !v_basis.is_empty() {
            let basis_m = v_basis.len();
            let coeff_m = ritz_vecs_small.nrows().min(basis_m);
            for j in 0..coeff_m {
                let y_jk = ritz_vecs_small[(j, k)];
                if y_jk.abs() < 1e-300 {
                    continue;
                }
                let vj = &v_basis[j];
                for i in 0..n {
                    x[i] += y_jk * vj[i];
                }
            }
            let mut mx = vec![0.0f64; n];
            m_mat.matvec(&x, &mut mx, comm);
            let norm_sq: f64 = x.iter().zip(mx.iter()).map(|(a, b)| a * b).sum();
            if norm_sq > 1e-300 {
                let s = 1.0 / norm_sq.sqrt();
                for xi in &mut x {
                    *xi *= s;
                }
            }
        }

        eigenpairs.push((freq_hz, x));
    }

    eigenpairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    Ok(EigenResult {
        frequencies_hz: eigenpairs.iter().map(|(f, _)| *f).collect(),
        eigenvectors: eigenpairs.into_iter().map(|(_, v)| v).collect(),
        q_factors: None,
        is_hcurl: true,
    })
}

fn assemble_hcurl_system<S: FESpace>(space: &S, domain_map: &DomainMap) -> (CsrMatrix, CsrMatrix) {
    let inv_mu = CtxFnCoeff(|ctx: &fem_assembly::coefficient::CoeffCtx<'_>| {
        domain_map.get(ctx.elem_tag as u32).reluctivity()
    });
    let eps = CtxFnCoeff(|ctx: &fem_assembly::coefficient::CoeffCtx<'_>| {
        domain_map.get(ctx.elem_tag as u32).epsilon_abs()
    });

    let curl_curl = CurlCurlIntegrator { mu: inv_mu };
    let mass = VectorMassIntegrator { alpha: eps };

    let k_fem = VectorAssembler::assemble_bilinear(space, &[&curl_curl], 4);
    let m_fem = VectorAssembler::assemble_bilinear(space, &[&mass], 4);

    (CsrMatrix::from_fem_csr(k_fem), CsrMatrix::from_fem_csr(m_fem))
}
