use std::collections::HashMap;
use std::f64::consts::PI;

use fem_assembly::coefficient::CtxFnCoeff;
use fem_assembly::standard::{CurlCurlIntegrator, VectorMassIntegrator};
use fem_assembly::{VectorAssembler, VectorBoundaryAssembler, TangentialMassIntegrator};
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

    let q_factors = compute_q_factors_hcurl(&m_mat, &eigenpairs, config, mesh, domain_map, n, comm)?;

    Ok(EigenResult {
        frequencies_hz: eigenpairs.iter().map(|(f, _)| *f).collect(),
        eigenvectors: eigenpairs.into_iter().map(|(_, v)| v).collect(),
        q_factors,
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

/// Build an HCurl bilinear form matrix with a scalar coefficient and return as CSR.
fn assemble_hcurl_form<C>(
    mesh: &RemMesh,
    order: u8,
    coeff: C,
) -> RemResult<CsrMatrix>
where
    C: Fn(&fem_assembly::coefficient::CoeffCtx<'_>) -> f64 + Send + Sync + 'static,
{
    let cf = CtxFnCoeff(coeff);
    let mass = VectorMassIntegrator { alpha: cf };
    let m_fem = if mesh.dim == 2 {
        let simplex = mesh.to_simplex_mesh_2d();
        let space = HCurlSpace::new(simplex, order);
        VectorAssembler::assemble_bilinear(&space, &[&mass], 4)
    } else if mesh.dim == 3 {
        let simplex = mesh.to_simplex_mesh();
        let space = HCurlSpace::new(simplex, order);
        VectorAssembler::assemble_bilinear(&space, &[&mass], 4)
    } else {
        return Err(RemError::Config(format!(
            "HCurl assemble_hcurl_form only supports 2-D/3-D, got dim={}", mesh.dim
        )));
    };
    Ok(CsrMatrix::from_fem_csr(m_fem))
}

/// Assemble the tangential mass matrix on boundary tags (for conductor loss).
fn assemble_hcurl_surface_matrix(
    mesh: &RemMesh,
    order: u8,
    boundary_tags: &[i32],
) -> RemResult<CsrMatrix> {
    let integ = TangentialMassIntegrator { gamma: 1.0 };
    let k_surf = if mesh.dim == 2 {
        let simplex = mesh.to_simplex_mesh_2d();
        let space = HCurlSpace::new(simplex, order);
        VectorBoundaryAssembler::assemble_boundary_bilinear(
            &space, &[&integ], boundary_tags, 4,
        )
    } else if mesh.dim == 3 {
        let simplex = mesh.to_simplex_mesh();
        let space = HCurlSpace::new(simplex, order);
        VectorBoundaryAssembler::assemble_boundary_bilinear(
            &space, &[&integ], boundary_tags, 4,
        )
    } else {
        return Err(RemError::Config(format!(
            "HCurl surface matrix only supports 2-D/3-D, got dim={}", mesh.dim
        )));
    };
    Ok(CsrMatrix::from_fem_csr(k_surf))
}

/// Compute Q-factors for HCurl eigenmodes using perturbation.
///
/// Dielectric: Q_d = (xᵀ M x) / (xᵀ M_loss x)  where M_loss uses ε·tanδ.
/// Conductor:  1/Q_c = R_s/(ωμ₀)·(xᵀ K_surf x)/(xᵀ M x),  R_s = √(ωμ₀/(2σ_wall)).
/// Combined:   1/Q_total = 1/Q_d + 1/Q_c
fn compute_q_factors_hcurl(
    m_mat: &CsrMatrix,
    eigenpairs: &[(f64, Vec<f64>)],
    config: &PalaceConfig,
    mesh: &RemMesh,
    domain_map: &DomainMap,
    n: usize,
    comm: &dyn Comm,
) -> RemResult<Option<Vec<f64>>> {
    let order = config.solver.eigenmode_hcurl_order().clamp(1, 2) as u8;

    // ── Dielectric loss: M_loss with ε·tanδ coefficient ──────────────────
    let has_lossy = mesh.domain_tags.keys().any(|&tag| domain_map.get(tag).is_lossy());
    let m_loss: Option<CsrMatrix> = if has_lossy {
        // Pre-collect (ε·tanδ) per element tag into an owned Vec (avoids capturing &DomainMap).
        let max_tag = mesh.domain_tags.keys().copied().max().unwrap_or(0) as usize;
        let mut loss_coeffs = vec![0.0f64; max_tag + 1];
        for (&tag, _) in &mesh.domain_tags {
            let mat = domain_map.get(tag);
            loss_coeffs[tag as usize] = mat.epsilon_abs() * mat.loss_tangent;
        }
        Some(assemble_hcurl_form(mesh, order, move |ctx| {
            let tid = ctx.elem_tag as usize;
            *loss_coeffs.get(tid).unwrap_or(&0.0)
        })?)
    } else {
        None
    };

    // ── Conductor loss: tangential mass on PEC/Ground boundaries ─────────
    let sigma_wall = config.solver.eigenmode.as_ref()
        .map(|e| e.wall_conductivity)
        .unwrap_or(0.0);

    let has_conductor = sigma_wall > 0.0
        && mesh.boundary_tags.values().any(|bc| {
            matches!(bc, BoundaryTag::Pec | BoundaryTag::Ground)
        });

    let k_surf: Option<CsrMatrix> = if has_conductor {
        let pec_tags: Vec<i32> = mesh.boundary_tags.iter()
            .filter_map(|(tag, bc)| {
                if matches!(bc, BoundaryTag::Pec | BoundaryTag::Ground) {
                    Some(*tag as i32)
                } else {
                    None
                }
            })
            .collect();
        let ks = assemble_hcurl_surface_matrix(mesh, order, &pec_tags)?;
        log::info!(
            "HCurl conductor Q: assembled {}×{} tangential mass matrix on {} PEC tags",
            ks.nrows, ks.ncols, pec_tags.len()
        );
        Some(ks)
    } else {
        None
    };

    if m_loss.is_none() && k_surf.is_none() {
        return Ok(None);
    }

    let mu0: f64 = 1.25663706212e-6;
    let freqs: Vec<f64> = eigenpairs.iter().map(|(f, _)| *f).collect();

    let qs: Vec<f64> = eigenpairs.iter().zip(freqs.iter()).map(|((freq_hz, phi), _)| {
        let omega = 2.0 * PI * freq_hz;

        let mut m_phi = vec![0.0f64; n];
        m_mat.matvec(phi, &mut m_phi, comm);
        let denom: f64 = phi.iter().zip(m_phi.iter()).map(|(a, b)| a * b).sum();
        let denom_safe = if denom.abs() > 1e-300 { denom } else { 1.0 };

        // Dielectric contribution: 1/Q_d = (xᵀ M_loss x) / (xᵀ M x)
        let inv_q_diel = m_loss.as_ref().map(|ml| {
            let mut ml_phi = vec![0.0f64; n];
            ml.matvec(phi, &mut ml_phi, comm);
            let num: f64 = phi.iter().zip(ml_phi.iter()).map(|(a, b)| a * b).sum();
            if num.abs() > 1e-300 { num / denom_safe } else { 0.0 }
        }).unwrap_or(0.0);

        // Conductor contribution: 1/Q_c = R_s/(ωμ₀)·(xᵀ K_surf x)/(xᵀ M x)
        let inv_q_cond = k_surf.as_ref().map(|ks| {
            let mut ks_phi = vec![0.0f64; n];
            ks.matvec(phi, &mut ks_phi, comm);
            let surf: f64 = phi.iter().zip(ks_phi.iter()).map(|(a, b)| a * b).sum();
            if omega > 0.0 && denom.abs() > 1e-300 && surf > 1e-300 {
                let r_s = (omega * mu0 / (2.0 * sigma_wall)).sqrt();
                (r_s * surf) / (omega * mu0 * denom_safe)
            } else {
                0.0
            }
        }).unwrap_or(0.0);

        let inv_total = inv_q_diel + inv_q_cond;
        if inv_total > 1e-300 { 1.0 / inv_total } else { f64::INFINITY }
    }).collect();

    if has_conductor {
        let f0 = eigenpairs.first().map(|(f, _)| *f).unwrap_or(1e9);
        let omega0 = 2.0 * PI * f0;
        let r_s0 = (omega0 * mu0 / (2.0 * sigma_wall)).sqrt();
        log::info!(
            "HCurl Q-factors: dielectric + conductor (σ_wall={:.3e} S/m, R_s={:.4} mΩ/□ @ {:.3} GHz)",
            sigma_wall, r_s0 * 1e3, f0 / 1e9,
        );
    } else if m_loss.is_some() {
        log::info!("HCurl Q-factors: dielectric loss only");
    }

    Ok(Some(qs))
}
