//! rem-bem — Laplace BEM solver with adaptive mesh refinement.
//!
//! Solves exterior Laplace problems using the BIE approach.
//!
//! # Multi-conductor capacitance extraction
//!
//! When `Boundaries.Terminal[]` is present, each entry defines a conductor.
//! For each terminal, Dirichlet BIE is solved (φ=1 on one terminal, φ=0 on
//! others), giving the surface charge.  Integrated charge → capacitance matrix.
//!
//! # Adaptive mesh refinement (AMR)
//!
//! When `Solver.Refinement.MaxIter > 0`, a Dörfler-marking loop refines the
//! surface mesh where the BIE residual is largest.

pub mod kernel;
pub mod assemble;
pub mod solve;
pub mod postprocess;

use rem_config::PalaceConfig;
use rem_core::{RemResult, EPS0};
use rem_surface::quadrature::TriQuad;
use rem_parallel::NoComm;

/// Run the BEM solver — single solve or AMR loop.
pub fn run(config: &PalaceConfig) -> RemResult<()> {
    log::info!("\n=== Boundary Element Method (BEM) solver ===\n");

    let output_dir = std::path::Path::new(config.problem.output_dir());
    #[cfg(not(target_arch = "wasm32"))]
    std::fs::create_dir_all(output_dir)?;

    let amr_cfg = &config.model.refinement;
    let max_amr = amr_cfg.max_iter;
    let theta = if amr_cfg.tol > 0.0 { amr_cfg.tol } else { 0.5 };

    // Detect terminals
    let n_terminals = config.boundaries.terminal.len();
    let terminals: Vec<(u32, Vec<u32>)> = if n_terminals > 0 {
        config.boundaries.terminal.iter().map(|t| (t.index, t.attributes.clone())).collect()
    } else {
        let pec_attrs: Vec<u32> = config.boundaries.pec.as_ref()
            .map(|p| p.attributes.clone()).unwrap_or_default();
        if pec_attrs.is_empty() {
            return Err(rem_core::RemError::Config(
                "BEM requires Boundaries.PEC.Attributes or Boundaries.Terminal[]".to_string()));
        }
        vec![(1, pec_attrs)]
    };
    let names: Vec<String> = (0..terminals.len()).map(|i| format!("C{}", i + 1)).collect();

    // Load initial mesh
    let mesh = rem_mesh::load_mesh(config, &NoComm)?;
    let mut cur_mesh = mesh;
    let mut c_matrix = None;

    for amr_iter in 0..=max_amr {
        if amr_iter > 0 {
            log::info!("AMR iteration {} ({} nodes, {} vol elements)...",
                amr_iter, cur_mesh.nodes.len(), cur_mesh.volume_elements.len());
        }

        // Extract surface mesh for conductors
        let all_attrs: Vec<u32> = terminals.iter().flat_map(|(_, a)| a.iter().copied()).collect();
        let surf = rem_surface::surface_mesh::SurfaceMesh::extract(&cur_mesh, &all_attrs)?;
        let n = surf.faces.len();

        // Build face masks
        let term_masks: Vec<Vec<bool>> = terminals.iter().map(|(_, attrs)|
            surf.face_attrs.iter().map(|tag| attrs.contains(tag)).collect()
        ).collect();

        if amr_iter == 0 {
            for (i, (_, attrs)) in terminals.iter().enumerate() {
                let count = term_masks[i].iter().filter(|&&m| m).count();
                log::info!("  Terminal {}: {} faces (attrs {:?})", i, count, attrs);
            }
        }

        // Assemble V and K
        let quad = TriQuad::new(3);
        let (v_mat, k_mat) = crate::assemble::assemble_laplace_p0(&surf, &quad, 4)?;

        // LU factorize V
        let v_lu = v_mat.clone().lu();

        // Solve for each terminal
        let n_cond = terminals.len();
        let mut cm = vec![vec![0.0_f64; n_cond]; n_cond];
        // Collect sigma vectors for error estimation (last solve)
        let mut last_sigma: Option<Vec<f64>> = None;
        let mut last_phi: Option<Vec<f64>> = None;

        for col in 0..n_cond {
            use nalgebra::DVector;
            let phi: Vec<f64> = term_masks[col].iter().map(|&m| if m { 1.0 } else { 0.0 }).collect();
            let mut rhs = DVector::<f64>::zeros(n);
            for i in 0..n {
                rhs[i] = (0..n).map(|j| k_mat[(i, j)] * phi[j]).sum();
            }
            let sigma_vec = v_lu.solve(&rhs)
                .ok_or_else(|| rem_core::RemError::Other(
                    format!("Dirichlet solve failed for terminal {}", col)))?;
            let sigma_slice: Vec<f64> = sigma_vec.iter().copied().collect();
            for row in 0..n_cond {
                let q_norm = postprocess::charge_on_mask(&sigma_slice, &surf, &term_masks[row]);
                cm[row][col] = EPS0 * q_norm;
            }
            if col == n_cond - 1 {
                last_sigma = Some(sigma_slice);
                last_phi = Some(phi);
            }
            log::info!("  Excitation {}: C_self = {:.6e} F ({:.3e} pF)",
                col, cm[col][col], cm[col][col] * 1e12);
        }
        // Symmetrize
        for i in 0..n_cond { for j in i+1..n_cond {
            let avg = (cm[i][j] + cm[j][i]) / 2.0;
            cm[i][j] = avg; cm[j][i] = avg;
        }}
        log::info!("  C matrix [pF]:");
        for row in 0..n_cond {
            let vals: Vec<String> = cm[row].iter().map(|v| format!("{:.3e}", v * 1e12)).collect();
            log::info!("    {}", vals.join("  "));
        }
        c_matrix = Some(cm);

        // AMR: estimate error and refine
        if amr_iter >= max_amr { break; }
        if n == 0 { break; }

        let phi_ref = last_phi.as_ref().unwrap();
        let sigma_ref = last_sigma.as_ref().unwrap();

        // Element-wise residual: r[m] = |Vσ - Kφ|[m] / max(|Vσ|)
        let mut res = vec![0.0_f64; n];
        let mut max_rhs = 0.0_f64;
        for i in 0..n {
            let mut vs = 0.0; let mut kp = 0.0;
            for j in 0..n {
                vs += v_mat[(i, j)] * sigma_ref[j];
                kp += k_mat[(i, j)] * phi_ref[j];
            }
            res[i] = (vs - kp).abs();
            max_rhs = max_rhs.max(vs.abs());
        }
        if max_rhs > 1e-30 {
            for r in &mut res { *r /= max_rhs; }
        }
        let total_err: f64 = res.iter().map(|r| r * r).sum::<f64>().sqrt();
        log::info!("  AMR error = {:.6e}", total_err);

        // Dörfler marking: elements with error > theta × max(error)
        let max_res = res.iter().cloned().fold(0.0_f64, f64::max);
        let marked: Vec<usize> = res.iter().enumerate()
            .filter(|&(_, &r)| r > theta * max_res)
            .map(|(i, _)| i)
            .collect();

        if marked.is_empty() || marked.len() >= n {
            log::info!("  AMR converged ({}/{} marked).", marked.len(), n);
            break;
        }
        log::info!("  Marked {}/{} elements (theta={:.2})", marked.len(), n, theta);

        // For surface-only meshes, volume_elements ARE the BEM faces.
        // Use marked indices directly as volume element indices.
        let (fine_mesh, _midpoints) = rem_mesh::amr::refine_marked(&cur_mesh, &marked);
        log::info!("  Refined: {}→{} nodes, {}→{} elems",
            cur_mesh.nodes.len(), fine_mesh.nodes.len(),
            cur_mesh.volume_elements.len(), fine_mesh.volume_elements.len());
        cur_mesh = fine_mesh;
    }

    // Write final C matrix
    if let Some(ref cm) = c_matrix {
        let csv_path = output_dir.join("bem-capacitance.csv");
        #[cfg(not(target_arch = "wasm32"))]
        {
            let mut wtr = csv::Writer::from_path(&csv_path)
                .map_err(|e| rem_core::RemError::Io(e.into()))?;
            let mut header = vec!["".to_string()];
            for name in &names { header.push(name.clone()); }
            wtr.write_record(&header).map_err(|e| rem_core::RemError::Io(e.into()))?;
            for (i, name) in names.iter().enumerate() {
                let mut row = vec![name.clone()];
                for j in 0..terminals.len() {
                    row.push(format!("{:.12e}", cm[i][j]));
                }
                wtr.write_record(&row).map_err(|e| rem_core::RemError::Io(e.into()))?;
            }
            wtr.flush().map_err(|e| rem_core::RemError::Io(e.into()))?;
        }
        log::info!("Capacitance matrix written to {}", csv_path.display());
    }

    log::info!("\nBEM solve completed.\n");
    Ok(())
}
