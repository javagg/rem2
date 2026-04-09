//! Electrostatic solver — Phase 3.
//!
//! Solves: −∇·(ε ∇φ) = ρ  with Dirichlet BCs on conductors.
//!
//! Pipeline:
//!   1. Load mesh + build domain/boundary maps
//!   2. Assemble stiffness matrix (P1, variable ε)
//!   3. Apply Dirichlet BCs (PEC → φ=0, Ground → φ=0, excitation → φ=V)
//!   4. Solve with PCG + Jacobi preconditioner
//!   5. Recover E = −∇φ (nodal average)
//!   6. Compute electrostatic energy
//!   7. Output CSV + VTK

pub mod assemble;
pub mod bc;
pub mod postprocess;
pub mod output;

use rem_config::PalaceConfig;
use rem_core::{RemResult, solve_pcg, report_peak_memory};
use rem_parallel::Comm;
use rem_materials::DomainMap;
use rem_mesh::{RemMesh, BoundaryTag, amr};
use rem_mesh::gmsh::read_msh_file;
use std::path::Path;

/// Entry point called from rem-cli.
pub fn run(config: &PalaceConfig, comm: &dyn Comm) -> RemResult<()> {
    log::info!("=== Electrostatic solver ===");

    if config.solver.order > 1 {
        log::warn!(
            "Solver.Order={} requested but only P1 (order=1) is implemented; \
             higher-order assembly is pending. Running P1.",
            config.solver.order
        );
    }

    // 1. Load mesh
    let mesh_path = Path::new(&config.model.mesh);
    log::info!("Loading mesh: {}", mesh_path.display());
    let raw = read_msh_file(mesh_path)?;
    let mut mesh = RemMesh::from_raw(raw, config)?;
    mesh.set_comm(comm.rank(), comm.size());
    mesh.partition(comm);
    log::info!(
        "Mesh: {} nodes, {} volume elements, {} boundary elements (dim={})",
        mesh.n_nodes(), mesh.n_volume_elements(), mesh.n_boundary_elements(), mesh.dim
    );

    // 2. Build material map
    let domain_map = DomainMap::from_config(config)?;

    // 3. Identify excited conductor from boundary tags
    let excited_port = find_excited_port(&mesh);
    let output_dir = Path::new(config.problem.output_dir());

    // 4. Solve (with optional AMR loop)
    let amr_cfg = &config.model.refinement;
    let max_amr_iter = if amr_cfg.max_iter > 0 { amr_cfg.max_iter } else { 0 };
    let amr_theta    = if amr_cfg.tol > 0.0 { amr_cfg.tol } else { 0.5 };

    if max_amr_iter > 0 {
        log::info!("AMR enabled: max_iter={}, θ={}", max_amr_iter, amr_theta);
        let mut cur_mesh = mesh;
        let mut phi = if let Some(exc_tag) = excited_port {
            solve_one(config, &cur_mesh, &domain_map, Some(exc_tag), 1.0, comm)?
        } else {
            solve_one(config, &cur_mesh, &domain_map, None, 0.0, comm)?
        };

        for amr_iter in 1..=max_amr_iter {
            let eta = amr::zz_estimator(&cur_mesh, &phi);
            let total_err: f64 = eta.iter().map(|&e| e * e).sum::<f64>().sqrt();
            log::info!("AMR iter {amr_iter}: nodes={}, |η|={total_err:.3e}", cur_mesh.n_nodes());

            let marked = amr::dorfler_mark(&eta, amr_theta);
            if marked.is_empty() {
                log::info!("AMR converged: no elements marked.");
                break;
            }

            let (fine_mesh, midpoints) = amr::refine_marked(&cur_mesh, &marked);
            // Prolongate solution as initial guess (not used for linear Poisson — re-solve)
            let _ = amr::prolongate_p1(&phi, fine_mesh.n_nodes(), &midpoints);

            // Re-solve on fine mesh
            let exc = if let Some(exc_tag) = excited_port { Some(exc_tag) } else { None };
            phi = if let Some(exc_tag) = exc {
                solve_one(config, &fine_mesh, &domain_map, Some(exc_tag), 1.0, comm)?
            } else {
                solve_one(config, &fine_mesh, &domain_map, None, 0.0, comm)?
            };
            cur_mesh = fine_mesh;
        }

        let eta_final = amr::zz_estimator(&cur_mesh, &phi);
        let total_err: f64 = eta_final.iter().map(|&e| e * e).sum::<f64>().sqrt();
        log::info!("AMR final: nodes={}, |η|={total_err:.3e}", cur_mesh.n_nodes());
        finalize(config, &cur_mesh, &domain_map, &phi, output_dir, None)?;
    } else if let Some(exc_tag) = excited_port {
        log::info!("Excited port tag: {}", exc_tag);
        let phi = solve_one(config, &mesh, &domain_map, Some(exc_tag), 1.0, comm)?;
        finalize(config, &mesh, &domain_map, &phi, output_dir, None)?;
    } else {
        let phi = solve_one(config, &mesh, &domain_map, None, 0.0, comm)?;
        finalize(config, &mesh, &domain_map, &phi, output_dir, None)?;
    }

    report_peak_memory("Electrostatic solver");
    Ok(())
}

/// Solve a single electrostatic problem with one conductor excited.
pub fn solve_one(
    config: &PalaceConfig,
    mesh: &RemMesh,
    domain_map: &DomainMap,
    excitation_tag: Option<u32>,
    excitation_val: f64,
    comm: &dyn Comm,
) -> RemResult<Vec<f64>> {
    let n = mesh.n_nodes();

    // Collect periodic node pairs (empty if no Periodic BCs configured)
    let periodic_pairs = bc::collect_periodic_node_pairs(mesh, config);

    // Assemble stiffness matrix — scalar or tensor path
    let mut triplet = if domain_map.any_anisotropic() {
        log::info!("Anisotropic material(s) detected — using tensor stiffness assembly.");
        let tensor_fn = |tag: u32| domain_map.get(tag).epsilon_tensor;
        assemble::assemble_stiffness_aniso(mesh, tensor_fn)?
    } else {
        let eps_fn = |tag: u32| domain_map.get(tag).epsilon_abs();
        assemble::assemble_stiffness(mesh, eps_fn)?
    };

    // Apply periodic remapping before converting to CSR
    if !periodic_pairs.is_empty() {
        triplet.remap_periodic_nodes(&periodic_pairs);
    }

    let mut mat = triplet.to_csr();
    let mut rhs = vec![0.0f64; n];

    // Apply Dirichlet BCs
    let mut dofs = bc::collect_dirichlet_dofs(mesh, excitation_tag, excitation_val);
    if !periodic_pairs.is_empty() {
        bc::apply_periodic(&mut dofs, &periodic_pairs);
    }
    log::info!("Dirichlet DOFs: {}", dofs.len());
    bc::apply_dirichlet(&mut mat, &mut rhs, &dofs);

    // Solve
    let lin = &config.solver.linear;
    let result = solve_pcg(&mat, &rhs, lin.tol, lin.max_iter, comm);
    if result.converged {
        log::info!("PCG converged in {} iterations (|r|={:.2e})", result.iterations, result.residual_norm);
    } else {
        log::warn!(
            "PCG did NOT converge after {} iterations (|r|={:.2e})",
            result.iterations, result.residual_norm
        );
    }

    let mut phi = result.solution;

    // Propagate periodic DOF values: φ[recv] = φ[donor]
    if !periodic_pairs.is_empty() {
        bc::propagate_periodic(&mut phi, &periodic_pairs);
    }

    Ok(phi)
}

/// Post-process and write output files.
fn finalize(
    _config: &PalaceConfig,
    mesh: &RemMesh,
    domain_map: &DomainMap,
    phi: &[f64],
    output_dir: &Path,
    c_matrix: Option<&[Vec<f64>]>,
) -> RemResult<()> {
    let eps_fn = |tag: u32| domain_map.get(tag).epsilon_abs();

    // E-field recovery
    let e_field = postprocess::gradient_recovery(phi, mesh);

    // Electrostatic energy
    let energy = postprocess::electrostatic_energy(phi, mesh, eps_fn);
    log::info!("Electrostatic energy: {:.6e} J", energy);

    // CSV outputs
    output::write_domain_energy(output_dir, energy)?;
    if let Some(c) = c_matrix {
        output::write_capacitance_matrix(output_dir, c)?;
    }

    // VTK
    output::write_vtk(output_dir, mesh, phi, &e_field)?;

    Ok(())
}

/// Return the INDEX of the first Terminal, LumpedPort, or WavePort boundary.
fn find_excited_port(mesh: &RemMesh) -> Option<u32> {
    for bc in mesh.boundary_tags.values() {
        match bc {
            BoundaryTag::Terminal { index } => return Some(*index),
            BoundaryTag::LumpedPort { index, .. } => return Some(*index),
            BoundaryTag::WavePort { index } => return Some(*index),
            _ => {}
        }
    }
    None
}

// ---------------------------------------------------------------------------
// Integration tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use rem_config::{load_config_from_str, ConfigFormat};
    use rem_mesh::{Node, Element, ElementKind};
    use rem_parallel::NoComm;
    use std::collections::HashMap;

    /// Unit square: 4 nodes, 2 triangles (tag=1).
    /// Bottom edge (y=0) nodes 0,1 → tag=10 (Ground)
    /// Top    edge (y=1) nodes 2,3 → tag=11 (LumpedPort index=1)
    fn unit_square_mesh() -> RemMesh {
        let nodes = vec![
            Node { id: 0, x: 0.0, y: 0.0, z: 0.0 },
            Node { id: 1, x: 1.0, y: 0.0, z: 0.0 },
            Node { id: 2, x: 1.0, y: 1.0, z: 0.0 },
            Node { id: 3, x: 0.0, y: 1.0, z: 0.0 },
        ];
        let volume_elements = vec![
            Element { id: 1, kind: ElementKind::Tri3, tag: 1, node_ids: vec![0, 1, 2] , rank: 0 },
            Element { id: 2, kind: ElementKind::Tri3, tag: 1, node_ids: vec![0, 2, 3] , rank: 0 },
        ];
        let boundary_elements = vec![
            Element { id: 3, kind: ElementKind::Line2, tag: 10, node_ids: vec![0, 1] , rank: 0 },
            Element { id: 4, kind: ElementKind::Line2, tag: 11, node_ids: vec![2, 3] , rank: 0 },
        ];
        let mut boundary_tags: HashMap<u32, BoundaryTag> = HashMap::new();
        boundary_tags.insert(10, BoundaryTag::Ground);
        boundary_tags.insert(11, BoundaryTag::LumpedPort { index: 1, r: 0.0 });

        RemMesh {
            nodes,
            volume_elements,
            boundary_elements,
            domain_tags: Default::default(),
            boundary_tags,
            dim: 2,
            rank: 0,
            size: 1,
        }
    }

    fn default_config() -> PalaceConfig {
        load_config_from_str(
            r#"{"Problem":{"Type":"Electrostatic"},"Model":{"Mesh":"x.msh"},
                "Solver":{"Linear":{"Tol":1e-12,"MaxIter":200}}}"#,
            ConfigFormat::Json,
        )
        .unwrap()
    }

    #[test]
    fn parallel_plate_linear_phi() {
        // Exact solution: φ(x,y) = y  (Poisson in 2D, linear E-field)
        let mesh = unit_square_mesh();
        let config = default_config();
        let domain_map = DomainMap::from_config(&config).unwrap();

        // Excite LumpedPort index=1 (top edge, physical tag 11) with V=1
        let phi = solve_one(&config, &mesh, &domain_map, Some(1), 1.0, &NoComm).unwrap();

        for (i, node) in mesh.nodes.iter().enumerate() {
            let exact = node.y;
            let err = (phi[i] - exact).abs();
            assert!(
                err < 1e-10,
                "node {} ({:.1},{:.1}): phi={:.6}, exact={:.6}, err={:.2e}",
                i, node.x, node.y, phi[i], exact, err
            );
        }
    }

    #[test]
    fn e_field_recovery_unit_square() {
        // φ = y → E = (0, -1, 0) everywhere
        let mesh = unit_square_mesh();
        let phi: Vec<f64> = mesh.nodes.iter().map(|n| n.y).collect();
        let e = postprocess::gradient_recovery(&phi, &mesh);
        for (i, ev) in e.iter().enumerate() {
            assert!(ev[0].abs() < 1e-12, "node {}: E_x={:.2e}", i, ev[0]);
            assert!((ev[1] + 1.0).abs() < 1e-12, "node {}: E_y={:.6}", i, ev[1]);
        }
    }

    #[test]
    fn energy_parallel_plate_vacuum() {
        // φ=y, ε=ε₀, unit square → U = ε₀/2  [J/m (2D per unit depth)]
        use rem_core::constants::EPS0;
        let mesh = unit_square_mesh();
        let phi: Vec<f64> = mesh.nodes.iter().map(|n| n.y).collect();
        let energy = postprocess::electrostatic_energy(&phi, &mesh, |_| EPS0);
        assert!(
            (energy - EPS0 / 2.0).abs() < 1e-30,
            "energy={:.6e}, expected={:.6e}", energy, EPS0 / 2.0
        );
    }

    #[test]
    fn energy_dielectric_medium() {
        // Same as above but ε = 4.5 ε₀ (SiO₂)
        use rem_core::constants::EPS0;
        let mesh = unit_square_mesh();
        let phi: Vec<f64> = mesh.nodes.iter().map(|n| n.y).collect();
        let eps_r = 4.5_f64;
        let energy = postprocess::electrostatic_energy(&phi, &mesh, |_| eps_r * EPS0);
        assert!(
            (energy - eps_r * EPS0 / 2.0).abs() < 1e-28,
            "energy={:.6e}, expected={:.6e}", energy, eps_r * EPS0 / 2.0
        );
    }

    #[test]
    fn amr_loop_refines_and_resolves() {
        // AMR loop: force refinement by marking all elements, verify fine mesh + re-solve
        use rem_config::{load_config_from_str, ConfigFormat};
        use rem_mesh::amr;

        let mesh = unit_square_mesh();
        let config = load_config_from_str(
            r#"{"Problem":{"Type":"Electrostatic"},"Model":{"Mesh":"x.msh"},
                "Solver":{"Linear":{"Tol":1e-12,"MaxIter":200}}}"#,
            ConfigFormat::Json,
        ).unwrap();
        let domain_map = DomainMap::from_config(&config).unwrap();

        // Mark all elements explicitly (simulate AMR trigger)
        let all_marked: Vec<usize> = (0..mesh.volume_elements.len()).collect();
        let (fine, midpoints) = amr::refine_marked(&mesh, &all_marked);
        assert!(fine.n_nodes() > mesh.n_nodes(),
            "red refinement should add midpoint nodes");
        assert!(fine.volume_elements.len() >= 4 * mesh.volume_elements.len(),
            "each Tri3 should produce 4 children");

        // Re-solve on fine mesh: φ=y should still be exact
        let phi1 = solve_one(&config, &fine, &domain_map, Some(1), 1.0, &NoComm).unwrap();
        for (i, node) in fine.nodes.iter().enumerate() {
            assert!((phi1[i] - node.y).abs() < 1e-8,
                "node {i} at y={:.3}: φ={:.6}, err={:.2e}",
                node.y, phi1[i], (phi1[i] - node.y).abs());
        }

        // Prolongation preserves linear fields
        let phi0: Vec<f64> = mesh.nodes.iter().map(|n| n.y).collect();
        let phi_prolonged = amr::prolongate_p1(&phi0, fine.n_nodes(), &midpoints);
        for (i, node) in fine.nodes.iter().enumerate() {
            assert!((phi_prolonged[i] - node.y).abs() < 1e-12,
                "prolongated node {i}: expected y={:.3}, got {:.6}", node.y, phi_prolonged[i]);
        }
    }
}
