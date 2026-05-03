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
pub mod assemble_fem;
pub mod bc;
pub mod postprocess;
pub mod output;

use rem_config::PalaceConfig;
use rem_core::{RemResult, solve_spd, report_peak_memory, CsrMatrix, TripletMatrix};
use rem_parallel::Comm;
use rem_materials::DomainMap;
use rem_mesh::{RemMesh, BoundaryTag, ElementKind, FemSubMesh2d, amr, extract_submesh_tri3, refine_marked_tri3};
use rem_mesh::gmsh::read_msh_file;
use std::path::Path;

/// Entry point called from rem-cli.
pub fn run(config: &PalaceConfig, comm: &dyn Comm) -> RemResult<()> {
    log::info!("\n=== Electrostatic solver ===\n");

    if config.solver.order > 1 {
        log::warn!(
            "Solver.Order={} requested but only P1 (order=1) is implemented; \
             higher-order assembly is pending. Running P1.",
            config.solver.order
        );
    }

    // 1. Load mesh
    let mesh_path = Path::new(&config.model.mesh);
    let raw = read_msh_file(mesh_path)?;
    let mut mesh = RemMesh::from_raw(raw, config)?;
    mesh.set_comm(comm.rank(), comm.size());
    mesh.partition(comm);
    log::info!("Mesh loaded:");
    log::info!("  {} nodes", mesh.n_nodes());
    log::info!("  {} volume elements", mesh.n_volume_elements());
    log::info!("  {} boundary elements", mesh.n_boundary_elements());
    log::info!("");

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
        log::info!("Adaptive mesh refinement (AMR):");
        log::info!("  Max iterations = {}", max_amr_iter);
        log::info!("  Dörfler marking = {:.1}%", amr_theta * 100.0);
        log::info!("");

        let mut cur_mesh = mesh;
        let mut phi = if let Some(exc_tag) = excited_port {
            solve_one(config, &cur_mesh, &domain_map, Some(exc_tag), 1.0, comm)?
        } else {
            solve_one(config, &cur_mesh, &domain_map, None, 0.0, comm)?
        };

        for amr_iter in 1..=max_amr_iter {
            let eta = amr::zz_estimator(&cur_mesh, &phi);
            let total_err: f64 = eta.iter().map(|&e| e * e).sum::<f64>().sqrt();
            log::info!("  Iteration {}: {} nodes, error = {:.3e}", amr_iter, cur_mesh.n_nodes(), total_err);

            let marked = amr::dorfler_mark(&eta, amr_theta);
            if marked.is_empty() {
                log::info!("  → Converged: no elements marked for refinement");
                break;
            }

            let (fine_mesh, midpoints) = refine_amr_mesh(&cur_mesh, &marked);
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
        log::info!("  → Final: {} nodes, error = {:.3e}\n", cur_mesh.n_nodes(), total_err);
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

fn refine_amr_mesh(
    mesh: &RemMesh,
    marked: &[usize],
) -> (RemMesh, std::collections::HashMap<(usize, usize), usize>) {
    if mesh.dim == 2
        && mesh.volume_elements.iter().all(|element| element.kind == ElementKind::Tri3)
        && mesh.boundary_elements.iter().all(|element| element.kind == ElementKind::Line2)
    {
        match refine_marked_tri3(mesh, marked) {
            Ok((fine_mesh, midpoint_map)) => {
                log::info!("AMR refine backend: fem-rs Tri3 bridge");
                return (fine_mesh, midpoint_map);
            }
            Err(err) => {
                log::warn!(
                    "fem-rs Tri3 bridge refinement failed ({}); falling back to legacy AMR",
                    err
                );
            }
        }
    }

    amr::refine_marked(mesh, marked)
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
    // Collect periodic node pairs (empty if no Periodic BCs configured)
    let periodic_pairs = bc::collect_periodic_node_pairs(mesh, config);

    let order     = if config.solver.order >= 2 { 2u8 } else { 1u8 };
    let quad_order = order * 2;

    log::info!(
        "Using fem-assembly backend (order=P{order}, dim={}, periodic_pairs={}).",
        mesh.dim,
        periodic_pairs.len()
    );
    let mut mat = if domain_map.any_anisotropic() {
        log::info!("Anisotropic material(s) detected — tensor assembly.");
        assemble_fem::assemble_stiffness_aniso_fem(mesh, domain_map, order, quad_order)
    } else {
        assemble_fem::assemble_stiffness_fem(mesh, domain_map, order, quad_order)
    };
    if !periodic_pairs.is_empty() {
        mat = remap_periodic_csr(&mat, &periodic_pairs);
    }

    let system_n = mat.nrows;
    let mesh_n = mesh.n_nodes();
    if system_n != mesh_n {
        log::warn!(
            "fem-assembly DOF count ({system_n}) differs from mesh node count ({mesh_n}); \
             using system size for solve and padding output to mesh size"
        );
    }

    let mut rhs = vec![0.0f64; system_n];

    // Apply Dirichlet BCs
    let mut dofs = bc::collect_dirichlet_dofs(mesh, excitation_tag, excitation_val);
    if !periodic_pairs.is_empty() {
        bc::apply_periodic(&mut dofs, &periodic_pairs);
    }
    let original_dof_count = dofs.len();
    dofs.retain(|&idx, _| idx < system_n);
    let dropped = original_dof_count.saturating_sub(dofs.len());
    if dropped > 0 {
        log::warn!(
            "Dropped {dropped} Dirichlet DOFs outside assembled system range [0, {}).",
            system_n
        );
    }
    log::info!("Dirichlet DOFs: {}", dofs.len());
    bc::apply_dirichlet(&mut mat, &mut rhs, &dofs);

    // Solve
    let lin = &config.solver.linear;
    let result = solve_spd(&mat, &rhs, lin.tol, lin.max_iter, comm);
    if result.converged {
        log::info!("PCG converged in {} iterations (|r|={:.2e})", result.iterations, result.residual_norm);
    } else {
        log::warn!(
            "PCG did NOT converge after {} iterations (|r|={:.2e})",
            result.iterations, result.residual_norm
        );
    }

    let mut phi = result.solution;
    if phi.len() != mesh_n {
        let mut phi_full = vec![0.0f64; mesh_n];
        let copy_n = phi.len().min(mesh_n);
        phi_full[..copy_n].copy_from_slice(&phi[..copy_n]);
        phi = phi_full;
    }

    // Propagate periodic DOF values: φ[recv] = φ[donor]
    if !periodic_pairs.is_empty() {
        bc::propagate_periodic(&mut phi, &periodic_pairs);
    }

    Ok(phi)
}

fn remap_periodic_csr(mat: &CsrMatrix, pairs: &[(usize, usize)]) -> CsrMatrix {
    if pairs.is_empty() {
        return mat.clone();
    }

    let mut triplet = TripletMatrix::with_capacity(mat.nrows, mat.ncols, mat.nnz());
    for row in 0..mat.nrows {
        for k in mat.row_ptr[row]..mat.row_ptr[row + 1] {
            triplet.add(row, mat.col_idx[k], mat.values[k]);
        }
    }
    triplet.remap_periodic_nodes(pairs);
    triplet.to_csr()
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
    let domain_energies = domain_energy_records(mesh, domain_map, phi, energy);
    for record in &domain_energies {
        log::info!(
            "Electrostatic energy [domain tag {}, mat {}]: {:.6e} J ({:.2}%)",
            record.domain_tag,
            record.material_index.unwrap_or(usize::MAX),
            record.energy,
            100.0 * record.fraction
        );
    }

    // CSV outputs
    output::write_domain_energy(output_dir, energy)?;
    output::write_domain_energy_by_tag(output_dir, &domain_energies)?;
    if let Some(c) = c_matrix {
        output::write_capacitance_matrix(output_dir, c)?;
    }

    // VTK
    output::write_vtk(output_dir, mesh, phi, &e_field)?;

    Ok(())
}

#[cfg(test)]
fn domain_energy_breakdown(mesh: &RemMesh, domain_map: &DomainMap, phi: &[f64]) -> Vec<(u32, f64)> {
    domain_energy_records(mesh, domain_map, phi, 0.0)
        .into_iter()
        .map(|record| (record.domain_tag, record.energy))
        .collect()
}

struct DomainEnergyRecord {
    domain_tag: u32,
    material_index: Option<usize>,
    energy: f64,
    fraction: f64,
}

fn domain_energy_records(
    mesh: &RemMesh,
    domain_map: &DomainMap,
    phi: &[f64],
    total_energy: f64,
) -> Vec<DomainEnergyRecord> {
    let mut domain_tags: Vec<u32> = mesh.domain_tags.keys().copied().collect();
    domain_tags.sort_unstable();

    if mesh.dim == 2
        && mesh.volume_elements.iter().all(|element| element.kind == ElementKind::Tri3)
        && mesh.boundary_elements.iter().all(|element| element.kind == ElementKind::Line2)
    {
        let mut energies = Vec::new();
        for tag in domain_tags {
            let (material_index, _) = domain_map.get_indexed(tag);
            match extract_domain_submesh(mesh, tag) {
                Some(submesh) => {
                    let sub_phi = submesh.transfer_from_parent(phi);
                    let energy = postprocess::electrostatic_energy(
                        &sub_phi,
                        &submesh.mesh,
                        |sub_tag| domain_map.get(sub_tag).epsilon_abs(),
                    );
                    energies.push(DomainEnergyRecord {
                        domain_tag: tag,
                        material_index: (material_index != usize::MAX).then_some(material_index),
                        energy,
                        fraction: if total_energy.abs() > 1e-300 { energy / total_energy } else { 0.0 },
                    });
                }
                None => energies.push(DomainEnergyRecord {
                    domain_tag: tag,
                    material_index: (material_index != usize::MAX).then_some(material_index),
                    energy: 0.0,
                    fraction: 0.0,
                }),
            }
        }
        return energies;
    }

    domain_tags
        .into_iter()
        .map(|tag| {
            let (material_index, _) = domain_map.get_indexed(tag);
            let energy = postprocess::electrostatic_energy(phi, mesh, |elem_tag| {
                if elem_tag == tag {
                    domain_map.get(elem_tag).epsilon_abs()
                } else {
                    0.0
                }
            });
            DomainEnergyRecord {
                domain_tag: tag,
                material_index: (material_index != usize::MAX).then_some(material_index),
                energy,
                fraction: if total_energy.abs() > 1e-300 { energy / total_energy } else { 0.0 },
            }
        })
        .collect()
}

impl output::DomainEnergyRow for DomainEnergyRecord {
    fn domain_tag(&self) -> u32 {
        self.domain_tag
    }

    fn material_index(&self) -> Option<usize> {
        self.material_index
    }

    fn energy(&self) -> f64 {
        self.energy
    }

    fn fraction(&self) -> f64 {
        self.fraction
    }
}

fn extract_domain_submesh(mesh: &RemMesh, domain_tag: u32) -> Option<FemSubMesh2d> {
    match extract_submesh_tri3(mesh, &[domain_tag]) {
        Ok(submesh) if !submesh.mesh.volume_elements.is_empty() => Some(submesh),
        Ok(_) => None,
        Err(err) => {
            log::warn!(
                "fem-rs Tri3 bridge submesh extraction failed for domain tag {} ({})",
                domain_tag,
                err
            );
            None
        }
    }
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
    use rem_materials::Material;
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
        boundary_tags.insert(11, BoundaryTag::LumpedPort { index: 1, r: 0.0, l: 0.0, c: 0.0 });

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
        let (fine, midpoints) = refine_amr_mesh(&mesh, &all_marked);
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

    #[test]
    fn amr_helper_uses_tri3_bridge_path() {
        let mesh = unit_square_mesh();
        let all_marked: Vec<usize> = (0..mesh.volume_elements.len()).collect();

        let (fine, midpoints) = refine_amr_mesh(&mesh, &all_marked);

        assert_eq!(fine.dim, 2);
        assert!(fine.volume_elements.iter().all(|element| element.kind == ElementKind::Tri3));
        assert!(fine.boundary_elements.iter().all(|element| element.kind == ElementKind::Line2));
        assert_eq!(midpoints.len(), 5, "unit square split should expose five coarse-edge midpoints");
    }

    #[test]
    fn domain_submesh_energy_breakdown_matches_total() {
        use rem_core::constants::EPS0;

        let mut mesh = unit_square_mesh();
        mesh.volume_elements[1].tag = 2;
        mesh.domain_tags = [(1u32, 0usize), (2u32, 0usize)].into_iter().collect();

        let domain_map = DomainMap::from_materials(vec![Material::default()], [(1u32, 0usize)]);
        let phi: Vec<f64> = mesh.nodes.iter().map(|n| n.y).collect();

        let total = postprocess::electrostatic_energy(&phi, &mesh, |_| EPS0);
        let simple_parts = domain_energy_breakdown(&mesh, &domain_map, &phi);
        let parts = domain_energy_records(&mesh, &domain_map, &phi, total);

        assert_eq!(parts.len(), 2);
        assert_eq!(simple_parts.len(), 2);
        let summed: f64 = parts.iter().map(|record| record.energy).sum();
        assert!((summed - total).abs() < 1e-30, "summed={summed:.6e}, total={total:.6e}");
        assert_eq!(parts[0].material_index, Some(0));
        assert_eq!(parts[1].material_index, None);

        let submesh = extract_domain_submesh(&mesh, 2).expect("domain 2 submesh should exist");
        assert!(submesh.mesh.volume_elements.iter().all(|element| element.tag == 2));
    }
}
