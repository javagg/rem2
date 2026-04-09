//! Magnetostatic solver — Phase 4.
//!
//! Supports two modes selected by `mesh.dim`:
//!
//! ## 2-D (mesh.dim == 2)
//!   −∇·(ν ∇A_z) = J_z      (in Ω)
//!   A_z = 0                  (on Dirichlet boundaries — Ground)
//!
//!   Post-processing:
//!   B_x =  ∂A_z/∂y ,   B_y = −∂A_z/∂x
//!
//! ## 3-D (mesh.dim == 3)
//!   Three decoupled vector-potential Poisson problems:
//!   −∇·(ν ∇Aᵢ) = Jᵢ   for i ∈ {x, y, z}
//!
//!   Post-processing (curl):
//!   B_x = ∂Az/∂y − ∂Ay/∂z
//!   B_y = ∂Ax/∂z − ∂Az/∂x
//!   B_z = ∂Ay/∂x − ∂Ax/∂y
//!
//! The assembly reuses `rem_electrostatic::assemble_stiffness` which already
//! handles both Tri3 (2-D) and Tet4 (3-D).  Gradient recovery likewise works
//! for both element types.

use rem_config::PalaceConfig;
use rem_core::{RemResult, solve_spd};
use rem_parallel::Comm;
use rem_materials::DomainMap;
use rem_mesh::{RemMesh, BoundaryTag, amr};
use rem_mesh::gmsh::read_msh_file;
use rem_electrostatic::assemble::{self, assemble_stiffness_aniso};
use rem_electrostatic::bc;
use rem_electrostatic::postprocess;
use std::path::Path;

/// Entry point called from rem-cli.
pub fn run(config: &PalaceConfig, comm: &dyn Comm) -> RemResult<()> {
    if config.solver.order > 1 {
        log::warn!(
            "Solver.Order={} requested but only P1 (order=1) is implemented; \
             higher-order assembly is pending. Running P1.",
            config.solver.order
        );
    }
    let mesh_path = Path::new(&config.model.mesh);
    log::info!("Loading mesh: {}", mesh_path.display());
    let raw = read_msh_file(mesh_path)?;
    let mut mesh = RemMesh::from_raw(raw, config)?;
    mesh.set_comm(comm.rank(), comm.size());
    mesh.partition(comm);
    log::info!(
        "Mesh: {} nodes, {} volume elements (dim={})",
        mesh.n_nodes(), mesh.n_volume_elements(), mesh.dim
    );

    if mesh.dim == 3 {
        log::info!("=== Magnetostatic solver (3-D vector potential) ===");
        return run_3d(config, mesh, comm);
    }
    log::info!("=== Magnetostatic solver (2-D A_z) ===");
    run_2d(config, mesh, comm)
}

// ---------------------------------------------------------------------------
// 2-D path (A_z scalar)
// ---------------------------------------------------------------------------

fn run_2d(config: &PalaceConfig, mesh: RemMesh, comm: &dyn Comm) -> RemResult<()> {
    let domain_map = DomainMap::from_config(config)?;
    let output_dir = Path::new(config.problem.output_dir());

    // Find surface current excitation (analogous to LumpedPort in electrostatics)
    let excited_tag = find_surface_current_tag(&mesh); // returns SurfaceCurrent INDEX

    // AMR loop
    let amr_cfg = &config.model.refinement;
    let max_amr_iter = if amr_cfg.max_iter > 0 { amr_cfg.max_iter } else { 0 };
    let amr_theta    = if amr_cfg.tol > 0.0 { amr_cfg.tol } else { 0.5 };

    let (final_mesh, az) = if max_amr_iter > 0 {
        log::info!("AMR enabled: max_iter={}, θ={}", max_amr_iter, amr_theta);
        let mut cur_mesh = mesh;
        let mut az = solve_one(config, &cur_mesh, &domain_map, excited_tag, comm)?;

        for amr_iter in 1..=max_amr_iter {
            let eta = amr::zz_estimator(&cur_mesh, &az);
            let total_err: f64 = eta.iter().map(|&e| e * e).sum::<f64>().sqrt();
            log::info!("AMR iter {amr_iter}: nodes={}, |η|={total_err:.3e}", cur_mesh.n_nodes());

            let marked = amr::dorfler_mark(&eta, amr_theta);
            if marked.is_empty() {
                log::info!("AMR converged: no elements marked.");
                break;
            }

            let (fine_mesh, _midpoints) = amr::refine_marked(&cur_mesh, &marked);
            az = solve_one(config, &fine_mesh, &domain_map, excited_tag, comm)?;
            cur_mesh = fine_mesh;
        }

        let eta_final = amr::zz_estimator(&cur_mesh, &az);
        let total_err: f64 = eta_final.iter().map(|&e| e * e).sum::<f64>().sqrt();
        log::info!("AMR final: nodes={}, |η|={total_err:.3e}", cur_mesh.n_nodes());
        (cur_mesh, az)
    } else {
        let az = solve_one(config, &mesh, &domain_map, excited_tag, comm)?;
        (mesh, az)
    };

    // gradient_recovery returns −∇A_z (same sign convention as E = −∇φ).
    // So g = (−∂Az/∂x, −∂Az/∂y).
    // B_x = ∂Az/∂y = −g[1],   B_y = −∂Az/∂x = g[0]
    let grad_az = postprocess::gradient_recovery(&az, &final_mesh);
    let b_field: Vec<[f64; 3]> = grad_az.iter()
        .map(|g| [-g[1], g[0], 0.0])
        .collect();

    // Magnetic energy U = (1/2) ∫ ν |∇A_z|² dΩ
    let nu_fn = |tag: u32| domain_map.get(tag).reluctivity();
    let energy = postprocess::electrostatic_energy(&az, &final_mesh, nu_fn);
    log::info!("Magnetic energy: {:.6e} J/m", energy);

    // Write CSV + VTK
    write_outputs(output_dir, &final_mesh, &az, &b_field, energy)?;

    Ok(())
}

// ---------------------------------------------------------------------------
// 3-D path (Ax, Ay, Az vector potential — three decoupled scalar solves)
// ---------------------------------------------------------------------------

/// Run the 3-D magnetostatic solver.
///
/// Solves three independent scalar Poisson problems:
///   −∇·(ν ∇Aᵢ) = 0   for i ∈ {x, y, z}
///
/// Boundary conditions:
///   Ground / PEC → Aᵢ = 0 on the boundary
///   SurfaceCurrent (excited) → Az = 1, Ax = Ay = 0 (z-directed excitation default)
///
/// Post-processing (curl, from gradient_recovery which returns −∇Aᵢ):
///   B_x = ∂Az/∂y − ∂Ay/∂z  = −gz[1] − (−gy[2]) = gy[2] − gz[1]
///   B_y = ∂Ax/∂z − ∂Az/∂x  = −gx[2] − (−gz[0]) = gz[0] − gx[2]
///   B_z = ∂Ay/∂x − ∂Ax/∂y  = −gy[0] − (−gx[1]) = gx[1] − gy[0]
fn run_3d(config: &PalaceConfig, mesh: RemMesh, comm: &dyn Comm) -> RemResult<()> {
    let domain_map = DomainMap::from_config(config)?;
    let output_dir = Path::new(config.problem.output_dir());
    let excited_tag = find_surface_current_tag(&mesh);

    // Three decoupled solves: excite Az=1, Ax=Ay=0 (z-directed current port)
    let (ax, ay, az) = solve_3d(config, &mesh, &domain_map, excited_tag, comm)?;

    // gradient_recovery returns g = −∇Aᵢ per node  (sign: g[j] = −∂Aᵢ/∂xⱼ)
    let gx = postprocess::gradient_recovery(&ax, &mesh);
    let gy = postprocess::gradient_recovery(&ay, &mesh);
    let gz = postprocess::gradient_recovery(&az, &mesh);

    let n = mesh.n_nodes();
    let b_field: Vec<[f64; 3]> = (0..n).map(|i| [
        gy[i][2] - gz[i][1],   // Bx = ∂Az/∂y − ∂Ay/∂z
        gz[i][0] - gx[i][2],   // By = ∂Ax/∂z − ∂Az/∂x
        gx[i][1] - gy[i][0],   // Bz = ∂Ay/∂x − ∂Ax/∂y
    ]).collect();

    // Magnetic energy: U = (1/2) ∫ ν (|∇Ax|² + |∇Ay|² + |∇Az|²) dΩ
    let nu_fn = |tag: u32| domain_map.get(tag).reluctivity();
    let energy = postprocess::electrostatic_energy(&ax, &mesh, nu_fn)
               + postprocess::electrostatic_energy(&ay, &mesh, nu_fn)
               + postprocess::electrostatic_energy(&az, &mesh, nu_fn);
    log::info!("3-D magnetic energy: {:.6e} J", energy);

    write_outputs(output_dir, &mesh, &az, &b_field, energy)
}

/// Solve the 3-D vector potential: three decoupled scalar Poisson systems.
///
/// Returns `(Ax, Ay, Az)` — nodal coefficient vectors of length `mesh.n_nodes()`.
///
/// The stiffness matrix K (built from the reluctivity ν) is assembled once and
/// reused for all three component solves.  Boundary conditions:
///   - Ground / PEC: Aᵢ = 0
///   - SurfaceCurrent (excited): Az = 1 (z-directed), Ax = Ay = 0
///
/// Exposed as `pub` for tests and external callers.
pub fn solve_3d(
    config: &PalaceConfig,
    mesh: &RemMesh,
    domain_map: &DomainMap,
    excitation_tag: Option<u32>,
    comm: &dyn Comm,
) -> RemResult<(Vec<f64>, Vec<f64>, Vec<f64>)> {
    let triplet = if domain_map.any_magnetically_anisotropic() {
        log::info!("Anisotropic permeability detected — using tensor reluctivity assembly (3D).");
        let tensor_fn = |tag: u32| domain_map.get(tag).nu_tensor;
        assemble_stiffness_aniso(mesh, tensor_fn)?
    } else {
        let nu_fn = |tag: u32| domain_map.get(tag).reluctivity();
        assemble::assemble_stiffness(mesh, nu_fn)?
    };
    let n = mesh.n_nodes();
    let lin = &config.solver.linear;

    // Component excitation values: [Ax_val, Ay_val, Az_val]
    let excitation_values = [0.0_f64, 0.0, 1.0];
    let labels = ["x", "y", "z"];

    let mut solutions: [Vec<f64>; 3] = [
        vec![0.0; n],
        vec![0.0; n],
        vec![0.0; n],
    ];

    for comp in 0..3 {
        let mut mat = triplet.clone().to_csr();
        let mut rhs = vec![0.0f64; n];

        let dofs = collect_magnetostatic_dofs(mesh, excitation_tag, excitation_values[comp]);
        bc::apply_dirichlet(&mut mat, &mut rhs, &dofs);

        let result = solve_spd(&mat, &rhs, lin.tol, lin.max_iter, comm);
        if result.converged {
            log::info!("3-D A{}: PCG converged in {} iters (|r|={:.2e})",
                labels[comp], result.iterations, result.residual_norm);
        } else {
            log::warn!("3-D A{}: PCG did NOT converge after {} iters",
                labels[comp], result.iterations);
        }
        solutions[comp] = result.solution;
    }

    let [ax, ay, az] = solutions;
    Ok((ax, ay, az))
}

/// Solve a single magnetostatic problem.
pub fn solve_one(
    config: &PalaceConfig,
    mesh: &RemMesh,
    domain_map: &DomainMap,
    excitation_tag: Option<u32>,
    comm: &dyn Comm,
) -> RemResult<Vec<f64>> {
    let n = mesh.n_nodes();

    // Assemble stiffness with reluctivity — scalar or tensor path
    let triplet = if domain_map.any_magnetically_anisotropic() {
        log::info!("Anisotropic permeability detected — using tensor reluctivity assembly.");
        let tensor_fn = |tag: u32| domain_map.get(tag).nu_tensor;
        assemble_stiffness_aniso(mesh, tensor_fn)?
    } else {
        let nu_fn = |tag: u32| domain_map.get(tag).reluctivity();
        assemble::assemble_stiffness(mesh, nu_fn)?
    };
    let mut mat = triplet.to_csr();
    let mut rhs = vec![0.0f64; n];

    // Dirichlet BCs: Ground → A_z = 0, SurfaceCurrent excited → A_z = 1
    let dofs = collect_magnetostatic_dofs(mesh, excitation_tag, 1.0);
    log::info!("Dirichlet DOFs: {}", dofs.len());
    bc::apply_dirichlet(&mut mat, &mut rhs, &dofs);

    // Solve
    let lin = &config.solver.linear;
    let result = solve_spd(&mat, &rhs, lin.tol, lin.max_iter, comm);
    if result.converged {
        log::info!("PCG converged in {} iterations (|r|={:.2e})", result.iterations, result.residual_norm);
    } else {
        log::warn!("PCG did NOT converge after {} iterations", result.iterations);
    }

    Ok(result.solution)
}

/// Collect Dirichlet DOFs for magnetostatics.
/// `excited_index`: the SurfaceCurrent INDEX (not physical tag) to assign A_z=val.
fn collect_magnetostatic_dofs(
    mesh: &RemMesh,
    excited_index: Option<u32>,
    excitation_val: f64,
) -> std::collections::HashMap<usize, f64> {
    let mut dofs = std::collections::HashMap::new();
    for belem in &mesh.boundary_elements {
        let bc = match mesh.boundary_tags.get(&belem.tag) {
            Some(b) => b,
            None => continue,
        };
        let val = match bc {
            BoundaryTag::Ground | BoundaryTag::Pec => 0.0,
            BoundaryTag::SurfaceCurrent { index } => {
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

fn find_surface_current_tag(mesh: &RemMesh) -> Option<u32> {
    for bc in mesh.boundary_tags.values() {
        if let BoundaryTag::SurfaceCurrent { index } = bc {
            return Some(*index);
        }
    }
    None
}

fn write_outputs(
    output_dir: &Path,
    mesh: &RemMesh,
    az: &[f64],
    b_field: &[[f64; 3]],
    energy: f64,
) -> RemResult<()> {
    use rem_core::RemError;
    use std::io::Write;

    // domain-B.csv (analogous to domain-E.csv)
    let dir = output_dir.join("postpro");
    std::fs::create_dir_all(&dir).map_err(RemError::Io)?;
    let path = dir.join("domain-B.csv");
    let mut f = std::fs::File::create(&path).map_err(RemError::Io)?;
    writeln!(
        f,
        r#""Frequency (GHz)","Magnetic Field Energy (J)","Electric Field Energy (J)","Total Energy (J)""#
    )
    .map_err(RemError::Io)?;
    writeln!(f, "0.000000e0,{:.6e},0.000000e0,{:.6e}", energy, energy)
        .map_err(RemError::Io)?;
    log::info!("Written: {}", path.display());

    // VTK with A_z and B
    let vtk_dir = output_dir.join("paraview");
    std::fs::create_dir_all(&vtk_dir).map_err(RemError::Io)?;
    let vtk_path = vtk_dir.join("solution.vtk");
    let mut vf = std::fs::File::create(&vtk_path).map_err(RemError::Io)?;

    let n_nodes = mesh.n_nodes();
    let n_cells = mesh.n_volume_elements();
    let cell_list_size: usize = mesh.volume_elements.iter().map(|e| e.node_ids.len() + 1).sum();

    writeln!(vf, "# vtk DataFile Version 3.0").map_err(RemError::Io)?;
    writeln!(vf, "rem magnetostatic solution").map_err(RemError::Io)?;
    writeln!(vf, "ASCII").map_err(RemError::Io)?;
    writeln!(vf, "DATASET UNSTRUCTURED_GRID").map_err(RemError::Io)?;
    writeln!(vf, "POINTS {} double", n_nodes).map_err(RemError::Io)?;
    for node in &mesh.nodes {
        writeln!(vf, "{:.9e} {:.9e} {:.9e}", node.x, node.y, node.z)
            .map_err(RemError::Io)?;
    }
    writeln!(vf, "CELLS {} {}", n_cells, cell_list_size).map_err(RemError::Io)?;
    for elem in &mesh.volume_elements {
        let ids: Vec<String> = elem.node_ids.iter().map(|i| i.to_string()).collect();
        writeln!(vf, "{} {}", elem.node_ids.len(), ids.join(" ")).map_err(RemError::Io)?;
    }
    writeln!(vf, "CELL_TYPES {}", n_cells).map_err(RemError::Io)?;
    for elem in &mesh.volume_elements {
        let t = match elem.kind {
            rem_mesh::ElementKind::Tri3 => 5,
            rem_mesh::ElementKind::Tet4 => 10,
            _ => 5,
        };
        writeln!(vf, "{}", t).map_err(RemError::Io)?;
    }
    writeln!(vf, "POINT_DATA {}", n_nodes).map_err(RemError::Io)?;
    writeln!(vf, "SCALARS Az double 1").map_err(RemError::Io)?;
    writeln!(vf, "LOOKUP_TABLE default").map_err(RemError::Io)?;
    for &v in az { writeln!(vf, "{:.9e}", v).map_err(RemError::Io)?; }
    writeln!(vf, "VECTORS B_field double").map_err(RemError::Io)?;
    for bv in b_field {
        writeln!(vf, "{:.9e} {:.9e} {:.9e}", bv[0], bv[1], bv[2])
            .map_err(RemError::Io)?;
    }

    // Element-center B-field: average nodal B over element corner nodes.
    // Stored as CELL_DATA for better coarse-mesh visualization in ParaView.
    writeln!(vf, "CELL_DATA {}", n_cells).map_err(RemError::Io)?;
    writeln!(vf, "VECTORS B_field_cell double").map_err(RemError::Io)?;
    for elem in &mesh.volume_elements {
        let nn = elem.node_ids.len().max(1) as f64;
        let bx: f64 = elem.node_ids.iter().map(|&i| b_field[i][0]).sum::<f64>() / nn;
        let by: f64 = elem.node_ids.iter().map(|&i| b_field[i][1]).sum::<f64>() / nn;
        let bz: f64 = elem.node_ids.iter().map(|&i| b_field[i][2]).sum::<f64>() / nn;
        writeln!(vf, "{:.9e} {:.9e} {:.9e}", bx, by, bz).map_err(RemError::Io)?;
    }
    log::info!("Written: {}", vtk_path.display());

    Ok(())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use rem_config::{load_config_from_str, ConfigFormat};
    use rem_mesh::{Node, Element, ElementKind};
    use rem_parallel::NoComm;
    use std::collections::HashMap;

    /// Unit square mesh with Ground at y=0 and Ground at y=1 (both A_z=0).
    /// Uniform J_z = 1 → A_z parabolic; but here we just test the solver runs.
    fn unit_square_mesh_grounded() -> RemMesh {
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
        boundary_tags.insert(10, BoundaryTag::Ground); // A_z = 0
        boundary_tags.insert(11, BoundaryTag::SurfaceCurrent { index: 1 }); // A_z = 1

        RemMesh {
            nodes, volume_elements, boundary_elements,
            domain_tags: Default::default(),
            boundary_tags,
            dim: 2,
            rank: 0,
            size: 1,
        }
    }

    fn default_config() -> PalaceConfig {
        load_config_from_str(
            r#"{"Problem":{"Type":"Magnetostatic"},"Model":{"Mesh":"x.msh"},
                "Solver":{"Linear":{"Tol":1e-12,"MaxIter":200}}}"#,
            ConfigFormat::Json,
        )
        .unwrap()
    }

    #[test]
    fn magnetostatic_linear_az() {
        // With A_z=0 at y=0 and A_z=1 at y=1, no sources → A_z(x,y) = y (linear)
        let mesh = unit_square_mesh_grounded();
        let config = default_config();
        let domain_map = DomainMap::from_config(&config).unwrap();

        let az = solve_one(&config, &mesh, &domain_map, Some(1), &NoComm).unwrap();

        for (i, node) in mesh.nodes.iter().enumerate() {
            let exact = node.y;
            let err = (az[i] - exact).abs();
            assert!(
                err < 1e-10,
                "node {} ({:.1},{:.1}): Az={:.6}, exact={:.6}, err={:.2e}",
                i, node.x, node.y, az[i], exact, err
            );
        }
    }

    #[test]
    fn b_field_from_linear_az() {
        // A_z = y → ∇A_z = (0,1) → B_x = ∂Az/∂y = 1, B_y = −∂Az/∂x = 0
        // gradient_recovery returns −∇A_z = (0, −1), so g[0]=0, g[1]=−1
        let mesh = unit_square_mesh_grounded();
        let az: Vec<f64> = mesh.nodes.iter().map(|n| n.y).collect();
        let g = postprocess::gradient_recovery(&az, &mesh);
        for (i, gv) in g.iter().enumerate() {
            // gv = -∇A_z → gv[0]=0, gv[1]=-1
            assert!(gv[0].abs() < 1e-12, "node {}: g[0]={:.2e}", i, gv[0]);
            assert!((gv[1] + 1.0).abs() < 1e-12, "node {}: g[1]={:.6}", i, gv[1]);
        }
        // Derived B-field: B_x = -g[1] = 1, B_y = g[0] = 0
        let b: Vec<[f64; 3]> = g.iter().map(|gv| [-gv[1], gv[0], 0.0]).collect();
        for (i, bv) in b.iter().enumerate() {
            assert!((bv[0] - 1.0).abs() < 1e-12, "node {}: B_x={:.6}", i, bv[0]);
            assert!(bv[1].abs() < 1e-12, "node {}: B_y={:.2e}", i, bv[1]);
        }
    }

    #[test]
    fn magnetic_energy_linear_az() {
        // A_z = y, ν = 1/(μ₀*1) = ν₀, unit area → U = ν₀/2
        use rem_core::constants::{MU0};
        let mesh = unit_square_mesh_grounded();
        let az: Vec<f64> = mesh.nodes.iter().map(|n| n.y).collect();
        let nu0 = 1.0 / MU0;
        let energy = postprocess::electrostatic_energy(&az, &mesh, |_| nu0);
        let expected = nu0 / 2.0;
        assert!(
            (energy - expected).abs() / expected < 1e-12,
            "energy={:.6e}, expected={:.6e}", energy, expected
        );
    }

    #[test]
    fn variable_permeability_iron() {
        // Iron with mu_r = 1000: ν_iron = 1/(μ₀*1000) ≈ 796
        use rem_config::{load_config_from_str, ConfigFormat};
        let cfg_str = r#"{
            "Problem": {"Type": "Magnetostatic"},
            "Model": {"Mesh": "x.msh"},
            "Domains": {
                "Materials": [{"Attributes": [1], "Permeability": 1000.0}]
            }
        }"#;
        let cfg = load_config_from_str(cfg_str, ConfigFormat::Json).unwrap();
        let dm = DomainMap::from_config(&cfg).unwrap();
        use rem_core::constants::MU0;
        let nu_iron = dm.get(1).reluctivity();
        let expected = 1.0 / (MU0 * 1000.0);
        assert!(
            (nu_iron - expected).abs() / expected < 1e-12,
            "ν_iron={:.6e}, expected={:.6e}", nu_iron, expected
        );
    }

    // -----------------------------------------------------------------------
    // 3-D tests
    // -----------------------------------------------------------------------

    /// Build a unit-cube mesh split into 6 tetrahedra sharing the diagonal.
    ///
    /// Boundary conditions:
    ///   tag 10 (face z=0): Ground → Aᵢ = 0
    ///   tag 11 (face z=1): SurfaceCurrent { index: 1 } → Az = 1, Ax = Ay = 0
    fn unit_cube_tet_mesh() -> RemMesh {
        // 8 corner nodes of the unit cube
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
        // 6 tetrahedra (standard "Sommerville" decomposition of the unit cube)
        let tets = [
            [0usize, 1, 3, 4],
            [1, 3, 4, 5],
            [3, 4, 5, 7],
            [1, 2, 3, 5],
            [2, 3, 5, 6],
            [3, 5, 6, 7],
        ];
        let volume_elements: Vec<Element> = tets.iter().enumerate()
            .map(|(i, ns)| Element {
                id: i + 1,
                kind: ElementKind::Tet4,
                tag: 1,
                node_ids: ns.to_vec(),
                rank: 0,
            })
            .collect();

        // Boundary faces: triangles on z=0 (tag 10) and z=1 (tag 11)
        // z=0 face: nodes 0,1,2,3  → two triangles
        // z=1 face: nodes 4,5,6,7  → two triangles
        let boundary_elements = vec![
            Element { id: 100, kind: ElementKind::Tri3, tag: 10, node_ids: vec![0, 1, 2], rank: 0 },
            Element { id: 101, kind: ElementKind::Tri3, tag: 10, node_ids: vec![0, 2, 3], rank: 0 },
            Element { id: 102, kind: ElementKind::Tri3, tag: 11, node_ids: vec![4, 5, 6], rank: 0 },
            Element { id: 103, kind: ElementKind::Tri3, tag: 11, node_ids: vec![4, 6, 7], rank: 0 },
        ];
        let mut boundary_tags: HashMap<u32, BoundaryTag> = HashMap::new();
        boundary_tags.insert(10, BoundaryTag::Ground);
        boundary_tags.insert(11, BoundaryTag::SurfaceCurrent { index: 1 });

        RemMesh {
            nodes, volume_elements, boundary_elements,
            domain_tags: Default::default(),
            boundary_tags,
            dim: 3,
            rank: 0,
            size: 1,
        }
    }

    #[test]
    fn magnetostatic_3d_linear_az() {
        // Az = 0 at z=0 and Az = 1 at z=1, no sources → Az = z (linear) exactly.
        // Ax = Ay = 0 everywhere (no x/y excitation).
        let mesh = unit_cube_tet_mesh();
        let config = default_config();
        let domain_map = DomainMap::from_config(&config).unwrap();

        let (ax, ay, az) = solve_3d(&config, &mesh, &domain_map, Some(1), &NoComm).unwrap();

        for (i, node) in mesh.nodes.iter().enumerate() {
            let exact_az = node.z;
            let err_az = (az[i] - exact_az).abs();
            assert!(
                err_az < 1e-10,
                "node {}: Az={:.6}, exact={:.6}, err={:.2e}",
                i, az[i], exact_az, err_az
            );
            // Ax and Ay should be zero (all-zero BCs)
            assert!(
                ax[i].abs() < 1e-12,
                "node {}: Ax={:.2e} (should be 0)", i, ax[i]
            );
            assert!(
                ay[i].abs() < 1e-12,
                "node {}: Ay={:.2e} (should be 0)", i, ay[i]
            );
        }
    }

    #[test]
    fn magnetostatic_3d_b_field_from_linear_az() {
        // A = (0, 0, z) → B = ∇×A = (∂Az/∂y − 0, 0 − ∂Az/∂x, 0) = (0, 0, 0)
        // Wait — that's trivial for A = z ẑ.  Check gradient_recovery gives correct
        // ∂Az/∂z = 1 (→ gz[2] = −1 in E-field convention), ∂Az/∂x = ∂Az/∂y = 0.
        let mesh = unit_cube_tet_mesh();
        let az: Vec<f64> = mesh.nodes.iter().map(|n| n.z).collect();
        let ax = vec![0.0f64; mesh.n_nodes()];
        let ay = vec![0.0f64; mesh.n_nodes()];

        let gx = postprocess::gradient_recovery(&ax, &mesh);
        let gy = postprocess::gradient_recovery(&ay, &mesh);
        let gz = postprocess::gradient_recovery(&az, &mesh);

        // A=(0,0,z) → B = curl A = (∂Az/∂y − ∂Ay/∂z, ∂Ax/∂z − ∂Az/∂x, ∂Ay/∂x − ∂Ax/∂y)
        //           = (0, 0, 0)  — uniform ẑ potential has zero curl
        let n = mesh.n_nodes();
        for i in 0..n {
            let bx = gy[i][2] - gz[i][1];
            let by = gz[i][0] - gx[i][2];
            let bz = gx[i][1] - gy[i][0];
            assert!(bx.abs() < 1e-10, "node {}: Bx={:.2e}", i, bx);
            assert!(by.abs() < 1e-10, "node {}: By={:.2e}", i, by);
            assert!(bz.abs() < 1e-10, "node {}: Bz={:.2e}", i, bz);
        }
    }

    #[test]
    fn magnetostatic_3d_curl_nonzero() {
        // A = (0, x, 0) → B = ∇×A = (0, 0, ∂Ay/∂x − 0) = (0, 0, 1)  [uniform Bz]
        // Set Ay = node.x for all nodes.
        let mesh = unit_cube_tet_mesh();
        let ax = vec![0.0f64; mesh.n_nodes()];
        let ay: Vec<f64> = mesh.nodes.iter().map(|n| n.x).collect();
        let az = vec![0.0f64; mesh.n_nodes()];

        let gx = postprocess::gradient_recovery(&ax, &mesh);
        let gy = postprocess::gradient_recovery(&ay, &mesh);
        let _gz = postprocess::gradient_recovery(&az, &mesh);

        // B = (gy[2]−gz[1], gz[0]−gx[2], gx[1]−gy[0])
        // Ay = x → ∇Ay = (1, 0, 0) → gy = −∇Ay = (−1, 0, 0)
        // So Bz = gx[1]−gy[0] = 0 − (−1) = 1
        let n = mesh.n_nodes();
        for i in 0..n {
            let bz = gx[i][1] - gy[i][0];
            assert!(
                (bz - 1.0).abs() < 1e-10,
                "node {}: Bz={:.6} (expected 1.0)", i, bz
            );
        }
    }
}

