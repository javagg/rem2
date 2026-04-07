//! Magnetostatic solver — Phase 4.
//!
//! Solves the 2-D vector potential formulation:
//!
//!   −∇·(ν ∇A_z) = J_z      (in Ω)
//!   A_z = 0                  (on Dirichlet boundaries — Ground)
//!
//! where ν = 1/(μ₀ μᵣ) is the magnetic reluctivity.
//!
//! Post-processing:
//!   B_x =  ∂A_z/∂y ,   B_y = −∂A_z/∂x   (from gradient of A_z)
//!   Inductance = 2 U_mag / I²  where U_mag = (1/2) ∫ ν |∇A_z|² dΩ
//!
//! The assembly is identical to the electrostatic solver — only the
//! coefficient (ν vs ε) and the interpretation of the solution differ.

use rem_config::PalaceConfig;
use rem_core::{RemResult, solve_pcg};
use rem_parallel::Comm;
use rem_materials::DomainMap;
use rem_mesh::{RemMesh, BoundaryTag};
use rem_mesh::gmsh::read_msh_file;
use rem_electrostatic::assemble;
use rem_electrostatic::bc;
use rem_electrostatic::postprocess;
use std::path::Path;

/// Entry point called from rem-cli.
pub fn run(config: &PalaceConfig, comm: &dyn Comm) -> RemResult<()> {
    log::info!("=== Magnetostatic solver (2-D A_z) ===");

    if config.solver.order > 1 {
        log::warn!(
            "Solver.Order={} requested but only P1 (order=1) is implemented; \
             higher-order assembly is pending. Running P1.",
            config.solver.order
        );
    }
    if config.model.refinement.max_iter > 0 {
        log::warn!(
            "Model.Refinement.MaxIter={} requested but AMR loop is not yet integrated \
             into the solver; running a single solve pass.",
            config.model.refinement.max_iter
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

    let domain_map = DomainMap::from_config(config)?;
    let output_dir = Path::new(config.problem.output_dir());

    // Find surface current excitation (analogous to LumpedPort in electrostatics)
    let excited_tag = find_surface_current_tag(&mesh); // returns SurfaceCurrent INDEX

    let az = solve_one(config, &mesh, &domain_map, excited_tag, comm)?;

    // gradient_recovery returns −∇A_z (same sign convention as E = −∇φ).
    // So g = (−∂Az/∂x, −∂Az/∂y).
    // B_x = ∂Az/∂y = −g[1],   B_y = −∂Az/∂x = g[0]
    let grad_az = postprocess::gradient_recovery(&az, &mesh);
    let b_field: Vec<[f64; 3]> = grad_az.iter()
        .map(|g| [-g[1], g[0], 0.0])
        .collect();

    // Magnetic energy U = (1/2) ∫ ν |∇A_z|² dΩ
    let nu_fn = |tag: u32| domain_map.get(tag).reluctivity();
    let energy = postprocess::electrostatic_energy(&az, &mesh, nu_fn);
    log::info!("Magnetic energy: {:.6e} J/m", energy);

    // Write CSV + VTK
    write_outputs(output_dir, &mesh, &az, &b_field, energy)?;

    Ok(())
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

    // Assemble stiffness with reluctivity coefficient
    let nu_fn = |tag: u32| domain_map.get(tag).reluctivity();
    let triplet = assemble::assemble_stiffness(mesh, nu_fn)?;
    let mut mat = triplet.to_csr();
    let mut rhs = vec![0.0f64; n];

    // Dirichlet BCs: Ground → A_z = 0, SurfaceCurrent excited → A_z = 1
    let dofs = collect_magnetostatic_dofs(mesh, excitation_tag, 1.0);
    log::info!("Dirichlet DOFs: {}", dofs.len());
    bc::apply_dirichlet(&mut mat, &mut rhs, &dofs);

    // Solve
    let lin = &config.solver.linear;
    let result = solve_pcg(&mat, &rhs, lin.tol, lin.max_iter, comm);
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
}
