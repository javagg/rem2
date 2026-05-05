//! rem-bem — Laplace/Helmholtz Boundary Element Method solver
//!
//! Solves exterior Laplace and Helmholtz problems using the boundary integral
//! equation approach with P0 (constant) and P1 (linear) basis functions.
//!
//! # Formulation
//!
//! For the exterior Laplace problem (electrostatics):
//! ```text
//! ½ φ(r) + ∫_S ∂G/∂n'(r,r') φ(r') dS' = ∫_S G(r,r') σ(r') dS'
//! ```
//! where G = 1/(4πR) is the Laplace Green function,
//! σ = ∂φ/∂n is the normal flux (surface charge / ε₀).
//!
//! # Architecture
//! ```text
//! SurfaceMesh (from rem-surface)
//!     ↓
//! assemble_laplace_bem (V + K matrices)
//!     ↓
//! solve (dense LU)
//!     ↓
//! postprocess (capacitance, potential)
//! ```

pub mod kernel;
pub mod assemble;
pub mod solve;
pub mod postprocess;

use rem_config::PalaceConfig;
use rem_core::RemResult;
use rem_surface::quadrature::TriQuad;
use rem_parallel::NoComm;

/// Run the BEM solver from a Palace config.
///
/// Loads the mesh specified in `config.model.mesh`, extracts PEC conductor
/// surfaces from `config.boundaries.pec`, assembles the Laplace BIE system,
/// solves the exterior Neumann problem, and writes capacitance results.
pub fn run(config: &PalaceConfig) -> RemResult<()> {
    log::info!("\n=== Boundary Element Method (BEM) solver ===\n");

    let output_dir = std::path::Path::new(config.problem.output_dir());
    #[cfg(not(target_arch = "wasm32"))]
    std::fs::create_dir_all(output_dir)?;

    // Load mesh
    let mesh = rem_mesh::load_mesh(config, &NoComm)?;
    log::info!("Mesh loaded: {} nodes, {} volume elements, {} boundary elements",
        mesh.nodes.len(), mesh.volume_elements.len(), mesh.boundary_elements.len());

    // Get conductor surface attributes from boundaries
    let pec_attrs: Vec<u32> = config.boundaries.pec
        .as_ref()
        .map(|p| p.attributes.clone())
        .unwrap_or_default();

    if pec_attrs.is_empty() {
        return Err(rem_core::RemError::Config(
            "BEM solver requires at least one PEC boundary (Boundaries.PEC.Attributes)".to_string()
        ));
    }

    // Extract surface mesh
    let surf = rem_surface::surface_mesh::SurfaceMesh::extract(&mesh, &pec_attrs)?;
    log::info!("Surface mesh: {} triangular faces", surf.faces.len());

    // Quadrature rules
    let quad = TriQuad::new(3);
    let n_duffy = 4;

    // Assemble Laplace BEM matrices (P0 basis)
    log::info!("Assembling BEM matrices...");
    let (v_mat, k_mat) = assemble::assemble_laplace_p0(&surf, &quad, n_duffy)?;
    log::info!("  V: {}×{}, K: {}×{}", v_mat.nrows(), v_mat.ncols(), k_mat.nrows(), k_mat.ncols());

    // Apply known normal flux σ = 1.0 on all conductor surfaces (unit excitation)
    let sigma: Vec<f64> = vec![1.0_f64; surf.faces.len()];

    // Solve Neumann problem: K φ = V σ
    log::info!("Solving Neumann problem...");
    let phi = assemble::solve_neumann(&v_mat, &k_mat, &sigma)?;

    // Post-process: compute capacitance
    let q = postprocess::total_charge(&sigma, &surf);
    let v0 = phi.iter().copied().fold(0.0_f64, f64::max);
    let c = postprocess::capacitance(&sigma, &surf, v0);

    log::info!("Total charge Q = {:.6e} C", q);
    log::info!("Max surface potential φ_max = {:.6e} V", v0);
    log::info!("Capacitance C = {:.6e} F ({:.3e} pF)", c, c * 1e12);

    // Write results CSV
    let csv_path = output_dir.join("bem-capacitance.csv");
    #[cfg(not(target_arch = "wasm32"))]
    {
        let mut wtr = csv::Writer::from_path(&csv_path)
            .map_err(|e| rem_core::RemError::Io(e.into()))?;
        wtr.write_record(&["quantity", "value", "unit"])
            .map_err(|e| rem_core::RemError::Io(e.into()))?;
        wtr.write_record(&["total_charge", &format!("{:.12e}", q), "C"])
            .map_err(|e| rem_core::RemError::Io(e.into()))?;
        wtr.write_record(&["max_potential", &format!("{:.12e}", v0), "V"])
            .map_err(|e| rem_core::RemError::Io(e.into()))?;
        wtr.write_record(&["capacitance", &format!("{:.12e}", c), "F"])
            .map_err(|e| rem_core::RemError::Io(e.into()))?;
        wtr.flush().map_err(|e| rem_core::RemError::Io(e.into()))?;
    }
    log::info!("Results written to {}", csv_path.display());

    log::info!("\nBEM solve completed.\n");
    Ok(())
}
