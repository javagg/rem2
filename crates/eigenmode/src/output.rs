//! Output routines for eigenmode results.

use fem_element::nedelec::{TetND1, TetND2, TriND1, TriND2};
use fem_element::reference::VectorReferenceElement;
use fem_mesh::topology::MeshTopology;
use fem_mesh::transformation::ElementTransformation;
use fem_space::{FESpace, HCurlSpace};
use rem_core::{RemError, RemResult};
use rem_mesh::RemMesh;
use std::io::Write;
use std::path::Path;

/// Write eigenfrequencies CSV in Palace format:
///   m, f (Hz)
pub fn write_eigenfrequencies(out_dir: &str, result: &super::EigenResult) -> RemResult<()> {
    let path = Path::new(out_dir).join("eigenfrequencies.csv");
    let mut f = std::fs::File::create(&path).map_err(RemError::Io)?;
    writeln!(f, "m,f (Hz)").map_err(RemError::Io)?;
    for (idx, &freq) in result.frequencies_hz.iter().enumerate() {
        writeln!(f, "{},{:.6e}", idx + 1, freq).map_err(RemError::Io)?;
    }
    log::info!("Wrote eigenfrequencies to {}", path.display());
    Ok(())
}

/// Write a single mode as a VTK legacy file.
pub fn write_mode_vtk(
    out_dir: &str,
    mesh: &RemMesh,
    phi: &[f64],
    mode_idx: usize,
) -> RemResult<()> {
    let path = Path::new(out_dir).join(format!("mode_{}.vtk", mode_idx));
    let mut f = std::fs::File::create(&path).map_err(RemError::Io)?;

    writeln!(f, "# vtk DataFile Version 2.0").map_err(RemError::Io)?;
    writeln!(f, "REM Eigenmode {}", mode_idx).map_err(RemError::Io)?;
    writeln!(f, "ASCII").map_err(RemError::Io)?;
    writeln!(f, "DATASET UNSTRUCTURED_GRID").map_err(RemError::Io)?;

    let n_nodes = mesh.nodes.len();
    writeln!(f, "POINTS {} double", n_nodes).map_err(RemError::Io)?;
    for node in &mesh.nodes {
        writeln!(f, "{:.6e} {:.6e} {:.6e}", node.x, node.y, node.z).map_err(RemError::Io)?;
    }

    let n_vols = mesh.volume_elements.len();
    let cells_size: usize = mesh.volume_elements.iter()
        .map(|e| 1 + e.node_ids.len())
        .sum();
    writeln!(f, "CELLS {} {}", n_vols, cells_size).map_err(RemError::Io)?;
    for elem in &mesh.volume_elements {
        write!(f, "{}", elem.node_ids.len()).map_err(RemError::Io)?;
        for &nid in &elem.node_ids {
            write!(f, " {}", nid).map_err(RemError::Io)?;
        }
        writeln!(f).map_err(RemError::Io)?;
    }

    writeln!(f, "CELL_TYPES {}", n_vols).map_err(RemError::Io)?;
    for elem in &mesh.volume_elements {
        let vtk_type = match elem.kind {
            rem_mesh::ElementKind::Tri3  => 5,
            rem_mesh::ElementKind::Tet4  => 10,
            rem_mesh::ElementKind::Hex8  => 12,
            rem_mesh::ElementKind::Quad4 => 9,
            rem_mesh::ElementKind::Tri6  => 22,
            rem_mesh::ElementKind::Tri10 => 22,  // VTK_QUADRATIC_TRIANGLE (approx)
            rem_mesh::ElementKind::Tet10 => 24,
            rem_mesh::ElementKind::Tet20 => 24,  // VTK_QUADRATIC_TETRA (approx)
            rem_mesh::ElementKind::Line2 => 3,
            rem_mesh::ElementKind::Line3 => 21,  // VTK_QUADRATIC_EDGE
        };
        writeln!(f, "{}", vtk_type).map_err(RemError::Io)?;
    }

    writeln!(f, "POINT_DATA {}", n_nodes).map_err(RemError::Io)?;
    writeln!(f, "SCALARS phi double 1").map_err(RemError::Io)?;
    writeln!(f, "LOOKUP_TABLE default").map_err(RemError::Io)?;
    for i in 0..n_nodes {
        let v = if i < phi.len() { phi[i] } else { 0.0 };
        writeln!(f, "{:.6e}", v).map_err(RemError::Io)?;
    }

    log::info!("Wrote mode {} VTK to {}", mode_idx, path.display());
    Ok(())
}

/// Write a single HCurl mode as a VTK file with vector E-field cell data.
///
/// Evaluates the Nedelec edge-element solution E = Σ x_e · φ_e at the centroid
/// of each element and writes `CELL_DATA VECTORS E_field`.
#[cfg(not(target_arch = "wasm32"))]
pub fn write_mode_vector_vtk(
    out_dir: &str,
    mesh: &RemMesh,
    phi: &[f64],
    mode_idx: usize,
    order: u8,
) -> RemResult<()> {
    let path = Path::new(out_dir).join(format!("mode_{}.vtk", mode_idx));
    let mut f = std::fs::File::create(&path).map_err(RemError::Io)?;

    writeln!(f, "# vtk DataFile Version 2.0").map_err(RemError::Io)?;
    writeln!(f, "REM HCurl Eigenmode {}", mode_idx).map_err(RemError::Io)?;
    writeln!(f, "ASCII").map_err(RemError::Io)?;
    writeln!(f, "DATASET UNSTRUCTURED_GRID").map_err(RemError::Io)?;

    match mesh.dim {
        2 => write_vector_vtk_2d(&mut f, mesh, phi, order)?,
        3 => write_vector_vtk_3d(&mut f, mesh, phi, order)?,
        d => return Err(RemError::Config(format!(
            "HCurl VTK output requires 2-D or 3-D mesh, got dim={}", d
        ))),
    }

    log::info!("Wrote HCurl mode {} VTK to {}", mode_idx, path.display());
    Ok(())
}

/// 2-D helper: TriND1/TriND2 field recovery + SimplexMesh VTK.
#[cfg(not(target_arch = "wasm32"))]
fn write_vector_vtk_2d(
    f: &mut std::fs::File,
    mesh: &RemMesh,
    phi: &[f64],
    order: u8,
) -> RemResult<()> {
    use fem_mesh::SimplexMesh;
    let simplex: SimplexMesh<2> = mesh.to_simplex_mesh_2d();
    let space = HCurlSpace::new(simplex, order);
    let smesh = space.mesh();

    let dim = 2usize;
    let ref_elem: Box<dyn VectorReferenceElement> = match order {
        1 => Box::new(TriND1),
        2 => Box::new(TriND2),
        o => return Err(RemError::Config(format!(
            "HCurl VTK supports order 1/2 only, got {o}"
        ))),
    };
    let n_ldofs = ref_elem.n_dofs();
    let xi: Vec<f64> = vec![1.0 / 3.0; 2]; // triangle centroid

    let mut ref_phi = vec![0.0; n_ldofs * dim];
    let mut phys_phi = vec![0.0; n_ldofs * dim];
    let n_elem = smesh.n_elements();
    let n_node = smesh.n_nodes();

    // ── Write POINTS ────────────────────────────────────────────────────────
    writeln!(f, "POINTS {} double", n_node).map_err(RemError::Io)?;
    for n in 0..n_node as u32 {
        let c = smesh.coords_of(n);
        writeln!(f, "{:.6e} {:.6e} 0.0", c[0], c[1]).map_err(RemError::Io)?;
    }

    // ── Write CELLS + compute E-field ───────────────────────────────────────
    let cells_size: usize = smesh.elem_iter()
        .map(|e| 1 + smesh.elem_nodes(e).len())
        .sum();
    writeln!(f, "CELLS {} {}", n_elem, cells_size).map_err(RemError::Io)?;
    let mut e_fields = Vec::with_capacity(n_elem);
    for e in smesh.elem_iter() {
        let nodes = smesh.elem_nodes(e);
        write!(f, "{}", nodes.len()).map_err(RemError::Io)?;
        for &nid in nodes {
            write!(f, " {}", nid).map_err(RemError::Io)?;
        }
        writeln!(f).map_err(RemError::Io)?;

        // Field recovery at centroid
        let elem_dofs = space.element_dofs(e);
        let signs = space.element_signs(e);
        let tr = ElementTransformation::from_simplex_nodes(smesh, nodes);
        let j_inv_t = tr.jacobian_inv_t().clone();

        ref_elem.eval_basis_vec(&xi, &mut ref_phi);
        for i in 0..n_ldofs {
            for r in 0..dim {
                let mut s = 0.0;
                for c in 0..dim {
                    s += j_inv_t[(r, c)] * ref_phi[i * dim + c];
                }
                phys_phi[i * dim + r] = s;
            }
        }
        for i in 0..n_ldofs {
            for c in 0..dim {
                phys_phi[i * dim + c] *= signs[i];
            }
        }

        let mut e = [0.0f64; 3];
        for i in 0..n_ldofs {
            let c = phi[elem_dofs[i] as usize];
            e[0] += c * phys_phi[i * dim];
            e[1] += c * phys_phi[i * dim + 1];
        }
        e_fields.push(e);
    }

    // ── Write CELL_TYPES ────────────────────────────────────────────────────
    writeln!(f, "CELL_TYPES {}", n_elem).map_err(RemError::Io)?;
    for e in smesh.elem_iter() {
        match smesh.element_type(e) {
            fem_mesh::ElementType::Tri6 => writeln!(f, "22")?,
            _ => writeln!(f, "5")?, // Tri3 or fallback
        }
    }

    // ── Write CELL_DATA VECTORS ─────────────────────────────────────────────
    writeln!(f, "CELL_DATA {}", n_elem).map_err(RemError::Io)?;
    writeln!(f, "VECTORS E_field double").map_err(RemError::Io)?;
    for e in &e_fields {
        writeln!(f, "{:.6e} {:.6e} {:.6e}", e[0], e[1], e[2]).map_err(RemError::Io)?;
    }

    Ok(())
}

/// 3-D helper: TetND1/TetND2 field recovery + SimplexMesh VTK.
#[cfg(not(target_arch = "wasm32"))]
fn write_vector_vtk_3d(
    f: &mut std::fs::File,
    mesh: &RemMesh,
    phi: &[f64],
    order: u8,
) -> RemResult<()> {
    use fem_mesh::SimplexMesh;
    let simplex: SimplexMesh<3> = mesh.to_simplex_mesh();
    let space = HCurlSpace::new(simplex, order);
    let smesh = space.mesh();

    let dim = 3usize;
    let ref_elem: Box<dyn VectorReferenceElement> = match order {
        1 => Box::new(TetND1),
        2 => Box::new(TetND2),
        o => return Err(RemError::Config(format!(
            "HCurl VTK supports order 1/2 only, got {o}"
        ))),
    };
    let n_ldofs = ref_elem.n_dofs();
    let xi: Vec<f64> = vec![1.0 / 4.0; 3]; // tet centroid

    let mut ref_phi = vec![0.0; n_ldofs * dim];
    let mut phys_phi = vec![0.0; n_ldofs * dim];
    let n_elem = smesh.n_elements();
    let n_node = smesh.n_nodes();

    // ── Write POINTS ────────────────────────────────────────────────────────
    writeln!(f, "POINTS {} double", n_node).map_err(RemError::Io)?;
    for n in 0..n_node as u32 {
        let c = smesh.coords_of(n);
        writeln!(f, "{:.6e} {:.6e} {:.6e}", c[0], c[1], c[2]).map_err(RemError::Io)?;
    }

    // ── Write CELLS + compute E-field ───────────────────────────────────────
    let cells_size: usize = smesh.elem_iter()
        .map(|e| 1 + smesh.elem_nodes(e).len())
        .sum();
    writeln!(f, "CELLS {} {}", n_elem, cells_size).map_err(RemError::Io)?;
    let mut e_fields = Vec::with_capacity(n_elem);
    for e in smesh.elem_iter() {
        let nodes = smesh.elem_nodes(e);
        write!(f, "{}", nodes.len()).map_err(RemError::Io)?;
        for &nid in nodes {
            write!(f, " {}", nid).map_err(RemError::Io)?;
        }
        writeln!(f).map_err(RemError::Io)?;

        // Field recovery at centroid
        let elem_dofs = space.element_dofs(e);
        let signs = space.element_signs(e);
        let tr = ElementTransformation::from_simplex_nodes(smesh, nodes);
        let j_inv_t = tr.jacobian_inv_t().clone();

        ref_elem.eval_basis_vec(&xi, &mut ref_phi);
        for i in 0..n_ldofs {
            for r in 0..dim {
                let mut s = 0.0;
                for c in 0..dim {
                    s += j_inv_t[(r, c)] * ref_phi[i * dim + c];
                }
                phys_phi[i * dim + r] = s;
            }
        }
        for i in 0..n_ldofs {
            for c in 0..dim {
                phys_phi[i * dim + c] *= signs[i];
            }
        }

        let mut e = [0.0f64; 3];
        for i in 0..n_ldofs {
            let c = phi[elem_dofs[i] as usize];
            for d in 0..dim {
                e[d] += c * phys_phi[i * dim + d];
            }
        }
        e_fields.push(e);
    }

    // ── Write CELL_TYPES ────────────────────────────────────────────────────
    writeln!(f, "CELL_TYPES {}", n_elem).map_err(RemError::Io)?;
    for e in smesh.elem_iter() {
        match smesh.element_type(e) {
            fem_mesh::ElementType::Tet10 => writeln!(f, "24")?,
            _ => writeln!(f, "10")?, // Tet4 or fallback
        }
    }

    // ── Write CELL_DATA VECTORS ─────────────────────────────────────────────
    writeln!(f, "CELL_DATA {}", n_elem).map_err(RemError::Io)?;
    writeln!(f, "VECTORS E_field double").map_err(RemError::Io)?;
    for e in &e_fields {
        writeln!(f, "{:.6e} {:.6e} {:.6e}", e[0], e[1], e[2]).map_err(RemError::Io)?;
    }

    Ok(())
}
