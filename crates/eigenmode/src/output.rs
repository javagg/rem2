//! Output routines for eigenmode results.

use rem_core::{RemError, RemResult};
use rem_mesh::RemMesh;
use std::path::Path;
use std::io::Write;

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
            rem_mesh::ElementKind::Tet10 => 24,
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
