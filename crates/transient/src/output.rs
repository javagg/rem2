//! Output routines for the transient solver.

use rem_core::{RemError, RemResult};
use rem_mesh::RemMesh;
use std::io::Write;
use std::path::Path;

/// Write a VTK file with the scalar field `v` at the given step index.
pub fn write_field_vtk(
    out_dir: &str,
    mesh: &RemMesh,
    v: &[f64],
    step: usize,
) -> RemResult<()> {
    let vtk_dir = Path::new(out_dir).join("paraview");
    std::fs::create_dir_all(&vtk_dir).map_err(RemError::Io)?;
    let path = vtk_dir.join(format!("transient_{:04}.vtk", step));
    let mut f = std::fs::File::create(&path).map_err(RemError::Io)?;

    let n_nodes = mesh.n_nodes();
    let n_cells = mesh.n_volume_elements();
    let cell_list_size: usize = mesh.volume_elements.iter().map(|e| e.node_ids.len() + 1).sum();

    writeln!(f, "# vtk DataFile Version 3.0").map_err(RemError::Io)?;
    writeln!(f, "rem transient step {}", step).map_err(RemError::Io)?;
    writeln!(f, "ASCII").map_err(RemError::Io)?;
    writeln!(f, "DATASET UNSTRUCTURED_GRID").map_err(RemError::Io)?;
    writeln!(f, "POINTS {} double", n_nodes).map_err(RemError::Io)?;
    for node in &mesh.nodes {
        writeln!(f, "{:.9e} {:.9e} {:.9e}", node.x, node.y, node.z).map_err(RemError::Io)?;
    }
    writeln!(f, "CELLS {} {}", n_cells, cell_list_size).map_err(RemError::Io)?;
    for elem in &mesh.volume_elements {
        let ids: Vec<String> = elem.node_ids.iter().map(|i| i.to_string()).collect();
        writeln!(f, "{} {}", elem.node_ids.len(), ids.join(" ")).map_err(RemError::Io)?;
    }
    writeln!(f, "CELL_TYPES {}", n_cells).map_err(RemError::Io)?;
    for elem in &mesh.volume_elements {
        let t = match elem.kind {
            rem_mesh::ElementKind::Tri3 => 5,
            rem_mesh::ElementKind::Tet4 => 10,
            _ => 5,
        };
        writeln!(f, "{}", t).map_err(RemError::Io)?;
    }
    writeln!(f, "POINT_DATA {}", n_nodes).map_err(RemError::Io)?;
    writeln!(f, "SCALARS phi double 1").map_err(RemError::Io)?;
    writeln!(f, "LOOKUP_TABLE default").map_err(RemError::Io)?;
    for &val in v {
        writeln!(f, "{:.9e}", val).map_err(RemError::Io)?;
    }

    log::debug!("Written: {}", path.display());
    Ok(())
}

/// Write a CSV time series of port voltage.
pub fn write_time_series(
    out_dir: &str,
    times: &[f64],
    voltages: &[f64],
) -> RemResult<()> {
    let dir = Path::new(out_dir).join("postpro");
    std::fs::create_dir_all(&dir).map_err(RemError::Io)?;
    let path = dir.join("port-t.csv");
    let mut f = std::fs::File::create(&path).map_err(RemError::Io)?;

    writeln!(f, r#""Time (s)","Port Voltage (V)""#).map_err(RemError::Io)?;
    for (t, v) in times.iter().zip(voltages.iter()) {
        writeln!(f, "{:.9e},{:.9e}", t, v).map_err(RemError::Io)?;
    }

    log::info!("Written: {}", path.display());
    Ok(())
}
