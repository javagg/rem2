//! Output routines for driven (frequency-domain) results.

use rem_core::{RemError, RemResult};
use rem_mesh::RemMesh;
use std::path::Path;
use std::io::Write;

/// Write S-parameter CSV in Palace format:
///   f (Hz), Re(S11), Im(S11), |S11| (dB)
pub(crate) fn write_s_params(out_dir: &str, results: &[super::FreqResult]) -> RemResult<()> {
    let path = Path::new(out_dir).join("port-S.csv");
    let mut f = std::fs::File::create(&path).map_err(RemError::Io)?;
    writeln!(f, "f (Hz),Re(S[1][1]),Im(S[1][1]),|S[1][1]| (dB)").map_err(RemError::Io)?;
    for r in results {
        let mag2 = r.s11_re * r.s11_re + r.s11_im * r.s11_im;
        let db = if mag2 > 1e-300 { 10.0 * mag2.log10() } else { -300.0 };
        writeln!(f, "{:.6e},{:.6e},{:.6e},{:.4}", r.freq_hz, r.s11_re, r.s11_im, db)
            .map_err(RemError::Io)?;
    }
    log::info!("Wrote S-parameters to {}", path.display());
    Ok(())
}

/// Write field solution as VTK legacy file.
pub fn write_field_vtk(
    out_dir: &str,
    mesh: &RemMesh,
    phi: &[f64],
    step: usize,
) -> RemResult<()> {
    let path = Path::new(out_dir).join(format!("driven_{:04}.vtk", step));
    let mut f = std::fs::File::create(&path).map_err(RemError::Io)?;

    writeln!(f, "# vtk DataFile Version 2.0").map_err(RemError::Io)?;
    writeln!(f, "REM Driven step {}", step).map_err(RemError::Io)?;
    writeln!(f, "ASCII").map_err(RemError::Io)?;
    writeln!(f, "DATASET UNSTRUCTURED_GRID").map_err(RemError::Io)?;

    let n_nodes = mesh.nodes.len();
    writeln!(f, "POINTS {} double", n_nodes).map_err(RemError::Io)?;
    for node in &mesh.nodes {
        writeln!(f, "{:.6e} {:.6e} {:.6e}", node.x, node.y, node.z).map_err(RemError::Io)?;
    }

    let n_vols = mesh.volume_elements.len();
    let cells_size: usize = mesh.volume_elements.iter().map(|e| 1 + e.node_ids.len()).sum();
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

    Ok(())
}
