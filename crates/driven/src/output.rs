//! Output routines for driven (frequency-domain) results.

use rem_core::{RemError, RemResult};
use rem_mesh::RemMesh;
use std::path::Path;
use std::io::Write;

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct DomainEnergyRecord {
    pub domain_tag: u32,
    pub material_index: Option<usize>,
    pub energy: f64,
    pub fraction: f64,
}

pub(crate) fn write_wave_port_support_regions(
    out_dir: &str,
    summaries: &[super::port_modal::PortSupportRegionSummary],
) -> RemResult<()> {
    let dir = Path::new(out_dir).join("postpro");
    std::fs::create_dir_all(&dir).map_err(RemError::Io)?;
    let path = dir.join("wave-port-support.csv");
    let mut f = std::fs::File::create(&path).map_err(RemError::Io)?;

    writeln!(
        f,
        "Port Index,Tri3 Elements,Nodes,Boundary Length,X Min,Y Min,X Max,Y Max,Domain Tags"
    )
    .map_err(RemError::Io)?;
    for summary in summaries {
        let tags = summary.domain_tags.iter().map(u32::to_string).collect::<Vec<_>>().join(";");
        writeln!(
            f,
            "{},{},{},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{}",
            summary.port_index,
            summary.n_volume_elements,
            summary.n_nodes,
            summary.boundary_length,
            summary.x_min,
            summary.y_min,
            summary.x_max,
            summary.y_max,
            tags
        )
        .map_err(RemError::Io)?;
    }

    log::info!("Wrote wave-port support regions to {}", path.display());
    Ok(())
}

pub(crate) fn write_peak_domain_energy(
    out_dir: &str,
    peak_freq_hz: f64,
    energies: &[DomainEnergyRecord],
) -> RemResult<()> {
    let dir = Path::new(out_dir).join("postpro");
    std::fs::create_dir_all(&dir).map_err(RemError::Io)?;
    let path = dir.join("domain-E-peak-by-tag.csv");
    let mut f = std::fs::File::create(&path).map_err(RemError::Io)?;

    writeln!(
        f,
        "Peak Frequency (Hz),Domain Tag,Material Index,Electric Field Energy (J),Energy Fraction"
    )
    .map_err(RemError::Io)?;
    for record in energies {
        let material_index = record.material_index.map(|idx| idx.to_string()).unwrap_or_default();
        writeln!(
            f,
            "{:.6e},{},{},{:.6e},{:.6e}",
            peak_freq_hz,
            record.domain_tag,
            material_index,
            record.energy,
            record.fraction
        )
        .map_err(RemError::Io)?;
    }

    log::info!("Wrote peak domain electric energy to {}", path.display());
    Ok(())
}

/// Write S-parameter CSV in Palace format.
///
/// For a single port: `f (Hz), Re(S[1][1]), Im(S[1][1]), |S[1][1]| (dB)`
/// For N ports: columns for every S[i][j] pair.
pub(crate) fn write_s_params(out_dir: &str, results: &[super::FreqResult]) -> RemResult<()> {
    let path = Path::new(out_dir).join("port-S.csv");
    let mut f = std::fs::File::create(&path).map_err(RemError::Io)?;

    // Determine port list from first result
    let port_list = results.first().map(|r| r.port_list.as_slice()).unwrap_or(&[]);
    let n_ports = port_list.len();

    if n_ports > 1 && !results.first().map(|r| r.s_matrix.is_empty()).unwrap_or(true) {
        // Multi-port header
        let mut header = String::from("f (Hz)");
        for &pi in port_list {
            for &pj in port_list {
                header.push_str(&format!(",Re(S[{pi}][{pj}]),Im(S[{pi}][{pj}]),|S[{pi}][{pj}]| (dB)"));
            }
        }
        writeln!(f, "{}", header).map_err(RemError::Io)?;
        for r in results {
            write!(f, "{:.6e}", r.freq_hz).map_err(RemError::Io)?;
            for i in 0..n_ports {
                for j in 0..n_ports {
                    let s = r.s_matrix.get(i).and_then(|row| row.get(j))
                        .copied()
                        .unwrap_or_default();
                    let mag2 = s.norm_sqr();
                    let db = if mag2 > 1e-300 { 10.0 * mag2.log10() } else { -300.0 };
                    write!(f, ",{:.6e},{:.6e},{:.4}", s.re, s.im, db).map_err(RemError::Io)?;
                }
            }
            writeln!(f).map_err(RemError::Io)?;
        }
    } else {
        // Single-port (backward compat) header
        writeln!(f, "f (Hz),Re(S[1][1]),Im(S[1][1]),|S[1][1]| (dB)").map_err(RemError::Io)?;
        for r in results {
            let mag2 = r.s11_re * r.s11_re + r.s11_im * r.s11_im;
            let db = if mag2 > 1e-300 { 10.0 * mag2.log10() } else { -300.0 };
            writeln!(f, "{:.6e},{:.6e},{:.6e},{:.4}", r.freq_hz, r.s11_re, r.s11_im, db)
                .map_err(RemError::Io)?;
        }
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

    Ok(())
}
