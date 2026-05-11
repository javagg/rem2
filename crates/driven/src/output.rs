//! Output routines for driven (frequency-domain) results.

use fem_element::nedelec::{TetND1, TetND2, TriND1, TriND2};
use fem_element::reference::VectorReferenceElement;
use fem_mesh::topology::MeshTopology;
use fem_mesh::transformation::ElementTransformation;
use fem_space::{FESpace, HCurlSpace};
use num_complex::Complex64;
use rem_core::{RemError, RemResult};
use rem_mesh::RemMesh;
use std::io::Write;
use std::path::Path;
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


/// Write Palace-compatible `postpro/port-VI.csv` with complex voltage, current,
/// and time-average power at each frequency.
///
/// Palace format (one row per frequency, columns repeated for each port):
/// ```text
/// f (Hz), Re(V[1]), Im(V[1]), Re(I[1]), Im(I[1]), Re(P[1]), Im(P[1]), ...
/// ```
pub(crate) fn write_port_vi_csv(
    out_dir: &str,
    results: &[super::FreqResult],
) -> RemResult<()> {
    // Collect the superset of port indices across all results
    let mut port_set: Vec<u32> = results
        .iter()
        .flat_map(|r| r.port_vi.iter().map(|v| v.port_index))
        .collect::<std::collections::BTreeSet<_>>()
        .into_iter()
        .collect();
    if port_set.is_empty() {
        return Ok(()); // nothing to write
    }
    port_set.sort();

    let dir = Path::new(out_dir).join("postpro");
    std::fs::create_dir_all(&dir).map_err(RemError::Io)?;
    let path = dir.join("port-VI.csv");
    let mut f = std::fs::File::create(&path).map_err(RemError::Io)?;

    // Header
    let mut header = String::from("f (Hz)");
    for &p in &port_set {
        header.push_str(&format!(
            ",Re(V[{p}]),Im(V[{p}]),Re(I[{p}]),Im(I[{p}]),Re(P[{p}]),Im(P[{p}])"
        ));
    }
    writeln!(f, "{}", header).map_err(RemError::Io)?;

    for r in results {
        write!(f, "{:.6e}", r.freq_hz).map_err(RemError::Io)?;
        for &p in &port_set {
            if let Some(vi) = r.port_vi.iter().find(|v| v.port_index == p) {
                write!(
                    f,
                    ",{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e}",
                    vi.v.re, vi.v.im,
                    vi.i.re, vi.i.im,
                    vi.p.re, vi.p.im
                )
                .map_err(RemError::Io)?;
            } else {
                write!(f, ",0,0,0,0,0,0").map_err(RemError::Io)?;
            }
        }
        writeln!(f).map_err(RemError::Io)?;
    }

    log::info!("Wrote port voltage/current/power to {}", path.display());
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

    Ok(())
}

/// Write HCurl driven field as VTK with vector E-field cell data (real + imag).
///
/// Evaluates the Nedelec edge-element solution E = Σ x_e · φ_e at each
/// element centroid and writes `CELL_DATA VECTORS E_real` / `E_imag`.
#[cfg(not(target_arch = "wasm32"))]
pub fn write_field_vector_vtk(
    out_dir: &str,
    mesh: &RemMesh,
    e_dofs: &[Complex64],
    step: usize,
    order: u8,
) -> RemResult<()> {
    let path = Path::new(out_dir).join(format!("driven_{:04}.vtk", step));
    let mut f = std::fs::File::create(&path).map_err(RemError::Io)?;

    writeln!(f, "# vtk DataFile Version 2.0").map_err(RemError::Io)?;
    writeln!(f, "REM HCurl Driven step {}", step).map_err(RemError::Io)?;
    writeln!(f, "ASCII").map_err(RemError::Io)?;
    writeln!(f, "DATASET UNSTRUCTURED_GRID").map_err(RemError::Io)?;

    match mesh.dim {
        2 => write_driven_vtk_2d(&mut f, mesh, e_dofs, order)?,
        3 => write_driven_vtk_3d(&mut f, mesh, e_dofs, order)?,
        d => return Err(RemError::Config(format!(
            "HCurl driven VTK requires 2-D or 3-D, got dim={}", d
        ))),
    }
    Ok(())
}

/// 2-D driven VTK helper (TriND1/TriND2).
#[cfg(not(target_arch = "wasm32"))]
fn write_driven_vtk_2d(
    f: &mut std::fs::File,
    mesh: &RemMesh,
    e_dofs: &[Complex64],
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
        o => return Err(RemError::Config(format!("Driven VTK supports order 1/2 only, got {o}"))),
    };
    let n_ldofs = ref_elem.n_dofs();
    let xi: Vec<f64> = vec![1.0 / 3.0; 2];
    let mut ref_phi = vec![0.0; n_ldofs * dim];
    let mut phys_phi = vec![0.0; n_ldofs * dim];
    let n_elem = smesh.n_elements();
    let n_node = smesh.n_nodes();

    writeln!(f, "POINTS {} double", n_node).map_err(RemError::Io)?;
    for n in 0..n_node as u32 {
        let c = smesh.coords_of(n);
        writeln!(f, "{:.6e} {:.6e} 0.0", c[0], c[1]).map_err(RemError::Io)?;
    }

    let cells_size: usize = smesh.elem_iter().map(|e| 1 + smesh.elem_nodes(e).len()).sum();
    writeln!(f, "CELLS {} {}", n_elem, cells_size).map_err(RemError::Io)?;
    let mut e_real_fields = Vec::with_capacity(n_elem);
    let mut e_imag_fields = Vec::with_capacity(n_elem);
    for e in smesh.elem_iter() {
        let nodes = smesh.elem_nodes(e);
        write!(f, "{}", nodes.len()).map_err(RemError::Io)?;
        for &nid in nodes { write!(f, " {}", nid).map_err(RemError::Io)?; }
        writeln!(f).map_err(RemError::Io)?;

        let elem_dofs = space.element_dofs(e);
        let signs = space.element_signs(e);
        let tr = ElementTransformation::from_simplex_nodes(smesh, nodes);
        let j_inv_t = tr.jacobian_inv_t().clone();
        ref_elem.eval_basis_vec(&xi, &mut ref_phi);
        for i in 0..n_ldofs {
            for r in 0..dim {
                let mut s = 0.0;
                for c in 0..dim { s += j_inv_t[(r, c)] * ref_phi[i * dim + c]; }
                phys_phi[i * dim + r] = s;
            }
        }
        for i in 0..n_ldofs {
            for c in 0..dim { phys_phi[i * dim + c] *= signs[i]; }
        }

        let mut er = [0.0f64; 3];
        let mut ei = [0.0f64; 3];
        for i in 0..n_ldofs {
            let c = e_dofs[elem_dofs[i] as usize];
            let base_re = c.re * phys_phi[i * dim];
            let base_im = c.im * phys_phi[i * dim];
            er[0] += base_re; er[1] += c.re * phys_phi[i * dim + 1];
            ei[0] += base_im; ei[1] += c.im * phys_phi[i * dim + 1];
        }
        e_real_fields.push(er);
        e_imag_fields.push(ei);
    }

    writeln!(f, "CELL_TYPES {}", n_elem).map_err(RemError::Io)?;
    for e in smesh.elem_iter() {
        match smesh.element_type(e) {
            fem_mesh::ElementType::Tri6 => writeln!(f, "22")?,
            _ => writeln!(f, "5")?,
        }
    }

    writeln!(f, "CELL_DATA {}", n_elem).map_err(RemError::Io)?;
    writeln!(f, "VECTORS E_real double").map_err(RemError::Io)?;
    for e in &e_real_fields { writeln!(f, "{:.6e} {:.6e} {:.6e}", e[0], e[1], e[2])?; }
    writeln!(f, "VECTORS E_imag double").map_err(RemError::Io)?;
    for e in &e_imag_fields { writeln!(f, "{:.6e} {:.6e} {:.6e}", e[0], e[1], e[2])?; }

    Ok(())
}

/// 3-D driven VTK helper (TetND1/TetND2).
#[cfg(not(target_arch = "wasm32"))]
fn write_driven_vtk_3d(
    f: &mut std::fs::File,
    mesh: &RemMesh,
    e_dofs: &[Complex64],
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
        o => return Err(RemError::Config(format!("Driven VTK supports order 1/2 only, got {o}"))),
    };
    let n_ldofs = ref_elem.n_dofs();
    let xi: Vec<f64> = vec![1.0 / 4.0; 3];
    let mut ref_phi = vec![0.0; n_ldofs * dim];
    let mut phys_phi = vec![0.0; n_ldofs * dim];
    let n_elem = smesh.n_elements();
    let n_node = smesh.n_nodes();

    writeln!(f, "POINTS {} double", n_node).map_err(RemError::Io)?;
    for n in 0..n_node as u32 {
        let c = smesh.coords_of(n);
        writeln!(f, "{:.6e} {:.6e} {:.6e}", c[0], c[1], c[2]).map_err(RemError::Io)?;
    }

    let cells_size: usize = smesh.elem_iter().map(|e| 1 + smesh.elem_nodes(e).len()).sum();
    writeln!(f, "CELLS {} {}", n_elem, cells_size).map_err(RemError::Io)?;
    let mut e_real_fields = Vec::with_capacity(n_elem);
    let mut e_imag_fields = Vec::with_capacity(n_elem);
    for e in smesh.elem_iter() {
        let nodes = smesh.elem_nodes(e);
        write!(f, "{}", nodes.len()).map_err(RemError::Io)?;
        for &nid in nodes { write!(f, " {}", nid).map_err(RemError::Io)?; }
        writeln!(f).map_err(RemError::Io)?;

        let elem_dofs = space.element_dofs(e);
        let signs = space.element_signs(e);
        let tr = ElementTransformation::from_simplex_nodes(smesh, nodes);
        let j_inv_t = tr.jacobian_inv_t().clone();
        ref_elem.eval_basis_vec(&xi, &mut ref_phi);
        for i in 0..n_ldofs {
            for r in 0..dim {
                let mut s = 0.0;
                for c in 0..dim { s += j_inv_t[(r, c)] * ref_phi[i * dim + c]; }
                phys_phi[i * dim + r] = s;
            }
        }
        for i in 0..n_ldofs {
            for c in 0..dim { phys_phi[i * dim + c] *= signs[i]; }
        }

        let mut er = [0.0f64; 3];
        let mut ei = [0.0f64; 3];
        for i in 0..n_ldofs {
            let c = e_dofs[elem_dofs[i] as usize];
            for d in 0..dim {
                er[d] += c.re * phys_phi[i * dim + d];
                ei[d] += c.im * phys_phi[i * dim + d];
            }
        }
        e_real_fields.push(er);
        e_imag_fields.push(ei);
    }

    writeln!(f, "CELL_TYPES {}", n_elem).map_err(RemError::Io)?;
    for e in smesh.elem_iter() {
        match smesh.element_type(e) {
            fem_mesh::ElementType::Tet10 => writeln!(f, "24")?,
            _ => writeln!(f, "10")?,
        }
    }

    writeln!(f, "CELL_DATA {}", n_elem).map_err(RemError::Io)?;
    writeln!(f, "VECTORS E_real double").map_err(RemError::Io)?;
    for e in &e_real_fields { writeln!(f, "{:.6e} {:.6e} {:.6e}", e[0], e[1], e[2])?; }
    writeln!(f, "VECTORS E_imag double").map_err(RemError::Io)?;
    for e in &e_imag_fields { writeln!(f, "{:.6e} {:.6e} {:.6e}", e[0], e[1], e[2])?; }

    Ok(())
}
