/// Palace-compatible CSV output for electrostatic results.
///
/// Reference: Palace postpro/domain-E.csv and terminal-C.csv formats.

use rem_config::{BoundaryPostprocessingSpec, EnergyPostSpec};
use rem_core::{RemError, RemResult};
use std::io::Write;
use std::path::Path;

/// Write the electrostatic energy to `<output_dir>/postpro/domain-E.csv`.
///
/// Palace format:
/// ```text
/// "Frequency (GHz)","Electric Field Energy (J)","Magnetic Field Energy (J)","Total Energy (J)"
/// 0.000000e0,1.234567e-12,0.000000e0,1.234567e-12
/// ```
pub fn write_domain_energy(output_dir: &Path, e_energy: f64) -> RemResult<()> {
    let dir = output_dir.join("postpro");
    std::fs::create_dir_all(&dir).map_err(RemError::Io)?;
    let path = dir.join("domain-E.csv");

    let mut f = std::fs::File::create(&path).map_err(RemError::Io)?;
    writeln!(
        f,
        r#""Frequency (GHz)","Electric Field Energy (J)","Magnetic Field Energy (J)","Total Energy (J)""#
    )
    .map_err(RemError::Io)?;
    writeln!(
        f,
        "0.000000e0,{:.6e},0.000000e0,{:.6e}",
        e_energy, e_energy
    )
    .map_err(RemError::Io)?;

    log::info!("Written: {}", path.display());
    Ok(())
}

/// Write per-domain electrostatic diagnostics to `<output_dir>/postpro/domain-E-by-tag.csv`.
pub fn write_domain_energy_by_tag<T>(output_dir: &Path, energies: &[T]) -> RemResult<()>
where
    T: DomainEnergyRow,
{
    let dir = output_dir.join("postpro");
    std::fs::create_dir_all(&dir).map_err(RemError::Io)?;
    let path = dir.join("domain-E-by-tag.csv");

    let mut f = std::fs::File::create(&path).map_err(RemError::Io)?;
    writeln!(
        f,
        r#""Domain Tag","Material Index","Electric Field Energy (J)","Energy Fraction""#
    ).map_err(RemError::Io)?;
    for energy in energies {
        let material_index = energy
            .material_index()
            .map(|idx| idx.to_string())
            .unwrap_or_default();
        writeln!(
            f,
            "{},{},{:.6e},{:.6e}",
            energy.domain_tag(),
            material_index,
            energy.energy(),
            energy.fraction()
        ).map_err(RemError::Io)?;
    }

    log::info!("Written: {}", path.display());
    Ok(())
}

pub trait DomainEnergyRow {
    fn domain_tag(&self) -> u32;
    fn material_index(&self) -> Option<usize>;
    fn energy(&self) -> f64;
    fn fraction(&self) -> f64;
}

/// Write per-group electrostatic energy to `<output_dir>/postpro/energy-E.csv`.
///
/// Each row corresponds to one `EnergyPostSpec`, summing energies over all
/// domain tags listed in `spec.attributes`.  Matches Palace's grouping output.
///
/// `per_tag_energies` is a slice of `(domain_tag, energy_J)` pairs.
pub fn write_energy_groups_csv(
    output_dir: &Path,
    specs: &[EnergyPostSpec],
    per_tag_energies: &[(u32, f64)],
) -> RemResult<()> {
    let dir = output_dir.join("postpro");
    std::fs::create_dir_all(&dir).map_err(RemError::Io)?;
    let path = dir.join("energy-E.csv");

    let mut f = std::fs::File::create(&path).map_err(RemError::Io)?;
    writeln!(f, r#""Index","Attributes","Electric Field Energy (J)""#)
        .map_err(RemError::Io)?;

    for spec in specs {
        let group_energy: f64 = per_tag_energies
            .iter()
            .filter(|(tag, _)| spec.attributes.contains(tag))
            .map(|(_, e)| e)
            .sum();
        let attrs_str = spec.attributes.iter()
            .map(|a| a.to_string())
            .collect::<Vec<_>>()
            .join(";");
        writeln!(f, "{},\"{}\",{:.6e}", spec.index, attrs_str, group_energy)
            .map_err(RemError::Io)?;
    }

    log::info!("Written: {}", path.display());
    Ok(())
}

/// Write magnetic energy groups to `<output_dir>/postpro/energy-B.csv`.
///
/// Mirrors `write_energy_groups_csv` but uses the "Magnetic Field Energy (J)" label
/// and writes to `energy-B.csv` (static / magnetostatic solver output).
pub fn write_energy_groups_csv_magnetic(
    output_dir: &Path,
    specs: &[EnergyPostSpec],
    per_tag_energies: &[(u32, f64)],
) -> RemResult<()> {
    let dir = output_dir.join("postpro");
    std::fs::create_dir_all(&dir).map_err(RemError::Io)?;
    let path = dir.join("energy-B.csv");

    let mut f = std::fs::File::create(&path).map_err(RemError::Io)?;
    writeln!(f, r#""Index","Attributes","Magnetic Field Energy (J)""#)
        .map_err(RemError::Io)?;

    for spec in specs {
        let group_energy: f64 = per_tag_energies
            .iter()
            .filter(|(tag, _)| spec.attributes.contains(tag))
            .map(|(_, e)| e)
            .sum();
        let attrs_str = spec.attributes.iter()
            .map(|a| a.to_string())
            .collect::<Vec<_>>()
            .join(";");
        writeln!(f, "{},\"{}\",{:.6e}", spec.index, attrs_str, group_energy)
            .map_err(RemError::Io)?;
    }

    log::info!("Written: {}", path.display());
    Ok(())
}

/// Write surface flux results to `<output_dir>/postpro/surface-flux.csv`.
///
/// Palace format:
/// ```text
/// "Frequency (GHz)","Flux[1] (C)","Flux[2] (C)",...
/// 0.000000e0,1.234e-12,...
/// ```
/// `specs` contains the `BoundaryPostprocessingSpec` entries (for headers/index).
/// `fluxes` is the precomputed flux value for each spec entry in order.
/// `unit` is the unit label, e.g. `"C"` (electric) or `"Wb"` (magnetic).
pub fn write_surface_flux_csv(
    output_dir: &Path,
    specs: &[&BoundaryPostprocessingSpec],
    fluxes: &[f64],
    unit: &str,
) -> RemResult<()> {
    if specs.is_empty() { return Ok(()); }
    let dir = output_dir.join("postpro");
    std::fs::create_dir_all(&dir).map_err(RemError::Io)?;
    let path = dir.join("surface-flux.csv");

    let mut f = std::fs::File::create(&path).map_err(RemError::Io)?;

    // Header
    let cols: Vec<String> = specs.iter()
        .map(|s| format!("\"Flux[{}] ({})\"", s.index, unit))
        .collect();
    writeln!(f, "\"Frequency (GHz)\",{}", cols.join(",")).map_err(RemError::Io)?;

    // Data row (static, frequency = 0 for static solvers)
    let vals: Vec<String> = fluxes.iter().map(|v| format!("{:.6e}", v)).collect();
    writeln!(f, "0.000000e0,{}", vals.join(",")).map_err(RemError::Io)?;

    log::info!("Written: {}", path.display());
    Ok(())
}

/// Write the capacitance matrix to `<output_dir>/postpro/terminal-C.csv`.
///
/// Palace format (N ports):
/// ```text
/// "Frequency (GHz)","C[1][1] (F)","C[1][2] (F)",...
/// 0.000000e0,1.234e-12,...
/// ```
pub fn write_capacitance_matrix(
    output_dir: &Path,
    c_matrix: &[Vec<f64>], // c_matrix[i][j] = C_ij
) -> RemResult<()> {
    let n = c_matrix.len();
    if n == 0 { return Ok(()); }

    let dir = output_dir.join("postpro");
    std::fs::create_dir_all(&dir).map_err(RemError::Io)?;
    let path = dir.join("terminal-C.csv");

    let mut f = std::fs::File::create(&path).map_err(RemError::Io)?;

    // Header
    let header_cols: Vec<String> = (0..n)
        .flat_map(|i| (0..n).map(move |j| format!("\"C[{}][{}] (F)\"", i + 1, j + 1)))
        .collect();
    writeln!(f, "\"Frequency (GHz)\",{}", header_cols.join(","))
        .map_err(RemError::Io)?;

    // Data row
    let vals: Vec<String> = c_matrix
        .iter()
        .flat_map(|row| row.iter().map(|&v| format!("{:.6e}", v)))
        .collect();
    writeln!(f, "0.000000e0,{}", vals.join(",")).map_err(RemError::Io)?;

    log::info!("Written: {}", path.display());
    Ok(())
}

/// Write VTU (legacy ASCII VTK) for scalar potential and vector E-field.
///
/// The output is a minimal but valid VTK UnstructuredGrid that ParaView can open.
pub fn write_vtk(
    output_dir: &Path,
    mesh: &rem_mesh::RemMesh,
    phi: &[f64],
    e_field: &[[f64; 3]],
) -> RemResult<()> {
    let dir = output_dir.join("paraview");
    std::fs::create_dir_all(&dir).map_err(RemError::Io)?;
    let path = dir.join("solution.vtk");
    let mut f = std::fs::File::create(&path).map_err(RemError::Io)?;

    let n_nodes = mesh.n_nodes();
    let n_cells = mesh.n_volume_elements();

    // --- Header ---
    writeln!(f, "# vtk DataFile Version 3.0").map_err(RemError::Io)?;
    writeln!(f, "rem electrostatic solution").map_err(RemError::Io)?;
    writeln!(f, "ASCII").map_err(RemError::Io)?;
    writeln!(f, "DATASET UNSTRUCTURED_GRID").map_err(RemError::Io)?;

    // --- Points ---
    writeln!(f, "POINTS {} double", n_nodes).map_err(RemError::Io)?;
    for node in &mesh.nodes {
        writeln!(f, "{:.9e} {:.9e} {:.9e}", node.x, node.y, node.z)
            .map_err(RemError::Io)?;
    }

    // --- Cells ---
    // Count total integers needed
    let cell_list_size: usize = mesh.volume_elements
        .iter()
        .map(|e| e.node_ids.len() + 1)
        .sum();
    writeln!(f, "CELLS {} {}", n_cells, cell_list_size).map_err(RemError::Io)?;
    for elem in &mesh.volume_elements {
        let ids_str: Vec<String> = elem.node_ids.iter().map(|i| i.to_string()).collect();
        writeln!(f, "{} {}", elem.node_ids.len(), ids_str.join(" "))
            .map_err(RemError::Io)?;
    }

    // --- Cell types ---
    writeln!(f, "CELL_TYPES {}", n_cells).map_err(RemError::Io)?;
    for elem in &mesh.volume_elements {
        let vtk_type = match elem.kind {
            rem_mesh::ElementKind::Tri3  => 5,
            rem_mesh::ElementKind::Quad4 => 9,
            rem_mesh::ElementKind::Tet4  => 10,
            rem_mesh::ElementKind::Hex8  => 12,
            _ => 5,
        };
        writeln!(f, "{}", vtk_type).map_err(RemError::Io)?;
    }

    // --- Point data: φ and E ---
    writeln!(f, "POINT_DATA {}", n_nodes).map_err(RemError::Io)?;

    writeln!(f, "SCALARS phi double 1").map_err(RemError::Io)?;
    writeln!(f, "LOOKUP_TABLE default").map_err(RemError::Io)?;
    for &v in phi {
        writeln!(f, "{:.9e}", v).map_err(RemError::Io)?;
    }

    writeln!(f, "VECTORS E_field double").map_err(RemError::Io)?;
    for ev in e_field {
        writeln!(f, "{:.9e} {:.9e} {:.9e}", ev[0], ev[1], ev[2])
            .map_err(RemError::Io)?;
    }

    log::info!("Written: {}", path.display());
    Ok(())
}
