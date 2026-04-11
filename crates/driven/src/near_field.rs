//! Near-field export and near-field source import for the driven solver.

use rem_config::NearFieldExportConfig;
use rem_core::{NearFieldPoint, write_near_field_csv, read_near_field_csv, interpolate_e_at};
use rem_mesh::RemMesh;
use rem_electrostatic::postprocess::gradient_recovery;
use std::collections::{HashMap, HashSet};
use std::path::Path;

// ---------------------------------------------------------------------------
// Near-field export
// ---------------------------------------------------------------------------

/// Export E-field on specified boundary elements to a near-field CSV file.
///
/// Computes E = -grad(phi) via nodal-averaged P1 gradient recovery, then
/// samples on boundary elements matching `cfg.attributes`.
pub fn export_near_field(
    mesh: &RemMesh,
    phi: &[f64],
    cfg: &NearFieldExportConfig,
) -> rem_core::RemResult<Vec<NearFieldPoint>> {
    if phi.is_empty() {
        return Ok(vec![]);
    }

    let e_field = gradient_recovery(phi, mesh);

    let mut points = Vec::new();
    let attr_set: std::collections::HashSet<u32> = cfg.attributes.iter().copied().collect();

    for belem in &mesh.boundary_elements {
        if !attr_set.is_empty() && !attr_set.contains(&belem.tag) {
            continue;
        }

        match belem.kind {
            rem_mesh::ElementKind::Tri3 => {
                if belem.node_ids.len() < 3 { continue; }
                let n0 = &mesh.nodes[belem.node_ids[0]];
                let n1 = &mesh.nodes[belem.node_ids[1]];
                let n2 = &mesh.nodes[belem.node_ids[2]];
                let cx = (n0.x + n1.x + n2.x) / 3.0;
                let cy = (n0.y + n1.y + n2.y) / 3.0;
                let cz = (n0.z + n1.z + n2.z) / 3.0;

                let e0 = e_field.get(belem.node_ids[0]).copied().unwrap_or([0.0; 3]);
                let e1 = e_field.get(belem.node_ids[1]).copied().unwrap_or([0.0; 3]);
                let e2 = e_field.get(belem.node_ids[2]).copied().unwrap_or([0.0; 3]);
                let e_avg = [
                    (e0[0] + e1[0] + e2[0]) / 3.0,
                    (e0[1] + e1[1] + e2[1]) / 3.0,
                    (e0[2] + e1[2] + e2[2]) / 3.0,
                ];
                // E = -grad(phi), negate
                let e = [-e_avg[0], -e_avg[1], -e_avg[2]];
                points.push(NearFieldPoint::from_real_e(cx, cy, cz, e));
            }
            rem_mesh::ElementKind::Line2 => {
                if belem.node_ids.len() < 2 { continue; }
                let n0 = &mesh.nodes[belem.node_ids[0]];
                let n1 = &mesh.nodes[belem.node_ids[1]];
                let cx = (n0.x + n1.x) * 0.5;
                let cy = (n0.y + n1.y) * 0.5;
                let cz = (n0.z + n1.z) * 0.5;

                let e0 = e_field.get(belem.node_ids[0]).copied().unwrap_or([0.0; 3]);
                let e1 = e_field.get(belem.node_ids[1]).copied().unwrap_or([0.0; 3]);
                let e_avg = [(e0[0] + e1[0]) * 0.5, (e0[1] + e1[1]) * 0.5, (e0[2] + e1[2]) * 0.5];
                let e = [-e_avg[0], -e_avg[1], -e_avg[2]];
                points.push(NearFieldPoint::from_real_e(cx, cy, cz, e));
            }
            _ => {}
        }
    }

    Ok(points)
}

/// Write near-field points to the configured output path.
pub fn write_near_field(
    output_dir: &Path,
    points: &[NearFieldPoint],
    cfg: &NearFieldExportConfig,
) -> rem_core::RemResult<()> {
    let filename = cfg.output_file.as_deref().unwrap_or("postpro/near_field.csv");
    let path = output_dir.join(filename);
    write_near_field_csv(&path, points)
}

// ---------------------------------------------------------------------------
// Near-field source import
// ---------------------------------------------------------------------------

/// Load near-field data from a CSV file and build Dirichlet values for
/// excited port boundaries.
///
/// For each node on an excited port boundary, the near-field E values at
/// the node position are interpolated from the CSV data.  The magnitude
/// |E| is used as the prescribed potential value (simplified coupling).
///
/// Returns a HashMap<node_index, prescribed_phi_value>.
pub fn build_near_field_dirichlet(
    mesh: &RemMesh,
    nf_path: &Path,
    port_attr_tags: &HashSet<u32>,
) -> rem_core::RemResult<HashMap<usize, f64>> {
    let nf_points = read_near_field_csv(nf_path)?;
    let mut dofs: HashMap<usize, f64> = HashMap::new();

    // Find nodes on port boundaries
    for (_elem_idx, belem) in mesh.boundary_elements.iter().enumerate() {
        if !port_attr_tags.contains(&belem.tag) {
            continue;
        }

        for &nid in &belem.node_ids {
            if dofs.contains_key(&nid) {
                continue;
            }
            let node = &mesh.nodes[nid];
            let pos = [node.x, node.y, node.z];
            let e_mag = interpolate_e_at(pos, &nf_points, 3);
            // Use the real part of |E| as the potential value
            dofs.insert(nid, e_mag.re.abs());
        }
    }

    log::info!("[REM] NearFieldSource: {} port boundary DOFs prescribed from near-field", dofs.len());
    Ok(dofs)
}
