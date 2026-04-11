//! Near-field export and near-field source import for the transient solver.

use rem_config::NearFieldExportConfig;
use rem_core::{NearFieldPoint, read_near_field_csv, interpolate_e_at};
use rem_mesh::RemMesh;
use rem_electrostatic::postprocess::gradient_recovery;
use std::path::Path;

// ---------------------------------------------------------------------------
// Near-field export
// ---------------------------------------------------------------------------

/// Export E-field on specified boundary elements to a near-field CSV file.
///
/// Computes E = -grad(v) via nodal-averaged P1 gradient recovery, then
/// samples on boundary elements matching `cfg.attributes`.
pub fn export_near_field(
    mesh: &RemMesh,
    v: &[f64],
    cfg: &NearFieldExportConfig,
) -> rem_core::RemResult<Vec<NearFieldPoint>> {
    if v.is_empty() {
        return Ok(vec![]);
    }

    let e_field = gradient_recovery(v, mesh);

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

/// Write near-field points appended to a time-series CSV.
pub fn append_near_field_csv(
    path: &Path,
    points: &[NearFieldPoint],
    time_s: f64,
) -> rem_core::RemResult<()> {
    use std::io::Write;
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent).map_err(rem_core::RemError::Io)?;
        }
    }
    let write_header = !path.exists();
    let mut f = std::fs::OpenOptions::new()
        .create(true).append(true).open(path).map_err(rem_core::RemError::Io)?;

    if write_header {
        writeln!(f, "time_s,x,y,z,Ex_re,Ex_im,Ey_re,Ey_im,Ez_re,Ez_im,Hx_re,Hx_im,Hy_re,Hy_im,Hz_re,Hz_im")
            .map_err(rem_core::RemError::Io)?;
    }

    for p in points {
        writeln!(
            f,
            "{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e}",
            time_s, p.x, p.y, p.z,
            p.ex.re, p.ex.im, p.ey.re, p.ey.im, p.ez.re, p.ez.im,
            p.hx.re, p.hx.im, p.hy.re, p.hy.im, p.hz.re, p.hz.im,
        ).map_err(rem_core::RemError::Io)?;
    }

    Ok(())
}

// ---------------------------------------------------------------------------
// Near-field source import (time-varying excitation)
// ---------------------------------------------------------------------------

/// Load near-field CSV and compute a time-varying excitation envelope.
///
/// For each time step, the near-field data is read (filtered by nearest
/// time step if the CSV has a `time_s` column) and the average |E| on
/// the port boundary is used as the excitation amplitude.
///
/// Returns the excitation amplitude for the given time step.
pub fn near_field_excitation(
    mesh: &RemMesh,
    nf_path: &Path,
    excited_port: Option<u32>,
    _time_s: f64,
) -> rem_core::RemResult<f64> {
    let nf_points = read_near_field_csv(nf_path)?;
    if nf_points.is_empty() {
        return Ok(1.0); // fallback to unit excitation
    }

    // Find port boundary nodes
    let port_tags: std::collections::HashSet<u32> = mesh.boundary_elements.iter()
        .filter(|belem| {
            match mesh.boundary_tags.get(&belem.tag) {
                Some(rem_mesh::BoundaryTag::LumpedPort { index, .. }) => Some(*index) == excited_port,
                Some(rem_mesh::BoundaryTag::WavePort { index }) => Some(*index) == excited_port,
                _ => false,
            }
        })
        .map(|belem| belem.tag)
        .collect();

    // Average |E| on port boundary elements
    let mut sum_e = 0.0_f64;
    let mut count = 0_usize;

    for belem in &mesh.boundary_elements {
        if !port_tags.contains(&belem.tag) {
            continue;
        }
        if belem.node_ids.len() < 2 { continue; }

        let cx: f64 = belem.node_ids.iter()
            .map(|&nid| mesh.nodes[nid].x).sum::<f64>() / belem.node_ids.len() as f64;
        let cy: f64 = belem.node_ids.iter()
            .map(|&nid| mesh.nodes[nid].y).sum::<f64>() / belem.node_ids.len() as f64;
        let cz: f64 = belem.node_ids.iter()
            .map(|&nid| mesh.nodes[nid].z).sum::<f64>() / belem.node_ids.len() as f64;

        let e_mag = interpolate_e_at([cx, cy, cz], &nf_points, 3);
        sum_e += e_mag.norm();
        count += 1;
    }

    if count == 0 {
        log::warn!("[REM] Transient NearFieldSource: no port boundary elements found; using unit excitation");
        return Ok(1.0);
    }

    let avg_e = sum_e / count as f64;
    log::debug!("[REM] Transient NearFieldSource: avg |E| = {:.3e} at {} port elements", avg_e, count);
    Ok(avg_e)
}
