//! Output routines for SBR+: VTK surface current and RCS CSV.
//!
//! Follows the same file layout as `rem_mom::postprocess` so results can be
//! compared directly.

use rem_core::RemResult;
use rem_mom::surface_mesh::SurfaceMesh;
use std::io::Write as IoWrite;
use std::path::Path;

use crate::po_integral::{CurrentMap, rcs_pattern, rcs_pattern_with_ptd};
use crate::ptd::BoundaryEdge;
use crate::excitation::PlaneWave;
use num_complex::Complex64;

// ---------------------------------------------------------------------------
// VTK surface current output
// ---------------------------------------------------------------------------

/// Write surface current distribution to VTK Legacy ASCII format.
///
/// Each triangular face carries scalars: `J_mag`, `J_x_re`, `J_y_re`, `J_z_re`.
/// Compatible with ParaView and VisIt.
pub fn write_surface_vtk(
    path: &Path,
    currents: &CurrentMap,
    surf: &SurfaceMesh,
) -> RemResult<()> {
    let mut f = std::fs::File::create(path)?;

    let n_pts   = surf.nodes.len();
    let n_cells = surf.faces.len();

    writeln!(f, "# vtk DataFile Version 3.0")?;
    writeln!(f, "SBR+ surface current")?;
    writeln!(f, "ASCII")?;
    writeln!(f, "DATASET UNSTRUCTURED_GRID")?;
    writeln!(f)?;

    writeln!(f, "POINTS {} float", n_pts)?;
    for &[x, y, z] in &surf.nodes {
        writeln!(f, "{:.8e} {:.8e} {:.8e}", x, y, z)?;
    }
    writeln!(f)?;

    writeln!(f, "CELLS {} {}", n_cells, n_cells * 4)?;
    for face in &surf.faces {
        writeln!(f, "3 {} {} {}", face.nodes[0], face.nodes[1], face.nodes[2])?;
    }
    writeln!(f)?;

    writeln!(f, "CELL_TYPES {}", n_cells)?;
    for _ in 0..n_cells {
        writeln!(f, "5")?; // VTK_TRIANGLE
    }
    writeln!(f)?;

    writeln!(f, "CELL_DATA {}", n_cells)?;

    // |J| magnitude
    writeln!(f, "SCALARS J_mag float 1")?;
    writeln!(f, "LOOKUP_TABLE default")?;
    for fc in currents {
        let jmag = (fc.j[0].norm_sqr() + fc.j[1].norm_sqr() + fc.j[2].norm_sqr()).sqrt();
        writeln!(f, "{:.8e}", jmag)?;
    }

    // J vector components (real parts)
    for (label, idx) in [("J_x_re", 0usize), ("J_y_re", 1), ("J_z_re", 2)] {
        writeln!(f, "SCALARS {} float 1", label)?;
        writeln!(f, "LOOKUP_TABLE default")?;
        for fc in currents {
            writeln!(f, "{:.8e}", fc.j[idx].re)?;
        }
    }

    Ok(())
}

// ---------------------------------------------------------------------------
// RCS CSV output
// ---------------------------------------------------------------------------

/// Append RCS results to `<output_dir>/postpro/rcs_sbr.csv`.
///
/// Header is written on first creation. Format matches `rem_mom::postprocess`.
/// If `ptd_edges` and `wave` are provided, PTD fringe correction is applied.
pub fn write_rcs(
    output_dir: &Path,
    freq: f64,
    currents: &CurrentMap,
    surf: &SurfaceMesh,
    k: f64,
    theta_deg: &[f64],
    phi_deg: &[f64],
) -> RemResult<()> {
    write_rcs_inner(output_dir, freq, currents, surf, k, theta_deg, phi_deg, None, None, None)
}

/// Append RCS results with PTD correction.
pub fn write_rcs_with_ptd(
    output_dir: &Path,
    freq: f64,
    currents: &CurrentMap,
    surf: &SurfaceMesh,
    k: f64,
    theta_deg: &[f64],
    phi_deg: &[f64],
    wave: &PlaneWave,
    ptd_edges: &[BoundaryEdge],
    e_inc_at: &dyn Fn(&[f64; 3]) -> [Complex64; 3],
) -> RemResult<()> {
    write_rcs_inner(output_dir, freq, currents, surf, k, theta_deg, phi_deg,
                    Some(wave), Some(ptd_edges), Some(e_inc_at))
}

fn write_rcs_inner(
    output_dir: &Path,
    freq: f64,
    currents: &CurrentMap,
    surf: &SurfaceMesh,
    k: f64,
    theta_deg: &[f64],
    phi_deg: &[f64],
    wave: Option<&PlaneWave>,
    ptd_edges: Option<&[BoundaryEdge]>,
    e_inc_at: Option<&dyn Fn(&[f64; 3]) -> [Complex64; 3]>,
) -> RemResult<()> {
    let path = output_dir.join("postpro").join("rcs_sbr.csv");
    let write_header = !path.exists();

    let mut file = std::fs::OpenOptions::new()
        .create(true).append(true).open(&path)?;

    if write_header {
        writeln!(file, "Freq (GHz),Theta (deg),Phi (deg),RCS (dBsm)")?;
    }

    let freq_ghz = freq / 1.0e9;

    let pattern = match (wave, ptd_edges, e_inc_at) {
        (Some(w), Some(e), Some(f)) => {
            rcs_pattern_with_ptd(currents, surf, k, theta_deg, phi_deg, w, e, f)
        }
        _ => rcs_pattern(currents, surf, k, theta_deg, phi_deg),
    };

    for (i_th, &th) in theta_deg.iter().enumerate() {
        for (i_ph, &ph) in phi_deg.iter().enumerate() {
            let sigma_m2 = pattern[i_th][i_ph];
            let rcs_dbsm = if sigma_m2 > 1e-40 {
                10.0 * sigma_m2.log10()
            } else {
                -999.9
            };
            writeln!(file, "{:.6e},{:.1},{:.1},{:.4}", freq_ghz, th, ph, rcs_dbsm)?;
        }
    }

    Ok(())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use rem_mom::surface_mesh::{SurfaceMesh, TriFace};
    use rem_mom::surface_mesh::tri_geometry;
    use crate::po_integral::zero_currents;

    fn flat_surf() -> SurfaceMesh {
        let nodes = vec![[0.0f64,0.0,0.0],[1.0,0.0,0.0],[0.0,1.0,0.0]];
        let (c,n,a) = tri_geometry(&nodes[0],&nodes[1],&nodes[2]);
        SurfaceMesh {
            nodes, faces: vec![TriFace{nodes:[0,1,2],centroid:c,normal:n,area:a}],
            edges: vec![], boundary_edges: vec![], face_attrs: vec![0], global_node_ids: vec![],
        }
    }

    #[test]
    fn vtk_writes_and_has_header() {
        let surf = flat_surf();
        let cur = zero_currents(&surf);
        let tmp = std::env::temp_dir().join("sbr_test_current.vtk");
        write_surface_vtk(&tmp, &cur, &surf).expect("VTK write failed");
        let content = std::fs::read_to_string(&tmp).unwrap();
        assert!(content.contains("SBR+ surface current"));
        assert!(content.contains("J_mag"));
        let _ = std::fs::remove_file(&tmp);
    }

    #[test]
    fn rcs_csv_writes() {
        let surf = flat_surf();
        let cur = zero_currents(&surf);
        let tmp_dir = std::env::temp_dir().join("sbr_rcs_test");
        std::fs::create_dir_all(tmp_dir.join("postpro")).ok();
        let path = tmp_dir.join("postpro").join("rcs_sbr.csv");
        let _ = std::fs::remove_file(&path);
        write_rcs(&tmp_dir, 1e9, &cur, &surf, 20.0, &[0.0, 90.0], &[0.0]).expect("CSV write failed");
        assert!(path.exists());
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("Freq (GHz)"));
        let _ = std::fs::remove_dir_all(&tmp_dir);
    }
}
