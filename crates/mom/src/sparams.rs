//! S-parameter matrix computation for MoM port-excited problems.
//!
//! Given N MomLumpedPort definitions and a pre-assembled Z matrix,
//! runs one solve per port and extracts the N×N S-matrix.

use crate::port::MomLumpedPort;
use crate::surface_mesh::SurfaceMesh;
use crate::basis::rwg::RwgBasis;
use crate::assemble::lu_solve;
use nalgebra::DMatrix;
use num_complex::Complex64;
use rem_core::RemResult;
use rem_touchstone::{write_snp, TsFormat, TsFreqUnit};
use std::path::Path;

/// N×N S-matrix at a single frequency.
#[derive(Debug, Clone)]
pub struct SMatrix {
    /// Number of ports.
    pub n_ports: usize,
    /// Frequency [Hz].
    pub freq_hz: f64,
    /// Row-major S-matrix: data[i * n_ports + j] = S_{i+1, j+1}.
    pub data: Vec<Complex64>,
}

impl SMatrix {
    /// S_{row+1, col+1} (0-based indices).
    pub fn get(&self, row: usize, col: usize) -> Complex64 {
        self.data[row * self.n_ports + col]
    }
}

/// Compute the N×N S-matrix for a set of lumped ports.
///
/// For each excitation port `p`:
/// 1. Build RHS from `port_p.excitation_rhs(v0=1V)`.
/// 2. Solve Z·I = V_rhs → current coefficients.
/// 3. For each observation port `q`: extract I_q = port_q.extract_current(coeffs).
/// 4. S_{qp} = V_q - Z0_q * I_q  (simplified wave-port formula, V_p_fwd = 1V)
///
/// S_{pp} = 1 - Z0_p * I_p  and  S_{qp} = -Z0_q * I_q for q≠p.
pub fn compute_s_matrix(
    surf: &SurfaceMesh,
    bases: &[RwgBasis],
    ports: &[MomLumpedPort],
    z_mat: &DMatrix<Complex64>,
    freq_hz: f64,
) -> RemResult<SMatrix> {
    let n_ports = ports.len();
    let n_rwg   = bases.len();
    let v0      = Complex64::new(1.0, 0.0);

    // Solve one system per excitation port
    let mut all_currents: Vec<Vec<Complex64>> = Vec::with_capacity(n_ports);
    for port_p in ports {
        let rhs = port_p.excitation_rhs(surf, bases, n_rwg, v0);
        let coeffs = lu_solve(z_mat, &rhs)?;
        all_currents.push(coeffs);
    }

    // Build N×N S-matrix
    let mut data = vec![Complex64::ZERO; n_ports * n_ports];
    for (p, currents_p) in all_currents.iter().enumerate() {
        for (q, port_q) in ports.iter().enumerate() {
            let z0_q = port_q.z0;
            let i_q  = port_q.extract_current(surf, bases, currents_p);
            // S_{qp} = V_q - Z0_q * I_q  (V_p_fwd = 1V)
            let v_q = if q == p { v0 } else { Complex64::ZERO };
            let s_qp = v_q - Complex64::new(z0_q, 0.0) * i_q;
            data[q * n_ports + p] = s_qp;
        }
    }

    Ok(SMatrix { n_ports, freq_hz, data })
}

/// Run a full S-parameter sweep over multiple frequencies.
///
/// Returns one `SMatrix` per frequency.
pub fn s_param_sweep(
    surf: &SurfaceMesh,
    bases: &[RwgBasis],
    ports: &[MomLumpedPort],
    freq_hz_list: &[f64],
    build_z: &dyn Fn(f64) -> RemResult<DMatrix<Complex64>>,
) -> RemResult<Vec<SMatrix>> {
    freq_hz_list.iter().map(|&f| {
        let z = build_z(f)?;
        compute_s_matrix(surf, bases, ports, &z, f)
    }).collect()
}

/// Write all frequency-sweep S-matrices to a Touchstone `.s{N}p` file.
pub fn write_touchstone(matrices: &[SMatrix], path: &Path, z0: f64) -> RemResult<()> {
    use std::io::Write;
    if matrices.is_empty() { return Ok(()); }
    let n_ports = matrices[0].n_ports;

    let freqs: Vec<f64>            = matrices.iter().map(|m| m.freq_hz).collect();
    let s_data: Vec<Vec<Complex64>> = matrices.iter().map(|m| m.data.clone()).collect();

    let content = write_snp(&freqs, &s_data, n_ports, z0, TsFormat::Ri, TsFreqUnit::Ghz);
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }
    let mut f = std::fs::File::create(path)?;
    f.write_all(content.as_bytes())?;
    Ok(())
}

/// Append S-parameter data to a Palace-compatible `port-S.csv`.
pub fn append_palace_csv(matrices: &[SMatrix], path: &Path) -> RemResult<()> {
    use std::io::Write;
    let write_header = !path.exists();
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }
    let mut f = std::fs::OpenOptions::new().create(true).append(true).open(path)?;
    if write_header {
        let mut hdr = "Freq (GHz)".to_string();
        if let Some(m) = matrices.first() {
            for i in 0..m.n_ports {
                for j in 0..m.n_ports {
                    hdr.push_str(&format!(",Re(S{}{}),Im(S{}{}),|S{}{}| (dB)",
                        i+1, j+1, i+1, j+1, i+1, j+1));
                }
            }
        }
        writeln!(f, "{hdr}")?;
    }
    for m in matrices {
        let mut line = format!("{:.9e}", m.freq_hz / 1e9);
        for &s in &m.data {
            let db = if s.norm() > 1e-300 { 20.0 * s.norm().log10() } else { -999.0 };
            line.push_str(&format!(",{:.8e},{:.8e},{:.4}", s.re, s.im, db));
        }
        writeln!(f, "{line}")?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::surface_mesh::{SurfaceMesh, TriFace, SharedEdge, tri_geometry};
    use crate::basis::rwg::generate_rwg_bases;
    use crate::port::MomLumpedPort;

    fn two_tri_surf_with_attrs(attr: u32) -> SurfaceMesh {
        let nodes = vec![
            [0.0_f64, 0.0, 0.0],
            [1.0,     0.0, 0.0],
            [0.5,     1.0, 0.0],
            [-0.5,    1.0, 0.0],
        ];
        let (c0, n0, a0) = tri_geometry(&nodes[0], &nodes[1], &nodes[2]);
        let (c1, n1, a1) = tri_geometry(&nodes[0], &nodes[2], &nodes[3]);
        let faces = vec![
            TriFace { nodes:[0,1,2], centroid:c0, normal:n0, area:a0 },
            TriFace { nodes:[0,2,3], centroid:c1, normal:n1, area:a1 },
        ];
        let edges = vec![SharedEdge {
            nodes: [0,2], plus_face: 0, minus_face: 1,
            length: (0.5_f64.powi(2) + 1.0_f64.powi(2)).sqrt(),
        }];
        SurfaceMesh {
            nodes, faces, edges,
            boundary_edges: vec![[0,1],[1,2],[2,3],[3,0]],
            face_attrs: vec![attr, attr],
        }
    }

    /// Build a trivial 1×1 Z-matrix (identity) for testing.
    fn identity_z(n: usize) -> DMatrix<Complex64> {
        let mut z = DMatrix::<Complex64>::zeros(n, n);
        for i in 0..n { z[(i,i)] = Complex64::new(1.0, 0.0); }
        z
    }

    #[test]
    fn s_matrix_shape_single_port() {
        let surf = two_tri_surf_with_attrs(1);
        let bases = generate_rwg_bases(&surf);
        let port = MomLumpedPort::from_surface(&surf, &bases, &[1], 1, "x", 50.0).unwrap();
        let ports = vec![port];
        let z = identity_z(bases.len());
        let sm = compute_s_matrix(&surf, &bases, &ports, &z, 1e9).unwrap();
        assert_eq!(sm.n_ports, 1);
        assert_eq!(sm.data.len(), 1);
        assert!(sm.data[0].re.is_finite());
    }

    #[test]
    fn write_touchstone_creates_file() {
        let matrices = vec![
            SMatrix { n_ports: 1, freq_hz: 1e9, data: vec![Complex64::new(0.5, -0.3)] },
            SMatrix { n_ports: 1, freq_hz: 2e9, data: vec![Complex64::new(0.4, -0.2)] },
        ];
        let tmp = std::env::temp_dir().join("test_mom_s1p.s1p");
        write_touchstone(&matrices, &tmp, 50.0).expect("write failed");
        let content = std::fs::read_to_string(&tmp).unwrap();
        assert!(content.contains("# GHz S RI"), "option line missing");
        let data_lines: usize = content.lines()
            .filter(|l| !l.starts_with('!') && !l.starts_with('#') && !l.trim().is_empty())
            .count();
        assert_eq!(data_lines, 2, "expected 2 data lines");
        let _ = std::fs::remove_file(&tmp);
    }

    #[test]
    fn append_palace_csv_creates_file() {
        let matrices = vec![
            SMatrix { n_ports: 1, freq_hz: 1e9, data: vec![Complex64::new(0.5, -0.3)] },
        ];
        let tmp = std::env::temp_dir().join("test_mom_port_s.csv");
        let _ = std::fs::remove_file(&tmp);
        append_palace_csv(&matrices, &tmp).expect("csv write failed");
        let content = std::fs::read_to_string(&tmp).unwrap();
        assert!(content.contains("Freq"), "CSV header missing");
        let _ = std::fs::remove_file(&tmp);
    }
}
