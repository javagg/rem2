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
use rem_core::{RemResult, C0};
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

/// Variant of `compute_s_matrix` that also returns the current coefficient
/// vectors (one `Vec<Complex64>` per excitation port).  The current for
/// excitation port `p` is at index `p` in the returned `Vec`.
pub fn compute_s_matrix_with_currents(
    surf: &SurfaceMesh,
    bases: &[RwgBasis],
    ports: &[MomLumpedPort],
    z_mat: &DMatrix<Complex64>,
    freq_hz: f64,
) -> RemResult<(SMatrix, Vec<Vec<Complex64>>)> {
    let n_ports = ports.len();
    let n_rwg   = bases.len();
    let v0      = Complex64::new(1.0, 0.0);

    let mut all_currents: Vec<Vec<Complex64>> = Vec::with_capacity(n_ports);
    for port_p in ports {
        let rhs = port_p.excitation_rhs(surf, bases, n_rwg, v0);
        let coeffs = lu_solve(z_mat, &rhs)?;
        all_currents.push(coeffs);
    }

    let mut data = vec![Complex64::ZERO; n_ports * n_ports];
    for (p, currents_p) in all_currents.iter().enumerate() {
        for (q, port_q) in ports.iter().enumerate() {
            let z0_q = port_q.z0;
            let i_q  = port_q.extract_current(surf, bases, currents_p);
            let v_q  = if q == p { v0 } else { Complex64::ZERO };
            let s_qp = v_q - Complex64::new(z0_q, 0.0) * i_q;
            data[q * n_ports + p] = s_qp;
        }
    }

    let sm = SMatrix { n_ports, freq_hz, data };
    Ok((sm, all_currents))
}

// ── Z / Y matrix conversions ───────────────────────────────────────────────

/// Convert an S-parameter matrix to a Z-parameter matrix.
///
/// Uses the relation Z = Z₀ · (I + S)(I - S)⁻¹  (uniform reference impedance).
pub fn s_to_z(s: &SMatrix, z0: f64) -> RemResult<SMatrix> {
    let n = s.n_ports;
    let z0c = Complex64::new(z0, 0.0);
    let i_mat = DMatrix::<Complex64>::identity(n, n);
    let s_mat = DMatrix::from_row_slice(n, n, &s.data);
    let ip_s  = &i_mat + &s_mat;
    let im_s  = &i_mat - &s_mat;
    let z_mat = ip_s * im_s.try_inverse()
        .ok_or_else(|| rem_core::RemError::Config("S→Z: (I−S) is singular".into()))?
        * z0c;
    let data: Vec<Complex64> = (0..n)
        .flat_map(|r| (0..n).map(|c| z_mat[(r, c)]).collect::<Vec<_>>())
        .collect();
    Ok(SMatrix { n_ports: n, freq_hz: s.freq_hz, data })
}

/// Convert a Z-parameter matrix to a Y-parameter matrix (Y = Z⁻¹).
pub fn z_to_y(z: &SMatrix) -> RemResult<SMatrix> {
    let n = z.n_ports;
    let z_mat = DMatrix::from_row_slice(n, n, &z.data);
    let y_mat = z_mat.try_inverse()
        .ok_or_else(|| rem_core::RemError::Config("Z→Y: Z is singular".into()))?;
    let data: Vec<Complex64> = (0..n)
        .flat_map(|r| (0..n).map(|c| y_mat[(r, c)]).collect::<Vec<_>>())
        .collect();
    Ok(SMatrix { n_ports: n, freq_hz: z.freq_hz, data })
}

/// Write a sequence of parameter matrices (Z or Y) to a Palace-compatible CSV.
///
/// `param_name` should be `"Z"` or `"Y"`.  Column format mirrors `port-S.csv`.
pub fn write_param_csv(matrices: &[SMatrix], path: &Path, param_name: &str) -> RemResult<()> {
    use std::io::Write;
    if matrices.is_empty() { return Ok(()); }
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() { std::fs::create_dir_all(parent)?; }
    }
    let write_header = !path.exists();
    let mut f = std::fs::OpenOptions::new().create(true).append(true).open(path)?;
    if write_header {
        let mut hdr = "Freq (GHz)".to_string();
        if let Some(m) = matrices.first() {
            for i in 0..m.n_ports {
                for j in 0..m.n_ports {
                    hdr.push_str(&format!(",Re({pn}{i}{j}),Im({pn}{i}{j})",
                        pn = param_name, i = i+1, j = j+1));
                }
            }
        }
        writeln!(f, "{hdr}")?;
    }
    for m in matrices {
        let mut line = format!("{:.9e}", m.freq_hz / 1e9);
        for &v in &m.data {
            line.push_str(&format!(",{:.8e},{:.8e}", v.re, v.im));
        }
        writeln!(f, "{line}")?;
    }
    Ok(())
}

/// Apply per-port reference-plane de-embedding to an S matrix.
///
/// For each entry, applies a two-sided shift:
///   S'_{ij} = S_{ij} · exp[(α + jβ)·(l_i + l_j)]
/// where β = 2πf/c0·sqrt(eps_eff).
pub fn apply_reference_plane_deembed(
    s: &SMatrix,
    lengths_m: &[f64],
    eps_eff: f64,
    alpha_np_per_m: f64,
) -> RemResult<SMatrix> {
    if lengths_m.len() != s.n_ports {
        return Err(rem_core::RemError::Config(format!(
            "Deembed length count {} does not match n_ports {}",
            lengths_m.len(), s.n_ports,
        )));
    }
    let mut out = s.clone();
    let beta = 2.0 * std::f64::consts::PI * s.freq_hz / C0 * eps_eff.max(1.0e-12).sqrt();
    let gamma = Complex64::new(alpha_np_per_m, beta);
    for i in 0..s.n_ports {
        for j in 0..s.n_ports {
            let idx = i * s.n_ports + j;
            let phase = gamma * Complex64::new(lengths_m[i] + lengths_m[j], 0.0);
            out.data[idx] = s.data[idx] * phase.exp();
        }
    }
    Ok(out)
}

/// Convert single-ended S-matrix to mixed-mode S-matrix for differential pairs.
///
/// `pairs` are 0-based single-ended port indices `(p, n)`.
/// The output ordering is `[d1, c1, d2, c2, ...]`.
pub fn single_ended_to_mixed_mode(
    s: &SMatrix,
    pairs: &[(usize, usize)],
) -> RemResult<SMatrix> {
    if pairs.is_empty() {
        return Err(rem_core::RemError::Config("Mixed-mode conversion requires at least one pair".into()));
    }
    if 2 * pairs.len() != s.n_ports {
        return Err(rem_core::RemError::Config(format!(
            "Mixed-mode conversion requires all ports paired: got {} ports, {} pairs",
            s.n_ports,
            pairs.len(),
        )));
    }

    let mut seen = vec![false; s.n_ports];
    for &(p, n) in pairs {
        if p >= s.n_ports || n >= s.n_ports || p == n {
            return Err(rem_core::RemError::Config("Invalid differential pair indices".into()));
        }
        if seen[p] || seen[n] {
            return Err(rem_core::RemError::Config("Differential pair indices overlap".into()));
        }
        seen[p] = true;
        seen[n] = true;
    }

    let n = s.n_ports;
    let inv_sqrt2 = 1.0 / 2.0_f64.sqrt();
    let mut m = DMatrix::<Complex64>::zeros(n, n);
    for (k, &(p, nneg)) in pairs.iter().enumerate() {
        // Differential mode row: (p - n)/sqrt(2)
        let rd = 2 * k;
        m[(rd, p)] = Complex64::new(inv_sqrt2, 0.0);
        m[(rd, nneg)] = Complex64::new(-inv_sqrt2, 0.0);
        // Common mode row: (p + n)/sqrt(2)
        let rc = rd + 1;
        m[(rc, p)] = Complex64::new(inv_sqrt2, 0.0);
        m[(rc, nneg)] = Complex64::new(inv_sqrt2, 0.0);
    }

    let s_se = DMatrix::from_row_slice(n, n, &s.data);
    let m_t = m.transpose();
    let s_mm = &m * s_se * m_t;

    let data: Vec<Complex64> = (0..n)
        .flat_map(|r| (0..n).map(|c| s_mm[(r, c)]).collect::<Vec<_>>())
        .collect();
    Ok(SMatrix { n_ports: n, freq_hz: s.freq_hz, data })
}

// ── Transmission-line RLGC extraction ─────────────────────────────────────

/// Per-unit-length RLGC parameters extracted from a 2-port S-matrix.
#[derive(Debug, Clone)]
pub struct TlineParams {
    pub freq_hz:  f64,
    /// Characteristic impedance [Ω] (complex).
    pub z0_tl:    num_complex::Complex64,
    /// Propagation constant γ = α + jβ [1/m].
    pub gamma:    num_complex::Complex64,
    /// R per unit length [Ω/m].
    pub r_per_m:  f64,
    /// L per unit length [H/m].
    pub l_per_m:  f64,
    /// G per unit length [S/m].
    pub g_per_m:  f64,
    /// C per unit length [F/m].
    pub c_per_m:  f64,
}

/// Extract RLGC per-unit-length parameters from 2-port S-matrices.
///
/// Uses the ABCD matrix formulation:
/// A = ((1+S11)(1-S22) + S12·S21) / (2·S21)
/// Z0_TL = √(B/C),  γ·ℓ = acosh(A)
///
/// Valid only for 2-port networks (returns empty `Vec` otherwise).
pub fn extract_tline_rlgc(s_matrices: &[SMatrix], z0_ref: f64, length_m: f64) -> Vec<TlineParams> {
    use std::f64::consts::PI;
    if length_m <= 0.0 { return vec![]; }

    s_matrices.iter().filter_map(|s| {
        if s.n_ports != 2 { return None; }
        let s11 = s.data[0];
        let s12 = s.data[1];
        let s21 = s.data[2];
        let s22 = s.data[3];
        if s21.norm() < 1e-30 { return None; }

        let z0r = Complex64::new(z0_ref, 0.0);
        let one = Complex64::new(1.0, 0.0);

        // ABCD matrix
        let two_s21 = Complex64::new(2.0, 0.0) * s21;
        let a_abcd = ((one + s11) * (one - s22) + s12 * s21) / two_s21;
        let b_abcd = z0r * ((one + s11) * (one + s22) - s12 * s21) / two_s21;
        let c_abcd = ((one - s11) * (one - s22) - s12 * s21) / (two_s21 * z0r);

        // Characteristic impedance Z0_TL = √(B/C)
        let bc = b_abcd / c_abcd;
        // sqrt for Complex64: use (re+j·im)^{1/2}
        let z0_tl = bc.sqrt();

        // Propagation constant: γ·ℓ = acosh(A) = ln(A + √(A²-1))
        let a_sq_m1 = a_abcd * a_abcd - one;
        let gamma_l = (a_abcd + a_sq_m1.sqrt()).ln();
        let gamma = gamma_l / Complex64::new(length_m, 0.0);

        let omega = 2.0 * PI * s.freq_hz;
        let rl = gamma * z0_tl; // (R + jωL)
        let gy = gamma / z0_tl; // (G + jωC)

        Some(TlineParams {
            freq_hz: s.freq_hz,
            z0_tl,
            gamma,
            r_per_m:  rl.re,
            l_per_m:  if omega > 0.0 { rl.im / omega } else { 0.0 },
            g_per_m:  gy.re,
            c_per_m:  if omega > 0.0 { gy.im / omega } else { 0.0 },
        })
    }).collect()
}

/// Write RLGC parameters to `tline_params.csv`.
pub fn write_tline_csv(params: &[TlineParams], path: &Path) -> RemResult<()> {
    use std::io::Write;
    if params.is_empty() { return Ok(()); }
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() { std::fs::create_dir_all(parent)?; }
    }
    let mut f = std::fs::File::create(path)?;
    writeln!(f, "Freq (GHz),Re(Z0_TL),Im(Z0_TL),Re(gamma) [1/m],Im(gamma) [1/m],R (Ohm/m),L (H/m),G (S/m),C (F/m)")?;
    for p in params {
        writeln!(f,
            "{:.9e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e}",
            p.freq_hz / 1e9,
            p.z0_tl.re, p.z0_tl.im,
            p.gamma.re,  p.gamma.im,
            p.r_per_m,   p.l_per_m,
            p.g_per_m,   p.c_per_m,
        )?;
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
            global_node_ids: vec![],
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

    // ── Phase 20 tests ────────────────────────────────────────────────────

    /// Round-trip S→Z→Y: for a matched load S11=0 (Z11=Z0, Y11=1/Z0).
    #[test]
    fn s_to_z_matched_load() {
        let z0 = 50.0;
        // S11 = 0 → Z11 = Z0 · (1+0)(1-0)^{-1} = Z0
        let s = SMatrix { n_ports: 1, freq_hz: 1e9, data: vec![Complex64::ZERO] };
        let z = s_to_z(&s, z0).unwrap();
        assert!((z.data[0].re - z0).abs() < 1e-10, "Z11 should equal Z0={z0}");
        assert!(z.data[0].im.abs() < 1e-10);
    }

    /// Round-trip Z→Y→Z recovers original Z.
    #[test]
    fn z_to_y_round_trip() {
        let z0 = 50.0;
        // S11 = -0.5 (reactive load)
        let s = SMatrix { n_ports: 1, freq_hz: 1e9, data: vec![Complex64::new(-0.5, 0.0)] };
        let z = s_to_z(&s, z0).unwrap();
        let y = z_to_y(&z).unwrap();
        // Y = 1/Z
        let z_re = z.data[0];
        let y_re = y.data[0];
        assert!((z_re * y_re - Complex64::new(1.0, 0.0)).norm() < 1e-10);
    }

    /// extract_tline_rlgc: lossless line (Z0=50Ω, β*ℓ=π/4) → L/C should give √(L/C)≈50.
    #[test]
    fn tline_rlgc_lossless_50ohm() {
        use std::f64::consts::PI;
        let z0_ref = 50.0_f64;
        let z0_tl  = 50.0_f64;
        let freq   = 1e9_f64;
        let length = 0.03; // 30 mm ≈ λ/10 at 1 GHz in free space
        let beta   = 2.0 * PI * freq / rem_core::C0;
        let theta  = beta * length;

        // S-parameters of a lossless 50Ω transmission line section
        // S11 = S22 = 0,  S21 = S12 = exp(-j θ)
        let s21 = Complex64::from_polar(1.0, -theta);
        let s = SMatrix {
            n_ports: 2, freq_hz: freq,
            data: vec![Complex64::ZERO, s21, s21, Complex64::ZERO],
        };
        let result = extract_tline_rlgc(&[s], z0_ref, length);
        assert_eq!(result.len(), 1);
        let p = &result[0];
        // Z0_TL ≈ 50 Ω (real)
        assert!((p.z0_tl.re - z0_tl).abs() / z0_tl < 0.01,
            "Z0_TL re = {:.2}, expected ≈ {z0_tl}", p.z0_tl.re);
        assert!(p.z0_tl.im.abs() < 1.0, "Z0_TL im should be ≈0");
        // R ≈ 0, G ≈ 0 (lossless)
        assert!(p.r_per_m.abs() < 0.1, "R = {:.4e} should be ≈0", p.r_per_m);
        assert!(p.g_per_m.abs() < 1e-6, "G = {:.4e} should be ≈0", p.g_per_m);
        // √(L/C) ≈ Z0_TL
        if p.l_per_m > 0.0 && p.c_per_m > 0.0 {
            let z0_lc = (p.l_per_m / p.c_per_m).sqrt();
            assert!((z0_lc - z0_tl).abs() / z0_tl < 0.01,
                "√(L/C) = {z0_lc:.2}, expected ≈ {z0_tl}");
        }
    }

    #[test]
    fn deembed_identity_when_zero_lengths() {
        let s = SMatrix {
            n_ports: 2,
            freq_hz: 1.0e9,
            data: vec![
                Complex64::new(0.1, 0.2), Complex64::new(0.3, -0.1),
                Complex64::new(0.4, 0.5), Complex64::new(-0.2, 0.05),
            ],
        };
        let out = apply_reference_plane_deembed(&s, &[0.0, 0.0], 1.0, 0.0).unwrap();
        for i in 0..s.data.len() {
            assert!((out.data[i] - s.data[i]).norm() < 1.0e-12);
        }
    }

    #[test]
    fn mixed_mode_identity_for_through_pair() {
        // Single-ended 2-port ideal through: S21=S12=1, S11=S22=0
        let s = SMatrix {
            n_ports: 2,
            freq_hz: 1.0e9,
            data: vec![
                Complex64::ZERO, Complex64::new(1.0, 0.0),
                Complex64::new(1.0, 0.0), Complex64::ZERO,
            ],
        };
        let mm = single_ended_to_mixed_mode(&s, &[(0, 1)]).unwrap();
        // dd and cc both through, dc/cd ~ 0 for symmetric pair
        assert!((mm.data[0] - Complex64::new(-1.0, 0.0)).norm() < 1.0e-12 ||
                (mm.data[0] - Complex64::new(1.0, 0.0)).norm() < 1.0e-12);
        assert!(mm.data[1].norm() < 1.0e-12);
        assert!(mm.data[2].norm() < 1.0e-12);
    }
}
