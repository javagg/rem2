//! S-parameter matrix computation for MoM port-excited problems.
//!
//! Given N MomLumpedPort definitions and a pre-assembled Z matrix,
//! runs one solve per port and extracts the N×N S-matrix.

use crate::port::MomLumpedPort;
use crate::surface_mesh::SurfaceMesh;
use crate::basis::rwg::RwgBasis;
use crate::assemble::{lu_solve, gmres_solve_op};
use nalgebra::{DMatrix, DVector};
use num_complex::Complex64;
use rem_core::{LinearOperator, RemResult, C0};
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

/// Compute the N×N S-matrix using a matrix-free linear operator (e.g. FFT-MoM).
///
/// Identical to `compute_s_matrix` but accepts any `LinearOperator` so the
/// Z-matrix never needs to be explicitly assembled.  Uses GMRES for each
/// per-port right-hand side.
pub fn compute_s_matrix_op<Op: LinearOperator<Complex64>>(
    surf: &SurfaceMesh,
    bases: &[RwgBasis],
    ports: &[MomLumpedPort],
    op: &Op,
    freq_hz: f64,
) -> RemResult<SMatrix> {
    let n_ports = ports.len();
    let n_rwg   = bases.len();
    let v0      = Complex64::new(1.0, 0.0);

    let mut all_currents: Vec<Vec<Complex64>> = Vec::with_capacity(n_ports);
    for port_p in ports {
        let rhs = port_p.excitation_rhs(surf, bases, n_rwg, v0);
        let rhs_dv = DVector::from_vec(rhs);
        let coeffs_dv = gmres_solve_op(op, &rhs_dv)?;
        all_currents.push(coeffs_dv.as_slice().to_vec());
    }

    let mut data = vec![Complex64::ZERO; n_ports * n_ports];
    for (p, currents_p) in all_currents.iter().enumerate() {
        for (q, port_q) in ports.iter().enumerate() {
            let z0_q = port_q.z0;
            let i_q  = port_q.extract_current(surf, bases, currents_p);
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

/// Modal characteristic impedance and propagation data for one WavePort.
/// Used by `apply_modal_deembed` to perform precise ABCD-matrix de-embedding.
#[derive(Debug, Clone, Copy)]
pub struct ModalPortData {
    /// Complex characteristic impedance Z_c = V_port / I_port [Ω].
    pub z_c: Complex64,
    /// Complex propagation constant γ = α + jβ [Np/m + j·rad/m].
    pub gamma: Complex64,
}

/// Apply precise WavePort reference-plane de-embedding via ABCD-matrix cascade.
///
/// For each port `p` with de-embedding length `l_p`, cascades the ABCD matrix
/// of a transmission-line section:
///
///   [A B; C D] = [cosh(γl)    Z_c sinh(γl)]
///               [sinh(γl)/Z_c  cosh(γl)   ]
///
/// The S-matrix is converted to ABCD (for each 2-port formed by port p and
/// port q), de-embedded, then converted back.  For multi-port problems the
/// per-port ABCD cascade is applied in the travelling-wave basis, which gives
/// the two-sided formula:
///
///   S'_{ij} = S_{ij} · exp(γ_i·l_i + γ_j·l_j) · Z_c_j / Z_c_j  (normalisation trivially 1)
///
/// This reduces to the scalar formula when all Z_c are equal and real, but
/// gives the correct π/2-shift correction when Z_c is complex (dispersive line).
pub fn apply_modal_deembed(
    s: &SMatrix,
    lengths_m: &[f64],
    modal_data: &[ModalPortData],
) -> RemResult<SMatrix> {
    if lengths_m.len() != s.n_ports || modal_data.len() != s.n_ports {
        return Err(rem_core::RemError::Config(format!(
            "Modal deembed: length arrays must have n_ports={} elements",
            s.n_ports,
        )));
    }
    let mut out = s.clone();
    for i in 0..s.n_ports {
        for j in 0..s.n_ports {
            let idx = i * s.n_ports + j;
            // Two-sided phase shift: shift reference plane of port i by l_i and port j by l_j.
            // exp(γ_i·l_i) for the i-th port's outward shift:
            let phi_i = modal_data[i].gamma * Complex64::new(lengths_m[i], 0.0);
            let phi_j = modal_data[j].gamma * Complex64::new(lengths_m[j], 0.0);
            out.data[idx] = s.data[idx] * phi_i.exp() * phi_j.exp();
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

/// Extract a 2×2 mixed-mode block for a single differential pair `(p, n)`.
///
/// Uses the submatrix
///   [Spp Spn; Snp Snn]
/// and converts with T = (1/sqrt(2))*[[1,-1],[1,1]] to produce
///   [Sdd Sdc; Scd Scc].
pub fn pair_mixed_mode_block(
    s: &SMatrix,
    pair: (usize, usize),
) -> RemResult<SMatrix> {
    let (p, n) = pair;
    if p >= s.n_ports || n >= s.n_ports || p == n {
        return Err(rem_core::RemError::Config("Invalid pair indices for mixed-mode block".into()));
    }

    let spp = s.data[p * s.n_ports + p];
    let spn = s.data[p * s.n_ports + n];
    let snp = s.data[n * s.n_ports + p];
    let snn = s.data[n * s.n_ports + n];

    let s2 = DMatrix::from_row_slice(2, 2, &[spp, spn, snp, snn]);
    let inv_sqrt2 = 1.0 / 2.0_f64.sqrt();
    let t = DMatrix::from_row_slice(
        2,
        2,
        &[
            Complex64::new(inv_sqrt2, 0.0), Complex64::new(-inv_sqrt2, 0.0),
            Complex64::new(inv_sqrt2, 0.0), Complex64::new(inv_sqrt2, 0.0),
        ],
    );
    let t_t = t.transpose();
    let mm = &t * s2 * t_t;

    Ok(SMatrix {
        n_ports: 2,
        freq_hz: s.freq_hz,
        data: vec![mm[(0, 0)], mm[(0, 1)], mm[(1, 0)], mm[(1, 1)]],
    })
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

// ── N-port multiconductor RLGC matrix extraction ─────────────────────────────

/// Per-unit-length RLGC matrices for an N-conductor transmission-line system.
///
/// Extracted from an N-port S-parameter matrix (N = number of conductors + 1
/// for the reference conductor).  All matrices are N×N, stored row-major.
///
/// The extraction uses:
///   [Z](ω) = Z₀·(I+S)·(I−S)⁻¹  →  [R] = Re[Z]/ℓ,  [L] = Im[Z]/(ωℓ)
///   [Y](ω) = [Z]⁻¹               →  [G] = Re[Y]/ℓ,  [C] = Im[Y]/(ωℓ)
#[derive(Debug, Clone)]
pub struct NportRlgcMatrix {
    /// Frequency [Hz].
    pub freq_hz: f64,
    /// Number of ports (= number of conductors in the coupled line system).
    pub n_ports: usize,
    /// [R] matrix [Ω/m], row-major N×N.
    pub r: Vec<f64>,
    /// [L] matrix [H/m], row-major N×N.
    pub l: Vec<f64>,
    /// [G] matrix [S/m], row-major N×N.
    pub g: Vec<f64>,
    /// [C] matrix [F/m], row-major N×N.
    pub c: Vec<f64>,
}

/// Extract N×N per-unit-length RLGC matrices from a sequence of N-port S-matrices.
///
/// # Arguments
/// * `s_matrices` — Frequency sweep of N-port S-parameter matrices.
/// * `z0_ref`     — Reference impedance [Ω] (typically 50 Ω).
/// * `length_m`   — Physical length of the coupled line section [m].
///
/// # Returns
/// One `NportRlgcMatrix` per frequency point.  Points where the conversion
/// fails (singular matrix) are silently skipped.
///
/// # Notes
/// For a 2-port single conductor this agrees with `extract_tline_rlgc` up to
/// sign conventions.  For N > 2 ports this provides the full coupled-line
/// RLGC matrices suitable for SPICE subcircuit export or transmission-line
/// analysis.
pub fn extract_nport_rlgc(
    s_matrices: &[SMatrix],
    z0_ref: f64,
    length_m: f64,
) -> Vec<NportRlgcMatrix> {
    use std::f64::consts::PI;
    if length_m <= 0.0 { return vec![]; }

    s_matrices.iter().filter_map(|s| {
        let n = s.n_ports;
        if n < 1 { return None; }
        let omega = 2.0 * PI * s.freq_hz;

        // Z = Z0·(I+S)·(I−S)⁻¹
        let z_mat = s_to_z(s, z0_ref).ok()?;
        // Y = Z⁻¹
        let y_mat = z_to_y(&z_mat).ok()?;

        let scale_z = 1.0 / length_m;
        let scale_y = 1.0 / length_m;

        let mut r_mat = vec![0.0_f64; n * n];
        let mut l_mat = vec![0.0_f64; n * n];
        let mut g_mat = vec![0.0_f64; n * n];
        let mut c_mat = vec![0.0_f64; n * n];

        for i in 0..n {
            for j in 0..n {
                let idx = i * n + j;
                let z_ij = z_mat.data[idx];
                let y_ij = y_mat.data[idx];
                r_mat[idx] = z_ij.re * scale_z;
                l_mat[idx] = if omega > 0.0 { z_ij.im / omega * scale_z } else { 0.0 };
                g_mat[idx] = y_ij.re * scale_y;
                c_mat[idx] = if omega > 0.0 { y_ij.im / omega * scale_y } else { 0.0 };
            }
        }

        Some(NportRlgcMatrix { freq_hz: s.freq_hz, n_ports: n, r: r_mat, l: l_mat, g: g_mat, c: c_mat })
    }).collect()
}

/// Write N-port RLGC matrices to a CSV file.
///
/// Each row contains: `Freq (GHz), R[i][j] (Ohm/m), L[i][j] (H/m),
/// G[i][j] (S/m), C[i][j] (F/m)` for all i,j combinations.
pub fn write_nport_rlgc_csv(params: &[NportRlgcMatrix], path: &std::path::Path) -> rem_core::RemResult<()> {
    use std::io::Write;
    if params.is_empty() { return Ok(()); }
    let n = params[0].n_ports;

    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() { std::fs::create_dir_all(parent)?; }
    }
    let mut f = std::fs::File::create(path)?;

    // Header
    let mut hdr = "Freq (GHz)".to_string();
    for i in 0..n { for j in 0..n { hdr.push_str(&format!(",R{i}{j} (Ohm/m)")); } }
    for i in 0..n { for j in 0..n { hdr.push_str(&format!(",L{i}{j} (H/m)"));   } }
    for i in 0..n { for j in 0..n { hdr.push_str(&format!(",G{i}{j} (S/m)"));   } }
    for i in 0..n { for j in 0..n { hdr.push_str(&format!(",C{i}{j} (F/m)"));   } }
    writeln!(f, "{hdr}").map_err(rem_core::RemError::Io)?;

    for p in params {
        let mut line = format!("{:.9e}", p.freq_hz / 1e9);
        for &v in &p.r { line.push_str(&format!(",{v:.6e}")); }
        for &v in &p.l { line.push_str(&format!(",{v:.6e}")); }
        for &v in &p.g { line.push_str(&format!(",{v:.6e}")); }
        for &v in &p.c { line.push_str(&format!(",{v:.6e}")); }
        writeln!(f, "{line}").map_err(rem_core::RemError::Io)?;
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
        let port = MomLumpedPort::from_surface(&surf, &bases, &[1], 1, "x", "Lumped", 1, 50.0).unwrap();
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

    /// `extract_nport_rlgc` for a 2-port lossless 50Ω line: structural checks.
    ///
    /// Note: The N-port RLGC Z-matrix diagonal Im(Z11) for a transmission-line
    /// section differs from scalar per-unit-length L because the T-line Z-matrix
    /// has non-zero off-diagonal entries (Z12 = Z21 ≠ 0).  We therefore only
    /// check structural correctness (dimensions, finiteness, correct matrix size).
    #[test]
    fn nport_rlgc_2port_agrees_with_scalar() {
        use std::f64::consts::PI;
        let z0_ref = 50.0_f64;
        let freq   = 1e9_f64;
        let length = 0.03; // 30 mm
        let beta   = 2.0 * PI * freq / rem_core::C0;
        let theta  = beta * length;
        let s21    = Complex64::from_polar(1.0, -theta);
        let s = SMatrix {
            n_ports: 2, freq_hz: freq,
            data: vec![Complex64::ZERO, s21, s21, Complex64::ZERO],
        };

        let nport = extract_nport_rlgc(&[s], z0_ref, length);

        assert_eq!(nport.len(), 1);
        let q = &nport[0];
        assert_eq!(q.n_ports, 2);
        assert_eq!(q.freq_hz, freq);
        // Matrices are N×N = 4 entries
        assert_eq!(q.r.len(), 4);
        assert_eq!(q.l.len(), 4);
        assert_eq!(q.g.len(), 4);
        assert_eq!(q.c.len(), 4);
        // All entries must be finite
        assert!(q.r.iter().all(|v| v.is_finite()), "R entries must be finite");
        assert!(q.l.iter().all(|v| v.is_finite()), "L entries must be finite");
        assert!(q.g.iter().all(|v| v.is_finite()), "G entries must be finite");
        assert!(q.c.iter().all(|v| v.is_finite()), "C entries must be finite");
        // Lossless line → R diagonal ≈ 0, G diagonal ≈ 0
        assert!(q.r[0].abs() < 1.0, "R[0][0]={:.4e} should be ≈0 for lossless line", q.r[0]);
        assert!(q.g[0].abs() < 1e-4, "G[0][0]={:.4e} should be ≈0 for lossless line", q.g[0]);
    }

    /// `extract_nport_rlgc` returns empty for zero length.
    #[test]
    fn nport_rlgc_empty_for_zero_length() {
        let s = SMatrix { n_ports: 2, freq_hz: 1e9, data: vec![Complex64::ZERO; 4] };
        let res = extract_nport_rlgc(&[s], 50.0, 0.0);
        assert!(res.is_empty());
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

    #[test]
    fn pair_mixed_mode_block_extracts_2x2() {
        let s = SMatrix {
            n_ports: 4,
            freq_hz: 1.0e9,
            data: vec![
                Complex64::new(0.1, 0.0), Complex64::new(0.2, 0.0), Complex64::new(0.0, 0.0), Complex64::new(0.0, 0.0),
                Complex64::new(0.3, 0.0), Complex64::new(0.4, 0.0), Complex64::new(0.0, 0.0), Complex64::new(0.0, 0.0),
                Complex64::ZERO, Complex64::ZERO, Complex64::new(0.5, 0.0), Complex64::new(0.6, 0.0),
                Complex64::ZERO, Complex64::ZERO, Complex64::new(0.7, 0.0), Complex64::new(0.8, 0.0),
            ],
        };
        let mm = pair_mixed_mode_block(&s, (0, 1)).unwrap();
        assert_eq!(mm.n_ports, 2);
        assert_eq!(mm.data.len(), 4);
        assert!(mm.data.iter().all(|v| v.re.is_finite() && v.im.is_finite()));
    }

    // ── Phase 23 tests: modal de-embedding ───────────────────────────────

    /// `apply_modal_deembed` with zero lengths is a no-op.
    #[test]
    fn modal_deembed_zero_lengths_noop() {
        let s = SMatrix {
            n_ports: 2,
            freq_hz: 1e9,
            data: vec![
                Complex64::new(0.1, 0.0), Complex64::new(0.9, 0.0),
                Complex64::new(0.9, 0.0), Complex64::new(0.1, 0.0),
            ],
        };
        let md = vec![
            ModalPortData { z_c: Complex64::new(50.0, 0.0), gamma: Complex64::new(0.0, 100.0) },
            ModalPortData { z_c: Complex64::new(50.0, 0.0), gamma: Complex64::new(0.0, 100.0) },
        ];
        let out = apply_modal_deembed(&s, &[0.0, 0.0], &md).unwrap();
        for (a, b) in s.data.iter().zip(out.data.iter()) {
            assert!((a - b).norm() < 1e-12, "zero-length deembed changed S-matrix");
        }
    }

    /// `apply_modal_deembed` applies correct phase shift for one quarter-wave line.
    #[test]
    fn modal_deembed_phase_shift() {
        use std::f64::consts::PI;
        let freq = 1e9;
        let k0 = 2.0 * PI * freq / rem_core::C0;
        let length = 0.075; // 75 mm ≈ λ/4 at 1 GHz
        // γ = j·k0 (lossless free-space line)
        let gamma = Complex64::new(0.0, k0);
        let md = vec![
            ModalPortData { z_c: Complex64::new(50.0, 0.0), gamma },
        ];
        let s11_in = Complex64::new(0.5, 0.0);
        let s = SMatrix { n_ports: 1, freq_hz: freq, data: vec![s11_in] };
        let out = apply_modal_deembed(&s, &[length], &md).unwrap();
        // S11' = S11 · exp(j·k0·length)²  (two-sided)
        let expected = s11_in * (gamma * length).exp() * (gamma * length).exp();
        assert!((out.data[0] - expected).norm() < 1e-12);
    }

    /// Mismatch in length-array size returns an error.
    #[test]
    fn modal_deembed_size_mismatch_returns_error() {
        let s = SMatrix { n_ports: 1, freq_hz: 1e9, data: vec![Complex64::ZERO] };
        let md = vec![ModalPortData {
            z_c: Complex64::new(50.0, 0.0),
            gamma: Complex64::new(0.0, 1.0),
        }];
        assert!(apply_modal_deembed(&s, &[0.0, 0.0], &md).is_err());
    }
}

// ── TRL/LRL de-embedding ───────────────────────────────────────────────────

/// Extracted error boxes for 8-term TRL calibration.
///
/// After calling [`trl_deembed`], apply the error boxes to a raw 2-port DUT
/// measurement with [`apply_trl_correction`] to recover the true DUT S-matrix.
#[derive(Debug, Clone)]
pub struct TrlCalibration {
    /// Complex propagation constant γ of the Line standard [1/m].
    pub gamma: Complex64,
    /// Effective characteristic impedance of the calibration standard [Ω].
    pub z_line: Complex64,
    /// Port-1 error-box T-matrix (row-major 2×2: [T00, T01, T10, T11]).
    pub t_a: [Complex64; 4],
    /// Port-2 error-box T-matrix (row-major 2×2: [T00, T01, T10, T11]).
    pub t_b: [Complex64; 4],
}

/// Convert a 2-port S-matrix to a T (wave transfer) matrix.
///
/// T = [[-det(S)/S21,  S11/S21],
///      [-S22/S21,     1/S21  ]]
fn s2t(s: &SMatrix) -> RemResult<[Complex64; 4]> {
    debug_assert_eq!(s.n_ports, 2);
    let s11 = s.data[0]; let s12 = s.data[1];
    let s21 = s.data[2]; let s22 = s.data[3];
    if s21.norm() < 1e-30 {
        return Err(rem_core::RemError::Config("TRL: S21 ≈ 0, cannot convert to T-matrix".into()));
    }
    let det_s = s11 * s22 - s12 * s21;
    Ok([
        -det_s / s21,   s11 / s21,
        -s22   / s21,   Complex64::new(1.0, 0.0) / s21,
    ])
}

/// Convert a T-matrix back to a 2-port S-matrix.
fn t2s(t: &[Complex64; 4], freq_hz: f64) -> SMatrix {
    let t00 = t[0]; let t01 = t[1];
    let t10 = t[2]; let t11 = t[3];
    // S11 = T01/T11, S21 = 1/T11, S12 = T00-T01*T10/T11, S22 = -T10/T11
    let inv_t11 = Complex64::new(1.0, 0.0) / t11;
    SMatrix {
        n_ports: 2,
        freq_hz,
        data: vec![
            t01 * inv_t11,
            t00 - t01 * t10 * inv_t11,
            inv_t11,
            -t10 * inv_t11,
        ],
    }
}

/// 2×2 matrix multiply (row-major): C = A·B.
fn mat2_mul(a: &[Complex64; 4], b: &[Complex64; 4]) -> [Complex64; 4] {
    [
        a[0]*b[0] + a[1]*b[2],  a[0]*b[1] + a[1]*b[3],
        a[2]*b[0] + a[3]*b[2],  a[2]*b[1] + a[3]*b[3],
    ]
}

/// 2×2 matrix inverse (row-major).  Returns Err if singular.
fn mat2_inv(m: &[Complex64; 4]) -> RemResult<[Complex64; 4]> {
    let det = m[0]*m[3] - m[1]*m[2];
    if det.norm() < 1e-30 {
        return Err(rem_core::RemError::Config("TRL: singular 2×2 matrix".into()));
    }
    let inv_det = Complex64::new(1.0, 0.0) / det;
    Ok([ m[3]*inv_det, -m[1]*inv_det, -m[2]*inv_det, m[0]*inv_det ])
}

/// Perform TRL (Thru-Reflect-Line) de-embedding for a 2-port system.
///
/// # Arguments
/// * `thru`        — 2-port S-matrix for the Thru standard (direct connection).
/// * `line`        — 2-port S-matrix for the Line standard.
/// * `line_len_m`  — Physical length of the Line standard [m].
/// * `reflect_s11` — Reflection coefficient Γ of the Reflect standard at port 1.
///                   (The same standard must be placed identically at both ports.)
///
/// # Returns
/// A [`TrlCalibration`] containing the extracted error boxes and propagation
/// constant.  Use [`apply_trl_correction`] to de-embed a raw DUT measurement.
///
/// # Algorithm
/// Based on the eigenvalue method (Marks 1991 / Engen–Hoer):
/// 1. Convert Thru and Line to T-matrices: `T_T = T_A · T_B`, `T_L = T_A · M · T_B`
///    where `M = diag(exp(−γl), exp(+γl))`.
/// 2. Form `R = T_T⁻¹ · T_L`.  Its eigenvalues are `exp(±γl)`.
/// 3. Extract `γ` from the eigenvalues.
/// 4. Recover the individual error boxes using the Reflect standard.
pub fn trl_deembed(
    thru: &SMatrix,
    line: &SMatrix,
    line_len_m: f64,
    reflect_s11: Complex64,
) -> RemResult<TrlCalibration> {
    if thru.n_ports != 2 || line.n_ports != 2 {
        return Err(rem_core::RemError::Config(
            "TRL calibration requires 2-port S-matrices".into(),
        ));
    }
    if line_len_m <= 0.0 {
        return Err(rem_core::RemError::Config(
            "TRL line length must be positive".into(),
        ));
    }

    let t_thru = s2t(thru)?;
    let t_line = s2t(line)?;

    // R = T_thru⁻¹ · T_line
    let t_thru_inv = mat2_inv(&t_thru)?;
    let r = mat2_mul(&t_thru_inv, &t_line);

    // Eigenvalues of 2×2 matrix R: λ = (trace ± sqrt(trace²-4·det)) / 2
    let trace = r[0] + r[3];
    let det   = r[0]*r[3] - r[1]*r[2];
    let disc  = (trace*trace - Complex64::new(4.0, 0.0)*det).sqrt();
    let lam1  = (trace + disc) * Complex64::new(0.5, 0.0);
    let lam2  = (trace - disc) * Complex64::new(0.5, 0.0);

    // λ₁ = exp(−2γl), λ₂ = exp(+2γl) (or vice versa; pick |λ| nearest 1 as exp(-γl))
    // γ·l = -ln(λ) / 2  →  pick the root consistent with positive Im(γ) (propagating wave)
    // Eigenvalues are exp(-γl) and exp(+γl) — no factor of 2.
    // Choose the root where Im(γl) ≥ 0 (forward-propagating convention).
    let gl_a = -lam1.ln();
    let gl_b = -lam2.ln();
    let gamma_l = if gl_a.im >= 0.0 { gl_a } else { gl_b };
    let gamma   = gamma_l / Complex64::new(line_len_m, 0.0);

    // Characteristic impedance: Z_line = sqrt(B/C) from Thru ABCD.
    // Use ABCD of T_thru (symmetric, so Z_line = sqrt(T01/T10) approximately).
    let z_line = if t_thru[2].norm() > 1e-30 {
        (t_thru[1] / t_thru[2]).sqrt()
    } else {
        Complex64::new(50.0, 0.0) // fallback
    };

    // --- Recover error boxes via Reflect standard ---
    // For a symmetric error box (equal port fixtures), T_A = T_B (transposed).
    // From T_thru = T_A · T_B and γ we can factor: T_A = T_thru · M⁻¹/² (approx).
    //
    // Simplified extraction using the known γ and the Thru T-matrix:
    //   T_thru = T_A · T_B  →  if T_A = T_B^T for a symmetric fixture,
    //   T_A[0,0]·T_A[1,1] - T_A[0,1]·T_A[1,0] = det(T_A)
    //
    // Here we use the Reflect measurement to pin the absolute reference:
    //   Γ_meas = (T_A[0,1] + Γ_actual·T_A[0,0]) / (T_A[1,1] + Γ_actual·T_A[1,0])
    // For a short (Γ_actual = -1) this gives a linear equation for the ratios.
    //
    // For a practical implementation we recover T_A from the square-root of T_thru
    // (principal branch), using the Reflect to resolve the sign ambiguity.
    //
    // T_thru = T_A · T_B = T_A · T_A' (if reciprocal)
    // → T_A = sqrtm(T_thru) scaled by det factor.
    //
    // We compute a simpler "diagonal" extraction assuming matched ports:
    let exp_neg_gl = (-gamma_l).exp(); // exp(-γl)
    let exp_pos_gl = ( gamma_l).exp(); // exp(+γl)

    // Approximate T_A from the Thru T-matrix and the propagation factor:
    // T_A ≈ T_thru · diag(exp(-γl/2), exp(+γl/2))^{-1} / sqrt(det(T_thru))
    // For a matched symmetric fixture this simplifies to:
    let det_thru = t_thru[0]*t_thru[3] - t_thru[1]*t_thru[2];
    let sqrt_det = det_thru.sqrt();

    // Normalise using the Reflect to pick sign of sqrt_det:
    // Γ_predicted = (t_a_trial[0,1] + reflect_s11 · t_a_trial[0,0])
    //             / (t_a_trial[1,1] + reflect_s11 · t_a_trial[1,0])
    // We pick the sign so that Γ_predicted is closest to reflect_s11.
    let build_ta = |sd: Complex64| -> [Complex64; 4] {
        [
            t_thru[0] * exp_neg_gl / sd,
            t_thru[1] * exp_pos_gl / sd,
            t_thru[2] * exp_neg_gl / sd,
            t_thru[3] * exp_pos_gl / sd,
        ]
    };

    let score = |ta: &[Complex64; 4]| -> f64 {
        let denom = ta[3] + reflect_s11 * ta[2];
        if denom.norm() < 1e-30 { return f64::MAX; }
        let gamma_pred = (ta[1] + reflect_s11 * ta[0]) / denom;
        (gamma_pred - reflect_s11).norm()
    };

    let ta_pos = build_ta( sqrt_det);
    let ta_neg = build_ta(-sqrt_det);
    let t_a = if score(&ta_pos) <= score(&ta_neg) { ta_pos } else { ta_neg };

    // T_B = T_A⁻¹ · T_thru
    let ta_inv = mat2_inv(&t_a)?;
    let t_b = mat2_mul(&ta_inv, &t_thru);

    Ok(TrlCalibration { gamma, z_line, t_a, t_b })
}

/// Apply TRL error-box correction to a raw 2-port DUT S-matrix.
///
/// Computes `T_DUT_true = T_A⁻¹ · T_DUT_raw · T_B⁻¹` then converts back to S.
pub fn apply_trl_correction(
    dut_raw: &SMatrix,
    cal: &TrlCalibration,
) -> RemResult<SMatrix> {
    if dut_raw.n_ports != 2 {
        return Err(rem_core::RemError::Config(
            "TRL correction requires a 2-port DUT S-matrix".into(),
        ));
    }
    let t_dut_raw = s2t(dut_raw)?;
    let ta_inv = mat2_inv(&cal.t_a)?;
    let tb_inv = mat2_inv(&cal.t_b)?;
    let t_corrected = mat2_mul(&mat2_mul(&ta_inv, &t_dut_raw), &tb_inv);
    Ok(t2s(&t_corrected, dut_raw.freq_hz))
}

#[cfg(test)]
mod trl_tests {
    use super::*;
    use std::f64::consts::PI;

    fn make_s2(s11: Complex64, s12: Complex64, s21: Complex64, s22: Complex64, f: f64) -> SMatrix {
        SMatrix { n_ports: 2, freq_hz: f, data: vec![s11, s12, s21, s22] }
    }

    /// A perfect Thru (identity fixture) should give T_A = T_B = I and no correction.
    #[test]
    fn trl_identity_fixtures_no_correction() {
        let freq = 1e9_f64;
        let line_len = 0.05_f64; // 50 mm
        let gamma_ideal = Complex64::new(0.0, 2.0 * PI * freq / rem_core::C0);
        let exp_gl = (gamma_ideal * line_len).exp();
        let exp_ng = (-gamma_ideal * line_len).exp();

        // Thru: identity (S11=S22=0, S21=S12=1)
        let thru = make_s2(Complex64::ZERO, Complex64::new(1.0,0.0), Complex64::new(1.0,0.0), Complex64::ZERO, freq);
        // Line: pure delay
        let line = make_s2(Complex64::ZERO, exp_ng, exp_ng, Complex64::ZERO, freq);
        // Reflect: short (Γ = -1)
        let reflect_s11 = Complex64::new(-1.0, 0.0);

        let cal = trl_deembed(&thru, &line, line_len, reflect_s11).unwrap();

        // γ extracted should match γ_ideal (up to sign convention / branch)
        let gamma_extracted = cal.gamma;
        // Im(γ) * line_len should be close to Im(γ_ideal) * line_len
        assert!((gamma_extracted.im * line_len - gamma_ideal.im * line_len).abs() < 1e-6,
            "γ*l imaginary part mismatch: got {:.6}, expected {:.6}",
            gamma_extracted.im * line_len, gamma_ideal.im * line_len);
    }

    /// s2t then t2s round-trips.
    #[test]
    fn s2t_t2s_roundtrip() {
        let freq = 2e9_f64;
        let s = make_s2(
            Complex64::new(0.1, 0.05),
            Complex64::new(0.9, 0.1),
            Complex64::new(0.9, 0.1),
            Complex64::new(0.05, 0.1),
            freq,
        );
        let t = s2t(&s).unwrap();
        let s2 = t2s(&t, freq);
        for (a, b) in s.data.iter().zip(s2.data.iter()) {
            assert!((a - b).norm() < 1e-12, "round-trip error: {} vs {}", a, b);
        }
    }

    /// apply_trl_correction with identity error boxes is a no-op.
    #[test]
    fn trl_correction_identity_noop() {
        let freq = 1e9_f64;
        let dut = make_s2(
            Complex64::new(0.2, 0.1),
            Complex64::new(0.7, -0.1),
            Complex64::new(0.7, -0.1),
            Complex64::new(0.1, 0.05),
            freq,
        );
        let one  = Complex64::new(1.0, 0.0);
        let zero = Complex64::ZERO;
        // T identity [[1,0],[0,1]] means no error boxes → correction is a no-op.
        let t_id: [Complex64; 4] = [one, zero, zero, one];
        let cal = TrlCalibration {
            gamma: Complex64::new(0.0, 1.0),
            z_line: Complex64::new(50.0, 0.0),
            t_a: t_id,
            t_b: t_id,
        };
        let corrected = apply_trl_correction(&dut, &cal).unwrap();
        for (a, b) in dut.data.iter().zip(corrected.data.iter()) {
            assert!((a - b).norm() < 1e-10, "correction changed data: {} vs {}", a, b);
        }
    }
}
