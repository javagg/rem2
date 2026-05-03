//! Snapshot-based Reduced Order Model (ROM) for MoM frequency sweeps.
//!
//! # Algorithm
//!
//! Given anchor frequencies f₁ … f_r, solve the full N-DOF MoM system at each
//! anchor to obtain current vectors **c**₁ … **c**_r.  Orthonormalize these
//! into a basis **V** ∈ ℂ^{N×r} via modified Gram-Schmidt.  At any target
//! frequency f, project:
//!
//! ```text
//! Z_r(f) = V^H Z(f) V   (r×r)
//! b_r(p) = V^H b_p       (r×1, per port excitation p)
//! x_p(f) = Z_r^{-1} b_r(p)   (r×r solve)
//! c_p(f) ≈ V x_p(f)           (reconstruct N-DOF current)
//! ```
//!
//! The resulting S-matrix is extracted from the approximate currents as usual.
//!
//! # Complexity
//! - Setup:   r × O(N²) for anchor solves + O(r N²) for basis projection
//! - Sweep:   O(r² N) per frequency point (matrix-matrix multiply for Z_r)
//!   versus O(N²) per point for the full system.
//!
//! # Limitations
//! The approximation accuracy degrades for frequency bands wider than the
//! "coherence length" of the Green function.  Always verify against at least
//! one full solve at the band edges.

use crate::basis::rwg::RwgBasis;
use crate::port::MomLumpedPort;
use crate::sparams::{SMatrix, compute_s_matrix};
use crate::surface_mesh::SurfaceMesh;
use nalgebra::DMatrix;
use num_complex::Complex64;
use rem_core::{RemError, RemResult};

// ─── Snapshot basis ──────────────────────────────────────────────────────────

/// Orthonormal snapshot basis for the MoM solution space.
///
/// Built from anchor-frequency current vectors and used to project the
/// system to a low-dimensional subspace.
pub struct MomRom {
    /// Orthonormal basis matrix V: N rows × r columns (one per anchor).
    pub basis: DMatrix<Complex64>,
    /// Anchor frequencies [Hz] in the order they were added.
    pub anchor_freqs: Vec<f64>,
}

impl MomRom {
    /// Build a ROM from current snapshots at the anchor frequencies.
    ///
    /// `snapshots[i]` is the solved current vector at `anchor_freqs[i]`.
    /// Vectors are orthonormalized via modified Gram-Schmidt; linearly
    /// dependent vectors (norm < `tol_rel` × max-norm) are discarded.
    pub fn build(
        anchor_freqs: Vec<f64>,
        snapshots: Vec<Vec<Complex64>>,
        tol_rel: f64,
    ) -> RemResult<Self> {
        if snapshots.is_empty() {
            return Err(RemError::Config("ROM: no snapshots provided".to_string()));
        }
        let n = snapshots[0].len();
        for (i, s) in snapshots.iter().enumerate() {
            if s.len() != n {
                return Err(RemError::Config(format!(
                    "ROM: snapshot {i} has length {}, expected {n}", s.len()
                )));
            }
        }

        // Modified Gram-Schmidt orthonormalization
        let mut qs: Vec<Vec<Complex64>> = Vec::with_capacity(snapshots.len());
        let mut kept_freqs: Vec<f64> = Vec::with_capacity(snapshots.len());

        // Find max-norm snapshot for relative tolerance
        let max_norm = snapshots.iter()
            .map(|s| vec_norm(s))
            .fold(0.0_f64, f64::max);
        let abs_tol = tol_rel * max_norm;

        for (snap, &freq) in snapshots.iter().zip(anchor_freqs.iter()) {
            let mut v: Vec<Complex64> = snap.clone();
            // Subtract projections onto existing basis vectors
            for q in &qs {
                let proj = dot(&v, q);
                for (vi, &qi) in v.iter_mut().zip(q.iter()) {
                    *vi -= proj * qi;
                }
            }
            let nrm = vec_norm(&v);
            if nrm < abs_tol {
                log::debug!("ROM: snapshot at {freq:.3e} Hz is linearly dependent — skipped");
                continue;
            }
            // Normalize
            let inv_n = Complex64::new(1.0 / nrm, 0.0);
            for vi in v.iter_mut() { *vi *= inv_n; }
            qs.push(v);
            kept_freqs.push(freq);
        }

        if qs.is_empty() {
            return Err(RemError::Config(
                "ROM: all snapshots were linearly dependent".to_string()
            ));
        }

        let r = qs.len();
        // Pack into N×r nalgebra matrix
        let basis = DMatrix::from_fn(n, r, |row, col| qs[col][row]);

        log::info!("ROM: built basis with r={r} vectors from {} snapshots (N={n})", snapshots.len());
        Ok(Self { basis, anchor_freqs: kept_freqs })
    }

    /// Project the N×N system matrix `z_full` to the r×r reduced system.
    ///
    /// Z_r = V^H Z V
    pub fn project_system(&self, z_full: &DMatrix<Complex64>) -> DMatrix<Complex64> {
        // Z_r = V^H * Z * V  where V is N×r
        let vt = self.basis.conjugate_transpose();
        &vt * z_full * &self.basis
    }

    /// Project a port RHS vector to the reduced subspace.
    ///
    /// b_r = V^H b
    pub fn project_rhs(&self, rhs: &[Complex64]) -> Vec<Complex64> {
        let r = self.basis.ncols();
        let n = self.basis.nrows();
        assert_eq!(rhs.len(), n);
        (0..r).map(|j| {
            (0..n).map(|i| self.basis[(i, j)].conj() * rhs[i]).sum()
        }).collect()
    }

    /// Reconstruct full-space vector from reduced coefficients.
    ///
    /// c ≈ V x
    pub fn reconstruct(&self, x: &[Complex64]) -> Vec<Complex64> {
        let n = self.basis.nrows();
        let r = self.basis.ncols();
        assert_eq!(x.len(), r);
        (0..n).map(|i| {
            (0..r).map(|j| self.basis[(i, j)] * x[j]).sum()
        }).collect()
    }
}

// ─── ROM frequency sweep ─────────────────────────────────────────────────────

/// Run a MoM S-parameter sweep accelerated by a snapshot ROM.
///
/// The workflow:
/// 1. Solve at `anchor_freqs` (subset of `freq_list`) to build the ROM basis.
/// 2. For every other frequency in `freq_list`, project Z and RHS to the
///    r×r subspace, solve cheaply, and reconstruct approximate currents.
/// 3. Extract S-matrix from approximate currents.
///
/// `build_z`: closure that assembles the N×N impedance matrix at a given frequency.
/// `anchor_freqs`: if `None`, automatically distribute r anchors uniformly over `freq_list`.
///
/// Returns one `SMatrix` per element of `freq_list` in order.
pub fn mom_rom_sweep(
    surf: &SurfaceMesh,
    bases: &[RwgBasis],
    ports: &[MomLumpedPort],
    freq_list: &[f64],
    anchor_count: usize,
    tol_rel: f64,
    build_z: &dyn Fn(f64) -> RemResult<DMatrix<Complex64>>,
) -> RemResult<Vec<SMatrix>> {
    if freq_list.is_empty() {
        return Ok(vec![]);
    }
    let r = anchor_count.max(1).min(freq_list.len());

    // ── 1. Choose anchor frequencies uniformly from freq_list ────────────
    let anchor_indices = uniform_subset(freq_list.len(), r);
    let anchor_freqs: Vec<f64> = anchor_indices.iter().map(|&i| freq_list[i]).collect();

    // ── 2. Full solve at each anchor ─────────────────────────────────────
    let mut snapshots: Vec<Vec<Complex64>> = Vec::with_capacity(r);
    let mut anchor_z_cache: Vec<(f64, DMatrix<Complex64>)> = Vec::with_capacity(r);

    for &freq in &anchor_freqs {
        let z = build_z(freq)?;
        // Solve for port 0 excitation (representative snapshot)
        let rhs = if !ports.is_empty() {
            ports[0].excitation_rhs(surf, bases, bases.len(), Complex64::new(1.0, 0.0))
        } else {
            vec![Complex64::new(1.0, 0.0); bases.len()]
        };
        let currents = crate::assemble::lu_solve(&z, &rhs)?;
        snapshots.push(currents);
        anchor_z_cache.push((freq, z));
    }

    // ── 3. Build ROM basis ────────────────────────────────────────────────
    let rom = MomRom::build(anchor_freqs.clone(), snapshots, tol_rel)?;

    // ── 4. Sweep over all frequencies ────────────────────────────────────
    let mut results = Vec::with_capacity(freq_list.len());
    let anchor_set: std::collections::HashSet<usize> = anchor_indices.iter().cloned().collect();

    for (fi, &freq) in freq_list.iter().enumerate() {
        let s_mat = if anchor_set.contains(&fi) {
            // Use cached full Z for anchor frequencies
            let z = &anchor_z_cache.iter().find(|(f, _)| (*f - freq).abs() < 1e-3).unwrap().1;
            compute_s_matrix(surf, bases, ports, z, freq)?
        } else {
            // Project to reduced subspace
            let z_full = build_z(freq)?;
            let z_r = rom.project_system(&z_full);

            let n_ports = ports.len();
            let mut data = vec![Complex64::ZERO; n_ports * n_ports];

            for (p, port_p) in ports.iter().enumerate() {
                let rhs_full = port_p.excitation_rhs(surf, bases, bases.len(),
                    Complex64::new(1.0, 0.0));
                let rhs_r = rom.project_rhs(&rhs_full);

                // Solve r×r system (small, use Gaussian elimination via nalgebra)
                let rhs_r_dv = nalgebra::DVector::from_vec(rhs_r);
                let z_r_dv = z_r.clone();
                let x_r = z_r_dv.lu().solve(&rhs_r_dv)
                    .ok_or_else(|| RemError::Config("ROM LU solve failed".to_string()))?;

                let currents: Vec<Complex64> = rom.reconstruct(x_r.as_slice());

                for (q, port_q) in ports.iter().enumerate() {
                    let z0_q = port_q.z0;
                    let i_q = port_q.extract_current(surf, bases, &currents);
                    let v_q = if q == p { Complex64::new(1.0, 0.0) } else { Complex64::ZERO };
                    data[q * n_ports + p] = v_q - Complex64::new(z0_q, 0.0) * i_q;
                }
            }
            SMatrix { n_ports, freq_hz: freq, data }
        };
        results.push(s_mat);
    }

    Ok(results)
}

// ─── Helpers ─────────────────────────────────────────────────────────────────

fn vec_norm(v: &[Complex64]) -> f64 {
    v.iter().map(|x| x.norm_sqr()).sum::<f64>().sqrt()
}

fn dot(a: &[Complex64], b: &[Complex64]) -> Complex64 {
    a.iter().zip(b.iter()).map(|(ai, bi)| bi.conj() * ai).sum()
}

/// Select `r` indices uniformly distributed over `[0, n)`.
fn uniform_subset(n: usize, r: usize) -> Vec<usize> {
    if r >= n {
        return (0..n).collect();
    }
    (0..r).map(|i| i * (n - 1) / (r - 1).max(1)).collect()
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use num_complex::Complex64;
    use nalgebra::DMatrix;

    fn make_diag_z(n: usize, scale: Complex64) -> DMatrix<Complex64> {
        let mut z = DMatrix::zeros(n, n);
        for i in 0..n { z[(i, i)] = scale; }
        z
    }

    #[test]
    fn build_rom_from_snapshots() {
        let n = 4;
        let freqs = vec![1e9, 2e9, 3e9];
        let snaps: Vec<Vec<Complex64>> = freqs.iter().enumerate().map(|(i, _)| {
            (0..n).map(|j| Complex64::new((i + 1) as f64 * (j + 1) as f64, 0.0)).collect()
        }).collect();
        let rom = MomRom::build(freqs, snaps, 1e-10).unwrap();
        assert!(rom.basis.ncols() >= 1);
        assert_eq!(rom.basis.nrows(), n);
    }

    #[test]
    fn project_reconstruct_roundtrip() {
        // Single-snapshot ROM: the only basis vector = normalised snapshot
        let n = 4;
        let snap = vec![
            Complex64::new(1.0, 0.0),
            Complex64::new(0.0, 1.0),
            Complex64::new(1.0, 1.0),
            Complex64::new(0.0, 0.0),
        ];
        let nrm: f64 = snap.iter().map(|x| x.norm_sqr()).sum::<f64>().sqrt();
        let rom = MomRom::build(vec![1e9], vec![snap.clone()], 1e-10).unwrap();
        assert_eq!(rom.basis.ncols(), 1);

        // Project then reconstruct should recover scaled version
        let x_r = rom.project_rhs(&snap);
        let rec = rom.reconstruct(&x_r);
        // ||rec||²  ≈  nrm²  (projection of snap onto its own normalised basis)
        let rec_norm: f64 = rec.iter().map(|x| x.norm_sqr()).sum::<f64>().sqrt();
        assert!((rec_norm - nrm).abs() < 1e-10 * nrm,
            "reconstruct norm mismatch: got {rec_norm}, expected {nrm}");
    }

    #[test]
    fn uniform_subset_endpoints() {
        let sub = uniform_subset(10, 3);
        assert_eq!(sub[0], 0);
        assert_eq!(*sub.last().unwrap(), 9);
        assert_eq!(sub.len(), 3);
    }

    #[test]
    fn rom_dependent_snapshots_reduced() {
        // Two identical snapshots → only 1 independent basis vector
        let n = 3;
        let snap = vec![Complex64::new(1.0, 0.0), Complex64::new(2.0, 0.0), Complex64::new(3.0, 0.0)];
        let rom = MomRom::build(
            vec![1e9, 2e9],
            vec![snap.clone(), snap],
            1e-10,
        ).unwrap();
        assert_eq!(rom.basis.ncols(), 1, "duplicate snapshots should yield 1 basis vector");
    }
}
