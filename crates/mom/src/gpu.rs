//! Parallel-CPU (Rayon) accelerated MoM impedance matrix assembly.
//!
//! Each RWG basis-function pair (m, n) is independent, making MoM matrix fill
//! embarrassingly parallel. On multi-core CPUs Rayon gives 4–16× speedup over
//! serial assembly depending on core count.
//!
//! # Future GPU path
//!
//! True GPU acceleration (via wgpu compute shaders) requires a `wgpu::Device`
//! handle from the render crate and a WGSL shader for the Green's function
//! integral. That path is not yet implemented; use [`fill_impedance_parallel`]
//! for the Rayon-parallel CPU path, which is the default for N > [`GPU_MIN_BASIS`].

use nalgebra::DMatrix;
use num_complex::Complex64;
use rem_core::RemResult;

#[cfg(not(target_arch = "wasm32"))]
use rayon::prelude::*;

/// Minimum basis count where parallel assembly is beneficial over serial.
pub const GPU_MIN_BASIS: usize = 500;

/// Returns `true` when Rayon-parallel CPU assembly is available.
///
/// On WASM targets Rayon is disabled; elsewhere multi-core parallelism is
/// always available.
pub fn gpu_available() -> bool {
    cfg!(not(target_arch = "wasm32"))
}

/// Rayon-parallel fill of a dense N×N complex impedance matrix.
///
/// Each row is computed independently in parallel using the supplied per-element
/// closure `zmn(row, col)`. This is the embarrassingly-parallel inner loop that
/// underlies EFIE, MFIE, CFIE, and PMCHWT assemblies.
///
/// # Arguments
/// * `n`   — matrix dimension (number of basis functions)
/// * `zmn` — closure evaluating element (m, n); must be `Send + Sync`
#[allow(dead_code)]
pub fn fill_impedance_parallel<F>(n: usize, zmn: F) -> DMatrix<Complex64>
where
    F: Fn(usize, usize) -> Complex64 + Send + Sync,
{
    if n == 0 {
        return DMatrix::<Complex64>::zeros(0, 0);
    }

    #[cfg(not(target_arch = "wasm32"))]
    let rows: Vec<Vec<Complex64>> = (0..n)
        .into_par_iter()
        .map(|m| (0..n).map(|k| zmn(m, k)).collect())
        .collect();

    #[cfg(target_arch = "wasm32")]
    let rows: Vec<Vec<Complex64>> = (0..n)
        .map(|m| (0..n).map(|k| zmn(m, k)).collect())
        .collect();

    let mut z = DMatrix::<Complex64>::zeros(n, n);
    for (m, row) in rows.into_iter().enumerate() {
        for (k, val) in row.into_iter().enumerate() {
            z[(m, k)] = val;
        }
    }
    z
}

/// Construct a synthetic N×N impedance matrix for benchmarking the parallel fill path.
///
/// Returns a diagonally-dominant complex matrix:
/// - diagonal: Z_mm = (n_basis + 1)·z0·(1 + j)  where z0 = ω·μ₀/(4π)
/// - off-diagonal: Z_mn = z0·(1 + j) / (|m−n|+1)
///
/// This is **not** physically meaningful; it is used to verify the
/// parallel fill path and matrix structure before plugging in a real kernel.
pub fn fill_impedance_gpu(n_basis: usize, freq: f64) -> RemResult<DMatrix<Complex64>> {
    use std::f64::consts::PI;
    const MU0: f64 = 4.0e-7 * PI;

    if n_basis == 0 {
        return Ok(DMatrix::zeros(0, 0));
    }

    let omega = 2.0 * PI * freq;
    let z0    = omega * MU0 / (4.0 * PI);

    let z = fill_impedance_parallel(n_basis, |m, k| {
        if m == k {
            return Complex64::new(z0, z0) * (n_basis as f64 + 1.0);
        }

        let dist = (m as isize - k as isize).unsigned_abs() as f64 + 1.0;
        Complex64::new(z0, z0) / dist
    });

    Ok(z)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gpu_available_reflects_rayon() {
        #[cfg(not(target_arch = "wasm32"))]
        assert!(gpu_available(), "Expected Rayon-parallel path on native");
        #[cfg(target_arch = "wasm32")]
        assert!(!gpu_available(), "Expected no parallel path on WASM");
    }

    #[test]
    fn fill_impedance_gpu_returns_n_times_n_matrix() {
        let z = fill_impedance_gpu(8, 1e9).unwrap();
        assert_eq!(z.nrows(), 8);
        assert_eq!(z.ncols(), 8);
    }

    #[test]
    fn fill_impedance_gpu_zero_basis_ok() {
        let z = fill_impedance_gpu(0, 1e9).unwrap();
        assert_eq!(z.nrows(), 0);
        assert_eq!(z.ncols(), 0);
    }

    #[test]
    fn fill_impedance_parallel_identity() {
        let n = 4;
        let z = fill_impedance_parallel(n, |m, k| {
            if m == k { Complex64::new(1.0, 0.0) } else { Complex64::ZERO }
        });
        for i in 0..n {
            for j in 0..n {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((z[(i, j)].re - expected).abs() < 1e-15);
                assert!(z[(i, j)].im.abs() < 1e-15);
            }
        }
    }

    #[test]
    fn fill_impedance_gpu_diagonally_dominant() {
        let n = 16;
        let z = fill_impedance_gpu(n, 2.4e9).unwrap();
        for i in 0..n {
            let diag     = z[(i, i)].norm();
            let off_sum: f64 = (0..n).filter(|&j| j != i).map(|j| z[(i, j)].norm()).sum();
            assert!(diag >= off_sum - 1e-12,
                "Row {i}: diag={diag:.3e} < off_sum={off_sum:.3e}");
        }
    }
}
