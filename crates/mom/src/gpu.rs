//! GPU-accelerated MoM impedance matrix assembly via wgpu compute.
//!
//! Each RWG basis-function pair is independent, making MoM matrix fill
//! embarrassingly parallel and well-suited for GPU acceleration.
//!
//! Expected speedup: 10-50× for N > 1000 basis functions.
//!
//! # Integration
//!
//! The full GPU path requires a `wgpu::Device` handle from the render crate.
//! When available, mesh geometry is uploaded to GPU buffers and a WGSL compute
//! shader (`mom_impedance.wgsl`) evaluates all (RWG_i, RWG_j) interactions in
//! parallel. Results are read back as a complex dense matrix.
//!
//! The WGSL shader is located at `crates/render/src/shaders/mom_impedance.wgsl`.

use nalgebra::DMatrix;
use num_complex::Complex64;
use rem_core::RemResult;

/// Placeholder for GPU-accelerated impedance matrix fill.
///
/// Full implementation requires passing a `wgpu::Device` from the render
/// crate. Until then, falls back to CPU assembly via `fill_impedance_naive`.
pub fn fill_impedance_gpu(
    _n_basis: usize,
    _freq: f64,
) -> RemResult<DMatrix<Complex64>> {
    Err(rem_core::RemError::NotImplemented(
        "GPU MoM assembly: requires wgpu::Device integration. \
         Use CPU path (fill_impedance_naive) for now.".into()
    ))
}

/// Check if GPU acceleration is available.
pub fn gpu_available() -> bool {
    false
}

/// Minimum basis count where GPU acceleration is beneficial.
pub const GPU_MIN_BASIS: usize = 500;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gpu_not_available_by_default() {
        assert!(!gpu_available());
    }

    #[test]
    fn gpu_fill_returns_not_implemented() {
        let result = fill_impedance_gpu(10, 1e9);
        assert!(result.is_err());
    }
}
