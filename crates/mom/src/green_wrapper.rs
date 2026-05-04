//! Wrapper to adapt rem-layered-green's GreenFunction trait for use in MoM assembly.
//! 
//! This module bridges the existing assemble.rs functions with the new trait-based
//! GreenFunction interface, enabling support for layered-media Green's functions.

use rem_layered_green::{GreenFunction, FreeSpaceGreen, LayeredGreen, DielectricLayer};
use std::f64::consts::PI;
use rem_core::C0;

/// Create a free-space Green's function from frequency.
pub fn free_space_from_freq(freq: f64) -> FreeSpaceGreen {
    FreeSpaceGreen::from_freq(freq)
}

/// Wrapper: obtain a `&dyn GreenFunction` trait object from a free-space wavenumber.
///
/// This is used internally by assemble.rs to keep the interface extensible without
/// immediately breaking existing code.
///
/// Returns a boxed trait object suitable for use in matrix assembly.
pub fn box_free_space_from_freq(freq: f64) -> Box<dyn GreenFunction> {
    Box::new(FreeSpaceGreen::from_freq(freq))
}

/// Shorthand: create free-space GreenFunction from k (wavenumber).
pub fn box_free_space_from_k(k: f64) -> Box<dyn GreenFunction> {
    let freq = k * C0 / (2.0 * PI);
    Box::new(FreeSpaceGreen::from_freq(freq))
}

/// Create a layered-media Green's function from frequency and substrate definition.
///
/// # Arguments
/// * `freq` - Frequency [Hz]
/// * `layers` - Vector of dielectric layers (from bottom to top)
///
/// # Returns
/// Boxed LayeredGreen trait object
pub fn box_layered_green_from_freq(freq: f64, layers: Vec<DielectricLayer>) -> Box<dyn GreenFunction> {
    let k0 = 2.0 * PI * freq / C0;
    Box::new(LayeredGreen::new(layers, k0))
}

/// Create a single-layer substrate over PEC-like region.
///
/// Common use: microstrip patch antenna or substrate-backed antenna.
/// Returns a simple single-layer DielectricLayer for testing.
pub fn single_layer_fr4() -> DielectricLayer {
    // FR-4 dielectric: εᵣ ≈ 4.2, tan(δ) ≈ 0.02, μᵣ = 1
    // Thickness ≈ 1 mm typical
    DielectricLayer {
        eps_r: 4.2,
        loss_tan: 0.02,
        mu_r: 1.0,
        thickness_m: 1.0e-3,
        eps_r_complex_override: None,
        eps_r_z: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wrapper_creates_green_function() {
        let green = free_space_from_freq(1.0e9);
        let r = [1.0, 0.0, 0.0];
        let rp = [0.0, 0.0, 0.0];
        let g = green.g(&r, &rp);
        assert!(g.norm() > 0.0);
    }

    #[test]
    fn layered_green_creation() {
        let layers = vec![single_layer_fr4()];
        let green = box_layered_green_from_freq(1.0e9, layers);
        let r = [1.0, 0.0, 0.0];
        let rp = [0.0, 0.0, 0.0];
        let g = green.g(&r, &rp);
        // Should produce a valid complex number
        assert!(!g.re.is_nan());
        assert!(!g.im.is_nan());
    }

    #[test]
    fn free_space_vs_layered_comparison() {
        // At moderate frequency with free-space-like substrate,
        // layered should produce non-zero values
        let layers = vec![DielectricLayer {
            eps_r: 1.0,
            loss_tan: 0.0,
            mu_r: 1.0,
            thickness_m: 1e10, // Very thick "air"
            eps_r_complex_override: None,
            eps_r_z: None,
        }];
        let freq = 1.0e9; // 1 GHz
        let r = [0.1, 0.0, 0.0]; // 10 cm
        let rp = [0.0, 0.0, 0.0];
        
        let free_space = free_space_from_freq(freq);
        let layered = box_layered_green_from_freq(freq, layers);
        
        let g_free = free_space.g(&r, &rp);
        let g_lay = layered.g(&r, &rp);
        
        // Free-space should definitely be non-zero
        assert!(g_free.norm() > 0.0, "Free-space Green's function should be non-zero");
        
        // Layered result should be valid (no NaN/Inf)
        assert!(!g_lay.re.is_nan());
        assert!(!g_lay.im.is_nan());
        assert!(!g_lay.re.is_infinite());
        assert!(!g_lay.im.is_infinite());
    }
}
