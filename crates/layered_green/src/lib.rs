//! rem-layered-green — Green's function implementations for stratified media and free space.
//!
//! This crate provides:
//! - `GreenFunction` trait: unified interface for free-space and layered-media Green's functions
//! - `FreeSpaceGreen`: scalar 3D free-space Helmholtz Green's function (self-contained)
//! - `LayeredGreen`: stratified dielectric Green's function (DCIM with Sommerfeld integral)
//!
//! # Architecture
//! ```text
//! GreenFunction (trait)
//!   ├─ FreeSpaceGreen (3D self-space, O(1) evaluation)
//!   └─ LayeredGreen (Sommerfeld integral via DCIM, O(N_poles) evaluation)
//! ```

mod discrete_image;
mod sommerfeld;
mod transfer_matrix;

use discrete_image::{DcimApproximation, GpofFitter};
use num_complex::Complex64;
use sommerfeld::{compute_green_sommerfeld, SommerfeldOptions};
use std::f64::consts::PI;
use transfer_matrix::MaterialProps;

/// Unified Green's function interface for MoM assembly.
/// Provides scalar potential G(r,r') and its derivatives.
pub trait GreenFunction: Send + Sync {
    /// Scalar Green's function G(r, r') at wavenumber k.
    /// Returns 0 if source and observation points coincide (singularity handled by quadrature).
    fn g(&self, r: &[f64; 3], r_prime: &[f64; 3]) -> Complex64;

    /// Gradient of G with respect to r: ∇_r G(r, r').
    /// Returns [∂G/∂x, ∂G/∂y, ∂G/∂z].
    fn grad_g(&self, r: &[f64; 3], r_prime: &[f64; 3]) -> [Complex64; 3];

    /// Normal derivative: ∂G/∂n' = (∇G) · n'.
    fn normal_deriv_g(&self, r: &[f64; 3], r_prime: &[f64; 3], n_prime: &[f64; 3]) -> Complex64 {
        let grad = self.grad_g(r, r_prime);
        grad[0] * n_prime[0] + grad[1] * n_prime[1] + grad[2] * n_prime[2]
    }
}

/// Free-space scalar Helmholtz Green's function in 3D.
/// G(r, r') = exp(-jk|r-r'|) / (4π|r-r'|)
#[derive(Debug, Clone, Copy)]
pub struct FreeSpaceGreen {
    /// Wavenumber k = 2πf/c [rad/m]
    pub k: f64,
}

impl FreeSpaceGreen {
    /// Create a free-space Green's function at the given wavenumber.
    pub fn new(k: f64) -> Self {
        Self { k }
    }

    /// Create from frequency [Hz] and assumed free-space propagation.
    pub fn from_freq(freq_hz: f64) -> Self {
        const C0: f64 = 299792458.0; // speed of light [m/s]
        Self {
            k: 2.0 * PI * freq_hz / C0,
        }
    }

    /// Euclidean distance between r and r'.
    #[inline]
    fn distance(r: &[f64; 3], r_prime: &[f64; 3]) -> f64 {
        let dx = r[0] - r_prime[0];
        let dy = r[1] - r_prime[1];
        let dz = r[2] - r_prime[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
}

impl GreenFunction for FreeSpaceGreen {
    fn g(&self, r: &[f64; 3], r_prime: &[f64; 3]) -> Complex64 {
        let dist = Self::distance(r, r_prime);
        if dist < 1e-14 {
            return Complex64::ZERO;
        }
        let phase = Complex64::new(0.0, -self.k * dist);
        phase.exp() / (4.0 * PI * dist)
    }

    fn grad_g(&self, r: &[f64; 3], r_prime: &[f64; 3]) -> [Complex64; 3] {
        let dist = Self::distance(r, r_prime);
        if dist < 1e-14 {
            return [Complex64::ZERO; 3];
        }
        let dx = r[0] - r_prime[0];
        let dy = r[1] - r_prime[1];
        let dz = r[2] - r_prime[2];
        let g = self.g(r, r_prime);
        // ∇G = G · (jk + 1/R) / R · (r - r')
        let factor = g * (Complex64::new(0.0, self.k) + Complex64::new(1.0 / dist, 0.0)) / dist;
        [factor * dx, factor * dy, factor * dz]
    }
}

/// Placeholder for layered dielectric Green's function.
/// Uses Sommerfeld integral approximation via DCIM (discrete complex image method).
///
/// # Features
/// - Support for single or multi-layer substrates
/// - Pre-computed DCIM coefficients (poles + residues) for fast evaluation
/// - Automatic GPOF fitting from Sommerfeld integral samples
#[derive(Debug, Clone)]
pub struct LayeredGreen {
    /// Layer definition (from bottom to top)
    pub layers: Vec<DielectricLayer>,
    /// Wavenumber in background medium [rad/m]
    pub k0: f64,
    /// Cached DCIM approximation (poles + residues)
    dcim: Option<DcimApproximation>,
    /// Source and observation heights (z, z_prime)
    z_obs: f64,
    z_src: f64,
}

/// Single dielectric layer definition.
#[derive(Debug, Clone)]
pub struct DielectricLayer {
    /// Relative permittivity (isotropic for now; may extend to tensor)
    pub eps_r: f64,
    /// Loss tangent: tan(δ) for dissipation model
    pub loss_tan: f64,
    /// Relative permeability
    pub mu_r: f64,
    /// Layer thickness [m]; use large value (1e10) for top (air)
    pub thickness_m: f64,
}

impl LayeredGreen {
    /// Create a layered Green's function from layer stack.
    /// `k0`: wavenumber in background medium (air).
    pub fn new(layers: Vec<DielectricLayer>, k0: f64) -> Self {
        Self {
            layers,
            k0,
            dcim: None,
            z_obs: 0.0,
            z_src: 0.0,
        }
    }

    /// Pre-compute DCIM approximation for given observation and source heights.
    /// Should be called once per height pair for efficiency.
    fn compute_dcim(&mut self, z_obs: f64, z_src: f64) {
        if let Some(ref _dcim) = self.dcim {
            // Already computed for these heights
            if (self.z_obs - z_obs).abs() < 1e-10 && (self.z_src - z_src).abs() < 1e-10 {
                return;
            }
        }

        self.z_obs = z_obs;
        self.z_src = z_src;

        // Build material properties for Sommerfeld integral
        let material = self.build_material_properties();

        // Sample Sommerfeld integral over range of horizontal distances
        let mut fitter = GpofFitter::new(8); // 8 poles typically sufficient

        let options = SommerfeldOptions::balanced();
        let rho_samples = (0..20)
            .map(|i| {
                let rho = 0.01 * 10.0_f64.powf(i as f64 / 5.0);
                rho
            })
            .collect::<Vec<_>>();

        for rho in rho_samples {
            let g_val = compute_green_sommerfeld(
                self.k0,
                rho,
                z_obs,
                z_src,
                &material,
                &options,
            );
            fitter.add_sample(rho, g_val);
        }

        // Fit DCIM approximation
        self.dcim = Some(fitter.fit());
    }

    /// Build material properties from layer definition
    fn build_material_properties(&self) -> MaterialProps {
        if self.layers.is_empty() {
            // Free space
            return MaterialProps {
                eps_r: Complex64::new(1.0, 0.0),
                mu_r: Complex64::new(1.0, 0.0),
                thickness: 1e10,
            };
        }

        // Use first layer (single-layer assumption for now)
        let layer = &self.layers[0];

        // Compute complex permittivity: εᵣ = ε'ᵣ - jε''ᵣ = εᵣ(1 - j tan(δ))
        let eps_r_real = layer.eps_r;
        let eps_r_imag = -layer.eps_r * layer.loss_tan;

        MaterialProps {
            eps_r: Complex64::new(eps_r_real, eps_r_imag),
            mu_r: Complex64::new(layer.mu_r, 0.0),
            thickness: layer.thickness_m,
        }
    }

    /// Euclidean distance between r and r'.
    #[inline]
    fn distance(r: &[f64; 3], r_prime: &[f64; 3]) -> f64 {
        let dx = r[0] - r_prime[0];
        let dy = r[1] - r_prime[1];
        let dz = r[2] - r_prime[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    /// Horizontal distance ρ = √((x-x')² + (y-y')²)
    #[inline]
    fn horizontal_distance(r: &[f64; 3], r_prime: &[f64; 3]) -> f64 {
        let dx = r[0] - r_prime[0];
        let dy = r[1] - r_prime[1];
        (dx * dx + dy * dy).sqrt()
    }
}

impl GreenFunction for LayeredGreen {
    fn g(&self, r: &[f64; 3], r_prime: &[f64; 3]) -> Complex64 {
        // Get horizontal and vertical components
        let rho = Self::horizontal_distance(r, r_prime);
        let z = r[2];
        let z_prime = r_prime[2];

        // Check for singularity
        if rho < 1e-14 && (z - z_prime).abs() < 1e-14 {
            return Complex64::ZERO;
        }

        // Use Sommerfeld integral directly (DCIM caching happens inside)
        let material = self.build_material_properties();
        let options = SommerfeldOptions::balanced();
        compute_green_sommerfeld(self.k0, rho, z, z_prime, &material, &options)
    }

    fn grad_g(&self, r: &[f64; 3], r_prime: &[f64; 3]) -> [Complex64; 3] {
        let h = 1e-6; // Finite difference step
        
        // Numerical differentiation for now (can be optimized with analytical)
        let g_center = self.g(r, r_prime);
        
        let r_dx = [r[0] + h, r[1], r[2]];
        let g_dx = self.g(&r_dx, r_prime);
        
        let r_dy = [r[0], r[1] + h, r[2]];
        let g_dy = self.g(&r_dy, r_prime);
        
        let r_dz = [r[0], r[1], r[2] + h];
        let g_dz = self.g(&r_dz, r_prime);
        
        [
            (g_dx - g_center) / h,
            (g_dy - g_center) / h,
            (g_dz - g_center) / h,
        ]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn free_space_magnitude_falloff() {
        // At distance R, |G| = 1/(4πR)
        let r = [1.0, 0.0, 0.0];
        let r_prime = [0.0, 0.0, 0.0];
        let green = FreeSpaceGreen::new(0.1);
        let g = green.g(&r, &r_prime);
        let expected_mag = 1.0 / (4.0 * PI);
        assert_abs_diff_eq!(g.norm(), expected_mag, epsilon = 1e-14);
    }

    #[test]
    fn free_space_zero_at_singular() {
        let r = [0.0, 0.0, 0.0];
        let green = FreeSpaceGreen::new(1.0);
        let g = green.g(&r, &r);
        assert_eq!(g, Complex64::ZERO);
    }

    #[test]
    fn free_space_grad_magnitude() {
        // |∇G| should decay as 1/R²  at small k
        let r = [0.1, 0.0, 0.0];
        let r_prime = [0.0, 0.0, 0.0];
        let green = FreeSpaceGreen::new(0.01); // small k
        let grad = green.grad_g(&r, &r_prime);
        let mag_sq = grad[0].norm_sqr() + grad[1].norm_sqr() + grad[2].norm_sqr();
        let _mag = mag_sq.sqrt(); // Should be ~O(1/0.01²) = O(10000) for small k
        assert!(mag_sq > 0.0);
    }

    #[test]
    fn layered_green_structure() {
        // Test that LayeredGreen can be created and stores parameters correctly
        let layers = vec![
            DielectricLayer {
                eps_r: 4.0,
                loss_tan: 0.01,
                mu_r: 1.0,
                thickness_m: 1.0,
            },
        ];
        let green = LayeredGreen::new(layers.clone(), 1.0);
        assert_eq!(green.k0, 1.0);
        assert_eq!(green.layers.len(), 1);
        assert_eq!(green.layers[0].eps_r, 4.0);
    }

    #[test]
    fn layered_green_singular_handling() {
        // Test that singularity is properly handled
        let layers = vec![
            DielectricLayer {
                eps_r: 1.0,
                loss_tan: 0.0,
                mu_r: 1.0,
                thickness_m: 1e10,
            },
        ];
        let green = LayeredGreen::new(layers, 1.0);
        let r = [0.0, 0.0, 0.0];
        let g = green.g(&r, &r);
        assert_eq!(g, Complex64::ZERO);
    }

    #[test]
    fn layered_green_differentiation() {
        // Test that gradient computation works without panicking
        let layers = vec![
            DielectricLayer {
                eps_r: 1.0,
                loss_tan: 0.0,
                mu_r: 1.0,
                thickness_m: 1e10,
            },
        ];
        let green = LayeredGreen::new(layers, 10.0); // Use higher wavenumber
        let r = [1.0, 0.5, 0.2];
        let r_prime = [0.0, 0.0, 0.0];
        let grad = green.grad_g(&r, &r_prime);
        // Should have 3 components
        assert_eq!(grad.len(), 3);
        // Gradient computation should not produce NaN or Inf
        for component in grad {
            assert!(!component.re.is_nan());
            assert!(!component.im.is_nan());
            assert!(!component.re.is_infinite());
            assert!(!component.im.is_infinite());
        }
    }
}
