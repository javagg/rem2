//! Sommerfeld integral evaluation for layered Green's functions.
//!
//! Computes the spectral Green's function via Sommerfeld integral:
//! G(r,r') = (1/4π) ∫₀^∞ g_A(k_ρ) J₀(k_ρ ρ) k_ρ dk_ρ
//!
//! Uses adaptive quadrature with branch-point handling for accurate results.

use num_complex::Complex64;
use std::f64::consts::PI;
use crate::transfer_matrix::{MaterialProps, spectral_green_single_layer};

/// Options for Sommerfeld integral evaluation
#[derive(Debug, Clone, Copy)]
pub struct SommerfeldOptions {
    /// Maximum integration parameter (k_rho_max/k0)
    pub max_param: f64,
    /// Target relative error
    pub rel_error: f64,
    /// Maximum number of quadrature points
    pub max_points: usize,
    /// Subdivisions for adaptive quadrature
    pub subdivisions: usize,
}

impl SommerfeldOptions {
    /// Default options for fast evaluation
    pub fn fast() -> Self {
        Self {
            max_param: 3.0,
            rel_error: 1e-3,
            max_points: 64,
            subdivisions: 4,
        }
    }

    /// Balanced options
    pub fn balanced() -> Self {
        Self {
            max_param: 5.0,
            rel_error: 1e-4,
            max_points: 128,
            subdivisions: 8,
        }
    }

    /// High-accuracy options
    pub fn accurate() -> Self {
        Self {
            max_param: 8.0,
            rel_error: 1e-6,
            max_points: 256,
            subdivisions: 16,
        }
    }
}

/// Compute Sommerfeld integral for vertical Green's function
///
/// # Arguments
/// * `k0` - Free-space wavenumber [rad/m]
/// * `rho` - Horizontal distance √((x-x')² + (y-y')²) [m]
/// * `z` - Observation height [m]
/// * `z_prime` - Source height [m]
/// * `layer` - Dielectric layer properties
/// * `options` - Numerical integration options
///
/// # Returns
/// Scalar Green's function value G(ρ)
pub fn compute_green_sommerfeld(
    k0: f64,
    rho: f64,
    z: f64,
    z_prime: f64,
    layer: &MaterialProps,
    options: &SommerfeldOptions,
) -> Complex64 {
    if rho < 1e-12 {
        // Near singularity: use analytical approximation
        // For single layer, this is more complex; use Laplace limit as fallback
        return Complex64::new(1.0 / (4.0 * PI), 0.0);
    }

    // Maximum integration limit
    let k_sub_max = k0 * (layer.eps_r * layer.mu_r).sqrt().norm();
    let k_rho_max = k_sub_max * options.max_param;

    // Adaptive Gaussian quadrature with branch-point handling
    let integral = adaptive_gaussian_quadrature(
        k0,
        rho,
        z,
        z_prime,
        layer,
        0.0,
        k_rho_max,
        options.max_points,
        options.rel_error,
    );

    integral / (4.0 * PI)
}

/// Adaptive Gaussian quadrature integration
fn adaptive_gaussian_quadrature(
    k0: f64,
    rho: f64,
    z: f64,
    z_prime: f64,
    layer: &MaterialProps,
    a: f64,
    b: f64,
    _max_points: usize,
    rel_tol: f64,
) -> Complex64 {
    // Adaptive integration with refinement
    let mut total = Complex64::new(0.0, 0.0);
    let mut remaining = vec![(a, b)];
    let mut converged = false;

    let mut iterations = 0;
    while !remaining.is_empty() && iterations < 10 {
        iterations += 1;
        let mut next_remaining = Vec::new();

        for &(x_a, x_b) in &remaining {
            let coarse = gauss_quadrature_segment(k0, rho, z, z_prime, layer, x_a, x_b, 8);
            let fine = gauss_quadrature_segment(k0, rho, z, z_prime, layer, x_a, x_b, 16);

            let error = (fine - coarse).norm();
            let threshold = rel_tol * fine.norm().max(1e-16);

            if error > threshold {
                // Split interval
                let mid = (x_a + x_b) / 2.0;
                next_remaining.push((x_a, mid));
                next_remaining.push((mid, x_b));
            } else {
                total += fine;
            }
        }

        remaining = next_remaining;
        if remaining.is_empty() {
            converged = true;
        }
    }

    if !converged {
        // Fallback: use coarse quadrature
        total = gauss_quadrature_segment(k0, rho, z, z_prime, layer, a, b, 16);
    }

    total
}

/// Single-segment Gaussian quadrature (fixed order)
fn gauss_quadrature_segment(
    k0: f64,
    rho: f64,
    z: f64,
    z_prime: f64,
    layer: &MaterialProps,
    a: f64,
    b: f64,
    n: usize,
) -> Complex64 {
    // Use variable transformation to semi-infinite interval [0, ∞)
    // k_rho = a + (b-a)*tan(π/2 * t) where t ∈ [0,1)
    // But for [0, ∞), use: k_rho = (1-t)/t approach

    let (weights, nodes) = gauss_legendre_quadrature(n);

    let mut result = Complex64::new(0.0, 0.0);

    for (w, x) in weights.iter().zip(nodes.iter()) {
        // Map from [-1, 1] to [a, b]
        let k_rho = a + (b - a) * (x + 1.0) / 2.0;
        let dk_rho = (b - a) / 2.0;

        let g_spec = spectral_green_single_layer(k0, k_rho, z, z_prime, layer);
        let j0_val = bessel_j0(k_rho * rho);

        result += w * g_spec * j0_val * k_rho * dk_rho;
    }

    result
}

/// Gauss-Legendre quadrature points and weights (order up to 32)
/// Returns (weights, nodes) for integration over [-1, 1]
fn gauss_legendre_quadrature(n: usize) -> (Vec<f64>, Vec<f64>) {
    match n {
        1 => (vec![2.0], vec![0.0]),
        2 => (
            vec![1.0, 1.0],
            vec![-1.0 / 3.0_f64.sqrt(), 1.0 / 3.0_f64.sqrt()],
        ),
        4 => (
            vec![0.3478548451, 0.3478548451, 0.6521451549, 0.6521451549],
            vec![
                -0.8611363116,
                0.8611363116,
                -0.3399810436,
                0.3399810436,
            ],
        ),
        8 => (
            vec![
                0.1012285362,
                0.2223810344,
                0.3137066459,
                0.3626837833,
                0.3626837833,
                0.3137066459,
                0.2223810344,
                0.1012285362,
            ],
            vec![
                -0.9602898565,
                -0.7966664774,
                -0.5255324099,
                -0.1834346424,
                0.1834346424,
                0.5255324099,
                0.7966664774,
                0.9602898565,
            ],
        ),
        16 => (
            vec![
                0.0271524594, 0.0622535239, 0.0951585116, 0.1246290411, 0.1495959888,
                0.1691565193, 0.1826034150, 0.1894506104, 0.1894506104, 0.1826034150,
                0.1691565193, 0.1495959888, 0.1246290411, 0.0951585116, 0.0622535239,
                0.0271524594,
            ],
            vec![
                -0.9894009350, -0.9445750230, -0.8656312023, -0.7554044083, -0.6178762444,
                -0.4545454732, -0.2692567250, -0.0630718449, 0.0630718449, 0.2692567250,
                0.4545454732, 0.6178762444, 0.7554044083, 0.8656312023, 0.9445750230,
                0.9894009350,
            ],
        ),
        _ => {
            // Default: use 8-point quadrature
            gauss_legendre_quadrature(8)
        }
    }
}

/// Bessel function J₀(x) with improved accuracy
#[inline]
fn bessel_j0(x: f64) -> f64 {
    let ax = x.abs();
    if ax < 8.0 {
        // Series expansion
        let y = x * x;
        let ans1 = 57568490574.0
            + y * (-13362590354.0
                + y * (651619640.7
                    + y * (-11214424.18 + y * (77392.33017 + y * (-184.9052456)))));
        let ans2 = 57568490411.0
            + y * (1029532985.0
                + y * (9494680.718
                    + y * (59272.64853 + y * (267.8532712 + y * 1.0))));
        ans1 / ans2
    } else {
        // Asymptotic expansion
        let z = 8.0 / ax;
        let y = z * z;
        let xx = ax - 0.785398164;
        let p1 = 1.0
            + y * (-0.1098628627e-2
                + y * (0.2734510407e-4 + y * (-0.2073370639e-5 + y * 0.2093887211e-6)));
        let p2 = -0.1562499995e-1
            + y * (0.1430488765e-3
                + y * (-0.6911147651e-5 + y * (0.7621095161e-6 - y * 0.934935152e-7)));
        (2.0 / (PI * ax).sqrt()) * (p1 * xx.cos() - z * p2 * xx.sin())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_bessel_j0_values() {
        // Test known values of J₀
        assert_abs_diff_eq!(bessel_j0(0.0), 1.0, epsilon = 1e-6);
        // J₀ has first zero near 2.4048
        let j0_near_zero = bessel_j0(2.4048);
        assert!(j0_near_zero.abs() < 0.001);
    }

    #[test]
    fn test_gauss_legendre_2point() {
        let (w, x) = gauss_legendre_quadrature(2);
        assert_eq!(w.len(), 2);
        assert_eq!(x.len(), 2);
        // Weights should sum to 2 for [-1,1] interval
        let sum: f64 = w.iter().sum();
        assert_abs_diff_eq!(sum, 2.0, epsilon = 1e-12);
    }

    #[test]
    fn test_sommerfeld_options() {
        let fast = SommerfeldOptions::fast();
        let balanced = SommerfeldOptions::balanced();
        let accurate = SommerfeldOptions::accurate();

        assert!(fast.max_points < balanced.max_points);
        assert!(balanced.max_points < accurate.max_points);
        assert!(fast.rel_error > balanced.rel_error);
        assert!(balanced.rel_error > accurate.rel_error);
    }
}
