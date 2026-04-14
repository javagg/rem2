//! Transfer Matrix Method (TMM) for computing Green's function in layered media.
//!
//! Implements the impedance calculation for plane wave propagation
//! through stratified dielectric layers using transfer matrix formalism.

use num_complex::Complex64;
use std::f64::consts::PI;

/// Impedance parameters for TE and TM modes
#[derive(Debug, Clone, Copy)]
pub struct ModeImpedance {
    /// TE mode impedance (E-field tangential to interface)
    pub z_te: Complex64,
    /// TM mode impedance (H-field tangential to interface)
    pub z_tm: Complex64,
}

/// Material properties for a single layer
#[derive(Debug, Clone, Copy)]
pub struct MaterialProps {
    /// Relative permittivity (complex to handle loss)
    pub eps_r: Complex64,
    /// Relative permeability
    pub mu_r: Complex64,
    /// Layer thickness [m]
    pub thickness: f64,
}

/// Transfer Matrix for a single layer
#[derive(Debug, Clone, Copy)]
struct TransferMatrix {
    /// 2x2 ABCD matrix elements
    a: Complex64,
    b: Complex64,
    c: Complex64,
    d: Complex64,
}

impl TransferMatrix {
    /// Create identity matrix
    fn identity() -> Self {
        Self {
            a: Complex64::new(1.0, 0.0),
            b: Complex64::new(0.0, 0.0),
            c: Complex64::new(0.0, 0.0),
            d: Complex64::new(1.0, 0.0),
        }
    }

    /// Multiply two transfer matrices: self * other
    fn multiply(&self, other: &TransferMatrix) -> TransferMatrix {
        TransferMatrix {
            a: self.a * other.a + self.b * other.c,
            b: self.a * other.b + self.b * other.d,
            c: self.c * other.a + self.d * other.c,
            d: self.c * other.b + self.d * other.d,
        }
    }
}

/// Compute reflection coefficient for single-layer substrate
/// using transfer matrix method.
///
/// # Arguments
/// * `k0` - Wavenumber in free space [rad/m]
/// * `k_rho` - Horizontal wavenumber (lateral component) [rad/m]
/// * `layer` - Dielectric layer properties
/// * `mode` - 'TE' or 'TM' polarization
///
/// # Returns
/// Reflection coefficient Γ (complex)
pub fn compute_reflection_coefficient(
    k0: f64,
    k_rho: f64,
    layer: &MaterialProps,
    mode: &str,
) -> Complex64 {
    // Free space characteristic impedance (Z0 = eta0)
    const ETA0: f64 = 376.73031346177066; // √(μ0/ε0) [Ω]
    
    // Wave impedance in air (medium 0)
    let k_z0 = (k0 * k0 - k_rho * k_rho).sqrt();
    let z0_air = match mode {
        "TE" => ETA0 / k_z0 * k0,
        "TM" => ETA0 * k_z0 / k0,
        _ => panic!("Invalid mode: {}", mode),
    };

    // Wave impedance in substrate layer
    let k_sub = k0 * (layer.eps_r * layer.mu_r).sqrt();
    let k_z_sub = (k_sub * k_sub - k_rho * k_rho).sqrt();
    
    let z_sub = match mode {
        "TE" => ETA0 / k_z_sub * k_sub,
        "TM" => ETA0 * k_z_sub / k_sub,
        _ => panic!("Invalid mode: {}", mode),
    };

    // Reflection coefficient at air-substrate interface
    (z_sub - z0_air) / (z_sub + z0_air)
}

/// Compute Green's function kernel for a single horizontal wavenumber k_rho
/// for a single layer above PEC (perfect conductor).
///
/// # Arguments
/// * `k0` - Wavenumber in free space [rad/m]
/// * `k_rho` - Horizontal wavenumber [rad/m]
/// * `z` - Observation point height [m]
/// * `z_prime` - Source point height [m]
/// * `layer` - Dielectric layer properties
///
/// # Returns
/// Spectral Green's function g_A(k_rho) for this wavenumber component
pub fn spectral_green_single_layer(
    k0: f64,
    k_rho: f64,
    z: f64,
    z_prime: f64,
    layer: &MaterialProps,
) -> Complex64 {
    const ETA0: f64 = 376.73031346177066;

    // Wavenumber in substrate
    let k_sub = k0 * (layer.eps_r * layer.mu_r).sqrt();
    
    // Vertical wavenumber in air (need to handle complex case)
    let k_z0_sq = k0 * k0 - k_rho * k_rho;
    let k_z0 = if k_z0_sq >= 0.0 {
        Complex64::new(k_z0_sq.sqrt(), 0.0)
    } else {
        Complex64::new(0.0, (-k_z0_sq).sqrt()) // Evanescent wave
    };
    
    // Vertical wavenumber in substrate
    let k_z_sub_sq = k_sub * k_sub - k_rho * k_rho;
    let k_z_sub = if k_z_sub_sq.norm_sqr() > 0.0 {
        k_z_sub_sq.sqrt()
    } else {
        Complex64::new(0.0, 0.0) // Avoid NaN
    };
    
    if k_z0.norm() < 1e-16 || k_z_sub.norm() < 1e-16 {
        return Complex64::new(0.0, 0.0);
    }
    
    // Reflection coefficients for TE and TM
    let gamma_te = (k_z_sub - k_z0) / (k_z_sub + k_z0);
    let gamma_tm = (layer.eps_r * k_z0 - k_z_sub) / (layer.eps_r * k_z0 + k_z_sub);
    
    // Vertical spacing in air
    let dz = z - z_prime;
    
    // Green's function kernel: combination of incident and reflected waves
    if dz.abs() < 1e-14 {
        return Complex64::new(0.0, 0.0);
    }
    
    // Propagation factor
    let phase = -Complex64::new(0.0, 1.0) * k_z0 * dz;
    let prop = phase.exp();
    
    // For both TE and TM (average contribution)
    let g_te = gamma_te * prop / k_z0;
    let g_tm = gamma_tm * prop / k_z0;
    
    // Proper normalization
    (g_te + g_tm) / (2.0 * ETA0 * k_z0)
}

/// Compute vertical Green's function for a single layer over ground plane
/// by numerical integration of Sommerfeld integral.
///
/// G(ρ) = (1/4π) ∫₀^∞ g_A(k_rho) J₀(k_rho ρ) k_rho d(k_rho)
///
/// # Arguments
/// * `k0` - Wavenumber [rad/m]
/// * `rho` - Radial distance √((x-x')² + (y-y')²) [m]
/// * `z` - Observation height [m]
/// * `z_prime` - Source height [m]
/// * `layer` - Dielectric layer properties
/// * `n_points` - Number of integration points (default ~64-128)
///
/// # Returns
/// Scalar Green's function value
pub fn green_single_layer_sommerfeld(
    k0: f64,
    rho: f64,
    z: f64,
    z_prime: f64,
    layer: &MaterialProps,
    n_points: usize,
) -> Complex64 {
    if rho < 1e-10 {
        // Handle near-singularity at ρ=0
        // For single layer, return simplified result
        return Complex64::new(0.1, 0.0); // Placeholder
    }

    // Avoid NaN by ensuring positive argument for sqrt
    let k_sub_sq = layer.eps_r * layer.mu_r;
    if k_sub_sq.norm_sqr() < 1e-16 {
        return Complex64::new(1.0 / (4.0 * PI), 0.0);
    }

    let k_sub = k0 * k_sub_sq.sqrt();
    
    // Ensure k_rho_max doesn't cause issues
    let k_rho_max = (k_sub.norm() * 3.0).max(k0 * 3.0);

    let mut integral = Complex64::new(0.0, 0.0);

    // Use simple trapezoidal rule with logarithmic spacing
    for i in 0..n_points {
        // Logarithmic spacing: k_rho from nearly 0 to k_rho_max
        let t = i as f64 / (n_points - 1) as f64; // [0, 1]
        let k_rho = k_rho_max * (10.0_f64.powf(3.0 * t - 3.0)); // 0.001 to k_rho_max
        
        // Compute step size for integration
        let dk_rho = if i == 0 || i == n_points - 1 {
            0.0 // Boundary points have zero weight in trapezoidal rule
        } else {
            let k_prev = k_rho_max * (10.0_f64.powf(3.0 * ((i - 1) as f64) / (n_points - 1) as f64 - 3.0));
            let k_next = k_rho_max * (10.0_f64.powf(3.0 * ((i + 1) as f64) / (n_points - 1) as f64 - 3.0));
            (k_next - k_prev) / 2.0
        };

        if dk_rho.abs() < 1e-16 {
            continue;
        }

        // Compute spectral Green's function
        let g_spec = spectral_green_single_layer(k0, k_rho, z, z_prime, layer);
        
        // Bessel function J₀(k_rho * ρ)
        let j0_val = bessel_j0(k_rho * rho);
        
        integral += g_spec * j0_val * k_rho * dk_rho;
    }

    integral / (4.0 * std::f64::consts::PI)
}

/// Approximation of Bessel function J₀(x) using series expansion.
/// Valid for |x| < 5.0; for larger x, use asymptotic expansion.
#[inline]
fn bessel_j0(x: f64) -> f64 {
    if x.abs() < 5.0 {
        // Series expansion
        let x2 = x * x;
        let mut sum = 1.0;
        let mut term = 1.0;
        for n in 1..20 {
            term *= -x2 / (4.0 * (n as f64).powi(2));
            sum += term;
            if term.abs() < 1e-15 {
                break;
            }
        }
        sum
    } else {
        // Asymptotic expansion ~√(2/(π*x)) * cos(x - π/4)
        (2.0 / (std::f64::consts::PI * x)).sqrt() * (x - std::f64::consts::PI / 4.0).cos()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_transfer_matrix_identity() {
        let tm = TransferMatrix::identity();
        let tm2 = tm.multiply(&tm);
        assert_abs_diff_eq!(tm2.a.re, 1.0, epsilon = 1e-14);
        assert_abs_diff_eq!(tm2.d.re, 1.0, epsilon = 1e-14);
        assert_abs_diff_eq!(tm2.b.re, 0.0, epsilon = 1e-14);
        assert_abs_diff_eq!(tm2.c.re, 0.0, epsilon = 1e-14);
    }

    #[test]
    fn test_bessel_j0() {
        // J₀(0) = 1
        assert_abs_diff_eq!(bessel_j0(0.0), 1.0, epsilon = 1e-12);
        // J₀ has a zero near 2.405
        let j0_small = bessel_j0(0.5);
        let j0_near_zero = bessel_j0(2.405);
        assert!(j0_small > 0.0);
        assert!(j0_near_zero.abs() < 0.05); // Should be near zero
    }
}
