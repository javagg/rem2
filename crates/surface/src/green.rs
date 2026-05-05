//! Free-space Green's function and its derivatives.

use num_complex::Complex64;
use std::f64::consts::PI;

/// 3D free-space scalar Green's function:
/// G(r, r') = exp(-jkR) / (4πR),  R = |r - r'|
///
/// Returns `Complex64::ZERO` if R < `eps` (singular point — handled by Duffy transform).
#[inline]
pub fn green3d(r: &[f64; 3], r_prime: &[f64; 3], k: f64) -> Complex64 {
    let dist = dist3(r, r_prime);
    if dist < 1e-14 { return Complex64::ZERO; }
    let phase = Complex64::new(0.0, -k * dist);
    phase.exp() / (4.0 * PI * dist)
}

/// Normal derivative of G w.r.t. r': ∂G/∂n' = G · (jkR + 1)/R² · (r-r')·n'
#[inline]
pub fn green3d_normal_deriv(
    r: &[f64; 3],
    r_prime: &[f64; 3],
    n_prime: &[f64; 3],
    k: f64,
) -> Complex64 {
    let dist = dist3(r, r_prime);
    if dist < 1e-14 { return Complex64::ZERO; }
    let g = green3d(r, r_prime, k);
    let dot = (r[0]-r_prime[0])*n_prime[0]
            + (r[1]-r_prime[1])*n_prime[1]
            + (r[2]-r_prime[2])*n_prime[2];
    let factor = Complex64::new(1.0, k * dist) / (dist * dist);
    g * factor * dot
}

/// Static (k=0) Laplace Green function: G_L(r,r') = 1/(4πR)
#[inline]
pub fn green_laplace(r: &[f64; 3], r_prime: &[f64; 3]) -> f64 {
    let dist = dist3(r, r_prime);
    if dist < 1e-14 { return 0.0; }
    1.0 / (4.0 * PI * dist)
}

#[inline]
fn dist3(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    let dx = a[0]-b[0]; let dy = a[1]-b[1]; let dz = a[2]-b[2];
    (dx*dx + dy*dy + dz*dz).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn green3d_magnitude_falls_off() {
        // At distance R, |G| = 1/(4πR)
        let r  = [1.0, 0.0, 0.0];
        let rp = [0.0, 0.0, 0.0];
        let k  = 0.1;
        let g  = green3d(&r, &rp, k);
        let expected_mag = 1.0 / (4.0 * std::f64::consts::PI * 1.0);
        assert_abs_diff_eq!(g.norm(), expected_mag, epsilon = 1e-14);
    }

    #[test]
    fn green3d_zero_at_singular() {
        let r = [0.0, 0.0, 0.0];
        assert_eq!(green3d(&r, &r, 1.0), Complex64::ZERO);
    }

    #[test]
    fn green_laplace_matches_static_limit() {
        // At k→0, |G(k)| → G_L
        let r  = [0.5, 0.0, 0.0];
        let rp = [0.0, 0.0, 0.0];
        let g_wave = green3d(&r, &rp, 0.0).re;
        let g_lap  = green_laplace(&r, &rp);
        assert_abs_diff_eq!(g_wave, g_lap, epsilon = 1e-14);
    }
}
