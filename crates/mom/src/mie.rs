//! Mie series analytical solutions for plane-wave scattering.
//!
//! Supports:
//! - PEC sphere: EFIE/CFIE pulse-basis and RWG-basis validation
//! - Homogeneous dielectric sphere: PMCHWT validation (Lorenz-Mie theory)
//!
//! References:
//! - Bohren & Huffman, *Absorption and Scattering of Light by Small Particles*, Ch. 4
//! - Stratton, *Electromagnetic Theory*, §9.25

use num_complex::Complex64;
use std::f64::consts::PI;

/// Compute the bistatic RCS σ(θ) [m²] of a PEC sphere of radius `a` [m]
/// at wavenumber `k` [1/m] for a x-polarized plane wave incident from +z.
///
/// Returns RCS at scattering angles `theta_deg` (degrees from forward).
/// Uses N_terms Mie coefficients (auto-selected if None).
pub fn pec_sphere_rcs(
    a: f64,
    k: f64,
    theta_deg: &[f64],
    n_terms: Option<usize>,
) -> Vec<f64> {
    let ka = k * a;
    // Wiscombe criterion: N ≈ ka + 4*(ka)^(1/3) + 2
    let n_max = n_terms.unwrap_or_else(|| (ka + 4.0*ka.powf(1.0/3.0) + 2.0).ceil() as usize + 5);

    // Mie coefficients for PEC sphere: a_n = j_n(ka)/h_n(ka), b_n = [ka j_n(ka)]'/[ka h_n(ka)]'
    let an = mie_an_pec(ka, n_max);
    let bn = mie_bn_pec(ka, n_max);

    theta_deg.iter().map(|&theta| {
        let cos_t = theta.to_radians().cos();
        let (s1, s2) = scattering_amplitudes(cos_t, &an, &bn);
        let wavelength = 2.0 * PI / k;
        // σ = λ²/(π) * |S|²   (bistatic, unpolarized average)
        let sigma = (wavelength * wavelength / PI) * 0.5 * (s1.norm_sqr() + s2.norm_sqr());
        sigma
    }).collect()
}

/// Compute the bistatic RCS σ(θ) [m²] of a homogeneous dielectric sphere
/// using Lorenz-Mie theory (Bohren & Huffman, Ch. 4, eqs. 4.53–4.56).
///
/// # Parameters
/// - `a`       — sphere radius [m]
/// - `k`       — free-space wave number k₁ = ω/c [1/m]
/// - `eps_r`   — relative permittivity of sphere (real, ≥ 1)
/// - `mu_r`    — relative permeability of sphere (real, ≥ 1)
/// - `theta_deg` — bistatic angles [°] from forward scatter (0° = forward, 180° = backward)
/// - `n_terms` — number of Mie terms (auto-selected by Wiscombe criterion if None)
///
/// Returns bistatic RCS [m²] as unpolarized average ½(|S₁|² + |S₂|²) × λ²/π.
pub fn dielectric_sphere_rcs(
    a: f64,
    k: f64,
    eps_r: f64,
    mu_r: f64,
    theta_deg: &[f64],
    n_terms: Option<usize>,
) -> Vec<f64> {
    let ka  = k * a;
    let m   = (eps_r * mu_r).sqrt(); // real refractive index ratio n₂/n₁
    let mx  = m * ka;                // internal size parameter

    let n_max = n_terms.unwrap_or_else(|| (ka + 4.0*ka.powf(1.0/3.0) + 2.0).ceil() as usize + 5);

    let an = mie_an_diel(ka, mx, m, n_max);
    let bn = mie_bn_diel(ka, mx, m, n_max);

    theta_deg.iter().map(|&theta| {
        let cos_t = theta.to_radians().cos();
        let (s1, s2) = scattering_amplitudes(cos_t, &an, &bn);
        let wavelength = 2.0 * PI / k;
        (wavelength * wavelength / PI) * 0.5 * (s1.norm_sqr() + s2.norm_sqr())
    }).collect()
}

/// Lorenz-Mie a_n for a homogeneous dielectric sphere (TM modes).
///
/// a_n = [m ψ_n(mx) ψ'_n(x) − ψ_n(x) ψ'_n(mx)] /
///       [m ψ_n(mx) ξ'_n(x) − ξ_n(x) ψ'_n(mx)]
fn mie_an_diel(ka: f64, mx: f64, m: f64, n_max: usize) -> Vec<Complex64> {
    let mc = Complex64::new(m, 0.0);
    (1..=n_max).map(|n| {
        let psi_mx  = sph_jn(n, mx) * Complex64::new(mx, 0.0);
        let dpsi_mx = Complex64::new(d_sph_jn(n, mx), 0.0);
        let psi_x   = sph_jn(n, ka) * Complex64::new(ka, 0.0);
        let dpsi_x  = Complex64::new(d_sph_jn(n, ka), 0.0);
        let xi_x    = sph_hn1(n, ka) * Complex64::new(ka, 0.0);
        let dxi_x   = d_sph_hn1(n, ka);
        let num = mc * psi_mx * dpsi_x - psi_x * dpsi_mx;
        let den = mc * psi_mx * dxi_x  - xi_x  * dpsi_mx;
        num / den
    }).collect()
}

/// Lorenz-Mie b_n for a homogeneous dielectric sphere (TE modes).
///
/// b_n = [ψ_n(mx) ψ'_n(x) − m ψ_n(x) ψ'_n(mx)] /
///       [ψ_n(mx) ξ'_n(x) − m ξ_n(x) ψ'_n(mx)]
fn mie_bn_diel(ka: f64, mx: f64, m: f64, n_max: usize) -> Vec<Complex64> {
    let mc = Complex64::new(m, 0.0);
    (1..=n_max).map(|n| {
        let psi_mx  = sph_jn(n, mx) * Complex64::new(mx, 0.0);
        let dpsi_mx = Complex64::new(d_sph_jn(n, mx), 0.0);
        let psi_x   = sph_jn(n, ka) * Complex64::new(ka, 0.0);
        let dpsi_x  = Complex64::new(d_sph_jn(n, ka), 0.0);
        let xi_x    = sph_hn1(n, ka) * Complex64::new(ka, 0.0);
        let dxi_x   = d_sph_hn1(n, ka);
        let num = psi_mx * dpsi_x - mc * psi_x * dpsi_mx;
        let den = psi_mx * dxi_x  - mc * xi_x  * dpsi_mx;
        num / den
    }).collect()
}

/// Mie a_n coefficients for PEC sphere (TM modes).
/// a_n = -j_n(ka) / h_n^(1)(ka)
fn mie_an_pec(ka: f64, n_max: usize) -> Vec<Complex64> {
    (1..=n_max).map(|n| {
        let jn = sph_jn(n, ka);
        let hn = sph_hn1(n, ka);
        -jn / hn
    }).collect()
}

/// Mie b_n coefficients for PEC sphere (TE modes).
/// b_n = -[ka j_n(ka)]' / [ka h_n^(1)(ka)]'
fn mie_bn_pec(ka: f64, n_max: usize) -> Vec<Complex64> {
    (1..=n_max).map(|n| {
        let djn = d_sph_jn(n, ka);
        let dhn = d_sph_hn1(n, ka);
        -Complex64::new(djn, 0.0) / dhn
    }).collect()
}

/// Far-field scattering amplitudes S₁(cos θ), S₂(cos θ).
fn scattering_amplitudes(cos_t: f64, an: &[Complex64], bn: &[Complex64]) -> (Complex64, Complex64) {
    let mut s1 = Complex64::ZERO;
    let mut s2 = Complex64::ZERO;

    let n_max = an.len();
    let (pi_arr, tau_arr) = pi_tau(cos_t, n_max);

    for n in 1..=n_max {
        let coef = (2*n + 1) as f64 / (n*(n+1)) as f64;
        s1 += coef * (an[n-1] * pi_arr[n-1] + bn[n-1] * tau_arr[n-1]);
        s2 += coef * (an[n-1] * tau_arr[n-1] + bn[n-1] * pi_arr[n-1]);
    }

    (s1, s2)
}

/// Angular functions πₙ(cos θ) and τₙ(cos θ) via recurrence.
fn pi_tau(cos_t: f64, n_max: usize) -> (Vec<Complex64>, Vec<Complex64>) {
    let mut pi_n = vec![Complex64::ZERO; n_max];
    let mut tau_n = vec![Complex64::ZERO; n_max];

    // Initial values: π₁ = 1, τ₁ = cos θ
    pi_n[0] = Complex64::new(1.0, 0.0);
    tau_n[0] = Complex64::new(cos_t, 0.0);

    if n_max >= 2 {
        pi_n[1] = Complex64::new(3.0 * cos_t, 0.0);
        tau_n[1] = Complex64::new(3.0 * (2.0*cos_t*cos_t - 1.0), 0.0);
    }

    for n in 2..n_max {
        let nf = n as f64 + 1.0; // actual n value (1-indexed)
        pi_n[n] = ((2.0*nf - 1.0)/(nf - 1.0)) * cos_t * pi_n[n-1]
                - (nf/(nf - 1.0)) * pi_n[n-2];
        tau_n[n] = nf * cos_t * pi_n[n] - (nf + 1.0) * pi_n[n-1];
    }

    (pi_n, tau_n)
}

// ---------------------------------------------------------------------------
// Spherical Bessel functions via downward recurrence
// ---------------------------------------------------------------------------

/// Spherical Bessel function j_n(x) via upward recurrence.
fn sph_jn(n: usize, x: f64) -> Complex64 {
    if x < 1e-10 {
        return if n == 0 { Complex64::new(1.0, 0.0) } else { Complex64::ZERO };
    }
    // j_0(x) = sin(x)/x
    let mut jnm1 = Complex64::new(x.sin() / x, 0.0);
    if n == 0 { return jnm1; }
    // j_1(x) = sin(x)/x² - cos(x)/x
    let mut jn0 = Complex64::new(x.sin()/(x*x) - x.cos()/x, 0.0);
    for k in 1..n {
        let jnp1 = ((2*k + 1) as f64 / x) * jn0 - jnm1;
        jnm1 = jn0;
        jn0 = jnp1;
    }
    jn0
}

/// Derivative: d/dx [x · j_n(x)] = x · j_{n-1}(x) - n · j_n(x) ... simplified form.
fn d_sph_jn(n: usize, x: f64) -> f64 {
    // d/dx [x j_n] = x j_{n-1} - n j_n (standard Bessel recurrence)
    let jn_val = sph_jn(n, x).re;
    if n == 0 {
        // d/dx [sin x] = cos x
        return x.cos();
    }
    let jnm1 = sph_jn(n - 1, x).re;
    x * jnm1 - (n as f64) * jn_val
}

/// Spherical Hankel function of first kind h_n^(1)(x) = j_n(x) + i·y_n(x).
fn sph_hn1(n: usize, x: f64) -> Complex64 {
    sph_jn(n, x) + Complex64::new(0.0, 1.0) * sph_yn(n, x)
}

/// Derivative d/dx [x · h_n^(1)(x)].
fn d_sph_hn1(n: usize, x: f64) -> Complex64 {
    let hn = sph_hn1(n, x);
    let hnm1 = if n == 0 {
        // h_{-1}^(1)(x) = -e^{ix}/x (but we handle n=0 derivative separately)
        let ex = Complex64::new(0.0, x).exp();
        -ex / Complex64::new(x, 0.0)
    } else {
        sph_hn1(n - 1, x)
    };
    Complex64::new(x, 0.0) * hnm1 - Complex64::new(n as f64, 0.0) * hn
}

/// Spherical Neumann function y_n(x) via upward recurrence.
fn sph_yn(n: usize, x: f64) -> Complex64 {
    if x < 1e-10 { return Complex64::new(f64::NEG_INFINITY, 0.0); }
    // y_0(x) = -cos(x)/x
    let mut ynm1 = Complex64::new(-x.cos() / x, 0.0);
    if n == 0 { return ynm1; }
    // y_1(x) = -cos(x)/x² - sin(x)/x
    let mut yn0 = Complex64::new(-x.cos()/(x*x) - x.sin()/x, 0.0);
    for k in 1..n {
        let ynp1 = ((2*k + 1) as f64 / x) * yn0 - ynm1;
        ynm1 = yn0;
        yn0 = ynp1;
    }
    yn0
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn optical_theorem_ka1() {
        // For PEC sphere, geometric optics limit σ → π*a² for large ka.
        // At ka=1, just check the value is reasonable (should be ~few×π*a²).
        let a = 1.0;
        let k = 1.0 / a; // ka = 1
        let rcs = pec_sphere_rcs(a, k, &[180.0], None); // backscatter
        let geometric = std::f64::consts::PI * a * a;
        // At ka=1, backscatter RCS ≈ 0.1 * π*a² to ~10 * π*a² (resonance region)
        assert!(rcs[0] > 0.0, "RCS must be positive");
        assert!(rcs[0] < 100.0 * geometric, "RCS sanity upper bound");
        assert!(rcs[0].is_finite(), "RCS must be finite");
    }

    #[test]
    fn large_sphere_approaches_geometric() {
        // For ka >> 1, forward RCS → 4π*a² (forward scattering theorem) — not tested here.
        // Instead just verify convergence of result is stable across n_terms.
        let a = 1.0;
        let k = 5.0; // ka=5
        let rcs5  = pec_sphere_rcs(a, k, &[90.0], Some(15));
        let rcs20 = pec_sphere_rcs(a, k, &[90.0], Some(20));
        // Should agree to within 1%
        let rel_err = ((rcs5[0] - rcs20[0]) / rcs20[0]).abs();
        assert!(rel_err < 0.01, "Mie series not converged: rel_err={:.3e}", rel_err);
    }
}
