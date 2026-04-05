//! Fresnel reflection and transmission coefficients for SBR+.
//!
//! Handles PEC (perfect conductor) and general dielectric interfaces.
//! All angles in radians; complex permittivity/permeability supported.

use num_complex::Complex64;
use rem_core::ETA0;
use crate::ray::{cross3, dot3, sub3, scale3, normalize3};

// ---------------------------------------------------------------------------
// Material interface
// ---------------------------------------------------------------------------

/// Material parameters at the interface (medium 2 = transmitted side).
#[derive(Debug, Clone)]
pub struct Interface {
    /// Relative permittivity of medium 2 (may be complex for lossy media)
    pub eps_r: Complex64,
    /// Relative permeability of medium 2
    pub mu_r: Complex64,
    /// Is this a PEC boundary? (overrides eps/mu)
    pub is_pec: bool,
}

impl Interface {
    /// Perfect electric conductor
    pub fn pec() -> Self {
        Self { eps_r: Complex64::new(1.0, 0.0), mu_r: Complex64::new(1.0, 0.0), is_pec: true }
    }

    /// Lossless dielectric
    pub fn dielectric(eps_r: f64, mu_r: f64) -> Self {
        Self { eps_r: Complex64::new(eps_r, 0.0), mu_r: Complex64::new(mu_r, 0.0), is_pec: false }
    }

    /// Lossy dielectric (complex permittivity)
    pub fn lossy(eps_r: Complex64, mu_r: Complex64) -> Self {
        Self { eps_r, mu_r, is_pec: false }
    }

    /// Wave impedance η₂ = η₀ √(μᵣ / εᵣ)
    pub fn eta(&self) -> Complex64 {
        let eta0 = Complex64::new(ETA0, 0.0);
        if self.is_pec { Complex64::ZERO } else { eta0 * (self.mu_r / self.eps_r).sqrt() }
    }
}

// ---------------------------------------------------------------------------
// Fresnel coefficients (TE and TM)
// ---------------------------------------------------------------------------

/// Fresnel reflection coefficients for a planar interface.
///
/// `cos_theta_i` = cosine of incidence angle (> 0 by convention).
/// Medium 1 = free space (η₁ = η₀).
///
/// Returns (Γ_TE, Γ_TM).
pub fn fresnel_refl(iface: &Interface, cos_theta_i: f64) -> (Complex64, Complex64) {
    if iface.is_pec {
        // PEC: total reflection with phase reversal
        return (Complex64::new(-1.0, 0.0), Complex64::new(1.0, 0.0));
    }

    let eta1 = Complex64::new(ETA0, 0.0);
    let eta2 = iface.eta();
    let n_rel_sq = iface.eps_r * iface.mu_r; // (n₂/n₁)²

    // Snell's law: cos θₜ via n₁² sin²θᵢ + cos²θₜ = 1 + (n₁/n₂)²... use:
    //   cos θₜ = sqrt(1 - sin²θₜ) = sqrt(1 - (1 - cos²θᵢ)/n_rel²)
    let sin2_i = 1.0 - cos_theta_i * cos_theta_i;
    let cos_theta_t = (Complex64::new(1.0, 0.0) - Complex64::new(sin2_i, 0.0) / n_rel_sq).sqrt();

    let ci = Complex64::new(cos_theta_i, 0.0);

    // TE: Γ = (η₂ cosθᵢ − η₁ cosθₜ) / (η₂ cosθᵢ + η₁ cosθₜ)
    let gamma_te = (eta2 * ci - eta1 * cos_theta_t) / (eta2 * ci + eta1 * cos_theta_t);

    // TM: Γ = (η₁ cosθᵢ − η₂ cosθₜ) / (η₁ cosθᵢ + η₂ cosθₜ)
    let gamma_tm = (eta1 * ci - eta2 * cos_theta_t) / (eta1 * ci + eta2 * cos_theta_t);

    (gamma_te, gamma_tm)
}

// ---------------------------------------------------------------------------
// Reflected field update
// ---------------------------------------------------------------------------

/// Compute the reflected electric field after a Fresnel bounce.
///
/// `e_inc` = incident E-field at the hit point
/// `h_inc` = incident H-field at the hit point
/// `dir_inc` = unit incident ray direction
/// `normal` = outward surface unit normal (pointing away from surface into incident medium)
///
/// Returns the reflected E-field vector and the new ray direction.
pub fn reflect_field(
    e_inc: &[Complex64; 3],
    _h_inc: &[Complex64; 3],
    dir_inc: &[f64; 3],
    normal: &[f64; 3],
    iface: &Interface,
) -> ([Complex64; 3], [f64; 3]) {
    // Ensure normal points toward incident side
    let n = if dot3(dir_inc, normal) < 0.0 {
        *normal
    } else {
        [-normal[0], -normal[1], -normal[2]]
    };

    let cos_i = (-dot3(dir_inc, &n)).max(0.0);
    let (gamma_te, gamma_tm) = fresnel_refl(iface, cos_i);

    // Decompose E_inc into TE (s) and TM (p) components
    // TE polarisation: ê_s = k̂ × n̂ / |k̂ × n̂|
    let k_cross_n = cross3(dir_inc, &n);
    let k_cross_n_len = dot3(&k_cross_n, &k_cross_n).sqrt();

    let e_refl = if k_cross_n_len < 1e-12 {
        // Normal incidence: both polarisations degenerate → Γ_TM applies
        scale_complex3(e_inc, gamma_tm)
    } else {
        let es_hat: [f64; 3] = normalize3(k_cross_n); // TE unit vector
        // ê_p = ê_s × k̂ (TM unit vector in incident plane)
        let ep_hat: [f64; 3] = normalize3(cross3(&es_hat, dir_inc));

        // Project E onto TE and TM
        let e_te: Complex64 = dot3_complex(e_inc, &es_hat);
        let e_tm: Complex64 = dot3_complex(e_inc, &ep_hat);

        // Reflected TM direction: ê_p_refl = ê_s × k̂_refl
        let new_dir = mirror3(dir_inc, &n);
        let ep_hat_refl = normalize3(cross3(&es_hat, &new_dir));

        // Reflected E = Γ_TE * E_TE * ê_s + Γ_TM * E_TM * ê_p_refl
        let mut e_r = [Complex64::ZERO; 3];
        for i in 0..3 {
            e_r[i] = gamma_te * e_te * es_hat[i] + gamma_tm * e_tm * ep_hat_refl[i];
        }
        e_r
    };

    let new_dir = mirror3(dir_inc, &n);
    (e_refl, new_dir)
}

/// Mirror a direction vector across a normal.
#[inline]
pub fn mirror3(dir: &[f64; 3], n: &[f64; 3]) -> [f64; 3] {
    let d2n = 2.0 * dot3(dir, n);
    normalize3(sub3(dir, &scale3(n, d2n)))
}

// ---------------------------------------------------------------------------
// PO surface current from H_inc (PEC only)
// ---------------------------------------------------------------------------

/// Physical Optics surface current on a PEC surface:
/// `J_PO = 2 n̂ × H_inc`
///
/// Returns `[0;3]` if the face is in shadow (n̂ · (−k̂) < 0).
pub fn po_current_pec(
    h_inc: &[Complex64; 3],
    normal: &[f64; 3],
    dir_inc: &[f64; 3],
) -> [Complex64; 3] {
    // Illumination check: face must face the incoming ray
    if dot3(normal, dir_inc) >= 0.0 {
        return [Complex64::ZERO; 3]; // facing away → shadow
    }
    cross3_complex(normal, h_inc, 2.0)
}

// ---------------------------------------------------------------------------
// Vector helpers for complex fields
// ---------------------------------------------------------------------------

#[inline]
fn dot3_complex(a: &[Complex64; 3], b: &[f64; 3]) -> Complex64 {
    a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
}

#[inline]
fn scale_complex3(v: &[Complex64; 3], s: Complex64) -> [Complex64; 3] {
    [v[0]*s, v[1]*s, v[2]*s]
}

/// `s * (n × v)` where n is real and v is complex.
#[inline]
fn cross3_complex(n: &[f64; 3], v: &[Complex64; 3], s: f64) -> [Complex64; 3] {
    let sc = Complex64::new(s, 0.0);
    [
        sc * (n[1]*v[2] - n[2]*v[1]),
        sc * (n[2]*v[0] - n[0]*v[2]),
        sc * (n[0]*v[1] - n[1]*v[0]),
    ]
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn pec_te_reflection_is_minus_one() {
        let (g_te, _) = fresnel_refl(&Interface::pec(), 0.5);
        assert_abs_diff_eq!(g_te.re, -1.0, epsilon = 1e-12);
        assert_abs_diff_eq!(g_te.im,  0.0, epsilon = 1e-12);
    }

    #[test]
    fn normal_incidence_on_dielectric() {
        // At normal incidence (cos_theta = 1) for ε=4 medium: Γ_TM = (1-2)/(1+2) = -1/3
        let iface = Interface::dielectric(4.0, 1.0);
        let (_, g_tm) = fresnel_refl(&iface, 1.0);
        // η2 = η0/2, Γ_TM = (η0 - η0/2)/(η0 + η0/2) = (1/2)/(3/2) = 1/3
        // Our sign convention: (η1 cosθi - η2 cosθt)/(η1 cosθi + η2 cosθt)
        // = (η0 - η0/2)/(η0 + η0/2) = 1/3
        assert!((g_tm.re - 1.0/3.0).abs() < 1e-6, "got {}", g_tm);
    }

    #[test]
    fn mirror_direction() {
        let d = normalize3([1.0, 0.0, -1.0]);
        let n = [0.0, 0.0, 1.0];
        let r = mirror3(&d, &n);
        // Should reflect z-component
        assert!((r[0] - d[0]).abs() < 1e-12);
        assert!((r[2] + d[2]).abs() < 1e-12);
    }

    #[test]
    fn po_current_zero_in_shadow() {
        let h_inc = [Complex64::new(1.0, 0.0); 3];
        let normal  = [0.0, 0.0, 1.0];
        let dir_inc = [0.0, 0.0, 1.0]; // same direction as normal → shadow
        let j = po_current_pec(&h_inc, &normal, &dir_inc);
        assert!(j.iter().all(|c| c.norm() < 1e-14));
    }

    #[test]
    fn po_current_nonzero_illuminated() {
        let h_inc = [Complex64::new(1.0, 0.0), Complex64::ZERO, Complex64::ZERO];
        let normal  = [0.0, 0.0, 1.0];
        let dir_inc = [0.0, 0.0, -1.0]; // incoming from +z
        let j = po_current_pec(&h_inc, &normal, &dir_inc);
        // J = 2 n̂ × H = 2 ẑ × x̂ = 2 ŷ
        assert!((j[1].re - 2.0).abs() < 1e-12, "J_y = {}", j[1]);
    }
}
