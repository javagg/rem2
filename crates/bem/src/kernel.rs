//! Laplace and Helmholtz BEM kernels.
//!
//! ## Laplace kernels
//!
//! Single-layer potential (V kernel):
//!   G(r,r') = 1 / (4π|r-r'|)
//!
//! Double-layer potential (K kernel):
//!   ∂G/∂n'(r,r') = -1/(4π) * (r-r')·n̂' / |r-r'|³
//!
//! ## Regularisation for near-singular points
//!
//! When |r-r'| < tol, kernels return 0.0 (caller uses Duffy/analytic treatment).

use std::f64::consts::PI;

/// Laplace single-layer kernel: G(r, r') = 1/(4π R).
#[inline]
#[allow(non_snake_case)]
pub fn laplace_G(r: &[f64; 3], rp: &[f64; 3]) -> f64 {
    let rx = r[0]-rp[0]; let ry = r[1]-rp[1]; let rz = r[2]-rp[2];
    let dist = (rx*rx + ry*ry + rz*rz).sqrt();
    if dist < 1e-14 { return 0.0; }
    1.0 / (4.0 * PI * dist)
}

/// Laplace double-layer kernel: ∂G/∂n'(r, r') = -(r-r')·n̂' / (4π R³).
///
/// `n_prime`: outward unit normal at source point r'.
#[inline]
#[allow(non_snake_case)]
pub fn laplace_dG_dn(r: &[f64; 3], rp: &[f64; 3], n_prime: &[f64; 3]) -> f64 {
    let rx = r[0]-rp[0]; let ry = r[1]-rp[1]; let rz = r[2]-rp[2];
    let dist2 = rx*rx + ry*ry + rz*rz;
    let dist = dist2.sqrt();
    if dist < 1e-14 { return 0.0; }
    let dot = rx*n_prime[0] + ry*n_prime[1] + rz*n_prime[2];
    -dot / (4.0 * PI * dist2 * dist)
}

/// Laplace hypersingular kernel T(r, r') = ∂G/∂n (derivative wrt observer normal).
///
/// `n`: outward normal at observer point r.
#[inline]
#[allow(non_snake_case)]
pub fn laplace_dG_dn_obs(r: &[f64; 3], rp: &[f64; 3], n: &[f64; 3]) -> f64 {
    let rx = r[0]-rp[0]; let ry = r[1]-rp[1]; let rz = r[2]-rp[2];
    let dist2 = rx*rx + ry*ry + rz*rz;
    let dist = dist2.sqrt();
    if dist < 1e-14 { return 0.0; }
    let dot = rx*n[0] + ry*n[1] + rz*n[2];
    dot / (4.0 * PI * dist2 * dist)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn laplace_G_falls_off_as_1_over_r() {
        let r  = [1.0, 0.0, 0.0];
        let rp = [0.0; 3];
        let g = laplace_G(&r, &rp);
        let expected = 1.0 / (4.0 * PI);
        assert!((g - expected).abs() < 1e-14, "G(1,0,0) = {g}, expected {expected}");
    }

    #[test]
    fn laplace_G_symmetric() {
        let r  = [1.0, 2.0, 3.0];
        let rp = [0.5, 0.1, 0.7];
        assert!((laplace_G(&r, &rp) - laplace_G(&rp, &r)).abs() < 1e-14);
    }

    #[test]
    fn laplace_dG_dn_antisymmetric_normals() {
        let r  = [1.0, 0.0, 0.0];
        let rp = [0.0; 3];
        let n1 = [1.0, 0.0, 0.0];
        let n2 = [-1.0, 0.0, 0.0];
        let d1 = laplace_dG_dn(&r, &rp, &n1);
        let d2 = laplace_dG_dn(&r, &rp, &n2);
        assert!((d1 + d2).abs() < 1e-14, "dG/dn not antisymmetric under n flip: {d1} + {d2}");
    }

    #[test]
    fn laplace_G_zero_at_coincident() {
        let r = [1.0, 2.0, 3.0];
        assert_eq!(laplace_G(&r, &r), 0.0);
        assert_eq!(laplace_dG_dn(&r, &r, &[0.0, 0.0, 1.0]), 0.0);
    }
}
