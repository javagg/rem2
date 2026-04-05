//! Excitation vectors for MoM (plane wave, lumped port, etc.)

use crate::surface_mesh::SurfaceMesh;
use num_complex::Complex64;

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Parameters for an incident plane wave.
///
/// Propagation direction: k̂ = (sin θ cos φ, sin θ sin φ, cos θ)
/// where θ (theta_inc) is the polar angle from +z and φ (phi_inc) from +x.
///
/// Polarization `pol` selects the electric field direction:
/// - `"theta"` (default): E ∝ θ̂ (theta unit vector, H-pol in xz-plane)
/// - `"phi"`:             E ∝ φ̂ (phi unit vector)
/// - `"x"`, `"y"`, `"z"`: Cartesian polarization (projected to transverse)
#[derive(Debug, Clone)]
pub struct PlaneWave {
    /// Polar angle of incidence [radians] (0 = +z direction)
    pub theta_inc: f64,
    /// Azimuth angle of incidence [radians] (0 = +x)
    pub phi_inc: f64,
    /// Polarization: "theta" | "phi" | "x" | "y" | "z"
    pub pol: String,
}

impl Default for PlaneWave {
    fn default() -> Self {
        // Default: +z incidence, x-polarised (matches original implementation)
        PlaneWave { theta_inc: 0.0, phi_inc: 0.0, pol: "x".to_string() }
    }
}

impl PlaneWave {
    /// Propagation unit vector k̂.
    pub fn k_hat(&self) -> [f64; 3] {
        let (st, ct) = (self.theta_inc.sin(), self.theta_inc.cos());
        let (sp, cp) = (self.phi_inc.sin(), self.phi_inc.cos());
        [st*cp, st*sp, ct]
    }

    /// Electric field polarisation unit vector ê.
    /// Transverse to k̂.
    pub fn e_hat(&self) -> [f64; 3] {
        let (st, ct) = (self.theta_inc.sin(), self.theta_inc.cos());
        let (sp, cp) = (self.phi_inc.sin(), self.phi_inc.cos());
        match self.pol.to_lowercase().as_str() {
            "theta" => [ct*cp, ct*sp, -st],   // θ̂
            "phi"   => [-sp, cp, 0.0],          // φ̂
            "x"     => {
                // x̂ projected transverse to k̂, then normalized
                let kh = self.k_hat();
                let dot = kh[0]; // k̂·x̂
                let e = [1.0 - dot*kh[0], -dot*kh[1], -dot*kh[2]];
                normalize3(e)
            }
            "y"     => {
                let kh = self.k_hat();
                let dot = kh[1];
                let e = [-dot*kh[0], 1.0 - dot*kh[1], -dot*kh[2]];
                normalize3(e)
            }
            "z"     => {
                let kh = self.k_hat();
                let dot = kh[2];
                let e = [-dot*kh[0], -dot*kh[1], 1.0 - dot*kh[2]];
                normalize3(e)
            }
            _ => [ct*cp, ct*sp, -st], // fallback: theta
        }
    }
}

fn normalize3(v: [f64; 3]) -> [f64; 3] {
    let len = (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]).sqrt();
    if len < 1e-14 { [1.0, 0.0, 0.0] } else { [v[0]/len, v[1]/len, v[2]/len] }
}

// ---------------------------------------------------------------------------
// Pulse basis RHS
// ---------------------------------------------------------------------------

/// Build RHS for pulse basis using a general PlaneWave.
///
/// V[m] = -∫_Tm ê · E_inc(r) dS ≈ -ê · E_inc(centroid_m) * area_m
///
/// E_inc(r) = ê * exp(-jk k̂·r)
pub fn plane_wave_rhs_general(surf: &SurfaceMesh, k: f64, wave: &PlaneWave, basis: &str) -> Vec<Complex64> {
    match basis.to_lowercase().as_str() {
        "pulse" => pulse_rhs(surf, k, wave),
        _       => rwg_rhs(surf, k, wave),
    }
}

/// Backward-compatible wrapper: +z incidence, x-polarised.
pub fn plane_wave_rhs(surf: &SurfaceMesh, k: f64, basis: &str) -> Vec<Complex64> {
    plane_wave_rhs_general(surf, k, &PlaneWave::default(), basis)
}

fn pulse_rhs(surf: &SurfaceMesh, k: f64, wave: &PlaneWave) -> Vec<Complex64> {
    let kh = wave.k_hat();
    let eh = wave.e_hat();
    surf.faces.iter().map(|face| {
        let r = &face.centroid;
        // Phase: exp(-jk k̂·r)
        let phase = k * (kh[0]*r[0] + kh[1]*r[1] + kh[2]*r[2]);
        let e_inc = Complex64::new(0.0, -phase).exp(); // e^{-jk k̂·r}
        // Galerkin test with pulse basis (constant = 1 on face):
        // V[m] = -∫ (ê · E_inc) f_m dS ≈ -(ê · ê * e_inc) * area
        // For pulse basis the testing function dotted with E_inc gives scalar
        // We project E_inc onto the face normal to get the "equivalent testing"
        // but for EFIE the testing is with the basis itself (ê direction).
        // Simplest: V[m] = -e_inc * face.area (x-component coupling, eh[0])
        // General: project onto tangential plane via n̂ × (E × n̂) = E_tan
        let e_x = eh[0] * e_inc;
        let e_y = eh[1] * e_inc;
        let e_z = eh[2] * e_inc;
        // Tangential component along the face (ê · ê = 1, already transverse to k̂)
        // For scalar pulse basis: use dominant Cartesian component
        // Full treatment: -∫ f_m · E_inc dS, f_m is the x̂ direction (assumed)
        // To generalise: dot ê with the "effective" testing direction.
        // For now we project E onto the face tangential plane using the normal:
        let nn = &face.normal;
        let ndote = nn[0]*e_x.re + nn[1]*e_y.re + nn[2]*e_z.re;
        let et_x = e_x - Complex64::new(ndote * nn[0], 0.0);
        let et_y = e_y - Complex64::new(ndote * nn[1], 0.0);
        let et_z = e_z - Complex64::new(ndote * nn[2], 0.0);
        // Scalar projection onto ê direction (tangential ê):
        let val = eh[0]*et_x + eh[1]*et_y + eh[2]*et_z;
        -val * face.area
    }).collect()
}

fn rwg_rhs(surf: &SurfaceMesh, k: f64, wave: &PlaneWave) -> Vec<Complex64> {
    use crate::basis::rwg::generate_rwg_bases;
    let kh = wave.k_hat();
    let eh = wave.e_hat();
    let bases = generate_rwg_bases(surf);

    bases.iter().map(|b| {
        let mut val = Complex64::ZERO;
        for &(face_idx, in_plus) in &[(b.plus_face, true), (b.minus_face, false)] {
            let face = &surf.faces[face_idx];
            let r = &face.centroid;
            let phase = k * (kh[0]*r[0] + kh[1]*r[1] + kh[2]*r[2]);
            let e_inc_scalar = Complex64::new(0.0, -phase).exp();
            let fn_ = b.eval(r, surf, in_plus);
            // f_n · E_inc = f_n · (ê * e_inc_scalar)
            let dot = eh[0]*fn_[0] + eh[1]*fn_[1] + eh[2]*fn_[2];
            val += Complex64::new(dot, 0.0) * e_inc_scalar * face.area;
        }
        -val
    }).collect()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_wave_khat_is_plus_z() {
        let w = PlaneWave::default();
        let k = w.k_hat();
        assert!((k[0]).abs() < 1e-14);
        assert!((k[1]).abs() < 1e-14);
        assert!((k[2] - 1.0).abs() < 1e-14);
    }

    #[test]
    fn theta90_phi0_wave_khat_is_plus_x() {
        let w = PlaneWave { theta_inc: std::f64::consts::FRAC_PI_2, phi_inc: 0.0, pol: "z".to_string() };
        let k = w.k_hat();
        assert!((k[0] - 1.0).abs() < 1e-14, "kx={}", k[0]);
        assert!((k[1]).abs() < 1e-14);
        assert!((k[2]).abs() < 1e-14);
    }

    #[test]
    fn e_hat_orthogonal_to_k_hat() {
        for (ti, pi) in [(0.0f64, 0.0f64), (1.0, 0.5), (0.7, 1.2)] {
            for pol in ["theta", "phi", "x", "y", "z"] {
                let w = PlaneWave { theta_inc: ti, phi_inc: pi, pol: pol.to_string() };
                let k = w.k_hat();
                let e = w.e_hat();
                let dot = k[0]*e[0] + k[1]*e[1] + k[2]*e[2];
                assert!(dot.abs() < 1e-13,
                    "k̂·ê = {:.2e} not zero for θ={ti},φ={pi},pol={pol}", dot);
            }
        }
    }
}
