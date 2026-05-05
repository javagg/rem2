//! Conductor surface roughness correction models for planar transmission lines.
//!
//! Real printed-circuit-board and MMIC conductors have rough surfaces that
//! increase the effective conductor resistance at high frequencies beyond the
//! smooth-surface skin-effect prediction.  This module provides two industry-
//! standard correction-factor models:
//!
//! # Hammerstad–Jensen (HJ) Model
//!
//! Hammerstad and Jensen (1980) proposed a semi-empirical correction factor:
//!
//! ```text
//!     K_HJ(f) = 1 + (2/π) · arctan[1.4 · (Δ / δₛ(f))²]
//! ```
//!
//! where:
//! - Δ  = RMS surface roughness [m]
//! - δₛ = skin depth = sqrt(ρ / (π f μ₀)) [m]
//! - K_HJ ∈ [1, 2]; K_HJ → 1 for δₛ ≫ Δ (low freq), K_HJ → 2 for δₛ ≪ Δ (high freq)
//!
//! The corrected attenuation constant and surface impedance are:
//! ```text
//!     α_c_corrected = α_c · K_HJ
//!     Z_s_corrected = Z_s · K_HJ
//! ```
//!
//! # Groisse Model
//!
//! Groisse et al. (1994) provide a slightly more accurate fit to measurement data:
//!
//! ```text
//!     K_G(f) = 1 / [1 − exp(−(δₛ / Δ)^1.6)]
//! ```
//!
//! # Cannonball–Huray Model
//!
//! The Cannonball model (Huray 2010) treats the surface roughness as a pile of
//! hemispherical bosses of radius *a* with tile area *A_flat*:
//!
//! ```text
//!     K_CB(f) = 1 + (A_sphere / A_flat) · [1 − exp(−2a/δₛ)]
//!     A_sphere = 2π a²  (half-sphere surface area)
//! ```
//!
//! This model better captures the saturation behaviour at high frequencies.
//!
//! # Usage
//!
//! ```rust
//! use rem_planar::roughness::{skin_depth, hammerstad_jensen_factor, apply_roughness_to_loss};
//!
//! let delta_rms = 1.5e-6;  // 1.5 µm RMS roughness
//! let rho_cu    = 1.72e-8; // copper resistivity [Ω·m]
//! let freq      = 10.0e9;  // 10 GHz
//!
//! let delta_s   = skin_depth(freq, rho_cu);
//! let k_hj      = hammerstad_jensen_factor(delta_rms, delta_s);
//! let alpha_raw = 0.05;    // smooth-conductor attenuation [Np/m]
//! let alpha_hj  = apply_roughness_to_loss(alpha_raw, k_hj);
//! ```
//!
//! # References
//! Hammerstad, E. & Jensen, O. (1980) "Accurate Models for Microstrip Computer-Aided
//! Design." *IEEE MTT-S International Microwave Symposium Digest*, pp. 407–409.
//!
//! Groisse, P. et al. (1994) "Parameters for the global RLCG transmission line model."
//! *Proc. 3rd IEEE Int. Conf. on Electromagnetic Compatibility*.
//!
//! Huray, P. G. et al. (2010) "Fundamentals of a 3-D snowball model for surface
//! roughness power losses." *IEEE Proceedings*, *SI8P*.

use std::f64::consts::PI;

/// Permeability of free space [H/m].
const MU0: f64 = 4.0e-7 * PI;

/// Compute the electromagnetic skin depth at frequency `freq_hz` for a conductor
/// with bulk resistivity `rho_ohm_m` [Ω·m].
///
/// ```text
///     δₛ = sqrt(ρ / (π f μ₀))
/// ```
///
/// Returns skin depth in metres.
pub fn skin_depth(freq_hz: f64, rho_ohm_m: f64) -> f64 {
    if freq_hz <= 0.0 || rho_ohm_m <= 0.0 {
        return f64::INFINITY;
    }
    (rho_ohm_m / (PI * freq_hz * MU0)).sqrt()
}

/// Compute the Hammerstad–Jensen roughness correction factor K_HJ.
///
/// # Arguments
/// * `rms_roughness_m` — RMS surface roughness Δ [m] (typical PCB Cu: 0.5–3 µm)
/// * `skin_depth_m`    — skin depth δₛ at the frequency of interest [m]
///
/// # Returns
/// K_HJ ∈ [1, 2]; multiply conductor loss or surface impedance by this factor.
pub fn hammerstad_jensen_factor(rms_roughness_m: f64, skin_depth_m: f64) -> f64 {
    if skin_depth_m <= 0.0 || rms_roughness_m <= 0.0 {
        return 1.0;
    }
    let ratio = rms_roughness_m / skin_depth_m;
    1.0 + (2.0 / PI) * (1.4 * ratio * ratio).atan()
}

/// Compute the Groisse roughness correction factor K_G.
///
/// # Arguments
/// * `rms_roughness_m` — RMS surface roughness Δ [m]
/// * `skin_depth_m`    — skin depth δₛ [m]
///
/// # Returns
/// K_G ≥ 1.0; multiply conductor loss or surface impedance by this factor.
pub fn groisse_factor(rms_roughness_m: f64, skin_depth_m: f64) -> f64 {
    if skin_depth_m <= 0.0 || rms_roughness_m <= 0.0 {
        return 1.0;
    }
    let ratio = skin_depth_m / rms_roughness_m;
    let exponent = -(ratio.powf(1.6));
    1.0 / (1.0 - exponent.exp())
}

/// Compute the Cannonball–Huray roughness correction factor K_CB.
///
/// # Arguments
/// * `boss_radius_m`   — radius of hemispherical bosses [m]
/// * `tile_area_m2`    — flat tile area per boss (1/boss_density) [m²]
/// * `skin_depth_m`    — skin depth δₛ [m]
///
/// # Returns
/// K_CB ≥ 1.0; multiply conductor loss by this factor.
pub fn cannonball_huray_factor(boss_radius_m: f64, tile_area_m2: f64, skin_depth_m: f64) -> f64 {
    if skin_depth_m <= 0.0 || boss_radius_m <= 0.0 || tile_area_m2 <= 0.0 {
        return 1.0;
    }
    let a_sphere = 2.0 * PI * boss_radius_m * boss_radius_m; // half-sphere area
    let corr     = 1.0 - (-2.0 * boss_radius_m / skin_depth_m).exp();
    1.0 + (a_sphere / tile_area_m2) * corr
}

/// Apply a roughness correction factor to a conductor loss (attenuation constant
/// or surface resistance).
///
/// # Arguments
/// * `smooth_loss` — loss for a perfectly smooth conductor (any consistent unit)
/// * `k_factor`    — correction factor (K_HJ, K_G, or K_CB)
///
/// # Returns
/// Corrected loss = `smooth_loss * k_factor`.
#[inline]
pub fn apply_roughness_to_loss(smooth_loss: f64, k_factor: f64) -> f64 {
    smooth_loss * k_factor
}

/// Frequency sweep: compute Hammerstad–Jensen factor at a set of frequencies.
///
/// Returns a `Vec<(freq_hz, skin_depth_m, k_hj)>` for each frequency point.
pub fn hj_sweep(
    freqs_hz: &[f64],
    rms_roughness_m: f64,
    rho_ohm_m: f64,
) -> Vec<(f64, f64, f64)> {
    freqs_hz.iter().map(|&f| {
        let ds = skin_depth(f, rho_ohm_m);
        let k  = hammerstad_jensen_factor(rms_roughness_m, ds);
        (f, ds, k)
    }).collect()
}

/// Write roughness sweep results to CSV.
pub fn write_roughness_csv(
    data: &[(f64, f64, f64)],
    output_dir: &std::path::Path,
) -> Result<std::path::PathBuf, std::io::Error> {
    use std::io::Write;
    let dir = output_dir.join("postpro");
    std::fs::create_dir_all(&dir)?;
    let path = dir.join("conductor_roughness.csv");
    let mut f = std::fs::File::create(&path)?;
    writeln!(f, "FreqHz,SkinDepth_m,K_HJ")?;
    for &(freq, ds, k) in data {
        writeln!(f, "{:.9e},{:.9e},{:.9e}", freq, ds, k)?;
    }
    Ok(path)
}

// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// At f → 0 (large δₛ), K_HJ → 1 (no roughness effect at low frequency).
    #[test]
    fn hj_low_freq_limit_approaches_one() {
        let delta = 1.5e-6; // 1.5 µm roughness
        let large_delta_s = 1.0; // 1 m skin depth ≫ roughness
        let k = hammerstad_jensen_factor(delta, large_delta_s);
        assert!((k - 1.0).abs() < 0.01, "K_HJ at low freq = {k:.6}, expected ≈ 1");
    }

    /// At f → ∞ (δₛ → 0), K_HJ → 2.
    #[test]
    fn hj_high_freq_limit_approaches_two() {
        let delta = 1.5e-6;
        let tiny_delta_s = 1e-12; // essentially zero
        let k = hammerstad_jensen_factor(delta, tiny_delta_s);
        assert!((k - 2.0).abs() < 0.01, "K_HJ at high freq = {k:.6}, expected ≈ 2");
    }

    /// K_HJ is always in [1, 2].
    #[test]
    fn hj_bounds_always_met() {
        let rho_cu = 1.72e-8;
        let delta   = 2.0e-6;
        for &f in &[1e6_f64, 1e8, 1e9, 5e9, 10e9, 100e9] {
            let ds = skin_depth(f, rho_cu);
            let k  = hammerstad_jensen_factor(delta, ds);
            assert!(k >= 1.0 && k <= 2.0,
                "K_HJ={k:.6} out of [1,2] at f={f:.2e} Hz");
        }
    }

    /// Skin depth for copper at 1 GHz should be ≈ 2.1 µm.
    #[test]
    fn copper_skin_depth_1ghz() {
        let rho_cu = 1.72e-8;
        let ds = skin_depth(1.0e9, rho_cu);
        // Expected: sqrt(1.72e-8 / (π · 1e9 · 4π×10⁻⁷)) ≈ 2.09 µm
        assert!((ds - 2.09e-6).abs() < 0.1e-6,
            "δₛ(Cu,1GHz) = {:.3e} m, expected ≈ 2.09e-6 m", ds);
    }

    /// Groisse factor: K_G ≥ 1 for all inputs; for large δₛ (low freq) → ≈ 1.
    #[test]
    fn groisse_low_freq_near_one() {
        let delta   = 1.5e-6;
        let large_ds = 1.0;
        let k = groisse_factor(delta, large_ds);
        // For δₛ/Δ = 1/1.5e-6 >> 1: K_G = 1/(1-exp(-very_large)) ≈ 1.
        assert!((k - 1.0).abs() < 0.01, "Groisse K at low freq = {k:.6}");
    }

    /// Cannonball factor should be ≥ 1 always.
    #[test]
    fn cannonball_always_gte_one() {
        let a  = 0.5e-6;     // 0.5 µm boss radius
        let at = 4.0e-12;    // tile area 2×2 µm²
        let rho_cu = 1.72e-8;
        for &f in &[1e8_f64, 1e9, 10e9, 100e9] {
            let ds = skin_depth(f, rho_cu);
            let k  = cannonball_huray_factor(a, at, ds);
            assert!(k >= 1.0, "Cannonball K={k:.4} < 1 at {f:.1e} Hz");
        }
    }

    /// Zero/negative inputs return K = 1 (graceful).
    #[test]
    fn degenerate_inputs_return_one() {
        assert_eq!(hammerstad_jensen_factor(0.0, 1e-6), 1.0);
        assert_eq!(hammerstad_jensen_factor(1e-6, 0.0), 1.0);
        assert_eq!(groisse_factor(0.0, 1e-6), 1.0);
        assert_eq!(cannonball_huray_factor(0.0, 1e-12, 1e-6), 1.0);
        assert_eq!(skin_depth(0.0, 1.72e-8), f64::INFINITY);
    }

    /// apply_roughness_to_loss is just multiplication.
    #[test]
    fn apply_loss_basic() {
        let a = apply_roughness_to_loss(1.5, 1.8);
        assert!((a - 2.7).abs() < 1e-10);
    }
}
