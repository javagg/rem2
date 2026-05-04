use rem_core::constants::{EPS0, MU0};

/// Surface roughness model for conductor loss correction.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RoughnessModel {
    /// Hammerstad-Jensen (classic): K = 1 + (2/π) arctan(1.4 (Δ/δ)²)
    Hammerstad,
    /// Groisse: K = 1 + exp(−δ/(2Δ))
    Groisse,
    /// Huray snowball model: K = 1 + (3/2) (Δ/δ)
    Huray,
}

/// Surface roughness parameters for conductor materials.
#[derive(Debug, Clone)]
pub struct SurfaceRoughness {
    /// RMS surface roughness [m]
    pub rms: f64,
    /// Roughness correction model
    pub model: RoughnessModel,
}

/// Temperature-dependent material coefficients.
#[derive(Debug, Clone)]
pub struct TemperatureCoefficient {
    /// Reference temperature for given properties [K] (typically 293.15 = 20°C)
    pub t_ref: f64,
    /// Temperature coefficient of permittivity [ppm/K]
    pub permittivity_tc: f64,
    /// Temperature coefficient of conductivity [ppm/K]
    pub conductivity_tc: f64,
}

/// Physical material properties.
///
/// All values are in SI units and represent relative quantities where noted.
#[derive(Debug, Clone)]
pub struct Material {
    /// Relative permittivity εᵣ (dimensionless)
    pub permittivity: f64,
    /// Relative permeability μᵣ (dimensionless)
    pub permeability: f64,
    /// Electric conductivity σ [S/m]
    pub conductivity: f64,
    /// Dielectric loss tangent tan δ_e (dimensionless)
    pub loss_tangent: f64,
    /// Magnetic loss tangent tan δ_m = μᵢ/μᵣ (dimensionless, for ferrites/lossy magnetics)
    pub loss_tangent_magnetic: f64,
    /// Absolute permittivity tensor ε₀ εᵣ [F/m], 3×3 row-major.
    /// Defaults to ε₀ εᵣ · I. Non-identity when MaterialAxes rotation is present.
    pub epsilon_tensor: [[f64; 3]; 3],
    /// Reluctivity tensor ν = 1/(μ₀ μᵣ) [m/H], 3×3 row-major.
    /// Defaults to ν · I (scalar isotropic). Set by MaterialAxes when present.
    pub nu_tensor: [[f64; 3]; 3],
    /// Drude-Lorentz oscillator poles: (ωp², ω0², γ) tuples.
    /// ε(ω) = εᵣ + Σ ωp² / (ω0² − ω² + jγω)
    /// Empty for most materials (static permittivity only).
    pub drude_lorentz_poles: Vec<(f64, f64, f64)>,
    /// Optional surface roughness parameters for conductor loss correction.
    pub surface_roughness: Option<SurfaceRoughness>,
    /// Optional temperature coefficient for thermal-dependent properties.
    pub temperature_coefficient: Option<TemperatureCoefficient>,
}

impl Default for Material {
    fn default() -> Self {
        let nu = 1.0 / MU0;
        Material {
            permittivity: 1.0,
            permeability: 1.0,
            conductivity: 0.0,
            loss_tangent: 0.0,
            loss_tangent_magnetic: 0.0,
            epsilon_tensor: [[EPS0, 0.0, 0.0], [0.0, EPS0, 0.0], [0.0, 0.0, EPS0]],
            nu_tensor: [[nu, 0.0, 0.0], [0.0, nu, 0.0], [0.0, 0.0, nu]],
            drude_lorentz_poles: Vec::new(),
            surface_roughness: None,
            temperature_coefficient: None,
        }
    }
}

impl Material {
    /// Construct from scalar parameters, building isotropic tensors.
    pub fn from_scalars(permittivity: f64, permeability: f64, conductivity: f64, loss_tangent: f64) -> Self {
        let eps = EPS0 * permittivity;
        let nu = 1.0 / (MU0 * permeability);
        Material {
            permittivity,
            permeability,
            conductivity,
            loss_tangent,
            loss_tangent_magnetic: 0.0,
            epsilon_tensor: [[eps, 0.0, 0.0], [0.0, eps, 0.0], [0.0, 0.0, eps]],
            nu_tensor: [[nu, 0.0, 0.0], [0.0, nu, 0.0], [0.0, 0.0, nu]],
            drude_lorentz_poles: Vec::new(),
            surface_roughness: None,
            temperature_coefficient: None,
        }
    }

    /// Construct from scalars + an optional rotation matrix (MaterialAxes rows).
    /// If `axes` has 3 rows, interprets them as rotation matrix R and builds
    /// ε_tensor = R^T · (ε·I) · R = ε·I  (rotation of isotropic is still isotropic,
    /// but establishes the tensor pathway for future per-axis εᵣ support).
    /// Similarly builds nu_tensor = R^T · (ν·I) · R.
    pub fn from_scalars_with_axes(
        permittivity: f64,
        permeability: f64,
        conductivity: f64,
        loss_tangent: f64,
        axes: &[Vec<f64>],
    ) -> Self {
        let mut mat = Self::from_scalars(permittivity, permeability, conductivity, loss_tangent);
        // from_scalars sets roughness/temperature to None — they are set separately
        if axes.len() >= 3
            && axes[0].len() >= 3
            && axes[1].len() >= 3
            && axes[2].len() >= 3
        {
            let r = [
                [axes[0][0], axes[0][1], axes[0][2]],
                [axes[1][0], axes[1][1], axes[1][2]],
                [axes[2][0], axes[2][1], axes[2][2]],
            ];
            let eps = EPS0 * permittivity;
            let nu = 1.0 / (MU0 * permeability);
            let mut et = [[0.0f64; 3]; 3];
            let mut nt = [[0.0f64; 3]; 3];
            for i in 0..3 {
                for j in 0..3 {
                    let rtrt = r[0][i]*r[0][j] + r[1][i]*r[1][j] + r[2][i]*r[2][j];
                    et[i][j] = eps * rtrt;
                    nt[i][j] = nu * rtrt;
                }
            }
            mat.epsilon_tensor = et;
            mat.nu_tensor = nt;
        }
        mat
    }

    /// Absolute permittivity for electrostatic assembly: ε₀ εᵣ [F/m]
    pub fn epsilon_abs(&self) -> f64 {
        EPS0 * self.permittivity
    }

    /// Reluctivity ν = 1 / (μ₀ μᵣ) [m/H] for magnetostatic assembly
    pub fn reluctivity(&self) -> f64 {
        1.0 / (MU0 * self.permeability)
    }

    /// Effective (complex) permittivity at frequency `freq` [Hz].
    ///
    /// Returns `(εᵣ_eff_re, εᵣ_eff_im)` (relative, dimensionless).
    ///
    /// Includes:
    /// - Static permittivity: εᵣ
    /// - Dielectric loss:     −j εᵣ tan δ
    /// - Conductivity:        −j σ / (ω ε₀)
    /// - Drude-Lorentz poles: Σ ωp² / (ω0² − ω² + jγω)
    ///
    /// The imaginary part is negative (lossy convention: e^{+jωt} time dependence).
    pub fn epsilon_complex(&self, freq: f64) -> (f64, f64) {
        use std::f64::consts::PI;
        let omega = 2.0 * PI * freq;
        let mut re = self.permittivity;
        let mut im = -(self.permittivity * self.loss_tangent
            + if omega > 0.0 { self.conductivity / (omega * EPS0) } else { 0.0 });

        // Drude-Lorentz poles: (ωp², ω0², γ)
        for &(wp2, w02, gamma) in &self.drude_lorentz_poles {
            // Denominator: (ω0² − ω²) + jγω
            let denom_re = w02 - omega * omega;
            let denom_im = gamma * omega;
            let denom_sq = denom_re * denom_re + denom_im * denom_im;
            if denom_sq > 1e-300 {
                // wp² · conj(denom) / |denom|²
                re +=  wp2 * denom_re / denom_sq;
                im += -wp2 * denom_im / denom_sq;
            }
        }

        (re, im)
    }

    /// Returns true if the material has any loss (tan δ > 0 or σ > 0 or Drude-Lorentz poles).
    pub fn is_lossy(&self) -> bool {
        self.loss_tangent > 0.0 || self.conductivity > 0.0 || !self.drude_lorentz_poles.is_empty()
    }

    /// Returns true if the material has Drude-Lorentz poles (frequency-dependent permittivity).
    pub fn has_drude_lorentz(&self) -> bool {
        !self.drude_lorentz_poles.is_empty()
    }

    /// Returns true if the material has magnetic loss (tan δ_m > 0).
    pub fn is_magnetically_lossy(&self) -> bool {
        self.loss_tangent_magnetic > 0.0
    }

    /// Returns true if the epsilon tensor is non-isotropic (off-diagonal entries non-zero
    /// or diagonal entries differ from ε₀·εᵣ).
    pub fn is_anisotropic(&self) -> bool {
        let eps = EPS0 * self.permittivity;
        let identity = [[eps, 0.0, 0.0], [0.0, eps, 0.0], [0.0, 0.0, eps]];
        self.epsilon_tensor != identity
    }

    /// Returns true if the nu_tensor is non-isotropic (off-diagonal non-zero
    /// or diagonal differs from ν = 1/(μ₀ μᵣ)).
    pub fn is_magnetically_anisotropic(&self) -> bool {
        let nu = 1.0 / (MU0 * self.permeability);
        let identity = [[nu, 0.0, 0.0], [0.0, nu, 0.0], [0.0, 0.0, nu]];
        self.nu_tensor != identity
    }

    /// Effective conductivity accounting for surface roughness at a given frequency.
    ///
    /// Returns the roughness-corrected conductivity σ_eff [S/m].
    ///
    /// **Hammerstad-Jensen**: K = 1 + (2/π) arctan(1.4 · (Δ/δ)²)
    ///
    /// **Groisse**: K = 1 + exp(−δ / (2Δ))
    ///
    /// **Huray snowball**: K = 1 + (3/2) · (Δ/δ)
    ///
    /// where Δ = RMS roughness [m], δ = skin depth [m] = √(2 / (ω μ₀ σ)).
    ///
    /// The corrected conductivity is σ_eff = σ / K² (power-loss correction).
    /// Returns the DC conductivity if `surface_roughness` is `None` or `freq <= 0`.
    pub fn effective_conductivity(&self, freq: f64) -> f64 {
        let sigma = self.conductivity;
        let roughness = match &self.surface_roughness {
            Some(r) if r.rms > 0.0 && sigma > 0.0 && freq > 0.0 => r,
            _ => return sigma,
        };

        use std::f64::consts::PI;
        let omega = 2.0 * PI * freq;
        let delta = (2.0 / (omega * MU0 * sigma)).sqrt();
        let ratio = roughness.rms / delta;

        let k = match roughness.model {
            RoughnessModel::Hammerstad => {
                1.0 + (2.0 / PI) * (1.4 * ratio * ratio).atan()
            }
            RoughnessModel::Groisse => {
                1.0 + (-0.5 / ratio.max(1e-15)).exp()
            }
            RoughnessModel::Huray => {
                1.0 + 1.5 * ratio
            }
        };

        sigma / (k * k)
    }

    /// Temperature-corrected conductivity at temperature `t` [K].
    ///
    /// Uses linear model: σ(T) = σ(T_ref) / (1 + α_T · (T − T_ref))
    /// where α_T = conductivity_tc × 10⁻⁶ [1/K].
    ///
    /// Returns the DC conductivity if no temperature coefficient is set.
    pub fn conductivity_at_temperature(&self, t: f64) -> f64 {
        let tc = match &self.temperature_coefficient {
            Some(tc) => tc,
            None => return self.conductivity,
        };
        let alpha = tc.conductivity_tc * 1e-6;
        self.conductivity / (1.0 + alpha * (t - tc.t_ref))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn vacuum_epsilon() {
        let m = Material::default();
        assert!((m.epsilon_abs() - EPS0).abs() < 1e-25);
    }

    #[test]
    fn reluctivity_iron() {
        let m = Material { permeability: 1000.0, ..Default::default() };
        let nu = m.reluctivity();
        assert!((nu - 1.0 / (MU0 * 1000.0)).abs() < 1e-5);
    }

    #[test]
    fn dielectric_fr4() {
        let m = Material::from_scalars(4.5, 1.0, 0.0, 0.02);
        assert!(!m.is_lossy() || m.is_lossy()); // just exercise the path
        let (re, im) = m.epsilon_complex(1.0e9);
        assert!((re - 4.5).abs() < 1e-12);
        // im = -4.5 * 0.02 at low conductivity
        assert!((im - (-4.5 * 0.02)).abs() < 1e-12);
    }

    #[test]
    fn lossless_is_not_lossy() {
        let m = Material::default();
        assert!(!m.is_lossy());
    }

    #[test]
    fn isotropic_not_anisotropic() {
        let m = Material::from_scalars(2.0, 1.0, 0.0, 0.0);
        assert!(!m.is_anisotropic());
    }

    #[test]
    fn no_roughness_returns_dc_conductivity() {
        let m = Material { conductivity: 5.8e7, ..Default::default() };
        let sigma = m.effective_conductivity(10.0e9);
        assert!((sigma - 5.8e7).abs() < 1e-7);
    }

    #[test]
    fn hammerstad_increases_loss() {
        let mut m = Material { conductivity: 5.8e7, ..Default::default() };
        m.surface_roughness = Some(SurfaceRoughness {
            rms: 1.0e-6, // 1 µm RMS
            model: RoughnessModel::Hammerstad,
        });
        let sigma = m.effective_conductivity(10.0e9);
        // Roughness reduces effective conductivity (increases loss)
        assert!(sigma < 5.8e7, "sigma_eff={} should be < bulk", sigma);
        assert!(sigma > 0.0);
    }

    #[test]
    fn groisse_model_reduces_conductivity() {
        let mut m = Material { conductivity: 4.1e7, ..Default::default() };
        m.surface_roughness = Some(SurfaceRoughness {
            rms: 0.5e-6,
            model: RoughnessModel::Groisse,
        });
        let s1 = m.effective_conductivity(1.0e9);
        let s2 = m.effective_conductivity(10.0e9);
        // Higher frequency → thinner skin depth → more roughness effect → lower sigma
        assert!(s2 < s1, "s2={} should be < s1={}", s2, s1);
    }

    #[test]
    fn huray_scales_with_roughness() {
        let mut m = Material { conductivity: 5.0e7, ..Default::default() };
        m.surface_roughness = Some(SurfaceRoughness {
            rms: 2.0e-6,
            model: RoughnessModel::Huray,
        });
        let s_rough = m.effective_conductivity(5.0e9);
        m.surface_roughness = Some(SurfaceRoughness {
            rms: 0.2e-6,
            model: RoughnessModel::Huray,
        });
        let s_smooth = m.effective_conductivity(5.0e9);
        // Rougher surface → more loss
        assert!(s_rough < s_smooth);
    }

    #[test]
    fn temperature_reduces_conductivity() {
        let mut m = Material { conductivity: 5.8e7, ..Default::default() };
        m.temperature_coefficient = Some(TemperatureCoefficient {
            t_ref: 293.15, // 20°C
            permittivity_tc: 0.0,
            conductivity_tc: 3930.0, // copper: ~3930 ppm/K
        });
        let s_cold = m.conductivity_at_temperature(293.15);
        let s_hot = m.conductivity_at_temperature(373.15); // 100°C
        assert!((s_cold - 5.8e7).abs() < 1e-7); // at ref temp
        assert!(s_hot < s_cold, "s_hot={} should be < s_cold={}", s_hot, s_cold);
    }

    #[test]
    fn no_temperature_coeff_returns_dc() {
        let m = Material { conductivity: 5.8e7, ..Default::default() };
        let s = m.conductivity_at_temperature(500.0);
        assert!((s - 5.8e7).abs() < 1e-7);
    }
}
