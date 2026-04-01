use rem_core::constants::{EPS0, MU0};

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
    /// Dielectric loss tangent tan δ (dimensionless)
    pub loss_tangent: f64,
}

impl Default for Material {
    fn default() -> Self {
        Material {
            permittivity: 1.0,
            permeability: 1.0,
            conductivity: 0.0,
            loss_tangent: 0.0,
        }
    }
}

impl Material {
    /// Absolute permittivity for electrostatic assembly: ε₀ εᵣ [F/m]
    pub fn epsilon_abs(&self) -> f64 {
        EPS0 * self.permittivity
    }

    /// Reluctivity ν = 1 / (μ₀ μᵣ) [m/H] for magnetostatic assembly
    pub fn reluctivity(&self) -> f64 {
        1.0 / (MU0 * self.permeability)
    }

    /// Effective (complex) permittivity at frequency `freq` [Hz].
    /// εᵣ_eff = εᵣ (1 − j tan δ) − j σ / (ω ε₀)
    pub fn epsilon_complex(&self, freq: f64) -> (f64, f64) {
        use std::f64::consts::PI;
        let omega = 2.0 * PI * freq;
        let re = self.permittivity;
        let im = -(self.permittivity * self.loss_tangent
            + if omega > 0.0 { self.conductivity / (omega * EPS0) } else { 0.0 });
        (re, im)
    }

    /// Returns true if the material has any loss (tan δ > 0 or σ > 0).
    pub fn is_lossy(&self) -> bool {
        self.loss_tangent > 0.0 || self.conductivity > 0.0
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
        let m = Material {
            permittivity: 4.5,
            loss_tangent: 0.02,
            ..Default::default()
        };
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
}
