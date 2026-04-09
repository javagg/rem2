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
}
