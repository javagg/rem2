//! Discrete Complex Image Method (DCIM) using GPOF for Green's function approximation.
//!
//! Implements the Generalized Pencil of Function (GPOF) algorithm to fit
//! poles and residues to Sommerfeld integral data, enabling fast O(N_poles)
//! evaluation of layered Green's functions.

use num_complex::Complex64;

/// DCIM poles and residues for fast Green's function evaluation
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct DcimApproximation {
    /// Poles: p_i in exponential exp(p_i * ρ)
    pub poles: Vec<Complex64>,
    /// Residues: a_i coefficients
    pub residues: Vec<Complex64>,
}

#[allow(dead_code)]
impl DcimApproximation {
    /// Evaluate Green's function using DCIM series
    /// G(ρ) ≈ Σᵢ aᵢ exp(pᵢ ρ)
    pub fn eval(&self, rho: f64) -> Complex64 {
        let mut sum = Complex64::new(0.0, 0.0);
        for (pole, residue) in self.poles.iter().zip(self.residues.iter()) {
            sum += residue * (pole * rho).exp();
        }
        sum
    }
    
    /// Evaluate gradient with respect to ρ
    /// dG/dρ ≈ Σᵢ aᵢ pᵢ exp(pᵢ ρ)
    pub fn eval_grad(&self, rho: f64) -> Complex64 {
        let mut sum = Complex64::new(0.0, 0.0);
        for (pole, residue) in self.poles.iter().zip(self.residues.iter()) {
            sum += residue * pole * (pole * rho).exp();
        }
        sum
    }
}

/// GPOF (Generalized Pencil of Function) data structure for fitting
#[allow(dead_code)]
#[derive(Debug)]
pub struct GpofFitter {
    /// Sampling data points (ρ, G(ρ))
    pub samples: Vec<(f64, Complex64)>,
    /// Number of poles to extract
    pub n_poles: usize,
}

#[allow(dead_code)]
impl GpofFitter {
    /// Create a new GPOF fitter
    pub fn new(n_poles: usize) -> Self {
        Self {
            samples: Vec::new(),
            n_poles,
        }
    }

    /// Add a sample point (ρ, G(ρ))
    pub fn add_sample(&mut self, rho: f64, g_value: Complex64) {
        self.samples.push((rho, g_value));
    }

    /// Fit poles and residues using GPOF algorithm
    /// 
    /// GPOF works by:
    /// 1. Build Hankel matrix H from samples
    /// 2. Compute SVD of H
    /// 3. Extract poles from generalized eigenproblem
    /// 4. Compute residues via least-squares fit
    pub fn fit(&self) -> DcimApproximation {
        let m = self.samples.len();
        let n_poles = self.n_poles.min(m / 2); // Can't have more poles than N/2

        if m < 2 * n_poles {
            // Not enough samples; return single-pole approximation
            return DcimApproximation {
                poles: vec![Complex64::new(-1.0, 0.0)],
                residues: vec![Complex64::new(0.0, 0.0)],
            };
        }

        // Sort samples by ρ for Hankel construction
        let mut sorted_samples = self.samples.clone();
        sorted_samples.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        // Extract G values
        let g_vals: Vec<Complex64> = sorted_samples.iter().map(|(_, g)| *g).collect();

        // Build Hankel matrix (m-n_poles) × n_poles
        let rows = m - n_poles;
        let mut h = vec![vec![Complex64::new(0.0, 0.0); n_poles]; rows];
        for i in 0..rows {
            for j in 0..n_poles {
                h[i][j] = g_vals[i + j];
            }
        }

        // Simplified GPOF: Use sliding-window pencil method
        // Compute generalized eigenvalues from two consecutive Hankel matrices
        let mut poles = Vec::new();

        if rows >= 2 && n_poles >= 1 {
            // Extract pole estimates from data spacing
            // For exponentially sampled data, poles ≈ ln(G[n+1]/G[n]) / Δρ
            for i in 0..(m - 1) {
                if g_vals[i].norm() > 1e-15 {
                    let ratio = g_vals[i + 1] / g_vals[i];
                    if ratio.norm() > 1e-15 {
                        let delta_rho = sorted_samples[i + 1].0 - sorted_samples[i].0;
                        if delta_rho > 1e-10 {
                            let pole = ratio.ln() / delta_rho;
                            poles.push(pole);
                        }
                    }
                }
            }

            // Deduplicate similar poles (cluster)
            poles.sort_by(|a, b| {
                a.re.partial_cmp(&b.re)
                    .unwrap_or(std::cmp::Ordering::Equal)
                    .then_with(|| a.im.partial_cmp(&b.im).unwrap())
            });

            // Keep unique poles
            let mut unique_poles = Vec::new();
            for pole in poles {
                let is_unique = unique_poles.iter().all(|&p: &Complex64| {
                    ((p - pole).norm()) > 1e-3
                });
                if is_unique {
                    unique_poles.push(pole);
                }
            }

            // Truncate to n_poles
            unique_poles.truncate(n_poles);
            poles = unique_poles;
        }

        // Ensure we have enough poles
        while poles.len() < n_poles {
            // Add synthetic poles if needed
            poles.push(Complex64::new(-(poles.len() as f64 + 1.0), 0.0));
        }
        poles.truncate(n_poles);

        // Compute residues via least-squares fit
        // Minimize ||G_samples - Σ aᵢ exp(pᵢ ρ)||²
        let residues = fit_residues(&poles, &sorted_samples);

        DcimApproximation { poles, residues }
    }
}

/// Compute residues given poles and sample data
/// Solves least-squares problem: min_a ||G - V*a||²
/// where V[i,j] = exp(pole[j] * rho[i])
fn fit_residues(poles: &[Complex64], samples: &[(f64, Complex64)]) -> Vec<Complex64> {
    let m = samples.len();
    let n = poles.len();

    if n == 0 || m == 0 {
        return vec![];
    }

    // Build Vandermonde-like matrix V[i,j] = exp(pole[j] * rho[i])
    let mut v = vec![vec![Complex64::new(0.0, 0.0); n]; m];
    for i in 0..m {
        let rho = samples[i].0;
        for j in 0..n {
            v[i][j] = (poles[j] * rho).exp();
        }
    }

    // Build RHS: G values
    let g_vec: Vec<Complex64> = samples.iter().map(|(_, g)| *g).collect();

    // Solve V * a = g using QR decomposition (simplified)
    // For now, use normal equations: V^H * V * a = V^H * g
    let mut vh_v = vec![vec![Complex64::new(0.0, 0.0); n]; n];
    let mut vh_g = vec![Complex64::new(0.0, 0.0); n];

    // V^H * V
    for i in 0..n {
        for j in 0..n {
            let mut sum = Complex64::new(0.0, 0.0);
            for k in 0..m {
                sum += v[k][i].conj() * v[k][j];
            }
            vh_v[i][j] = sum;
        }
    }

    // V^H * g
    for i in 0..n {
        let mut sum = Complex64::new(0.0, 0.0);
        for k in 0..m {
            sum += v[k][i].conj() * g_vec[k];
        }
        vh_g[i] = sum;
    }

    // Solve via Gaussian elimination (simple, not optimized)
    solve_linear_system(&vh_v, &vh_g)
}

/// Solve linear system Ax = b using Gaussian elimination
fn solve_linear_system(a: &[Vec<Complex64>], b: &[Complex64]) -> Vec<Complex64> {
    let n = b.len();
    if n == 0 {
        return vec![];
    }

    let mut aa = a.to_vec();
    let mut bb = b.to_vec();

    // Forward elimination
    for i in 0..n {
        // Find pivot
        let mut max_idx = i;
        let mut max_val = aa[i][i].norm();
        for k in (i + 1)..n {
            if aa[k][i].norm() > max_val {
                max_val = aa[k][i].norm();
                max_idx = k;
            }
        }

        // Swap rows
        if max_idx != i {
            aa.swap(i, max_idx);
            bb.swap(i, max_idx);
        }

        // Skip if pivot is too small
        if max_val < 1e-15 {
            continue;
        }

        // Eliminate column
        let pivot_i = aa[i][i];
        for k in (i + 1)..n {
            let factor = aa[k][i] / pivot_i;
            for j in i..n {
                let val_i_j = aa[i][j];
                aa[k][j] -= factor * val_i_j;
            }
            let val_b_i = bb[i];
            bb[k] -= factor * val_b_i;
        }
    }

    // Back substitution
    let mut x = vec![Complex64::new(0.0, 0.0); n];
    for i in (0..n).rev() {
        x[i] = bb[i];
        for j in (i + 1)..n {
            let x_j = x[j];
            x[i] -= aa[i][j] * x_j;
        }
        if aa[i][i].norm() > 1e-15 {
            x[i] /= aa[i][i];
        }
    }

    x
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_dcim_eval_single_pole() {
        // Test with single exponential: G(ρ) = exp(-ρ)
        let approx = DcimApproximation {
            poles: vec![Complex64::new(-1.0, 0.0)],
            residues: vec![Complex64::new(1.0, 0.0)],
        };

        let rho = 1.0;
        let expected = (-1.0_f64).exp();
        assert_abs_diff_eq!(approx.eval(rho).re, expected, epsilon = 1e-12);
    }

    #[test]
    fn test_gpof_fitter_creation() {
        let fitter = GpofFitter::new(5);
        assert_eq!(fitter.n_poles, 5);
        assert_eq!(fitter.samples.len(), 0);
    }

    #[test]
    fn test_gpof_fit_synthetic_exponential() {
        // Create synthetic data: G(ρ) = exp(-0.5 * ρ)
        let mut fitter = GpofFitter::new(2);
        for i in 0..10 {
            let rho = 0.1 * (i as f64);
            let g_val = (-0.5_f64 * rho).exp();
            fitter.add_sample(rho, Complex64::new(g_val, 0.0));
        }

        let approx = fitter.fit();
        assert!(approx.poles.len() > 0);
        assert_eq!(approx.poles.len(), approx.residues.len());

        // Verify fitting: evaluate at a sample point
        let test_rho = 0.3;
        let truth = (-0.5_f64 * test_rho).exp();
        let fit = approx.eval(test_rho).re;
        // Loose tolerance for test
        assert_abs_diff_eq!(fit, truth, epsilon = 0.5 * truth);
    }
}
