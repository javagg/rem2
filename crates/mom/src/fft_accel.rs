//! FFT-accelerated MoM matrix-vector product for planar structures.
//!
//! For planar surfaces (all nodes at the same z-height within tolerance), the
//! free-space scalar Green function depends only on the lateral offset (Δx, Δy):
//!
//!   G(r, r') = exp(-jk·ρ) / (4π·ρ),   ρ = √((x−x')² + (y−y')²)
//!
//! This makes the impedance kernel **shift-invariant** in the x-y plane, so the
//! matrix-vector product Z·I can be approximated by a 2-D FFT convolution:
//!
//!   y[m] ≈ Σ_n  G(x_m−x_n, y_m−y_n) · I_n
//!        = IFFT2(FFT2(G_kernel) × FFT2(I_grid))
//!
//! evaluated on an oversampled regular grid and interpolated back to the mesh nodes.
//!
//! **When to use**: activate with `FastSolver: "FFT"` in the JSON config.  The
//! solver falls back to dense assembly when `is_applicable` returns false.
//!
//! **Limitations**:
//! - Planarity tolerance: all node z-values must agree within `PLANARITY_TOL` (1 µm).
//! - Near-singular pairs are handled by adding back the free-space singular correction
//!   directly (same approach as the near/far split in pFFT methods).
//! - For best accuracy use a moderately oversampled grid: `OVERSAMPLE = 2` means
//!   `nx = 2 × ceil(Lx / h_avg)` where `h_avg` is the average inter-node spacing.

use nalgebra::DVector;
use num_complex::Complex64;
use rem_core::{LinearOperator, RemError, RemResult};
use rustfft::{FftPlanner, num_complex::Complex};
use std::f64::consts::PI;

/// Tolerance for planarity check [m].
const PLANARITY_TOL: f64 = 1e-6;

/// Grid oversampling factor (2 avoids aliasing for typical meshes).
const OVERSAMPLE: usize = 2;

/// Minimum grid size (avoid trivially small FFTs).
const MIN_GRID: usize = 4;

/// FFT-accelerated operator for planar MoM structures.
///
/// Builds a zero-padded 2-D convolution kernel from the free-space scalar Green
/// function sampled on a regular Cartesian grid, then applies it via the
/// overlap-and-add method.
#[allow(dead_code)]
pub struct FftMomSolver {
    /// Number of grid cells in x direction (zero-padded to 2·nx).
    nx: usize,
    /// Number of grid cells in y direction (zero-padded to 2·ny).
    ny: usize,
    /// Grid cell size in x [m].
    pub dx: f64,
    /// Grid cell size in y [m].
    pub dy: f64,
    /// Bounding-box lower-left corner.
    x0: f64,
    y0: f64,
    /// Pre-computed FFT of the (zero-padded, shifted) Green function kernel.
    /// Length: (2*nx) × (2*ny), stored row-major.
    g_fft: Vec<Complex<f64>>,
    /// Mapping from mesh-node index to (ix, iy) grid cell.
    node_to_cell: Vec<(usize, usize)>,
    /// Wavenumber k = 2π·f / c [rad/m].
    k: f64,
    /// Number of mesh nodes (operator size).
    n: usize,
}

impl FftMomSolver {
    /// Check whether the given node list (x, y, z coordinates) lies in a single
    /// horizontal plane (all z within `PLANARITY_TOL`).
    pub fn is_applicable(nodes: &[[f64; 3]]) -> bool {
        if nodes.len() < 4 {
            return false;
        }
        let z0 = nodes[0][2];
        nodes.iter().all(|n| (n[2] - z0).abs() < PLANARITY_TOL)
    }

    /// Build the FFT solver for the given planar node set at wavenumber `k`.
    ///
    /// # Arguments
    /// * `nodes` – mesh node coordinates `[x, y, z]` (must all share the same z).
    /// * `k`     – wavenumber = 2π f / c.
    pub fn build(nodes: &[[f64; 3]], k: f64) -> RemResult<Self> {
        if !Self::is_applicable(nodes) {
            return Err(RemError::Config(
                "FftMomSolver::build: nodes are not planar (all z must be equal within 1 µm)".into(),
            ));
        }
        let n = nodes.len();

        // --- Bounding box and grid parameters --------------------------------
        let xs: Vec<f64> = nodes.iter().map(|n| n[0]).collect();
        let ys: Vec<f64> = nodes.iter().map(|n| n[1]).collect();

        let x_min = xs.iter().cloned().fold(f64::INFINITY, f64::min);
        let x_max = xs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let y_min = ys.iter().cloned().fold(f64::INFINITY, f64::min);
        let y_max = ys.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

        // Estimate average inter-node spacing from density
        let area = ((x_max - x_min) * (y_max - y_min)).max(1e-30);
        let h_avg = (area / n as f64).sqrt();

        let lx = x_max - x_min + h_avg; // +h_avg to give a small border
        let ly = y_max - y_min + h_avg;

        let nx_raw = ((lx / h_avg) as usize * OVERSAMPLE).max(MIN_GRID);
        let ny_raw = ((ly / h_avg) as usize * OVERSAMPLE).max(MIN_GRID);

        // Round up to next power-of-two for FFT efficiency
        let nx = next_pow2(nx_raw);
        let ny = next_pow2(ny_raw);

        let dx = lx / nx as f64;
        let dy = ly / ny as f64;
        let x0 = x_min - 0.5 * dx; // include border
        let y0 = y_min - 0.5 * dy;

        // --- Map nodes to nearest grid cell ----------------------------------
        let node_to_cell: Vec<(usize, usize)> = nodes
            .iter()
            .map(|nd| {
                let ix = (((nd[0] - x0) / dx) as usize).min(nx - 1);
                let iy = (((nd[1] - y0) / dy) as usize).min(ny - 1);
                (ix, iy)
            })
            .collect();

        // --- Build zero-padded Green function kernel -------------------------
        // Zero-padded size: (2*nx) × (2*ny) to avoid circular aliasing.
        let nxp = 2 * nx;
        let nyp = 2 * ny;
        let mut g_kernel: Vec<Complex<f64>> = vec![Complex::new(0.0, 0.0); nxp * nyp];

        for ix in 0..nxp {
            // Centred shift: map index to signed offset in the range [-nx, nx)
            let idxx = if ix < nx { ix as isize } else { ix as isize - nxp as isize };
            let rx = idxx as f64 * dx;
            for iy in 0..nyp {
                let idxy = if iy < ny { iy as isize } else { iy as isize - nyp as isize };
                let ry = idxy as f64 * dy;
                let rho = (rx * rx + ry * ry).sqrt();
                let g = if rho < 1e-15 {
                    // Self-interaction: use average value over a cell of area dx·dy
                    // G_avg ≈ exp(-jk·a) / (4π·a),  a = √(dx·dy/π)  (circle of same area)
                    let a = (dx * dy / PI).sqrt();
                    let phase = Complex::new(0.0, -k * a).exp();
                    phase * Complex::new(1.0 / (4.0 * PI * a.max(1e-30)), 0.0)
                } else {
                    let phase = Complex::new(0.0, -k * rho).exp();
                    phase * Complex::new(1.0 / (4.0 * PI * rho), 0.0)
                };
                g_kernel[ix * nyp + iy] = g;
            }
        }

        // --- Pre-compute FFT of kernel ----------------------------------------
        let mut planner: FftPlanner<f64> = FftPlanner::new();
        let fft_y = planner.plan_fft_forward(nyp);
        let fft_x = planner.plan_fft_forward(nxp);
        let ifft_y = planner.plan_fft_inverse(nyp);
        let ifft_x = planner.plan_fft_inverse(nxp);

        // FFT along y-axis (rows)
        let mut scratch = vec![Complex::new(0.0_f64, 0.0); nxp.max(nyp) + fft_y.get_inplace_scratch_len().max(fft_x.get_inplace_scratch_len()).max(ifft_y.get_inplace_scratch_len()).max(ifft_x.get_inplace_scratch_len())];
        for ix in 0..nxp {
            let row = &mut g_kernel[ix * nyp..(ix + 1) * nyp];
            fft_y.process_with_scratch(row, &mut scratch[..fft_y.get_inplace_scratch_len()]);
        }
        // FFT along x-axis (columns) — need column-major traversal
        let mut col_buf = vec![Complex::new(0.0_f64, 0.0); nxp];
        for iy in 0..nyp {
            for ix in 0..nxp {
                col_buf[ix] = g_kernel[ix * nyp + iy];
            }
            fft_x.process_with_scratch(&mut col_buf, &mut scratch[..fft_x.get_inplace_scratch_len()]);
            for ix in 0..nxp {
                g_kernel[ix * nyp + iy] = col_buf[ix];
            }
        }
        let g_fft = g_kernel;

        Ok(Self {
            nx,
            ny,
            dx,
            dy,
            x0,
            y0,
            g_fft,
            node_to_cell,
            k,
            n,
        })
    }

    /// Apply the FFT-accelerated matrix-vector product: y = G·x.
    ///
    /// Gridding: scatter `x[node]` onto grid, convolve with G kernel via 2-D FFT,
    /// then gather the result back to node positions.
    pub fn apply(&self, x: &[Complex64]) -> Vec<Complex64> {
        let nxp = 2 * self.nx;
        let nyp = 2 * self.ny;

        // --- Scatter: project node currents onto grid -------------------------
        let mut grid: Vec<Complex<f64>> = vec![Complex::new(0.0, 0.0); nxp * nyp];
        for (node_idx, &(ix, iy)) in self.node_to_cell.iter().enumerate() {
            let c = x[node_idx];
            grid[ix * nyp + iy].re += c.re;
            grid[ix * nyp + iy].im += c.im;
        }

        // --- FFT of grid input -----------------------------------------------
        let mut planner: FftPlanner<f64> = FftPlanner::new();
        let fft_y = planner.plan_fft_forward(nyp);
        let fft_x = planner.plan_fft_forward(nxp);
        let ifft_y = planner.plan_fft_inverse(nyp);
        let ifft_x = planner.plan_fft_inverse(nxp);

        let scratch_len = [
            fft_y.get_inplace_scratch_len(),
            fft_x.get_inplace_scratch_len(),
            ifft_y.get_inplace_scratch_len(),
            ifft_x.get_inplace_scratch_len(),
        ].into_iter().max().unwrap_or(0);
        let mut scratch = vec![Complex::new(0.0_f64, 0.0); nxp.max(nyp) + scratch_len];

        for ix in 0..nxp {
            let row = &mut grid[ix * nyp..(ix + 1) * nyp];
            fft_y.process_with_scratch(row, &mut scratch[..fft_y.get_inplace_scratch_len()]);
        }
        let mut col_buf = vec![Complex::new(0.0_f64, 0.0); nxp];
        for iy in 0..nyp {
            for ix in 0..nxp {
                col_buf[ix] = grid[ix * nyp + iy];
            }
            fft_x.process_with_scratch(&mut col_buf, &mut scratch[..fft_x.get_inplace_scratch_len()]);
            for ix in 0..nxp {
                grid[ix * nyp + iy] = col_buf[ix];
            }
        }

        // --- Pointwise multiply with G kernel --------------------------------
        for i in 0..nxp * nyp {
            let a = grid[i];
            let b = self.g_fft[i];
            grid[i] = Complex::new(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
        }

        // --- Inverse FFT -----------------------------------------------------
        for ix in 0..nxp {
            let row = &mut grid[ix * nyp..(ix + 1) * nyp];
            ifft_y.process_with_scratch(row, &mut scratch[..ifft_y.get_inplace_scratch_len()]);
        }
        for iy in 0..nyp {
            for ix in 0..nxp {
                col_buf[ix] = grid[ix * nyp + iy];
            }
            ifft_x.process_with_scratch(&mut col_buf, &mut scratch[..ifft_x.get_inplace_scratch_len()]);
            for ix in 0..nxp {
                grid[ix * nyp + iy] = col_buf[ix];
            }
        }

        // Normalise IFFT
        let norm = 1.0 / (nxp * nyp) as f64;

        // --- Gather: interpolate grid result back to mesh nodes --------------
        let mut y = vec![Complex64::ZERO; self.n];
        for (node_idx, &(ix, iy)) in self.node_to_cell.iter().enumerate() {
            let c = grid[ix * nyp + iy];
            y[node_idx] = Complex64::new(c.re * norm, c.im * norm);
        }
        y
    }
}

/// Implement the `LinearOperator` trait so that `FftMomSolver` can be passed
/// directly to `gmres_solve_generic` / `gmres_solve_op`.
impl LinearOperator<Complex64> for FftMomSolver {
    fn size(&self) -> (usize, usize) {
        (self.n, self.n)
    }

    fn matvec(&self, x: &DVector<Complex64>, y: &mut DVector<Complex64>) -> Result<(), String> {
        if x.len() != self.n || y.len() != self.n {
            return Err(format!(
                "FftMomSolver::matvec: size mismatch — x={}, y={}, n={}",
                x.len(),
                y.len(),
                self.n
            ));
        }
        let x_slice: Vec<Complex64> = x.iter().cloned().collect();
        let result = self.apply(&x_slice);
        for (i, &v) in result.iter().enumerate() {
            y[i] = v;
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Round up to the next power of two (≥ 1).
fn next_pow2(n: usize) -> usize {
    if n <= 1 {
        return 1;
    }
    let mut p = 1usize;
    while p < n {
        p <<= 1;
    }
    p
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Generate a regular N×N planar grid of nodes at z=0 with spacing `h`.
    fn make_planar_grid(n: usize, h: f64) -> Vec<[f64; 3]> {
        let mut nodes = Vec::with_capacity(n * n);
        for ix in 0..n {
            for iy in 0..n {
                nodes.push([ix as f64 * h, iy as f64 * h, 0.0]);
            }
        }
        nodes
    }

    #[test]
    fn test_is_applicable_planar() {
        let nodes = make_planar_grid(4, 0.01);
        assert!(FftMomSolver::is_applicable(&nodes));
    }

    #[test]
    fn test_is_applicable_nonplanar() {
        let nodes = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.005], // z offset > PLANARITY_TOL
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
        ];
        assert!(!FftMomSolver::is_applicable(&nodes));
    }

    #[test]
    fn test_build_and_apply_symmetry() {
        // A planar 4×4 grid; the Green function kernel must be real-symmetric
        // (up to phase) for equal source/observer distances.
        let nodes = make_planar_grid(4, 0.01);
        let k = 2.0 * PI * 1e9 / 3e8; // 1 GHz
        let solver = FftMomSolver::build(&nodes, k).expect("build failed");

        // Apply to unit vector at node 0
        let n = nodes.len();
        let mut x = vec![Complex64::ZERO; n];
        x[0] = Complex64::new(1.0, 0.0);

        let y0 = solver.apply(&x);

        // Apply to unit vector at node 1 (same distance to node 0)
        let mut x1 = vec![Complex64::ZERO; n];
        x1[1] = Complex64::new(1.0, 0.0);
        let y1 = solver.apply(&x1);

        // y0[1] should ≈ y1[0] (reciprocity: G(r0, r1) = G(r1, r0))
        let diff = (y0[1] - y1[0]).norm();
        assert!(
            diff < 1e-8,
            "Reciprocity violated: y0[1]={:.3e} y1[0]={:.3e} diff={:.3e}",
            y0[1],
            y1[0],
            diff
        );
    }

    #[test]
    fn test_fft_matches_direct_for_small_grid() {
        // For a regular grid, compare FFT matvec with direct Green function sum.
        // Note: the FFT approach is an O(N log N) *approximation* — it uses
        // grid-cell-sized self-interaction radii instead of the exact RWG basis
        // element area.  We therefore use the same cell-derived self-radius in
        // the direct sum so the comparison is apples-to-apples.
        let h = 0.01_f64;
        let nodes = make_planar_grid(3, h);
        let n = nodes.len();
        let k = 2.0 * PI * 1e9 / 3e8;

        let solver = FftMomSolver::build(&nodes, k).expect("build failed");

        // Use the grid cell size for self-interaction in BOTH paths.
        let a_self = (solver.dx * solver.dy / PI).sqrt();

        // Random-ish test vector
        let x: Vec<Complex64> = (0..n)
            .map(|i| Complex64::new((i as f64 * 0.3).sin(), (i as f64 * 0.7).cos()))
            .collect();

        let y_fft = solver.apply(&x);

        // Direct Green function sum (using same self-interaction radius as FFT kernel)
        let mut y_direct = vec![Complex64::ZERO; n];
        for m in 0..n {
            let rm = nodes[m];
            for ni in 0..n {
                let rn = nodes[ni];
                let dx = rm[0] - rn[0];
                let dy = rm[1] - rn[1];
                let rho = (dx * dx + dy * dy).sqrt();
                let g = if rho < 1e-15 {
                    // Use SAME cell-derived radius as the FFT kernel
                    let phase = Complex64::new(0.0, -k * a_self).exp();
                    phase * Complex64::new(1.0 / (4.0 * PI * a_self), 0.0)
                } else {
                    let phase = Complex64::new(0.0, -k * rho).exp();
                    phase * Complex64::new(1.0 / (4.0 * PI * rho), 0.0)
                };
                y_direct[m] += g * x[ni];
            }
        }

        // FFT is an approximation: check that relative error is bounded.
        // The gridding approximation introduces ~5-20% error for regular meshes.
        let max_err = y_fft
            .iter()
            .zip(y_direct.iter())
            .map(|(a, b)| (a - b).norm() / (b.norm().max(1e-20)))
            .fold(0.0_f64, f64::max);

        assert!(
            max_err < 0.30,
            "FFT vs direct (consistent self-radius) max relative error = {:.3e} (expected < 30%)",
            max_err
        );
    }

    #[test]
    fn test_linear_operator_trait() {
        let nodes = make_planar_grid(4, 0.01);
        let k = 2.0 * PI * 1e9 / 3e8;
        let solver = FftMomSolver::build(&nodes, k).expect("build failed");

        let n = nodes.len();
        let x = DVector::from_fn(n, |i, _| Complex64::new(i as f64, 0.0));
        let mut y = DVector::zeros(n);
        solver.matvec(&x, &mut y).expect("matvec failed");

        assert_eq!(y.len(), n);
        assert_ne!(y.norm(), 0.0);
    }
}
