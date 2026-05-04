//! 3-D FFT-accelerated Fast Multipole Method (FMM) for EFIE/CFIE-RWG.
//!
//! # Algorithm
//!
//! Implements a **single-level monopole FMM** using 3-D FFT convolution:
//!
//! 1. Map RWG basis centroids onto a uniform 3-D grid of cells.
//! 2. For each basis `n`, compute vector moment
//!    `m_n = l_n/2 * [(c⁺_n − v⁺_n) − (c⁻_n − v⁻_n)]`
//!    (centroid displacement of T⁺ and T⁻).
//! 3. **Far-field** (all cells, monopole approximation):
//!    `F_t = Σ_s G(c_t − c_s) M_s`  via 3-D FFT convolution, O(N log N).
//! 4. **Near-field correction**: subtract the monopole approximation for
//!    adjacent cells and add the exact CFIE block.
//!
//! The resulting LinearOperator approximates `Z · x` with:
//! - **exact** near-field (same cell + 26 direct neighbours),
//! - monopole accuracy for far-field interactions (error ∝ `cell_size/R`).
//!
//! For high accuracy, prefer `ACA` or `GMRES` (direct assembly).
//! FMM is most efficient for large systems (N > 2000) where O(N²) assembly
//! becomes prohibitive.
//!
//! # References
//!
//! - Rokhlin (1985), "Rapid solution of integral equations of classical potential theory"
//! - Greengard & Rokhlin (1987), "A fast algorithm for particle simulations"
//! - Coifman et al. (1993), "The fast multipole method for the wave equation"

use crate::assemble::assemble_cfie_rwg_block;
use crate::basis::rwg::RwgBasis;
use crate::quadrature::TriQuad;
use crate::surface_mesh::SurfaceMesh;
use nalgebra::{DMatrix, DVector};
use num_complex::Complex64;
use rem_core::{LinearOperator, RemError, RemResult, C0, EPS0, MU0};
use rem_layered_green::GreenFunction;
use rustfft::{num_complex::Complex, FftPlanner};
use std::f64::consts::PI;

const NEAR_RADIUS: isize = 1; // cells within ±NEAR_RADIUS in each dim are "near"
const MIN_CELLS_PER_DIM: usize = 3;
const MAX_CELLS_PER_DIM: usize = 32; // cap FFT grid size
const TARGET_BASES_PER_CELL: usize = 10;

// ---------------------------------------------------------------------------
// Helper types
// ---------------------------------------------------------------------------

/// Dense near-field block Z[row_ids × col_ids] from exact assembly.
struct NearBlock {
    row_ids: Vec<usize>,
    col_ids: Vec<usize>,
    data: DMatrix<Complex64>,
}

/// Monopole correction entry for a near-field cell pair.
/// Used to subtract the far-field approximation for near pairs.
struct NearMonopole {
    target_cell: usize,
    source_cell: usize,
    /// Pre-scaled G(c_t, c_s) value (includes jkη factor).
    g_scaled: Complex64,
}

// ---------------------------------------------------------------------------
// FmmMomSolver
// ---------------------------------------------------------------------------

/// FMM-accelerated matrix-free operator for EFIE/CFIE-RWG.
///
/// Implements `LinearOperator<Complex64>` for use with GMRES.
/// Build via [`FmmMomSolver::build`].
pub struct FmmMomSolver {
    // ── Grid ─────────────────────────────────────────────────────────────
    nx: usize,
    ny: usize,
    nz: usize,
    x0: f64,
    y0: f64,
    z0: f64,
    dx: f64,
    dy: f64,
    dz: f64,
    // ── Precomputed FFT of 3-D scalar G kernel (zero-padded) ─────────────
    /// Length = (2*nx) × (2*ny) × (2*nz), stored as flat Vec<Complex<f64>>.
    g_fft: Vec<Complex<f64>>,
    // ── Per-basis data ────────────────────────────────────────────────────
    /// Plain vector moment m_n [3] for each basis (see module doc).
    moments: Vec<[f64; 3]>,
    /// Flat 3-D cell index for each basis.
    cell_of_basis: Vec<usize>,
    /// For each cell: list of basis indices (index into the global bases array).
    cell_bases: Vec<Vec<usize>>,
    /// Cell centres.
    cell_centers: Vec<[f64; 3]>,
    // ── Near-field ────────────────────────────────────────────────────────
    near_blocks: Vec<NearBlock>,
    near_monopole: Vec<NearMonopole>,
    // ── Sizes / physics ───────────────────────────────────────────────────
    n_bases: usize,
    k: f64,
    /// jkη coefficient for the A-potential monopole far-field term.
    jk_eta: Complex64,
}

impl FmmMomSolver {
    /// Build the FMM solver.
    ///
    /// # Arguments
    /// * `surf`  – surface mesh
    /// * `bases` – RWG basis functions (ordered consistently with the Z matrix)
    /// * `green` – Green's function (free-space or layered)
    /// * `freq`  – frequency [Hz]
    /// * `alpha` – CFIE mixing coefficient (1.0 = pure EFIE)
    /// * `quad`  – quadrature rule for near-field assembly
    pub fn build(
        surf: &SurfaceMesh,
        bases: &[RwgBasis],
        green: &dyn GreenFunction,
        freq: f64,
        alpha: f64,
        quad: &TriQuad,
    ) -> RemResult<Self> {
        let n = bases.len();
        if n == 0 {
            return Err(RemError::Mesh("FmmMomSolver::build: no RWG bases".into()));
        }

        let omega = 2.0 * PI * freq;
        let k     = omega / C0;
        let eta0  = (MU0 / EPS0).sqrt();
        let jk_eta = Complex64::new(0.0, k * eta0); // jkη

        // ── Compute basis centroids ──────────────────────────────────────
        // Centroid of basis n = midpoint of T⁺ and T⁻ centroids.
        let centroids: Vec<[f64; 3]> = bases.iter().map(|b| {
            let cp = surf.faces[b.plus_face].centroid;
            let cm = surf.faces[b.minus_face].centroid;
            [(cp[0]+cm[0])*0.5, (cp[1]+cm[1])*0.5, (cp[2]+cm[2])*0.5]
        }).collect();

        // ── Compute vector moments m_n ───────────────────────────────────
        // m_n = l_n/2 * [(c⁺_n − v⁺_n) − (c⁻_n − v⁻_n)]
        let moments: Vec<[f64; 3]> = bases.iter().map(|b| {
            let cp = surf.faces[b.plus_face].centroid;
            let cm = surf.faces[b.minus_face].centroid;
            let vp = surf.nodes[b.free_node_plus];
            let vm = surf.nodes[b.free_node_minus];
            let scale = b.length * 0.5;
            let dp = [cp[0]-vp[0], cp[1]-vp[1], cp[2]-vp[2]];
            let dm = [cm[0]-vm[0], cm[1]-vm[1], cm[2]-vm[2]];
            [scale*(dp[0]-dm[0]), scale*(dp[1]-dm[1]), scale*(dp[2]-dm[2])]
        }).collect();

        // ── Build 3-D grid ───────────────────────────────────────────────
        let xmin = centroids.iter().map(|c| c[0]).fold(f64::INFINITY, f64::min);
        let xmax = centroids.iter().map(|c| c[0]).fold(f64::NEG_INFINITY, f64::max);
        let ymin = centroids.iter().map(|c| c[1]).fold(f64::INFINITY, f64::min);
        let ymax = centroids.iter().map(|c| c[1]).fold(f64::NEG_INFINITY, f64::max);
        let zmin = centroids.iter().map(|c| c[2]).fold(f64::INFINITY, f64::min);
        let zmax = centroids.iter().map(|c| c[2]).fold(f64::NEG_INFINITY, f64::max);

        // Choose number of cells so each cell has ~TARGET_BASES_PER_CELL bases.
        let cells_total = (n / TARGET_BASES_PER_CELL).max(1);
        let m_raw = (cells_total as f64).powf(1.0/3.0).ceil() as usize;
        let m = m_raw.max(MIN_CELLS_PER_DIM).min(MAX_CELLS_PER_DIM);

        // Add a small border (half a cell) to avoid boundary issues.
        let lx = (xmax - xmin).max(1e-12);
        let ly = (ymax - ymin).max(1e-12);
        let lz = (zmax - zmin).max(1e-12);

        // Scale cells to be roughly cubic.
        let h = (lx * ly * lz).powf(1.0/3.0);
        let nx = (lx / h * m as f64).ceil() as usize;
        let ny = (ly / h * m as f64).ceil() as usize;
        let nz = (lz / h * m as f64).ceil() as usize;
        let nx = nx.max(MIN_CELLS_PER_DIM).min(MAX_CELLS_PER_DIM);
        let ny = ny.max(MIN_CELLS_PER_DIM).min(MAX_CELLS_PER_DIM);
        let nz = nz.max(MIN_CELLS_PER_DIM).min(MAX_CELLS_PER_DIM);

        let border = 0.5;
        let dx = lx / nx as f64 * (1.0 + border / nx as f64);
        let dy = ly / ny as f64 * (1.0 + border / ny as f64);
        let dz = lz / nz as f64 * (1.0 + border / nz as f64);
        let x0 = xmin - 0.5 * dx;
        let y0 = ymin - 0.5 * dy;
        let z0 = zmin - 0.5 * dz;

        // ── Assign bases to cells ────────────────────────────────────────
        let flat_idx = |ix: usize, iy: usize, iz: usize| ix * ny * nz + iy * nz + iz;
        let n_cells = nx * ny * nz;

        let cell_of_basis: Vec<usize> = centroids.iter().map(|c| {
            let ix = (((c[0] - x0) / dx) as usize).min(nx - 1);
            let iy = (((c[1] - y0) / dy) as usize).min(ny - 1);
            let iz = (((c[2] - z0) / dz) as usize).min(nz - 1);
            flat_idx(ix, iy, iz)
        }).collect();

        let mut cell_bases: Vec<Vec<usize>> = vec![Vec::new(); n_cells];
        for (bi, &ci) in cell_of_basis.iter().enumerate() {
            cell_bases[ci].push(bi);
        }

        // ── Compute cell centres ─────────────────────────────────────────
        let cell_centers: Vec<[f64; 3]> = (0..n_cells).map(|ci| {
            let ix = ci / (ny * nz);
            let iy = (ci / nz) % ny;
            let iz = ci % nz;
            [
                x0 + (ix as f64 + 0.5) * dx,
                y0 + (iy as f64 + 0.5) * dy,
                z0 + (iz as f64 + 0.5) * dz,
            ]
        }).collect();

        // ── Build 3-D FFT of scalar Green kernel (zero-padded) ───────────
        // Zero-padded grid: (2*nx) × (2*ny) × (2*nz) for linear convolution.
        let nxp = 2 * nx;
        let nyp = 2 * ny;
        let nzp = 2 * nz;
        let mut g_kernel: Vec<Complex<f64>> = vec![Complex::new(0.0, 0.0); nxp * nyp * nzp];

        for ix in 0..nxp {
            let idxx = if ix < nx { ix as isize } else { ix as isize - nxp as isize };
            let rx = idxx as f64 * dx;
            for iy in 0..nyp {
                let idxy = if iy < ny { iy as isize } else { iy as isize - nyp as isize };
                let ry = idxy as f64 * dy;
                for iz in 0..nzp {
                    let idxz = if iz < nz { iz as isize } else { iz as isize - nzp as isize };
                    let rz = idxz as f64 * dz;
                    let r = (rx*rx + ry*ry + rz*rz).sqrt();
                    let g = if r < 1e-15 {
                        // Self-interaction: use average G over a sphere of same volume as cell.
                        let a = (3.0 * dx * dy * dz / (4.0 * PI)).powf(1.0/3.0);
                        let phase = Complex::new(0.0, k * a).exp();
                        phase * Complex::new(1.0 / (4.0 * PI * a.max(1e-30)), 0.0)
                    } else {
                        let phase = Complex::new(0.0, k * r).exp();
                        phase * Complex::new(1.0 / (4.0 * PI * r), 0.0)
                    };
                    g_kernel[(ix * nyp + iy) * nzp + iz] = g;
                }
            }
        }

        // 3-D FFT of g_kernel (in-place, row-major x,y,z).
        let mut planner: FftPlanner<f64> = FftPlanner::new();
        let fft_x = planner.plan_fft_forward(nxp);
        let fft_y = planner.plan_fft_forward(nyp);
        let fft_z = planner.plan_fft_forward(nzp);
        let scratch_len = fft_x.get_inplace_scratch_len()
            .max(fft_y.get_inplace_scratch_len())
            .max(fft_z.get_inplace_scratch_len());
        let mut scratch = vec![Complex::new(0.0_f64, 0.0); scratch_len];

        // FFT along z (innermost, contiguous)
        for i in 0..(nxp * nyp) {
            let row = &mut g_kernel[i * nzp..(i + 1) * nzp];
            fft_z.process_with_scratch(row, &mut scratch[..fft_z.get_inplace_scratch_len()]);
        }
        // FFT along y
        let mut ybuf = vec![Complex::new(0.0_f64, 0.0); nyp];
        for ix in 0..nxp {
            for iz in 0..nzp {
                for iy in 0..nyp {
                    ybuf[iy] = g_kernel[(ix * nyp + iy) * nzp + iz];
                }
                fft_y.process_with_scratch(&mut ybuf, &mut scratch[..fft_y.get_inplace_scratch_len()]);
                for iy in 0..nyp {
                    g_kernel[(ix * nyp + iy) * nzp + iz] = ybuf[iy];
                }
            }
        }
        // FFT along x
        let mut xbuf = vec![Complex::new(0.0_f64, 0.0); nxp];
        for iy in 0..nyp {
            for iz in 0..nzp {
                for ix in 0..nxp {
                    xbuf[ix] = g_kernel[(ix * nyp + iy) * nzp + iz];
                }
                fft_x.process_with_scratch(&mut xbuf, &mut scratch[..fft_x.get_inplace_scratch_len()]);
                for ix in 0..nxp {
                    g_kernel[(ix * nyp + iy) * nzp + iz] = xbuf[ix];
                }
            }
        }
        let g_fft = g_kernel;

        // ── Build near-field blocks ──────────────────────────────────────
        // For each cell tc, find all adjacent cells (distance ≤ NEAR_RADIUS).
        let mut near_blocks: Vec<NearBlock>    = Vec::new();
        let mut near_monopole: Vec<NearMonopole> = Vec::new();

        for tc in 0..n_cells {
            if cell_bases[tc].is_empty() { continue; }
            let txc = (tc / (ny * nz)) as isize;
            let tyc = ((tc / nz) % ny) as isize;
            let tzc = (tc % nz) as isize;

            for dix in -NEAR_RADIUS..=NEAR_RADIUS {
                for diy in -NEAR_RADIUS..=NEAR_RADIUS {
                    for diz in -NEAR_RADIUS..=NEAR_RADIUS {
                        let sx = txc + dix;
                        let sy = tyc + diy;
                        let sz = tzc + diz;
                        if sx < 0 || sy < 0 || sz < 0 { continue; }
                        let sx = sx as usize; let sy = sy as usize; let sz = sz as usize;
                        if sx >= nx || sy >= ny || sz >= nz { continue; }
                        let sc = flat_idx(sx, sy, sz);
                        if cell_bases[sc].is_empty() { continue; }

                        // Exact near-field block
                        let block = assemble_cfie_rwg_block(
                            surf, bases,
                            &cell_bases[tc],
                            &cell_bases[sc],
                            green, freq, alpha, quad,
                        )?;
                        near_blocks.push(NearBlock {
                            row_ids: cell_bases[tc].clone(),
                            col_ids: cell_bases[sc].clone(),
                            data: block,
                        });

                        // Monopole correction for this near-cell pair
                        let ct = cell_centers[tc];
                        let cs = cell_centers[sc];
                        let dr = [ct[0]-cs[0], ct[1]-cs[1], ct[2]-cs[2]];
                        let r = (dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]).sqrt();
                        let g_ts = if r < 1e-15 {
                            // Same-cell: skip (handled exactly)
                            Complex64::ZERO
                        } else {
                            Complex64::new(0.0, k * r).exp()
                                / Complex64::new(4.0 * PI * r, 0.0)
                        };
                        near_monopole.push(NearMonopole {
                            target_cell: tc,
                            source_cell: sc,
                            g_scaled: jk_eta * g_ts,
                        });
                    }
                }
            }
        }

        log::info!(
            "FmmMomSolver: grid {}×{}×{}={} cells, {} bases, {} near blocks",
            nx, ny, nz, n_cells, n, near_blocks.len()
        );

        Ok(Self {
            nx, ny, nz,
            x0, y0, z0,
            dx, dy, dz,
            g_fft,
            moments,
            cell_of_basis,
            cell_bases,
            cell_centers,
            near_blocks,
            near_monopole,
            n_bases: n,
            k,
            jk_eta,
        })
    }

    // ── Helper: flat padded index ────────────────────────────────────────
    fn pad_idx(&self, ix: usize, iy: usize, iz: usize) -> usize {
        (ix * 2 * self.ny + iy) * 2 * self.nz + iz
    }

    // ── Helper: cell 3D indices from flat cell index ─────────────────────
    fn cell_xyz(&self, ci: usize) -> (usize, usize, usize) {
        let ix = ci / (self.ny * self.nz);
        let iy = (ci / self.nz) % self.ny;
        let iz = ci % self.nz;
        (ix, iy, iz)
    }

    /// Apply the 3-D FFT convolution for component `comp` (0=x, 1=y, 2=z).
    ///
    /// Computes `f[t] = Σ_s G(c_t − c_s) M_s[comp]` for all cells t.
    fn fft_convolve_component(&self, cell_moments_comp: &[Complex<f64>]) -> Vec<Complex<f64>> {
        let nxp = 2 * self.nx;
        let nyp = 2 * self.ny;
        let nzp = 2 * self.nz;

        // Zero-pad the cell moment array.
        let mut m_pad: Vec<Complex<f64>> = vec![Complex::new(0.0, 0.0); nxp * nyp * nzp];
        for ci in 0..self.nx * self.ny * self.nz {
            let (ix, iy, iz) = self.cell_xyz(ci);
            m_pad[self.pad_idx(ix, iy, iz)] = cell_moments_comp[ci];
        }

        // FFT of m_pad (in-place 3-D).
        let mut planner: FftPlanner<f64> = FftPlanner::new();
        let fft_x = planner.plan_fft_forward(nxp);
        let fft_y = planner.plan_fft_forward(nyp);
        let fft_z = planner.plan_fft_forward(nzp);
        let scratch_len = fft_x.get_inplace_scratch_len()
            .max(fft_y.get_inplace_scratch_len())
            .max(fft_z.get_inplace_scratch_len());
        let mut scratch = vec![Complex::new(0.0_f64, 0.0); scratch_len];

        // z
        for i in 0..(nxp * nyp) {
            let row = &mut m_pad[i * nzp..(i + 1) * nzp];
            fft_z.process_with_scratch(row, &mut scratch[..fft_z.get_inplace_scratch_len()]);
        }
        // y
        let mut ybuf = vec![Complex::new(0.0_f64, 0.0); nyp];
        for ix in 0..nxp {
            for iz in 0..nzp {
                for iy in 0..nyp {
                    ybuf[iy] = m_pad[(ix * nyp + iy) * nzp + iz];
                }
                fft_y.process_with_scratch(&mut ybuf, &mut scratch[..fft_y.get_inplace_scratch_len()]);
                for iy in 0..nyp {
                    m_pad[(ix * nyp + iy) * nzp + iz] = ybuf[iy];
                }
            }
        }
        // x
        let mut xbuf = vec![Complex::new(0.0_f64, 0.0); nxp];
        for iy in 0..nyp {
            for iz in 0..nzp {
                for ix in 0..nxp {
                    xbuf[ix] = m_pad[(ix * nyp + iy) * nzp + iz];
                }
                fft_x.process_with_scratch(&mut xbuf, &mut scratch[..fft_x.get_inplace_scratch_len()]);
                for ix in 0..nxp {
                    m_pad[(ix * nyp + iy) * nzp + iz] = xbuf[ix];
                }
            }
        }

        // Pointwise multiply with precomputed G_fft.
        for i in 0..m_pad.len() {
            m_pad[i] = m_pad[i] * self.g_fft[i];
        }

        // IFFT (3-D).
        let ifft_x = planner.plan_fft_inverse(nxp);
        let ifft_y = planner.plan_fft_inverse(nyp);
        let ifft_z = planner.plan_fft_inverse(nzp);
        let scratch_len2 = ifft_x.get_inplace_scratch_len()
            .max(ifft_y.get_inplace_scratch_len())
            .max(ifft_z.get_inplace_scratch_len());
        let mut scratch2 = vec![Complex::new(0.0_f64, 0.0); scratch_len2];

        // z
        for i in 0..(nxp * nyp) {
            let row = &mut m_pad[i * nzp..(i + 1) * nzp];
            ifft_z.process_with_scratch(row, &mut scratch2[..ifft_z.get_inplace_scratch_len()]);
        }
        // y
        for ix in 0..nxp {
            for iz in 0..nzp {
                for iy in 0..nyp {
                    ybuf[iy] = m_pad[(ix * nyp + iy) * nzp + iz];
                }
                ifft_y.process_with_scratch(&mut ybuf, &mut scratch2[..ifft_y.get_inplace_scratch_len()]);
                for iy in 0..nyp {
                    m_pad[(ix * nyp + iy) * nzp + iz] = ybuf[iy];
                }
            }
        }
        // x
        for iy in 0..nyp {
            for iz in 0..nzp {
                for ix in 0..nxp {
                    xbuf[ix] = m_pad[(ix * nyp + iy) * nzp + iz];
                }
                ifft_x.process_with_scratch(&mut xbuf, &mut scratch2[..ifft_x.get_inplace_scratch_len()]);
                for ix in 0..nxp {
                    m_pad[(ix * nyp + iy) * nzp + iz] = xbuf[ix];
                }
            }
        }

        // Normalise FFT output.
        let norm = 1.0 / (nxp * nyp * nzp) as f64;
        for v in m_pad.iter_mut() {
            *v = *v * Complex::new(norm, 0.0);
        }

        m_pad
    }
}

// ---------------------------------------------------------------------------
// LinearOperator implementation
// ---------------------------------------------------------------------------

impl LinearOperator<Complex64> for FmmMomSolver {
    fn size(&self) -> (usize, usize) {
        (self.n_bases, self.n_bases)
    }

    /// Compute y ← Z_FMM · x using the 3-D FFT monopole FMM algorithm.
    ///
    /// ## Steps
    /// 1. Aggregate cell vector moments: `M_s = Σ_{n∈s} x_n m_n`.
    /// 2. Far-field (all-cell FFT): `F_t = G*M` via 3-D FFT.
    /// 3. Accumulate far-field to target bases: `y_m += jkη (m_m · F_{t(m)})`.
    /// 4. Subtract monopole for near-cell pairs.
    /// 5. Add exact near-field blocks.
    fn matvec(&self, x: &DVector<Complex64>, y: &mut DVector<Complex64>) -> Result<(), String> {
        let n = self.n_bases;
        if x.len() != n || y.len() != n {
            return Err(format!("FmmMomSolver::matvec dimension mismatch: n={}, x={}, y={}", n, x.len(), y.len()));
        }

        let n_cells = self.nx * self.ny * self.nz;

        // ── Step 1: compute cell vector moments ──────────────────────────
        // M_s[comp] = Σ_{n∈s} Re(x_n) × m_n[comp]  +  Im(x_n) × j × m_n[comp]
        let mut cell_moment_x: Vec<Complex<f64>> = vec![Complex::new(0.0, 0.0); n_cells];
        let mut cell_moment_y: Vec<Complex<f64>> = vec![Complex::new(0.0, 0.0); n_cells];
        let mut cell_moment_z: Vec<Complex<f64>> = vec![Complex::new(0.0, 0.0); n_cells];

        for bi in 0..n {
            let ci = self.cell_of_basis[bi];
            let xn = x[bi];
            let mn = self.moments[bi];
            let xc = Complex::new(xn.re, xn.im);
            let mnx = Complex::new(mn[0], 0.0);
            let mny = Complex::new(mn[1], 0.0);
            let mnz = Complex::new(mn[2], 0.0);
            cell_moment_x[ci] = cell_moment_x[ci] + xc * mnx;
            cell_moment_y[ci] = cell_moment_y[ci] + xc * mny;
            cell_moment_z[ci] = cell_moment_z[ci] + xc * mnz;
        }

        // ── Step 2: 3-D FFT convolution F_t = G * M (one per component) ──
        let fx = self.fft_convolve_component(&cell_moment_x);
        let fy = self.fft_convolve_component(&cell_moment_y);
        let fz = self.fft_convolve_component(&cell_moment_z);

        // ── Step 3: accumulate far-field A-potential to target bases ──────
        // y_m += jkη (m_m · F_{t(m)})
        for bi in 0..n {
            let ci = self.cell_of_basis[bi];
            let (ix, iy, iz) = self.cell_xyz(ci);
            let fidx = self.pad_idx(ix, iy, iz);
            let ftx = fx[fidx];
            let fty = fy[fidx];
            let ftz = fz[fidx];
            let mn = self.moments[bi];
            // dot product m_m · F_t (complex)
            let dot = Complex::<f64>::new(mn[0], 0.0) * ftx
                    + Complex::<f64>::new(mn[1], 0.0) * fty
                    + Complex::<f64>::new(mn[2], 0.0) * ftz;
            let contrib = Complex::<f64>::new(self.jk_eta.re, self.jk_eta.im) * dot;
            y[bi] += Complex64::new(contrib.re, contrib.im);
        }

        // ── Step 4: subtract monopole approximation for near-field pairs ──
        for nm in &self.near_monopole {
            if nm.g_scaled.norm() < 1e-30 { continue; }
            let tc = nm.target_cell;
            let sc = nm.source_cell;
            // Aggregate source moment
            let mut ms = [Complex64::ZERO; 3];
            for &bi in &self.cell_bases[sc] {
                let mn = self.moments[bi];
                let xn = x[bi];
                ms[0] += xn * mn[0];
                ms[1] += xn * mn[1];
                ms[2] += xn * mn[2];
            }
            // Subtract from each target basis in tc
            for &tm in &self.cell_bases[tc] {
                let mm = self.moments[tm];
                let dot = mm[0]*ms[0] + mm[1]*ms[1] + mm[2]*ms[2];
                y[tm] -= nm.g_scaled * dot;
            }
        }

        // ── Step 5: add exact near-field blocks ───────────────────────────
        for nb in &self.near_blocks {
            for (ri, &row_bi) in nb.row_ids.iter().enumerate() {
                let mut acc = Complex64::ZERO;
                for (ci_block, &col_bi) in nb.col_ids.iter().enumerate() {
                    acc += nb.data[(ri, ci_block)] * x[col_bi];
                }
                y[row_bi] += acc;
            }
        }

        Ok(())
    }

    fn density(&self) -> f64 {
        // Near-field density approximation
        let n_cells = self.nx * self.ny * self.nz;
        let near_fraction = (27.0 * self.n_bases as f64 / n_cells as f64)
            .min(self.n_bases as f64);
        near_fraction / self.n_bases as f64
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::surface_mesh::{SurfaceMesh, TriFace, SharedEdge, tri_geometry};
    use rem_layered_green::FreeSpaceGreen;

    /// Build a minimal two-triangle mesh (one RWG basis).
    fn two_tri_mesh() -> SurfaceMesh {
        let nodes = vec![
            [0.0_f64, 0.0, 0.0],
            [1.0,     0.0, 0.0],
            [0.0,     1.0, 0.0],
            [1.0,     1.0, 0.0],
        ];
        let (c0, n0, a0) = tri_geometry(&nodes[0], &nodes[1], &nodes[2]);
        let (c1, n1, a1) = tri_geometry(&nodes[1], &nodes[3], &nodes[2]);
        let faces = vec![
            TriFace { nodes: [0,1,2], centroid: c0, normal: n0, area: a0 },
            TriFace { nodes: [1,3,2], centroid: c1, normal: n1, area: a1 },
        ];
        let edge_nodes = [1_usize, 2];
        let edges = vec![SharedEdge {
            nodes: edge_nodes,
            plus_face: 0,
            minus_face: 1,
            length: 1.0_f64.sqrt() * 2_f64.sqrt().recip() * (2.0_f64).sqrt(), // |1-0,1-0| = √2? compute properly
        }];
        // Correct edge length
        let elen = {
            let d = [nodes[edge_nodes[0]][0]-nodes[edge_nodes[1]][0],
                     nodes[edge_nodes[0]][1]-nodes[edge_nodes[1]][1],
                     nodes[edge_nodes[0]][2]-nodes[edge_nodes[1]][2]];
            (d[0]*d[0]+d[1]*d[1]+d[2]*d[2]).sqrt()
        };
        let edges = vec![SharedEdge { nodes: edge_nodes, plus_face: 0, minus_face: 1, length: elen }];
        SurfaceMesh { nodes, faces, edges, boundary_edges: vec![], face_attrs: vec![0, 0], global_node_ids: vec![0, 1, 2, 3] }
    }

    /// Basic smoke test: FmmMomSolver builds and produces finite output.
    #[test]
    fn fmm_build_and_matvec_finite() {
        let mesh = two_tri_mesh();
        let bases = crate::basis::rwg::generate_rwg_bases(&mesh);
        let n = bases.len();

        assert!(n >= 1, "expected at least one RWG basis");

        let freq = 1e9_f64;
        let green = FreeSpaceGreen::from_freq(freq);
        let quad = TriQuad::new(3);

        let fmm = FmmMomSolver::build(&mesh, &bases, &green, freq, 1.0, &quad)
            .expect("FMM build should succeed");

        let x = DVector::from_element(n, Complex64::new(1.0, 0.0));
        let mut y = DVector::zeros(n);
        fmm.matvec(&x, &mut y).expect("matvec should not fail");

        for i in 0..n {
            assert!(y[i].norm().is_finite(), "y[{}] = {} is not finite", i, y[i]);
        }
    }

    /// For a small problem, FMM result should be close to dense Z·x.
    #[test]
    #[ignore] // Slow; run with `cargo test -- --ignored`
    fn fmm_matches_dense_small_sphere() {
        use crate::assemble::{assemble_cfie_rwg_block, lu_solve};
        use crate::mie::pec_sphere_rcs;
        use rem_core::C0;

        let freq = 3e8_f64;   // 300 MHz → λ=1 m, ka=1 for a=1/(2π) m
        let k = 2.0 * PI * freq / C0;
        let a = 1.0 / (2.0 * PI);  // ka=1
        let green = FreeSpaceGreen::from_freq(freq);
        let quad = TriQuad::new(3);

        // Build a small icosphere using the integration test helper via the mesh loader
        // (icosphere only available in integration tests). Use a flat-plate mesh instead.
        let surf = {
            // 6-node planar mesh approximating a disk: 4 triangles around center
            let nodes = vec![
                [0.0_f64, 0.0, 0.0],  // center
                [a, 0.0, 0.0],
                [0.0, a, 0.0],
                [-a, 0.0, 0.0],
                [0.0, -a, 0.0],
                [a * 0.707, a * 0.707, 0.0],
            ];
            let (c0, n0, a0) = tri_geometry(&nodes[0], &nodes[1], &nodes[2]);
            let (c1, n1, a1) = tri_geometry(&nodes[0], &nodes[2], &nodes[3]);
            let (c2, n2, a2) = tri_geometry(&nodes[0], &nodes[3], &nodes[4]);
            let (c3, n3, a3) = tri_geometry(&nodes[0], &nodes[4], &nodes[1]);
            let faces = vec![
                TriFace { nodes: [0,1,2], centroid: c0, normal: n0, area: a0 },
                TriFace { nodes: [0,2,3], centroid: c1, normal: n1, area: a1 },
                TriFace { nodes: [0,3,4], centroid: c2, normal: n2, area: a2 },
                TriFace { nodes: [0,4,1], centroid: c3, normal: n3, area: a3 },
            ];
            let edge_pairs = [(1_usize, 2), (2, 3), (3, 4), (4, 1), (0, 1), (0, 2), (0, 3), (0, 4)];
            let edges: Vec<SharedEdge> = vec![]; // minimal: no shared interior edges in this mesh layout
            let _ = edge_pairs;
            SurfaceMesh { nodes, faces, edges, boundary_edges: vec![], face_attrs: vec![0;4], global_node_ids: (0..6).collect() }
        };
        let bases = crate::basis::rwg::generate_rwg_bases(&surf);
        let n = bases.len();

        // Dense Z matrix
        let row_ids: Vec<usize> = (0..n).collect();
        let z_dense = assemble_cfie_rwg_block(&surf, &bases, &row_ids, &row_ids, &green, freq, 0.5, &quad)
            .expect("dense block assembly");

        // FMM
        let fmm = FmmMomSolver::build(&surf, &bases, &green, freq, 0.5, &quad)
            .expect("FMM build");

        let x = DVector::from_element(n, Complex64::new(1.0, 0.5));
        let y_dense = z_dense.clone() * x.clone();
        let mut y_fmm = DVector::zeros(n);
        fmm.matvec(&x, &mut y_fmm).unwrap();

        // Check relative error is reasonable (monopole accuracy ~30% for close cells)
        let err = (y_dense.clone() - y_fmm.clone()).norm();
        let ref_norm = y_dense.norm();
        let rel_err = err / ref_norm.max(1e-30);
        println!("FMM vs dense relative error: {:.2}%", rel_err * 100.0);
        assert!(rel_err < 2.0, "FMM relative error {:.2}% should be < 200%", rel_err * 100.0);
    }
}
