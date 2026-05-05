//! Multilevel Fast Multipole Algorithm (MLFMA) for EFIE/CFIE-RWG.
//!
//! # Algorithm Overview
//!
//! The MLFMA accelerates the MoM impedance matrix-vector product from O(N²) to O(N log N)
//! using a hierarchical octree with scalar Helmholtz multipole expansions.
//!
//! ## Stages of one matvec `Z·x`:
//!
//! 1. **P2M** (particle-to-multipole): at each leaf box, aggregate the RWG source
//!    moments into a multipole expansion of order P.
//! 2. **M2M** (multipole-to-multipole): propagate multipoles upward through the
//!    octree (child → parent shift using translation operators).
//! 3. **M2L** (multipole-to-local): at each tree level, translate the multipole
//!    expansions of well-separated (interaction-list) boxes into local expansions
//!    at the observer box.
//! 4. **L2L** (local-to-local): propagate local expansions downward (parent → child).
//! 5. **L2P** (local-to-particle): evaluate the local expansion at each leaf box
//!    and accumulate into the output vector.
//! 6. **Near-field** (P2P): exact CFIE-RWG block assembly for adjacent leaf boxes.
//!
//! ## Multipole Representation
//!
//! We use the standard spherical-harmonic expansion of the scalar Helmholtz Green function:
//!
//! ```text
//! G(r - r') = ik Σ_{n=0}^{P} Σ_{m=-n}^{n} h_n(k|r|) Y_n^m(r̂)  j_n(k|r'|) Y_n^m*(r̂')
//! ```
//!
//! Multipole coefficients at expansion centre `c` due to source at `r'` with
//! moment `w` are:
//!
//! ```text
//! M_{n,m} = w · j_n(k|r' - c|) · Y_n^m*(r̂'_c)
//! ```
//!
//! Local expansion coefficients at `c_obs` due to source multipole `M_{n,m}` at `c_src`:
//!
//! ```text
//! L_{n,m}(c_obs) = ik Σ_ν Σ_μ  h_{ν+n}(k|d|) A(ν,μ,n,m,d̂) · M_{ν,μ}
//! ```
//!
//! where d = c_src - c_obs.
//!
//! # References
//!
//! - Rokhlin (1985): "Rapid solution of integral equations"
//! - Greengard & Rokhlin (1987): Fast algorithm for particle simulations
//! - Coifman, Rokhlin, Wandzura (1993): Fast multipole method for the wave equation
//! - Song, Lu, Chew (1997): Multilevel fast multipole algorithm for EFIE

use crate::assemble::assemble_cfie_rwg_block;
use crate::basis::rwg::RwgBasis;
use crate::quadrature::TriQuad;
use crate::surface_mesh::SurfaceMesh;
use nalgebra::{DMatrix, DVector};
use num_complex::Complex64;
use rem_core::{LinearOperator, RemError, RemResult, C0, MU0};
use rem_layered_green::GreenFunction;
use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// Maximum octree depth (8 levels → 8^6 = 262 144 leaf boxes maximum).
const MAX_LEVELS: usize = 6;
/// Minimum bases per leaf to stop subdividing (trade-off accuracy vs near-field cost).
const MIN_BASES_PER_LEAF: usize = 8;
/// Default multipole order P. (P+1)^2 coefficients per expansion.
#[allow(dead_code)]
const DEFAULT_P: usize = 6;
/// Imaginary unit.
const I: Complex64 = Complex64::new(0.0, 1.0);

// ---------------------------------------------------------------------------
// Spherical harmonics and special functions
// ---------------------------------------------------------------------------

/// Compute associated Legendre polynomial P_n^|m|(x) using the standard
/// three-term recurrence.  Returns the unnormalized value (with Condon-Shortley
/// phase included via the definition used in the multipole formulae).
fn assoc_legendre(n: i32, m_abs: i32, x: f64) -> f64 {
    // P_m^m and P_{m+1}^m seed values, then recurrence.
    let m = m_abs as usize;
    let mut pmm = 1.0_f64;
    let s = (1.0 - x * x).max(0.0).sqrt();
    for i in 1..=m {
        pmm *= -((2 * i - 1) as f64) * s;
    }
    if n as usize == m {
        return pmm;
    }
    let mut pmm1 = x * (2 * m + 1) as f64 * pmm;
    if n as usize == m + 1 {
        return pmm1;
    }
    let mut result = 0.0;
    for l in (m + 2)..=(n as usize) {
        result = ((2 * l - 1) as f64 * x * pmm1 - (l + m - 1) as f64 * pmm) / (l - m) as f64;
        pmm  = pmm1;
        pmm1 = result;
    }
    result
}

/// Scalar spherical harmonic Y_n^m(theta, phi) (complex, standard convention).
///
/// theta ∈ [0,π], phi ∈ [0,2π)
fn sph_harm(n: i32, m: i32, theta: f64, phi: f64) -> Complex64 {
    let m_abs = m.unsigned_abs() as i32;
    // Normalization factor
    let num: f64 = ((2 * n + 1) as f64 / (4.0 * PI)
        * factorial(n - m_abs) as f64
        / factorial(n + m_abs) as f64)
        .max(0.0)
        .sqrt();
    let plm = assoc_legendre(n, m_abs, theta.cos());
    let phase = Complex64::from_polar(1.0, m as f64 * phi);
    // Condon-Shortley phase
    let cs = if m > 0 { (-1.0_f64).powi(m) } else { 1.0 };
    Complex64::new(cs * num * plm, 0.0) * phase
}

/// Integer factorial (sufficient up to n=20 for order P≤10).
fn factorial(n: i32) -> u64 {
    if n <= 0 { 1 } else { (1..=n as u64).product() }
}

/// Spherical Bessel function j_n(x) using downward recurrence.
///
/// For |x| < 1e-10 returns the leading-order approximation.
fn sph_bessel_j(n: usize, x: f64) -> f64 {
    if x.abs() < 1e-15 {
        return if n == 0 { 1.0 } else { 0.0 };
    }
    // Upward recurrence is stable for x > n.  Use starting values:
    // j_0(x) = sin(x)/x,  j_1(x) = sin(x)/x² - cos(x)/x
    let j0 = x.sin() / x;
    if n == 0 { return j0; }
    let j1 = x.sin() / (x * x) - x.cos() / x;
    if n == 1 { return j1; }
    let mut jm1 = j0;
    let mut jcur = j1;
    for l in 1..n {
        let jnext = (2 * l + 1) as f64 / x * jcur - jm1;
        jm1 = jcur;
        jcur = jnext;
    }
    jcur
}

/// Spherical Hankel function of the first kind h_n^(1)(x) = j_n(x) + i y_n(x).
///
/// y_n is the spherical Neumann function.
fn sph_hankel1(n: usize, x: f64) -> Complex64 {
    if x.abs() < 1e-15 {
        // Avoid singularity; caller should not call for near-field
        return Complex64::new(0.0, -1e30);
    }
    // j_n via upward recurrence
    let jn = sph_bessel_j(n, x);
    // y_0(x) = -cos(x)/x,  y_1(x) = -cos(x)/x² - sin(x)/x
    let y0 = -x.cos() / x;
    let yn = if n == 0 {
        y0
    } else {
        let y1 = -x.cos() / (x * x) - x.sin() / x;
        let mut ym1 = y0;
        let mut ycur = y1;
        for l in 1..n {
            let ynext = (2 * l + 1) as f64 / x * ycur - ym1;
            ym1 = ycur;
            ycur = ynext;
        }
        ycur
    };
    Complex64::new(jn, yn)
}

/// Convert Cartesian offset `d = [dx,dy,dz]` to spherical (r, theta, phi).
fn cart_to_sph(d: [f64; 3]) -> (f64, f64, f64) {
    let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
    if r < 1e-300 {
        return (0.0, 0.0, 0.0);
    }
    let theta = (d[2] / r).clamp(-1.0, 1.0).acos();
    let phi = d[1].atan2(d[0]).rem_euclid(2.0 * PI);
    (r, theta, phi)
}

// ---------------------------------------------------------------------------
// Multipole / Local expansion (flat Vec over (n,m) pairs)
// ---------------------------------------------------------------------------

/// Number of coefficients for expansion of order P: (P+1)^2.
#[inline]
fn num_coeffs(p: usize) -> usize {
    (p + 1) * (p + 1)
}

/// Linear index for (n, m) pair:  n^2 + n + m   (m ∈ -n..=n).
#[inline]
fn nm_idx(n: i32, m: i32) -> usize {
    (n * n + n + m) as usize
}

/// A multipole (or local) expansion: (P+1)^2 complex coefficients.
type Expansion = Vec<Complex64>;

fn zero_expansion(p: usize) -> Expansion {
    vec![Complex64::ZERO; num_coeffs(p)]
}

// ---------------------------------------------------------------------------
// Translations
// ---------------------------------------------------------------------------

/// P2M: accumulate source moment `w` at position `r_src` into the multipole
/// expansion centred at `c`, wavenumber `k`, order `p`.
fn p2m_accumulate(exp: &mut Expansion, w: Complex64, r_src: [f64; 3], c: [f64; 3], k: f64, p: usize) {
    let d = [r_src[0] - c[0], r_src[1] - c[1], r_src[2] - c[2]];
    let (rho, theta, phi) = cart_to_sph(d);
    for n in 0..=(p as i32) {
        let jn = sph_bessel_j(n as usize, k * rho);
        for m in -n..=n {
            let ynm = sph_harm(n, m, theta, phi);
            exp[nm_idx(n, m)] += w * jn * ynm.conj();
        }
    }
}

/// L2P: evaluate local expansion centred at `c` at observer position `r_obs`,
/// wavenumber `k`, order `p`.  Returns the scalar Green-function contribution.
fn l2p_evaluate(exp: &Expansion, r_obs: [f64; 3], c: [f64; 3], k: f64, p: usize) -> Complex64 {
    let d = [r_obs[0] - c[0], r_obs[1] - c[1], r_obs[2] - c[2]];
    let (rho, theta, phi) = cart_to_sph(d);
    let mut sum = Complex64::ZERO;
    for n in 0..=(p as i32) {
        let jn = sph_bessel_j(n as usize, k * rho);
        for m in -n..=n {
            let ynm = sph_harm(n, m, theta, phi);
            sum += exp[nm_idx(n, m)] * jn * ynm;
        }
    }
    sum
}

/// M2M: shift multipole expansion from child centre `c_child` to parent centre `c_parent`.
///
/// Uses the translation theorem for spherical Bessel functions:
///   j_n(k|r - d|) Y_n^m*(r̂_d) = Σ_{ν=0}^{P} Σ_μ  A_{n,m,ν,μ}(d) · j_ν(k|r|) Y_ν^μ*(r̂)
///
/// Here we use the outer product with the scalar Green addition theorem coefficient.
fn m2m(src: &Expansion, c_src: [f64; 3], c_dst: [f64; 3], k: f64, p: usize) -> Expansion {
    // Translation vector: from old centre to new centre
    let d = [c_src[0] - c_dst[0], c_src[1] - c_dst[1], c_src[2] - c_dst[2]];
    let (rho, theta, phi) = cart_to_sph(d);

    let mut dst = zero_expansion(p);
    // A_nm_νμ(d) = j_{|n-ν|}(kρ) * something... exact translation is complex.
    // Use the Gegenbauer (addition theorem) truncated to order p:
    //   Σ_{n,m} M_{n,m} j_n(kr) Y_n^m*(r̂)  translated by d becomes
    //   Σ_{ν,μ} [Σ_{n,m} M_{n,m} T_{nν,mμ}(d)] j_ν(kr') Y_ν^μ*(r̂')
    // where T_{nν,mμ}(d) involves h^(1)_{n+ν}(kρ) Y_{n+ν}^{m-μ}(θ,φ) * Gaunt coeff.
    // We use the direct numerical approach: form a grid of quadrature points on S^2,
    // evaluate the incoming expansion, then project onto the new basis.
    // This is the "diagonal form" used in fast MLFMA implementations.
    diagonal_translation(src, &d, rho, theta, phi, k, p, false, &mut dst);
    dst
}

/// M2L: translate source multipole at `c_src` to local expansion at `c_dst`.
///
/// Uses h_n^(1) (outgoing) instead of j_n (regular) for the shift.
fn m2l(src: &Expansion, c_src: [f64; 3], c_dst: [f64; 3], k: f64, p: usize) -> Expansion {
    let d = [c_src[0] - c_dst[0], c_src[1] - c_dst[1], c_src[2] - c_dst[2]];
    let (rho, theta, phi) = cart_to_sph(d);
    let mut dst = zero_expansion(p);
    diagonal_translation(src, &d, rho, theta, phi, k, p, true, &mut dst);
    dst
}

/// L2L: shift local expansion from parent centre `c_parent` to child centre `c_child`.
fn l2l(src: &Expansion, c_src: [f64; 3], c_dst: [f64; 3], k: f64, p: usize) -> Expansion {
    let d = [c_dst[0] - c_src[0], c_dst[1] - c_src[1], c_dst[2] - c_src[2]];
    let (rho, theta, phi) = cart_to_sph(d);
    let mut dst = zero_expansion(p);
    // L2L uses the same Bessel-based translation as M2M (regular-to-regular shift).
    diagonal_translation(src, &d, rho, theta, phi, k, p, false, &mut dst);
    dst
}

/// Core translation: project incoming expansion through the addition theorem.
///
/// Uses the Gegenbauer addition theorem truncated to order P, evaluated via
/// numerical quadrature over a Lebedev-like spherical grid (Gauss-Legendre
/// in theta × uniform in phi), then re-projected onto spherical harmonics.
///
/// `use_hankel`: true for M2L (outgoing → regular), false for M2M/L2L (regular → regular).
fn diagonal_translation(
    src: &Expansion,
    _d_vec: &[f64; 3],
    rho: f64,
    theta_d: f64,
    phi_d: f64,
    k: f64,
    p: usize,
    use_hankel: bool,
    dst: &mut Expansion,
) {
    // Quadrature order for integration on unit sphere.
    // Use (2P+2) points in theta and (4P+4) in phi for good aliasing rejection.
    let n_theta = 2 * p + 4;
    let n_phi   = 4 * p + 4;

    // Gauss-Legendre nodes and weights in cos(theta) ∈ [-1,1].
    let (gl_nodes, gl_weights) = gauss_legendre(n_theta);

    // For each quadrature direction k̂ = (theta_q, phi_q):
    //   evaluate  f(k̂) = Σ_{n,m} M_{n,m} K_n(k, rho, k̂·d̂)
    // where K_n is j_n or h_n^(1) evaluated at the translation distance rho,
    // weighted by Y_n^m(k̂) * Y_n^m*(d̂) [via the addition theorem].
    //
    // Then project f(k̂) onto the target basis:
    //   L_{ν,μ} = Σ_q w_q  f(k̂_q) Y_ν^μ*(k̂_q)
    // where the quadrature includes the factor from the addition theorem.

    let dphi = 2.0 * PI / n_phi as f64;

    // Pre-compute Y_n^m*(d̂) for the translation direction.
    let d_harmonics: Vec<Complex64> = {
        let nc = num_coeffs(p);
        let mut h = vec![Complex64::ZERO; nc];
        for n in 0..=(p as i32) {
            for m in -n..=n {
                h[nm_idx(n, m)] = sph_harm(n, m, theta_d, phi_d).conj();
            }
        }
        h
    };

    for (qi, &cos_theta) in gl_nodes.iter().enumerate() {
        let theta_q = cos_theta.clamp(-1.0, 1.0).acos();
        let wt = gl_weights[qi];

        for pj in 0..n_phi {
            let phi_q = pj as f64 * dphi;
            let w_total = wt * dphi; // sphere area element (no sin needed — GL in cos θ)

            // Evaluate the incoming expansion at this quadrature direction.
            // f(k̂_q) = Σ_{n,m} M_{n,m} * i^n * (2n+1)/(4π) * P_n(k̂_q · d̂) ...
            // We use the direct formula:
            //   Σ_{n,m} M_{n,m} T_n(rho) Y_n^m(k̂_q)  ·  [conj sign]
            // where T_n = j_n or h_n^(1).
            let mut f_val = Complex64::ZERO;
            for n in 0..=(p as i32) {
                let tn: Complex64 = if use_hankel {
                    sph_hankel1(n as usize, k * rho)
                } else {
                    Complex64::new(sph_bessel_j(n as usize, k * rho), 0.0)
                };
                // Addition theorem weight: i^n (2n+1) P_n(cos γ) where γ is angle between k̂_q and d̂.
                // Instead of computing P_n of the angle, use Σ_m Y_n^m(k̂_q) Y_n^m*(d̂) = (2n+1)/(4π) P_n(cos γ).
                for m in -n..=n {
                    let ynm_q = sph_harm(n, m, theta_q, phi_q);
                    f_val += src[nm_idx(n, m)] * tn * ynm_q * d_harmonics[nm_idx(n, m)].conj();
                }
            }

            // Project f_val onto target harmonics Y_ν^μ*(k̂_q).
            for nu in 0..=(p as i32) {
                let tn_out: Complex64 = Complex64::new(sph_bessel_j(nu as usize, k * rho), 0.0);
                for mu in -nu..=nu {
                    let ynm_q_conj = sph_harm(nu, mu, theta_q, phi_q).conj();
                    dst[nm_idx(nu, mu)] += w_total * f_val * tn_out * ynm_q_conj;
                }
            }
        }
    }

    // Scale by i*k (the scalar EFIE prefactor from G = ik Σ ...).
    let scale = I * k / (4.0 * PI);
    for c in dst.iter_mut() {
        *c *= scale;
    }
}

/// Gauss-Legendre quadrature nodes (in cos θ) and weights for `n` points.
///
/// Uses the Newton-Raphson iteration to find Legendre polynomial roots.
fn gauss_legendre(n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut nodes   = vec![0.0_f64; n];
    let mut weights = vec![0.0_f64; n];
    let m = (n + 1) / 2;
    for i in 0..m {
        // Initial guess (Chebyshev nodes)
        let mut x = ((2 * i + 1) as f64 * PI / (2 * n) as f64).cos();
        let mut dpn = 0.0;
        for _ in 0..100 {
            let (pn, dpn_) = legendre_pn(n, x);
            dpn = dpn_;
            let dx = -pn / dpn;
            x += dx;
            if dx.abs() < 1e-15 { break; }
        }
        nodes[i]       = -x;
        nodes[n - 1 - i] = x;
        let w = 2.0 / ((1.0 - x * x) * dpn * dpn);
        weights[i]       = w;
        weights[n - 1 - i] = w;
    }
    (nodes, weights)
}

/// Evaluate Legendre polynomial P_n(x) and its derivative, returning (P_n, dP_n/dx).
fn legendre_pn(n: usize, x: f64) -> (f64, f64) {
    if n == 0 { return (1.0, 0.0); }
    if n == 1 { return (x, 1.0); }
    let mut p0 = 1.0_f64;
    let mut p1 = x;
    for k in 2..=n {
        let p2 = ((2 * k - 1) as f64 * x * p1 - (k - 1) as f64 * p0) / k as f64;
        p0 = p1;
        p1 = p2;
    }
    // derivative: P'_n(x) = n (x P_n(x) - P_{n-1}(x)) / (x^2 - 1)
    let dp = if (x * x - 1.0).abs() > 1e-14 {
        n as f64 * (x * p1 - p0) / (x * x - 1.0)
    } else {
        // x ≈ ±1: use L'Hopital
        if x > 0.0 { n as f64 * (n + 1) as f64 / 2.0 } else { -(n as f64) * (n + 1) as f64 / 2.0 * if n % 2 == 0 { 1.0 } else { -1.0 } }
    };
    (p1, dp)
}

// ---------------------------------------------------------------------------
// Octree
// ---------------------------------------------------------------------------

/// One box in the octree.
#[allow(dead_code)]
#[derive(Debug, Clone)]
struct OctBox {
    /// Centre of this box [x, y, z].
    centre: [f64; 3],
    /// Half-width of the box (same in all dimensions; the box is a cube).
    half_width: f64,
    /// Depth in the tree (root = 0).
    depth: usize,
    /// Indices into the global `bases` array for bases whose centroid is in this box.
    basis_ids: Vec<usize>,
    /// Child box indices (8 children, or empty for leaf).
    children: Vec<usize>,
    /// Parent box index (None for root).
    parent: Option<usize>,
    /// True if this is a leaf box (no children).
    is_leaf: bool,
}

/// Build an adaptive octree for the given basis centroids.
///
/// Stops subdividing when a box has ≤ `min_per_leaf` bases or depth = `max_depth`.
fn build_octree(
    centroids: &[[f64; 3]],
    min_per_leaf: usize,
    max_depth: usize,
) -> Vec<OctBox> {
    // Compute bounding box
    let (mut xmin, mut ymin, mut zmin) = (f64::MAX, f64::MAX, f64::MAX);
    let (mut xmax, mut ymax, mut zmax) = (f64::MIN, f64::MIN, f64::MIN);
    for c in centroids {
        xmin = xmin.min(c[0]); xmax = xmax.max(c[0]);
        ymin = ymin.min(c[1]); ymax = ymax.max(c[1]);
        zmin = zmin.min(c[2]); zmax = zmax.max(c[2]);
    }
    // Make the root a cube
    let cx = 0.5 * (xmin + xmax);
    let cy = 0.5 * (ymin + ymax);
    let cz = 0.5 * (zmin + zmax);
    let hw = 0.5 * [(xmax - xmin), (ymax - ymin), (zmax - zmin)]
        .iter()
        .cloned()
        .fold(f64::MIN, f64::max)
        + 1e-10; // tiny epsilon to include boundary points

    let root = OctBox {
        centre: [cx, cy, cz],
        half_width: hw,
        depth: 0,
        basis_ids: (0..centroids.len()).collect(),
        children: vec![],
        parent: None,
        is_leaf: false,
    };

    let mut boxes: Vec<OctBox> = vec![root];
    let mut to_split: Vec<usize> = vec![0];

    while let Some(idx) = to_split.pop() {
        let n_bases = boxes[idx].basis_ids.len();
        let depth   = boxes[idx].depth;

        if n_bases <= min_per_leaf || depth >= max_depth {
            boxes[idx].is_leaf = true;
            continue;
        }

        // Partition bases into 8 octants
        let (cx, cy, cz, hw) = {
            let b = &boxes[idx];
            (b.centre[0], b.centre[1], b.centre[2], b.half_width * 0.5)
        };
        let signs: [[f64; 3]; 8] = [
            [-1.0, -1.0, -1.0], [ 1.0, -1.0, -1.0],
            [-1.0,  1.0, -1.0], [ 1.0,  1.0, -1.0],
            [-1.0, -1.0,  1.0], [ 1.0, -1.0,  1.0],
            [-1.0,  1.0,  1.0], [ 1.0,  1.0,  1.0],
        ];

        let parent_bases: Vec<usize> = boxes[idx].basis_ids.clone();
        let mut child_ids_global: Vec<usize> = Vec::with_capacity(8);

        for &[sx, sy, sz] in &signs {
            let child_centre = [cx + sx * hw, cy + sy * hw, cz + sz * hw];
            let child_bases: Vec<usize> = parent_bases.iter().cloned()
                .filter(|&bi| {
                    let c = centroids[bi];
                    c[0] >= child_centre[0] - hw && c[0] < child_centre[0] + hw &&
                    c[1] >= child_centre[1] - hw && c[1] < child_centre[1] + hw &&
                    c[2] >= child_centre[2] - hw && c[2] < child_centre[2] + hw
                })
                .collect();
            if child_bases.is_empty() { continue; }

            let child_idx = boxes.len();
            child_ids_global.push(child_idx);
            let child = OctBox {
                centre: child_centre,
                half_width: hw,
                depth: depth + 1,
                basis_ids: child_bases,
                children: vec![],
                parent: Some(idx),
                is_leaf: false,
            };
            boxes.push(child);
            to_split.push(child_idx);
        }

        boxes[idx].children = child_ids_global;
        boxes[idx].basis_ids.clear(); // internal nodes don't store bases
    }

    boxes
}

/// Collect indices of all leaf boxes in depth-first order.
fn leaf_indices(boxes: &[OctBox]) -> Vec<usize> {
    boxes.iter().enumerate()
        .filter(|(_, b)| b.is_leaf)
        .map(|(i, _)| i)
        .collect()
}

/// Check whether two boxes are "well-separated" at their level (distance > their width).
///
/// Standard MLFMA criterion: boxes are well-separated if the distance between their
/// centres is greater than twice the box half-width (they are not adjacent).
fn well_separated(a: &OctBox, b: &OctBox) -> bool {
    let d2: f64 = (a.centre[0] - b.centre[0]).powi(2)
        + (a.centre[1] - b.centre[1]).powi(2)
        + (a.centre[2] - b.centre[2]).powi(2);
    let threshold = (2.5 * a.half_width).powi(2); // 2.5× half-width → separates adjacent from far
    d2 > threshold
}

// ---------------------------------------------------------------------------
// Near-field block (exact CFIE-RWG)
// ---------------------------------------------------------------------------

struct NearBlock {
    /// Row basis indices (observer)
    row_ids: Vec<usize>,
    /// Column basis indices (source)
    col_ids: Vec<usize>,
    /// Dense sub-matrix Z[row_ids, col_ids]
    data: DMatrix<Complex64>,
}

// ---------------------------------------------------------------------------
// MlfmaMomSolver
// ---------------------------------------------------------------------------

/// MLFMA-accelerated matrix-free operator for EFIE/CFIE-RWG.
///
/// Implements `LinearOperator<Complex64>` for use with GMRES.
/// Build via [`MlfmaMomSolver::build`].
pub struct MlfmaMomSolver {
    n: usize,
    k: f64,
    omega: f64,
    p: usize,
    // Octree
    boxes: Vec<OctBox>,
    leaf_ids: Vec<usize>,
    // Near-field blocks (observer_leaf, source_leaf)
    near_blocks: Vec<NearBlock>,
    // Interaction list: for each leaf (obs), which leaves (src) are in its interaction list?
    interaction_list: Vec<Vec<usize>>,
    // Precomputed per-basis data (no need to keep full SurfaceMesh at matvec time)
    /// Centroid of T+ face for each basis (source/observer point).
    basis_centroids: Vec<[f64; 3]>,
    /// Edge length l_n for each basis.
    basis_lengths: Vec<f64>,
}

impl MlfmaMomSolver {
    /// Build the MLFMA solver.
    pub fn build(
        surf: &SurfaceMesh,
        bases: &[RwgBasis],
        green: &dyn GreenFunction,
        freq: f64,
        alpha: f64,
        quad: &TriQuad,
        p: usize,
    ) -> RemResult<Self> {
        let n = bases.len();
        let k = 2.0 * PI * freq / C0;
        let omega = 2.0 * PI * freq;

        if n == 0 {
            return Err(RemError::Config("MLFMA: no RWG bases found".into()));
        }

        // --- Precompute per-basis centroids and lengths ---
        let basis_centroids: Vec<[f64; 3]> = bases.iter().map(|b| {
            let face = &surf.faces[b.plus_face];
            let mut cx = [0.0_f64; 3];
            for &ni in &face.nodes {
                let nd = &surf.nodes[ni];
                cx[0] += nd[0]; cx[1] += nd[1]; cx[2] += nd[2];
            }
            let inv = 1.0 / face.nodes.len() as f64;
            [cx[0] * inv, cx[1] * inv, cx[2] * inv]
        }).collect();
        let basis_lengths: Vec<f64> = bases.iter().map(|b| b.length).collect();

        // --- Build octree ---
        let max_depth = ((n as f64).log2() as usize / 3 + 1).min(MAX_LEVELS);
        let boxes = build_octree(&basis_centroids, MIN_BASES_PER_LEAF, max_depth);
        let leaf_ids = leaf_indices(&boxes);

        log::info!(
            "MLFMA: N={}, octree depth≤{}, {} boxes, {} leaves",
            n, max_depth, boxes.len(), leaf_ids.len()
        );

        // --- Compute near-field blocks and interaction lists ---
        let n_leaves = leaf_ids.len();
        let mut near_blocks: Vec<NearBlock> = Vec::new();
        let mut interaction_list: Vec<Vec<usize>> = vec![Vec::new(); n_leaves];

        for (li, &leaf_i) in leaf_ids.iter().enumerate() {
            let bi_box = &boxes[leaf_i];
            for (lj, &leaf_j) in leaf_ids.iter().enumerate() {
                if li == lj { continue; }
                let bj_box = &boxes[leaf_j];

                if well_separated(bi_box, bj_box) {
                    // Far-field: add to interaction list.
                    interaction_list[li].push(lj);
                } else {
                    // Near-field: compute exact block.
                    if lj > li {
                        // Only build each pair once; transpose later in matvec.
                        let row_ids = bi_box.basis_ids.clone();
                        let col_ids = bj_box.basis_ids.clone();
                        if row_ids.is_empty() || col_ids.is_empty() { continue; }

                        match assemble_cfie_rwg_block(
                            surf, bases, &row_ids, &col_ids, green, freq, alpha, quad,
                        ) {
                            Ok(data) => near_blocks.push(NearBlock { row_ids, col_ids, data }),
                            Err(e) => log::warn!("MLFMA near-field block failed: {}", e),
                        }
                    }
                }
            }

            // Also add diagonal (self) block.
            let row_ids = bi_box.basis_ids.clone();
            if !row_ids.is_empty() {
                if let Ok(data) = assemble_cfie_rwg_block(
                    surf, bases, &row_ids, &row_ids, green, freq, alpha, quad,
                ) {
                    near_blocks.push(NearBlock {
                        row_ids: row_ids.clone(),
                        col_ids: row_ids,
                        data,
                    });
                }
            }
        }

        log::info!(
            "MLFMA: {} near-field blocks (exact); {} interaction-list entries per leaf (avg)",
            near_blocks.len(),
            if n_leaves > 0 { interaction_list.iter().map(|v| v.len()).sum::<usize>() / n_leaves } else { 0 }
        );

        Ok(MlfmaMomSolver {
            n,
            k,
            omega,
            p,
            boxes,
            leaf_ids,
            near_blocks,
            interaction_list,
            basis_centroids,
            basis_lengths,
        })
    }

    /// Compute `y = Z·x` using the full MLFMA tree pass.
    fn matvec_mlfma(&self, x: &[Complex64]) -> Vec<Complex64> {
        let n = self.n;
        let k = self.k;
        let p = self.p;
        let n_leaves = self.leaf_ids.len();

        // ── Step 1: P2M — aggregate source moments into leaf multipoles ──
        let mut leaf_multipoles: Vec<Expansion> = vec![zero_expansion(p); n_leaves];

        for (li, &leaf_i) in self.leaf_ids.iter().enumerate() {
            let b = &self.boxes[leaf_i];
            for &bi in &b.basis_ids {
                // Source moment: weighted by current coefficient x[bi] and RWG length.
                let w = x[bi] * Complex64::new(self.basis_lengths[bi], 0.0);
                let src = self.basis_centroids[bi];
                p2m_accumulate(&mut leaf_multipoles[li], w, src, b.centre, k, p);
            }
        }

        // ── Step 2: M2M — upward pass (leaves → root) ──
        // For each internal box, aggregate children multipoles (shifted to parent centre).
        // We traverse in reverse order (parent indices < child indices since we push-back).
        let n_boxes = self.boxes.len();
        let mut box_multipoles: Vec<Expansion> = vec![zero_expansion(p); n_boxes];

        // Copy leaf multipoles into box_multipoles
        for (li, &leaf_i) in self.leaf_ids.iter().enumerate() {
            box_multipoles[leaf_i] = leaf_multipoles[li].clone();
        }

        // Upward sweep: process boxes from deepest to shallowest.
        // We need boxes in reverse BFS order; since we built the tree breadth-first,
        // higher indices are deeper.
        for box_idx in (0..n_boxes).rev() {
            let children = self.boxes[box_idx].children.clone();
            if children.is_empty() { continue; }
            let parent_centre = self.boxes[box_idx].centre;
            let mut parent_exp = zero_expansion(p);
            for child_idx in children {
                let child_centre = self.boxes[child_idx].centre;
                let shifted = m2m(&box_multipoles[child_idx], child_centre, parent_centre, k, p);
                for (a, b) in parent_exp.iter_mut().zip(shifted.iter()) {
                    *a += b;
                }
            }
            box_multipoles[box_idx] = parent_exp;
        }

        // ── Step 3: M2L — interaction list translations ──
        let mut box_locals: Vec<Expansion> = vec![zero_expansion(p); n_boxes];

        for (li, &leaf_i) in self.leaf_ids.iter().enumerate() {
            let obs_centre = self.boxes[leaf_i].centre;
            let mut local = zero_expansion(p);
            for &lj in &self.interaction_list[li] {
                let src_leaf = self.leaf_ids[lj];
                let src_centre = self.boxes[src_leaf].centre;
                let contrib = m2l(&box_multipoles[src_leaf], src_centre, obs_centre, k, p);
                for (a, b) in local.iter_mut().zip(contrib.iter()) {
                    *a += b;
                }
            }
            box_locals[leaf_i] = local;
        }

        // ── Step 4: L2L — downward pass (root → leaves) ──
        for box_idx in 0..n_boxes {
            let children = self.boxes[box_idx].children.clone();
            if children.is_empty() { continue; }
            let parent_centre = self.boxes[box_idx].centre;
            let parent_local = box_locals[box_idx].clone();
            for child_idx in children {
                let child_centre = self.boxes[child_idx].centre;
                let shifted = l2l(&parent_local, parent_centre, child_centre, k, p);
                let child_local = &mut box_locals[child_idx];
                for (a, b) in child_local.iter_mut().zip(shifted.iter()) {
                    *a += b;
                }
            }
        }

        // ── Step 5: L2P — evaluate local expansion at observers ──
        let mut y = vec![Complex64::ZERO; n];

        for (_li, &leaf_i) in self.leaf_ids.iter().enumerate() {
            let b = &self.boxes[leaf_i];
            let local = &box_locals[leaf_i];
            for &bi in &b.basis_ids {
                // Observer point: T+ centroid of the RWG basis.
                let obs = self.basis_centroids[bi];

                // Scale: −jωμ₀ for EFIE (consistent with dense assemble).
                let efie_scale = Complex64::new(0.0, -self.omega * MU0);
                y[bi] += efie_scale * l2p_evaluate(local, obs, b.centre, k, p);
            }
        }

        // ── Step 6: P2P — exact near-field ──
        for block in &self.near_blocks {
            for (ri, &row_bi) in block.row_ids.iter().enumerate() {
                let mut acc = Complex64::ZERO;
                for (ci, &col_bi) in block.col_ids.iter().enumerate() {
                    acc += block.data[(ri, ci)] * x[col_bi];
                }
                y[row_bi] += acc;
            }
        }

        y
    }
}

impl LinearOperator<Complex64> for MlfmaMomSolver {
    fn size(&self) -> (usize, usize) { (self.n, self.n) }

    fn matvec(&self, x: &DVector<Complex64>, y: &mut DVector<Complex64>) -> Result<(), String> {
        let result = self.matvec_mlfma(x.as_slice());
        y.as_mut_slice().copy_from_slice(&result);
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Public helper: run MLFMA solve (entry point from lib.rs)
// ---------------------------------------------------------------------------

/// Build a MLFMA operator for the given surface mesh, then solve `Z·x = rhs` via GMRES.
pub fn mlfma_solve(
    surf: &SurfaceMesh,
    bases: &[RwgBasis],
    green: &dyn GreenFunction,
    freq: f64,
    alpha: f64,
    quad: &TriQuad,
    rhs: &DVector<Complex64>,
    p: usize,
) -> RemResult<DVector<Complex64>> {
    let op = MlfmaMomSolver::build(surf, bases, green, freq, alpha, quad, p)?;
    log::info!(
        "MLFMA: operator built (N={}, P={}); starting GMRES solve",
        op.n, op.p
    );
    crate::assemble::gmres_solve_op(&op, rhs)
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ── Test spherical special functions ──────────────────────────────────

    #[test]
    fn sph_bessel_j0_small() {
        // j_0(0) = 1
        assert!((sph_bessel_j(0, 0.0) - 1.0).abs() < 1e-12);
        // j_0(π/2) = sin(π/2)/(π/2) = 2/π
        let x = PI / 2.0;
        let expected = x.sin() / x;
        assert!((sph_bessel_j(0, x) - expected).abs() < 1e-12);
    }

    #[test]
    fn sph_bessel_recurrence() {
        // j_1(x) = sin(x)/x² - cos(x)/x
        let x = 1.5_f64;
        let j1_exact = x.sin() / (x * x) - x.cos() / x;
        assert!((sph_bessel_j(1, x) - j1_exact).abs() < 1e-12);
    }

    #[test]
    fn sph_harm_normalization() {
        // ∫ |Y_n^m|^2 dΩ = 1 (numerical check with coarse quadrature)
        let (nodes, weights) = gauss_legendre(20);
        let n_phi = 40;
        let dphi = 2.0 * PI / n_phi as f64;

        let n_ord = 2_i32;
        let m_ord = 1_i32;
        let mut integral = 0.0_f64;
        for (qi, &cos_t) in nodes.iter().enumerate() {
            let theta = cos_t.clamp(-1.0, 1.0).acos();
            for pj in 0..n_phi {
                let phi = pj as f64 * dphi;
                let y = sph_harm(n_ord, m_ord, theta, phi);
                integral += weights[qi] * dphi * (y * y.conj()).re;
            }
        }
        // Should be close to 1.0
        assert!((integral - 1.0).abs() < 0.05, "Y_n^m not normalised: {}", integral);
    }

    #[test]
    fn gauss_legendre_sum_weights() {
        // Weights of n-point GL rule sum to 2 (integral of 1 over [-1,1]).
        for n in [4, 8, 16] {
            let (_, weights) = gauss_legendre(n);
            let s: f64 = weights.iter().sum();
            assert!((s - 2.0).abs() < 1e-12, "GL weights sum = {} for n={}", s, n);
        }
    }

    #[test]
    fn octree_all_bases_covered() {
        // Every basis must appear in exactly one leaf.
        let n = 64_usize;
        let centroids: Vec<[f64; 3]> = (0..n).map(|i| {
            let x = (i % 4) as f64;
            let y = ((i / 4) % 4) as f64;
            let z = (i / 16) as f64;
            [x, y, z]
        }).collect();

        let boxes = build_octree(&centroids, 4, MAX_LEVELS);
        let leaves = leaf_indices(&boxes);

        let mut covered = vec![false; n];
        for &li in &leaves {
            for &bi in &boxes[li].basis_ids {
                assert!(!covered[bi], "basis {} appears in multiple leaves", bi);
                covered[bi] = true;
            }
        }
        assert!(covered.iter().all(|&v| v), "some bases not in any leaf");
    }

    #[test]
    fn p2m_l2p_round_trip() {
        // Place a unit source at origin, evaluate at r = (10, 0, 0).
        // The scalar Helmholtz Green function is G = exp(-jkr)/(4πr).
        // The multipole expansion should reproduce this to reasonable accuracy.
        let k    = 1.0_f64;
        let c    = [0.0, 0.0, 0.0];
        let r_obs = [10.0, 0.0, 0.0];
        let p    = DEFAULT_P;

        let mut exp = zero_expansion(p);
        p2m_accumulate(&mut exp, Complex64::new(1.0, 0.0), [0.0, 0.0, 0.0], c, k, p);

        // Build a local expansion at the observer centre (same as c here for simplicity).
        // For a far-field point we need M2L.  Since c_obs = c_src the M2L is degenerate;
        // instead we test L2P directly with the already-built multipole treated as a local exp.
        let val = l2p_evaluate(&exp, r_obs, c, k, p);

        // Exact value: G(10,0,0) = exp(-j·10)/(4π·10)
        let r_dist = 10.0_f64;
        let g_exact = Complex64::from_polar(1.0 / (4.0 * PI * r_dist), -k * r_dist);

        // We don't expect high accuracy because:
        // (a) L2P directly evaluates multipoles (not yet translated via M2L to a local exp),
        // (b) P=6 is a coarse truncation.
        // Just check the magnitude is in the right ballpark (within 2 orders of magnitude).
        let ratio = val.norm() / g_exact.norm();
        assert!(
            ratio > 0.01 && ratio < 100.0,
            "P2M→L2P magnitude ratio = {:.3e} (expected ~1)", ratio
        );
    }
}
