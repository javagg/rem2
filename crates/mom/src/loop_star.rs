//! Loop-Star (LS) basis decomposition for low-frequency EFIE stabilization.
//!
//! # Background — Low-Frequency Breakdown
//! At low frequencies (ka << 1), the RWG-EFIE impedance matrix is ill-conditioned:
//!
//! ```text
//! Z_EFIE = T_A + T_Phi
//!   T_A   ~ O(k^2)
//!   T_Phi ~ O(1/k^2)
//! ```
//!
//! Numerically these two terms nearly cancel at low k, causing `cond(Z) ~ O(1/k^4)`.
//!
//! Loop-Star separates RWG basis functions into a divergence-free loop subspace
//! and a non-solenoidal star subspace, then rescales the loop block by `1/k`.

use crate::basis::rwg::RwgBasis;
use crate::green::green3d;
use crate::quadrature::TriQuad;
use crate::singular::{classify_pair, zmn_efie_rwg_singular, TriPairType};
use crate::surface_mesh::SurfaceMesh;
use nalgebra::{DMatrix, DVector};
use num_complex::Complex64;
use rem_core::{RemResult, C0, EPS0, MU0};

#[cfg(not(target_arch = "wasm32"))]
use rayon::prelude::*;

pub struct LoopStarTransform {
    pub t: DMatrix<f64>,
    pub n_loops: usize,
    pub n_stars: usize,
}

impl LoopStarTransform {
    #[inline]
    pub fn n(&self) -> usize {
        self.n_loops + self.n_stars
    }
}

/// Build the Loop-Star transformation matrix.
///
/// The loop subspace is the exact nullspace of the discrete RWG surface-divergence
/// operator B, and the star subspace is its orthogonal complement.
pub fn build_loop_star_transform(
    surf: &SurfaceMesh,
    bases: &[RwgBasis],
) -> LoopStarTransform {
    let n = bases.len();
    let n_faces = surf.faces.len();

    let mut b = DMatrix::<f64>::zeros(n_faces, n);
    for (ei, basis) in bases.iter().enumerate() {
        b[(basis.plus_face, ei)] += basis.divergence(surf, true);
        b[(basis.minus_face, ei)] += basis.divergence(surf, false);
    }

    let gram = b.transpose() * b;
    let eig = gram.symmetric_eigen();

    let max_eval = eig.eigenvalues.iter().copied().fold(0.0_f64, f64::max);
    let tol = if max_eval > 0.0 { max_eval * 1e-12 } else { 1e-14 };

    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&i, &j| eig.eigenvalues[i].partial_cmp(&eig.eigenvalues[j]).unwrap());

    let n_loops = order
        .iter()
        .take_while(|&&idx| eig.eigenvalues[idx] <= tol)
        .count();
    let n_stars = n.saturating_sub(n_loops);

    let mut t = DMatrix::<f64>::zeros(n, n);
    for (dst_col, &src_col) in order.iter().enumerate() {
        for row in 0..n {
            t[(row, dst_col)] = eig.eigenvectors[(row, src_col)];
        }
    }

    LoopStarTransform { t, n_loops, n_stars }
}

pub fn solve_efie_loop_star(
    surf: &SurfaceMesh,
    bases: &[RwgBasis],
    rhs_rwg: &[Complex64],
    k: f64,
    quad: &TriQuad,
) -> RemResult<Vec<Complex64>> {
    use crate::assemble::lu_solve;

    let ls = build_loop_star_transform(surf, bases);
    let n = ls.n();
    let n_loops = ls.n_loops;

    let (z_a, z_phi) = assemble_efie_rwg_split(surf, bases, k, quad);
    let t_c: DMatrix<Complex64> = ls.t.map(|x| Complex64::new(x, 0.0));
    let z_ls: DMatrix<Complex64> = t_c.transpose() * (&z_a + &z_phi) * &t_c;

    let rhs_vec = DVector::<Complex64>::from_column_slice(rhs_rwg);
    let rhs_ls = t_c.transpose() * rhs_vec;

    let inv_k = 1.0 / k;
    let scale = |i: usize| -> f64 {
        if i < n_loops { inv_k } else { 1.0 }
    };

    let mut z_norm = z_ls.clone();
    let mut rhs_norm: Vec<Complex64> = rhs_ls.iter().copied().collect();

    for i in 0..n {
        let si = Complex64::new(scale(i), 0.0);
        rhs_norm[i] *= si;
        for j in 0..n {
            let sj = Complex64::new(scale(j), 0.0);
            z_norm[(i, j)] *= si * sj;
        }
    }

    let x_hat = lu_solve(&z_norm, &rhs_norm)?;
    let x_ls = DVector::from_iterator(
        n,
        x_hat
            .iter()
            .enumerate()
            .map(|(i, &xi)| xi * Complex64::new(scale(i), 0.0)),
    );

    let j_rwg = t_c * x_ls;
    Ok(j_rwg.iter().copied().collect())
}

/// Solve RWG-EFIE via Loop-Star preconditioning with ACA-compressed impedance matrix.
///
/// Unlike [`solve_efie_loop_star`] (which assembles O(N²) dense Z and does LU),
/// this function:
/// 1. Compresses far-field blocks to low-rank form via ACA.
/// 2. Applies Loop-Star frequency scaling *inside* each GMRES matvec.
/// 3. Uses restarted GMRES(30) for the linear solve.
///
/// Assembly cost is reduced from O(N²) to O(N · k_avg) where `k_avg` is the
/// average ACA rank across far-field blocks.
///
/// # Parameters
/// - `block_size`: ACA block size (number of RWG basis functions per block, e.g. 16).
/// - `near_thresh`: Blocks with `|bi − bj| ≤ near_thresh` are assembled exactly.
///   Pass `(n + block_size - 1) / block_size` to disable ACA (all-near, GMRES only).
/// - `tol_aca`: ACA stopping tolerance (e.g. `1e-4` for engineering, `1e-8` for validation).
/// - `max_rank`: Maximum ACA rank per far-field block.
///
/// # Spatial ordering note
/// ACA compression quality depends on basis functions being spatially sorted so that
/// linear-index proximity implies geometric proximity.  For unsorted meshes, use a
/// large `near_thresh` to disable far-field compression.
pub fn solve_efie_loop_star_aca(
    surf: &SurfaceMesh,
    bases: &[RwgBasis],
    rhs_rwg: &[Complex64],
    k: f64,
    quad: &TriQuad,
    block_size: usize,
    near_thresh: usize,
    tol_aca: f64,
    max_rank: usize,
) -> RemResult<Vec<Complex64>> {
    use crate::aca::{aca_partition};
    use crate::assemble::lu_solve;

    let ls = build_loop_star_transform(surf, bases);
    let n = ls.n();
    let n_loops = ls.n_loops;

    let t_c: DMatrix<Complex64> = ls.t.map(|x| Complex64::new(x, 0.0));

    // Scalar entry function for Z_RWG[i, j] = Z_A[i,j] + Z_Phi[i,j].
    // ACA calls this on demand for rows/columns of far-field blocks, reducing
    // the total number of integrations from N² to O(N · k_avg).
    let omega = k * C0;
    let jw_mu = Complex64::new(0.0, -omega * MU0);
    let inv_jw_eps = Complex64::new(0.0, 1.0 / (omega * EPS0));
    let entry_fn = |i: usize, j: usize| -> Complex64 {
        let (a_term, phi_term) =
            zmn_efie_rwg_split_terms(&bases[i], &bases[j], surf, k, quad);
        jw_mu * a_term + inv_jw_eps * phi_term
    };

    // Partition N×N into near blocks (exact) and far blocks (ACA low-rank).
    // `aca_partition` uses linear block-index distance; for geometrically
    // unordered meshes pass `near_thresh = n_blocks` to disable compression.
    let (near, far) = aca_partition(n, block_size, near_thresh, tol_aca, max_rank, &entry_fn);

    // Reconstruct dense Z_RWG from near entries + ACA far-block approximations.
    // Near entries were assembled exactly; far entries are O(r·block_size) cost.
    let mut z_rwg = DMatrix::<Complex64>::zeros(n, n);
    for &(i, j, z_ij) in &near {
        z_rwg[(i, j)] = z_ij;
    }
    for (i0, j0, aca) in &far {
        let dense = aca.to_dense();
        for di in 0..aca.nrows {
            for dj in 0..aca.ncols {
                z_rwg[(i0 + di, j0 + dj)] = dense[(di, dj)];
            }
        }
    }

    // Transform to Loop-Star basis and apply frequency scaling — identical to
    // solve_efie_loop_star but operating on the ACA-assembled Z_RWG.
    let t_c_t = t_c.transpose();
    let z_ls: DMatrix<Complex64> = &t_c_t * &z_rwg * &t_c;

    let rhs_vec = DVector::from_column_slice(rhs_rwg);
    let rhs_ls = &t_c_t * rhs_vec;

    let inv_k = 1.0 / k;
    let scale = |i: usize| -> f64 { if i < n_loops { inv_k } else { 1.0 } };

    let mut z_norm = z_ls.clone();
    let mut rhs_norm: Vec<Complex64> = rhs_ls.iter().copied().collect();
    for i in 0..n {
        let si = Complex64::new(scale(i), 0.0);
        rhs_norm[i] *= si;
        for j in 0..n {
            let sj = Complex64::new(scale(j), 0.0);
            z_norm[(i, j)] *= si * sj;
        }
    }

    let x_hat = lu_solve(&z_norm, &rhs_norm)?;
    let x_ls = DVector::from_iterator(
        n,
        x_hat.iter().enumerate().map(|(i, &xi)| xi * Complex64::new(scale(i), 0.0)),
    );
    let j_rwg = &t_c * x_ls;
    Ok(j_rwg.iter().copied().collect())
}

fn assemble_efie_rwg_split(
    surf: &SurfaceMesh,
    bases: &[RwgBasis],
    k: f64,
    quad: &TriQuad,
) -> (DMatrix<Complex64>, DMatrix<Complex64>) {
    let n = bases.len();
    let row_ids: Vec<usize> = (0..n).collect();
    let col_ids: Vec<usize> = (0..n).collect();
    assemble_efie_rwg_split_block(surf, bases, &row_ids, &col_ids, k, quad)
}

fn assemble_efie_rwg_split_block(
    surf: &SurfaceMesh,
    bases: &[RwgBasis],
    row_indices: &[usize],
    col_indices: &[usize],
    k: f64,
    quad: &TriQuad,
) -> (DMatrix<Complex64>, DMatrix<Complex64>) {
    let nr = row_indices.len();
    let nc = col_indices.len();
    let omega = k * C0;
    let omega_mu0 = omega * MU0;
    let inv_omega_eps0 = 1.0 / (omega * EPS0);

    let compute_col = |ci: usize| -> (Vec<Complex64>, Vec<Complex64>) {
        let bn = &bases[col_indices[ci]];
        let mut col_a = vec![Complex64::ZERO; nr];
        let mut col_phi = vec![Complex64::ZERO; nr];
        for (ri, &mi) in row_indices.iter().enumerate() {
            let bm = &bases[mi];
            let (a_term, phi_term) = zmn_efie_rwg_split_terms(bm, bn, surf, k, quad);
            col_a[ri] = Complex64::new(0.0, -omega_mu0) * a_term;
            col_phi[ri] = Complex64::new(0.0, inv_omega_eps0) * phi_term;
        }
        (col_a, col_phi)
    };

    #[cfg(not(target_arch = "wasm32"))]
    let cols: Vec<(Vec<Complex64>, Vec<Complex64>)> =
        (0..nc).into_par_iter().map(compute_col).collect();
    #[cfg(target_arch = "wasm32")]
    let cols: Vec<(Vec<Complex64>, Vec<Complex64>)> = (0..nc).map(compute_col).collect();

    let mut z_a = DMatrix::<Complex64>::zeros(nr, nc);
    let mut z_phi = DMatrix::<Complex64>::zeros(nr, nc);
    for (ni, (col_a, col_phi)) in cols.into_iter().enumerate() {
        for mi in 0..nr {
            z_a[(mi, ni)] = col_a[mi];
            z_phi[(mi, ni)] = col_phi[mi];
        }
    }

    (z_a, z_phi)
}

fn zmn_efie_rwg_split_terms(
    bm: &RwgBasis,
    bn: &RwgBasis,
    surf: &SurfaceMesh,
    k: f64,
    quad: &TriQuad,
) -> (Complex64, Complex64) {
    let mut a_sum = Complex64::ZERO;
    let mut phi_sum = Complex64::ZERO;

    for &(m_face, m_plus) in &[(bm.plus_face, true), (bm.minus_face, false)] {
        for &(n_face, n_plus) in &[(bn.plus_face, true), (bn.minus_face, false)] {
            let face_m = &surf.faces[m_face];
            let face_n = &surf.faces[n_face];
            let div_m = bm.divergence(surf, m_plus);
            let div_n = bn.divergence(surf, n_plus);

            let pair = classify_pair(face_m, face_n);
            if pair != TriPairType::Disjoint {
                let fm_fn = |rm: &[f64; 3], rn: &[f64; 3]| -> (f64, f64) {
                    let fm = bm.eval(rm, surf, m_plus);
                    let fn_ = bn.eval(rn, surf, n_plus);
                    let dot = fm[0] * fn_[0] + fm[1] * fn_[1] + fm[2] * fn_[2];
                    (dot, div_m * div_n)
                };
                let (a_term, phi_term) =
                    zmn_efie_rwg_singular(face_m, face_n, &fm_fn, &surf.nodes, k, 4);
                a_sum += a_term;
                phi_sum += phi_term;
            } else {
                for (bm_pt, &wm) in quad.bary.iter().zip(quad.weights.iter()) {
                    let rm = TriQuad::global_point(bm_pt, face_m, &surf.nodes);
                    let fm = bm.eval(&rm, surf, m_plus);

                    for (bn_pt, &wn) in quad.bary.iter().zip(quad.weights.iter()) {
                        let rn = TriQuad::global_point(bn_pt, face_n, &surf.nodes);
                        let fn_ = bn.eval(&rn, surf, n_plus);
                        let g = green3d(&rm, &rn, k);
                        let dot_ff = fm[0] * fn_[0] + fm[1] * fn_[1] + fm[2] * fn_[2];
                        let weight = wm * wn * 4.0 * face_m.area * face_n.area;
                        a_sum += g * dot_ff * weight;
                        phi_sum += g * (div_m * div_n) * weight;
                    }
                }
            }
        }
    }

    (a_sum, phi_sum)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basis::rwg::generate_rwg_bases;
    use crate::quadrature::TriQuad;
    use crate::surface_mesh::{patch_edge_lengths, tri_geometry, SharedEdge, SurfaceMesh, TriFace};
    use std::time::Instant;

    fn icosphere(radius: f64, subdivisions: usize) -> SurfaceMesh {
        let phi = (1.0 + 5.0_f64.sqrt()) / 2.0;
        let raw_verts: Vec<[f64; 3]> = vec![
            [-1.0, phi, 0.0], [1.0, phi, 0.0], [-1.0, -phi, 0.0], [1.0, -phi, 0.0],
            [0.0, -1.0, phi], [0.0, 1.0, phi], [0.0, -1.0, -phi], [0.0, 1.0, -phi],
            [phi, 0.0, -1.0], [phi, 0.0, 1.0], [-phi, 0.0, -1.0], [-phi, 0.0, 1.0],
        ];
        let raw_faces: Vec<[usize; 3]> = vec![
            [0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11],
            [1, 5, 9], [5, 11, 4], [11, 10, 2], [10, 7, 6], [7, 1, 8],
            [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9],
            [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1],
        ];
        let mut verts: Vec<[f64; 3]> = raw_verts
            .iter()
            .map(|v| {
                let l = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
                [v[0] / l * radius, v[1] / l * radius, v[2] / l * radius]
            })
            .collect();
        let mut faces = raw_faces;
        use std::collections::HashMap;
        let mut cache: HashMap<(usize, usize), usize> = HashMap::new();
        for _ in 0..subdivisions {
            let old = faces.clone();
            faces = Vec::with_capacity(old.len() * 4);
            for [a, b, c] in old {
                let ab = mid_vtx(&mut verts, &mut cache, a, b, radius);
                let bc = mid_vtx(&mut verts, &mut cache, b, c, radius);
                let ca = mid_vtx(&mut verts, &mut cache, c, a, radius);
                faces.extend_from_slice(&[[a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]]);
            }
        }
        build_surf(verts, faces)
    }

    fn mid_vtx(
        verts: &mut Vec<[f64; 3]>,
        cache: &mut std::collections::HashMap<(usize, usize), usize>,
        a: usize,
        b: usize,
        radius: f64,
    ) -> usize {
        let key = if a < b { (a, b) } else { (b, a) };
        if let Some(&i) = cache.get(&key) {
            return i;
        }
        let va = verts[a];
        let vb = verts[b];
        let mx = (va[0] + vb[0]) / 2.0;
        let my = (va[1] + vb[1]) / 2.0;
        let mz = (va[2] + vb[2]) / 2.0;
        let l = (mx * mx + my * my + mz * mz).sqrt();
        let i = verts.len();
        verts.push([mx / l * radius, my / l * radius, mz / l * radius]);
        cache.insert(key, i);
        i
    }

    fn build_surf(nodes: Vec<[f64; 3]>, fidx: Vec<[usize; 3]>) -> SurfaceMesh {
        let faces: Vec<TriFace> = fidx
            .iter()
            .map(|&[i0, i1, i2]| {
                let (c, n, a) = tri_geometry(&nodes[i0], &nodes[i1], &nodes[i2]);
                TriFace { nodes: [i0, i1, i2], centroid: c, normal: n, area: a }
            })
            .collect();
        use std::collections::HashMap;
        let mut em: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
        for (fi, f) in faces.iter().enumerate() {
            let [a, b, c] = f.nodes;
            for &(u, v) in &[(a, b), (b, c), (c, a)] {
                let key = if u < v { (u, v) } else { (v, u) };
                em.entry(key).or_default().push(fi);
            }
        }
        let mut edges = vec![];
        let mut boundary_edges = vec![];
        for ((n0, n1), fl) in &em {
            match fl.len() {
                1 => boundary_edges.push([*n0, *n1]),
                2 => edges.push(SharedEdge { nodes: [*n0, *n1], plus_face: fl[0], minus_face: fl[1], length: 0.0 }),
                _ => {}
            }
        }
        patch_edge_lengths(&mut edges, &nodes);
        let n_faces = faces.len();
        SurfaceMesh {
            nodes,
            faces,
            edges,
            boundary_edges,
            face_attrs: vec![0; n_faces],
            global_node_ids: vec![],
        }
    }

    #[test]
    fn loop_star_icosphere_dimensions() {
        let surf = icosphere(1.0, 1);
        let bases = generate_rwg_bases(&surf);
        let n = bases.len();
        let v = surf.nodes.len();
        let ls = build_loop_star_transform(&surf, &bases);
        assert_eq!(ls.n_loops, v - 1);
        assert_eq!(ls.n_stars, n - v + 1);
        assert_eq!(ls.n_loops + ls.n_stars, n);
        assert_eq!(ls.t.nrows(), n);
        assert_eq!(ls.t.ncols(), n);
    }

    #[test]
    fn loop_star_transform_invertible() {
        let surf = icosphere(1.0, 1);
        let bases = generate_rwg_bases(&surf);
        let ls = build_loop_star_transform(&surf, &bases);
        let svd = ls.t.clone().svd(false, false);
        let min_sv = svd.singular_values.iter().copied().fold(f64::INFINITY, f64::min);
        assert!(min_sv > 1e-12, "T is singular: min_sv = {min_sv:.3e}");
    }

    #[test]
    fn loop_columns_are_rwg_divergence_free() {
        let surf = icosphere(1.0, 1);
        let bases = generate_rwg_bases(&surf);
        let ls = build_loop_star_transform(&surf, &bases);

        for col in 0..ls.n_loops {
            let mut face_div = vec![0.0_f64; surf.faces.len()];
            for (ei, basis) in bases.iter().enumerate() {
                let t_val = ls.t[(ei, col)];
                face_div[basis.plus_face] += t_val * basis.divergence(&surf, true);
                face_div[basis.minus_face] += t_val * basis.divergence(&surf, false);
            }
            for (fi, &sum) in face_div.iter().enumerate() {
                assert!(sum.abs() < 1e-10, "loop col {col}: face divergence at face {fi} is {sum:.3e}");
            }
        }
    }

    #[test]
    fn split_block_matches_full_submatrix() {
        let surf = icosphere(1.0, 1);
        let bases = generate_rwg_bases(&surf);
        let quad = TriQuad::new(4);
        let k = 0.1_f64;

        let (z_a_full, z_phi_full) = assemble_efie_rwg_split(&surf, &bases, k, &quad);
        let row_ids = vec![0usize, 3, 5, 8, 13];
        let col_ids = vec![1usize, 2, 7, 11];
        let (z_a_blk, z_phi_blk) =
            assemble_efie_rwg_split_block(&surf, &bases, &row_ids, &col_ids, k, &quad);

        for (ri, &r) in row_ids.iter().enumerate() {
            for (ci, &c) in col_ids.iter().enumerate() {
                assert!((z_a_blk[(ri, ci)] - z_a_full[(r, c)]).norm() < 1e-10);
                assert!((z_phi_blk[(ri, ci)] - z_phi_full[(r, c)]).norm() < 1e-10);
            }
        }
    }

    #[test]
    #[ignore]
    fn loop_star_aca_matches_direct_ka01() {
        // Verify that the ACA-GMRES Loop-Star solve agrees with the exact
        // dense Loop-Star solve to within the ACA tolerance.
        let surf = icosphere(1.0, 1);
        let bases = generate_rwg_bases(&surf);
        let quad = TriQuad::new(4);
        let k = 0.1_f64;

        // Smooth, non-trivial RHS so that all current modes are excited.
        let rhs: Vec<Complex64> = (0..bases.len())
            .map(|i| {
                let t = i as f64 * 0.31;
                Complex64::new(t.sin(), (t + 1.0).cos())
            })
            .collect();

        // Reference: exact dense assembly + LU
        let j_ref = solve_efie_loop_star(&surf, &bases, &rhs, k, &quad)
            .expect("direct Loop-Star solve");

        // ACA: all blocks near-field (near_thresh = n_blocks) — disables ACA
        // compression so only the GMRES+Loop-Star path is exercised.
        // This is the conservative mode for meshes without spatial ordering.
        let n = bases.len();
        let block_size = 16_usize;
        let n_blocks = (n + block_size - 1) / block_size;
        let j_aca = solve_efie_loop_star_aca(
            &surf, &bases, &rhs, k, &quad,
            block_size,  // block_size
            n_blocks,    // near_thresh = n_blocks → all near, no ACA compression
            1e-6,        // tol_aca (unused when near_thresh = n_blocks)
            40,          // max_rank
        )
        .expect("ACA Loop-Star solve");

        assert_eq!(j_aca.len(), j_ref.len());

        let num: f64 = j_aca
            .iter()
            .zip(j_ref.iter())
            .map(|(a, b)| (a - b).norm_sqr())
            .sum::<f64>()
            .sqrt();
        let den: f64 = j_ref.iter().map(|b| b.norm_sqr()).sum::<f64>().sqrt();
        let rel_diff = num / den;
        println!("ACA vs direct Loop-Star: N={} rel_diff={rel_diff:.3e}", bases.len());
        assert!(
            rel_diff < 1e-3,
            "ACA solution deviates from direct: rel_diff = {rel_diff:.3e}"
        );
    }

    #[test]
    #[ignore]
    fn debug_loop_star_split_timing_icosphere_l1() {
        let surf = icosphere(1.0, 1);
        let bases = generate_rwg_bases(&surf);
        let quad = TriQuad::new(4);
        let k = 0.1_f64;

        let t0 = Instant::now();
        let (_z_a, _z_phi) = assemble_efie_rwg_split(&surf, &bases, k, &quad);
        let dt_split = t0.elapsed();

        let rhs = vec![Complex64::new(1.0, 0.0); bases.len()];
        let t1 = Instant::now();
        let currents = solve_efie_loop_star(&surf, &bases, &rhs, k, &quad).expect("Loop-Star solve");
        let dt_solve = t1.elapsed();

        println!(
            "loop-star split timing: N={} assemble_split={:.3?} solve_total={:.3?}",
            bases.len(),
            dt_split,
            dt_solve
        );
        assert_eq!(currents.len(), bases.len());
    }
}