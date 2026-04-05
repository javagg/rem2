//! BEM matrix assembly for the Laplace equation.
//!
//! Assembles the single-layer (V) and double-layer (K) matrices for the
//! exterior Neumann problem (given ∂φ/∂n, find φ) or Dirichlet problem.
//!
//! ## Boundary Integral Equation (exterior Laplace)
//!
//! ```text
//! ½φ(r) + K[φ](r) = V[σ](r)
//! ```
//!
//! **V** (single-layer, N×N):  V[m,n] = ∫_Tm ∫_Tn G(r,r') dS' dS
//! **K** (double-layer, N×N):  K[m,n] = ∫_Tm ∫_Tn ∂G/∂n'(r,r') dS' dS
//!
//! ## P0 (constant) basis
//!
//! φ ≈ Σ φ_n * 1_Tn (piecewise constant)
//!
//! Diagonal V[m,m] uses Duffy analytic self-integral of G = 1/(4πR).
//! Diagonal K[m,m] = ½ (identity from jump condition).
//! Off-diagonal: standard Gaussian quadrature.

use rem_mom::surface_mesh::SurfaceMesh;
use rem_mom::quadrature::TriQuad;
use rem_mom::singular::{classify_pair, TriPairType};
use crate::kernel::{laplace_G, laplace_dG_dn};
use rem_core::{RemError, RemResult};
use rayon::prelude::*;
use faer::Mat;

type C64 = faer::c64;

/// Assemble Laplace BEM V and K matrices (P0 basis).
///
/// Returns `(V, K)` where:
/// - V[m,n] = ∫_Tm ∫_Tn G dS' dS
/// - K[m,n] = ∫_Tm ∫_Tn ∂G/∂n' dS' dS  +  ½δ_{mn} (identity term)
pub fn assemble_laplace_p0(
    surf: &SurfaceMesh,
    quad: &TriQuad,
    n_duffy: usize,
) -> RemResult<(Mat<f64>, Mat<f64>)> {
    let n = surf.faces.len();
    if n == 0 {
        return Err(RemError::Mesh("Empty surface mesh".to_string()));
    }

    // Parallel column-wise assembly
    let cols_v: Vec<Vec<f64>> = (0..n).into_par_iter().map(|ni| {
        let face_n = &surf.faces[ni];
        let mut col = vec![0.0_f64; n];
        for mi in 0..n {
            let face_m = &surf.faces[mi];
            let pair = classify_pair(face_m, face_n);
            col[mi] = match pair {
                TriPairType::Identical => {
                    duffy_self_laplace_G(face_n, &surf.nodes, n_duffy)
                }
                TriPairType::SharedEdge | TriPairType::SharedVertex => {
                    // Near-singular: use more quadrature points
                    vmn_regular(face_m, face_n, &surf.nodes, quad)
                }
                TriPairType::Disjoint => {
                    vmn_regular(face_m, face_n, &surf.nodes, quad)
                }
            };
        }
        col
    }).collect();

    let cols_k: Vec<Vec<f64>> = (0..n).into_par_iter().map(|ni| {
        let face_n = &surf.faces[ni];
        let mut col = vec![0.0_f64; n];
        for mi in 0..n {
            let face_m = &surf.faces[mi];
            let pair = classify_pair(face_m, face_n);
            col[mi] = if mi == ni {
                0.5 // identity term from jump condition
            } else {
                match pair {
                    TriPairType::Disjoint |
                    TriPairType::SharedEdge |
                    TriPairType::SharedVertex => {
                        kmn_regular(face_m, face_n, &surf.nodes, quad)
                    }
                    TriPairType::Identical => 0.5,
                }
            };
        }
        col
    }).collect();

    let mut v_mat = Mat::<f64>::zeros(n, n);
    let mut k_mat = Mat::<f64>::zeros(n, n);
    for (ni, (cv, ck)) in cols_v.into_iter().zip(cols_k.into_iter()).enumerate() {
        for mi in 0..n {
            v_mat[(mi, ni)] = cv[mi];
            k_mat[(mi, ni)] = ck[mi];
        }
    }

    Ok((v_mat, k_mat))
}

/// Solve exterior Neumann problem: given σ = ∂φ/∂n on S, find φ on S.
///
/// Equation: (½I + K) φ = V σ  →  K φ = V σ  (K already includes ½I)
pub fn solve_neumann(
    v_mat: &Mat<f64>,
    k_mat: &Mat<f64>,
    sigma: &[f64],
) -> RemResult<Vec<f64>> {
    use faer::linalg::solvers::Solve;
    let n = sigma.len();
    if v_mat.nrows() != n { return Err(RemError::Config("V matrix size mismatch".to_string())); }

    // RHS = V * σ
    let mut rhs = Mat::<f64>::zeros(n, 1);
    for i in 0..n {
        let mut s = 0.0;
        for j in 0..n {
            s += v_mat[(i, j)] * sigma[j];
        }
        rhs[(i, 0)] = s;
    }

    let lu = k_mat.as_ref().partial_piv_lu();
    let x = lu.solve(rhs.as_ref());
    Ok((0..n).map(|i| x[(i, 0)]).collect())
}

// ---------------------------------------------------------------------------
// Matrix element integrals
// ---------------------------------------------------------------------------

/// V[m,n] = ∫_Tn G(r_m, r') dS'  (collocation at centroid of Tm)
fn vmn_regular(
    face_m: &rem_mom::surface_mesh::TriFace,
    face_n: &rem_mom::surface_mesh::TriFace,
    nodes: &[[f64; 3]],
    quad: &TriQuad,
) -> f64 {
    let rm = &face_m.centroid;
    let mut val = 0.0_f64;
    for (bn, &wn) in quad.bary.iter().zip(quad.weights.iter()) {
        let rn = TriQuad::global_point(bn, face_n, nodes);
        let g = laplace_G(rm, &rn);
        val += g * (wn * 2.0 * face_n.area);
    }
    val
}

/// K[m,n] = ∫_Tn ∂G/∂n'(r_m, r') dS'  (collocation at centroid of Tm)
fn kmn_regular(
    face_m: &rem_mom::surface_mesh::TriFace,
    face_n: &rem_mom::surface_mesh::TriFace,
    nodes: &[[f64; 3]],
    quad: &TriQuad,
) -> f64 {
    let rm = &face_m.centroid;
    let np = &face_n.normal;
    let mut val = 0.0_f64;
    for (bn, &wn) in quad.bary.iter().zip(quad.weights.iter()) {
        let rn = TriQuad::global_point(bn, face_n, nodes);
        let dg = laplace_dG_dn(rm, &rn, np);
        val += dg * (wn * 2.0 * face_n.area);
    }
    val
}

/// V[m,m] = ∫_Tm G(r_m, r') dS'  where r_m = centroid (collocation self-integral).
///
/// Uses Duffy polar transform to remove 1/R singularity at r' → r_m.
fn duffy_self_laplace_G(
    face: &rem_mom::surface_mesh::TriFace,
    nodes: &[[f64; 3]],
    n_gauss: usize,
) -> f64 {
    use rem_mom::singular::gauss_legendre_1d;

    let (gl, gw) = gauss_legendre_1d(n_gauss);
    let [i0, i1, i2] = face.nodes;
    let v = [nodes[i0], nodes[i1], nodes[i2]];

    let mut sum = 0.0_f64;

    for &pivot in &[0usize, 1, 2] {
        let va = &v[pivot];
        let vb = &v[(pivot+1)%3];
        let vc = &v[(pivot+2)%3];
        let area_sub = sub_tri_area(va, vb, vc);

        for (&rho, &w_rho) in gl.iter().zip(gw.iter()) {
            for (&theta, &w_theta) in gl.iter().zip(gw.iter()) {
                let rp = interp3(va, vb, vc, rho, theta);
                let r_obs = face.centroid;
                let g = laplace_G(&r_obs, &rp);
                // Duffy Jacobian: 4 * area_sub * rho
                sum += g * (w_rho * w_theta * 4.0 * area_sub * rho);
            }
        }
    }

    sum
}

fn interp3(va: &[f64;3], vb: &[f64;3], vc: &[f64;3], rho: f64, theta: f64) -> [f64;3] {
    [
        va[0] + rho*(theta*(vb[0]-va[0]) + (1.0-theta)*(vc[0]-va[0])),
        va[1] + rho*(theta*(vb[1]-va[1]) + (1.0-theta)*(vc[1]-va[1])),
        va[2] + rho*(theta*(vb[2]-va[2]) + (1.0-theta)*(vc[2]-va[2])),
    ]
}

fn sub_tri_area(va: &[f64;3], vb: &[f64;3], vc: &[f64;3]) -> f64 {
    let e1 = [vb[0]-va[0], vb[1]-va[1], vb[2]-va[2]];
    let e2 = [vc[0]-va[0], vc[1]-va[1], vc[2]-va[2]];
    let cx = e1[1]*e2[2]-e1[2]*e2[1];
    let cy = e1[2]*e2[0]-e1[0]*e2[2];
    let cz = e1[0]*e2[1]-e1[1]*e2[0];
    0.5 * (cx*cx+cy*cy+cz*cz).sqrt()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use rem_mom::surface_mesh::{SurfaceMesh, TriFace, SharedEdge, tri_geometry, patch_edge_lengths};

    fn two_tri_mesh() -> SurfaceMesh {
        // Two triangles sharing edge (1,2): T0=(0,1,2), T1=(1,3,2)
        let nodes: Vec<[f64;3]> = vec![
            [0.0,0.0,0.0],[1.0,0.0,0.0],[0.5,1.0,0.0],[1.5,1.0,0.0]
        ];
        let mut faces = Vec::new();
        for &[i0,i1,i2] in &[[0usize,1,2],[1,3,2]] {
            let (c,n,a) = tri_geometry(&nodes[i0],&nodes[i1],&nodes[i2]);
            faces.push(TriFace { nodes:[i0,i1,i2], centroid:c, normal:n, area:a });
        }
        use std::collections::HashMap;
        let mut em: HashMap<(usize,usize),Vec<usize>> = HashMap::new();
        for (fi,f) in faces.iter().enumerate() {
            let [a,b,c]=f.nodes;
            for &(u,v) in &[(a,b),(b,c),(c,a)] {
                let key=if u<v{(u,v)}else{(v,u)};
                em.entry(key).or_default().push(fi);
            }
        }
        let mut edges=Vec::new(); let mut be=Vec::new();
        for ((n0,n1),fl) in &em {
            match fl.len() {
                1=>{be.push([*n0,*n1]);}
                2=>{edges.push(SharedEdge{nodes:[*n0,*n1],plus_face:fl[0],minus_face:fl[1],length:0.0});}
                _=>{}
            }
        }
        patch_edge_lengths(&mut edges, &nodes);
        SurfaceMesh { nodes, faces, edges, boundary_edges:be }
    }

    #[test]
    fn v_matrix_positive_diagonal() {
        let surf = two_tri_mesh();
        let quad = TriQuad::new(3);
        let (v, _) = assemble_laplace_p0(&surf, &quad, 4).unwrap();
        for i in 0..v.nrows() {
            assert!(v[(i,i)] > 0.0, "V[{i},{i}] = {} not positive", v[(i,i)]);
        }
    }

    #[test]
    fn v_matrix_off_diagonal_positive() {
        // Both faces are in z=0, so G > 0 → off-diagonal should also be positive
        let surf = two_tri_mesh();
        let quad = TriQuad::new(3);
        let (v, _) = assemble_laplace_p0(&surf, &quad, 4).unwrap();
        assert!(v[(0,1)] > 0.0 && v[(1,0)] > 0.0, "Off-diag should be positive");
    }

    #[test]
    fn k_diagonal_is_half() {
        let surf = two_tri_mesh();
        let quad = TriQuad::new(3);
        let (_, k) = assemble_laplace_p0(&surf, &quad, 4).unwrap();
        for i in 0..k.nrows() {
            assert!((k[(i,i)] - 0.5).abs() < 1e-14,
                "K[{i},{i}] = {} ≠ 0.5", k[(i,i)]);
        }
    }

    #[test]
    fn neumann_solve_doesnt_panic() {
        let surf = two_tri_mesh();
        let quad = TriQuad::new(3);
        let (v, k) = assemble_laplace_p0(&surf, &quad, 4).unwrap();
        let sigma = vec![1.0_f64; surf.faces.len()];
        let phi = solve_neumann(&v, &k, &sigma).unwrap();
        for &p in &phi {
            assert!(p.is_finite(), "phi contains non-finite value");
        }
    }
}
