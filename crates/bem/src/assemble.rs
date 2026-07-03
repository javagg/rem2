//! BEM matrix assembly for the Laplace equation.
//!
//! Assembles V (single-layer) and K (double-layer) matrices for P0 basis:
//!
//! V[m,n] = ∫_Tm ∫_Tn G(r,r') dS' dS    G = 1/(4πR)
//! K[m,n] = ∫_Tm ∫_Tn ∂G/∂n' dS' dS  +  ½δ_{mn}  (jump term on diagonal)

use rem_surface::surface_mesh::SurfaceMesh;
use rem_surface::quadrature::TriQuad;
use rem_surface::singular::{classify_pair, TriPairType};
use crate::kernel::{laplace_G, laplace_dG_dn};
use nalgebra::{DMatrix, DVector};
use rem_core::{RemError, RemResult};
use rayon::prelude::*;

pub fn assemble_laplace_p0(
    surf: &SurfaceMesh, quad: &TriQuad, n_duffy: usize,
) -> RemResult<(DMatrix<f64>, DMatrix<f64>)> {
    let n = surf.faces.len();
    if n == 0 { return Err(RemError::Mesh("Empty surface mesh".to_string())); }

    let cols_v: Vec<Vec<f64>> = (0..n).into_par_iter().map(|ni| {
        let fn_ = &surf.faces[ni];
        let mut col = vec![0.0_f64; n];
        for mi in 0..n {
            let fm = &surf.faces[mi];
            col[mi] = match classify_pair(fm, fn_) {
                TriPairType::Identical => duffy_self_G(fn_, &surf.nodes, n_duffy),
                _ => vmn(fm, fn_, &surf.nodes, quad),
            };
        }
        col
    }).collect();

    let cols_k: Vec<Vec<f64>> = (0..n).into_par_iter().map(|ni| {
        let fn_ = &surf.faces[ni];
        let mut col = vec![0.0_f64; n];
        for mi in 0..n {
            col[mi] = if mi == ni { 0.5 } else {
                kmn(&surf.faces[mi], fn_, &surf.nodes, quad)
            };
        }
        col
    }).collect();

    let mut v_mat = DMatrix::<f64>::zeros(n, n);
    let mut k_mat = DMatrix::<f64>::zeros(n, n);
    for (ni, (cv, ck)) in cols_v.into_iter().zip(cols_k.into_iter()).enumerate() {
        for mi in 0..n { v_mat[(mi, ni)] = cv[mi]; k_mat[(mi, ni)] = ck[mi]; }
    }
    Ok((v_mat, k_mat))
}

/// Solve Neumann problem: K φ = V σ  (K includes ½I)
pub fn solve_neumann(v: &DMatrix<f64>, k: &DMatrix<f64>, sigma: &[f64]) -> RemResult<Vec<f64>> {
    let n = sigma.len();
    let rhs = DVector::from_iterator(n, (0..n).map(|i| {
        (0..n).map(|j| v[(i, j)] * sigma[j]).sum::<f64>()
    }));
    let lu = k.clone().lu();
    let x = lu.solve(&rhs)
        .ok_or_else(|| RemError::Other("solve_neumann failed: singular matrix".to_string()))?;
    Ok(x.iter().copied().collect())
}

// ── Integral helpers ────────────────────────────────────────────────

fn vmn(fm: &rem_surface::surface_mesh::TriFace,
       fn_: &rem_surface::surface_mesh::TriFace,
       nodes: &[[f64; 3]], quad: &TriQuad) -> f64 {
    let rm = &fm.centroid;
    let mut val = 0.0;
    for (bn, &wn) in quad.bary.iter().zip(quad.weights.iter()) {
        let rn = TriQuad::global_point(bn, fn_, nodes);
        val += laplace_G(rm, &rn) * (wn * 2.0 * fn_.area);
    }
    val
}

fn kmn(fm: &rem_surface::surface_mesh::TriFace,
       fn_: &rem_surface::surface_mesh::TriFace,
       nodes: &[[f64; 3]], quad: &TriQuad) -> f64 {
    let rm = &fm.centroid;
    let np = &fn_.normal;
    let mut val = 0.0;
    for (bn, &wn) in quad.bary.iter().zip(quad.weights.iter()) {
        let rn = TriQuad::global_point(bn, fn_, nodes);
        val += laplace_dG_dn(rm, &rn, np) * (wn * 2.0 * fn_.area);
    }
    val
}

/// Self-term V[m,m] via Duffy transform for 1/(4πR) singularity.
fn duffy_self_G(face: &rem_surface::surface_mesh::TriFace,
                nodes: &[[f64; 3]], n_gauss: usize) -> f64 {
    use rem_surface::singular::gauss_legendre_1d;
    let (gl, gw) = gauss_legendre_1d(n_gauss);
    let [i0, i1, i2] = face.nodes;
    let v = [nodes[i0], nodes[i1], nodes[i2]];
    let mut sum = 0.0;
    for &pivot in &[0usize, 1, 2] {
        let (va, vb, vc) = (&v[pivot], &v[(pivot+1)%3], &v[(pivot+2)%3]);
        let area_sub = sub_tri_area(va, vb, vc);
        for (&rho, &w_rho) in gl.iter().zip(gw.iter()) {
            for (&theta, &w_theta) in gl.iter().zip(gw.iter()) {
                let rp = interp3(va, vb, vc, rho, theta);
                sum += laplace_G(&face.centroid, &rp)
                    * (w_rho * w_theta * 4.0 * area_sub * rho);
            }
        }
    }
    sum
}

fn interp3(va: &[f64;3], vb: &[f64;3], vc: &[f64;3], rho: f64, theta: f64) -> [f64;3] {
    [va[0]+rho*(theta*(vb[0]-va[0])+(1.0-theta)*(vc[0]-va[0])),
     va[1]+rho*(theta*(vb[1]-va[1])+(1.0-theta)*(vc[1]-va[1])),
     va[2]+rho*(theta*(vb[2]-va[2])+(1.0-theta)*(vc[2]-va[2]))]
}

fn sub_tri_area(va: &[f64;3], vb: &[f64;3], vc: &[f64;3]) -> f64 {
    let (e1, e2) = ([vb[0]-va[0],vb[1]-va[1],vb[2]-va[2]],
                   [vc[0]-va[0],vc[1]-va[1],vc[2]-va[2]]);
    let cx = e1[1]*e2[2]-e1[2]*e2[1];
    let cy = e1[2]*e2[0]-e1[0]*e2[2];
    let cz = e1[0]*e2[1]-e1[1]*e2[0];
    0.5 * (cx*cx+cy*cy+cz*cz).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use rem_surface::surface_mesh::{SurfaceMesh, TriFace, SharedEdge, tri_geometry, patch_edge_lengths};

    fn two_tri_mesh() -> SurfaceMesh {
        let nodes: Vec<[f64;3]> = vec![[0.0,0.0,0.0],[1.0,0.0,0.0],[0.5,1.0,0.0],[1.5,1.0,0.0]];
        let mut faces = Vec::new();
        for &[i0,i1,i2] in &[[0usize,1,2],[1,3,2]] {
            let (c,n,a) = tri_geometry(&nodes[i0],&nodes[i1],&nodes[i2]);
            faces.push(TriFace{nodes:[i0,i1,i2],centroid:c,normal:n,area:a});
        }
        let mut em = std::collections::HashMap::new();
        for (fi,f) in faces.iter().enumerate() {
            let [a,b,c]=f.nodes;
            for &(u,v) in &[(a,b),(b,c),(c,a)] { let key=if u<v{(u,v)}else{(v,u)}; em.entry(key).or_default().push(fi); }
        }
        let mut edges=Vec::new(); let mut be=Vec::new();
        for ((n0,n1),fl) in &em { match fl.len() { 1=>{be.push([*n0,*n1]);} 2=>{edges.push(SharedEdge{nodes:[*n0,*n1],plus_face:fl[0],minus_face:fl[1],length:0.0});} _=>{} } }
        patch_edge_lengths(&mut edges, &nodes);
        SurfaceMesh{nodes,faces,edges,boundary_edges:be,face_attrs:vec![0;2],global_node_ids:vec![]}
    }

    #[test] fn v_positive_diag() { let s=two_tri_mesh(); let (v,_)=assemble_laplace_p0(&s,&TriQuad::new(3),4).unwrap();
        for i in 0..v.nrows() { assert!(v[(i,i)]>0.0); } }
    #[test] fn k_diag_half() { let s=two_tri_mesh(); let (_,k)=assemble_laplace_p0(&s,&TriQuad::new(3),4).unwrap();
        for i in 0..k.nrows() { assert!((k[(i,i)]-0.5).abs()<1e-14); } }
    #[test] fn neumann_solve() { let s=two_tri_mesh(); let (v,k)=assemble_laplace_p0(&s,&TriQuad::new(3),4).unwrap();
        let phi=solve_neumann(&v,&k,&vec![1.0;2]).unwrap();
        for &p in &phi { assert!(p.is_finite()); } }
}
