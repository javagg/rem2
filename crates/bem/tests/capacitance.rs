//! BEM integration test: sphere capacitance via Laplace P0 BEM.
//!
//! Analytical: C = 4πε₀a  for a sphere of radius a.
//! We set φ = 1 V on the sphere surface, solve for σ = ∂φ/∂n,
//! compute Q = ∫σ dS, then C = ε₀ Q.
//!
//! With a coarse icosphere and P0 BEM the expected error is ~5-20%.

use rem_bem::{
    assemble::assemble_laplace_p0,
    assemble::solve_neumann,
    postprocess::capacitance,
};
use rem_mom::{
    surface_mesh::{SurfaceMesh, TriFace, SharedEdge, tri_geometry, patch_edge_lengths},
    quadrature::TriQuad,
};
use rem_core::EPS0;
use std::f64::consts::PI;
use nalgebra::DVector;

fn icosphere(radius: f64, subdivisions: usize) -> SurfaceMesh {
    let phi = (1.0 + 5.0_f64.sqrt()) / 2.0;
    let raw_verts: Vec<[f64;3]> = vec![
        [-1.0,phi,0.0],[1.0,phi,0.0],[-1.0,-phi,0.0],[1.0,-phi,0.0],
        [0.0,-1.0,phi],[0.0,1.0,phi],[0.0,-1.0,-phi],[0.0,1.0,-phi],
        [phi,0.0,-1.0],[phi,0.0,1.0],[-phi,0.0,-1.0],[-phi,0.0,1.0],
    ];
    let raw_faces: Vec<[usize;3]> = vec![
        [0,11,5],[0,5,1],[0,1,7],[0,7,10],[0,10,11],
        [1,5,9],[5,11,4],[11,10,2],[10,7,6],[7,1,8],
        [3,9,4],[3,4,2],[3,2,6],[3,6,8],[3,8,9],
        [4,9,5],[2,4,11],[6,2,10],[8,6,7],[9,8,1],
    ];
    let mut verts: Vec<[f64;3]> = raw_verts.iter().map(|v| {
        let l = (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]).sqrt();
        [v[0]/l*radius,v[1]/l*radius,v[2]/l*radius]
    }).collect();
    let mut faces = raw_faces;
    use std::collections::HashMap;
    let mut cache: HashMap<(usize,usize),usize> = HashMap::new();
    for _ in 0..subdivisions {
        let old = faces.clone(); faces = Vec::with_capacity(old.len()*4);
        for [a,b,c] in old {
            let ab = mid(&mut verts,&mut cache,a,b,radius);
            let bc = mid(&mut verts,&mut cache,b,c,radius);
            let ca = mid(&mut verts,&mut cache,c,a,radius);
            faces.extend_from_slice(&[[a,ab,ca],[b,bc,ab],[c,ca,bc],[ab,bc,ca]]);
        }
    }
    build_mesh(verts, faces)
}

fn mid(v: &mut Vec<[f64;3]>, cache: &mut std::collections::HashMap<(usize,usize),usize>,
       a: usize, b: usize, r: f64) -> usize {
    let key = if a<b {(a,b)} else {(b,a)};
    if let Some(&i) = cache.get(&key) { return i; }
    let va=v[a]; let vb=v[b];
    let mx=(va[0]+vb[0])/2.0; let my=(va[1]+vb[1])/2.0; let mz=(va[2]+vb[2])/2.0;
    let l=(mx*mx+my*my+mz*mz).sqrt(); let i=v.len();
    v.push([mx/l*r,my/l*r,mz/l*r]); cache.insert(key,i); i
}

fn build_mesh(nodes: Vec<[f64;3]>, fidx: Vec<[usize;3]>) -> SurfaceMesh {
    let faces: Vec<TriFace> = fidx.iter().map(|&[i0,i1,i2]| {
        let (c,n,a) = tri_geometry(&nodes[i0],&nodes[i1],&nodes[i2]);
        TriFace { nodes:[i0,i1,i2], centroid:c, normal:n, area:a }
    }).collect();
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

/// Sphere capacitance test using first-kind BIE.
///
/// Represent φ as single-layer potential: φ(r) = ∫_S G(r,r') σ(r') dS'
/// On the boundary: V σ = φ_D   (first-kind equation)
/// Then Q = ∫σ dS, C = ε₀ Q / φ_D.
///
/// Analytic: C = 4πε₀a
#[test]
fn sphere_capacitance_vs_analytic() {
    let a = 1.0_f64;
    let v0 = 1.0_f64;
    let c_analytic = 4.0 * PI * EPS0 * a;

    let surf = icosphere(a, 1); // 80 faces
    let n = surf.faces.len();
    let quad = TriQuad::new(5);

    let (v_mat, _k_mat) = assemble_laplace_p0(&surf, &quad, 5).unwrap();

    // Solve V σ = φ_D (first-kind Dirichlet)
    let phi_d = vec![v0; n];
    let b = DVector::<f64>::from_iterator(n, phi_d.iter().copied());
    let lu = v_mat.clone().lu();
    let x = lu
        .solve(&b)
        .expect("Dirichlet solve failed: V matrix may be singular");
    let sigma: Vec<f64> = x.iter().copied().collect();

    let c_bem = capacitance(&sigma, &surf, v0);

    let rel_err = ((c_bem - c_analytic) / c_analytic).abs();
    println!("Sphere capacitance: BEM={:.4e} F, analytic={:.4e} F, rel_err={:.1}%",
        c_bem, c_analytic, rel_err*100.0);

    assert!(c_bem > 0.0, "Capacitance must be positive");
    // First-kind collocation on coarse icosphere (80 faces): expect < 40%
    assert!(rel_err < 0.40,
        "Sphere capacitance error {:.1}% exceeds 40%", rel_err*100.0);
}
