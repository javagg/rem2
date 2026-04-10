//! Integration tests for MoM correctness.
//!
//! ## EFIE-pulse flat plate test
//!
//! For a flat PEC plate in the z=0 plane, normal-incidence plane wave:
//! - Surface current J_x ≈ 2 H_inc_y = 2/η₀ (physical optics limit)
//! - EFIE-pulse should reproduce this in the limit of large plate / fine mesh
//!
//! This is a sanity check that:
//!   1. The impedance matrix assembles without NaN/Inf
//!   2. The solve converges to a physically reasonable current magnitude
//!   3. Singular integrals (shared-edge/vertex) are handled without blowup
//!
//! ## Mie convergence test (qualitative)
//!
//! EFIE-pulse is a scalar approximation that models only J_x.
//! For a PEC sphere it underestimates back-scatter (which requires J_z)
//! but the forward-scatter magnitude should be in the right order of magnitude.
//! We test that the result is within a factor of 10 of Mie, confirming the
//! RCS formula and solver are not grossly wrong.

use rem_mom::{
    surface_mesh::{SurfaceMesh, TriFace, SharedEdge, tri_geometry, patch_edge_lengths},
    quadrature::TriQuad,
    assemble::{assemble_efie_pulse, lu_solve},
    excitation::plane_wave_rhs,
    postprocess::rcs_pattern,
    mie::pec_sphere_rcs,
};
use rem_core::{C0, ETA0};
use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Mesh helpers
// ---------------------------------------------------------------------------

fn flat_plate_mesh(nx: usize, ny: usize, lx: f64, ly: f64) -> SurfaceMesh {
    // nx × ny grid of quads → 2*nx*ny triangles in z=0 plane
    let mut nodes: Vec<[f64;3]> = Vec::new();
    for j in 0..=ny {
        for i in 0..=nx {
            let x = (i as f64 / nx as f64) * lx - lx/2.0;
            let y = (j as f64 / ny as f64) * ly - ly/2.0;
            nodes.push([x, y, 0.0]);
        }
    }
    let idx = |i: usize, j: usize| j*(nx+1) + i;
    let mut faces: Vec<TriFace> = Vec::new();
    for j in 0..ny {
        for i in 0..nx {
            let a = idx(i,j); let b = idx(i+1,j);
            let c = idx(i,j+1); let d = idx(i+1,j+1);
            for &[n0,n1,n2] in &[[a,b,d],[a,d,c]] {
                let (centroid,normal,area) = tri_geometry(&nodes[n0],&nodes[n1],&nodes[n2]);
                faces.push(TriFace { nodes:[n0,n1,n2], centroid, normal, area });
            }
        }
    }
    use std::collections::HashMap;
    let mut edge_map: HashMap<(usize,usize),Vec<usize>> = HashMap::new();
    for (fi,f) in faces.iter().enumerate() {
        let [a,b,c] = f.nodes;
        for &(u,v) in &[(a,b),(b,c),(c,a)] {
            let key = if u<v {(u,v)} else {(v,u)};
            edge_map.entry(key).or_default().push(fi);
        }
    }
    let mut edges: Vec<SharedEdge> = Vec::new();
    let mut boundary_edges: Vec<[usize;2]> = Vec::new();
    for ((n0,n1),fl) in &edge_map {
        match fl.len() {
            1 => { boundary_edges.push([*n0,*n1]); }
            2 => { edges.push(SharedEdge { nodes:[*n0,*n1], plus_face:fl[0], minus_face:fl[1], length:0.0 }); }
            _ => {}
        }
    }
    patch_edge_lengths(&mut edges, &nodes);
    let n_faces = faces.len();
    SurfaceMesh { nodes, faces, edges, boundary_edges, face_attrs: vec![0; n_faces] }
}

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
            let key=if u<v {(u,v)} else {(v,u)};
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
    let n_faces = faces.len();
    SurfaceMesh { nodes, faces, edges, boundary_edges:be, face_attrs: vec![0; n_faces] }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

/// Flat PEC plate: check that EFIE-pulse gives physically reasonable current.
///
/// For normal incidence (+z → -z), x-polarised plane wave on a plate in z=0:
/// Physical Optics (PO) current: J_x = 2 * H_inc_y = 2/η₀ [A/m]
/// (at ka >> 1; at low ka the MoM result will differ from PO)
///
/// We just verify the result is finite, nonzero, and within 2 orders of 2/η₀.
#[test]
fn flat_plate_current_finite() {
    let lx = 1.0_f64; let ly = 1.0_f64;
    let freq = 3e8; // 300 MHz → λ = 1 m → ka ≈ 1 (plate ~ 1λ × 1λ)
    let k = 2.0*PI*freq/C0;
    let surf = flat_plate_mesh(4, 4, lx, ly); // 32 triangles

    let quad = TriQuad::new(3);
    let z = assemble_efie_pulse(&surf, freq, &quad, 1e-6).expect("assemble");
    let rhs = plane_wave_rhs(&surf, k, "Pulse");
    let currents = lu_solve(&z, &rhs).expect("solve");

    // Check all currents are finite
    for (i,&j) in currents.iter().enumerate() {
        assert!(j.norm().is_finite(), "current[{}] = {:?} is not finite", i, j);
    }

    // Average current magnitude
    let avg_j = currents.iter().map(|j| j.norm()).sum::<f64>() / currents.len() as f64;
    let po_j = 2.0 / ETA0; // Physical Optics reference [A/m]

    println!("Flat plate: avg |J| = {:.3e} A/m, PO ref = {:.3e} A/m", avg_j, po_j);
    assert!(avg_j > po_j * 0.01, "Current too small: {:.3e}", avg_j);
    assert!(avg_j < po_j * 100.0, "Current too large: {:.3e}", avg_j);
}

/// EFIE-pulse on PEC sphere at ka=1: verify forward RCS is in the right order
/// of magnitude vs Mie (within a factor of 5).
///
/// EFIE-pulse is a scalar approximation (only J_x). Forward scatter is dominated
/// by the forward lit face J_x, so it should be qualitatively correct.
#[test]
fn sphere_forward_rcs_order_of_magnitude() {
    let a = 1.0_f64;
    let freq = C0 / (2.0*PI*a); // ka=1
    let k = 2.0*PI*freq/C0;
    let surf = icosphere(a, 1); // 80 faces

    let quad = TriQuad::new(3);
    let z = assemble_efie_pulse(&surf, freq, &quad, 1e-6).expect("assemble");
    let rhs = plane_wave_rhs(&surf, k, "Pulse");
    let currents = lu_solve(&z, &rhs).expect("solve");

    let theta = vec![0.0_f64]; // forward scatter only
    let rcs_mom = &rcs_pattern(&currents, &surf, k, &theta, &[0.0])[0][0];
    let rcs_mie = pec_sphere_rcs(a, k, &theta, None)[0];

    println!("Sphere ka=1 forward: RCS_MoM={:.3e}, RCS_Mie={:.3e}", rcs_mom, rcs_mie);

    // Require same order of magnitude (within factor 5)
    let ratio = rcs_mom / rcs_mie;
    assert!(ratio > 0.2 && ratio < 5.0,
        "Forward RCS ratio MoM/Mie = {:.2} outside [0.2, 5.0]", ratio);
}
