/// P1/P2 finite element stiffness matrix assembly.
///
/// Supports:
///  - 2-D: Tri3 (P1 linear triangle), Tri6 (P2 quadratic triangle)
///  - 3-D: Tet4 (P1 linear tetrahedron), Tet10 (P2 quadratic tetrahedron), Hex8 (trilinear)
///
/// The assembled system is:
///   K_ij = Σ_e  ε_e · ∫_Ωe  ∇φ_i · ∇φ_j  dΩ
///
/// No external dependencies (WASM-compatible).

use rem_core::{TripletMatrix, RemError, RemResult};
use rem_mesh::{RemMesh, ElementKind};
use rem_materials::DomainMap;

/// Assemble the global stiffness (diffusion) matrix for Poisson's equation.
///
/// `coeff_fn` maps physical group tag → diffusion coefficient (e.g. ε_abs).
pub fn assemble_stiffness(
    mesh: &RemMesh,
    coeff_fn: impl Fn(u32) -> f64,
) -> RemResult<TripletMatrix> {
    let n = mesh.n_nodes();
    // Upper bound on nnz: each Tri3 contributes 9 entries, Tet4 contributes 16.
    let cap = mesh.n_volume_elements() * 16;
    let mut triplet = TripletMatrix::with_capacity(n, n, cap);

    for elem in &mesh.volume_elements {
        if mesh.size > 1 && elem.rank != mesh.rank {
            continue;
        }
        let eps = coeff_fn(elem.tag);
        match elem.kind {
            ElementKind::Tri3 => {
                assemble_tri3(mesh, elem, eps, &mut triplet)?;
            }
            ElementKind::Tet4 => {
                assemble_tet4(mesh, elem, eps, &mut triplet)?;
            }
            ElementKind::Tri6 => {
                assemble_tri6(mesh, elem, eps, &mut triplet)?;
            }
            ElementKind::Tet10 => {
                assemble_tet10(mesh, elem, eps, &mut triplet)?;
            }
            ElementKind::Hex8 => {
                assemble_hex8(mesh, elem, eps, &mut triplet)?;
            }
            other => {
                log::warn!(
                    "Element kind {:?} not supported in P1 assembly — skipping",
                    other
                );
            }
        }
    }

    Ok(triplet)
}

/// Local stiffness for a linear triangle (Tri3).
///
/// Basis function gradients are constant; area integration is exact.
fn assemble_tri3(
    mesh: &RemMesh,
    elem: &rem_mesh::Element,
    eps: f64,
    triplet: &mut TripletMatrix,
) -> RemResult<()> {
    debug_assert_eq!(elem.node_ids.len(), 3);
    let [n0, n1, n2] = [elem.node_ids[0], elem.node_ids[1], elem.node_ids[2]];
    let (x0, y0) = (mesh.nodes[n0].x, mesh.nodes[n0].y);
    let (x1, y1) = (mesh.nodes[n1].x, mesh.nodes[n1].y);
    let (x2, y2) = (mesh.nodes[n2].x, mesh.nodes[n2].y);

    // det(J) = (x1-x0)(y2-y0) - (x2-x0)(y1-y0)
    let det_j = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
    let area = 0.5 * det_j.abs();
    if area < 1e-300 {
        return Err(RemError::Mesh(format!(
            "Degenerate Tri3 element {} (area ≈ 0)", elem.id
        )));
    }

    // Gradients of the three P1 basis functions (constant per element)
    // ∇λ_i = (1/(2A)) * [y_{j} - y_{k}, x_{k} - x_{j}]
    let inv2a = 1.0 / (2.0 * area);
    let grads = [
        [(y1 - y2) * inv2a, (x2 - x1) * inv2a],
        [(y2 - y0) * inv2a, (x0 - x2) * inv2a],
        [(y0 - y1) * inv2a, (x1 - x0) * inv2a],
    ];
    let nodes = [n0, n1, n2];

    // K_ij^e = eps * area * (∇λ_i · ∇λ_j)
    for i in 0..3 {
        for j in 0..3 {
            let k_ij = eps * area * (grads[i][0] * grads[j][0] + grads[i][1] * grads[j][1]);
            triplet.add(nodes[i], nodes[j], k_ij);
        }
    }
    Ok(())
}

/// Local stiffness for a linear tetrahedron (Tet4).
fn assemble_tet4(
    mesh: &RemMesh,
    elem: &rem_mesh::Element,
    eps: f64,
    triplet: &mut TripletMatrix,
) -> RemResult<()> {
    debug_assert_eq!(elem.node_ids.len(), 4);
    assemble_tet4_by_nodes(mesh, &elem.node_ids, elem.id, eps, triplet)
}

/// Assemble Tet4 stiffness from an explicit node slice (used for Tet10 corner approximation).
fn assemble_tet4_by_nodes(
    mesh: &RemMesh,
    node_ids: &[usize],
    elem_id: usize,
    eps: f64,
    triplet: &mut TripletMatrix,
) -> RemResult<()> {
    let [n0, n1, n2, n3] = [node_ids[0], node_ids[1], node_ids[2], node_ids[3]];
    let nodes = [n0, n1, n2, n3];
    let x = [
        mesh.nodes[n0].x, mesh.nodes[n1].x, mesh.nodes[n2].x, mesh.nodes[n3].x,
    ];
    let y = [
        mesh.nodes[n0].y, mesh.nodes[n1].y, mesh.nodes[n2].y, mesh.nodes[n3].y,
    ];
    let z = [
        mesh.nodes[n0].z, mesh.nodes[n1].z, mesh.nodes[n2].z, mesh.nodes[n3].z,
    ];

    // Jacobian columns: J = [x1-x0, x2-x0, x3-x0; y1-y0, ...; z1-z0, ...]
    let j = [
        [x[1]-x[0], x[2]-x[0], x[3]-x[0]],
        [y[1]-y[0], y[2]-y[0], y[3]-y[0]],
        [z[1]-z[0], z[2]-z[0], z[3]-z[0]],
    ];

    let det = det3(&j);
    let vol = det.abs() / 6.0;
    if vol < 1e-300 {
        return Err(RemError::Mesh(format!(
            "Degenerate Tet4 element {} (volume ≈ 0)", elem_id
        )));
    }

    // J^{-T} rows: the gradients of ξ, η, ζ w.r.t. x, y, z
    // Each column of J^{-1} gives ∇_x ξ, ∇_x η, ∇_x ζ
    let j_inv = inv3(&j, det);

    // Gradients in reference coordinates:
    // ∇_ξ λ_0 = (-1,-1,-1),  ∇_ξ λ_1 = (1,0,0),
    // ∇_ξ λ_2 = (0,1,0),     ∇_ξ λ_3 = (0,0,1)
    let ref_grads = [
        [-1.0f64, -1.0, -1.0],
        [ 1.0,     0.0,  0.0],
        [ 0.0,     1.0,  0.0],
        [ 0.0,     0.0,  1.0],
    ];

    // Physical gradients: ∇_x λ_i = J^{-T} * ∇_ξ λ_i
    let mut grads = [[0.0f64; 3]; 4];
    for i in 0..4 {
        for row in 0..3 {
            for col in 0..3 {
                grads[i][row] += j_inv[col][row] * ref_grads[i][col];
            }
        }
    }

    // K_ij^e = eps * vol * (∇λ_i · ∇λ_j)
    for i in 0..4 {
        for j in 0..4 {
            let k_ij = eps * vol * dot3(&grads[i], &grads[j]);
            triplet.add(nodes[i], nodes[j], k_ij);
        }
    }
    Ok(())
}

/// Local stiffness for a trilinear hexahedron (Hex8).
///
/// Uses 2×2×2 Gauss quadrature (exact for trilinear hex on parallelepipeds,
/// approximate for general hex). Each Gauss point carries weight 1.0.
fn assemble_hex8(
    mesh: &RemMesh,
    elem: &rem_mesh::Element,
    eps: f64,
    triplet: &mut TripletMatrix,
) -> RemResult<()> {
    debug_assert!(elem.node_ids.len() >= 8);
    let nids: [usize; 8] = [
        elem.node_ids[0], elem.node_ids[1], elem.node_ids[2], elem.node_ids[3],
        elem.node_ids[4], elem.node_ids[5], elem.node_ids[6], elem.node_ids[7],
    ];
    // GMSH Hex8 corner node order: same as standard right-hand hex
    // Nodes 0-3: bottom face (z-), nodes 4-7: top face (z+), each face in CCW order
    let xi_ref  = [-1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0_f64];
    let eta_ref = [-1.0,-1.0, 1.0,  1.0, -1.0,-1.0, 1.0,  1.0_f64];
    let zet_ref = [-1.0,-1.0,-1.0, -1.0,  1.0, 1.0, 1.0,  1.0_f64];

    let coords: [[f64; 3]; 8] = {
        let mut c = [[0.0f64; 3]; 8];
        for (i, &n) in nids.iter().enumerate() {
            c[i] = [mesh.nodes[n].x, mesh.nodes[n].y, mesh.nodes[n].z];
        }
        c
    };

    // 2-point Gauss: pts = ±1/√3, weight = 1.0
    let gp = 1.0_f64 / 3.0_f64.sqrt();
    let gauss_pts = [-gp, gp];

    let mut ke = [[0.0f64; 8]; 8];

    for &xi in &gauss_pts {
        for &eta in &gauss_pts {
            for &zet in &gauss_pts {
                // Shape functions and their reference-space gradients
                let mut n_val = [0.0f64; 8];
                let mut dn_dxi  = [0.0f64; 8];
                let mut dn_deta = [0.0f64; 8];
                let mut dn_dzet = [0.0f64; 8];
                for i in 0..8 {
                    let a = xi_ref[i];
                    let b = eta_ref[i];
                    let c = zet_ref[i];
                    n_val[i]  = 0.125 * (1.0 + a*xi) * (1.0 + b*eta) * (1.0 + c*zet);
                    dn_dxi[i]  = 0.125 * a * (1.0 + b*eta) * (1.0 + c*zet);
                    dn_deta[i] = 0.125 * b * (1.0 + a*xi)  * (1.0 + c*zet);
                    dn_dzet[i] = 0.125 * c * (1.0 + a*xi)  * (1.0 + b*eta);
                }

                // Jacobian J_kl = d(x_k)/d(xi_l)
                let mut jac = [[0.0f64; 3]; 3];
                for i in 0..8 {
                    jac[0][0] += dn_dxi[i]  * coords[i][0];
                    jac[0][1] += dn_deta[i] * coords[i][0];
                    jac[0][2] += dn_dzet[i] * coords[i][0];
                    jac[1][0] += dn_dxi[i]  * coords[i][1];
                    jac[1][1] += dn_deta[i] * coords[i][1];
                    jac[1][2] += dn_dzet[i] * coords[i][1];
                    jac[2][0] += dn_dxi[i]  * coords[i][2];
                    jac[2][1] += dn_deta[i] * coords[i][2];
                    jac[2][2] += dn_dzet[i] * coords[i][2];
                }
                let det_j = det3(&jac);
                if det_j.abs() < 1e-300 { continue; }
                let j_inv = inv3(&jac, det_j);

                // Physical gradients: ∇φ_i = J^{-T} * [dN/dxi, dN/deta, dN/dzet]
                let mut grad = [[0.0f64; 3]; 8];
                for i in 0..8 {
                    let ref_g = [dn_dxi[i], dn_deta[i], dn_dzet[i]];
                    for row in 0..3 {
                        for col in 0..3 {
                            grad[i][row] += j_inv[col][row] * ref_g[col];
                        }
                    }
                }

                // Accumulate K_ij += eps * w * det(J) * ∇φ_i · ∇φ_j  (w=1 for 2-pt Gauss)
                let wdet = eps * det_j.abs();
                for i in 0..8 {
                    for j in 0..8 {
                        ke[i][j] += wdet * dot3(&grad[i], &grad[j]);
                    }
                }
            }
        }
    }

    for i in 0..8 {
        for j in 0..8 {
            triplet.add(nids[i], nids[j], ke[i][j]);
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// P2 quadratic triangle (Tri6) stiffness assembly
// ---------------------------------------------------------------------------

/// Local stiffness for a quadratic triangle (Tri6) using P2 basis functions.
///
/// Reference element: ξ∈[0,1], η∈[0,1-ξ].  Barycentric: λ1=1-ξ-η, λ2=ξ, λ3=η.
///
/// GMSH node order:
///   0:(0,0), 1:(1,0), 2:(0,1), 3:(1/2,0), 4:(1/2,1/2), 5:(0,1/2)
///
/// P2 basis: vertex φi = λi(2λi-1), edge φ3=4λ1λ2, φ4=4λ2λ3, φ5=4λ1λ3.
///
/// Gauss quadrature: 3-point degree-2 rule (exact for ∇φi·∇φj which is degree 2).
fn assemble_tri6(
    mesh: &RemMesh,
    elem: &rem_mesh::Element,
    eps: f64,
    triplet: &mut TripletMatrix,
) -> RemResult<()> {
    debug_assert_eq!(elem.node_ids.len(), 6);
    let nids: [usize; 6] = [
        elem.node_ids[0], elem.node_ids[1], elem.node_ids[2],
        elem.node_ids[3], elem.node_ids[4], elem.node_ids[5],
    ];

    // Nodal coordinates
    let xy: [[f64; 2]; 6] = {
        let mut c = [[0.0f64; 2]; 6];
        for (i, &n) in nids.iter().enumerate() {
            c[i] = [mesh.nodes[n].x, mesh.nodes[n].y];
        }
        c
    };

    // 3-point Gauss rule on reference triangle (exact for degree ≤ 2)
    // Points: (1/6,1/6), (2/3,1/6), (1/6,2/3);  weight = 1/6 each
    let gauss_pts = [[1.0/6.0, 1.0/6.0], [2.0/3.0, 1.0/6.0], [1.0/6.0, 2.0/3.0]];
    let gauss_w = 1.0 / 6.0;

    let mut ke = [[0.0f64; 6]; 6];

    for &[xi, eta] in &gauss_pts {
        let l1 = 1.0 - xi - eta;
        let l2 = xi;
        let l3 = eta;

        // P2 basis gradients in reference coords (∂φ/∂ξ, ∂φ/∂η)
        // Corner nodes: ∂(λ(2λ-1))/∂α = (4λ-1) ∂λ/∂α
        // ∂l1/∂ξ=-1, ∂l1/∂η=-1; ∂l2/∂ξ=1,∂l2/∂η=0; ∂l3/∂ξ=0,∂l3/∂η=1
        let dndxi  = [
            -(4.0*l1 - 1.0),          // φ0 = l1(2l1-1)
             4.0*l2 - 1.0,            // φ1 = l2(2l2-1)
             0.0,                     // φ2 = l3(2l3-1)
             4.0*(l1 - l2),           // φ3 = 4l1l2: ∂/∂ξ = 4(∂l1/∂ξ·l2 + l1·∂l2/∂ξ) = 4(-l2+l1)
             4.0*l3,                  // φ4 = 4l2l3: 4l3
            -4.0*l3,                  // φ5 = 4l1l3: 4(∂l1/∂ξ·l3) = -4l3
        ];
        let dndeta = [
            -(4.0*l1 - 1.0),          // φ0
             0.0,                     // φ1
             4.0*l3 - 1.0,            // φ2
            -4.0*l2,                  // φ3 = 4l1l2: 4·∂l1/∂η·l2 = -4l2
             4.0*l2,                  // φ4 = 4l2l3: 4l2
             4.0*(l1 - l3),           // φ5 = 4l1l3: 4(-l3+l1)
        ];

        // Jacobian J = [∂x/∂ξ, ∂x/∂η; ∂y/∂ξ, ∂y/∂η]
        let mut jac = [[0.0f64; 2]; 2];
        for i in 0..6 {
            jac[0][0] += dndxi[i]  * xy[i][0]; // ∂x/∂ξ
            jac[0][1] += dndeta[i] * xy[i][0]; // ∂x/∂η
            jac[1][0] += dndxi[i]  * xy[i][1]; // ∂y/∂ξ
            jac[1][1] += dndeta[i] * xy[i][1]; // ∂y/∂η
        }
        let det_j = jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0];
        if det_j.abs() < 1e-300 { continue; }

        // J^{-1}: inv([[a,b],[c,d]]) = (1/det) [[d,-b],[-c,a]]
        let inv_det = 1.0 / det_j;
        let ji = [
            [ jac[1][1] * inv_det, -jac[0][1] * inv_det],
            [-jac[1][0] * inv_det,  jac[0][0] * inv_det],
        ];

        // Physical gradients: [∂φ/∂x, ∂φ/∂y] = J^{-T} [∂φ/∂ξ, ∂φ/∂η]
        let mut grad = [[0.0f64; 2]; 6];
        for i in 0..6 {
            grad[i][0] = ji[0][0] * dndxi[i] + ji[1][0] * dndeta[i];
            grad[i][1] = ji[0][1] * dndxi[i] + ji[1][1] * dndeta[i];
        }

        let wdet = eps * gauss_w * det_j.abs();
        for i in 0..6 {
            for j in 0..6 {
                ke[i][j] += wdet * (grad[i][0]*grad[j][0] + grad[i][1]*grad[j][1]);
            }
        }
    }

    for i in 0..6 {
        for j in 0..6 {
            triplet.add(nids[i], nids[j], ke[i][j]);
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// P2 quadratic tetrahedron (Tet10) stiffness assembly
// ---------------------------------------------------------------------------

/// Local stiffness for a quadratic tetrahedron (Tet10) using P2 basis functions.
///
/// Reference element: ξ,η,ζ≥0, ξ+η+ζ≤1.  Barycentric: λ1=1-ξ-η-ζ, λ2=ξ, λ3=η, λ4=ζ.
///
/// GMSH node order (0-indexed):
///   0:(0,0,0), 1:(1,0,0), 2:(0,1,0), 3:(0,0,1),
///   4:mid(0-1), 5:mid(1-2), 6:mid(0-2), 7:mid(0-3), 8:mid(1-3), 9:mid(2-3)
///
/// P2 basis: vertex φi = λi(2λi-1), edge φij = 4λiλj.
///
/// Gauss quadrature: Keast 4-point degree-2 rule (exact for ∇φi·∇φj which is degree 2).
fn assemble_tet10(
    mesh: &RemMesh,
    elem: &rem_mesh::Element,
    eps: f64,
    triplet: &mut TripletMatrix,
) -> RemResult<()> {
    debug_assert_eq!(elem.node_ids.len(), 10);
    let nids: [usize; 10] = {
        let ids = &elem.node_ids;
        [ids[0],ids[1],ids[2],ids[3],ids[4],ids[5],ids[6],ids[7],ids[8],ids[9]]
    };

    let xyz: [[f64; 3]; 10] = {
        let mut c = [[0.0f64; 3]; 10];
        for (i, &n) in nids.iter().enumerate() {
            c[i] = [mesh.nodes[n].x, mesh.nodes[n].y, mesh.nodes[n].z];
        }
        c
    };

    // Keast 4-point rule on reference tet (exact for degree 2).
    // a ≈ 0.5854101966, b ≈ 0.1381966011; weight = 1/24 each.
    const A: f64 = 0.5854101966249685;
    const B: f64 = 0.1381966011250105;
    const W: f64 = 1.0 / 24.0;
    let gauss_pts = [[B,B,B],[A,B,B],[B,A,B],[B,B,A]];

    let mut ke = [[0.0f64; 10]; 10];

    for &[xi, eta, zet] in &gauss_pts {
        let l1 = 1.0 - xi - eta - zet; // λ1
        // l2 = xi, l3 = eta, l4 = zet

        // P2 basis gradients in reference coords (∂φ/∂ξ, ∂φ/∂η, ∂φ/∂ζ).
        // Corner nodes: ∂(λ(2λ-1))/∂α = (4λ-1)·∂λ/∂α
        //   ∂l1/∂ξ=∂l1/∂η=∂l1/∂ζ=-1
        //   ∂l2/∂ξ=1; ∂l3/∂η=1; ∂l4/∂ζ=1; all others 0
        // Edge nodes φ_ab = 4λaλb: ∂φ/∂α = 4(∂λa/∂α·λb + λa·∂λb/∂α)
        let dndxi = [
            -(4.0*l1 - 1.0),              // φ0: (4l1-1)*(-1)
             4.0*xi - 1.0,                // φ1: (4l2-1)*1
             0.0,                         // φ2
             0.0,                         // φ3
             4.0*(l1 - xi),               // φ4=4l1l2: 4(-l2+l1)
             4.0*eta,                     // φ5=4l2l3: 4l3
            -4.0*eta,                     // φ6=4l1l3: 4(-l3)
            -4.0*zet,                     // φ7=4l1l4: 4(-l4)
             4.0*zet,                     // φ8=4l2l4: 4l4
             0.0,                         // φ9=4l3l4: 0
        ];
        let dndeta = [
            -(4.0*l1 - 1.0),              // φ0
             0.0,                         // φ1
             4.0*eta - 1.0,               // φ2: (4l3-1)*1
             0.0,                         // φ3
            -4.0*xi,                      // φ4=4l1l2: 4(-l2)
             4.0*xi,                      // φ5=4l2l3: 4l2
             4.0*(l1 - eta),              // φ6=4l1l3: 4(-l3+l1)
            -4.0*zet,                     // φ7=4l1l4: 4(-l4)
             0.0,                         // φ8=4l2l4: 0
             4.0*zet,                     // φ9=4l3l4: 4l4
        ];
        let dndzet = [
            -(4.0*l1 - 1.0),              // φ0
             0.0,                         // φ1
             0.0,                         // φ2
             4.0*zet - 1.0,               // φ3: (4l4-1)*1
            -4.0*xi,                      // φ4=4l1l2: 4(-l2) [no l4 term]... actually
             0.0,                         // φ5=4l2l3: 0
            -4.0*eta,                     // φ6=4l1l3: 4(-l3)
             4.0*(l1 - zet),              // φ7=4l1l4: 4(-l4+l1)
             4.0*xi,                      // φ8=4l2l4: 4l2
             4.0*eta,                     // φ9=4l3l4: 4l3
        ];

        // Jacobian J_kl = ∂x_k/∂ξ_l  (k=row, l=col: ξ,η,ζ)
        let mut jac = [[0.0f64; 3]; 3];
        for i in 0..10 {
            let dref = [dndxi[i], dndeta[i], dndzet[i]];
            for k in 0..3 {
                jac[k][0] += dref[0] * xyz[i][k];
                jac[k][1] += dref[1] * xyz[i][k];
                jac[k][2] += dref[2] * xyz[i][k];
            }
        }
        let det_j = det3(&jac);
        if det_j.abs() < 1e-300 { continue; }
        let j_inv = inv3(&jac, det_j);

        // Physical gradients: ∇φ_i = J^{-T} ∇_ξ φ_i
        let mut grad = [[0.0f64; 3]; 10];
        for i in 0..10 {
            let ref_g = [dndxi[i], dndeta[i], dndzet[i]];
            for row in 0..3 {
                for col in 0..3 {
                    grad[i][row] += j_inv[col][row] * ref_g[col];
                }
            }
        }

        let wdet = eps * W * det_j.abs();
        for i in 0..10 {
            for j in 0..10 {
                ke[i][j] += wdet * dot3(&grad[i], &grad[j]);
            }
        }
    }

    for i in 0..10 {
        for j in 0..10 {
            triplet.add(nids[i], nids[j], ke[i][j]);
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// 3×3 matrix helpers (pub for use in postprocess.rs)
// ---------------------------------------------------------------------------

pub fn det3_pub(m: &[[f64; 3]; 3]) -> f64 { det3(m) }
pub fn inv3_pub(m: &[[f64; 3]; 3], det: f64) -> [[f64; 3]; 3] { inv3(m, det) }

fn det3(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
  - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
  + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

fn inv3(m: &[[f64; 3]; 3], det: f64) -> [[f64; 3]; 3] {
    let inv_det = 1.0 / det;
    [
        [
            (m[1][1]*m[2][2] - m[1][2]*m[2][1]) * inv_det,
            (m[0][2]*m[2][1] - m[0][1]*m[2][2]) * inv_det,
            (m[0][1]*m[1][2] - m[0][2]*m[1][1]) * inv_det,
        ],
        [
            (m[1][2]*m[2][0] - m[1][0]*m[2][2]) * inv_det,
            (m[0][0]*m[2][2] - m[0][2]*m[2][0]) * inv_det,
            (m[0][2]*m[1][0] - m[0][0]*m[1][2]) * inv_det,
        ],
        [
            (m[1][0]*m[2][1] - m[1][1]*m[2][0]) * inv_det,
            (m[0][1]*m[2][0] - m[0][0]*m[2][1]) * inv_det,
            (m[0][0]*m[1][1] - m[0][1]*m[1][0]) * inv_det,
        ],
    ]
}

#[inline]
fn dot3(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
}

/// Assemble the global stiffness matrix with a full anisotropic (tensor) coefficient.
///
/// `tensor_fn` maps physical group tag → absolute permittivity tensor [F/m], 3×3 row-major.
/// Element integral: K_ij^e = ∫_Ωe  ∇λ_i^T · A · ∇λ_j  dΩ
///
/// Falls back to the 2-D xy submatrix for Tri3 elements.
pub fn assemble_stiffness_aniso(
    mesh: &RemMesh,
    tensor_fn: impl Fn(u32) -> [[f64; 3]; 3],
) -> RemResult<TripletMatrix> {
    let n = mesh.n_nodes();
    let cap = mesh.n_volume_elements() * 16;
    let mut triplet = TripletMatrix::with_capacity(n, n, cap);

    for elem in &mesh.volume_elements {
        if mesh.size > 1 && elem.rank != mesh.rank {
            continue;
        }
        let a = tensor_fn(elem.tag);
        match elem.kind {
            ElementKind::Tri3 => {
                assemble_tri3_aniso(mesh, elem, &a, &mut triplet)?;
            }
            ElementKind::Tri6 => {
                assemble_tri6_aniso(mesh, elem, &a, &mut triplet)?;
            }
            ElementKind::Tet4 => {
                assemble_tet4_aniso_by_nodes(mesh, &elem.node_ids, elem.id, &a, &mut triplet)?;
            }
            ElementKind::Tet10 => {
                assemble_tet10_aniso(mesh, elem, &a, &mut triplet)?;
            }
            other => {
                log::warn!(
                    "Element kind {:?} not supported in anisotropic P1 assembly — skipping",
                    other
                );
            }
        }
    }

    Ok(triplet)
}

/// Anisotropic Tri3: uses the 2×2 upper-left submatrix of `a`.
fn assemble_tri3_aniso(
    mesh: &RemMesh,
    elem: &rem_mesh::Element,
    a: &[[f64; 3]; 3],
    triplet: &mut TripletMatrix,
) -> RemResult<()> {
    debug_assert_eq!(elem.node_ids.len(), 3);
    let [n0, n1, n2] = [elem.node_ids[0], elem.node_ids[1], elem.node_ids[2]];
    let (x0, y0) = (mesh.nodes[n0].x, mesh.nodes[n0].y);
    let (x1, y1) = (mesh.nodes[n1].x, mesh.nodes[n1].y);
    let (x2, y2) = (mesh.nodes[n2].x, mesh.nodes[n2].y);

    let det_j = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
    let area = 0.5 * det_j.abs();
    if area < 1e-300 {
        return Err(RemError::Mesh(format!("Degenerate Tri3 element {} (area ≈ 0)", elem.id)));
    }

    let inv2a = 1.0 / (2.0 * area);
    // Only x,y components for 2-D
    let grads: [[f64; 2]; 3] = [
        [(y1 - y2) * inv2a, (x2 - x1) * inv2a],
        [(y2 - y0) * inv2a, (x0 - x2) * inv2a],
        [(y0 - y1) * inv2a, (x1 - x0) * inv2a],
    ];
    let nodes = [n0, n1, n2];

    // K_ij = area * g_i^T · A_2d · g_j
    for i in 0..3 {
        for j in 0..3 {
            let k_ij = area * (
                grads[i][0] * (a[0][0]*grads[j][0] + a[0][1]*grads[j][1])
              + grads[i][1] * (a[1][0]*grads[j][0] + a[1][1]*grads[j][1])
            );
            triplet.add(nodes[i], nodes[j], k_ij);
        }
    }
    Ok(())
}

/// Anisotropic Tet4 stiffness from explicit node slice.
fn assemble_tet4_aniso_by_nodes(
    mesh: &RemMesh,
    node_ids: &[usize],
    elem_id: usize,
    a: &[[f64; 3]; 3],
    triplet: &mut TripletMatrix,
) -> RemResult<()> {
    let [n0, n1, n2, n3] = [node_ids[0], node_ids[1], node_ids[2], node_ids[3]];
    let nodes = [n0, n1, n2, n3];
    let x = [mesh.nodes[n0].x, mesh.nodes[n1].x, mesh.nodes[n2].x, mesh.nodes[n3].x];
    let y = [mesh.nodes[n0].y, mesh.nodes[n1].y, mesh.nodes[n2].y, mesh.nodes[n3].y];
    let z = [mesh.nodes[n0].z, mesh.nodes[n1].z, mesh.nodes[n2].z, mesh.nodes[n3].z];

    let jac = [
        [x[1]-x[0], x[2]-x[0], x[3]-x[0]],
        [y[1]-y[0], y[2]-y[0], y[3]-y[0]],
        [z[1]-z[0], z[2]-z[0], z[3]-z[0]],
    ];
    let det = det3(&jac);
    let vol = det.abs() / 6.0;
    if vol < 1e-300 {
        return Err(RemError::Mesh(format!("Degenerate Tet4 element {} (volume ≈ 0)", elem_id)));
    }

    let j_inv = inv3(&jac, det);
    let ref_grads = [[-1.0f64,-1.0,-1.0],[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]];
    let mut grads = [[0.0f64; 3]; 4];
    for i in 0..4 {
        for row in 0..3 {
            for col in 0..3 {
                grads[i][row] += j_inv[col][row] * ref_grads[i][col];
            }
        }
    }

    // K_ij = vol * g_i^T · A · g_j
    for i in 0..4 {
        for j in 0..4 {
            let mut k_ij = 0.0;
            for r in 0..3 {
                let ag_j_r: f64 = a[r][0]*grads[j][0] + a[r][1]*grads[j][1] + a[r][2]*grads[j][2];
                k_ij += grads[i][r] * ag_j_r;
            }
            triplet.add(nodes[i], nodes[j], vol * k_ij);
        }
    }
    Ok(())
}

/// Anisotropic Tri6 P2 stiffness: uses 3-point Gauss quadrature with tensor coefficient.
fn assemble_tri6_aniso(
    mesh: &RemMesh,
    elem: &rem_mesh::Element,
    a: &[[f64; 3]; 3],
    triplet: &mut TripletMatrix,
) -> RemResult<()> {
    debug_assert_eq!(elem.node_ids.len(), 6);
    let nids: [usize; 6] = {
        let ids = &elem.node_ids;
        [ids[0],ids[1],ids[2],ids[3],ids[4],ids[5]]
    };
    let xy: [[f64; 2]; 6] = {
        let mut c = [[0.0f64; 2]; 6];
        for (i, &n) in nids.iter().enumerate() {
            c[i] = [mesh.nodes[n].x, mesh.nodes[n].y];
        }
        c
    };

    let gauss_pts = [[1.0/6.0, 1.0/6.0], [2.0/3.0, 1.0/6.0], [1.0/6.0, 2.0/3.0]];
    let gauss_w = 1.0 / 6.0;
    let mut ke = [[0.0f64; 6]; 6];

    for &[xi, eta] in &gauss_pts {
        let l1 = 1.0 - xi - eta;
        let l2 = xi;
        let l3 = eta;
        let dndxi  = [-(4.0*l1-1.0), 4.0*l2-1.0, 0.0, 4.0*(l1-l2), 4.0*l3, -4.0*l3];
        let dndeta = [-(4.0*l1-1.0), 0.0, 4.0*l3-1.0, -4.0*l2, 4.0*l2, 4.0*(l1-l3)];

        let mut jac = [[0.0f64; 2]; 2];
        for i in 0..6 {
            jac[0][0] += dndxi[i]*xy[i][0]; jac[0][1] += dndeta[i]*xy[i][0];
            jac[1][0] += dndxi[i]*xy[i][1]; jac[1][1] += dndeta[i]*xy[i][1];
        }
        let det_j = jac[0][0]*jac[1][1] - jac[0][1]*jac[1][0];
        if det_j.abs() < 1e-300 { continue; }
        let inv_det = 1.0 / det_j;
        let ji = [
            [ jac[1][1]*inv_det, -jac[0][1]*inv_det],
            [-jac[1][0]*inv_det,  jac[0][0]*inv_det],
        ];
        let mut grad = [[0.0f64; 3]; 6];
        for i in 0..6 {
            grad[i][0] = ji[0][0]*dndxi[i] + ji[1][0]*dndeta[i];
            grad[i][1] = ji[0][1]*dndxi[i] + ji[1][1]*dndeta[i];
            // z-component is 0 (2D element)
        }
        let wdet = gauss_w * det_j.abs();
        for i in 0..6 {
            for j in 0..6 {
                let mut k_ij = 0.0;
                for r in 0..3 {
                    let ag_r = a[r][0]*grad[j][0] + a[r][1]*grad[j][1] + a[r][2]*grad[j][2];
                    k_ij += grad[i][r] * ag_r;
                }
                ke[i][j] += wdet * k_ij;
            }
        }
    }
    for i in 0..6 {
        for j in 0..6 {
            triplet.add(nids[i], nids[j], ke[i][j]);
        }
    }
    Ok(())
}

/// Anisotropic Tet10 P2 stiffness: uses Keast 4-point Gauss quadrature with tensor coefficient.
fn assemble_tet10_aniso(
    mesh: &RemMesh,
    elem: &rem_mesh::Element,
    a: &[[f64; 3]; 3],
    triplet: &mut TripletMatrix,
) -> RemResult<()> {
    debug_assert_eq!(elem.node_ids.len(), 10);
    let nids: [usize; 10] = {
        let ids = &elem.node_ids;
        [ids[0],ids[1],ids[2],ids[3],ids[4],ids[5],ids[6],ids[7],ids[8],ids[9]]
    };
    let xyz: [[f64; 3]; 10] = {
        let mut c = [[0.0f64; 3]; 10];
        for (i, &n) in nids.iter().enumerate() {
            c[i] = [mesh.nodes[n].x, mesh.nodes[n].y, mesh.nodes[n].z];
        }
        c
    };

    const AA: f64 = 0.5854101966249685;
    const BB: f64 = 0.1381966011250105;
    const WW: f64 = 1.0 / 24.0;
    let gauss_pts = [[BB,BB,BB],[AA,BB,BB],[BB,AA,BB],[BB,BB,AA]];

    let mut ke = [[0.0f64; 10]; 10];

    for &[xi, eta, zet] in &gauss_pts {
        let l1 = 1.0 - xi - eta - zet;
        let dndxi = [-(4.0*l1-1.0), 4.0*xi-1.0, 0.0, 0.0,
                      4.0*(l1-xi), 4.0*eta, -4.0*eta, -4.0*zet, 4.0*zet, 0.0];
        let dndeta = [-(4.0*l1-1.0), 0.0, 4.0*eta-1.0, 0.0,
                       -4.0*xi, 4.0*xi, 4.0*(l1-eta), -4.0*zet, 0.0, 4.0*zet];
        let dndzet = [-(4.0*l1-1.0), 0.0, 0.0, 4.0*zet-1.0,
                       -4.0*xi, 0.0, -4.0*eta, 4.0*(l1-zet), 4.0*xi, 4.0*eta];

        let mut jac = [[0.0f64; 3]; 3];
        for i in 0..10 {
            let dref = [dndxi[i], dndeta[i], dndzet[i]];
            for k in 0..3 {
                jac[k][0] += dref[0]*xyz[i][k];
                jac[k][1] += dref[1]*xyz[i][k];
                jac[k][2] += dref[2]*xyz[i][k];
            }
        }
        let det_j = det3(&jac);
        if det_j.abs() < 1e-300 { continue; }
        let j_inv = inv3(&jac, det_j);

        let mut grad = [[0.0f64; 3]; 10];
        for i in 0..10 {
            let ref_g = [dndxi[i], dndeta[i], dndzet[i]];
            for row in 0..3 {
                for col in 0..3 {
                    grad[i][row] += j_inv[col][row] * ref_g[col];
                }
            }
        }
        let wdet = WW * det_j.abs();
        for i in 0..10 {
            for j in 0..10 {
                let mut k_ij = 0.0;
                for r in 0..3 {
                    let ag_r = a[r][0]*grad[j][0] + a[r][1]*grad[j][1] + a[r][2]*grad[j][2];
                    k_ij += grad[i][r] * ag_r;
                }
                ke[i][j] += wdet * k_ij;
            }
        }
    }
    for i in 0..10 {
        for j in 0..10 {
            triplet.add(nids[i], nids[j], ke[i][j]);
        }
    }
    Ok(())
}

/// Build per-element epsilon from the domain map.
pub fn element_epsilon(mesh: &RemMesh, domain_map: &DomainMap) -> Vec<f64> {
    mesh.volume_elements
        .iter()
        .map(|e| domain_map.get(e.tag).epsilon_abs())
        .collect()
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
pub mod tests {
    use super::*;
    use rem_core::CsrMatrix;

    /// Build a trivial 2-node, 1-triangle mesh for testing.
    /// Triangle: (0,0), (1,0), (0,1) — tag 1
    pub fn unit_triangle_mesh() -> RemMesh {
        use rem_mesh::{Node, Element};
        RemMesh {
            nodes: vec![
                Node { id: 0, x: 0.0, y: 0.0, z: 0.0 },
                Node { id: 1, x: 1.0, y: 0.0, z: 0.0 },
                Node { id: 2, x: 0.0, y: 1.0, z: 0.0 },
            ],
            volume_elements: vec![
                Element { id: 1, kind: ElementKind::Tri3, tag: 1, node_ids: vec![0, 1, 2] , rank: 0 },
            ],
            boundary_elements: vec![],
            domain_tags: Default::default(),
            boundary_tags: Default::default(),
            dim: 2,
            rank: 0,
            size: 1,
        }
    }

    #[test]
    fn tri3_stiffness_row_sum_zero() {
        // For a constant-coefficient problem, the row sums of K must be zero
        // (partition of unity: Σ_j K_ij = 0 for interior nodes).
        let mesh = unit_triangle_mesh();
        let triplet = assemble_stiffness(&mesh, |_| 1.0).unwrap();
        let csr = triplet.to_csr();
        let n = mesh.n_nodes();
        let mut x = vec![1.0; n];
        let mut y = vec![0.0; n];
        csr.matvec(&x, &mut y, &rem_parallel::NoComm);
        for &yi in &y {
            assert!(yi.abs() < 1e-13, "row sum = {}", yi);
        }
    }

    #[test]
    fn tri3_stiffness_symmetry() {
        let mesh = unit_triangle_mesh();
        let triplet = assemble_stiffness(&mesh, |_| 2.5).unwrap();
        let csr = triplet.to_csr();
        // K[i,j] == K[j,i] for all i,j
        let n = csr.nrows;
        for i in 0..n {
            for k in csr.row_ptr[i]..csr.row_ptr[i + 1] {
                let j = csr.col_idx[k];
                let kij = csr.values[k];
                // find K[j,i]
                let kji = {
                    let mut v = None;
                    for kk in csr.row_ptr[j]..csr.row_ptr[j + 1] {
                        if csr.col_idx[kk] == i { v = Some(csr.values[kk]); break; }
                    }
                    v.unwrap_or(0.0)
                };
                assert!((kij - kji).abs() < 1e-13, "K[{},{}]={} != K[{},{}]={}", i, j, kij, j, i, kji);
            }
        }
    }

    /// Anisotropic stiffness with identity tensor must equal scalar stiffness.
    #[test]
    fn aniso_identity_tensor_matches_scalar() {
        let mesh = unit_triangle_mesh();
        let eps = 3.0_f64;
        let tensor = [[eps, 0.0, 0.0], [0.0, eps, 0.0], [0.0, 0.0, eps]];
        let aniso = assemble_stiffness_aniso(&mesh, |_| tensor).unwrap().to_csr();
        let scalar = assemble_stiffness(&mesh, |_| eps).unwrap().to_csr();

        assert_eq!(aniso.nrows, scalar.nrows);
        assert_eq!(aniso.nnz(), scalar.nnz());
        for (a, b) in aniso.values.iter().zip(scalar.values.iter()) {
            assert!((a - b).abs() < 1e-12, "aniso={} scalar={}", a, b);
        }
    }

    /// Anisotropic row sum must be zero (constant fields in null space).
    #[test]
    fn aniso_row_sum_zero() {
        let mesh = unit_triangle_mesh();
        // Off-diagonal tensor to test genuine anisotropy
        let tensor = [[2.0, 0.5, 0.0], [0.5, 3.0, 0.0], [0.0, 0.0, 1.0]];
        let triplet = assemble_stiffness_aniso(&mesh, |_| tensor).unwrap();
        let csr = triplet.to_csr();
        let n = mesh.n_nodes();
        let x = vec![1.0; n];
        let mut y = vec![0.0; n];
        csr.matvec(&x, &mut y, &rem_parallel::NoComm);
        for &yi in &y {
            assert!(yi.abs() < 1e-12, "aniso row sum = {}", yi);
        }
    }

    /// Anisotropic stiffness must be symmetric for symmetric tensor.
    #[test]
    fn aniso_symmetry() {
        let mesh = unit_triangle_mesh();
        let tensor = [[4.0, 1.0, 0.0], [1.0, 2.0, 0.0], [0.0, 0.0, 3.0]];
        let csr = assemble_stiffness_aniso(&mesh, |_| tensor).unwrap().to_csr();
        let n = csr.nrows;
        for i in 0..n {
            for k in csr.row_ptr[i]..csr.row_ptr[i + 1] {
                let j = csr.col_idx[k];
                let kij = csr.values[k];
                let kji = {
                    let mut v = None;
                    for kk in csr.row_ptr[j]..csr.row_ptr[j + 1] {
                        if csr.col_idx[kk] == i { v = Some(csr.values[kk]); break; }
                    }
                    v.unwrap_or(0.0)
                };
                assert!((kij - kji).abs() < 1e-12, "K_aniso[{},{}]={} != K_aniso[{},{}]={}", i, j, kij, j, i, kji);
            }
        }
    }

    // -----------------------------------------------------------------------
    // P2 element tests
    // -----------------------------------------------------------------------

    /// Build a single Tri6 element mesh from the reference triangle.
    /// Nodes: 0:(0,0), 1:(1,0), 2:(0,1), 3:(0.5,0), 4:(0.5,0.5), 5:(0,0.5)
    fn unit_tri6_mesh() -> RemMesh {
        use rem_mesh::{Node, Element};
        RemMesh {
            nodes: vec![
                Node { id: 0, x: 0.0, y: 0.0, z: 0.0 },
                Node { id: 1, x: 1.0, y: 0.0, z: 0.0 },
                Node { id: 2, x: 0.0, y: 1.0, z: 0.0 },
                Node { id: 3, x: 0.5, y: 0.0, z: 0.0 }, // mid 0-1
                Node { id: 4, x: 0.5, y: 0.5, z: 0.0 }, // mid 1-2
                Node { id: 5, x: 0.0, y: 0.5, z: 0.0 }, // mid 0-2
            ],
            volume_elements: vec![
                Element { id: 1, kind: ElementKind::Tri6, tag: 1, node_ids: vec![0,1,2,3,4,5], rank: 0 },
            ],
            boundary_elements: vec![],
            domain_tags: Default::default(),
            boundary_tags: Default::default(),
            dim: 2, rank: 0, size: 1,
        }
    }

    /// Build a single Tet10 element mesh from the reference tetrahedron.
    /// Nodes 0-3: corners; 4-9: edge midpoints (GMSH order).
    fn unit_tet10_mesh() -> RemMesh {
        use rem_mesh::{Node, Element};
        RemMesh {
            nodes: vec![
                Node { id: 0, x: 0.0, y: 0.0, z: 0.0 }, // corner 0
                Node { id: 1, x: 1.0, y: 0.0, z: 0.0 }, // corner 1
                Node { id: 2, x: 0.0, y: 1.0, z: 0.0 }, // corner 2
                Node { id: 3, x: 0.0, y: 0.0, z: 1.0 }, // corner 3
                Node { id: 4, x: 0.5, y: 0.0, z: 0.0 }, // mid(0-1)
                Node { id: 5, x: 0.5, y: 0.5, z: 0.0 }, // mid(1-2)
                Node { id: 6, x: 0.0, y: 0.5, z: 0.0 }, // mid(0-2)
                Node { id: 7, x: 0.0, y: 0.0, z: 0.5 }, // mid(0-3)
                Node { id: 8, x: 0.5, y: 0.0, z: 0.5 }, // mid(1-3)
                Node { id: 9, x: 0.0, y: 0.5, z: 0.5 }, // mid(2-3)
            ],
            volume_elements: vec![
                Element { id: 1, kind: ElementKind::Tet10, tag: 1,
                    node_ids: vec![0,1,2,3,4,5,6,7,8,9], rank: 0 },
            ],
            boundary_elements: vec![],
            domain_tags: Default::default(),
            boundary_tags: Default::default(),
            dim: 3, rank: 0, size: 1,
        }
    }

    #[test]
    fn tri6_stiffness_row_sum_zero() {
        let mesh = unit_tri6_mesh();
        let triplet = assemble_stiffness(&mesh, |_| 1.0).unwrap();
        let csr = triplet.to_csr();
        let n = mesh.n_nodes();
        let x = vec![1.0; n];
        let mut y = vec![0.0; n];
        csr.matvec(&x, &mut y, &rem_parallel::NoComm);
        for (i, &yi) in y.iter().enumerate() {
            assert!(yi.abs() < 1e-12, "Tri6 row {} sum = {:.3e}", i, yi);
        }
    }

    #[test]
    fn tri6_stiffness_symmetry() {
        let mesh = unit_tri6_mesh();
        let csr = assemble_stiffness(&mesh, |_| 1.0).unwrap().to_csr();
        let n = csr.nrows;
        for i in 0..n {
            for k in csr.row_ptr[i]..csr.row_ptr[i+1] {
                let j = csr.col_idx[k];
                let kij = csr.values[k];
                let kji = (csr.row_ptr[j]..csr.row_ptr[j+1])
                    .find(|&kk| csr.col_idx[kk] == i)
                    .map(|kk| csr.values[kk]).unwrap_or(0.0);
                assert!((kij - kji).abs() < 1e-12,
                    "Tri6 K[{},{}]={:.4e} != K[{},{}]={:.4e}", i, j, kij, j, i, kji);
            }
        }
    }

    #[test]
    fn tet10_stiffness_row_sum_zero() {
        let mesh = unit_tet10_mesh();
        let triplet = assemble_stiffness(&mesh, |_| 1.0).unwrap();
        let csr = triplet.to_csr();
        let n = mesh.n_nodes();
        let x = vec![1.0; n];
        let mut y = vec![0.0; n];
        csr.matvec(&x, &mut y, &rem_parallel::NoComm);
        for (i, &yi) in y.iter().enumerate() {
            assert!(yi.abs() < 1e-12, "Tet10 row {} sum = {:.3e}", i, yi);
        }
    }

    #[test]
    fn tet10_stiffness_symmetry() {
        let mesh = unit_tet10_mesh();
        let csr = assemble_stiffness(&mesh, |_| 2.0).unwrap().to_csr();
        let n = csr.nrows;
        for i in 0..n {
            for k in csr.row_ptr[i]..csr.row_ptr[i+1] {
                let j = csr.col_idx[k];
                let kij = csr.values[k];
                let kji = (csr.row_ptr[j]..csr.row_ptr[j+1])
                    .find(|&kk| csr.col_idx[kk] == i)
                    .map(|kk| csr.values[kk]).unwrap_or(0.0);
                assert!((kij - kji).abs() < 1e-12,
                    "Tet10 K[{},{}]={:.4e} != K[{},{}]={:.4e}", i, j, kij, j, i, kji);
            }
        }
    }

    /// P2 Tet10 trace(K) should equal 3 * eps * (surface integral via stiffness).
    /// For the reference tet (vol=1/6) and eps=1, trace(K) should be positive.
    #[test]
    fn tet10_stiffness_positive_diagonal() {
        let mesh = unit_tet10_mesh();
        let csr = assemble_stiffness(&mesh, |_| 1.0).unwrap().to_csr();
        let n = csr.nrows;
        for i in 0..n {
            let kii = (csr.row_ptr[i]..csr.row_ptr[i+1])
                .find(|&k| csr.col_idx[k] == i)
                .map(|k| csr.values[k]).unwrap_or(0.0);
            assert!(kii > 0.0, "Tet10 diagonal K[{},{}] = {:.4e} not positive", i, i, kii);
        }
    }

    /// P2 Tri6 trace: all diagonal entries must be positive.
    #[test]
    fn tri6_stiffness_positive_diagonal() {
        let mesh = unit_tri6_mesh();
        let csr = assemble_stiffness(&mesh, |_| 1.0).unwrap().to_csr();
        let n = csr.nrows;
        for i in 0..n {
            let kii = (csr.row_ptr[i]..csr.row_ptr[i+1])
                .find(|&k| csr.col_idx[k] == i)
                .map(|k| csr.values[k]).unwrap_or(0.0);
            assert!(kii > 0.0, "Tri6 diagonal K[{},{}] = {:.4e} not positive", i, i, kii);
        }
    }
}
