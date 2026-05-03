//! P1/P2 finite element mass matrix assembly.
//!
//! M_ij = Σ_e ε_e ∫_Ωe φ_i φ_j dΩ
//!
//! Supports: Tri3, Tri6, Tet4, Tet10, Hex8.

use rem_core::{TripletMatrix, RemError, RemResult};
use rem_mesh::{RemMesh, ElementKind};

/// Assemble the global mass matrix with a scalar coefficient per element.
///
/// `coeff_fn` maps physical group tag → permittivity ε (or similar coefficient).
pub fn assemble_mass(
    mesh: &RemMesh,
    coeff_fn: impl Fn(u32) -> f64,
) -> RemResult<TripletMatrix> {
    let n = mesh.n_nodes();
    let cap = mesh.n_volume_elements() * 16;
    let mut triplet = TripletMatrix::with_capacity(n, n, cap);

    for elem in &mesh.volume_elements {
        if mesh.size > 1 && elem.rank != mesh.rank {
            continue;
        }
        let eps = coeff_fn(elem.tag);
        match elem.kind {
            ElementKind::Tri3 => mass_tri3(mesh, elem, eps, &mut triplet)?,
            ElementKind::Tri6 => mass_tri6(mesh, elem, eps, &mut triplet)?,
            ElementKind::Tet4 => mass_tet4(mesh, elem, eps, &mut triplet)?,
            ElementKind::Tet10 => mass_tet10(mesh, elem, eps, &mut triplet)?,
            ElementKind::Hex8 => mass_hex8(mesh, elem, eps, &mut triplet)?,
            other => {
                log::warn!("Mass matrix: element {:?} not supported — skipping", other);
            }
        }
    }

    Ok(triplet)
}

/// Assemble the global mass matrix with a full anisotropic tensor coefficient.
///
/// The mass integral `M_ij = ∫ ε̄ φ_i φ_j dΩ` uses a scalar effective permittivity
/// derived from the tensor as `ε̄ = trace(A) / 3`.  For isotropic rotation this equals
/// ε₀εᵣ exactly; for genuinely anisotropic media it gives the isotropic equivalent.
///
/// `tensor_fn` maps physical group tag → absolute permittivity tensor [F/m], 3×3 row-major.
pub fn assemble_mass_aniso(
    mesh: &RemMesh,
    tensor_fn: impl Fn(u32) -> [[f64; 3]; 3],
) -> RemResult<TripletMatrix> {
    assemble_mass(mesh, |tag| {
        let a = tensor_fn(tag);
        (a[0][0] + a[1][1] + a[2][2]) / 3.0
    })
}

fn mass_tri3(
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

    let det_j = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
    let area = 0.5 * det_j.abs();
    if area < 1e-300 {
        return Err(RemError::Mesh(format!(
            "Degenerate Tri3 element {} in mass assembly", elem.id
        )));
    }

    // ∫_Tri φ_i φ_j dA = area * (1 + δ_{ij}) / 12   (exact for linear triangle)
    let nodes = [n0, n1, n2];
    let diag = eps * area / 6.0;
    let off  = eps * area / 12.0;

    for i in 0..3 {
        for j in 0..3 {
            let val = if i == j { diag } else { off };
            triplet.add(nodes[i], nodes[j], val);
        }
    }
    Ok(())
}

fn mass_tet4(
    mesh: &RemMesh,
    elem: &rem_mesh::Element,
    eps: f64,
    triplet: &mut TripletMatrix,
) -> RemResult<()> {
    mass_tet4_first4(mesh, &elem.node_ids, elem.id, eps, triplet)
}

/// Use first 4 (corner) nodes of a Tet10 element for P1 approximation.
fn mass_tet4_corners(
    mesh: &RemMesh,
    elem: &rem_mesh::Element,
    eps: f64,
    triplet: &mut TripletMatrix,
) -> RemResult<()> {
    mass_tet4_first4(mesh, &elem.node_ids[..4], elem.id, eps, triplet)
}

// ---------------------------------------------------------------------------
// P2 quadratic triangle (Tri6) mass assembly
// ---------------------------------------------------------------------------

/// Mass matrix for a quadratic triangle (Tri6) using P2 basis functions.
///
/// Uses 6-point Gauss rule on reference triangle (exact for degree 4).
/// Reference element: ξ∈[0,1], η∈[0,1-ξ].
fn mass_tri6(
    mesh: &RemMesh,
    elem: &rem_mesh::Element,
    eps: f64,
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

    // 6-point Gauss rule on triangle (exact for degree 4).
    // Derived from two sets of permutations of (a1,b1,b1) and (a2,b2,b2).
    let a1: f64 = 0.816847572980459;
    let b1: f64 = 0.091576213509771;
    let w1: f64 = 0.109951743655322 / 2.0;
    let a2: f64 = 0.108103018168070;
    let b2: f64 = 0.445948490915965;
    let w2: f64 = 0.223381589678011 / 2.0;
    // Points (ξ, η) and weights
    let gauss: [[f64; 3]; 6] = [
        [b1, b1, w1], [a1, b1, w1], [b1, a1, w1],
        [b2, b2, w2], [a2, b2, w2], [b2, a2, w2],
    ];

    let mut me = [[0.0f64; 6]; 6];

    for &[xi, eta, w] in &gauss {
        let l1 = 1.0 - xi - eta;
        let l2 = xi;
        let l3 = eta;

        // P2 shape function values
        let n_val = [
            l1*(2.0*l1-1.0), l2*(2.0*l2-1.0), l3*(2.0*l3-1.0),
            4.0*l1*l2, 4.0*l2*l3, 4.0*l1*l3,
        ];
        // P2 gradients in reference coords (∂φ/∂ξ, ∂φ/∂η)
        let dndxi = [
            -(4.0*l1-1.0), 4.0*l2-1.0, 0.0,
            4.0*(l1-l2), 4.0*l3, -4.0*l3,
        ];
        let dndeta = [
            -(4.0*l1-1.0), 0.0, 4.0*l3-1.0,
            -4.0*l2, 4.0*l2, 4.0*(l1-l3),
        ];

        // Jacobian
        let mut jac = [[0.0f64; 2]; 2];
        for i in 0..6 {
            jac[0][0] += dndxi[i]  * xy[i][0];
            jac[0][1] += dndeta[i] * xy[i][0];
            jac[1][0] += dndxi[i]  * xy[i][1];
            jac[1][1] += dndeta[i] * xy[i][1];
        }
        let det_j = jac[0][0]*jac[1][1] - jac[0][1]*jac[1][0];
        if det_j.abs() < 1e-300 { continue; }

        let wdet = eps * w * det_j.abs();
        for i in 0..6 {
            for j in 0..6 {
                me[i][j] += wdet * n_val[i] * n_val[j];
            }
        }
    }

    for i in 0..6 {
        for j in 0..6 {
            triplet.add(nids[i], nids[j], me[i][j]);
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// P2 quadratic tetrahedron (Tet10) mass assembly
// ---------------------------------------------------------------------------

/// Mass matrix for a quadratic tetrahedron (Tet10) using P2 basis functions.
///
/// Uses Keast 5-point rule on reference tet (exact for degree 3; sufficient in practice
/// for Tet10 meshes where spatial accuracy is O(h³)).
fn mass_tet10(
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

    // Keast 5-point rule (exact for degree 3):
    // p1 = centroid (1/4,1/4,1/4), w1 = -4/5 * (1/6) = -2/15
    // p2..p5 = (1/6,1/6,1/6), (1/2,1/6,1/6), (1/6,1/2,1/6), (1/6,1/6,1/2), w2..5 = 9/20*(1/6) = 3/40
    let gauss: [[f64; 4]; 5] = [
        [0.25,      0.25,      0.25,      -2.0/15.0],
        [1.0/6.0,   1.0/6.0,   1.0/6.0,   3.0/40.0 ],
        [1.0/2.0,   1.0/6.0,   1.0/6.0,   3.0/40.0 ],
        [1.0/6.0,   1.0/2.0,   1.0/6.0,   3.0/40.0 ],
        [1.0/6.0,   1.0/6.0,   1.0/2.0,   3.0/40.0 ],
    ];

    let mut me = [[0.0f64; 10]; 10];

    for &[xi, eta, zet, w] in &gauss {
        let l1 = 1.0 - xi - eta - zet;

        // P2 shape function values at (ξ,η,ζ)
        let n_val = [
            l1*(2.0*l1-1.0),
            xi*(2.0*xi-1.0),
            eta*(2.0*eta-1.0),
            zet*(2.0*zet-1.0),
            4.0*l1*xi,
            4.0*xi*eta,
            4.0*l1*eta,
            4.0*l1*zet,
            4.0*xi*zet,
            4.0*eta*zet,
        ];
        // Reference-space gradients for Jacobian computation
        let dndxi = [
            -(4.0*l1-1.0), 4.0*xi-1.0, 0.0, 0.0,
            4.0*(l1-xi), 4.0*eta, -4.0*eta, -4.0*zet, 4.0*zet, 0.0,
        ];
        let dndeta = [
            -(4.0*l1-1.0), 0.0, 4.0*eta-1.0, 0.0,
            -4.0*xi, 4.0*xi, 4.0*(l1-eta), -4.0*zet, 0.0, 4.0*zet,
        ];
        let dndzet = [
            -(4.0*l1-1.0), 0.0, 0.0, 4.0*zet-1.0,
            -4.0*xi, 0.0, -4.0*eta, 4.0*(l1-zet), 4.0*xi, 4.0*eta,
        ];

        // Jacobian
        let mut jac = [[0.0f64; 3]; 3];
        for i in 0..10 {
            let dref = [dndxi[i], dndeta[i], dndzet[i]];
            for k in 0..3 {
                jac[k][0] += dref[0] * xyz[i][k];
                jac[k][1] += dref[1] * xyz[i][k];
                jac[k][2] += dref[2] * xyz[i][k];
            }
        }
        let det_j = jac[0][0]*(jac[1][1]*jac[2][2]-jac[1][2]*jac[2][1])
                   -jac[0][1]*(jac[1][0]*jac[2][2]-jac[1][2]*jac[2][0])
                   +jac[0][2]*(jac[1][0]*jac[2][1]-jac[1][1]*jac[2][0]);
        if det_j.abs() < 1e-300 { continue; }

        let wdet = eps * w * det_j.abs();
        for i in 0..10 {
            for j in 0..10 {
                me[i][j] += wdet * n_val[i] * n_val[j];
            }
        }
    }

    for i in 0..10 {
        for j in 0..10 {
            triplet.add(nids[i], nids[j], me[i][j]);
        }
    }
    Ok(())
}
fn mass_tet4_first4(
    mesh: &RemMesh,
    node_ids: &[usize],
    elem_id: usize,
    eps: f64,
    triplet: &mut TripletMatrix,
) -> RemResult<()> {
    let [n0, n1, n2, n3] = [node_ids[0], node_ids[1], node_ids[2], node_ids[3]];
    let nodes = [n0, n1, n2, n3];

    let x = [mesh.nodes[n0].x, mesh.nodes[n1].x, mesh.nodes[n2].x, mesh.nodes[n3].x];
    let y = [mesh.nodes[n0].y, mesh.nodes[n1].y, mesh.nodes[n2].y, mesh.nodes[n3].y];
    let z = [mesh.nodes[n0].z, mesh.nodes[n1].z, mesh.nodes[n2].z, mesh.nodes[n3].z];

    // det(J) = 6 * vol
    let j = [
        [x[1]-x[0], x[2]-x[0], x[3]-x[0]],
        [y[1]-y[0], y[2]-y[0], y[3]-y[0]],
        [z[1]-z[0], z[2]-z[0], z[3]-z[0]],
    ];
    let det = j[0][0]*(j[1][1]*j[2][2]-j[1][2]*j[2][1])
            - j[0][1]*(j[1][0]*j[2][2]-j[1][2]*j[2][0])
            + j[0][2]*(j[1][0]*j[2][1]-j[1][1]*j[2][0]);
    let vol = det.abs() / 6.0;
    if vol < 1e-300 {
        return Err(RemError::Mesh(format!(
            "Degenerate Tet4 element {} in mass assembly", elem_id
        )));
    }

    // ∫_Tet φ_i φ_j dV = vol * (1 + δ_{ij}) / 20   (exact for linear tetrahedron)
    let diag = eps * vol / 10.0;
    let off  = eps * vol / 20.0;

    for i in 0..4 {
        for j in 0..4 {
            let val = if i == j { diag } else { off };
            triplet.add(nodes[i], nodes[j], val);
        }
    }
    Ok(())
}

/// Mass matrix for a trilinear hexahedron (Hex8) using 2×2×2 Gauss quadrature.
fn mass_hex8(
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

    let gp = 1.0_f64 / 3.0_f64.sqrt();
    let gauss_pts = [-gp, gp];
    let mut me = [[0.0f64; 8]; 8];

    for &xi in &gauss_pts {
        for &eta in &gauss_pts {
            for &zet in &gauss_pts {
                let mut n_val   = [0.0f64; 8];
                let mut dn_dxi  = [0.0f64; 8];
                let mut dn_deta = [0.0f64; 8];
                let mut dn_dzet = [0.0f64; 8];
                for i in 0..8 {
                    let (a, b, c) = (xi_ref[i], eta_ref[i], zet_ref[i]);
                    n_val[i]   = 0.125 * (1.0 + a*xi) * (1.0 + b*eta) * (1.0 + c*zet);
                    dn_dxi[i]  = 0.125 * a * (1.0 + b*eta) * (1.0 + c*zet);
                    dn_deta[i] = 0.125 * b * (1.0 + a*xi)  * (1.0 + c*zet);
                    dn_dzet[i] = 0.125 * c * (1.0 + a*xi)  * (1.0 + b*eta);
                }
                // Jacobian determinant
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
                let det_j = jac[0][0]*(jac[1][1]*jac[2][2]-jac[1][2]*jac[2][1])
                           -jac[0][1]*(jac[1][0]*jac[2][2]-jac[1][2]*jac[2][0])
                           +jac[0][2]*(jac[1][0]*jac[2][1]-jac[1][1]*jac[2][0]);
                if det_j.abs() < 1e-300 { continue; }

                let wdet = eps * det_j.abs();
                for i in 0..8 {
                    for j in 0..8 {
                        me[i][j] += wdet * n_val[i] * n_val[j];
                    }
                }
            }
        }
    }

    for i in 0..8 {
        for j in 0..8 {
            triplet.add(nids[i], nids[j], me[i][j]);
        }
    }
    Ok(())
}
