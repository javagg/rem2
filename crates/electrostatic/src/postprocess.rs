/// Post-processing: gradient (E-field) recovery and capacitance extraction.

use rem_mesh::{RemMesh, ElementKind};

// ---------------------------------------------------------------------------
// Gradient recovery (E = -∇φ)
// ---------------------------------------------------------------------------

/// Compute the electrostatic field **E = −∇φ** by nodal averaging of
/// element-constant gradients (Zienkiewicz–Zhu style recovery).
///
/// For P1 elements (Tri3, Tet4): gradient is constant per element.
/// For P2 elements (Tri6, Tet10): gradient is linear per element; we evaluate it
/// at the element centroid (a good representative sample for ZZ recovery).
///
/// Returns a `Vec<[f64; 3]>` — one 3-component vector per node.
pub fn gradient_recovery(phi: &[f64], mesh: &RemMesh) -> Vec<[f64; 3]> {
    let n = mesh.n_nodes();
    let mut e_sum = vec![[0.0f64; 3]; n];
    let mut count = vec![0usize; n];

    for elem in &mesh.volume_elements {
        let grad_opt = match elem.kind {
            ElementKind::Tri3 => tri3_grad(phi, elem, mesh),
            ElementKind::Tet4 => tet4_grad(phi, elem, mesh),
            ElementKind::Tri6 => tri6_grad_centroid(phi, elem, mesh),
            ElementKind::Tet10 => tet10_grad_centroid(phi, elem, mesh),
            _ => None,
        };
        if let Some(grad) = grad_opt {
            for &nid in &elem.node_ids {
                e_sum[nid][0] -= grad[0]; // E = -∇φ
                e_sum[nid][1] -= grad[1];
                e_sum[nid][2] -= grad[2];
                count[nid] += 1;
            }
        }
    }

    // Nodal average
    let mut e_field = vec![[0.0f64; 3]; n];
    for i in 0..n {
        if count[i] > 0 {
            let c = count[i] as f64;
            e_field[i] = [e_sum[i][0] / c, e_sum[i][1] / c, e_sum[i][2] / c];
        }
    }
    e_field
}

/// Element-constant gradient for a Tri3 element.
fn tri3_grad(phi: &[f64], elem: &rem_mesh::Element, mesh: &RemMesh) -> Option<[f64; 3]> {
    let [n0, n1, n2] = [elem.node_ids[0], elem.node_ids[1], elem.node_ids[2]];
    let (x0, y0) = (mesh.nodes[n0].x, mesh.nodes[n0].y);
    let (x1, y1) = (mesh.nodes[n1].x, mesh.nodes[n1].y);
    let (x2, y2) = (mesh.nodes[n2].x, mesh.nodes[n2].y);
    let det_j = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
    let area = 0.5 * det_j.abs();
    if area < 1e-300 { return None; }
    let inv2a = 1.0 / (2.0 * area);

    let grads = [
        [(y1 - y2) * inv2a, (x2 - x1) * inv2a],
        [(y2 - y0) * inv2a, (x0 - x2) * inv2a],
        [(y0 - y1) * inv2a, (x1 - x0) * inv2a],
    ];
    let phis = [phi[n0], phi[n1], phi[n2]];

    let mut grad = [0.0f64; 3];
    for i in 0..3 {
        grad[0] += grads[i][0] * phis[i];
        grad[1] += grads[i][1] * phis[i];
    }
    Some(grad)
}

/// Element-constant gradient for a Tet4 element.
fn tet4_grad(phi: &[f64], elem: &rem_mesh::Element, mesh: &RemMesh) -> Option<[f64; 3]> {
    let ids = &elem.node_ids;
    let (x0,y0,z0) = (mesh.nodes[ids[0]].x, mesh.nodes[ids[0]].y, mesh.nodes[ids[0]].z);
    let (x1,y1,z1) = (mesh.nodes[ids[1]].x, mesh.nodes[ids[1]].y, mesh.nodes[ids[1]].z);
    let (x2,y2,z2) = (mesh.nodes[ids[2]].x, mesh.nodes[ids[2]].y, mesh.nodes[ids[2]].z);
    let (x3,y3,z3) = (mesh.nodes[ids[3]].x, mesh.nodes[ids[3]].y, mesh.nodes[ids[3]].z);

    let j = [
        [x1-x0, x2-x0, x3-x0],
        [y1-y0, y2-y0, y3-y0],
        [z1-z0, z2-z0, z3-z0],
    ];
    let det = crate::assemble::det3_pub(&j);
    if det.abs() < 1e-300 { return None; }
    let j_inv = crate::assemble::inv3_pub(&j, det);

    let ref_grads: [[f64;3]; 4] = [
        [-1.0,-1.0,-1.0], [1.0,0.0,0.0], [0.0,1.0,0.0], [0.0,0.0,1.0],
    ];
    let mut phys_grads = [[0.0f64;3]; 4];
    for i in 0..4 {
        for row in 0..3 {
            for col in 0..3 {
                phys_grads[i][row] += j_inv[col][row] * ref_grads[i][col];
            }
        }
    }

    let mut grad = [0.0f64; 3];
    for i in 0..4 {
        let pi = phi[ids[i]];
        grad[0] += phys_grads[i][0] * pi;
        grad[1] += phys_grads[i][1] * pi;
        grad[2] += phys_grads[i][2] * pi;
    }
    Some(grad)
}

/// Element centroid gradient for a Tri6 element (P2 quadratic triangle).
///
/// Evaluates ∇φ at the centroid (ξ=η=1/3) using P2 shape function gradients.
fn tri6_grad_centroid(phi: &[f64], elem: &rem_mesh::Element, mesh: &RemMesh) -> Option<[f64; 3]> {
    debug_assert_eq!(elem.node_ids.len(), 6);
    let ids = &elem.node_ids;
    let xy: [[f64; 2]; 6] = {
        let mut c = [[0.0f64; 2]; 6];
        for (i, &n) in ids.iter().enumerate() {
            c[i] = [mesh.nodes[n].x, mesh.nodes[n].y];
        }
        c
    };

    // Centroid of reference triangle: ξ=η=1/3, λ1=λ2=λ3=1/3
    let xi = 1.0 / 3.0;
    let eta = 1.0 / 3.0;
    let l1 = 1.0 - xi - eta;
    let l2 = xi;
    let l3 = eta;

    let dndxi  = [-(4.0*l1-1.0), 4.0*l2-1.0, 0.0, 4.0*(l1-l2), 4.0*l3, -4.0*l3];
    let dndeta = [-(4.0*l1-1.0), 0.0, 4.0*l3-1.0, -4.0*l2, 4.0*l2, 4.0*(l1-l3)];

    let mut jac = [[0.0f64; 2]; 2];
    for i in 0..6 {
        jac[0][0] += dndxi[i]  * xy[i][0];
        jac[0][1] += dndeta[i] * xy[i][0];
        jac[1][0] += dndxi[i]  * xy[i][1];
        jac[1][1] += dndeta[i] * xy[i][1];
    }
    let det_j = jac[0][0]*jac[1][1] - jac[0][1]*jac[1][0];
    if det_j.abs() < 1e-300 { return None; }
    let inv_det = 1.0 / det_j;
    let ji = [
        [ jac[1][1]*inv_det, -jac[0][1]*inv_det],
        [-jac[1][0]*inv_det,  jac[0][0]*inv_det],
    ];

    // Physical gradient of φ: ∇φ = Σ φ_i ∇N_i
    let mut dphidxi = 0.0;
    let mut dphideta = 0.0;
    for i in 0..6 {
        dphidxi  += dndxi[i]  * phi[ids[i]];
        dphideta += dndeta[i] * phi[ids[i]];
    }
    Some([
        ji[0][0] * dphidxi + ji[1][0] * dphideta,
        ji[0][1] * dphidxi + ji[1][1] * dphideta,
        0.0,
    ])
}

/// Element centroid gradient for a Tet10 element (P2 quadratic tetrahedron).
///
/// Evaluates ∇φ at the centroid (ξ=η=ζ=1/4) using P2 shape function gradients.
fn tet10_grad_centroid(phi: &[f64], elem: &rem_mesh::Element, mesh: &RemMesh) -> Option<[f64; 3]> {
    debug_assert_eq!(elem.node_ids.len(), 10);
    let ids = &elem.node_ids;
    let xyz: [[f64; 3]; 10] = {
        let mut c = [[0.0f64; 3]; 10];
        for (i, &n) in ids.iter().enumerate() {
            c[i] = [mesh.nodes[n].x, mesh.nodes[n].y, mesh.nodes[n].z];
        }
        c
    };

    // Centroid: ξ=η=ζ=1/4, λ1=1/4
    let xi = 0.25;
    let eta = 0.25;
    let zet = 0.25;
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
            jac[k][0] += dref[0] * xyz[i][k];
            jac[k][1] += dref[1] * xyz[i][k];
            jac[k][2] += dref[2] * xyz[i][k];
        }
    }
    let det_j = crate::assemble::det3_pub(&jac);
    if det_j.abs() < 1e-300 { return None; }
    let j_inv = crate::assemble::inv3_pub(&jac, det_j);

    // Physical gradient of φ: ∇φ = Σ φ_i J^{-T} ∇_ξ N_i
    let mut ref_grad = [0.0f64; 3];
    for i in 0..10 {
        ref_grad[0] += dndxi[i]  * phi[ids[i]];
        ref_grad[1] += dndeta[i] * phi[ids[i]];
        ref_grad[2] += dndzet[i] * phi[ids[i]];
    }
    let mut grad = [0.0f64; 3];
    for row in 0..3 {
        for col in 0..3 {
            grad[row] += j_inv[col][row] * ref_grad[col];
        }
    }
    Some(grad)
}

// ---------------------------------------------------------------------------
// Electrostatic energy
// ---------------------------------------------------------------------------

/// Total electrostatic energy: U = (1/2) ∫ ε |∇φ|² dΩ
///
/// For P1 elements (Tri3, Tet4): gradient is constant → U_e = ½ ε |∇φ|² Vol_e.
/// For P2 elements (Tri6, Tet10): gradient is linear → Gauss quadrature is used.
pub fn electrostatic_energy(
    phi: &[f64],
    mesh: &RemMesh,
    coeff_fn: impl Fn(u32) -> f64,
) -> f64 {
    let mut energy = 0.0;
    for elem in &mesh.volume_elements {
        let eps = coeff_fn(elem.tag);
        match elem.kind {
            ElementKind::Tri3 => {
                if let Some(grad) = tri3_grad(phi, elem, mesh) {
                    let n0 = elem.node_ids[0];
                    let n1 = elem.node_ids[1];
                    let n2 = elem.node_ids[2];
                    let area = tri3_area(
                        mesh.nodes[n0].x, mesh.nodes[n0].y,
                        mesh.nodes[n1].x, mesh.nodes[n1].y,
                        mesh.nodes[n2].x, mesh.nodes[n2].y,
                    );
                    let grad_sq = grad[0]*grad[0] + grad[1]*grad[1] + grad[2]*grad[2];
                    energy += 0.5 * eps * area * grad_sq;
                }
            }
            ElementKind::Tet4 => {
                if let Some(grad) = tet4_grad(phi, elem, mesh) {
                    let ids = &elem.node_ids;
                    let vol = tet4_volume(mesh, ids);
                    let grad_sq = grad[0]*grad[0] + grad[1]*grad[1] + grad[2]*grad[2];
                    energy += 0.5 * eps * vol * grad_sq;
                }
            }
            ElementKind::Tri6 => {
                energy += tri6_energy(phi, elem, mesh, eps);
            }
            ElementKind::Tet10 => {
                energy += tet10_energy(phi, elem, mesh, eps);
            }
            _ => {}
        }
    }
    energy
}

fn tri3_area(x0: f64, y0: f64, x1: f64, y1: f64, x2: f64, y2: f64) -> f64 {
    0.5 * ((x1-x0)*(y2-y0) - (x2-x0)*(y1-y0)).abs()
}

fn tet4_volume(mesh: &RemMesh, ids: &[usize]) -> f64 {
    let (x0,y0,z0) = (mesh.nodes[ids[0]].x, mesh.nodes[ids[0]].y, mesh.nodes[ids[0]].z);
    let (x1,y1,z1) = (mesh.nodes[ids[1]].x, mesh.nodes[ids[1]].y, mesh.nodes[ids[1]].z);
    let (x2,y2,z2) = (mesh.nodes[ids[2]].x, mesh.nodes[ids[2]].y, mesh.nodes[ids[2]].z);
    let (x3,y3,z3) = (mesh.nodes[ids[3]].x, mesh.nodes[ids[3]].y, mesh.nodes[ids[3]].z);
    let j = [
        [x1-x0, x2-x0, x3-x0],
        [y1-y0, y2-y0, y3-y0],
        [z1-z0, z2-z0, z3-z0],
    ];
    crate::assemble::det3_pub(&j).abs() / 6.0
}

/// Energy contribution from a Tri6 element: ½ ε ∫ |∇φ|² dΩ
/// Uses 3-point Gauss rule (degree-2 exact, sufficient for |∇φ|² of P2 field).
fn tri6_energy(phi: &[f64], elem: &rem_mesh::Element, mesh: &RemMesh, eps: f64) -> f64 {
    const GP: [[f64; 2]; 3] = [
        [1.0/6.0, 1.0/6.0],
        [2.0/3.0, 1.0/6.0],
        [1.0/6.0, 2.0/3.0],
    ];
    const W: [f64; 3] = [1.0/3.0, 1.0/3.0, 1.0/3.0];

    let ids = &elem.node_ids;
    let xy: [[f64; 2]; 6] = {
        let mut c = [[0.0f64; 2]; 6];
        for (i, &n) in ids.iter().enumerate() {
            c[i] = [mesh.nodes[n].x, mesh.nodes[n].y];
        }
        c
    };

    let mut energy = 0.0;
    for (p, &w) in GP.iter().zip(W.iter()) {
        let xi = p[0]; let eta = p[1];
        let l1 = 1.0 - xi - eta;
        let l2 = xi; let l3 = eta;

        let dndxi  = [-(4.0*l1-1.0), 4.0*l2-1.0, 0.0, 4.0*(l1-l2), 4.0*l3, -4.0*l3];
        let dndeta = [-(4.0*l1-1.0), 0.0, 4.0*l3-1.0, -4.0*l2, 4.0*l2, 4.0*(l1-l3)];

        let mut jac = [[0.0f64; 2]; 2];
        for i in 0..6 {
            jac[0][0] += dndxi[i]  * xy[i][0];
            jac[0][1] += dndeta[i] * xy[i][0];
            jac[1][0] += dndxi[i]  * xy[i][1];
            jac[1][1] += dndeta[i] * xy[i][1];
        }
        let det_j = jac[0][0]*jac[1][1] - jac[0][1]*jac[1][0];
        if det_j.abs() < 1e-300 { continue; }
        let inv_det = 1.0 / det_j;
        let ji = [
            [ jac[1][1]*inv_det, -jac[0][1]*inv_det],
            [-jac[1][0]*inv_det,  jac[0][0]*inv_det],
        ];

        let mut dphidxi = 0.0; let mut dphideta = 0.0;
        for i in 0..6 {
            dphidxi  += dndxi[i]  * phi[ids[i]];
            dphideta += dndeta[i] * phi[ids[i]];
        }
        let gx = ji[0][0]*dphidxi + ji[1][0]*dphideta;
        let gy = ji[0][1]*dphidxi + ji[1][1]*dphideta;
        let grad_sq = gx*gx + gy*gy;
        // Reference area is 0.5; weight includes det_j*area_ref = 0.5*|det_j|
        energy += 0.5 * eps * w * grad_sq * det_j.abs() * 0.5;
    }
    energy
}

/// Energy contribution from a Tet10 element: ½ ε ∫ |∇φ|² dΩ
/// Uses Keast 4-point Gauss rule (degree-2 exact).
fn tet10_energy(phi: &[f64], elem: &rem_mesh::Element, mesh: &RemMesh, eps: f64) -> f64 {
    const A: f64 = 0.585_410_196_624_97;
    const B: f64 = 0.138_196_601_125_01;
    const W4: f64 = 1.0 / 4.0;
    let gp: [[f64; 3]; 4] = [
        [A, B, B],
        [B, A, B],
        [B, B, A],
        [B, B, B],
    ];

    let ids = &elem.node_ids;
    let xyz: [[f64; 3]; 10] = {
        let mut c = [[0.0f64; 3]; 10];
        for (i, &n) in ids.iter().enumerate() {
            c[i] = [mesh.nodes[n].x, mesh.nodes[n].y, mesh.nodes[n].z];
        }
        c
    };

    let mut energy = 0.0;
    for p in &gp {
        let xi = p[0]; let eta = p[1]; let zet = p[2];
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
                jac[k][0] += dref[0] * xyz[i][k];
                jac[k][1] += dref[1] * xyz[i][k];
                jac[k][2] += dref[2] * xyz[i][k];
            }
        }
        let det_j = crate::assemble::det3_pub(&jac);
        if det_j.abs() < 1e-300 { continue; }
        let j_inv = crate::assemble::inv3_pub(&jac, det_j);

        let mut ref_grad = [0.0f64; 3];
        for i in 0..10 {
            ref_grad[0] += dndxi[i]  * phi[ids[i]];
            ref_grad[1] += dndeta[i] * phi[ids[i]];
            ref_grad[2] += dndzet[i] * phi[ids[i]];
        }
        let mut grad = [0.0f64; 3];
        for row in 0..3 {
            for col in 0..3 {
                grad[row] += j_inv[col][row] * ref_grad[col];
            }
        }
        let grad_sq = grad[0]*grad[0] + grad[1]*grad[1] + grad[2]*grad[2];
        // Reference tet volume = 1/6; weight = W4/6 * |det_j|
        energy += 0.5 * eps * W4 * grad_sq * det_j.abs() / 6.0;
    }
    energy
}

// ---------------------------------------------------------------------------
// Capacitance extraction
// ---------------------------------------------------------------------------

/// Extract the capacitance between the excited electrode and ground.
///
/// Uses the energy method: **C = 2U / V²**
///
/// This is exact for a linear FEM solution with Dirichlet BC φ = `v_applied`
/// on the excited conductor and φ = 0 on ground.
pub fn capacitance(
    phi: &[f64],
    mesh: &RemMesh,
    coeff_fn: impl Fn(u32) -> f64,
    v_applied: f64,
) -> f64 {
    if v_applied.abs() < 1e-300 {
        return 0.0;
    }
    let energy = electrostatic_energy(phi, mesh, coeff_fn);
    2.0 * energy / (v_applied * v_applied)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assemble::tests::unit_triangle_mesh;

    #[test]
    fn gradient_linear_phi() {
        // φ(x, y) = x → ∇φ = (1, 0) → E = (-1, 0)
        let mesh = unit_triangle_mesh();
        let phi = vec![
            mesh.nodes[0].x,
            mesh.nodes[1].x,
            mesh.nodes[2].x,
        ];
        let e = gradient_recovery(&phi, &mesh);
        for ev in &e {
            assert!((ev[0] - (-1.0)).abs() < 1e-13, "E_x={}", ev[0]);
            assert!(ev[1].abs() < 1e-13, "E_y={}", ev[1]);
        }
    }

    #[test]
    fn gradient_linear_phi_y() {
        // φ(x, y) = y → E = (0, -1)
        let mesh = unit_triangle_mesh();
        let phi = vec![
            mesh.nodes[0].y,
            mesh.nodes[1].y,
            mesh.nodes[2].y,
        ];
        let e = gradient_recovery(&phi, &mesh);
        for ev in &e {
            assert!(ev[0].abs() < 1e-13);
            assert!((ev[1] - (-1.0)).abs() < 1e-13);
        }
    }
}
