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
// Field probe sampling
// ---------------------------------------------------------------------------

/// Result of evaluating a field probe at a specific point.
#[derive(Debug, Clone)]
pub struct ProbeValue {
    /// Probe index (from config).
    pub index: u32,
    /// Probe coordinates [x, y, z].
    pub center: [f64; 3],
    /// Interpolated scalar potential φ.
    pub phi: f64,
    /// Interpolated E-field components (−∇φ).
    pub e_field: [f64; 3],
}

/// Evaluate the scalar field `phi` and its gradient at a list of probe points.
///
/// Each probe is described by `(index, center)` where `center` is `[x, y, z]`.
/// A brute-force element search is used (O(N_elements × N_probes)).
/// Returns one `ProbeValue` per probe; un-found probes get NaN values.
pub fn evaluate_probes(
    phi: &[f64],
    mesh: &RemMesh,
    probes: &[(u32, [f64; 3])],
) -> Vec<ProbeValue> {
    probes.iter().map(|&(index, center)| {
        let mut result = ProbeValue { index, center, phi: f64::NAN, e_field: [f64::NAN; 3] };

        'search: for elem in &mesh.volume_elements {
            match elem.kind {
                ElementKind::Tri3 => {
                    if let Some(bary) = bary_tri3(center, elem, mesh) {
                        result.phi = bary[0]*phi[elem.node_ids[0]]
                                   + bary[1]*phi[elem.node_ids[1]]
                                   + bary[2]*phi[elem.node_ids[2]];
                        if let Some(g) = tri3_grad(phi, elem, mesh) {
                            result.e_field = [-g[0], -g[1], -g[2]];
                        }
                        break 'search;
                    }
                }
                ElementKind::Tet4 => {
                    if let Some(bary) = bary_tet4(center, elem, mesh) {
                        result.phi = bary[0]*phi[elem.node_ids[0]]
                                   + bary[1]*phi[elem.node_ids[1]]
                                   + bary[2]*phi[elem.node_ids[2]]
                                   + bary[3]*phi[elem.node_ids[3]];
                        if let Some(g) = tet4_grad(phi, elem, mesh) {
                            result.e_field = [-g[0], -g[1], -g[2]];
                        }
                        break 'search;
                    }
                }
                ElementKind::Tri6 => {
                    if let Some(bary) = bary_tri3_corners(center, elem, mesh) {
                        // Use P1 interpolation on corner nodes for φ; centroid gradient
                        let ids = &elem.node_ids;
                        result.phi = bary[0]*phi[ids[0]] + bary[1]*phi[ids[1]] + bary[2]*phi[ids[2]];
                        if let Some(g) = tri6_grad_at_bary(phi, elem, mesh, bary) {
                            result.e_field = [-g[0], -g[1], -g[2]];
                        }
                        break 'search;
                    }
                }
                ElementKind::Tet10 => {
                    if let Some(bary) = bary_tet4_corners(center, elem, mesh) {
                        // Use P1 interpolation on corner nodes for φ; centroid gradient
                        let ids = &elem.node_ids;
                        result.phi = bary[0]*phi[ids[0]] + bary[1]*phi[ids[1]]
                                   + bary[2]*phi[ids[2]] + bary[3]*phi[ids[3]];
                        if let Some(g) = tet10_grad_at_bary(phi, elem, mesh, bary) {
                            result.e_field = [-g[0], -g[1], -g[2]];
                        }
                        break 'search;
                    }
                }
                _ => {}
            }
        }

        if result.phi.is_nan() {
            log::warn!("Probe {} at ({:.4}, {:.4}, {:.4}): no containing element found \
                        (point outside mesh or unsupported element type).",
                index, center[0], center[1], center[2]);
        }
        result
    }).collect()
}

// ---------------------------------------------------------------------------
// Probe output writers
// ---------------------------------------------------------------------------

/// Write probe potential values to `<output_dir>/postpro/probe-phi.csv`.
pub fn write_probe_phi_csv(
    output_dir: &std::path::Path,
    probes: &[ProbeValue],
) -> std::io::Result<()> {
    use std::io::Write;
    let dir = output_dir.join("postpro");
    std::fs::create_dir_all(&dir)?;
    let path = dir.join("probe-phi.csv");
    let mut f = std::fs::File::create(&path)?;
    writeln!(f, r#""Probe Index","x (m)","y (m)","z (m)","Phi (V)""#)?;
    for p in probes {
        writeln!(f, "{},{:.9e},{:.9e},{:.9e},{:.9e}",
            p.index, p.center[0], p.center[1], p.center[2], p.phi)?;
    }
    log::info!("Written: {}", path.display());
    Ok(())
}

/// Write probe E-field values to `<output_dir>/postpro/probe-E.csv`.
pub fn write_probe_e_csv(
    output_dir: &std::path::Path,
    probes: &[ProbeValue],
) -> std::io::Result<()> {
    use std::io::Write;
    let dir = output_dir.join("postpro");
    std::fs::create_dir_all(&dir)?;
    let path = dir.join("probe-E.csv");
    let mut f = std::fs::File::create(&path)?;
    writeln!(f, r#""Probe Index","x (m)","y (m)","z (m)","Ex (V/m)","Ey (V/m)","Ez (V/m)","||E|| (V/m)""#)?;
    for p in probes {
        let enorm = (p.e_field[0].powi(2) + p.e_field[1].powi(2) + p.e_field[2].powi(2)).sqrt();
        writeln!(f, "{},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e}",
            p.index, p.center[0], p.center[1], p.center[2],
            p.e_field[0], p.e_field[1], p.e_field[2], enorm)?;
    }
    log::info!("Written: {}", path.display());
    Ok(())
}

/// Write multi-mode probe potential values to `<output_dir>/postpro/probe-phi-modes.csv`.
///
/// Each row holds one probe sample from one mode: `mode_index, probe_index, x, y, z, phi`.
/// `mode_probes` is a slice of `(mode_index, probe_values)` pairs (1-based mode_index).
pub fn write_probe_modal_csv(
    output_dir: &std::path::Path,
    mode_probes: &[(usize, Vec<ProbeValue>)],
) -> std::io::Result<()> {
    use std::io::Write;
    let dir = output_dir.join("postpro");
    std::fs::create_dir_all(&dir)?;
    let path = dir.join("probe-phi-modes.csv");
    let mut f = std::fs::File::create(&path)?;
    writeln!(f, r#""Mode","Probe Index","x (m)","y (m)","z (m)","Phi""#)?;
    for (mode_idx, probes) in mode_probes {
        for p in probes {
            writeln!(f, "{},{},{:.9e},{:.9e},{:.9e},{:.9e}",
                mode_idx, p.index, p.center[0], p.center[1], p.center[2], p.phi)?;
        }
    }
    log::info!("Written: {}", path.display());
    Ok(())
}

// ---------------------------------------------------------------------------
// Barycentric coordinate helpers
// ---------------------------------------------------------------------------

const BARY_TOL: f64 = -1e-10;

/// Barycentric coordinates of `pt` in triangle (corner nodes of elem).
/// Returns Some([λ0, λ1, λ2]) if inside (λi ≥ BARY_TOL and Σλi ≈ 1).
fn bary_tri3(pt: [f64; 3], elem: &rem_mesh::Element, mesh: &RemMesh) -> Option<[f64; 3]> {
    let [n0, n1, n2] = [elem.node_ids[0], elem.node_ids[1], elem.node_ids[2]];
    bary_tri3_nodes(pt, &mesh.nodes[n0], &mesh.nodes[n1], &mesh.nodes[n2])
}

/// Same but uses only first 3 node IDs (for Tri6 corner-only containment test).
fn bary_tri3_corners(pt: [f64; 3], elem: &rem_mesh::Element, mesh: &RemMesh) -> Option<[f64; 3]> {
    let [n0, n1, n2] = [elem.node_ids[0], elem.node_ids[1], elem.node_ids[2]];
    bary_tri3_nodes(pt, &mesh.nodes[n0], &mesh.nodes[n1], &mesh.nodes[n2])
}

fn bary_tri3_nodes(pt: [f64; 3], n0: &rem_mesh::Node, n1: &rem_mesh::Node, n2: &rem_mesh::Node) -> Option<[f64; 3]> {
    let (x0, y0) = (n0.x, n0.y);
    let (x1, y1) = (n1.x, n1.y);
    let (x2, y2) = (n2.x, n2.y);
    let t = (y1-y2)*(x0-x2) + (x2-x1)*(y0-y2);
    if t.abs() < 1e-300 { return None; }
    let l0 = ((y1-y2)*(pt[0]-x2) + (x2-x1)*(pt[1]-y2)) / t;
    let l1 = ((y2-y0)*(pt[0]-x2) + (x0-x2)*(pt[1]-y2)) / t;
    let l2 = 1.0 - l0 - l1;
    if l0 >= BARY_TOL && l1 >= BARY_TOL && l2 >= BARY_TOL {
        Some([l0, l1, l2])
    } else {
        None
    }
}

/// Barycentric coordinates of `pt` in Tet4.
fn bary_tet4(pt: [f64; 3], elem: &rem_mesh::Element, mesh: &RemMesh) -> Option<[f64; 4]> {
    let ids = &elem.node_ids;
    bary_tet4_nodes(pt, &mesh.nodes[ids[0]], &mesh.nodes[ids[1]], &mesh.nodes[ids[2]], &mesh.nodes[ids[3]])
}

/// Same but uses only first 4 corner node IDs (for Tet10 containment test).
fn bary_tet4_corners(pt: [f64; 3], elem: &rem_mesh::Element, mesh: &RemMesh) -> Option<[f64; 4]> {
    let ids = &elem.node_ids;
    bary_tet4_nodes(pt, &mesh.nodes[ids[0]], &mesh.nodes[ids[1]], &mesh.nodes[ids[2]], &mesh.nodes[ids[3]])
}

fn bary_tet4_nodes(
    pt: [f64; 3],
    n0: &rem_mesh::Node, n1: &rem_mesh::Node,
    n2: &rem_mesh::Node, n3: &rem_mesh::Node,
) -> Option<[f64; 4]> {
    let (x0,y0,z0) = (n0.x, n0.y, n0.z);
    let (x1,y1,z1) = (n1.x, n1.y, n1.z);
    let (x2,y2,z2) = (n2.x, n2.y, n2.z);
    let (x3,y3,z3) = (n3.x, n3.y, n3.z);
    let j = [
        [x1-x0, x2-x0, x3-x0],
        [y1-y0, y2-y0, y3-y0],
        [z1-z0, z2-z0, z3-z0],
    ];
    let det = crate::assemble::det3_pub(&j);
    if det.abs() < 1e-300 { return None; }
    let j_inv = crate::assemble::inv3_pub(&j, det);
    let dp = [pt[0]-x0, pt[1]-y0, pt[2]-z0];
    let l1 = j_inv[0][0]*dp[0] + j_inv[0][1]*dp[1] + j_inv[0][2]*dp[2];
    let l2 = j_inv[1][0]*dp[0] + j_inv[1][1]*dp[1] + j_inv[1][2]*dp[2];
    let l3 = j_inv[2][0]*dp[0] + j_inv[2][1]*dp[1] + j_inv[2][2]*dp[2];
    let l0 = 1.0 - l1 - l2 - l3;
    if l0 >= BARY_TOL && l1 >= BARY_TOL && l2 >= BARY_TOL && l3 >= BARY_TOL {
        Some([l0, l1, l2, l3])
    } else {
        None
    }
}

// ---------------------------------------------------------------------------
// Gradient at barycentric coordinates for P2 elements
// ---------------------------------------------------------------------------

/// Gradient of φ at reference point (ξ=l1, η=l2) in a Tri6 element.
fn tri6_grad_at_bary(phi: &[f64], elem: &rem_mesh::Element, mesh: &RemMesh, bary: [f64; 3]) -> Option<[f64; 3]> {
    let ids = &elem.node_ids;
    let xy: [[f64; 2]; 6] = {
        let mut c = [[0.0f64; 2]; 6];
        for (i, &n) in ids.iter().enumerate() { c[i] = [mesh.nodes[n].x, mesh.nodes[n].y]; }
        c
    };
    let xi = bary[1]; let eta = bary[2];
    let l1 = 1.0 - xi - eta; let l2 = xi; let l3 = eta;
    let dndxi  = [-(4.0*l1-1.0), 4.0*l2-1.0, 0.0, 4.0*(l1-l2), 4.0*l3, -4.0*l3];
    let dndeta = [-(4.0*l1-1.0), 0.0, 4.0*l3-1.0, -4.0*l2, 4.0*l2, 4.0*(l1-l3)];
    let mut jac = [[0.0f64; 2]; 2];
    for i in 0..6 {
        jac[0][0] += dndxi[i]*xy[i][0]; jac[0][1] += dndeta[i]*xy[i][0];
        jac[1][0] += dndxi[i]*xy[i][1]; jac[1][1] += dndeta[i]*xy[i][1];
    }
    let det_j = jac[0][0]*jac[1][1] - jac[0][1]*jac[1][0];
    if det_j.abs() < 1e-300 { return None; }
    let inv_det = 1.0 / det_j;
    let ji = [[ jac[1][1]*inv_det, -jac[0][1]*inv_det], [-jac[1][0]*inv_det, jac[0][0]*inv_det]];
    let mut dphidxi = 0.0; let mut dphideta = 0.0;
    for i in 0..6 { dphidxi += dndxi[i]*phi[ids[i]]; dphideta += dndeta[i]*phi[ids[i]]; }
    Some([ji[0][0]*dphidxi + ji[1][0]*dphideta, ji[0][1]*dphidxi + ji[1][1]*dphideta, 0.0])
}

/// Gradient of φ at reference point (ξ=l1, η=l2, ζ=l3) in a Tet10 element.
fn tet10_grad_at_bary(phi: &[f64], elem: &rem_mesh::Element, mesh: &RemMesh, bary: [f64; 4]) -> Option<[f64; 3]> {
    let ids = &elem.node_ids;
    let xyz: [[f64; 3]; 10] = {
        let mut c = [[0.0f64; 3]; 10];
        for (i, &n) in ids.iter().enumerate() { c[i] = [mesh.nodes[n].x, mesh.nodes[n].y, mesh.nodes[n].z]; }
        c
    };
    let xi = bary[1]; let eta = bary[2]; let zet = bary[3];
    let l1 = 1.0 - xi - eta - zet;
    let dndxi = [-(4.0*l1-1.0), 4.0*xi-1.0, 0.0, 0.0, 4.0*(l1-xi), 4.0*eta, -4.0*eta, -4.0*zet, 4.0*zet, 0.0];
    let dndeta = [-(4.0*l1-1.0), 0.0, 4.0*eta-1.0, 0.0, -4.0*xi, 4.0*xi, 4.0*(l1-eta), -4.0*zet, 0.0, 4.0*zet];
    let dndzet = [-(4.0*l1-1.0), 0.0, 0.0, 4.0*zet-1.0, -4.0*xi, 0.0, -4.0*eta, 4.0*(l1-zet), 4.0*xi, 4.0*eta];
    let mut jac = [[0.0f64; 3]; 3];
    for i in 0..10 {
        for k in 0..3 {
            jac[k][0] += dndxi[i]*xyz[i][k];
            jac[k][1] += dndeta[i]*xyz[i][k];
            jac[k][2] += dndzet[i]*xyz[i][k];
        }
    }
    let det_j = crate::assemble::det3_pub(&jac);
    if det_j.abs() < 1e-300 { return None; }
    let j_inv = crate::assemble::inv3_pub(&jac, det_j);
    let mut ref_grad = [0.0f64; 3];
    for i in 0..10 { ref_grad[0] += dndxi[i]*phi[ids[i]]; ref_grad[1] += dndeta[i]*phi[ids[i]]; ref_grad[2] += dndzet[i]*phi[ids[i]]; }
    let mut grad = [0.0f64; 3];
    for row in 0..3 { for col in 0..3 { grad[row] += j_inv[col][row]*ref_grad[col]; } }
    Some(grad)
}

// ---------------------------------------------------------------------------
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
