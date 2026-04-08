//! Singular integral handlers for coincident and near-singular triangle pairs.
//!
//! ## Cases handled
//!
//! | Relationship      | Method                        | Used by         |
//! |-------------------|-------------------------------|-----------------|
//! | Identical (self)  | Duffy 2D polar transform      | EFIE pulse diag |
//! | Shared edge       | Sauter-Schwab 4-D rule (edge) | EFIE/MFIE off-diag |
//! | Shared vertex     | Sauter-Schwab 4-D rule (vtx)  | EFIE/MFIE off-diag |
//! | Well-separated    | Standard Gauss quadrature     | (caller handles) |
//!
//! Reference: Sauter & Schwab, *Boundary Element Methods* (2011), §5.2–5.3;
//!            Rao, Wilton, Glisson (1982), Appendix B.

use crate::surface_mesh::TriFace;
use crate::green::green3d;
use num_complex::Complex64;

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Classify the geometric relationship between two triangles.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TriPairType {
    /// Same triangle (m == n index)
    Identical,
    /// Share exactly one edge (2 common nodes)
    SharedEdge,
    /// Share exactly one vertex (1 common node)
    SharedVertex,
    /// No shared nodes
    Disjoint,
}

/// Determine the geometric relationship and return shared node local indices.
pub fn classify_pair(face_m: &TriFace, face_n: &TriFace) -> TriPairType {
    let mut common = 0u32;
    for &nm in &face_m.nodes {
        for &nn in &face_n.nodes {
            if nm == nn { common += 1; }
        }
    }
    match common {
        3 => TriPairType::Identical,
        2 => TriPairType::SharedEdge,
        1 => TriPairType::SharedVertex,
        _ => TriPairType::Disjoint,
    }
}

/// Compute ∫_Tm ∫_Tn G(r,r') dS' dS for a coincident or near-singular pair.
///
/// Multiplied by (-jωμ₀) externally (caller's responsibility).
/// Returns the raw double-surface integral of G (scalar).
pub fn zmn_singular_pulse(
    face_m: &TriFace,
    face_n: &TriFace,
    nodes: &[[f64; 3]],
    k: f64,
    n_gauss: usize,
) -> Complex64 {
    let pair = classify_pair(face_m, face_n);
    match pair {
        TriPairType::Identical => {
            duffy_self_integral(face_m, nodes, k, n_gauss)
        }
        TriPairType::SharedEdge => {
            sauter_schwab_shared_edge(face_m, face_n, nodes, k, n_gauss)
        }
        TriPairType::SharedVertex => {
            sauter_schwab_shared_vertex(face_m, face_n, nodes, k, n_gauss)
        }
        TriPairType::Disjoint => {
            panic!("zmn_singular_pulse called on disjoint pair — use regular quadrature");
        }
    }
}

// ---------------------------------------------------------------------------
// RWG singular integrals
// ---------------------------------------------------------------------------

/// Compute the EFIE integrand ⟨f_m, L(f_n)⟩ using Sauter-Schwab/Duffy quadrature
/// for a near-singular pair of RWG half-support faces.
///
/// The EFIE vector potential kernel is:
///   A_term = ∫∫ G(r,r') f_m(r) · f_n(r') dS' dS
///   Phi_term = ∫∫ G(r,r') divS_fm · divS_fn dS' dS  (scalar)
///
/// Both are integrated with singular-safe quadrature when r and r' are on
/// coincident, shared-edge, or shared-vertex triangles.
///
/// Returns (a_term, phi_term) as Complex64 values.
/// The caller assembles: Z[m,n] += -jωμ (a_term - inv_omega_eps/omega_mu * phi_term)
pub fn zmn_efie_rwg_singular(
    face_m: &TriFace,
    face_n: &TriFace,
    fm_fn: &dyn Fn(&[f64; 3], &[f64; 3]) -> (f64, f64),  // (f_m·f_n, div_m*div_n) at (r,r')
    nodes: &[[f64; 3]],
    k: f64,
    n_gauss: usize,
) -> (Complex64, Complex64) {
    let pair = classify_pair(face_m, face_n);
    match pair {
        TriPairType::Identical => {
            rwg_efie_duffy_self(face_m, fm_fn, nodes, k, n_gauss)
        }
        TriPairType::SharedEdge => {
            rwg_efie_sauter_schwab_edge(face_m, face_n, fm_fn, nodes, k, n_gauss)
        }
        TriPairType::SharedVertex => {
            rwg_efie_sauter_schwab_vertex(face_m, face_n, fm_fn, nodes, k, n_gauss)
        }
        TriPairType::Disjoint => {
            panic!("zmn_efie_rwg_singular called on disjoint pair");
        }
    }
}

/// Duffy self-integral for RWG EFIE (identical faces).
fn rwg_efie_duffy_self(
    face: &TriFace,
    fm_fn: &dyn Fn(&[f64; 3], &[f64; 3]) -> (f64, f64),
    nodes: &[[f64; 3]],
    k: f64,
    n_gauss: usize,
) -> (Complex64, Complex64) {
    let (gl_pts, gl_wts) = gauss_legendre_1d(n_gauss);
    let [i0, i1, i2] = face.nodes;
    let v = [nodes[i0], nodes[i1], nodes[i2]];
    let area = face.area;

    let mut a_sum = Complex64::ZERO;
    let mut phi_sum = Complex64::ZERO;

    for &pivot in &[0usize, 1, 2] {
        let va = &v[pivot];
        let vb = &v[(pivot + 1) % 3];
        let vc = &v[(pivot + 2) % 3];
        let area_sub = sub_triangle_area(va, vb, vc);

        // Outer loop: observation point r (standard Gauss on face)
        for (bm_pt_idx, &wm) in gl_pts.iter().enumerate() {
            // Use face-level Gauss for outer (observation)
            for (bm_pt_jdx, &wm2) in gl_pts.iter().enumerate() {
                let bary_m = [(1.0 - wm) * (1.0 - wm2), (1.0 - wm) * wm2, wm];
                let rm = bary_to_cart(&bary_m, &v);
                let wm_combined = wm * wm2 * area * 4.0;
                let _ = wm_combined; // used below

                // Inner Duffy loop: source point r' (Duffy-transformed on same face)
                for (&rho, &w_rho) in gl_pts.iter().zip(gl_wts.iter()) {
                    for (&theta, &w_theta) in gl_pts.iter().zip(gl_wts.iter()) {
                        let rp = interp3(va, vb, vc, rho, theta);
                        let g = green3d(&rm, &rp, k);
                        let jac_inner = 4.0 * area_sub * rho;
                        let (dot_ff, div_prod) = fm_fn(&rm, &rp);
                        let weight = wm * wm2 * w_rho * w_theta * jac_inner * (area * 4.0);
                        let _ = (bm_pt_idx, bm_pt_jdx);
                        a_sum   += g * dot_ff   * weight;
                        phi_sum += g * div_prod * weight;
                    }
                }
            }
        }
    }

    (a_sum, phi_sum)
}

/// Sauter-Schwab shared-edge for RWG EFIE.
fn rwg_efie_sauter_schwab_edge(
    face_m: &TriFace,
    face_n: &TriFace,
    fm_fn: &dyn Fn(&[f64; 3], &[f64; 3]) -> (f64, f64),
    nodes: &[[f64; 3]],
    k: f64,
    n_gauss: usize,
) -> (Complex64, Complex64) {
    let shared: Vec<usize> = face_m.nodes.iter()
        .filter(|&&nm| face_n.nodes.contains(&nm))
        .copied().collect();
    let unshared_m = face_m.nodes.iter().find(|&&nm| !shared.contains(&nm)).copied().unwrap();
    let unshared_n = face_n.nodes.iter().find(|&&nn| !shared.contains(&nn)).copied().unwrap();

    let vm = [nodes[shared[0]], nodes[shared[1]], nodes[unshared_m]];
    let vn = [nodes[shared[0]], nodes[shared[1]], nodes[unshared_n]];
    let area_m = sub_triangle_area(&vm[0], &vm[1], &vm[2]);
    let area_n = sub_triangle_area(&vn[0], &vn[1], &vn[2]);

    let (gl, gw) = gauss_legendre_1d(n_gauss);
    let mut a_sum = Complex64::ZERO;
    let mut phi_sum = Complex64::ZERO;

    for region in 0usize..5 {
        for (&x1, &w1) in gl.iter().zip(gw.iter()) {
            for (&x2, &w2) in gl.iter().zip(gw.iter()) {
                for (&x3, &w3) in gl.iter().zip(gw.iter()) {
                    for (&x4, &w4) in gl.iter().zip(gw.iter()) {
                        let (xi1, xi2, eta1, eta2, jac_extra) = match region {
                            0 => (x1, x1*x2, x1*x3, x1*x3*x4, x1*x1*x1*x3),
                            1 => (x1, x1*x2, x1*x3*x4, x1*x3, x1*x1*x1*x3),
                            2 => (x1*x2, x1, x1*x3, x1*x3*x4, x1*x1*x1*x3),
                            3 => (x1*x2, x1, x1*x3*x4, x1*x3, x1*x1*x1*x3),
                            _ => (x1, x1*x2*x3, x1*x2, x1*x4, x1*x1*x1*x2),
                        };
                        let bm = [1.0-xi1, xi1-xi2, xi2];
                        let bn = [1.0-eta1, eta1-eta2, eta2];
                        if bm[0] < 0.0 || bm[1] < 0.0 || bm[2] < 0.0 { continue; }
                        if bn[0] < 0.0 || bn[1] < 0.0 || bn[2] < 0.0 { continue; }

                        let rm = bary_to_cart(&bm, &vm);
                        let rn = bary_to_cart(&bn, &vn);
                        let g = green3d(&rm, &rn, k);
                        let jac = jac_extra * 4.0 * area_m * 4.0 * area_n;
                        let (dot_ff, div_prod) = fm_fn(&rm, &rn);
                        let weight = w1 * w2 * w3 * w4 * jac;
                        a_sum   += g * dot_ff   * weight;
                        phi_sum += g * div_prod * weight;
                    }
                }
            }
        }
    }

    (a_sum, phi_sum)
}

/// Sauter-Schwab shared-vertex for RWG EFIE.
fn rwg_efie_sauter_schwab_vertex(
    face_m: &TriFace,
    face_n: &TriFace,
    fm_fn: &dyn Fn(&[f64; 3], &[f64; 3]) -> (f64, f64),
    nodes: &[[f64; 3]],
    k: f64,
    n_gauss: usize,
) -> (Complex64, Complex64) {
    let shared = face_m.nodes.iter()
        .find(|&&nm| face_n.nodes.contains(&nm))
        .copied().unwrap();
    let vm = reorder_with_first(face_m.nodes, shared, nodes);
    let vn = reorder_with_first(face_n.nodes, shared, nodes);
    let area_m = sub_triangle_area(&vm[0], &vm[1], &vm[2]);
    let area_n = sub_triangle_area(&vn[0], &vn[1], &vn[2]);

    let (gl, gw) = gauss_legendre_1d(n_gauss);
    let mut a_sum = Complex64::ZERO;
    let mut phi_sum = Complex64::ZERO;

    for region in 0usize..2 {
        for (&x1, &w1) in gl.iter().zip(gw.iter()) {
            for (&x2, &w2) in gl.iter().zip(gw.iter()) {
                for (&x3, &w3) in gl.iter().zip(gw.iter()) {
                    for (&x4, &w4) in gl.iter().zip(gw.iter()) {
                        let (xi1, xi2, eta1, eta2, jac_extra) = match region {
                            0 => (x1, x1*x2, x1*x3, x1*x3*x4, x1*x1*x1*x3),
                            _ => (x1*x2, x1*x2*x3, x1, x1*x4, x1*x1*x2),
                        };
                        let bm = [1.0-xi1, xi1-xi2, xi2];
                        let bn = [1.0-eta1, eta1-eta2, eta2];
                        if bm[0] < 0.0 || bm[1] < 0.0 || bm[2] < 0.0 { continue; }
                        if bn[0] < 0.0 || bn[1] < 0.0 || bn[2] < 0.0 { continue; }

                        let rm = bary_to_cart(&bm, &vm);
                        let rn = bary_to_cart(&bn, &vn);
                        let g = green3d(&rm, &rn, k);
                        let jac = jac_extra * 4.0 * area_m * 4.0 * area_n;
                        let (dot_ff, div_prod) = fm_fn(&rm, &rn);
                        let weight = w1 * w2 * w3 * w4 * jac;
                        a_sum   += g * dot_ff   * weight;
                        phi_sum += g * div_prod * weight;
                    }
                }
            }
        }
    }

    (a_sum, phi_sum)
}

/// Convenience wrapper used by the assembler for the diagonal self-term.
/// Returns Z_self = -jωμ₀ * ∫∫ G dS' dS * (nothing — caller multiplies).
pub fn zmn_self_duffy_pulse(
    face: &TriFace,
    nodes: &[[f64; 3]],
    k: f64,
    omega_mu0: f64,
    n_gauss: usize,
) -> Complex64 {
    let integral = duffy_self_integral(face, nodes, k, n_gauss);
    Complex64::new(0.0, -omega_mu0) * integral * face.area
}

// ---------------------------------------------------------------------------
// Duffy self-integral  (identical triangles)
// ---------------------------------------------------------------------------
//
// Split the reference triangle T into 3 sub-triangles anchored at each vertex.
// For each sub-triangle apply the Duffy polar substitution:
//   (ρ, θ) ∈ [0,1]² ,  r' = v0 + ρ(θ·e1 + (1-θ)·e2)
//   Jacobian includes factor ρ that cancels the 1/R singularity.
// The observation point r is taken at the face centroid (sufficient for pulse basis
// where f_m = constant; accuracy improved by outer Gauss loop over m).

fn duffy_self_integral(
    face: &TriFace,
    nodes: &[[f64; 3]],
    k: f64,
    n_gauss: usize,
) -> Complex64 {
    let (gl_pts, gl_wts) = gauss_legendre_1d(n_gauss);
    let [i0, i1, i2] = face.nodes;
    let v = [nodes[i0], nodes[i1], nodes[i2]];
    let area = face.area;

    let mut sum = Complex64::ZERO;

    for &pivot in &[0usize, 1, 2] {
        let va = &v[pivot];
        let vb = &v[(pivot + 1) % 3];
        let vc = &v[(pivot + 2) % 3];
        let area_sub = sub_triangle_area(va, vb, vc);

        for (&rho, &w_rho) in gl_pts.iter().zip(gl_wts.iter()) {
            for (&theta, &w_theta) in gl_pts.iter().zip(gl_wts.iter()) {
                // Source point r' via Duffy coords
                let rp = interp3(va, vb, vc, rho, theta);
                // Observation: face centroid (pulse basis — constant over face)
                let r_obs = face.centroid;
                let g = green3d(&r_obs, &rp, k);
                // Jacobian: 4 * area_sub * rho
                let jac = 4.0 * area_sub * rho;
                sum += g * (w_rho * w_theta * jac);
            }
        }
    }

    // Return ∫∫ G dS' dS; the area factor for the observation (pulse basis) is applied by caller
    sum * area
}

// ---------------------------------------------------------------------------
// Sauter-Schwab shared-edge integral
// ---------------------------------------------------------------------------
//
// Two triangles T_m, T_n sharing exactly one edge.
// Reorder vertices so shared edge is (v0,v1) in both triangles:
//   T_m: v0, v1, v2_m    T_n: v0, v1, v2_n
//
// Apply Sauter-Schwab 4-D transformation (Sauter & Schwab §5.3.2):
// 5 sub-regions; combined with n_gauss^4 tensor product quadrature.
//
// For the scalar Green function G = exp(-jkR)/(4πR):
//   I = ∫_Tm ∫_Tn G(r,r') dS' dS
// The 1/R singularity is integrable in 4-D (bounded after the SS transform).

fn sauter_schwab_shared_edge(
    face_m: &TriFace,
    face_n: &TriFace,
    nodes: &[[f64; 3]],
    k: f64,
    n_gauss: usize,
) -> Complex64 {
    // Find the two shared nodes
    let shared: Vec<usize> = face_m.nodes.iter()
        .filter(|&&nm| face_n.nodes.contains(&nm))
        .copied().collect();
    let unshared_m = face_m.nodes.iter().find(|&&nm| !shared.contains(&nm)).copied().unwrap();
    let unshared_n = face_n.nodes.iter().find(|&&nn| !shared.contains(&nn)).copied().unwrap();

    // Local vertex arrays: v[0,1] = shared edge, v[2] = free vertex
    let vm = [nodes[shared[0]], nodes[shared[1]], nodes[unshared_m]];
    let vn = [nodes[shared[0]], nodes[shared[1]], nodes[unshared_n]];

    let (gl, gw) = gauss_legendre_1d(n_gauss);
    let mut sum = Complex64::ZERO;

    // Sauter-Schwab edge rule: 5 sub-regions, each mapped to [0,1]^4
    // We use the simplified "direct evaluation" form from Mesh-based quad.
    // Reference: Graglia & Lombardi (2008), Table I; or Eibert & Hansen (1995).
    //
    // Variables: (x1,x2,x3,x4) ∈ [0,1]^4
    // Region I:   ξ₁=x1, ξ₂=x1·x2, η₁=x1·x3, η₂=x1·x3·x4  ; J = x1³·x3
    // Region II:  ξ₁=x1, ξ₂=x1·x2, η₁=x1·x3·x4, η₂=x1·x3   ; J = x1³·x3
    // Region III: ξ₁=x1·x2, ξ₂=x1, η₁=x1·x3, η₂=x1·x3·x4   ; J = x1³·x3
    // Region IV:  ξ₁=x1·x2, ξ₂=x1, η₁=x1·x3·x4, η₂=x1·x3   ; J = x1³·x3
    // Region V:   ξ₁=x1, ξ₂=x1·x2·x3, η₁=x1·x2, η₂=x1·x4   ; J = x1³·x2
    //
    // Barycentric coords for T_m: (1-ξ₁, ξ₁-ξ₂, ξ₂) → vertices v0,v1,v2_m
    // Barycentric coords for T_n: (1-η₁, η₁-η₂, η₂) → vertices v0,v1,v2_n

    for region in 0usize..5 {
        for (&x1, &w1) in gl.iter().zip(gw.iter()) {
            for (&x2, &w2) in gl.iter().zip(gw.iter()) {
                for (&x3, &w3) in gl.iter().zip(gw.iter()) {
                    for (&x4, &w4) in gl.iter().zip(gw.iter()) {
                        let (xi1, xi2, eta1, eta2, jac_extra) = match region {
                            0 => (x1, x1*x2, x1*x3, x1*x3*x4, x1*x1*x1*x3),
                            1 => (x1, x1*x2, x1*x3*x4, x1*x3, x1*x1*x1*x3),
                            2 => (x1*x2, x1, x1*x3, x1*x3*x4, x1*x1*x1*x3),
                            3 => (x1*x2, x1, x1*x3*x4, x1*x3, x1*x1*x1*x3),
                            _ => (x1, x1*x2*x3, x1*x2, x1*x4, x1*x1*x1*x2),
                        };

                        let bm = [1.0-xi1, xi1-xi2, xi2];
                        let bn = [1.0-eta1, eta1-eta2, eta2];

                        if bm[0] < 0.0 || bm[1] < 0.0 || bm[2] < 0.0 { continue; }
                        if bn[0] < 0.0 || bn[1] < 0.0 || bn[2] < 0.0 { continue; }

                        let rm = bary_to_cart(&bm, &vm);
                        let rn = bary_to_cart(&bn, &vn);

                        let g = green3d(&rm, &rn, k);
                        let area_m = sub_triangle_area(&vm[0], &vm[1], &vm[2]);
                        let area_n = sub_triangle_area(&vn[0], &vn[1], &vn[2]);
                        let jac = jac_extra * 4.0 * area_m * 4.0 * area_n;

                        sum += g * (w1 * w2 * w3 * w4 * jac);
                    }
                }
            }
        }
    }

    sum
}

// ---------------------------------------------------------------------------
// Sauter-Schwab shared-vertex integral
// ---------------------------------------------------------------------------
//
// Two triangles sharing exactly one vertex.
// Reorder: shared vertex = v0 in both.
//   T_m: v0, v1_m, v2_m    T_n: v0, v1_n, v2_n
//
// Sauter-Schwab vertex rule: 2 sub-regions, [0,1]^4 each.
// Region I:  ξ₁=x1, ξ₂=x1·x2, η₁=x1·x3·x4, η₂=x1·x3    ; J=x1³·x3
// Region II: ξ₁=x1·x2·x3, ξ₂=x1·x2, η₁=x1, η₂=x1·x4    ; J=x1²·x2
//
// Barycentric for T_m: (1-ξ₁, ξ₁-ξ₂, ξ₂)
// Barycentric for T_n: (1-η₁, η₁-η₂, η₂)

fn sauter_schwab_shared_vertex(
    face_m: &TriFace,
    face_n: &TriFace,
    nodes: &[[f64; 3]],
    k: f64,
    n_gauss: usize,
) -> Complex64 {
    // Find shared node
    let shared = face_m.nodes.iter()
        .find(|&&nm| face_n.nodes.contains(&nm))
        .copied().unwrap();

    // Reorder so shared vertex is first
    let vm = reorder_with_first(face_m.nodes, shared, nodes);
    let vn = reorder_with_first(face_n.nodes, shared, nodes);

    let (gl, gw) = gauss_legendre_1d(n_gauss);
    let mut sum = Complex64::ZERO;

    for region in 0usize..2 {
        for (&x1, &w1) in gl.iter().zip(gw.iter()) {
            for (&x2, &w2) in gl.iter().zip(gw.iter()) {
                for (&x3, &w3) in gl.iter().zip(gw.iter()) {
                    for (&x4, &w4) in gl.iter().zip(gw.iter()) {
                        // Sauter-Schwab vertex rule (Sauter&Schwab §5.3.3)
                        // Bary for T_m: (1-ξ₁, ξ₁-ξ₂, ξ₂)
                        // Bary for T_n: (1-η₁, η₁-η₂, η₂)
                        // Both must satisfy 0 ≤ ξ₂ ≤ ξ₁ ≤ 1 (and same for η)
                        let (xi1, xi2, eta1, eta2, jac_extra) = match region {
                            // Region I:  ξ₁=x1, ξ₂=x1·x2, η₁=x1·x3, η₂=x1·x3·x4
                            0 => (x1, x1*x2, x1*x3, x1*x3*x4, x1*x1*x1*x3),
                            // Region II: ξ₁=x1·x2, ξ₂=x1·x2·x3, η₁=x1, η₂=x1·x4
                            _ => (x1*x2, x1*x2*x3, x1, x1*x4, x1*x1*x2),
                        };

                        let bm = [1.0-xi1, xi1-xi2, xi2];
                        let bn = [1.0-eta1, eta1-eta2, eta2];

                        if bm[0] < 0.0 || bm[1] < 0.0 || bm[2] < 0.0 { continue; }
                        if bn[0] < 0.0 || bn[1] < 0.0 || bn[2] < 0.0 { continue; }

                        let rm = bary_to_cart(&bm, &vm);
                        let rn = bary_to_cart(&bn, &vn);

                        let g = green3d(&rm, &rn, k);
                        let area_m = sub_triangle_area(&vm[0], &vm[1], &vm[2]);
                        let area_n = sub_triangle_area(&vn[0], &vn[1], &vn[2]);
                        let jac = jac_extra * 4.0 * area_m * 4.0 * area_n;

                        sum += g * (w1 * w2 * w3 * w4 * jac);
                    }
                }
            }
        }
    }

    sum
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn bary_to_cart(b: &[f64; 3], v: &[[f64; 3]; 3]) -> [f64; 3] {
    [
        b[0]*v[0][0] + b[1]*v[1][0] + b[2]*v[2][0],
        b[0]*v[0][1] + b[1]*v[1][1] + b[2]*v[2][1],
        b[0]*v[0][2] + b[1]*v[1][2] + b[2]*v[2][2],
    ]
}

fn reorder_with_first(nodes_idx: [usize; 3], first: usize, nodes: &[[f64; 3]]) -> [[f64; 3]; 3] {
    let mut ordered = nodes_idx;
    if ordered[1] == first { ordered.swap(0, 1); }
    else if ordered[2] == first { ordered.swap(0, 2); }
    [nodes[ordered[0]], nodes[ordered[1]], nodes[ordered[2]]]
}

/// Interpolation: r = va + rho*(theta*(vb-va) + (1-theta)*(vc-va))
fn interp3(va: &[f64; 3], vb: &[f64; 3], vc: &[f64; 3], rho: f64, theta: f64) -> [f64; 3] {
    [
        va[0] + rho*(theta*(vb[0]-va[0]) + (1.0-theta)*(vc[0]-va[0])),
        va[1] + rho*(theta*(vb[1]-va[1]) + (1.0-theta)*(vc[1]-va[1])),
        va[2] + rho*(theta*(vb[2]-va[2]) + (1.0-theta)*(vc[2]-va[2])),
    ]
}

fn sub_triangle_area(va: &[f64; 3], vb: &[f64; 3], vc: &[f64; 3]) -> f64 {
    let e1 = [vb[0]-va[0], vb[1]-va[1], vb[2]-va[2]];
    let e2 = [vc[0]-va[0], vc[1]-va[1], vc[2]-va[2]];
    let cx = e1[1]*e2[2] - e1[2]*e2[1];
    let cy = e1[2]*e2[0] - e1[0]*e2[2];
    let cz = e1[0]*e2[1] - e1[1]*e2[0];
    0.5 * (cx*cx + cy*cy + cz*cz).sqrt()
}

/// Gauss-Legendre quadrature points and weights on [0,1], n ∈ {1..5}.
pub fn gauss_legendre_1d(n: usize) -> (Vec<f64>, Vec<f64>) {
    let (pts_m1, wts_m1): (Vec<f64>, Vec<f64>) = match n {
        1 => (vec![0.0], vec![2.0]),
        2 => (
            vec![-0.577350269189626, 0.577350269189626],
            vec![1.0, 1.0],
        ),
        3 => (
            vec![-0.774596669241483, 0.0, 0.774596669241483],
            vec![0.555555555555556, 0.888888888888889, 0.555555555555556],
        ),
        4 => (
            vec![-0.861136311594953, -0.339981043584856,
                  0.339981043584856,  0.861136311594953],
            vec![ 0.347854845137454,  0.652145154862626,
                  0.652145154862626,  0.347854845137454],
        ),
        5 => (
            vec![-0.906179845938664, -0.538469310105683, 0.0,
                  0.538469310105683,  0.906179845938664],
            vec![ 0.236926885056189,  0.478628670499366, 0.568888888888889,
                  0.478628670499366,  0.236926885056189],
        ),
        n => panic!("gauss_legendre_1d: unsupported n={}", n),
    };
    let pts: Vec<f64> = pts_m1.iter().map(|&t| (1.0 + t) / 2.0).collect();
    let wts: Vec<f64> = wts_m1.iter().map(|&w| w / 2.0).collect();
    (pts, wts)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::surface_mesh::tri_geometry;

    fn make_face(p0: [f64;3], p1: [f64;3], p2: [f64;3]) -> (TriFace, Vec<[f64;3]>) {
        let nodes = vec![p0, p1, p2];
        let (centroid, normal, area) = tri_geometry(&p0, &p1, &p2);
        let face = TriFace { nodes: [0,1,2], centroid, normal, area };
        (face, nodes)
    }

    #[test]
    fn classify_identical() {
        let (f, _) = make_face([0.0,0.0,0.0],[1.0,0.0,0.0],[0.0,1.0,0.0]);
        assert_eq!(classify_pair(&f, &f), TriPairType::Identical);
    }

    #[test]
    fn classify_shared_edge() {
        // T_m: nodes 0,1,2   T_n: nodes 0,1,3
        let nm = TriFace { nodes:[0,1,2], centroid:[0.0;3], normal:[0.0,0.0,1.0], area:0.5 };
        let nn = TriFace { nodes:[0,1,3], centroid:[0.0;3], normal:[0.0,0.0,1.0], area:0.5 };
        assert_eq!(classify_pair(&nm, &nn), TriPairType::SharedEdge);
    }

    #[test]
    fn classify_shared_vertex() {
        let nm = TriFace { nodes:[0,1,2], centroid:[0.0;3], normal:[0.0,0.0,1.0], area:0.5 };
        let nn = TriFace { nodes:[0,3,4], centroid:[0.0;3], normal:[0.0,0.0,1.0], area:0.5 };
        assert_eq!(classify_pair(&nm, &nn), TriPairType::SharedVertex);
    }

    #[test]
    fn classify_disjoint() {
        let nm = TriFace { nodes:[0,1,2], centroid:[0.0;3], normal:[0.0,0.0,1.0], area:0.5 };
        let nn = TriFace { nodes:[3,4,5], centroid:[0.0;3], normal:[0.0,0.0,1.0], area:0.5 };
        assert_eq!(classify_pair(&nm, &nn), TriPairType::Disjoint);
    }

    #[test]
    fn duffy_self_is_finite_and_nonzero() {
        let (face, nodes) = make_face([0.0,0.0,0.0],[1.0,0.0,0.0],[0.0,1.0,0.0]);
        let k = 1.0;
        let omega_mu0 = 1.0;
        let z = zmn_self_duffy_pulse(&face, &nodes, k, omega_mu0, 4);
        assert!(z.norm() > 0.0, "self-integral should be nonzero");
        assert!(z.norm().is_finite(), "self-integral should be finite");
    }

    #[test]
    fn shared_edge_integral_finite() {
        // Two triangles sharing edge (0,1)
        let nodes: Vec<[f64;3]> = vec![
            [0.0,0.0,0.0],[1.0,0.0,0.0],[0.0,1.0,0.0],[1.0,1.0,0.0]
        ];
        let (c0,n0,a0) = tri_geometry(&nodes[0],&nodes[1],&nodes[2]);
        let (c1,n1,a1) = tri_geometry(&nodes[0],&nodes[1],&nodes[3]);
        let fm = TriFace { nodes:[0,1,2], centroid:c0, normal:n0, area:a0 };
        let fn_ = TriFace { nodes:[0,1,3], centroid:c1, normal:n1, area:a1 };

        let val = sauter_schwab_shared_edge(&fm, &fn_, &nodes, 1.0, 3);
        assert!(val.norm().is_finite(), "shared-edge integral finite");
        assert!(val.norm() > 0.0, "shared-edge integral nonzero");
    }

    #[test]
    fn shared_vertex_integral_finite() {
        let nodes: Vec<[f64;3]> = vec![
            [0.0,0.0,0.0],[1.0,0.0,0.0],[0.0,1.0,0.0],
            [-1.0,0.0,0.0],[0.0,-1.0,0.0],
        ];
        let (c0,n0,a0) = tri_geometry(&nodes[0],&nodes[1],&nodes[2]);
        let (c1,n1,a1) = tri_geometry(&nodes[0],&nodes[3],&nodes[4]);
        let fm = TriFace { nodes:[0,1,2], centroid:c0, normal:n0, area:a0 };
        let fn_ = TriFace { nodes:[0,3,4], centroid:c1, normal:n1, area:a1 };

        let val = sauter_schwab_shared_vertex(&fm, &fn_, &nodes, 1.0, 3);
        assert!(val.norm().is_finite(), "shared-vertex integral finite");
        assert!(val.norm() > 0.0, "shared-vertex integral nonzero");
    }
}
