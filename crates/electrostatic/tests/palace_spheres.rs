//! Palace "spheres" example — 2-D coaxial cross-section.
//!
//! Reference:  https://awslabs.github.io/palace/stable/examples/spheres/
//!
//! Geometry: annular domain, inner conductor r_i = 1 mm (Terminal, φ = 1 V),
//!           outer conductor r_o = 4 mm (Ground, φ = 0 V), ε_r = 1 (vacuum).
//!
//! Analytical solution:
//!   φ(r) = ln(r_o / r) / ln(r_o / r_i)   where r = √(x²+y²)
//!
//!   C/L  = 2π ε₀ / ln(r_o / r_i)
//!        = 2π × 8.854187817e-12 / ln(4)
//!        ≈ 40.12 pF/m
//!
//! Palace "spheres" result (for reference, 3-D):
//!   C_11 =  +1.2305  pF     (self-capacitance, inner sphere)
//!   C_22 =  +2.4315  pF     (self-capacitance, outer shell)
//!   C_12 =  −0.4946  pF     (mutual capacitance)
//!
//! This 2-D test verifies:
//!   1. φ L2-error < 0.5 % vs. analytical
//!   2. C/L within 1 % of 2πε₀/ln(4) ≈ 40.12 pF/m
//!   3. |E|-field recovered at representative nodes within 2 % of −dφ/dr

use rem_config::{load_config_from_str, ConfigFormat};
use rem_core::constants::EPS0;
use rem_electrostatic::{postprocess, solve_one};
use rem_materials::DomainMap;
use rem_mesh::{gen::annular_msh, gmsh::read_msh_str, RemMesh};
use std::f64::consts::PI;

// ── mesh / config parameters ─────────────────────────────────────────────────

const R_I_MM: f64 = 1.0; // inner radius [mesh units = mm]
const R_O_MM: f64 = 4.0; // outer radius
const L0: f64 = 1.0e-3;  // config L0: mm → m

const R_I: f64 = R_I_MM * L0; // 1e-3 m
const R_O: f64 = R_O_MM * L0; // 4e-3 m

// Physical group tags
const TAG_INNER: u32 = 1; // Terminal index 1 → φ = 1 V
const TAG_OUTER: u32 = 2; // Ground → φ = 0 V
const TAG_VOL:   u32 = 10;

// Mesh resolution
const N_R: usize     = 10;
const N_THETA: usize = 32;

fn config() -> rem_config::PalaceConfig {
    let s = format!(
        r#"{{
            "Problem": {{"Type": "Electrostatic"}},
            "Model":   {{"Mesh": "coaxial_2d.msh", "L0": {}}},
            "Domains": {{
                "Materials": [{{"Attributes": [{}], "Permittivity": 1.0}}]
            }},
            "Boundaries": {{
                "Ground":   {{"Attributes": [{}]}},
                "Terminal": [{{"Index": 1, "Attributes": [{}]}}]
            }},
            "Solver": {{"Linear": {{"Tol": 1e-10, "MaxIter": 500}}}}
        }}"#,
        L0, TAG_VOL, TAG_OUTER, TAG_INNER
    );
    load_config_from_str(&s, ConfigFormat::Json).unwrap()
}

fn solve() -> (rem_mesh::RemMesh, Vec<f64>) {
    let msh = annular_msh(R_I_MM, R_O_MM, N_R, N_THETA, TAG_INNER, TAG_OUTER, TAG_VOL);
    let raw = read_msh_str(&msh).expect("annular mesh should parse");
    let cfg = config();
    let mesh = RemMesh::from_raw(raw, &cfg).expect("RemMesh::from_raw failed");
    let dm   = DomainMap::from_config(&cfg).unwrap();
    let phi  = solve_one(&cfg, &mesh, &dm, Some(1), 1.0, &rem_parallel::NoComm).expect("solve failed");
    (mesh, phi)
}

/// Analytical potential φ(r) for the coaxial problem.
#[inline]
fn phi_exact(r: f64) -> f64 {
    (R_O / r).ln() / (R_O / R_I).ln()
}

/// Analytical capacitance per unit depth [F/m].
fn cap_analytical() -> f64 {
    2.0 * PI * EPS0 / (R_O / R_I).ln()
}

// ── Test 1: L2 relative error of φ ──────────────────────────────────────────

/// φ is accurate to < 0.5 % relative L2 error on the N_R=10, N_THETA=32 mesh.
#[test]
fn potential_l2_error_below_half_percent() {
    let (mesh, phi) = solve();

    let mut err_sq = 0.0_f64;
    let mut ref_sq = 0.0_f64;
    for (i, node) in mesh.nodes.iter().enumerate() {
        let r = (node.x * node.x + node.y * node.y).sqrt();
        let exact = phi_exact(r);
        let diff  = phi[i] - exact;
        err_sq += diff * diff;
        ref_sq += exact * exact;
    }
    let rel_err = (err_sq / ref_sq).sqrt();
    assert!(
        rel_err < 0.005,
        "L2 relative error = {:.4}%  (limit 0.5%)",
        rel_err * 100.0
    );
}

// ── Test 2: capacitance per unit depth ──────────────────────────────────────

/// C/L extracted from field energy agrees with analytical within 1 %.
/// Energy U = ½ C V² → C/L = 2U (V = 1 V).
#[test]
fn capacitance_per_unit_depth_within_1pct() {
    let (mesh, phi) = solve();
    let cfg = config();
    let dm  = DomainMap::from_config(&cfg).unwrap();

    let eps_fn = |tag: u32| dm.get(tag).epsilon_abs();
    let energy = postprocess::electrostatic_energy(&phi, &mesh, eps_fn);
    let c_fem  = 2.0 * energy; // C/L = 2U/V²  (V = 1)
    let c_ana  = cap_analytical();

    let rel_err = (c_fem - c_ana).abs() / c_ana;
    assert!(
        rel_err < 0.01,
        "C/L: FEM = {:.4} pF/m,  analytical = {:.4} pF/m,  err = {:.3}%",
        c_fem * 1e12,
        c_ana * 1e12,
        rel_err * 100.0
    );
}

// ── Test 3: Electric field magnitude ────────────────────────────────────────

/// |E(r)| = 1 / (r ln(r_o/r_i)).
/// Check recovered nodal E-field at inner-ring nodes within 2 %.
#[test]
fn e_field_radial_magnitude() {
    let (mesh, phi) = solve();

    let e_field = postprocess::gradient_recovery(&phi, &mesh);

    let ln_ratio = (R_O / R_I).ln();
    let mut max_err_pct = 0.0_f64;

    for (i, node) in mesh.nodes.iter().enumerate() {
        let r = (node.x * node.x + node.y * node.y).sqrt();
        // Only check nodes clearly away from both boundaries (0.3 < r/R_I < 3.5)
        if r < 1.1 * R_I || r > 0.9 * R_O {
            continue;
        }
        let e_exact = 1.0 / (r * ln_ratio);
        let e_fem   = (e_field[i][0].powi(2) + e_field[i][1].powi(2)).sqrt();
        let err_pct = (e_fem - e_exact).abs() / e_exact * 100.0;
        if err_pct > max_err_pct {
            max_err_pct = err_pct;
        }
    }

    assert!(
        max_err_pct < 2.0,
        "Max |E| error at interior nodes: {:.3}%  (limit 2%)",
        max_err_pct
    );
}

// ── Test 4: Boundary conditions exactly satisfied ───────────────────────────

/// All inner-ring nodes have φ = 1 V, all outer-ring nodes have φ = 0 V.
#[test]
fn boundary_conditions_exactly_satisfied() {
    let (mesh, phi) = solve();

    for (i, node) in mesh.nodes.iter().enumerate() {
        let r = (node.x * node.x + node.y * node.y).sqrt();
        if (r - R_I).abs() < 1e-8 {
            assert!(
                (phi[i] - 1.0).abs() < 1e-9,
                "inner BC: node {} r={:.4e} φ={:.6}",
                i, r, phi[i]
            );
        }
        if (r - R_O).abs() < 1e-8 {
            assert!(
                phi[i].abs() < 1e-9,
                "outer BC: node {} r={:.4e} φ={:.6}",
                i, r, phi[i]
            );
        }
    }
}

// ── Test 5: config round-trips through Palace JSON format ───────────────────

/// Terminal + Ground boundary are correctly parsed from the config string.
#[test]
fn terminal_bc_parses_in_palace_format() {
    use rem_config::{load_config_from_str, ConfigFormat, ProblemType};
    let s = r#"{
        "Problem": {"Type": "Electrostatic"},
        "Model":   {"Mesh": "x.msh"},
        "Boundaries": {
            "Ground":   {"Attributes": [2]},
            "Terminal": [{"Index": 1, "Attributes": [1]}]
        }
    }"#;
    let cfg = load_config_from_str(s, ConfigFormat::Json).unwrap();
    assert_eq!(cfg.problem.problem_type, ProblemType::Electrostatic);
    assert_eq!(cfg.boundaries.terminal.len(), 1);
    assert_eq!(cfg.boundaries.terminal[0].index, 1);
    assert_eq!(cfg.boundaries.terminal[0].attributes, vec![1]);
    assert!(cfg.boundaries.ground.is_some());
}
