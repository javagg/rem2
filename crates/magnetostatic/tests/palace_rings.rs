//! Palace "rings" example — 2-D bimaterial magnetostatic slab.
//!
//! Reference:  https://awslabs.github.io/palace/v0.11/examples/rings/
//!
//! Geometry: 1 mm × 1 mm square domain.
//!   Lower half  (y ∈ [0, 0.5 mm])  :  iron,  μ_r = 1000  (tag 10)
//!   Upper half  (y ∈ [0.5, 1 mm])  :  air,   μ_r = 1     (tag 20)
//!   Bottom (y = 0)   : Ground         → A_z = 0
//!   Top    (y = H)   : SurfaceCurrent → A_z = 1  (index 1)
//!   Left/Right       : Neumann (∂A_z/∂n = 0)
//!
//! Analytical solution (series-reluctance):
//!   ν × ∂A_z/∂y = const  throughout the domain.
//!   A_z at interface (y = H/2):
//!     A_z_mid = μ_r_iron / (μ_r_iron + μ_r_air)
//!             = 1000 / 1001  ≈ 0.999001
//!
//!   Gradient in iron (y ∈ [0, H/2]):
//!     ∂A_z/∂y = 2 μ_r_iron / ((μ_r_iron + μ_r_air) H)
//!             = 2000 / (1001 H)
//!
//!   Gradient in air (y ∈ [H/2, H]):
//!     ∂A_z/∂y = 2 / (1001 H)
//!
//! Palace "rings" result (for reference, 3-D):
//!   L_11 = 4.2585e-11 H  (self-inductance, inner ring)
//!   L_22 = 7.1269e-10 H  (self-inductance, outer ring)
//!   L_12 = 1.9585e-12 H  (mutual inductance)
//!
//! This 2-D test verifies:
//!   1. A_z at the iron/air interface equals analytical value to < 1e-8
//!   2. ∂A_z/∂y in each layer matches the analytical gradient within 0.1 %
//!   3. Magnetic energy computed from the field agrees with analytical within 1 %

use rem_config::{load_config_from_str, ConfigFormat};
use rem_core::constants::MU0;
use rem_electrostatic::postprocess;
use rem_magnetostatic::solve_one;
use rem_parallel::NoComm;
use rem_materials::DomainMap;
use rem_mesh::{gen::rect_bimaterial_msh, gmsh::read_msh_str, RemMesh};

// ── parameters ───────────────────────────────────────────────────────────────

const W_MM: f64 = 1.0;   // width  [mm]
const H_MM: f64 = 1.0;   // height [mm]
const L0:   f64 = 1.0e-3; // mm → m

const H:    f64 = H_MM * L0; // 1e-3 m
const MU_R_IRON: f64 = 1000.0;
const MU_R_AIR:  f64 = 1.0;

const TAG_BOTTOM:     u32 = 1;  // Ground
const TAG_TOP:        u32 = 2;  // SurfaceCurrent index 1
const TAG_IRON:       u32 = 10;
const TAG_AIR:        u32 = 20;

const N_X: usize = 20;
const N_Y: usize = 20; // must be even — 10 cells per layer

/// Analytical A_z at the iron/air interface (y = H/2).
fn az_interface_exact() -> f64 {
    MU_R_IRON / (MU_R_IRON + MU_R_AIR)
}

fn config() -> rem_config::PalaceConfig {
    let s = format!(
        r#"{{
            "Problem": {{"Type": "Magnetostatic"}},
            "Model":   {{"Mesh": "slab_2d.msh", "L0": {l0}}},
            "Domains": {{
                "Materials": [
                    {{"Attributes": [{iron}], "Permeability": {mu_iron}}},
                    {{"Attributes": [{air}],  "Permeability": {mu_air}}}
                ]
            }},
            "Boundaries": {{
                "Ground":         {{"Attributes": [{bot}]}},
                "SurfaceCurrent": [{{"Index": 1, "Attributes": [{top}], "Direction": "+Y"}}]
            }},
            "Solver": {{"Linear": {{"Tol": 1e-12, "MaxIter": 500}}}}
        }}"#,
        l0 = L0,
        iron = TAG_IRON, mu_iron = MU_R_IRON,
        air  = TAG_AIR,  mu_air  = MU_R_AIR,
        bot  = TAG_BOTTOM,
        top  = TAG_TOP,
    );
    load_config_from_str(&s, ConfigFormat::Json).unwrap()
}

fn solve() -> (RemMesh, Vec<f64>) {
    let msh = rect_bimaterial_msh(W_MM, H_MM, N_X, N_Y, TAG_BOTTOM, TAG_TOP, TAG_IRON, TAG_AIR);
    let raw = read_msh_str(&msh).expect("bimaterial mesh should parse");
    let cfg = config();
    let mesh = RemMesh::from_raw(raw, &cfg).expect("RemMesh::from_raw failed");
    let dm   = DomainMap::from_config(&cfg).unwrap();
    let az   = solve_one(&cfg, &mesh, &dm, Some(1), &rem_parallel::NoComm).expect("solve failed");
    (mesh, az)
}

// ── Test 1: Interface value ──────────────────────────────────────────────────

/// A_z at y = H/2 equals 1000/1001 to numerical precision.
/// This verifies the FEM correctly handles the μ-contrast interface condition.
#[test]
fn az_at_interface_equals_analytical() {
    let (mesh, az) = solve();
    let az_mid = az_interface_exact();
    let h_mid  = H / 2.0; // 0.5e-3 m

    let mut max_err = 0.0_f64;
    let mut count   = 0usize;

    for (i, node) in mesh.nodes.iter().enumerate() {
        if (node.y - h_mid).abs() < 1e-10 {
            let err = (az[i] - az_mid).abs();
            if err > max_err { max_err = err; }
            count += 1;
        }
    }

    assert!(count > 0, "No nodes found at y = {:.2e}", h_mid);
    assert!(
        max_err < 1e-8,
        "A_z at interface: max error = {:.2e}  (analytical = {:.8})",
        max_err,
        az_mid
    );
}

// ── Test 2: Gradient in each layer ──────────────────────────────────────────

/// ∂A_z/∂y in the iron and air layers matches analytical gradients within 0.1 %.
#[test]
fn az_gradient_per_layer() {
    let (mesh, az) = solve();

    // gradient_recovery returns −∇A_z (same sign convention as E = −∇φ).
    // So recovered[i][1] = −∂Az/∂y.
    let recovered = postprocess::gradient_recovery(&az, &mesh);

    // Analytical gradients (positive, since A_z increases with y)
    let g_iron_exact = 2.0 * MU_R_IRON / ((MU_R_IRON + MU_R_AIR) * H);
    let g_air_exact  = 2.0 * MU_R_AIR  / ((MU_R_IRON + MU_R_AIR) * H);

    let h_mid = H / 2.0;
    let mut iron_err_pct = 0.0_f64;
    let mut air_err_pct  = 0.0_f64;

    for (i, node) in mesh.nodes.iter().enumerate() {
        // Skip interface and boundary nodes (gradient discontinuous there)
        if node.y < 0.05 * H || node.y > 0.95 * H
            || (node.y - h_mid).abs() < 0.05 * H
        {
            continue;
        }
        // Recovered gradient in y: −recovered[i][1]
        let g_fem = -recovered[i][1];
        if node.y < h_mid {
            // Iron region
            let err = (g_fem - g_iron_exact).abs() / g_iron_exact * 100.0;
            if err > iron_err_pct { iron_err_pct = err; }
        } else {
            // Air region
            let err = (g_fem - g_air_exact).abs() / g_air_exact * 100.0;
            if err > air_err_pct { air_err_pct = err; }
        }
    }

    assert!(
        iron_err_pct < 0.1,
        "Iron ∂Az/∂y error: {:.4}%  (limit 0.1%)",
        iron_err_pct
    );
    assert!(
        air_err_pct < 0.1,
        "Air  ∂Az/∂y error: {:.4}%  (limit 0.1%)",
        air_err_pct
    );
}

// ── Test 3: Magnetic energy ──────────────────────────────────────────────────

/// Total magnetic energy U = ½ ∫ ν |∇A_z|² dΩ.
///
/// Analytical:
///   U = ½ [ ν_iron × g_iron² × (W×H/2) + ν_air × g_air² × (W×H/2) ]
///
///   = ½ W H / (2 (μ_r_iron + μ_r_air)² μ₀)
///     × (4 μ_r_iron² / μ_r_iron + 4 μ_r_air² / μ_r_air)  ... simplifies to:
///
///   U = W H / (μ₀ (μ_r_iron + μ_r_air))
///     = 1e-6 / (4π×10⁻⁷ × 1001)  ≈ 0.7955 J/m
#[test]
fn magnetic_energy_within_1pct() {
    let (mesh, az) = solve();
    let cfg = config();
    let dm  = DomainMap::from_config(&cfg).unwrap();

    let nu_fn  = |tag: u32| dm.get(tag).reluctivity();
    let u_fem  = postprocess::electrostatic_energy(&az, &mesh, nu_fn);

    // Analytical energy per unit depth [J/m]:
    //   Each layer (area = W×H/2) with uniform gradient g_i and reluctivity ν_i:
    //   U_i = ½ ν_i g_i² (W H/2)
    let w  = W_MM * L0;
    let nu_iron = 1.0 / (MU0 * MU_R_IRON);
    let nu_air  = 1.0 / (MU0 * MU_R_AIR);
    let g_iron  = 2.0 * MU_R_IRON / ((MU_R_IRON + MU_R_AIR) * H);
    let g_air   = 2.0 * MU_R_AIR  / ((MU_R_IRON + MU_R_AIR) * H);
    let u_iron  = 0.5 * nu_iron * g_iron.powi(2) * w * H / 2.0;
    let u_air   = 0.5 * nu_air  * g_air.powi(2)  * w * H / 2.0;
    let u_ana   = u_iron + u_air;

    let rel_err = (u_fem - u_ana).abs() / u_ana;
    assert!(
        rel_err < 0.01,
        "Magnetic energy: FEM = {:.6e} J/m,  analytical = {:.6e} J/m,  err = {:.4}%",
        u_fem, u_ana, rel_err * 100.0
    );
}

// ── Test 4: Boundary conditions exactly satisfied ───────────────────────────

/// All bottom nodes have A_z = 0, all top nodes have A_z = 1.
#[test]
fn boundary_conditions_exactly_satisfied() {
    let (mesh, az) = solve();

    for (i, node) in mesh.nodes.iter().enumerate() {
        if node.y < 1e-12 {
            assert!(
                az[i].abs() < 1e-10,
                "bottom BC: node {} y={:.2e}  A_z={:.4e}", i, node.y, az[i]
            );
        }
        if (node.y - H).abs() < 1e-12 {
            assert!(
                (az[i] - 1.0).abs() < 1e-10,
                "top BC: node {} y={:.4e}  A_z={:.6}", i, node.y, az[i]
            );
        }
    }
}

// ── Test 5: Single-material limit (μ_r = 1 everywhere) ─────────────────────

/// With μ_r = 1 in both layers the solution degenerates to A_z(y) = y/H (linear).
#[test]
fn single_material_limit_uniform_az() {
    let s = format!(
        r#"{{
            "Problem": {{"Type": "Magnetostatic"}},
            "Model":   {{"Mesh": "x.msh", "L0": {}}},
            "Domains": {{
                "Materials": [
                    {{"Attributes": [10, 20], "Permeability": 1.0}}
                ]
            }},
            "Boundaries": {{
                "Ground":         {{"Attributes": [1]}},
                "SurfaceCurrent": [{{"Index": 1, "Attributes": [2], "Direction": "+Y"}}]
            }},
            "Solver": {{"Linear": {{"Tol": 1e-12, "MaxIter": 500}}}}
        }}"#,
        L0
    );
    let cfg  = load_config_from_str(&s, ConfigFormat::Json).unwrap();
    let msh  = rect_bimaterial_msh(W_MM, H_MM, N_X, N_Y, TAG_BOTTOM, TAG_TOP, TAG_IRON, TAG_AIR);
    let raw  = read_msh_str(&msh).unwrap();
    let mesh = RemMesh::from_raw(raw, &cfg).unwrap();
    let dm   = DomainMap::from_config(&cfg).unwrap();
    let az   = solve_one(&cfg, &mesh, &dm, Some(1), &NoComm).unwrap();

    for (i, node) in mesh.nodes.iter().enumerate() {
        let exact = node.y / H;
        let err   = (az[i] - exact).abs();
        assert!(
            err < 1e-9,
            "uniform-μ: node {} y={:.4e}  A_z={:.8}  exact={:.8}  err={:.2e}",
            i, node.y, az[i], exact, err
        );
    }
}
