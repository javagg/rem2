//! Palace "parallel_plate" example — 2-D parallel plate capacitor.
//!
//! Geometry: 1 mm x 1 mm square.
//! Bottom (y=0): Ground (0 V)
//! Top    (y=1): Terminal index 1 (1 V)
//! Left/Right: Neumann
//!
//! Analytical solution:
//! phi(y) = y / H
//! C = eps0 * W / H = eps0 (for unit square)

use rem_config::{load_config_from_str, ConfigFormat};
use rem_core::constants::EPS0;
use rem_electrostatic::{postprocess, solve_one};
use rem_parallel::NoComm;
use rem_materials::DomainMap;
use rem_mesh::{gen::rect_msh, gmsh::read_msh_str, RemMesh};

#[test]
fn palace_parallel_plate() {
    let w_mm = 1.0;
    let h_mm = 1.0;
    let l0 = 1.0e-3;
    let n_x = 10;
    let n_y = 10;

    let tag_bottom = 1;
    let tag_top = 2;
    let tag_left = 3;
    let tag_right = 4;
    let tag_vol = 10;

    let s = format!(
        r#"{{
            "Problem": {{"Type": "Electrostatic"}},
            "Model":   {{"Mesh": "plate_2d.msh", "L0": {}}},
            "Domains": {{
                "Materials": [{{"Attributes": [{}], "Permittivity": 1.0}}]
            }},
            "Boundaries": {{
                "Ground":   {{"Attributes": [{}]}},
                "Terminal": [{{"Index": 1, "Attributes": [{}]}}]
            }},
            "Solver": {{"Linear": {{"Tol": 1e-12, "MaxIter": 500}}}}
        }}"#,
        l0, tag_vol, tag_bottom, tag_top
    );

    let cfg = load_config_from_str(&s, ConfigFormat::Json).unwrap();
    let msh = rect_msh(w_mm, h_mm, n_x, n_y, tag_bottom, tag_top, tag_left, tag_right, tag_vol);
    let raw = read_msh_str(&msh).unwrap();
    let mesh = RemMesh::from_raw(raw, &cfg).unwrap();
    let dm = DomainMap::from_config(&cfg).unwrap();

    let phi = solve_one(&cfg, &mesh, &dm, Some(1), 1.0, &NoComm).unwrap();

    // Verify energy: U = 0.5 * eps0 * V^2 * (W/H) = 0.5 * eps0
    let energy = postprocess::electrostatic_energy(&phi, &mesh, |_| EPS0);
    let expected_energy = 0.5 * EPS0;
    assert!((energy - expected_energy).abs() / expected_energy < 1e-3);
}
