//! HCurl eigenmode integration test: 2-D PEC cavity resonance.
//!
//! Verifies that the Nedelec edge-element eigenmode solver recovers the
//! TE₁₀ mode frequency for a rectangular PEC cavity to within 3 %.

use rem_config::{load_config_from_str, ConfigFormat};
use rem_materials::DomainMap;
use rem_mesh::{gen::rect_msh, gmsh::read_msh_str, RemMesh};
use rem_parallel::NoComm;

const C0: f64 = 2.997_924_58e8;

#[test]
fn hcurl_eigenmode_pec_cavity_te10_mode_2d() {
    // 2-D rectangular PEC cavity, a × b.
    // TE₁₀ (fundamental): f = c₀ / (2a)
    let a_m = 1.0_f64;
    let b_m = 0.5_f64;
    let f_te10 = C0 / (2.0 * a_m); // 150 MHz

    // Generate 2-D rectangle mesh with all four sides tagged separately.
    let msh = rect_msh(a_m, b_m, 16, 8, 1, 2, 3, 4, 10);
    let raw = read_msh_str(&msh).expect("rect_msh should parse");

    let json = format!(
        r#"{{
            "Problem": {{"Type": "Eigenmode"}},
            "Model":   {{"Mesh": "cavity.msh", "L0": 1.0}},
            "Domains": {{
                "Materials": [{{"Attributes": [10], "Permittivity": 1.0, "Permeability": 1.0}}]
            }},
            "Boundaries": {{
                "PEC": {{"Attributes": [1, 2, 3, 4]}}
            }},
            "Solver": {{
                "Order": 1,
                "Eigenmode": {{"N": 5, "Target": {tgt:.1e}, "Tol": 1e-8}},
                "Linear": {{"Tol": 1e-10, "MaxIter": 2000}}
            }}
        }}"#,
        tgt = f_te10
    );
    let config = load_config_from_str(&json, ConfigFormat::Json).expect("config should parse");
    let mesh = RemMesh::from_raw(raw, config.model.l0).expect("RemMesh::from_raw failed");
    let domain_map = DomainMap::from_config(&config).expect("DomainMap::from_config failed");

    let result =
        rem_eigenmode::solve(&config, &mesh, &domain_map, &NoComm)
            .expect("HCurl eigenmode solve failed");

    assert!(
        !result.frequencies_hz.is_empty(),
        "expected at least one eigenmode"
    );
    assert!(result.is_hcurl, "solver should report HCurl path");

    let f0 = result.frequencies_hz[0];
    let rel_err = (f0 - f_te10).abs() / f_te10;
    assert!(
        rel_err < 0.03,
        "TE₁₀ frequency: got {:.3e} Hz, expected {:.3e} Hz, rel err {:.4}%",
        f0,
        f_te10,
        rel_err * 100.0
    );
}
