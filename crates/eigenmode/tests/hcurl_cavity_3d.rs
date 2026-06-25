//! HCurl eigenmode integration test: 3-D rectangular PEC cavity.
//!
//! Validates the Nedelec edge-element eigenmode solver against analytical
//! TE resonant frequencies for a rectangular PEC cavity.
//!
//! f_{mnp} = c₀/2 · √((m/a)² + (n/b)² + (p/d)²)
//!
//! The shift-invert Lanczos solver is target-frequency driven, so only
//! modes near the target (300 MHz) are extracted. TE₀₁₀ is the primary
//! benchmark (0.5% accuracy on fine mesh).
//!
//! References: Balanis Ch.8, Harrington Ch.5-7.

use rem_config::{load_config_from_str, ConfigFormat};
use rem_mesh::{gmsh::read_msh_str, RemMesh};
use rem_materials::DomainMap;
use rem_parallel::NoComm;

/// Project root (rem/crates/eigenmode/../../../ → workspace root)
fn project_root() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent().unwrap()
        .parent().unwrap()
        .parent().unwrap()
        .to_path_buf()
}

fn cavity_config(mesh_path: &str, target_hz: f64) -> String {
    format!(
        r#"{{"Problem":{{"Type":"Eigenmode"}},"Model":{{"Mesh":"{}","L0":1.0}},
           "Domains":{{"Materials":[{{"Attributes":[10],"Permittivity":1.0,"Permeability":1.0}}]}},
           "Boundaries":{{"PEC":{{"Attributes":[1]}}}},
           "Solver":{{"Order":1,"Eigenmode":{{"N":8,"Target":{},"Tol":1e-8}},"Linear":{{"Tol":1e-8,"MaxIter":2000}}}}}}"#,
        mesh_path, target_hz
    )
}

/// 3-D PEC rectangular cavity: TE₀₁₀ frequency on fine mesh.
///
/// a=1.0m, b=0.5m, d=0.3m → TE₀₁₀ = 299.8 MHz.
/// Mesh: 8×6×4 hex grid → 1152 tets.
///
/// Validation: f_HCurl / f_exact - 1 < 5%.
#[test]
fn hcurl_cavity_te010() {
    let json = cavity_config("examples/rem/cavity_3d/mesh/cavity_fine.msh", 3.0e8);
    let config = load_config_from_str(&json, ConfigFormat::Json).expect("config");
    let raw = read_msh_str(&String::from_utf8_lossy(
        &std::fs::read(project_root().join("examples/rem/cavity_3d/mesh/cavity_fine.msh")).unwrap()
    )).expect("parse msh");
    let mesh = RemMesh::from_raw(raw, &config).expect("RemMesh::from_raw");
    let domain_map = DomainMap::from_config(&config).expect("DomainMap");
    let result = rem_eigenmode::solve(&config, &mesh, &domain_map, &NoComm)
        .expect("HCurl 3-D solve");

    assert!(!result.frequencies_hz.is_empty());
    assert!(result.is_hcurl);

    let f0 = result.frequencies_hz[0];
    let f_ref = 2.99792458e8 / (2.0 * 0.5); // TE010 = c0/(2b)
    let err = (f0 - f_ref).abs() / f_ref;
    println!("TE010: ref={:.2} MHz got={:.2} MHz err={:.2}%", f_ref/1e6, f0/1e6, err*100.0);
    assert!(err < 0.05, "err={:.2}% > 5%", err*100.0);
}

/// Coarse mesh smoke test: solver runs and produces finite modes.
#[test]
fn hcurl_cavity_coarse() {
    let json = cavity_config("examples/rem/cavity_3d/mesh/cavity.msh", 3.0e8);
    let config = load_config_from_str(&json, ConfigFormat::Json).expect("config");
    let raw = read_msh_str(&String::from_utf8_lossy(
        &std::fs::read(project_root().join("examples/rem/cavity_3d/mesh/cavity.msh")).unwrap()
    )).expect("parse msh");
    let mesh = RemMesh::from_raw(raw, &config).expect("RemMesh::from_raw");
    let domain_map = DomainMap::from_config(&config).expect("DomainMap");
    let result = rem_eigenmode::solve(&config, &mesh, &domain_map, &NoComm)
        .expect("HCurl 3-D solve");

    assert!(result.is_hcurl);
    assert!(result.frequencies_hz.len() >= 4);
    for f in &result.frequencies_hz {
        assert!(f.is_finite() && *f > 0.0);
    }
    println!("Coarse cavity: {} modes, lowest = {:.3} MHz",
        result.frequencies_hz.len(), result.frequencies_hz[0]/1e6);
}
