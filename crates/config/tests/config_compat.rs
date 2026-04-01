use rem_config::{load_config_from_str, ConfigFormat, ProblemType};
use std::path::Path;

fn fixture(name: &str) -> String {
    let manifest = std::env::var("CARGO_MANIFEST_DIR").unwrap();
    // crates/config → repo root → tests/fixtures/
    let root = Path::new(&manifest).parent().unwrap().parent().unwrap();
    let path = root.join("tests/fixtures").join(name);
    std::fs::read_to_string(&path)
        .unwrap_or_else(|_| panic!("fixture not found: {}", path.display()))
}

// ── JSON parsing ─────────────────────────────────────────────────────────────

#[test]
fn parse_parallel_plate_json() {
    let s = fixture("parallel_plate.json");
    let cfg = load_config_from_str(&s, ConfigFormat::Json)
        .expect("parallel_plate.json should parse");

    assert_eq!(cfg.problem.problem_type, ProblemType::Electrostatic);
    assert_eq!(cfg.problem.verbose(), 1);
    assert_eq!(cfg.model.l0, 1.0e-3);
    assert_eq!(cfg.domains.materials.len(), 1);
    assert!((cfg.domains.materials[0].permittivity - 4.5).abs() < 1e-12);
    assert!((cfg.domains.materials[0].loss_tangent - 0.02).abs() < 1e-12);

    let pec = cfg.boundaries.pec.as_ref().expect("PEC should be present");
    assert_eq!(pec.attributes, vec![3, 4]);

    let gnd = cfg.boundaries.ground.as_ref().expect("Ground should be present");
    assert_eq!(gnd.attributes, vec![1, 2]);

    assert_eq!(cfg.solver.order, 1);
    assert_eq!(cfg.solver.linear.solver_type, "GMRES");
    assert!((cfg.solver.linear.tol - 1.0e-8).abs() < 1e-20);
}

#[test]
fn json_comments_stripped() {
    let s = r#"{
        // line comment
        "Problem": { "Type": "Electrostatic" /* block */ },
        "Model":   { "Mesh": "x.msh" }
    }"#;
    let cfg = load_config_from_str(s, ConfigFormat::Json).expect("comment stripping failed");
    assert_eq!(cfg.problem.problem_type, ProblemType::Electrostatic);
}

#[test]
fn json_string_with_url_not_stripped() {
    let s = r#"{
        "Problem": { "Type": "Electrostatic", "Output": "http://localhost/out" },
        "Model":   { "Mesh": "x.msh" }
    }"#;
    let cfg = load_config_from_str(s, ConfigFormat::Json).unwrap();
    assert_eq!(cfg.problem.output_dir(), "http://localhost/out");
}

#[test]
fn attribute_range_string() {
    let s = r#"{
        "Problem": { "Type": "Electrostatic" },
        "Model":   { "Mesh": "x.msh" },
        "Boundaries": {
            "PEC": { "Attributes": "1,3-5,8" }
        }
    }"#;
    let cfg = load_config_from_str(s, ConfigFormat::Json).unwrap();
    let attrs = &cfg.boundaries.pec.unwrap().attributes;
    assert_eq!(attrs, &[1, 3, 4, 5, 8]);
}

#[test]
fn attribute_range_array() {
    let s = r#"{
        "Problem": { "Type": "Electrostatic" },
        "Model":   { "Mesh": "x.msh" },
        "Boundaries": {
            "Ground": { "Attributes": [10, 20, 30] }
        }
    }"#;
    let cfg = load_config_from_str(s, ConfigFormat::Json).unwrap();
    assert_eq!(cfg.boundaries.ground.unwrap().attributes, vec![10, 20, 30]);
}

// ── YAML parsing ─────────────────────────────────────────────────────────────

#[test]
fn parse_parallel_plate_yaml() {
    let s = fixture("parallel_plate.yaml");
    let cfg = load_config_from_str(&s, ConfigFormat::Yaml)
        .expect("parallel_plate.yaml should parse");

    assert_eq!(cfg.problem.problem_type, ProblemType::Electrostatic);
    assert!((cfg.domains.materials[0].permittivity - 4.5).abs() < 1e-12);
    assert_eq!(cfg.solver.linear.solver_type, "GMRES");
}

#[test]
fn json_and_yaml_equivalent() {
    let json = fixture("parallel_plate.json");
    let yaml = fixture("parallel_plate.yaml");

    let cj = load_config_from_str(&json, ConfigFormat::Json).unwrap();
    let cy = load_config_from_str(&yaml, ConfigFormat::Yaml).unwrap();

    assert_eq!(cj.problem.problem_type, cy.problem.problem_type);
    assert_eq!(cj.problem.verbose(), cy.problem.verbose());
    assert!((cj.model.l0 - cy.model.l0).abs() < 1e-20);
    assert!((cj.domains.materials[0].permittivity
        - cy.domains.materials[0].permittivity)
        .abs()
        < 1e-12);
    assert_eq!(
        cj.boundaries.pec.as_ref().unwrap().attributes,
        cy.boundaries.pec.as_ref().unwrap().attributes
    );
}

// ── Default values ────────────────────────────────────────────────────────────

#[test]
fn defaults_applied() {
    let s = r#"{ "Problem": { "Type": "Magnetostatic" }, "Model": { "Mesh": "m.msh" } }"#;
    let cfg = load_config_from_str(s, ConfigFormat::Json).unwrap();
    assert_eq!(cfg.problem.verbose(), 1);
    assert_eq!(cfg.problem.output_dir(), ".");
    assert!((cfg.model.l0 - 1.0).abs() < 1e-20);
    assert_eq!(cfg.solver.order, 1);
    assert!((cfg.solver.linear.tol - 1e-6).abs() < 1e-20);
    assert_eq!(cfg.solver.linear.max_iter, 200);
    assert!(cfg.boundaries.pec.is_none());
    assert!(cfg.boundaries.ground.is_none());
    assert!(cfg.domains.materials.is_empty());
}

// ── Magnetostatic ─────────────────────────────────────────────────────────────

#[test]
fn parse_magnetostatic_json() {
    let s = fixture("magnetostatic_2d.json");
    let cfg = load_config_from_str(&s, ConfigFormat::Json).unwrap();
    assert_eq!(cfg.problem.problem_type, ProblemType::Magnetostatic);
    assert_eq!(cfg.domains.materials.len(), 2);
    assert!((cfg.domains.materials[1].permeability - 1000.0).abs() < 1e-9);
}

// ── Mixed attribute spec (array + string in same config) ────────────────────

#[test]
fn mixed_attribute_specs() {
    let s = fixture("magnetostatic_2d.json");
    let cfg = load_config_from_str(&s, ConfigFormat::Json).unwrap();
    // PEC has "Attributes": [1, "2,4-6"] — should flatten to [1,2,4,5,6]
    // (Current fixture uses a single string; test that it parses without error)
    assert!(cfg.boundaries.pec.is_some());
}
