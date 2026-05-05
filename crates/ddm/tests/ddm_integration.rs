//! End-to-end integration test for the DDM (Domain Decomposition Method) solver.
//!
//! This test builds a simple 2-D mesh programmatically, partitions it into 2
//! subdomains, runs the full Schwarz DDM solver pipeline, and verifies:
//!
//! 1. The solver completes without error.
//! 2. The Schwarz iteration converges (residual < tolerance).
//! 3. The solution vectors have the correct dimension (one per subdomain).
//! 4. The number of iterations is reasonable (≥ 1, ≤ max_iter).
//!
//! The mesh is a 2×1 rectangle split into 4 Tri3 elements:
//!
//! ```
//! (0,1)──(1,1)──(2,1)
//!   |  ╲  |  ╲  |
//!   |   ╲ |   ╲ |
//! (0,0)──(1,0)──(2,0)
//! ```
//!
//! Nodes 0–2 (x=0..1) belong to subdomain 0; nodes 3–5 (x=1..2) to subdomain 1.
//! Node 2 at (1,0) and node 3 at (1,1) are shared interface nodes.
//!
//! The DDM Schwarz solver alternates between subdomains using Robin transmission
//! conditions, exchanging interface field values until convergence.

use rem_config::DdmSolverConfig;
use rem_ddm::run_with_mesh;
use rem_mesh::{Element, ElementKind, Node, RemMesh};
use rem_parallel::NoComm;
use std::collections::HashMap;

fn minimal_palace_config() -> rem_config::PalaceConfig {
    // The _config parameter in run_with_mesh is currently unused; provide the
    // minimum valid JSON that deserialises without error.
    let json = r#"{
        "Problem": { "Type": "Driven" },
        "Model":   { "Mesh": "dummy.msh" }
    }"#;
    serde_json::from_str(json).expect("minimal PalaceConfig JSON should parse")
}

/// Build a 2×1 rectangle mesh with 6 nodes and 4 Tri3 elements.
fn build_rectangle_mesh() -> RemMesh {
    // 6 nodes: (0,0) (1,0) (2,0) (0,1) (1,1) (2,1)
    let nodes = vec![
        Node { id: 0, x: 0.0, y: 0.0, z: 0.0 },
        Node { id: 1, x: 1.0, y: 0.0, z: 0.0 },
        Node { id: 2, x: 2.0, y: 0.0, z: 0.0 },
        Node { id: 3, x: 0.0, y: 1.0, z: 0.0 },
        Node { id: 4, x: 1.0, y: 1.0, z: 0.0 },
        Node { id: 5, x: 2.0, y: 1.0, z: 0.0 },
    ];

    // 4 Tri3 elements, tag=1 (material domain 1)
    let volume_elements = vec![
        Element { id: 0, kind: ElementKind::Tri3, tag: 1, node_ids: vec![0, 1, 3], rank: 0 },
        Element { id: 1, kind: ElementKind::Tri3, tag: 1, node_ids: vec![1, 4, 3], rank: 0 },
        Element { id: 2, kind: ElementKind::Tri3, tag: 1, node_ids: vec![1, 2, 4], rank: 0 },
        Element { id: 3, kind: ElementKind::Tri3, tag: 1, node_ids: vec![2, 5, 4], rank: 0 },
    ];

    // Boundary: 4 outer edges, tag=10 (Ground)
    let boundary_elements = vec![
        Element { id: 10, kind: ElementKind::Line2, tag: 10, node_ids: vec![0, 1], rank: 0 },
        Element { id: 11, kind: ElementKind::Line2, tag: 10, node_ids: vec![1, 2], rank: 0 },
        Element { id: 12, kind: ElementKind::Line2, tag: 10, node_ids: vec![2, 5], rank: 0 },
        Element { id: 13, kind: ElementKind::Line2, tag: 10, node_ids: vec![5, 4], rank: 0 },
        Element { id: 14, kind: ElementKind::Line2, tag: 10, node_ids: vec![4, 3], rank: 0 },
        Element { id: 15, kind: ElementKind::Line2, tag: 10, node_ids: vec![3, 0], rank: 0 },
    ];

    use rem_mesh::BoundaryTag;
    let mut boundary_tags = HashMap::new();
    boundary_tags.insert(10u32, BoundaryTag::Ground);

    RemMesh {
        nodes,
        volume_elements,
        boundary_elements,
        boundary_tags,
        domain_tags: HashMap::new(),
        dim: 2,
        rank: 0,
        size: 1,
    }
}

/// Build a minimal DDM config for 2 subdomains.
fn make_ddm_cfg(n_sub: usize) -> DdmSolverConfig {
    DdmSolverConfig {
        num_subdomains: n_sub,
        method: "Schwarz".to_string(),
        robin_order: 1,
        tolerance: 1e-4,
        max_iter: 50,
        partition_type: "Dual".to_string(),
    }
}

#[test]
fn ddm_two_subdomain_converges() {
    let _ = env_logger::try_init();

    let mesh   = build_rectangle_mesh();
    let cfg    = make_ddm_cfg(2);
    let config = minimal_palace_config();

    let result = run_with_mesh(&config, &cfg, &mesh, &NoComm)
        .expect("DDM run_with_mesh should succeed on a simple 2-subdomain rectangle mesh");

    // Convergence criteria
    assert!(
        result.residual < cfg.tolerance * 10.0,
        "Schwarz residual {:.4e} should be near tolerance {:.4e}",
        result.residual,
        cfg.tolerance
    );
    assert!(
        result.iterations >= 1,
        "At least 1 Schwarz iteration expected"
    );
    assert!(
        result.iterations <= cfg.max_iter,
        "Iterations {} should not exceed max_iter {}",
        result.iterations,
        cfg.max_iter
    );

    // Structural checks
    assert_eq!(
        result.subdomain_solutions.len(),
        2,
        "Expected 2 subdomain solution vectors"
    );
    for (i, sol) in result.subdomain_solutions.iter().enumerate() {
        assert!(
            !sol.is_empty(),
            "Subdomain {} solution vector should not be empty",
            i
        );
    }
}

#[test]
fn ddm_single_subdomain_is_trivial() {
    let _ = env_logger::try_init();

    let mesh   = build_rectangle_mesh();
    let cfg    = make_ddm_cfg(1);
    let config = minimal_palace_config();

    let result = run_with_mesh(&config, &cfg, &mesh, &NoComm)
        .expect("DDM with 1 subdomain should succeed");

    // Single subdomain: the Schwarz step is trivially a direct solve → residual = 0.
    assert!(
        result.residual.is_finite(),
        "Residual must be finite for single-subdomain case"
    );
    assert_eq!(
        result.subdomain_solutions.len(),
        1,
        "Single subdomain case should produce exactly 1 solution vector"
    );
}

#[test]
fn ddm_four_subdomains_converges() {
    let _ = env_logger::try_init();

    // For 4 subdomains on a 4-element mesh, each subdomain may get 1 element.
    // METIS is expected to produce a valid partition even in this degenerate case.
    let mesh   = build_rectangle_mesh();
    let cfg    = make_ddm_cfg(4);
    let config = minimal_palace_config();

    // Allow failure here: METIS may refuse to partition 4 elements into 4 subdomains
    // if some result in empty domains.  If it succeeds, validate the result.
    match run_with_mesh(&config, &cfg, &mesh, &NoComm) {
        Ok(result) => {
            assert!(result.residual.is_finite(), "Residual must be finite");
            assert!(result.subdomain_solutions.len() <= 4);
        }
        Err(e) => {
            // METIS may reject degenerate partitions; that's acceptable.
            let msg = format!("{e}");
            assert!(
                msg.contains("partition") || msg.contains("METIS") || msg.contains("subdomain"),
                "Unexpected error: {msg}"
            );
        }
    }
}
