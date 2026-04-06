use rem_mesh::gen::{annular_msh, rect_bimaterial_msh};

#[derive(Clone, Copy, PartialEq)]
pub enum ExampleStatus {
    Ready,
    Unimplemented,
}

pub struct ExampleMeta {
    pub key: &'static str,
    pub label: &'static str,
    pub problem_type: &'static str,
    pub status: ExampleStatus,
    pub config_json: &'static str,
    pub source_code: &'static str,
}

pub const EXAMPLES: &[ExampleMeta] = &[
    ExampleMeta {
        key: "spheres",
        label: "Spheres (Electrostatic)",
        problem_type: "Electrostatic",
        status: ExampleStatus::Ready,
        config_json: r#"{
  "Problem": { "Type": "Electrostatic", "Output": "." },
  "Model": { "Mesh": "spheres.msh", "L0": 0.001 },
  "Domains": {
    "Materials": [{ "Attributes": [10], "Permittivity": 1.0, "Permeability": 1.0 }]
  },
  "Boundaries": {
    "Ground": { "Attributes": [2] },
    "Terminal": [{ "Index": 1, "Attributes": [1] }]
  },
  "Solver": { "Order": 1, "Linear": { "Type": "GMRES", "Tol": 1e-10, "MaxIter": 500 } }
}"#,
        source_code: "// Simplified view of palace_spheres.rs\n#[test]\nfn solve_spheres() {\n    let mesh = annular_msh(1.0, 4.0, 10, 32, 1, 2, 10);\n    let phi = solve_one(&cfg, &mesh, &dm, Some(1), 1.0, &comm).unwrap();\n    // ... verification\n}",
    },
    ExampleMeta {
        key: "rings",
        label: "Rings (Magnetostatic)",
        problem_type: "Magnetostatic",
        status: ExampleStatus::Ready,
        config_json: r#"{
  "Problem": { "Type": "Magnetostatic", "Output": "." },
  "Model": { "Mesh": "rings.msh", "L0": 0.001 },
  "Domains": {
    "Materials": [
      { "Attributes": [10], "Permittivity": 1.0, "Permeability": 1000.0 },
      { "Attributes": [20], "Permittivity": 1.0, "Permeability": 1.0 }
    ]
  },
  "Boundaries": {
    "Ground": { "Attributes": [1] },
    "SurfaceCurrent": [{ "Index": 1, "Attributes": [2], "Direction": "+Y" }]
  },
  "Solver": { "Order": 1, "Linear": { "Type": "GMRES", "Tol": 1e-10, "MaxIter": 500 } }
}"#,
        source_code: "// Simplified view of palace_rings.rs\n#[test]\nfn solve_rings() {\n    let mesh = rect_bimaterial_msh(1.0, 1.0, 20, 20, 1, 2, 10, 20);\n    let az = solve_one(&cfg, &mesh, &dm, Some(1), &comm).unwrap();\n    // ... verification\n}",
    },
    ExampleMeta {
        key: "adapter",
        label: "Adapter (Driven, v0.2)",
        problem_type: "Driven",
      status: ExampleStatus::Ready,
        config_json: r#"{
  "Problem": { "Type": "Driven", "Output": "." },
  "Model": { "Mesh": "adapter.msh", "L0": 1.0 },
  "Domains": { "Materials": [{ "Attributes": [1], "Permittivity": 1.0, "Permeability": 1.0 }] },
  "Boundaries": {
    "PEC": { "Attributes": [2] },
    "LumpedPort": [{ "Index": 1, "Attributes": [3], "R": 50.0 }]
  },
  "Solver": { "Order": 1, "Driven": { "MinFreq": 1e9, "MaxFreq": 10e9, "FreqStep": 0.1e9 } }
}"#,
        source_code: "// Driven frequency sweep example (wasm-enabled)",
    },
    ExampleMeta {
        key: "antenna",
        label: "Antenna (Driven, v0.2)",
        problem_type: "Driven",
      status: ExampleStatus::Ready,
        config_json: r#"{
  "Problem": { "Type": "Driven", "Output": "." },
  "Model": { "Mesh": "antenna.msh", "L0": 0.001 },
  "Domains": { "Materials": [{ "Attributes": [1], "Permittivity": 1.0, "Permeability": 1.0 }] },
  "Boundaries": {
    "Absorbing": { "Attributes": [2] },
    "LumpedPort": [{ "Index": 1, "Attributes": [3], "R": 50.0 }]
  },
  "Solver": { "Order": 2, "Driven": { "MinFreq": 2e9, "MaxFreq": 3e9, "FreqStep": 0.05e9 } }
}"#,
        source_code: "// Driven frequency sweep example (wasm-enabled)",
    },
    ExampleMeta {
        key: "coaxial",
        label: "Coaxial (Electrostatic)",
        problem_type: "Electrostatic",
        status: ExampleStatus::Ready,
        config_json: r#"{
  "Problem": { "Type": "Electrostatic", "Output": "." },
  "Model": { "Mesh": "coaxial.msh", "L0": 0.001 },
  "Domains": { "Materials": [{ "Attributes": [1], "Permittivity": 2.1 }] },
  "Boundaries": {
    "Ground": { "Attributes": [2] },
    "Terminal": [{ "Index": 1, "Attributes": [3] }]
  },
  "Solver": { "Order": 1 }
}"#,
        source_code: "// Electrostatic coaxial example",
    },
    ExampleMeta {
        key: "cpw",
        label: "CPW (Driven, v0.2)",
        problem_type: "Driven",
      status: ExampleStatus::Ready,
        config_json: r#"{
  "Problem": { "Type": "Driven", "Output": "." },
  "Model": { "Mesh": "cpw.msh", "L0": 1e-6 },
  "Domains": { "Materials": [{ "Attributes": [1], "Permittivity": 11.7 }] },
  "Boundaries": {
    "PEC": { "Attributes": [2] },
    "LumpedPort": [{ "Index": 1, "Attributes": [3], "R": 50.0 }]
  },
  "Solver": { "Order": 2, "Driven": { "MinFreq": 4e9, "MaxFreq": 8e9, "FreqStep": 0.1e9 } }
}"#,
        source_code: "// Driven frequency sweep example (wasm-enabled)",
    },
    ExampleMeta {
        key: "cylinder",
      label: "Cylinder (Magnetostatic)",
        problem_type: "Magnetostatic",
      status: ExampleStatus::Ready,
        config_json: r#"{
  "Problem": { "Type": "Magnetostatic", "Output": "." },
  "Model": { "Mesh": "cylinder.msh", "L0": 0.001 },
  "Domains": { "Materials": [{ "Attributes": [1], "Permeability": 1000.0 }] },
  "Boundaries": {
    "Ground": { "Attributes": [4] },
    "SurfaceCurrent": [{ "Index": 1, "Attributes": [2], "Direction": "+Z" }]
  },
  "Solver": { "Order": 1 }
}"#,
        source_code: "// Magnetostatic cylinder example\n// Physical groups: 1=cylinder(vol), 2=top, 3=bottom, 4=exterior, 5=symmetry",
    },
    ExampleMeta {
        key: "parallel_plate",
        label: "Parallel Plate (Electrostatic)",
        problem_type: "Electrostatic",
        status: ExampleStatus::Ready,
        config_json: r#"{
  "Problem": { "Type": "Electrostatic", "Output": "." },
  "Model": { "Mesh": "plate_2d.msh", "L0": 0.001 },
  "Domains": {
    "Materials": [{ "Attributes": [10], "Permittivity": 1.0 }]
  },
  "Boundaries": {
    "Ground": { "Attributes": [1] },
    "Terminal": [{ "Index": 1, "Attributes": [2] }]
  },
  "Solver": { "Order": 1, "Linear": { "Type": "GMRES", "Tol": 1e-10, "MaxIter": 500 } }
}"#,
        source_code: "// Electrostatic parallel plate example",
    },
    ExampleMeta {
        key: "sbr_sphere",
        label: "Sphere (SBR+)",
        problem_type: "SBR",
        status: ExampleStatus::Ready,
        config_json: r#"{
  "Problem": { "Type": "SBR", "Output": "." },
  "Model": { "Mesh": "sphere.msh", "L0": 1.0 },
  "Boundaries": { "PEC": { "Attributes": [1] } },
  "Solver": {
    "SBR": {
      "FreqMin": 1.0e9,
      "FreqMax": 1.0e9,
      "FreqStep": 0.0,
      "RayDensity": 3000.0,
      "MaxBounces": 2,
      "WeightThresh": 1.0e-4,
      "TargetType": "PEC",
      "ThetaInc": 0.0,
      "PhiInc": 0.0,
      "Polarization": "theta"
    }
  },
  "Postprocessing": {
    "RCS": {
      "ThetaDeg": [180.0],
      "PhiDeg": [0.0]
    }
  }
}"#,
        source_code: "// SBR+ PEC sphere example (Mie-comparison setup)",
    },
    ExampleMeta {
        key: "mom_sphere",
        label: "Sphere (MoM)",
        problem_type: "MoM",
        status: ExampleStatus::Ready,
        config_json: r#"{
  "Problem": { "Type": "MoM", "Output": "." },
  "Model": { "Mesh": "sphere.msh", "L0": 1.0 },
  "Boundaries": { "PEC": { "Attributes": [1] } },
  "Solver": {
    "MoM": {
      "Equation": "CFIE",
      "Basis": "RWG",
      "FreqMin": 1.0e9,
      "FreqMax": 1.0e9,
      "FreqStep": 1.0,
      "Alpha": 0.5,
      "SingularTol": 1.0e-6,
      "FastSolver": "Direct",
      "ThetaInc": 0.0,
      "PhiInc": 0.0,
      "Polarization": "theta"
    }
  }
}"#,
        source_code: "// MoM PEC sphere example with CFIE + RWG basis",
    },
    ExampleMeta {
        key: "transmon",
        label: "Transmon (Eigenmode, v0.2)",
        problem_type: "Eigenmode",
        status: ExampleStatus::Unimplemented,
        config_json: r#"{
  "Problem": { "Type": "Eigenmode", "Output": "." },
  "Model": { "Mesh": "transmon.msh", "L0": 1e-6 },
  "Domains": { "Materials": [{ "Attributes": [1], "Permittivity": 11.4 }] },
  "Boundaries": { "PEC": { "Attributes": [2] } },
  "Solver": { "Order": 2, "Eigenmode": { "N": 5, "Tol": 1e-8, "Target": 5e9 } }
}"#,
        source_code: "// Eigenmode solver not yet implemented",
    },
];

pub fn find_example(key: &str) -> Option<&'static ExampleMeta> {
    EXAMPLES.iter().find(|e| e.key == key)
}

pub fn get_mesh_bytes(key: &str) -> Vec<u8> {
    match key {
        "spheres" => annular_msh(1.0, 4.0, 10, 32, 1, 2, 10).into_bytes(),
        "rings" => rect_bimaterial_msh(1.0, 1.0, 20, 20, 1, 2, 10, 20).into_bytes(),
        "adapter" => include_bytes!("../../../examples/adapter/mesh/adapter.msh").to_vec(),
        "antenna" => include_bytes!("../../../examples/antenna/mesh/antenna.msh").to_vec(),
        "coaxial" => include_bytes!("../../../examples/coaxial/mesh/coaxial.msh").to_vec(),
        "cpw" => include_bytes!("../../../examples/cpw/mesh/cpw_coax.msh").to_vec(),
        "cylinder" => include_bytes!("../../../examples/cylinder/mesh/cylinder_hex.msh").to_vec(),
        "parallel_plate" => include_bytes!("../../../examples/parallel_plate/mesh/plate_2d.msh").to_vec(),
        "sbr_sphere" | "mom_sphere" => include_bytes!("../../../examples/sbr_sphere/mesh/sphere.msh").to_vec(),
        "transmon" => include_bytes!("../../../examples/transmon/mesh/transmon.msh2").to_vec(),
        _ => vec![],
    }
}
