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
    // ── REM 独立示例 ──────────────────────────────────────────────────────────

    ExampleMeta {
        key: "spheres",
        label: "Spheres (Electrostatic)",
        problem_type: "Electrostatic",
        status: ExampleStatus::Ready,
        config_json: r#"{
  "Problem": { "Type": "Electrostatic", "Verbose": 1, "Output": "./output/spheres" },
  "Model": { "Mesh": "spheres.msh", "L0": 1.0e-3 },
  "Domains": {
    "Materials": [{ "Attributes": [10], "Permittivity": 1.0, "Permeability": 1.0 }]
  },
  "Boundaries": {
    "Ground": { "Attributes": [2] },
    "Terminal": [{ "Index": 1, "Attributes": [1] }]
  },
  "Solver": { "Order": 1, "Linear": { "Type": "GMRES", "Tol": 1.0e-10, "MaxIter": 500 } }
}"#,
        source_code: "// Simplified view of palace_spheres.rs\n#[test]\nfn solve_spheres() {\n    let mesh = annular_msh(1.0, 4.0, 10, 32, 1, 2, 10);\n    let phi = solve_one(&cfg, &mesh, &dm, Some(1), 1.0, &comm).unwrap();\n    // ... verification\n}",
    },

    ExampleMeta {
        key: "rings",
        label: "Rings (Magnetostatic)",
        problem_type: "Magnetostatic",
        status: ExampleStatus::Ready,
        config_json: r#"{
  "Problem": { "Type": "Magnetostatic", "Verbose": 1, "Output": "./output/rings" },
  "Model": { "Mesh": "rings.msh", "L0": 1.0e-3 },
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
  "Solver": { "Order": 1, "Linear": { "Type": "GMRES", "Tol": 1.0e-10, "MaxIter": 500 } }
}"#,
        source_code: "// Simplified view of palace_rings.rs\n#[test]\nfn solve_rings() {\n    let mesh = rect_bimaterial_msh(1.0, 1.0, 20, 20, 1, 2, 10, 20);\n    let az = solve_one(&cfg, &mesh, &dm, Some(1), &comm).unwrap();\n    // ... verification\n}",
    },

    ExampleMeta {
        key: "parallel_plate",
        label: "Parallel Plate (Electrostatic)",
        problem_type: "Electrostatic",
        status: ExampleStatus::Ready,
        config_json: r#"{
  "Problem": { "Type": "Electrostatic", "Verbose": 1, "Output": "./output/parallel_plate" },
  "Model": { "Mesh": "plate_2d.msh", "L0": 1.0e-3 },
  "Domains": {
    "Materials": [{ "Attributes": [10], "Permittivity": 1.0 }]
  },
  "Boundaries": {
    "Ground": { "Attributes": [1] },
    "Terminal": [{ "Index": 1, "Attributes": [2] }]
  },
  "Solver": { "Order": 1, "Linear": { "Type": "GMRES", "Tol": 1.0e-10, "MaxIter": 500 } }
}"#,
        source_code: "// Electrostatic parallel plate example",
    },

    ExampleMeta {
        key: "rem_es_fast",
        label: "REM Quick Check / Parallel Plate (Electrostatic)",
        problem_type: "Electrostatic",
        status: ExampleStatus::Ready,
        config_json: r#"{
  "Problem": { "Type": "Electrostatic", "Verbose": 1, "Output": "./output/rem_es_parallel_plate_fast" },
  "Model": { "Mesh": "plate_2d.msh", "L0": 1.0e-3 },
  "Domains": {
    "Materials": [{ "Attributes": [10], "Permittivity": 1.0 }]
  },
  "Boundaries": {
    "Ground": { "Attributes": [1] },
    "Terminal": [{ "Index": 1, "Attributes": [2] }]
  },
  "Solver": { "Order": 1, "Linear": { "Type": "GMRES", "Tol": 1.0e-10, "MaxIter": 200 } }
}"#,
        source_code: "// Fast electrostatic sanity check: small mesh + closed-form C = eps0*A/d",
    },

    ExampleMeta {
        key: "rem_ms_fast",
        label: "REM Quick Check / Surface Current Strip (Magnetostatic)",
        problem_type: "Magnetostatic",
        status: ExampleStatus::Ready,
        config_json: r#"{
  "Problem": { "Type": "Magnetostatic", "Verbose": 1, "Output": "./output/rem_ms_parallel_plate_fast" },
  "Model": { "Mesh": "plate_2d.msh", "L0": 1.0e-3 },
  "Domains": {
    "Materials": [{ "Attributes": [10], "Permeability": 1.0 }]
  },
  "Boundaries": {
    "Ground": { "Attributes": [1] },
    "SurfaceCurrent": [{ "Index": 1, "Attributes": [2], "Direction": "+Y" }]
  },
  "Solver": { "Order": 1, "Linear": { "Type": "GMRES", "Tol": 1.0e-10, "MaxIter": 200 } }
}"#,
        source_code: "// Fast magnetostatic sanity check: tiny mesh + one current boundary",
    },

    ExampleMeta {
        key: "rem_driven_fast",
        label: "REM Quick Check / CPW Single-Port (Driven)",
        problem_type: "Driven",
        status: ExampleStatus::Ready,
        config_json: r#"{
  "Problem": { "Type": "Driven", "Verbose": 1, "Output": "./output/rem_driven_cpw_fast" },
  "Model": { "Mesh": "cpw_coax.msh", "L0": 1.0e-6 },
  "Domains": { "Materials": [{ "Attributes": [1], "Permittivity": 11.7 }] },
  "Boundaries": {
    "PEC": { "Attributes": [2] },
    "LumpedPort": [{ "Index": 1, "Attributes": [3], "R": 50.0, "Excitation": true }]
  },
  "Solver": {
    "Order": 1,
    "Driven": { "MinFreq": 6.0e9, "MaxFreq": 6.2e9, "FreqStep": 0.2e9, "SaveStep": 1 }
  }
}"#,
        source_code: "// Fast driven sanity check: single lumped port + 2-point sweep",
    },

    ExampleMeta {
        key: "rem_eigen_fast",
        label: "REM Quick Check / Cylinder (Eigenmode)",
        problem_type: "Eigenmode",
        status: ExampleStatus::Ready,
        config_json: r#"{
  "Problem": { "Type": "Eigenmode", "Verbose": 1, "Output": "./output/rem_eigen_cylinder_fast" },
  "Model": { "Mesh": "cylinder_tet.msh", "L0": 1.0e-2 },
  "Domains": {
    "Materials": [{ "Attributes": [1], "Permeability": 1.0, "Permittivity": 2.08, "LossTan": 4.0e-4 }]
  },
  "Boundaries": { "PEC": { "Attributes": [4] } },
  "Solver": {
    "Order": 2,
    "Eigenmode": { "N": 2, "Tol": 1.0e-8, "Target": 2.0e9, "Save": 1 },
    "Linear": { "Tol": 1.0e-8, "MaxIter": 80 }
  }
}"#,
        source_code: "// Fast eigenmode sanity check: 2 modes on tet cylinder mesh",
    },

    ExampleMeta {
        key: "rem_transient_fast",
        label: "REM Quick Check / Coaxial Pulse (Transient)",
        problem_type: "Transient",
        status: ExampleStatus::Ready,
        config_json: r#"{
  "Problem": { "Type": "Transient", "Verbose": 1, "Output": "./output/rem_transient_coax_fast" },
  "Model": { "Mesh": "coaxial.msh", "L0": 1.0e-3 },
  "Domains": {
    "Materials": [{ "Attributes": [1], "Permeability": 1.0, "Permittivity": 2.08, "Conductivity": 4.629e-2 }]
  },
  "Boundaries": {
    "PEC": { "Attributes": [2] },
    "LumpedPort": [
      { "Index": 1, "Attributes": [3], "R": 50.0, "Direction": "+R", "Excitation": true },
      { "Index": 2, "Attributes": [4], "R": 50.0, "Direction": "+R" }
    ]
  },
  "Solver": {
    "Order": 2,
    "Transient": {
      "Type": "GeneralizedAlpha", "Excitation": "GaussianPulse",
      "ExcitationWidth": 0.03, "MaxTime": 0.20, "TimeStep": 0.01, "SaveStep": 5
    },
    "Linear": { "Tol": 1.0e-8, "MaxIter": 80 }
  }
}"#,
        source_code: "// Fast transient sanity check: short Gaussian pulse in matched coax",
    },

    ExampleMeta {
        key: "rem_mom_fast",
        label: "REM Quick Check / PEC Sphere (MoM)",
        problem_type: "MoM",
        status: ExampleStatus::Ready,
        config_json: r#"{
  "Problem": { "Type": "MoM", "Verbose": 1, "Output": "./output/rem_mom_sphere_fast" },
  "Model": { "Mesh": "sphere.msh", "L0": 1.0 },
  "Boundaries": { "PEC": { "Attributes": [1] } },
  "Solver": {
    "MoM": {
      "Equation": "EFIE", "Basis": "Pulse",
      "FreqMin": 1.0e9, "FreqMax": 1.0e9, "FreqStep": 1.0,
      "Alpha": 0.5, "SingularTol": 5.0e-4, "FastSolver": "GMRES",
      "ThetaInc": 0.0, "PhiInc": 0.0, "Polarization": "theta"
    }
  },
  "Postprocessing": { "RCS": { "ThetaDeg": "0:30:180", "PhiDeg": [0.0] } }
}"#,
        source_code: "// Fast MoM sanity check: single-frequency PEC sphere",
    },

    ExampleMeta {
        key: "rem_sbr_fast",
        label: "REM Quick Check / PEC Sphere (SBR)",
        problem_type: "SBR",
        status: ExampleStatus::Ready,
        config_json: r#"{
  "Problem": { "Type": "SBR", "Verbose": 1, "Output": "./output/rem_sbr_sphere_fast" },
  "Model": { "Mesh": "sphere.msh", "L0": 1.0 },
  "Boundaries": { "PEC": { "Attributes": [1] } },
  "Solver": {
    "SBR": {
      "FreqMin": 3.0e9, "FreqMax": 3.0e9, "FreqStep": 0.0,
      "RayDensity": 1200.0, "MaxBounces": 1, "WeightThresh": 1.0e-3,
      "TargetType": "PEC", "ThetaInc": 0.0, "PhiInc": 0.0, "Polarization": "theta"
    }
  },
  "Postprocessing": { "RCS": { "ThetaDeg": "0:10:180", "PhiDeg": [0.0] } }
}"#,
        source_code: "// Fast SBR sanity check: low ray density + 1 bounce",
    },

    ExampleMeta {
        key: "sbr_sphere",
        label: "Sphere (SBR+)",
        problem_type: "SBR",
        status: ExampleStatus::Ready,
        config_json: r#"{
  "Problem": { "Type": "SBR", "Verbose": 1, "Output": "output/sbr_sphere" },
  "Model": { "Mesh": "sphere.msh", "L0": 1.0 },
  "Boundaries": { "PEC": { "Attributes": [1] } },
  "Solver": {
    "SBR": {
      "FreqMin": 3.0e9, "FreqMax": 3.0e9, "FreqStep": 0.0,
      "RayDensity": 5000.0, "MaxBounces": 3, "WeightThresh": 1.0e-4,
      "TargetType": "PEC", "ThetaInc": 0.0, "PhiInc": 0.0, "Polarization": "theta"
    }
  },
  "Postprocessing": { "RCS": { "ThetaDeg": "0:5:180", "PhiDeg": [0.0] } }
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
      "Equation": "CFIE", "Basis": "RWG",
      "FreqMin": 1.0e9, "FreqMax": 1.0e9, "FreqStep": 1.0,
      "Alpha": 0.5, "SingularTol": 1.0e-6, "FastSolver": "Direct",
      "ThetaInc": 0.0, "PhiInc": 0.0, "Polarization": "theta"
    }
  }
}"#,
        source_code: "// MoM PEC sphere example with CFIE + RWG basis",
    },

    // ── Palace 对齐示例：adapter ──────────────────────────────────────────────

    ExampleMeta {
        key: "adapter",
        label: "Adapter / hybrid (Eigenmode)",
        problem_type: "Eigenmode",
        status: ExampleStatus::Ready,
        config_json: r#"{
  "Problem": { "Type": "Eigenmode", "Verbose": 2, "Output": "output/hybrid" },
  "Model": { "Mesh": "adapter.msh", "L0": 1.0e-2 },
  "Domains": { "Materials": [{ "Attributes": [1], "Permeability": 1.0, "Permittivity": 1.0, "LossTan": 0.0 }] },
  "Boundaries": {
    "WavePort": [
      { "Index": 1, "Attributes": [2], "Mode": 1 },
      { "Index": 2, "Attributes": [3], "Mode": 1 }
    ],
    "PEC": { "Attributes": [4] }
  },
  "Solver": {
    "Order": 2,
    "Eigenmode": { "N": 3, "Tol": 1.0e-6, "Target": 6.6e9, "Save": 3 },
    "Linear": { "Tol": 1.0e-8, "MaxIter": 100 }
  }
}"#,
        source_code: "// Eigenmode adapter example (hybrid.json)",
    },

    // ── Palace 对齐示例：antenna ─────────────────────────────────────────────

    ExampleMeta {
        key: "antenna_halfwave_dipole",
        label: "Antenna / halfwave dipole (Driven)",
        problem_type: "Driven",
        status: ExampleStatus::Ready,
        config_json: r#"{
  "Problem": { "Type": "Driven", "Verbose": 2, "Output": "postpro/antenna_halfwave_dipole" },
  "Model": { "Mesh": "antenna.msh", "L0": 1.0 },
  "Domains": { "Materials": [{ "Attributes": [7] }] },
  "Boundaries": {
    "Absorbing": { "Attributes": [4] },
    "PEC": { "Attributes": [1, 2] },
    "LumpedPort": [{ "Index": 1, "R": 50.0, "Excitation": true, "Attributes": [3], "Direction": "+Z" }]
  },
  "Solver": {
    "Order": 2,
    "Driven": { "Samples": [{ "Type": "Point", "Freq": [0.0749], "SaveStep": 1 }] },
    "Linear": { "Type": "Default", "KSPType": "GMRES", "Tol": 1.0e-10, "MaxIts": 100 }
  }
}"#,
        source_code: "// Halfwave dipole antenna: LumpedPort feed at gap, FarField postprocessing",
    },

    ExampleMeta {
        key: "antenna_short_dipole",
        label: "Antenna / short dipole (Driven)",
        problem_type: "Driven",
        status: ExampleStatus::Ready,
        config_json: r#"{
  "Problem": { "Type": "Driven", "Verbose": 2, "Output": "postpro/antenna_short_dipole" },
  "Model": { "Mesh": "antenna.msh", "L0": 1.0 },
  "Domains": {
    "Materials": [{ "Attributes": [5, 6, 7] }],
    "CurrentDipole": [{ "Index": 1, "Moment": 1.0, "Center": [0.0, 0.0, 0.0], "Direction": [0, 0, 1] }]
  },
  "Boundaries": {
    "Absorbing": { "Attributes": [4], "Order": 2 }
  },
  "Solver": {
    "Order": 2,
    "Driven": { "Samples": [{ "Type": "Point", "Freq": [0.0749], "SaveStep": 1 }] },
    "Linear": { "Type": "Default", "KSPType": "GMRES", "Tol": 1.0e-10, "MaxIts": 100 }
  }
}"#,
        source_code: "// Short dipole: CurrentDipole source (no physical port), Absorbing order 2",
    },

    // ── Palace 对齐示例：coaxial ─────────────────────────────────────────────

    ExampleMeta {
        key: "coaxial",
        label: "Coaxial (Electrostatic)",
        problem_type: "Electrostatic",
        status: ExampleStatus::Ready,
        config_json: r#"{
  "Problem": { "Type": "Electrostatic", "Verbose": 1, "Output": "./output/coaxial" },
  "Model": { "Mesh": "coaxial.msh", "L0": 1.0e-3 },
  "Domains": { "Materials": [{ "Attributes": [10], "Permittivity": 2.1 }] },
  "Boundaries": {
    "Ground": { "Attributes": [2] },
    "Terminal": [{ "Index": 1, "Attributes": [1] }]
  },
  "Solver": { "Order": 1 }
}"#,
        source_code: "// Electrostatic coaxial example",
    },

    ExampleMeta {
        key: "coaxial_matched",
        label: "Coaxial / matched (Transient)",
        problem_type: "Transient",
        status: ExampleStatus::Ready,
        config_json: r#"{
  "Problem": { "Type": "Transient", "Verbose": 2, "Output": "output/coaxial_matched" },
  "Model": { "Mesh": "coaxial.msh", "L0": 1.0e-3 },
  "Domains": { "Materials": [{ "Attributes": [1], "Permeability": 1.0, "Permittivity": 2.08, "Conductivity": 4.629e-2 }] },
  "Boundaries": {
    "PEC": { "Attributes": [2] },
    "LumpedPort": [
      { "Index": 1, "Attributes": [3], "R": 50.0, "Direction": "+R", "Excitation": true },
      { "Index": 2, "Attributes": [4], "R": 50.0, "Direction": "+R" }
    ]
  },
  "Solver": {
    "Order": 3,
    "Transient": {
      "Type": "GeneralizedAlpha", "Excitation": "ModulatedGaussian",
      "ExcitationFreq": 10.0, "ExcitationWidth": 0.05,
      "MaxTime": 1.0, "TimeStep": 0.005, "SaveStep": 10
    },
    "Linear": { "Tol": 1.0e-8, "MaxIter": 100 }
  }
}"#,
        source_code: "// Coaxial matched: ModulatedGaussian excitation, matched 50Ω load at far end",
    },

    ExampleMeta {
        key: "coaxial_open",
        label: "Coaxial / open (Transient)",
        problem_type: "Transient",
        status: ExampleStatus::Ready,
        config_json: r#"{
  "Problem": { "Type": "Transient", "Verbose": 2, "Output": "output/coaxial_open" },
  "Model": { "Mesh": "coaxial.msh", "L0": 1.0e-3 },
  "Domains": { "Materials": [{ "Attributes": [1], "Permeability": 1.0, "Permittivity": 2.08, "Conductivity": 4.629e-2 }] },
  "Boundaries": {
    "PEC": { "Attributes": [2] },
    "PMC": { "Attributes": [4] },
    "LumpedPort": [{ "Index": 1, "Attributes": [3], "R": 50.0, "Direction": "+R", "Excitation": true }]
  },
  "Solver": {
    "Order": 3,
    "Transient": {
      "Type": "GeneralizedAlpha", "Excitation": "ModulatedGaussian",
      "ExcitationFreq": 10.0, "ExcitationWidth": 0.05,
      "MaxTime": 1.0, "TimeStep": 0.005, "SaveStep": 10
    },
    "Linear": { "Tol": 1.0e-8, "MaxIter": 100 }
  }
}"#,
        source_code: "// Coaxial open: PMC at far end �?full reflection",
    },

    ExampleMeta {
        key: "coaxial_short",
        label: "Coaxial / short (Transient)",
        problem_type: "Transient",
        status: ExampleStatus::Ready,
        config_json: r#"{
  "Problem": { "Type": "Transient", "Verbose": 2, "Output": "output/coaxial_short" },
  "Model": { "Mesh": "coaxial.msh", "L0": 1.0e-3 },
  "Domains": { "Materials": [{ "Attributes": [1], "Permeability": 1.0, "Permittivity": 2.08, "Conductivity": 4.629e-2 }] },
  "Boundaries": {
    "PEC": { "Attributes": [2, 4] },
    "LumpedPort": [{ "Index": 1, "Attributes": [3], "R": 50.0, "Direction": "+R", "Excitation": true }]
  },
  "Solver": {
    "Order": 3,
    "Transient": {
      "Type": "GeneralizedAlpha", "Excitation": "ModulatedGaussian",
      "ExcitationFreq": 10.0, "ExcitationWidth": 0.05,
      "MaxTime": 1.0, "TimeStep": 0.005, "SaveStep": 10
    },
    "Linear": { "Tol": 1.0e-8, "MaxIter": 100 }
  }
}"#,
        source_code: "// Coaxial short: PEC at far end �?phase-reversed reflection",
    },

    // ── Palace 对齐示例：cpw ─────────────────────────────────────────────────

    ExampleMeta {
        key: "cpw_coax_adaptive",
        label: "CPW / coax adaptive (Driven)",
        problem_type: "Driven",
        status: ExampleStatus::Ready,
        config_json: r#"{
  "Problem": { "Type": "Driven", "Verbose": 2, "Output": "output/cpw_coax_adaptive" },
  "Model": { "Mesh": "cpw_coax_0.msh", "L0": 1.0e-6, "Refinement": {} },
  "Domains": {
    "Materials": [
      { "Attributes": [1], "Permeability": 1.0, "Permittivity": 1.0, "LossTan": 0.0 },
      { "Attributes": [2], "Permeability": 0.99999975, "Permittivity": 9.3, "LossTan": 3.0e-5 }
    ]
  },
  "Boundaries": {
    "PEC": { "Attributes": [8, 9, 11] },
    "Absorbing": { "Attributes": [10], "Order": 1 },
    "LumpedPort": [
      { "Index": 1, "Attributes": [4], "Direction": "+R", "R": 56.02, "Excitation": true },
      { "Index": 2, "Attributes": [5], "Direction": "+R", "R": 56.02 },
      { "Index": 3, "Attributes": [6], "Direction": "+R", "R": 56.02 },
      { "Index": 4, "Attributes": [7], "Direction": "+R", "R": 56.02 }
    ]
  },
  "Solver": {
    "Order": 2,
    "Driven": { "MinFreq": 2.0e9, "MaxFreq": 30.0e9, "FreqStep": 0.1e9, "SaveStep": 40, "AdaptiveTol": 1.0e-3 },
    "Linear": { "Tol": 1.0e-8, "MaxIter": 200 }
  }
}"#,
        source_code: "// CPW coax port, adaptive frequency sweep 2�?0 GHz, R=56.02Ω",
    },

    ExampleMeta {
        key: "cpw_coax_uniform",
        label: "CPW / coax uniform (Driven)",
        problem_type: "Driven",
        status: ExampleStatus::Ready,
        config_json: r#"{
  "Problem": { "Type": "Driven", "Verbose": 2, "Output": "output/cpw_coax_uniform" },
  "Model": { "Mesh": "cpw_coax_0.msh", "L0": 1.0e-6, "Refinement": {} },
  "Domains": {
    "Materials": [
      { "Attributes": [1], "Permeability": 1.0, "Permittivity": 1.0, "LossTan": 0.0 },
      { "Attributes": [2], "Permeability": 0.99999975, "Permittivity": 9.3, "LossTan": 3.0e-5 }
    ]
  },
  "Boundaries": {
    "PEC": { "Attributes": [8, 9, 11] },
    "Absorbing": { "Attributes": [10], "Order": 1 },
    "LumpedPort": [
      { "Index": 1, "Attributes": [4], "Direction": "+R", "R": 56.02, "Excitation": true },
      { "Index": 2, "Attributes": [5], "Direction": "+R", "R": 56.02 },
      { "Index": 3, "Attributes": [6], "Direction": "-R", "R": 56.02 },
      { "Index": 4, "Attributes": [7], "Direction": "-R", "R": 56.02 }
    ]
  },
  "Solver": {
    "Order": 2,
    "Driven": { "MinFreq": 2.0e9, "MaxFreq": 30.0e9, "FreqStep": 2.0e9, "SaveStep": 2 },
    "Linear": { "Tol": 1.0e-8, "MaxIter": 200 }
  }
}"#,
        source_code: "// CPW coax port, uniform frequency sweep 2�?0 GHz step 2 GHz",
    },

    ExampleMeta {
        key: "cpw_lumped_adaptive",
        label: "CPW / lumped adaptive (Driven)",
        problem_type: "Driven",
        status: ExampleStatus::Ready,
        config_json: r#"{
  "Problem": { "Type": "Driven", "Verbose": 2, "Output": "output/cpw_lumped_adaptive" },
  "Model": { "Mesh": "cpw_lumped_0.msh", "L0": 1.0e-6, "Refinement": {} },
  "Domains": {
    "Materials": [
      { "Attributes": [1], "Permeability": 1.0, "Permittivity": 1.0, "LossTan": 0.0 },
      { "Attributes": [2], "Permeability": 0.99999975, "Permittivity": 9.3, "LossTan": 3.0e-5 }
    ]
  },
  "Boundaries": {
    "PEC": { "Attributes": [13] },
    "Absorbing": { "Attributes": [4], "Order": 1 },
    "LumpedPort": [
      { "Index": 1, "Attributes": [5], "Direction": "+Y", "R": 56.02, "Excitation": true },
      { "Index": 2, "Attributes": [6], "Direction": "+Y", "R": 56.02 },
      { "Index": 3, "Attributes": [7], "Direction": "+Y", "R": 56.02 },
      { "Index": 4, "Attributes": [8], "Direction": "+Y", "R": 56.02 }
    ]
  },
  "Solver": {
    "Order": 2,
    "Driven": { "MinFreq": 2.0e9, "MaxFreq": 32.0e9, "FreqStep": 0.1e9, "SaveStep": 1, "AdaptiveTol": 1.0e-3 },
    "Linear": { "Tol": 1.0e-8, "MaxIter": 200 }
  }
}"#,
        source_code: "// CPW lumped port +Y, adaptive sweep 2�?2 GHz",
    },

    ExampleMeta {
        key: "cpw_lumped_uniform",
        label: "CPW / lumped uniform (Driven)",
        problem_type: "Driven",
        status: ExampleStatus::Ready,
        config_json: r#"{
  "Problem": { "Type": "Driven", "Verbose": 2, "Output": "output/cpw_lumped_uniform" },
  "Model": { "Mesh": "cpw_lumped_0.msh", "L0": 1.0e-6, "Refinement": {} },
  "Domains": {
    "Materials": [
      { "Attributes": [1], "Permeability": 1.0, "Permittivity": 1.0, "LossTan": 0.0 },
      { "Attributes": [2], "Permeability": 0.99999975, "Permittivity": 9.3, "LossTan": 3.0e-5 }
    ]
  },
  "Boundaries": {
    "PEC": { "Attributes": [13] },
    "Absorbing": { "Attributes": [4], "Order": 1 },
    "LumpedPort": [
      { "Index": 1, "Attributes": [5], "Direction": "+Y", "R": 56.02, "Excitation": true },
      { "Index": 2, "Attributes": [6], "Direction": "+Y", "R": 56.02 },
      { "Index": 3, "Attributes": [7], "Direction": "+Y", "R": 56.02 },
      { "Index": 4, "Attributes": [8], "Direction": "+Y", "R": 56.02 }
    ]
  },
  "Solver": {
    "Order": 2,
    "Driven": { "MinFreq": 2.0e9, "MaxFreq": 32.0e9, "FreqStep": 6.0e9, "SaveStep": 1 },
    "Linear": { "Tol": 1.0e-8, "MaxIter": 200 }
  }
}"#,
        source_code: "// CPW lumped port +Y, uniform sweep 2�?2 GHz step 6 GHz",
    },

    ExampleMeta {
        key: "cpw_lumped_eigen",
        label: "CPW / lumped (Eigenmode)",
        problem_type: "Eigenmode",
        status: ExampleStatus::Ready,
        config_json: r#"{
  "Problem": { "Type": "Eigenmode", "Verbose": 2, "Output": "output/cpw_lumped_eigen" },
  "Model": { "Mesh": "cpw_lumped_0.msh", "L0": 1.0e-6, "Refinement": {} },
  "Domains": {
    "Materials": [
      { "Attributes": [1], "Permeability": 1.0, "Permittivity": 1.0, "LossTan": 0.0 },
      { "Attributes": [2], "Permeability": 0.99999975, "Permittivity": 9.3, "LossTan": 3.0e-5 }
    ]
  },
  "Boundaries": {
    "PEC": { "Attributes": [13] },
    "Absorbing": { "Attributes": [4], "Order": 2 },
    "LumpedPort": [
      { "Index": 1, "Attributes": [5], "Direction": "+Y", "R": 56.02 },
      { "Index": 2, "Attributes": [6], "Direction": "+Y", "R": 56.02 },
      { "Index": 3, "Attributes": [7], "Direction": "+Y", "R": 56.02 },
      { "Index": 4, "Attributes": [8], "Direction": "+Y", "R": 56.02 }
    ]
  },
  "Solver": {
    "Order": 2,
    "Eigenmode": { "N": 1, "Tol": 1.0e-6, "Target": 16.0e9, "Save": 1 },
    "Linear": { "Tol": 1.0e-10, "MaxIter": 1000 }
  }
}"#,
        source_code: "// CPW lumped port, N=1 eigenmode near 16 GHz",
    },

    ExampleMeta {
        key: "cpw_wave_adaptive",
        label: "CPW / wave adaptive (Driven)",
        problem_type: "Driven",
        status: ExampleStatus::Ready,
        config_json: r#"{
  "Problem": { "Type": "Driven", "Verbose": 2, "Output": "output/cpw_wave_adaptive" },
  "Model": { "Mesh": "cpw_wave_0.msh", "L0": 1.0e-6, "Refinement": {} },
  "Domains": {
    "Materials": [
      { "Attributes": [1], "Permeability": 1.0, "Permittivity": 1.0, "LossTan": 0.0 },
      { "Attributes": [2], "Permeability": 0.99999975, "Permittivity": 9.3, "LossTan": 3.0e-5 }
    ]
  },
  "Boundaries": {
    "PEC": { "Attributes": [8, 9, 11] },
    "Absorbing": { "Attributes": [10], "Order": 1 },
    "WavePort": [
      { "Index": 1, "Attributes": [4], "Mode": 1, "Excitation": true },
      { "Index": 2, "Attributes": [5], "Mode": 1 },
      { "Index": 3, "Attributes": [6], "Mode": 1 },
      { "Index": 4, "Attributes": [7], "Mode": 1 }
    ]
  },
  "Solver": {
    "Order": 2,
    "Driven": { "MinFreq": 2.0e9, "MaxFreq": 32.0e9, "FreqStep": 0.1e9, "SaveStep": 1, "AdaptiveTol": 1.0e-3 },
    "Linear": { "Tol": 1.0e-8, "MaxIter": 200 }
  }
}"#,
        source_code: "// CPW WavePort mode 1, adaptive sweep 2�?2 GHz",
    },

    ExampleMeta {
        key: "cpw_wave_uniform",
        label: "CPW / wave uniform (Driven)",
        problem_type: "Driven",
        status: ExampleStatus::Ready,
        config_json: r#"{
  "Problem": { "Type": "Driven", "Verbose": 2, "Output": "postpro/wave_uniform" },
  "Model": { "Mesh": "cpw_wave_0.msh", "L0": 1.0e-6, "Refinement": {} },
  "Domains": {
    "Materials": [
      { "Attributes": [1], "Permeability": 1.0, "Permittivity": 1.0, "LossTan": 0.0 },
      {
        "Attributes": [2],
        "Permeability": [0.99999975, 0.99999975, 0.99999979],
        "Permittivity": [9.3, 9.3, 11.5],
        "LossTan": [3.0e-5, 3.0e-5, 8.6e-5],
        "MaterialAxes": [[0.8, 0.6, 0.0], [-0.6, 0.8, 0.0], [0.0, 0.0, 1.0]]
      }
    ]
  },
  "Boundaries": {
    "PEC": { "Attributes": [8, 9, 11] },
    "Absorbing": { "Attributes": [10], "Order": 1 },
    "WavePort": [
      { "Index": 1, "Attributes": [4], "Mode": 1, "Excitation": 1 },
      { "Index": 2, "Attributes": [5], "Mode": 1, "Excitation": 2 },
      { "Index": 3, "Attributes": [6], "Mode": 1 },
      { "Index": 4, "Attributes": [7], "Mode": 1 }
    ]
  },
  "Solver": {
    "Order": 2,
    "Driven": { "Samples": [{ "Type": "Linear", "MinFreq": 2.0, "MaxFreq": 32.0, "FreqStep": 6.0 }] },
    "Linear": { "Type": "Default", "KSPType": "GMRES", "Tol": 1.0e-8, "MaxIts": 200 }
  }
}"#,
        source_code: "// Anisotropic sapphire substrate, SurfaceFlux + Dielectric interface loss postprocessing",
    },

    ExampleMeta {
        key: "cpw_wave_eigen",
        label: "CPW / wave (Eigenmode)",
        problem_type: "Eigenmode",
        status: ExampleStatus::Ready,
        config_json: r#"{
  "Problem": { "Type": "Eigenmode", "Verbose": 2, "Output": "output/cpw_wave_eigen" },
  "Model": { "Mesh": "cpw_wave_0.msh", "L0": 1.0e-6, "Refinement": {} },
  "Domains": {
    "Materials": [
      { "Attributes": [1], "Permeability": 1.0, "Permittivity": 1.0, "LossTan": 0.0 },
      { "Attributes": [2], "Permeability": 0.99999975, "Permittivity": 9.3, "LossTan": 3.0e-5 }
    ]
  },
  "Boundaries": {
    "PEC": { "Attributes": [8, 9, 11] },
    "Absorbing": { "Attributes": [10], "Order": 1 },
    "WavePort": [
      { "Index": 1, "Attributes": [4], "Mode": 1 },
      { "Index": 2, "Attributes": [5], "Mode": 1 },
      { "Index": 3, "Attributes": [6], "Mode": 1 },
      { "Index": 4, "Attributes": [7], "Mode": 1 }
    ]
  },
  "Solver": {
    "Order": 2,
    "Eigenmode": { "N": 1, "Tol": 1.0e-6, "Target": 10.0e9, "Save": 1 },
    "Linear": { "Tol": 1.0e-8, "MaxIter": 200 }
  }
}"#,
        source_code: "// CPW WavePort, N=1 eigenmode near 10 GHz",
    },

    // ── Palace 对齐示例：cylinder ────────────────────────────────────────────

    ExampleMeta {
        key: "cavity_pec",
        label: "Cylinder / cavity PEC (Eigenmode)",
        problem_type: "Eigenmode",
        status: ExampleStatus::Ready,
        config_json: r#"{
  "Problem": { "Type": "Eigenmode", "Verbose": 2, "Output": "output/cavity_pec" },
  "Model": { "Mesh": "cylinder_hex.msh", "L0": 1.0e-2 },
  "Domains": { "Materials": [{ "Attributes": [1], "Permeability": 1.0, "Permittivity": 2.08, "LossTan": 4.0e-4 }] },
  "Boundaries": { "PEC": { "Attributes": [2, 3, 4] } },
  "Solver": {
    "Order": 4,
    "Eigenmode": { "N": 15, "Tol": 1.0e-8, "Target": 2.0e9, "Save": 15 },
    "Linear": { "Tol": 1.0e-8, "MaxIter": 100 }
  }
}"#,
        source_code: "// PEC cylindrical cavity, Order 4, N=15 modes near 2 GHz",
    },

    ExampleMeta {
        key: "cavity_impedance",
        label: "Cylinder / cavity impedance (Eigenmode)",
        problem_type: "Eigenmode",
        status: ExampleStatus::Ready,
        config_json: r#"{
  "Problem": { "Type": "Eigenmode", "Verbose": 2, "Output": "output/cavity_impedance" },
  "Model": { "Mesh": "cylinder_prism.msh", "L0": 1.0e-2 },
  "Domains": { "Materials": [{ "Attributes": [1], "Permeability": 1.0, "Permittivity": 2.08, "LossTan": 4.0e-4 }] },
  "Boundaries": {
    "Impedance": [{ "Attributes": [2, 3, 4], "Rs": 0.0184 }]
  },
  "Solver": {
    "Order": 4,
    "Eigenmode": { "N": 15, "Tol": 1.0e-8, "Target": 2.0e9, "Save": 15 },
    "Linear": { "Tol": 1.0e-8, "MaxIter": 100 }
  }
}"#,
        source_code: "// Impedance BC (Rs=0.0184 Ω) cavity, verifies Q factor vs PEC",
    },

    ExampleMeta {
        key: "driven_wave",
        label: "Cylinder / driven wave (Driven)",
        problem_type: "Driven",
        status: ExampleStatus::Ready,
        config_json: r#"{
  "Problem": { "Type": "Driven", "Verbose": 2, "Output": "output/driven_wave" },
  "Model": { "Mesh": "cylinder_hex.msh", "L0": 1.0e-2 },
  "Domains": { "Materials": [{ "Attributes": [1], "Permeability": 1.0, "Permittivity": 2.08, "LossTan": 4.0e-4 }] },
  "Boundaries": {
    "WavePort": [{ "Index": 1, "Attributes": [2], "Mode": 1, "Excitation": true }],
    "PEC": { "Attributes": [3, 4] }
  },
  "Solver": {
    "Order": 4,
    "Driven": { "MinFreq": 2.5e9, "MaxFreq": 5.0e9, "FreqStep": 0.5e9, "SaveStep": 1 },
    "Linear": { "Tol": 1.0e-9, "MaxIter": 100 }
  }
}"#,
        source_code: "// Cylindrical waveguide driven by WavePort mode 1, 2.5�? GHz",
    },

    ExampleMeta {
        key: "waveguide",
        label: "Cylinder / waveguide (Eigenmode)",
        problem_type: "Eigenmode",
        status: ExampleStatus::Ready,
        config_json: r#"{
  "Problem": { "Type": "Eigenmode", "Verbose": 2, "Output": "output/waveguide" },
  "Model": { "Mesh": "cylinder_tet.msh", "L0": 1.0e-2 },
  "Domains": { "Materials": [{ "Attributes": [1], "Permeability": 1.0, "Permittivity": 2.08, "LossTan": 4.0e-4 }] },
  "Boundaries": { "PEC": { "Attributes": [4] } },
  "Solver": {
    "Order": 4,
    "Eigenmode": { "N": 15, "Tol": 1.0e-9, "Target": 2.0e9, "Save": 15 },
    "Linear": { "Tol": 1.0e-9, "MaxIter": 100 }
  }
}"#,
        source_code: "// Cylindrical waveguide, Order 4, N=15 modes, tet mesh",
    },

    ExampleMeta {
        key: "floquet",
        label: "Cylinder / Floquet (Eigenmode)",
        problem_type: "Eigenmode",
        status: ExampleStatus::Ready,
        config_json: r#"{
  "Problem": { "Type": "Eigenmode", "Verbose": 2, "Output": "output/floquet" },
  "Model": { "Mesh": "cylinder_tet.msh", "L0": 1.0e-2 },
  "Domains": { "Materials": [{ "Attributes": [1], "Permeability": 1.0, "Permittivity": 2.08, "LossTan": 4.0e-4 }] },
  "Boundaries": { "PEC": { "Attributes": [4] } },
  "Solver": {
    "Order": 4,
    "Eigenmode": { "N": 15, "Tol": 1.0e-8, "Target": 2.0e9, "Save": 15 },
    "Linear": { "Tol": 1.0e-8, "MaxIter": 200 }
  }
}"#,
        source_code: "// Floquet Eigenmode on tet mesh (same geometry as waveguide, k·L=0 / Γ-point)",
    },

    // ── Palace 对齐示例：cylinder (REM magnetostatic) ─────────────────────────

    ExampleMeta {
        key: "cylinder",
        label: "Cylinder (Magnetostatic)",
        problem_type: "Magnetostatic",
        status: ExampleStatus::Ready,
        config_json: r#"{
  "Problem": { "Type": "Magnetostatic", "Verbose": 1, "Output": "./output/cylinder" },
  "Model": { "Mesh": "cylinder.msh", "L0": 1.0e-3 },
  "Domains": { "Materials": [{ "Attributes": [1], "Permeability": 1000.0 }] },
  "Boundaries": {
    "Ground": { "Attributes": [3] },
    "SurfaceCurrent": [{ "Index": 1, "Attributes": [2] }]
  },
  "Solver": { "Order": 1 }
}"#,
        source_code: "// Magnetostatic cylinder example\n// Physical groups: 1=cylinder(vol), 2=top, 3=bottom, 4=exterior, 5=symmetry",
    },

    // ── Palace 对齐示例：cpw (REM baseline) ──────────────────────────────────

    ExampleMeta {
        key: "cpw",
        label: "CPW baseline (Driven)",
        problem_type: "Driven",
        status: ExampleStatus::Ready,
        config_json: r#"{
  "Problem": { "Type": "Driven", "Verbose": 1, "Output": "./output/cpw" },
  "Model": { "Mesh": "cpw_coax.msh", "L0": 1.0e-6 },
  "Domains": { "Materials": [{ "Attributes": [1], "Permittivity": 11.7 }] },
  "Boundaries": {
    "PEC": { "Attributes": [2] },
    "LumpedPort": [{ "Index": 1, "Attributes": [3], "R": 50.0 }]
  },
  "Solver": { "Order": 2, "Driven": { "MinFreq": 4.0e9, "MaxFreq": 8.0e9, "FreqStep": 0.1e9 } }
}"#,
        source_code: "// REM CPW baseline: single LumpedPort, Silicon ε=11.7, 4�? GHz",
    },

    // ── Palace 对齐示例：transmon ─────────────────────────────────────────────

    ExampleMeta {
        key: "transmon_coarse",
        label: "Transmon / coarse (Eigenmode)",
        problem_type: "Eigenmode",
        status: ExampleStatus::Ready,
        config_json: r#"{
  "Problem": { "Type": "Eigenmode", "Verbose": 2, "Output": "output/transmon_coarse" },
  "Model": { "Mesh": "transmon.msh2", "L0": 1.0e-6, "Refinement": { "MaxIter": 0 } },
  "Domains": {
    "Materials": [
      { "Attributes": [2], "Permittivity": 1.0, "Permeability": 1.0 },
      { "Attributes": [1], "Permittivity": 9.3, "Permeability": 0.99999975, "LossTan": 3.0e-5 }
    ]
  },
  "Boundaries": {
    "PEC": { "Attributes": [5] },
    "Absorbing": { "Attributes": [3], "Order": 1 },
    "LumpedPort": [
      { "Index": 1, "Attributes": [6], "R": 50.0, "Direction": "+X" },
      { "Index": 2, "Attributes": [7], "R": 50.0, "Direction": "+X" },
      { "Index": 3, "Attributes": [4], "C": 5.5e-15, "L": 1.486e-8, "Direction": "+Y" }
    ]
  },
  "Solver": {
    "Order": 2,
    "Eigenmode": { "N": 2, "Save": 2, "Tol": 1.0e-8, "Target": 4.0e9 },
    "Linear": { "Tol": 1.0e-12, "MaxIter": 500 }
  }
}"#,
        source_code: "// Transmon coarse: JJ as C=5.5fF / L=14.86nH lumped port, no AMR",
    },

    ExampleMeta {
        key: "transmon_amr",
        label: "Transmon / AMR (Eigenmode)",
        problem_type: "Eigenmode",
        status: ExampleStatus::Ready,
        config_json: r#"{
  "Problem": { "Type": "Eigenmode", "Verbose": 2, "Output": "output/transmon_amr" },
  "Model": { "Mesh": "transmon.msh2", "L0": 1.0e-6, "Refinement": { "MaxIter": 2 } },
  "Domains": {
    "Materials": [
      { "Attributes": [2], "Permittivity": 1.0, "Permeability": 1.0 },
      { "Attributes": [1], "Permittivity": 9.3, "Permeability": 0.99999975, "LossTan": 3.0e-5 }
    ]
  },
  "Boundaries": {
    "PEC": { "Attributes": [5] },
    "Absorbing": { "Attributes": [3], "Order": 1 },
    "LumpedPort": [
      { "Index": 1, "Attributes": [6], "R": 50.0, "Direction": "+X" },
      { "Index": 2, "Attributes": [7], "R": 50.0, "Direction": "+X" },
      { "Index": 3, "Attributes": [4], "C": 5.5e-15, "L": 1.486e-8, "Direction": "+Y" }
    ]
  },
  "Solver": {
    "Order": 3,
    "Eigenmode": { "N": 2, "Save": 2, "Tol": 1.0e-8, "Target": 4.0e9 },
    "Linear": { "Tol": 1.0e-12, "MaxIter": 500 }
  }
}"#,
        source_code: "// Transmon AMR: 2 rounds adaptive mesh refinement, Order 3",
    },

    ExampleMeta {
        key: "transmon",
        label: "Transmon REM (Eigenmode)",
        problem_type: "Eigenmode",
        status: ExampleStatus::Ready,
        config_json: r#"{
  "Problem": { "Type": "Eigenmode", "Verbose": 1, "Output": "./output/transmon" },
  "Model": { "Mesh": "transmon.msh", "L0": 1.0e-6 },
  "Domains": { "Materials": [{ "Attributes": [1], "Permittivity": 11.4 }] },
  "Boundaries": { "PEC": { "Attributes": [2] } },
  "Solver": { "Order": 2, "Eigenmode": { "N": 5, "Tol": 1.0e-8, "Target": 5.0e9 } }
}"#,
        source_code: "// REM simplified transmon: PEC cavity, Si ε=11.4, N=5 modes",
    },
];

pub fn find_example(key: &str) -> Option<&'static ExampleMeta> {
    EXAMPLES.iter().find(|e| e.key == key)
}

pub fn get_config_json(key: &str) -> &'static str {
  match key {
    // REM + Palace-aligned example config files
    "spheres" => include_str!("../../../examples/palace/spheres/spheres.json"),
    "rings" => include_str!("../../../examples/palace/rings/rings.json"),
    "parallel_plate" => include_str!("../../../examples/rem/parallel_plate/parallel_plate.json"),
    "rem_es_fast" => include_str!("../../../examples/rem/es_parallel_plate_fast/es_parallel_plate_fast.json"),
    "rem_ms_fast" => include_str!("../../../examples/rem/ms_parallel_plate_fast/ms_parallel_plate_fast.json"),
    "rem_driven_fast" => include_str!("../../../examples/rem/driven_cpw_fast/driven_cpw_fast.json"),
    "rem_eigen_fast" => include_str!("../../../examples/rem/eigen_cylinder_fast/eigen_cylinder_fast.json"),
    "rem_transient_fast" => include_str!("../../../examples/rem/transient_coax_fast/transient_coax_fast.json"),
    "rem_mom_fast" => include_str!("../../../examples/rem/mom_sphere_fast/mom_sphere_fast.json"),
    "rem_sbr_fast" => include_str!("../../../examples/rem/sbr_sphere_fast/sbr_sphere_fast.json"),

    "adapter" => include_str!("../../../examples/palace/adapter/hybrid.json"),

    "antenna_halfwave_dipole" => {
      include_str!("../../../examples/palace/antenna/antenna_halfwave_dipole.json")
    }
    "antenna_short_dipole" => include_str!("../../../examples/palace/antenna/antenna_short_dipole.json"),

    "coaxial" => include_str!("../../../examples/palace/coaxial/coaxial.json"),
    "coaxial_matched" => include_str!("../../../examples/palace/coaxial/coaxial_matched.json"),
    "coaxial_open" => include_str!("../../../examples/palace/coaxial/coaxial_open.json"),
    "coaxial_short" => include_str!("../../../examples/palace/coaxial/coaxial_short.json"),

    "cpw" => include_str!("../../../examples/palace/cpw/cpw.json"),
    "cpw_coax_adaptive" => include_str!("../../../examples/palace/cpw/cpw_coax_adaptive.json"),
    "cpw_coax_uniform" => include_str!("../../../examples/palace/cpw/cpw_coax_uniform.json"),
    "cpw_lumped_adaptive" => include_str!("../../../examples/palace/cpw/cpw_lumped_adaptive.json"),
    "cpw_lumped_uniform" => include_str!("../../../examples/palace/cpw/cpw_lumped_uniform.json"),
    "cpw_lumped_eigen" => include_str!("../../../examples/palace/cpw/cpw_lumped_eigen.json"),
    "cpw_wave_adaptive" => include_str!("../../../examples/palace/cpw/cpw_wave_adaptive.json"),
    "cpw_wave_uniform" => include_str!("../../../examples/palace/cpw/cpw_wave_uniform.json"),
    "cpw_wave_eigen" => include_str!("../../../examples/palace/cpw/cpw_wave_eigen.json"),

    "cavity_pec" => include_str!("../../../examples/palace/cylinder/cavity_pec.json"),
    "cavity_impedance" => include_str!("../../../examples/palace/cylinder/cavity_impedance.json"),
    "driven_wave" => include_str!("../../../examples/palace/cylinder/driven_wave.json"),
    "waveguide" => include_str!("../../../examples/palace/cylinder/waveguide.json"),
    "floquet" => include_str!("../../../examples/palace/cylinder/floquet.json"),
    "cylinder" => include_str!("../../../examples/palace/cylinder/cylinder.json"),

    "sbr_sphere" => find_example(key).map(|e| e.config_json).unwrap_or("{}"),

    "transmon" => include_str!("../../../examples/palace/transmon/transmon.json"),
    "transmon_coarse" => include_str!("../../../examples/palace/transmon/transmon_coarse.json"),
    "transmon_amr" => include_str!("../../../examples/palace/transmon/transmon_amr.json"),

    // No dedicated file currently; keep embedded demo config.
    "mom_sphere" => find_example(key).map(|e| e.config_json).unwrap_or("{}"),

    _ => find_example(key).map(|e| e.config_json).unwrap_or("{}"),
  }
}

pub fn get_mesh_bytes(key: &str) -> Vec<u8> {
    match key {
    // Palace-aligned meshes
    "spheres" => include_bytes!("../../../examples/palace/spheres/mesh/spheres.msh").to_vec(),
    "rings"   => include_bytes!("../../../examples/palace/rings/mesh/rings.msh").to_vec(),

        // Palace: adapter
        "adapter" => include_bytes!("../../../examples/palace/adapter/mesh/adapter.msh").to_vec(),

        // Palace: antenna
        "antenna_halfwave_dipole" | "antenna_short_dipole"
            => include_bytes!("../../../examples/palace/antenna/mesh/antenna.msh").to_vec(),

        // Palace: coaxial
        "coaxial" | "coaxial_matched" | "coaxial_open" | "coaxial_short" | "rem_transient_fast"
            => include_bytes!("../../../examples/palace/coaxial/mesh/coaxial.msh").to_vec(),

        // Palace: cpw �?coax port variants
        "cpw_coax_adaptive" | "cpw_coax_uniform"
            => include_bytes!("../../../examples/palace/cpw/mesh/cpw_coax_0.msh").to_vec(),

        // Palace: cpw �?lumped port variants
        "cpw_lumped_adaptive" | "cpw_lumped_uniform" | "cpw_lumped_eigen"
            => include_bytes!("../../../examples/palace/cpw/mesh/cpw_lumped_0.msh").to_vec(),

        // Palace: cpw �?wave port variants
        "cpw_wave_adaptive" | "cpw_wave_uniform" | "cpw_wave_eigen"
            => include_bytes!("../../../examples/palace/cpw/mesh/cpw_wave_0.msh").to_vec(),

        // REM CPW baseline
        "cpw" | "rem_driven_fast"
          => include_bytes!("../../../examples/palace/cpw/mesh/cpw_coax.msh").to_vec(),

        // Palace: cylinder �?hex (cavity_pec, driven_wave, magnetostatic)
        "cavity_pec" | "driven_wave" | "cylinder"
            => include_bytes!("../../../examples/palace/cylinder/mesh/cylinder_hex.msh").to_vec(),

        // Palace: cylinder �?prism (cavity_impedance)
        "cavity_impedance"
            => include_bytes!("../../../examples/palace/cylinder/mesh/cylinder_prism.msh").to_vec(),

        // Palace: cylinder �?tet (waveguide, floquet)
        "waveguide" | "floquet" | "rem_eigen_fast"
            => include_bytes!("../../../examples/palace/cylinder/mesh/cylinder_tet.msh").to_vec(),

        // REM: parallel_plate
        "parallel_plate" | "rem_es_fast" | "rem_ms_fast"
            => include_bytes!("../../../examples/rem/parallel_plate/mesh/plate_2d.msh").to_vec(),

        // SBR / MoM
        "sbr_sphere" | "mom_sphere" | "rem_sbr_fast" | "rem_mom_fast"
            => include_bytes!("../../../examples/rem/sbr_sphere/mesh/sphere.msh").to_vec(),

        // Palace: transmon
        "transmon" | "transmon_coarse" | "transmon_amr"
            => include_bytes!("../../../examples/palace/transmon/mesh/transmon.msh2").to_vec(),

        _ => vec![],
    }
}
