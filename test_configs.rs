use rem_config::{load_config_from_str, ConfigFormat};

fn main() {
    println!("Testing example configs for parsing issues...\n");
    
    let examples = vec![
        ("spheres", r#"{"Problem":{"Type":"Electrostatic","Output":"."},"Model":{"Mesh":"spheres.msh","L0":0.001},"Domains":{"Materials":[{"Attributes":[10],"Permittivity":1.0,"Permeability":1.0}]},"Boundaries":{"Ground":{"Attributes":[2]},"Terminal":[{"Index":1,"Attributes":[1]}]},"Solver":{"Order":1,"Linear":{"Type":"GMRES","Tol":1e-10,"MaxIter":500}}}"#),
        ("rings", r#"{"Problem":{"Type":"Magnetostatic","Output":"."},"Model":{"Mesh":"rings.msh","L0":0.001},"Domains":{"Materials":[{"Attributes":[10],"Permittivity":1.0,"Permeability":1000.0},{"Attributes":[20],"Permittivity":1.0,"Permeability":1.0}]},"Boundaries":{"Ground":{"Attributes":[1]},"SurfaceCurrent":[{"Index":1,"Attributes":[2],"Direction":"+Y"}]},"Solver":{"Order":1,"Linear":{"Type":"GMRES","Tol":1e-10,"MaxIter":500}}}"#),
        ("coaxial", r#"{"Problem":{"Type":"Electrostatic","Output":"."},"Model":{"Mesh":"coaxial.msh","L0":0.001},"Domains":{"Materials":[{"Attributes":[1],"Permittivity":2.1}]},"Boundaries":{"Ground":{"Attributes":[2]},"Terminal":[{"Index":1,"Attributes":[3]}]},"Solver":{"Order":1}}"#),
        ("cylinder", r#"{"Problem":{"Type":"Magnetostatic","Output":"."},"Model":{"Mesh":"cylinder.msh","L0":0.001},"Domains":{"Materials":[{"Attributes":[1],"Permeability":1000.0}]},"Boundaries":{"Source":{"Current":[{"Index":1,"Attributes":[2],"Value":1.0}]}},"Solver":{"Order":1}}"#),
        ("adapter", r#"{"Problem":{"Type":"Driven","Output":"."},"Model":{"Mesh":"adapter.msh","L0":1.0},"Domains":{"Materials":[{"Attributes":[1],"Permittivity":1.0,"Permeability":1.0}]},"Boundaries":{"PEC":{"Attributes":[2]},"LumpedPort":[{"Index":1,"Attributes":[3],"R":50.0}]},"Solver":{"Order":1,"Driven":{"FreqStart":1e9,"FreqEnd":10e9,"FreqStep":0.1e9}}}"#),
        ("antenna", r#"{"Problem":{"Type":"Driven","Output":"."},"Model":{"Mesh":"antenna.msh","L0":0.001},"Domains":{"Materials":[{"Attributes":[1],"Permittivity":1.0,"Permeability":1.0}]},"Boundaries":{"Absorbing":{"Attributes":[2]},"LumpedPort":[{"Index":1,"Attributes":[3],"R":50.0}]},"Solver":{"Order":2,"Driven":{"FreqStart":2e9,"FreqEnd":3e9,"FreqStep":0.05e9}}}"#),
        ("cpw", r#"{"Problem":{"Type":"Driven","Output":"."},"Model":{"Mesh":"cpw.msh","L0":1e-6},"Domains":{"Materials":[{"Attributes":[1],"Permittivity":11.7}]},"Boundaries":{"PEC":{"Attributes":[2]},"LumpedPort":[{"Index":1,"Attributes":[3],"R":50.0}]},"Solver":{"Order":2,"Driven":{"MinFreq":4e9,"MaxFreq":8e9,"FreqStep":0.1e9}}}"#),
        ("transmon", r#"{"Problem":{"Type":"Eigenmode","Output":"."},"Model":{"Mesh":"transmon.msh","L0":1e-6},"Domains":{"Materials":[{"Attributes":[1],"Permittivity":11.4}]},"Boundaries":{"PEC":{"Attributes":[2]}},"Solver":{"Order":2,"Eigenmode":{"N":5,"Tol":1e-8,"Target":5e9}}}"#),
    ];

    for (name, config_json) in examples {
        match load_config_from_str(config_json, ConfigFormat::Json) {
            Ok(_cfg) => println!("✓ {}", name),
            Err(e) => println!("✗ {} - Error: {}", name, e),
        }
    }
}
