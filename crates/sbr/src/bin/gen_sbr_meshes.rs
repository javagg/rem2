//! Generate example SBR+ meshes.
//!
//! Usage:
//!   cargo run -p rem-sbr --bin gen_sbr_meshes

fn main() {
    use rem_mesh::gen::pec_sphere_msh;

    // PEC sphere: radius 0.5 m, 24 lat × 48 lon, tag=1
    let msh = pec_sphere_msh(0.5, 24, 48, 1);
    let out = "examples/sbr_sphere/mesh/sphere.msh";
    std::fs::write(out, &msh).expect("failed to write sphere.msh");
    println!("Wrote {} ({} bytes)", out, msh.len());
}
