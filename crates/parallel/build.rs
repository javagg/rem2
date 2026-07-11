fn main() {
    let have_mpi = std::env::var("CARGO_FEATURE_MPI").is_ok();
    if have_mpi {
        println!("cargo:rustc-link-lib=dylib=mpi");
        println!("cargo:rerun-if-env-changed=MPI_HOME");
    }
    println!("cargo:rerun-if-changed=build.rs");
}
