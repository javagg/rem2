fn main() {
    // Link against system MPI when native-mpi feature is active.
    let native_mpi = std::env::var("CARGO_FEATURE_NATIVE_MPI").is_ok();
    if native_mpi {
        println!("cargo:rustc-link-lib=dylib=mpi");
        println!("cargo:rerun-if-env-changed=MPI_HOME");
    }
    println!("cargo:rerun-if-changed=build.rs");
}
