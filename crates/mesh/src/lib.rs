pub mod gen;
pub mod gmsh;
pub mod mesh_data;

pub use mesh_data::{BoundaryTag, Node, Element, ElementKind, RemMesh};

use rem_config::PalaceConfig;
use rem_core::RemResult;
use std::path::Path;

/// Load a mesh and bind it to the material/BC configuration.
use rem_parallel::Comm;

/// Load a mesh and bind it to the material/BC configuration.
pub fn load_mesh(config: &PalaceConfig, comm: &impl Comm) -> RemResult<RemMesh> {
    // Resolve mesh path relative to cwd (or absolute)
    let path = Path::new(&config.model.mesh);
    let raw = gmsh::read_msh_file(path)?;
    let mut mesh = RemMesh::from_raw(raw, config)?;
    mesh.set_comm(comm.rank(), comm.size());
    Ok(mesh)
}

/// Load a mesh from bytes (GMSH v4.1 ASCII). Used in WASM where there is no filesystem.
/// Load a mesh from bytes (GMSH v4.1 ASCII). Used in WASM where there is no filesystem.
pub fn load_mesh_from_bytes(config: &PalaceConfig, data: &[u8], comm: &impl Comm) -> RemResult<RemMesh> {
    let text = std::str::from_utf8(data)
        .map_err(|_| rem_core::RemError::Mesh("mesh bytes are not valid UTF-8".into()))?;
    let raw = gmsh::read_msh_str(text)?;
    let mut mesh = RemMesh::from_raw(raw, config)?;
    mesh.set_comm(comm.rank(), comm.size());
    Ok(mesh)
}
