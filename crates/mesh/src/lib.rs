pub mod gen;
pub mod gmsh;
pub mod mesh_data;
pub mod amr;
pub mod fem_bridge;
pub mod p_refine;

pub use mesh_data::{BoundaryTag, Node, Element, ElementKind, RemMesh};
pub use fem_bridge::{FemSubMesh2d, extract_submesh_tri3, extract_submesh_by_element_ids_tri3, refine_marked_tri3};
pub use p_refine::{p_refine_mesh, p3_refine_mesh};

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

/// Load a mesh from bytes (GMSH v2.2 or v4.1, ASCII or binary). Used in WASM where there is no filesystem.
pub fn load_mesh_from_bytes(config: &PalaceConfig, data: &[u8], comm: &impl Comm) -> RemResult<RemMesh> {
    let raw = gmsh::read_msh_bytes(data)?;
    let mut mesh = RemMesh::from_raw(raw, config)?;
    mesh.set_comm(comm.rank(), comm.size());
    Ok(mesh)
}
