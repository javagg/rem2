pub mod element;
pub mod mesh;
pub mod node;
pub mod topology;

pub use element::{Element, ElementType};
pub use mesh::Mesh;
pub use node::Node;
pub use topology::{GEdge, GFace, GModel, GRegion, GSelection, GVertex};

// Backward compatibility aliases
pub use topology::{TopoEdge, TopoFace, TopoSelection, TopoVertex, TopoVolume, Topology};

/// 3D point type alias
pub type Point3 = nalgebra::Point3<f64>;
/// 3D vector type alias
pub type Vector3 = nalgebra::Vector3<f64>;
