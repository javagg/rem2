use serde::{Deserialize, Serialize};

/// A geometric vertex — a corner point where 3+ edges meet.
/// Dimension: 0 (point)
/// Belongs to a GModel (B-Rep topology).
///
/// Relationship with mesh entities:
/// - Contains a single node ID (typically Point1 element nodes, but can be any node)
/// - Appears as an endpoint of GEdges
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GVertex {
    pub id: usize,
    /// Mesh node ID at this geometric point.
    pub node_id: u64,
}

impl GVertex {
    /// Topological dimension: 0 (point)
    pub fn dimension(&self) -> u8 {
        0
    }
}

/// A geometric edge — a curve along a sharp feature boundary.
/// Dimension: 1 (curve)
/// Belongs to a GModel (B-Rep topology).
///
/// Relationship with mesh entities:
/// - Contains ordered node IDs (typically from Line2 elements)
/// - Bounded by two GVertices (except for closed loops where vertex_ids = [None, None])
/// - Appears as a bounding edge of GFaces
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GEdge {
    pub id: usize,
    /// Start and end geometric vertex IDs. `None` for closed loops.
    pub vertex_ids: [Option<usize>; 2],
    /// Ordered mesh node IDs along this edge.
    pub node_ids: Vec<u64>,
}

impl GEdge {
    /// Topological dimension: 1 (curve)
    pub fn dimension(&self) -> u8 {
        1
    }
}

/// A geometric face — a smooth surface region bounded by sharp edges.
/// Dimension: 2 (surface)
/// Belongs to a GModel (B-Rep topology).
///
/// Relationship with mesh entities:
/// - Contains mesh face polygons (each is a node ID sequence)
/// - Polygons correspond to Triangle3, Quad4, or boundary faces of volume elements
/// - Bounded by GEdges (forming loops)
/// - Appears as a bounding face of GRegions
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GFace {
    pub id: usize,
    /// Bounding geometric edge IDs.
    pub edge_ids: Vec<usize>,
    /// Mesh face polygons making up this face (each polygon is a node ID sequence).
    /// Typically 3-node or 4-node polygons (Triangle3 or Quad4 element faces).
    pub mesh_faces: Vec<Vec<u64>>,
}

impl GFace {
    /// Topological dimension: 2 (surface)
    pub fn dimension(&self) -> u8 {
        2
    }
}

/// A geometric region — a solid region of connected volume elements.
/// Dimension: 3 (volume)
/// Belongs to a GModel (B-Rep topology).
///
/// Relationship with mesh entities:
/// - Contains volume element IDs (Tetrahedron4, Hexahedron8, Prism6, Pyramid5)
/// - Bounded by GFaces (which form the boundary surfaces)
/// - All elements in the region are topologically connected via face adjacency
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GRegion {
    pub id: usize,
    /// Bounding geometric face IDs.
    pub face_ids: Vec<usize>,
    /// Mesh element IDs belonging to this region.
    /// Typically 3D elements: Tetrahedron4, Hexahedron8, Prism6, Pyramid5
    pub element_ids: Vec<u64>,
}

impl GRegion {
    /// Topological dimension: 3 (volume)
    pub fn dimension(&self) -> u8 {
        3
    }
}

/// Complete B-Rep-style geometric model of a mesh, built via dihedral-angle classification
/// following the approach used in gmsh's `classifyFaces`.
///
/// A GModel contains geometric entities extracted from or inferred from the mesh:
/// - GVertices: geometric points
/// - GEdges: geometric curves
/// - GFaces: geometric surfaces
/// - GRegions: geometric solids
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GModel {
    pub vertices: Vec<GVertex>,
    pub edges: Vec<GEdge>,
    pub faces: Vec<GFace>,
    pub regions: Vec<GRegion>,
    /// Dihedral angle threshold used for classification (degrees).
    pub angle_threshold_deg: f64,
}

impl Default for GModel {
    fn default() -> Self {
        Self {
            vertices: Vec::new(),
            edges: Vec::new(),
            faces: Vec::new(),
            regions: Vec::new(),
            angle_threshold_deg: 40.0,
        }
    }
}

/// Selection of a geometric entity for highlighting.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GSelection {
    Region(usize),
    Face(usize),
    Edge(usize),
    Vertex(usize),
}

// Backward-compatibility aliases (will be deprecated)
pub type Topology = GModel;
pub type TopoVertex = GVertex;
pub type TopoEdge = GEdge;
pub type TopoFace = GFace;
pub type TopoVolume = GRegion;
pub type TopoSelection = GSelection;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn geometric_entities_have_correct_dimensions() {
        // 0-Dimensional: GVertex
        let gvert = GVertex {
            id: 1,
            node_id: 100,
        };
        assert_eq!(gvert.dimension(), 0);

        // 1-Dimensional: GEdge
        let gedge = GEdge {
            id: 2,
            vertex_ids: [Some(1), Some(2)],
            node_ids: vec![100, 101, 102],
        };
        assert_eq!(gedge.dimension(), 1);

        // 2-Dimensional: GFace
        let gface = GFace {
            id: 3,
            edge_ids: vec![1, 2, 3],
            mesh_faces: vec![vec![100, 101, 102], vec![102, 103, 104]],
        };
        assert_eq!(gface.dimension(), 2);

        // 3-Dimensional: GRegion
        let gregion = GRegion {
            id: 4,
            face_ids: vec![1, 2, 3, 4, 5, 6],
            element_ids: vec![1000, 1001, 1002],
        };
        assert_eq!(gregion.dimension(), 3);
    }

    #[test]
    fn geometric_model_contains_all_dimensions() {
        let mut gmodel = GModel::default();

        // Add vertices (0D)
        gmodel.vertices.push(GVertex { id: 1, node_id: 1 });
        gmodel.vertices.push(GVertex { id: 2, node_id: 2 });

        // Add edges (1D)
        gmodel.edges.push(GEdge {
            id: 1,
            vertex_ids: [Some(1), Some(2)],
            node_ids: vec![1, 2],
        });

        // Add faces (2D)
        gmodel.faces.push(GFace {
            id: 1,
            edge_ids: vec![1],
            mesh_faces: vec![vec![1, 2, 3]],
        });

        // Add regions (3D)
        gmodel.regions.push(GRegion {
            id: 1,
            face_ids: vec![1],
            element_ids: vec![100],
        });

        // Verify model structure
        assert_eq!(gmodel.vertices.len(), 2);
        assert_eq!(gmodel.edges.len(), 1);
        assert_eq!(gmodel.faces.len(), 1);
        assert_eq!(gmodel.regions.len(), 1);

        // Verify dimensions
        assert_eq!(gmodel.vertices[0].dimension(), 0);
        assert_eq!(gmodel.edges[0].dimension(), 1);
        assert_eq!(gmodel.faces[0].dimension(), 2);
        assert_eq!(gmodel.regions[0].dimension(), 3);
    }

    #[test]
    fn dimension_containment_relationships() {
        // Verify typical containment patterns:
        // GVertex (0D) contains 1 node
        let gvert = GVertex { id: 1, node_id: 5 };
        assert_eq!(gvert.dimension(), 0);
        // Single node represents a 0D point

        // GEdge (1D) contains ordered nodes from Line2 chain
        let gedge = GEdge {
            id: 1,
            vertex_ids: [Some(1), Some(3)],
            node_ids: vec![5, 10, 11, 12], // nodes from Line2 chain
        };
        assert_eq!(gedge.dimension(), 1);
        assert!(gedge.node_ids.len() >= 2);

        // GFace (2D) contains Triangle3/Quad4 node polygons
        let gface = GFace {
            id: 1,
            edge_ids: vec![1, 2, 3],
            mesh_faces: vec![
                vec![10, 11, 12],     // Triangle3
                vec![12, 13, 14, 15], // Quad4
            ],
        };
        assert_eq!(gface.dimension(), 2);
        assert!(gface.mesh_faces.iter().all(|f| f.len() >= 3)); // All faces have 3+ nodes

        // GRegion (3D) contains volume element IDs
        let gregion = GRegion {
            id: 1,
            face_ids: vec![1, 2, 3, 4, 5, 6],
            element_ids: vec![100, 101, 102], // Tet4, Hex8, etc.
        };
        assert_eq!(gregion.dimension(), 3);
        assert!(!gregion.element_ids.is_empty());
    }
}
