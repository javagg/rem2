use crate::gmsh::RawMesh;
use rem_config::{Boundaries, PalaceConfig};
use rem_core::{RemError, RemResult};
use std::collections::HashMap;

// ---------------------------------------------------------------------------
// Geometry primitives
// ---------------------------------------------------------------------------

/// A mesh node with coordinates scaled by `L0`.
#[derive(Debug, Clone)]
pub struct Node {
    pub id: usize,
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

/// Element type (subset of GMSH types we actually use).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ElementKind {
    Line2,
    Tri3,
    Tri6,
    Quad4,
    Tet4,
    Tet10,
    Hex8,
}

impl ElementKind {
    pub fn n_nodes(self) -> usize {
        match self {
            ElementKind::Line2 => 2,
            ElementKind::Tri3  => 3,
            ElementKind::Tri6  => 6,
            ElementKind::Quad4 => 4,
            ElementKind::Tet4  => 4,
            ElementKind::Tet10 => 10,
            ElementKind::Hex8  => 8,
        }
    }

    pub fn dim(self) -> u8 {
        match self {
            ElementKind::Line2 => 1,
            ElementKind::Tri3 | ElementKind::Tri6 | ElementKind::Quad4 => 2,
            ElementKind::Tet4 | ElementKind::Tet10 | ElementKind::Hex8 => 3,
        }
    }

    /// Try to parse from GMSH element type integer.
    pub fn from_gmsh_type(t: u32) -> Option<Self> {
        match t {
            1  => Some(ElementKind::Line2),
            2  => Some(ElementKind::Tri3),
            3  => Some(ElementKind::Quad4),
            4  => Some(ElementKind::Tet4),
            5  => Some(ElementKind::Hex8),
            9  => Some(ElementKind::Tri6),
            11 => Some(ElementKind::Tet10),
            _  => None,
        }
    }
}

/// A mesh element (volume or boundary face/edge).
#[derive(Debug, Clone)]
pub struct Element {
    pub id: usize,
    pub kind: ElementKind,
    /// Physical group tag (from GMSH)
    pub tag: u32,
    /// Node indices (0-based into `RemMesh::nodes`)
    pub node_ids: Vec<usize>,
    /// MPI rank that owns this element (used during partitioning)
    pub rank: i32,
}

// ---------------------------------------------------------------------------
// Boundary condition tag (resolved from config)
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq)]
pub enum BoundaryTag {
    Pec,
    Pmc,
    Ground,
    ZeroCharge,
    Impedance { rs: f64, ls: f64, cs: f64 },
    /// Electrostatic terminal (Palace "Terminal"): equipotential conductor surface.
    Terminal { index: u32 },
    LumpedPort { index: u32, r: f64 },
    WavePort   { index: u32 },
    Absorbing  { order: u8 },
    SurfaceCurrent { index: u32 },
}

// ---------------------------------------------------------------------------
// RemMesh
// ---------------------------------------------------------------------------

/// The central mesh structure used by all solvers.
pub struct RemMesh {
    /// All nodes, coordinates in SI metres (after L0 scaling).
    pub nodes: Vec<Node>,

    /// All volume elements (highest-dimension elements).
    pub volume_elements: Vec<Element>,

    /// All boundary elements (dim = mesh_dim − 1).
    pub boundary_elements: Vec<Element>,

    /// Physical tag → material index in `PalaceConfig::domains.materials`.
    pub domain_tags: HashMap<u32, usize>,

    /// Physical tag → boundary condition.
    pub boundary_tags: HashMap<u32, BoundaryTag>,

    /// Spatial dimension of the mesh (2 or 3).
    pub dim: u8,
    /// MPI rank (for distributed solver)
    pub rank: i32,
    /// Total MPI processes
    pub size: i32,
}

impl RemMesh {
    /// Build from raw GMSH data + Palace config.
    pub fn from_raw(raw: RawMesh, config: &PalaceConfig) -> RemResult<Self> {
        let l0 = config.model.l0;

        // Scale coordinates
        let nodes: Vec<Node> = raw
            .nodes
            .into_iter()
            .map(|(id, x, y, z)| Node { id, x: x * l0, y: y * l0, z: z * l0 })
            .collect();

        // Classify elements by dimension
        let mesh_dim = raw
            .elements
            .iter()
            .filter_map(|e| ElementKind::from_gmsh_type(e.elem_type))
            .map(|k| k.dim())
            .max()
            .unwrap_or(2);

        let mut volume_elements = Vec::new();
        let mut boundary_elements = Vec::new();

        for re in &raw.elements {
            let kind = match ElementKind::from_gmsh_type(re.elem_type) {
                Some(k) => k,
                None => {
                    log::warn!("Skipping unsupported GMSH element type {}", re.elem_type);
                    continue;
                }
            };
            let elem = Element {
                id:       re.id,
                kind,
                tag:      re.phys_tag,
                node_ids: re.node_ids.iter().map(|&n| n - 1).collect(), // GMSH 1-based → 0-based
                rank:     0, // default to rank 0
            };
            if kind.dim() == mesh_dim {
                volume_elements.push(elem);
            } else if kind.dim() == mesh_dim - 1 {
                boundary_elements.push(elem);
            }
        }

        // Build domain tag map (physical group → material index)
        let mut domain_tags: HashMap<u32, usize> = HashMap::new();
        for (mat_idx, mat) in config.domains.materials.iter().enumerate() {
            for &attr in &mat.attributes {
                domain_tags.insert(attr, mat_idx);
            }
        }

        // Build boundary tag map
        let boundary_tags = build_boundary_tags(&config.boundaries)?;

        Ok(RemMesh {
            nodes,
            volume_elements,
            boundary_elements,
            domain_tags,
            boundary_tags,
            dim: mesh_dim,
            rank: 0,
            size: 1,
        })
    }

    pub fn set_comm(&mut self, rank: i32, size: i32) {
        self.rank = rank;
        self.size = size;
    }

    pub fn n_nodes(&self) -> usize { self.nodes.len() }
    pub fn n_volume_elements(&self) -> usize { self.volume_elements.len() }
    pub fn n_boundary_elements(&self) -> usize { self.boundary_elements.len() }

    /// Partition the volume elements among ranks using a simple geometric split (along X axis).
    /// This is a temporary placeholder for METIS partitioning.
    pub fn partition_geometric(&mut self) {
        if self.size <= 1 { return; }

        // Find bounding box for volume element centroids
        let mut x_min = f64::INFINITY;
        let mut x_max = f64::NEG_INFINITY;

        let centroids: Vec<f64> = self.volume_elements.iter().map(|e| {
            let mut cx = 0.0;
            for &nid in &e.node_ids {
                cx += self.nodes[nid].x;
            }
            cx /= e.node_ids.len() as f64;
            if cx < x_min { x_min = cx; }
            if cx > x_max { x_max = cx; }
            cx
        }).collect();

        let dx = x_max - x_min;
        if dx < 1e-12 {
            // Fallback: split by index if no spatial variation
            let n = self.volume_elements.len();
            for (i, e) in self.volume_elements.iter_mut().enumerate() {
                e.rank = (i * self.size as usize / n) as i32;
            }
        } else {
            for (e, cx) in self.volume_elements.iter_mut().zip(centroids) {
                let mut r = ((cx - x_min) / dx * self.size as f64) as i32;
                if r >= self.size { r = self.size - 1; }
                e.rank = r;
            }
        }

        // Boundary elements inherit rank from adjacent volume element (simplified)
        // In a real partitioner, we'd find which volume element it belongs to.
        // For now, use the rank of the first node's first adjacent volume element if we had that map.
        // Simplified: just use first node's position.
        for e in &mut self.boundary_elements {
             let mut cx = 0.0;
             for &nid in &e.node_ids {
                 cx += self.nodes[nid].x;
             }
             cx /= e.node_ids.len() as f64;
             let mut r = ((cx - x_min) / dx * self.size as f64) as i32;
             if r < 0 { r = 0; }
             if r >= self.size { r = self.size - 1; }
             e.rank = r;
        }
    }

    /// Return all boundary element tags that map to the given BoundaryTag variant.
    pub fn boundary_element_tags_of<F>(&self, predicate: F) -> Vec<u32>
    where
        F: Fn(&BoundaryTag) -> bool,
    {
        self.boundary_tags
            .iter()
            .filter_map(|(tag, bc)| if predicate(bc) { Some(*tag) } else { None })
            .collect()
    }
}

// ---------------------------------------------------------------------------
// Helper: build boundary tag map from Palace config
// ---------------------------------------------------------------------------

fn build_boundary_tags(b: &Boundaries) -> RemResult<HashMap<u32, BoundaryTag>> {
    let mut map: HashMap<u32, BoundaryTag> = HashMap::new();

    macro_rules! insert_unique {
        ($tag:expr, $bc:expr) => {{
            if let Some(existing) = map.insert($tag, $bc) {
                return Err(RemError::Config(format!(
                    "physical group {} is assigned to multiple boundary conditions (e.g. {:?})",
                    $tag, existing
                )));
            }
        }};
    }

    if let Some(pec) = &b.pec {
        for &t in &pec.attributes { insert_unique!(t, BoundaryTag::Pec); }
    }
    if let Some(pmc) = &b.pmc {
        for &t in &pmc.attributes { insert_unique!(t, BoundaryTag::Pmc); }
    }
    if let Some(gnd) = &b.ground {
        for &t in &gnd.attributes { insert_unique!(t, BoundaryTag::Ground); }
    }
    if let Some(zc) = &b.zero_charge {
        for &t in &zc.attributes { insert_unique!(t, BoundaryTag::ZeroCharge); }
    }
    for term in &b.terminal {
        let bc = BoundaryTag::Terminal { index: term.index };
        for &t in &term.attributes { insert_unique!(t, bc.clone()); }
    }
    for imp in &b.impedance {
        let bc = BoundaryTag::Impedance { rs: imp.rs, ls: imp.ls, cs: imp.cs };
        for &t in &imp.attributes { insert_unique!(t, bc.clone()); }
    }
    for port in &b.lumped_port {
        let bc = BoundaryTag::LumpedPort { index: port.index, r: port.r };
        for &t in &port.attributes { insert_unique!(t, bc.clone()); }
    }
    for port in &b.wave_port {
        let bc = BoundaryTag::WavePort { index: port.index };
        for &t in &port.attributes { insert_unique!(t, bc.clone()); }
    }
    if let Some(abs) = &b.absorbing {
        let bc = BoundaryTag::Absorbing { order: abs.order };
        for &t in &abs.attributes { insert_unique!(t, bc.clone()); }
    }
    for sc in &b.surface_current {
        let bc = BoundaryTag::SurfaceCurrent { index: sc.index };
        for &t in &sc.attributes { insert_unique!(t, bc.clone()); }
    }

    Ok(map)
}
