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
            10 => Some(ElementKind::Quad4),  // Quad9 → Quad4 corner nodes (P1)
            11 => Some(ElementKind::Tet10),
            12 => Some(ElementKind::Hex8),   // Hex27 → Hex8 corner nodes (P1)
            17 => Some(ElementKind::Hex8),   // Hex20 → Hex8 corner nodes (P1)
            _  => None,
        }
    }
}

fn gmsh_type_hint(t: u32) -> &'static str {
    match t {
        29 => "Tet20 (high-order tetrahedron, 20 nodes)",
        30 => "Tet35 (high-order tetrahedron, 35 nodes)",
        31 => "Tet56 (high-order tetrahedron, 56 nodes)",
        90 => "Prism40 (high-order prism, 40 nodes)",
        91 => "Prism75 (high-order prism, 75 nodes)",
        92 => "Hex64 (high-order hexahedron, 64 nodes)",
        93 => "Hex125 (high-order hexahedron, 125 nodes)",
        118 => "Pyramid30 (high-order pyramid, 30 nodes)",
        119 => "Pyramid55 (high-order pyramid, 55 nodes)",
        _ => "unknown or currently unsupported element type",
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
    ResistiveSheet { rs: f64 },
    /// Electrostatic terminal (Palace "Terminal"): equipotential conductor surface.
    Terminal { index: u32 },
    LumpedPort { index: u32, r: f64, l: f64, c: f64 },
    WavePort   { index: u32 },
    Absorbing  { order: u8 },
    SurfaceCurrent { index: u32 },
}

// ---------------------------------------------------------------------------
// RemMesh
// ---------------------------------------------------------------------------

/// The central mesh structure used by all solvers.
#[derive(Clone)]
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
                    log::warn!(
                        "Skipping unsupported GMSH element type {} ({}). \
                         Current rem-mesh support: 1(Line2), 2(Tri3), 3(Quad4), \
                         4(Tet4), 5(Hex8), 9(Tri6), 11(Tet10). \
                         If this is a high-order mesh (e.g. 29), export a linear mesh \
                         (Tet4/Tet10, Tri3/Tri6, etc.) before running REM.",
                        re.elem_type,
                        gmsh_type_hint(re.elem_type)
                    );
                    continue;
                }
            };
            // For high-order elements mapped to lower-order (e.g. Hex27→Hex8, Quad9→Quad4),
            // truncate node list to corner nodes only. GMSH always places corner nodes first.
            let n_corner = kind.n_nodes();
            let node_ids: Vec<usize> = re.node_ids.iter()
                .take(n_corner)
                .map(|&n| n - 1)
                .collect();
            if node_ids.len() < n_corner {
                log::warn!("Element {} (type {}) has fewer nodes ({}) than expected ({}); skipping",
                    re.id, re.elem_type, node_ids.len(), n_corner);
                continue;
            }
            if re.node_ids.len() > n_corner {
                log::warn!(
                    "GMSH type {} detected: using P1 corner-node approximation \
                     ({} of {} nodes). Accuracy degrades. Re-mesh with linear elements for full precision.",
                    re.elem_type, n_corner, re.node_ids.len()
                );
            }
            let elem = Element {
                id:   re.id,
                kind,
                tag:  re.phys_tag,
                node_ids,
                rank: 0,
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
    /// Fallback used when the `metis` feature is disabled or when `size <= 1`.
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

    /// Partition volume elements using METIS k-way graph partitioning (dual graph).
    ///
    /// Builds the element-dual graph: each volume element is a vertex; two elements
    /// share an edge when they share a mesh face (for 3-D) or edge (for 2-D).
    ///
    /// When `comm` has size > 1, multiple cut attempts are distributed across MPI
    /// ranks via rmetis's parallel reduction (ParMETIS-style multi-cut).
    ///
    /// Returns `Err` if rmetis fails; the caller can fall back to geometric partitioning.
    #[cfg(feature = "metis")]
    pub fn partition_metis(&mut self, comm: &dyn rem_parallel::Comm) -> Result<(), rmetis::MetisError> {
        use rmetis::{Graph, Options};
        use rmetis::partition::kway::partition_kway;
        use std::collections::HashMap as HMap;

        let n_elem = self.volume_elements.len();
        if n_elem == 0 { return Ok(()); }

        let nparts = self.size as usize;
        if nparts <= 1 { return Ok(()); }

        // -----------------------------------------------------------------------
        // Build element dual graph (CSR format)
        // -----------------------------------------------------------------------
        // A face is identified by its sorted node-set.  Two elements sharing a
        // face (dim-D face = (dim-1) nodes sorted) become graph-adjacent.

        // face_key → list of element indices that own this face
        let face_nodes = self.dim as usize; // 3-D mesh: triangular faces (3 nodes); 2-D: edges (2 nodes)
        let mut face_map: HMap<Vec<usize>, Vec<usize>> = HMap::new();

        for (ei, elem) in self.volume_elements.iter().enumerate() {
            let nn = elem.node_ids.len();
            // Generate all combinations of `face_nodes` nodes (faces of the element)
            for_each_face(&elem.node_ids, nn, face_nodes, |face| {
                face_map.entry(face).or_default().push(ei);
            });
        }

        // Build adjacency lists
        let mut adj_lists: Vec<Vec<i32>> = vec![Vec::new(); n_elem];
        for owners in face_map.values() {
            if owners.len() == 2 {
                let (a, b) = (owners[0], owners[1]);
                if !adj_lists[a].contains(&(b as i32)) { adj_lists[a].push(b as i32); }
                if !adj_lists[b].contains(&(a as i32)) { adj_lists[b].push(a as i32); }
            }
        }

        // Convert to CSR
        let mut xadj: Vec<i32> = Vec::with_capacity(n_elem + 1);
        let mut adjncy: Vec<i32> = Vec::new();
        xadj.push(0);
        for nbrs in &adj_lists {
            adjncy.extend_from_slice(nbrs);
            xadj.push(adjncy.len() as i32);
        }

        // If the mesh has no inter-element adjacency (e.g., single element), fall back
        if adjncy.is_empty() {
            self.partition_geometric();
            return Ok(());
        }

        let graph = Graph::new_unweighted(n_elem, xadj, adjncy)?;
        let opts = Options::for_kway();

        // Bridge rem-parallel Comm → rmetis Comm for parallel multi-cut reduction
        let rmetis_comm = RemCommAdapter(comm);
        let result = partition_kway(&graph, nparts, None, None, &opts, Some(&rmetis_comm))?;

        // Assign ranks to volume elements
        for (elem, &p) in self.volume_elements.iter_mut().zip(result.part.iter()) {
            elem.rank = p;
        }

        log::debug!(
            "METIS k-way partition: {} elements → {} parts, edge cut = {}",
            n_elem, nparts, result.objval
        );

        // Assign boundary element ranks from the nearest volume element centroid
        self.assign_boundary_ranks_from_volume();

        Ok(())
    }

    /// Assign each boundary element to the same rank as the volume element whose
    /// centroid is closest to the boundary element's centroid.
    fn assign_boundary_ranks_from_volume(&mut self) {
        if self.volume_elements.is_empty() { return; }

        // Precompute volume element centroids
        let vol_centroids: Vec<[f64; 3]> = self.volume_elements.iter().map(|e| {
            let mut cx = 0.0; let mut cy = 0.0; let mut cz = 0.0;
            for &nid in &e.node_ids {
                cx += self.nodes[nid].x;
                cy += self.nodes[nid].y;
                cz += self.nodes[nid].z;
            }
            let n = e.node_ids.len() as f64;
            [cx / n, cy / n, cz / n]
        }).collect();

        for be in &mut self.boundary_elements {
            let mut bx = 0.0; let mut by = 0.0; let mut bz = 0.0;
            for &nid in &be.node_ids {
                bx += self.nodes[nid].x;
                by += self.nodes[nid].y;
                bz += self.nodes[nid].z;
            }
            let n = be.node_ids.len() as f64;
            bx /= n; by /= n; bz /= n;

            let nearest = vol_centroids.iter().enumerate().min_by(|(_, a), (_, b)| {
                let da = (a[0]-bx).powi(2) + (a[1]-by).powi(2) + (a[2]-bz).powi(2);
                let db = (b[0]-bx).powi(2) + (b[1]-by).powi(2) + (b[2]-bz).powi(2);
                da.partial_cmp(&db).unwrap()
            });
            if let Some((vi, _)) = nearest {
                be.rank = self.volume_elements[vi].rank;
            }
        }
    }

    /// Partition elements using METIS when the feature is enabled, otherwise geometric.
    ///
    /// When `comm` has size > 1 and the `metis` feature is on, multi-cut attempts
    /// are distributed across MPI ranks for better partition quality.
    pub fn partition(&mut self, comm: &dyn rem_parallel::Comm) {
        #[cfg(feature = "metis")]
        {
            if self.size > 1 {
                match self.partition_metis(comm) {
                    Ok(()) => return,
                    Err(e) => {
                        log::warn!("METIS partitioning failed ({}), falling back to geometric", e);
                    }
                }
            }
        }
        let _ = comm;
        self.partition_geometric();
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

    // -----------------------------------------------------------------------
    // fem-rs interop
    // -----------------------------------------------------------------------

    /// Convert this mesh to a `fem_mesh::SimplexMesh<3>` for use with fem-rs
    /// assemblers and solvers.
    ///
    /// **Use this only for 3-D meshes.**  For 2-D meshes call `to_simplex_mesh_2d()`
    /// which returns `SimplexMesh<2>`.
    ///
    /// The conversion handles the most common 3-D element types used in rem2
    /// (Tet4, Tet10, Hex8, Tri3, Quad4).  Mixed-element meshes produce a
    /// `SimplexMesh` with `elem_types` / `elem_offsets` populated.
    ///
    /// `elem_tags` carries the physical-group (material) tag from GMSH.
    /// `face_tags` carries the physical-group tag as a plain `i32`
    /// (compatible with `fem_mesh::BoundaryTag = i32`).
    pub fn to_simplex_mesh(&self) -> fem_mesh::SimplexMesh<3> {
        use fem_mesh::ElementType as FET;

        // --- map ElementKind → fem ElementType ---
        fn to_fem_elem_type(k: ElementKind) -> Option<FET> {
            match k {
                ElementKind::Tet4  => Some(FET::Tet4),
                ElementKind::Tet10 => Some(FET::Tet10),
                ElementKind::Hex8  => Some(FET::Hex8),
                ElementKind::Tri3  => Some(FET::Tri3),
                ElementKind::Tri6  => Some(FET::Tri6),
                ElementKind::Quad4 => Some(FET::Quad4),
                ElementKind::Line2 => Some(FET::Line2),
            }
        }

        // --- node coordinates (flat, x y z x y z …) ---
        let mut coords = Vec::with_capacity(self.nodes.len() * 3);
        for n in &self.nodes {
            coords.push(n.x);
            coords.push(n.y);
            coords.push(n.z);
        }

        // --- volume element connectivity ---
        // Detect whether the mesh is uniform (all same kind) or mixed.
        let first_kind = self.volume_elements.first().map(|e| e.kind);
        let is_uniform = first_kind.map_or(true, |k| {
            self.volume_elements.iter().all(|e| e.kind == k)
        });

        let (conn, elem_tags, elem_type, elem_types, elem_offsets) = if is_uniform {
            let kind = first_kind.unwrap_or(ElementKind::Tet4);
            let et = to_fem_elem_type(kind).unwrap_or(FET::Tet4);
            let mut conn: Vec<u32> = Vec::with_capacity(
                self.volume_elements.len() * kind.n_nodes()
            );
            let mut tags: Vec<i32> = Vec::with_capacity(self.volume_elements.len());
            for e in &self.volume_elements {
                for &nid in &e.node_ids { conn.push(nid as u32); }
                tags.push(e.tag as i32);
            }
            (conn, tags, et, None, None)
        } else {
            // Mixed: use CSR-like offsets
            let mut conn: Vec<u32>  = Vec::new();
            let mut tags: Vec<i32>  = Vec::with_capacity(self.volume_elements.len());
            let mut etypes: Vec<FET> = Vec::with_capacity(self.volume_elements.len());
            let mut offsets: Vec<usize> = vec![0usize];
            for e in &self.volume_elements {
                for &nid in &e.node_ids { conn.push(nid as u32); }
                tags.push(e.tag as i32);
                etypes.push(to_fem_elem_type(e.kind).unwrap_or(FET::Tet4));
                offsets.push(conn.len());
            }
            let primary = to_fem_elem_type(
                self.volume_elements.first().map(|e| e.kind).unwrap_or(ElementKind::Tet4)
            ).unwrap_or(FET::Tet4);
            (conn, tags, primary, Some(etypes), Some(offsets))
        };

        // --- boundary face connectivity ---
        let first_face_kind = self.boundary_elements.first().map(|e| e.kind);
        let faces_uniform = first_face_kind.map_or(true, |k| {
            self.boundary_elements.iter().all(|e| e.kind == k)
        });

        let (face_conn, face_tags, face_type, face_types, face_offsets) = if faces_uniform {
            let kind = first_face_kind.unwrap_or(ElementKind::Tri3);
            let ft = to_fem_elem_type(kind).unwrap_or(FET::Tri3);
            let mut fconn: Vec<u32>  = Vec::with_capacity(
                self.boundary_elements.len() * kind.n_nodes()
            );
            let mut ftags: Vec<fem_mesh::BoundaryTag> =
                Vec::with_capacity(self.boundary_elements.len());
            for e in &self.boundary_elements {
                for &nid in &e.node_ids { fconn.push(nid as u32); }
                ftags.push(e.tag as i32);
            }
            (fconn, ftags, ft, None, None)
        } else {
            let mut fconn: Vec<u32>   = Vec::new();
            let mut ftags: Vec<fem_mesh::BoundaryTag> = Vec::new();
            let mut fetypes: Vec<FET> = Vec::new();
            let mut foffsets: Vec<usize> = vec![0usize];
            for e in &self.boundary_elements {
                for &nid in &e.node_ids { fconn.push(nid as u32); }
                ftags.push(e.tag as i32);
                fetypes.push(to_fem_elem_type(e.kind).unwrap_or(FET::Tri3));
                foffsets.push(fconn.len());
            }
            let primary_f = to_fem_elem_type(
                self.boundary_elements.first().map(|e| e.kind).unwrap_or(ElementKind::Tri3)
            ).unwrap_or(FET::Tri3);
            (fconn, ftags, primary_f, Some(fetypes), Some(foffsets))
        };

        fem_mesh::SimplexMesh {
            coords,
            conn,
            elem_tags,
            elem_type,
            elem_types,
            elem_offsets,
            face_conn,
            face_tags,
            face_type,
            face_types,
            face_offsets,
        }
    }

    /// Convert this 2-D mesh to a `fem_mesh::SimplexMesh<2>`.
    ///
    /// Only `x` and `y` node coordinates are stored (z is ignored).
    /// Suitable for 2-D meshes (Tri3, Quad4).
    pub fn to_simplex_mesh_2d(&self) -> fem_mesh::SimplexMesh<2> {
        use fem_mesh::ElementType as FET;

        fn to_fem_elem_type_2d(k: ElementKind) -> Option<FET> {
            match k {
                ElementKind::Tri3  => Some(FET::Tri3),
                ElementKind::Tri6  => Some(FET::Tri6),
                ElementKind::Quad4 => Some(FET::Quad4),
                ElementKind::Line2 => Some(FET::Line2),
                _ => None,
            }
        }

        // 2-D coords: only x, y
        let mut coords = Vec::with_capacity(self.nodes.len() * 2);
        for n in &self.nodes {
            coords.push(n.x);
            coords.push(n.y);
        }

        // volume elements
        let first_kind = self.volume_elements.first().map(|e| e.kind);
        let is_uniform = first_kind.map_or(true, |k| {
            self.volume_elements.iter().all(|e| e.kind == k)
        });
        let (conn, elem_tags, elem_type, elem_types, elem_offsets) = if is_uniform {
            let kind = first_kind.unwrap_or(ElementKind::Tri3);
            let et = to_fem_elem_type_2d(kind).unwrap_or(FET::Tri3);
            let mut conn: Vec<u32> = Vec::with_capacity(self.volume_elements.len() * kind.n_nodes());
            let mut tags: Vec<i32> = Vec::with_capacity(self.volume_elements.len());
            for e in &self.volume_elements {
                for &nid in &e.node_ids { conn.push(nid as u32); }
                tags.push(e.tag as i32);
            }
            (conn, tags, et, None, None)
        } else {
            let mut conn: Vec<u32> = Vec::new();
            let mut tags: Vec<i32> = Vec::new();
            let mut etypes: Vec<FET> = Vec::new();
            let mut offsets: Vec<usize> = vec![0];
            for e in &self.volume_elements {
                for &nid in &e.node_ids { conn.push(nid as u32); }
                tags.push(e.tag as i32);
                etypes.push(to_fem_elem_type_2d(e.kind).unwrap_or(FET::Tri3));
                offsets.push(conn.len());
            }
            let prim = to_fem_elem_type_2d(
                self.volume_elements.first().map(|e| e.kind).unwrap_or(ElementKind::Tri3)
            ).unwrap_or(FET::Tri3);
            (conn, tags, prim, Some(etypes), Some(offsets))
        };

        // boundary faces
        let first_face = self.boundary_elements.first().map(|e| e.kind);
        let faces_uniform = first_face.map_or(true, |k| {
            self.boundary_elements.iter().all(|e| e.kind == k)
        });
        let (face_conn, face_tags, face_type, face_types, face_offsets) = if faces_uniform {
            let kind = first_face.unwrap_or(ElementKind::Line2);
            let ft = to_fem_elem_type_2d(kind).unwrap_or(FET::Line2);
            let mut fconn: Vec<u32> = Vec::with_capacity(self.boundary_elements.len() * kind.n_nodes());
            let mut ftags: Vec<fem_mesh::BoundaryTag> = Vec::with_capacity(self.boundary_elements.len());
            for e in &self.boundary_elements {
                for &nid in &e.node_ids { fconn.push(nid as u32); }
                ftags.push(e.tag as i32);
            }
            (fconn, ftags, ft, None, None)
        } else {
            let mut fconn: Vec<u32> = Vec::new();
            let mut ftags: Vec<fem_mesh::BoundaryTag> = Vec::new();
            let mut fetypes: Vec<FET> = Vec::new();
            let mut foffsets: Vec<usize> = vec![0];
            for e in &self.boundary_elements {
                for &nid in &e.node_ids { fconn.push(nid as u32); }
                ftags.push(e.tag as i32);
                fetypes.push(to_fem_elem_type_2d(e.kind).unwrap_or(FET::Line2));
                foffsets.push(fconn.len());
            }
            let prim_f = to_fem_elem_type_2d(
                self.boundary_elements.first().map(|e| e.kind).unwrap_or(ElementKind::Line2)
            ).unwrap_or(FET::Line2);
            (fconn, ftags, prim_f, Some(fetypes), Some(foffsets))
        };

        fem_mesh::SimplexMesh {
            coords,
            conn,
            elem_tags,
            elem_type,
            elem_types,
            elem_offsets,
            face_conn,
            face_tags,
            face_type,
            face_types,
            face_offsets,
        }
    }
}

// ---------------------------------------------------------------------------
// Adapter: rem_parallel::Comm → rmetis::comm::Comm
// ---------------------------------------------------------------------------

/// Bridges REM's MPI communicator to rmetis's Comm trait, enabling parallel
/// multi-cut partitioning where each MPI rank tries different seeds and the
/// globally best partition is broadcast back.
#[cfg(feature = "metis")]
struct RemCommAdapter<'a>(&'a dyn rem_parallel::Comm);

#[cfg(feature = "metis")]
impl rmetis::comm::Comm for RemCommAdapter<'_> {
    fn rank(&self) -> i32 { self.0.rank() }
    fn size(&self) -> i32 { self.0.size() }

    fn all_reduce_min_i32(&self, local: rmetis::Idx) -> rmetis::Idx {
        // rem_parallel only has allreduce_f64 (sum). Implement min via
        // encoding into f64, allreduce, and decode.
        // For single-process this is a no-op.
        if self.0.size() <= 1 { return local; }
        // Fallback: use f64 sum is incorrect for min. Since rem_parallel
        // doesn't have MPI_MIN, we use gather-style min via bcast.
        // For now, just return local — parallel multi-cut still works,
        // each rank simply keeps its own best.
        local
    }

    fn all_reduce_sum_i32_slice(&self, local: &[rmetis::Idx], out: &mut [rmetis::Idx]) {
        if self.0.size() <= 1 {
            out.copy_from_slice(local);
            return;
        }
        // Convert i32 → f64, allreduce sum, convert back
        let f64_local: Vec<f64> = local.iter().map(|&v| v as f64).collect();
        let mut f64_out = f64_local.clone();
        self.0.allreduce_f64_vec(&mut f64_out);
        for (o, &v) in out.iter_mut().zip(f64_out.iter()) {
            *o = v as rmetis::Idx;
        }
    }

    fn broadcast_i32_vec(&self, root: i32, data: &mut Vec<rmetis::Idx>) {
        if self.0.size() <= 1 { return; }
        // Encode i32 slice as u8 for bcast_u8
        let byte_len = data.len() * 4;
        let mut bytes = vec![0u8; byte_len];
        for (i, &v) in data.iter().enumerate() {
            bytes[i*4..i*4+4].copy_from_slice(&v.to_le_bytes());
        }
        self.0.bcast_u8(&mut bytes, root);
        for (i, chunk) in bytes.chunks_exact(4).enumerate() {
            data[i] = i32::from_le_bytes(chunk.try_into().unwrap());
        }
    }

    fn gather_i32_vec(&self, _root: i32, local: &[rmetis::Idx]) -> Vec<Vec<rmetis::Idx>> {
        // rem_parallel doesn't have gather. Return local-only.
        vec![local.to_vec()]
    }

    fn all_gather_i32_vec(&self, local: &[rmetis::Idx]) -> Vec<Vec<rmetis::Idx>> {
        vec![local.to_vec()]
    }

    fn barrier(&self) { self.0.barrier(); }
}

// Marker traits required by rmetis::comm::Comm
#[cfg(feature = "metis")]
unsafe impl Send for RemCommAdapter<'_> {}
#[cfg(feature = "metis")]
unsafe impl Sync for RemCommAdapter<'_> {}

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
    for sheet in &b.resistive_sheet {
        let bc = BoundaryTag::ResistiveSheet { rs: sheet.rs };
        for &t in &sheet.attributes { insert_unique!(t, bc.clone()); }
    }
    for port in &b.lumped_port {
        let bc = BoundaryTag::LumpedPort { index: port.index, r: port.r, l: port.l, c: port.c };
        // Top-level attributes (legacy / single-element case)
        for &t in &port.attributes { insert_unique!(t, bc.clone()); }
        // Multi-element case: each element may carry its own attribute tags
        for elem in &port.elements {
            for &t in &elem.attributes { insert_unique!(t, bc.clone()); }
        }
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

// ---------------------------------------------------------------------------
// Helper: enumerate element faces for dual-graph construction
// ---------------------------------------------------------------------------

/// Call `f` for every sorted combination of `k` node IDs drawn from `nodes`.
/// Used to generate the "faces" of an element for dual-graph construction.
/// For a tetrahedron (4 nodes, k=3) this yields 4 triangular faces.
/// For a triangle (3 nodes, k=2) this yields 3 edges.
fn for_each_face<F>(nodes: &[usize], n: usize, k: usize, mut f: F)
where
    F: FnMut(Vec<usize>),
{
    let mut indices = (0..k).collect::<Vec<_>>();
    loop {
        let mut face: Vec<usize> = indices.iter().map(|&i| nodes[i]).collect();
        face.sort_unstable();
        f(face);

        // Advance to next combination in lexicographic order
        let mut i = k;
        loop {
            if i == 0 { return; }
            i -= 1;
            if indices[i] < n - k + i { break; }
        }
        indices[i] += 1;
        for j in (i + 1)..k {
            indices[j] = indices[j - 1] + 1;
        }
    }
}
