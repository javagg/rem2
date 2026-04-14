//! Wave-port modal analysis — TE/TM 2-D scalar cross-section eigenvalue.
//!
//! For a wave port whose cross-section is a set of **Line2** boundary
//! elements, this module
//!
//! 1. Builds the 1-D stiffness K_p and mass M_p for those line segments.
//! 2. Applies homogeneous Dirichlet at the two endpoints (metal walls).
//! 3. Solves the small tridiagonal eigenvalue problem K_p x = λ M_p x.
//! 4. Returns the fundamental mode shape (eigenvector) and cutoff
//!    wavenumber k_c = √λ₁.
//!
//! The caller can then
//! - use `mode.shape` as the Dirichlet excitation profile (φ_n = shape[n])
//!   instead of the uniform φ=V TEM approximation, and
//! - use `mode.te_impedance(freq_hz)` as the reference impedance Z₀.
//!
//! For a rectangular waveguide of width W the first eigenvalue is
//! λ₁ = (π/W)² and the eigenvector is sin(πx/W), matching TE₁₀.

use rem_mesh::{RemMesh, BoundaryTag, ElementKind, extract_submesh_by_element_ids_tri3};
use rem_core::TripletMatrix;
use nalgebra::DMatrix;
use std::collections::HashMap;

const MU0: f64 = 1.256_637_061_4e-6; // H/m
const C0:  f64 = 2.997_924_58e8;     // m/s

/// Whether the cross-section mode is TE or TM.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ModeType { Te, Tm }

/// Result of a wave-port cross-section eigenvalue solve.
#[derive(Debug, Clone)]
pub struct PortMode {
    /// Cutoff wavenumber k_c = √λ₁  [rad/m]
    pub kc: f64,
    /// Mode shape: global node index → normalised excitation value ∈ [-1, 1].
    /// Only WavePort nodes appear; all others default to 0.
    pub shape: HashMap<usize, f64>,
    /// TE or TM mode classification.
    pub mode_type: ModeType,
}

#[derive(Debug, Clone, PartialEq)]
pub struct PortSupportRegionSummary {
    pub port_index: u32,
    pub n_volume_elements: usize,
    pub n_nodes: usize,
    pub domain_tags: Vec<u32>,
    pub boundary_length: f64,
    pub x_min: f64,
    pub y_min: f64,
    pub x_max: f64,
    pub y_max: f64,
}

impl PortMode {
    /// TE-mode wave impedance Z_TE = ω μ₀ / k_z at `freq_hz` [Ω].
    ///
    /// Returns `f64::INFINITY` when below cutoff (evanescent mode).
    pub fn te_impedance(&self, freq_hz: f64) -> f64 {
        let k = 2.0 * std::f64::consts::PI * freq_hz / C0;
        let kz2 = k * k - self.kc * self.kc;
        if kz2 <= 0.0 {
            return f64::INFINITY; // below cutoff
        }
        let omega = 2.0 * std::f64::consts::PI * freq_hz;
        omega * MU0 / kz2.sqrt()
    }

    /// TM-mode wave impedance Z_TM = k_z / (ω ε₀) at `freq_hz` [Ω].
    ///
    /// Returns `0.0` when below cutoff (evanescent mode).
    pub fn tm_impedance(&self, freq_hz: f64) -> f64 {
        let k = 2.0 * std::f64::consts::PI * freq_hz / C0;
        let kz2 = k * k - self.kc * self.kc;
        if kz2 <= 0.0 {
            return 0.0; // below cutoff — evanescent
        }
        let eps0: f64 = 8.854_187_817e-12;
        let omega = 2.0 * std::f64::consts::PI * freq_hz;
        kz2.sqrt() / (omega * eps0)
    }

    /// Wave impedance appropriate for this mode's type.
    pub fn impedance(&self, freq_hz: f64) -> f64 {
        match self.mode_type {
            ModeType::Te => self.te_impedance(freq_hz),
            ModeType::Tm => self.tm_impedance(freq_hz),
        }
    }

    /// `true` when the port is above cutoff at `freq_hz`.
    pub fn is_propagating(&self, freq_hz: f64) -> bool {
        let k = 2.0 * std::f64::consts::PI * freq_hz / C0;
        k > self.kc
    }
}

/// Compute the wave-port mode for `port_idx` on `mesh`.
///
/// `mode_number` = 1 selects the fundamental (lowest k_c) TE mode.
/// `mode_number` > 1 selects higher-order modes in eigenvalue order.
/// When `mode_number` is 0 it is treated as 1.
///
/// Returns `None` when:
/// - No Line2 elements are found for this port tag.
/// - The cross-section has fewer than 3 free DOFs (degenerate geometry).
pub fn compute_wave_port_mode(mesh: &RemMesh, port_idx: u32) -> Option<PortMode> {
    compute_wave_port_mode_n(mesh, port_idx, 1)
}

/// As `compute_wave_port_mode` but selects the `mode_n`-th eigenvalue (1-based).
pub fn compute_wave_port_mode_n(mesh: &RemMesh, port_idx: u32, mode_n: u32) -> Option<PortMode> {
    let mode_n = mode_n.max(1) as usize;
    // 1. Find the WavePort physical tag for this index
    let port_tag = mesh.boundary_tags.iter()
        .find_map(|(&tag, bc)| {
            if let BoundaryTag::WavePort { index } = bc {
                if *index == port_idx { Some(tag) } else { None }
            } else {
                None
            }
        })?;

    // 2. Collect Line2 elements on this port
    let port_elems: Vec<_> = mesh.boundary_elements.iter()
        .filter(|e| e.tag == port_tag && e.kind == ElementKind::Line2)
        .collect();

    if port_elems.is_empty() {
        return None;
    }

    log_port_support_region(mesh, port_idx, &port_elems);

    // 3. Build a local node numbering (global id → local 0-based)
    let mut global_to_local: HashMap<usize, usize> = HashMap::new();
    let mut local_to_global: Vec<usize> = Vec::new();
    for elem in &port_elems {
        for &nid in &elem.node_ids {
            if !global_to_local.contains_key(&nid) {
                global_to_local.insert(nid, local_to_global.len());
                local_to_global.push(nid);
            }
        }
    }
    let n_local = local_to_global.len();

    // 4. Find endpoint nodes (appear in exactly one Line2 element)
    let mut node_count: HashMap<usize, usize> = HashMap::new();
    for elem in &port_elems {
        for &nid in &elem.node_ids {
            *node_count.entry(nid).or_insert(0) += 1;
        }
    }
    let endpoint_locals: Vec<usize> = node_count.iter()
        .filter_map(|(&nid, &cnt)| if cnt == 1 { global_to_local.get(&nid).copied() } else { None })
        .collect();

    // 5. Assemble 1-D K_p and M_p for line elements
    let mut k_trip = TripletMatrix::with_capacity(n_local, n_local, port_elems.len() * 4);
    let mut m_trip = TripletMatrix::with_capacity(n_local, n_local, port_elems.len() * 4);

    for elem in &port_elems {
        let n0 = global_to_local[&elem.node_ids[0]];
        let n1 = global_to_local[&elem.node_ids[1]];
        let x0 = &mesh.nodes[elem.node_ids[0]];
        let x1 = &mesh.nodes[elem.node_ids[1]];
        let dx = x1.x - x0.x;
        let dy = x1.y - x0.y;
        let dz = x1.z - x0.z;
        let len = (dx * dx + dy * dy + dz * dz).sqrt();
        if len < 1e-300 { continue; }

        // Stiffness: [1/L, -1/L; -1/L, 1/L]
        k_trip.add(n0, n0,  1.0 / len);
        k_trip.add(n0, n1, -1.0 / len);
        k_trip.add(n1, n0, -1.0 / len);
        k_trip.add(n1, n1,  1.0 / len);

        // Consistent mass: [L/3, L/6; L/6, L/3]
        m_trip.add(n0, n0, len / 3.0);
        m_trip.add(n0, n1, len / 6.0);
        m_trip.add(n1, n0, len / 6.0);
        m_trip.add(n1, n1, len / 3.0);
    }

    // 6. Apply Dirichlet at endpoints — symmetric elimination with unit diagonal
    let k_csr = k_trip.to_csr();
    let m_csr = m_trip.to_csr();

    // Determine free DOFs (all except endpoints)
    let is_constrained: Vec<bool> = (0..n_local)
        .map(|i| endpoint_locals.contains(&i))
        .collect();
    let free_dofs: Vec<usize> = (0..n_local).filter(|&i| !is_constrained[i]).collect();
    let n_free = free_dofs.len();

    if n_free < 2 {
        // Too few free DOFs for a meaningful mode — fall back to TEM
        return None;
    }

    let free_map: HashMap<usize, usize> = free_dofs.iter().enumerate().map(|(li, &gi)| (gi, li)).collect();

    // 7. Build dense reduced K̃ and M̃ (free × free sub-matrices)
    let mut kf = DMatrix::<f64>::zeros(n_free, n_free);
    let mut mf = DMatrix::<f64>::zeros(n_free, n_free);

    for i in 0..k_csr.nrows {
        let Some(&ri) = free_map.get(&i) else { continue };
        for ptr in k_csr.row_ptr[i]..k_csr.row_ptr[i+1] {
            let j = k_csr.col_idx[ptr];
            if let Some(&rj) = free_map.get(&j) {
                kf[(ri, rj)] += k_csr.values[ptr];
            }
        }
    }
    for i in 0..m_csr.nrows {
        let Some(&ri) = free_map.get(&i) else { continue };
        for ptr in m_csr.row_ptr[i]..m_csr.row_ptr[i+1] {
            let j = m_csr.col_idx[ptr];
            if let Some(&rj) = free_map.get(&j) {
                mf[(ri, rj)] += m_csr.values[ptr];
            }
        }
    }

    // 8. Solve K̃ x = λ M̃ x — use Cholesky of M̃ to reduce to standard form
    //    L L^T = M̃  →  (L^{-1} K̃ L^{-T}) y = λ y
    let mf_chol = match mf.clone().cholesky() {
        Some(c) => c,
        None => {
            log::warn!("WavePort {port_idx}: M_port Cholesky failed; using TEM approx");
            return None;
        }
    };
    let l = mf_chol.l();
    // Solve L A_red = K̃  →  A_red = L^{-1} K̃
    let linv_kf = l.solve_lower_triangular(&kf).unwrap_or(kf.clone());
    // A_sym = L^{-1} K̃ L^{-T}  (symmetric)
    let lt = l.transpose();
    let a_sym = &linv_kf * lt.solve_upper_triangular(&DMatrix::identity(n_free, n_free)).unwrap_or(DMatrix::identity(n_free, n_free));

    // Symmetric eigendecomposition
    let eig = a_sym.symmetric_eigen();
    let mut pairs: Vec<(f64, Vec<f64>)> = eig.eigenvalues.iter().copied()
        .zip(eig.eigenvectors.column_iter().map(|col| col.iter().copied().collect::<Vec<_>>()))
        .filter(|(lambda, _)| *lambda > 1e-14)
        .collect();
    pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    // Select the mode_n-th mode (1-based index).
    if pairs.len() < mode_n {
        log::warn!("WavePort {port_idx}: only {} positive modes found, requested mode {}", pairs.len(), mode_n);
        return None;
    }
    let (lambda_n, y_n) = pairs.remove(mode_n - 1);
    let kc = lambda_n.sqrt();

    // 9. Back-transform to get x = L^{-T} y
    let y1_vec = nalgebra::DVector::from_vec(y_n);
    let x_free = match lt.solve_upper_triangular(&y1_vec) {
        Some(v) => v,
        None => y1_vec,
    };

    // 10. Map back to global node indices; normalize so max abs = 1
    let mut raw: Vec<f64> = vec![0.0; n_local];
    for (fi, &li) in free_dofs.iter().enumerate() {
        raw[li] = x_free[fi];
    }
    let max_abs = raw.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
    if max_abs < 1e-300 {
        return None;
    }
    let shape_map: HashMap<usize, f64> = local_to_global.iter()
        .enumerate()
        .filter_map(|(li, &gi)| {
            let v = raw[li] / max_abs;
            if v.abs() > 1e-12 { Some((gi, v)) } else { None }
        })
        .collect();

    log::info!(
        "WavePort {port_idx} mode {mode_n}: k_c = {kc:.4e} rad/m, f_cutoff = {:.4e} Hz, {} free DOFs",
        kc * C0 / (2.0 * std::f64::consts::PI),
        n_free
    );

    // All modes solved with Dirichlet BCs at endpoints are TE modes.
    // (TM modes require Neumann BCs and a separate solve — not yet implemented.)
    Some(PortMode { kc, shape: shape_map, mode_type: ModeType::Te })
}

pub fn collect_port_support_region(mesh: &RemMesh, port_idx: u32) -> Option<PortSupportRegionSummary> {
    let port_tag = mesh.boundary_tags.iter()
        .find_map(|(&tag, bc)| match bc {
            BoundaryTag::WavePort { index } if *index == port_idx => Some(tag),
            _ => None,
        })?;
    let port_elems: Vec<_> = mesh.boundary_elements.iter()
        .filter(|element| element.tag == port_tag && element.kind == ElementKind::Line2)
        .collect();
    summarize_port_support_region(mesh, port_idx, &port_elems)
}

fn log_port_support_region(mesh: &RemMesh, port_idx: u32, port_elems: &[&rem_mesh::Element]) {
    if let Some(summary) = summarize_port_support_region(mesh, port_idx, port_elems) {
        log::info!(
            "WavePort {} support region: {} Tri3 elems, {} nodes, domain tags {:?}",
            summary.port_index,
            summary.n_volume_elements,
            summary.n_nodes,
            summary.domain_tags
        );
    }
}

fn summarize_port_support_region(
    mesh: &RemMesh,
    port_idx: u32,
    port_elems: &[&rem_mesh::Element],
) -> Option<PortSupportRegionSummary> {
    if mesh.dim != 2
        || !mesh.volume_elements.iter().all(|element| element.kind == ElementKind::Tri3)
        || !mesh.boundary_elements.iter().all(|element| element.kind == ElementKind::Line2)
    {
        return None;
    }

    let mut adjacent_volume_ids = Vec::new();
    for (elem_id, element) in mesh.volume_elements.iter().enumerate() {
        if element.kind != ElementKind::Tri3 {
            continue;
        }
        let contains_port_edge = port_elems.iter().any(|port_elem| {
            let a = port_elem.node_ids[0];
            let b = port_elem.node_ids[1];
            triangle_has_edge(&element.node_ids, a, b)
        });
        if contains_port_edge {
            adjacent_volume_ids.push(elem_id);
        }
    }

    if adjacent_volume_ids.is_empty() {
        return None;
    }

    match extract_submesh_by_element_ids_tri3(mesh, &adjacent_volume_ids) {
        Ok(submesh) => {
            let mut tags: Vec<u32> = adjacent_volume_ids
                .iter()
                .map(|&element_id| mesh.volume_elements[element_id].tag)
                .collect();
            tags.sort_unstable();
            tags.dedup();
            let boundary_length: f64 = port_elems.iter().map(|element| line_length(mesh, element)).sum();
            let mut unique_nodes: Vec<usize> = port_elems
                .iter()
                .flat_map(|element| element.node_ids.iter().copied())
                .collect();
            unique_nodes.sort_unstable();
            unique_nodes.dedup();
            let x_min = unique_nodes.iter().map(|&node_id| mesh.nodes[node_id].x).fold(f64::INFINITY, f64::min);
            let y_min = unique_nodes.iter().map(|&node_id| mesh.nodes[node_id].y).fold(f64::INFINITY, f64::min);
            let x_max = unique_nodes.iter().map(|&node_id| mesh.nodes[node_id].x).fold(f64::NEG_INFINITY, f64::max);
            let y_max = unique_nodes.iter().map(|&node_id| mesh.nodes[node_id].y).fold(f64::NEG_INFINITY, f64::max);
            Some(PortSupportRegionSummary {
                port_index: port_idx,
                n_volume_elements: submesh.mesh.n_volume_elements(),
                n_nodes: submesh.mesh.n_nodes(),
                domain_tags: tags,
                boundary_length,
                x_min,
                y_min,
                x_max,
                y_max,
            })
        }
        Err(err) => {
            log::warn!(
                "WavePort {} support-region extraction via fem bridge failed ({})",
                port_idx,
                err
            );
            None
        }
    }
}

fn triangle_has_edge(node_ids: &[usize], a: usize, b: usize) -> bool {
    let has_a = node_ids.contains(&a);
    let has_b = node_ids.contains(&b);
    has_a && has_b
}

fn line_length(mesh: &RemMesh, element: &rem_mesh::Element) -> f64 {
    let a = &mesh.nodes[element.node_ids[0]];
    let b = &mesh.nodes[element.node_ids[1]];
    let dx = b.x - a.x;
    let dy = b.y - a.y;
    let dz = b.z - a.z;
    (dx * dx + dy * dy + dz * dz).sqrt()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use rem_mesh::{Node, Element, ElementKind, RemMesh};

    /// Rectangular waveguide cross-section of width W.
    /// 5 nodes, 4 Line2 elements, endpoints constrained.
    /// Expected k_c = π / W, mode shape ≈ sin(πx/W).
    fn rectangular_port_mesh(w: f64, n_seg: usize) -> RemMesh {
        // n_seg equal-length segments across width W
        let dx = w / n_seg as f64;
        let nodes: Vec<Node> = (0..=n_seg)
            .map(|i| Node { id: i, x: i as f64 * dx, y: 0.0, z: 0.0 })
            .collect();
        let boundary_elements: Vec<Element> = (0..n_seg)
            .map(|i| Element {
                id: i + 1,
                kind: ElementKind::Line2,
                tag: 10,
                node_ids: vec![i, i + 1],
                rank: 0,
            })
            .collect();
        let mut boundary_tags = HashMap::new();
        boundary_tags.insert(10u32, BoundaryTag::WavePort { index: 1 });
        RemMesh {
            nodes,
            volume_elements: vec![],
            boundary_elements,
            domain_tags: Default::default(),
            boundary_tags,
            dim: 2,
            rank: 0,
            size: 1,
        }
    }

    fn tri3_support_mesh() -> RemMesh {
        let nodes = vec![
            Node { id: 0, x: 0.0, y: 0.0, z: 0.0 },
            Node { id: 1, x: 1.0, y: 0.0, z: 0.0 },
            Node { id: 2, x: 1.0, y: 1.0, z: 0.0 },
            Node { id: 3, x: 0.0, y: 1.0, z: 0.0 },
        ];
        let volume_elements = vec![
            Element { id: 1, kind: ElementKind::Tri3, tag: 7, node_ids: vec![0, 1, 2], rank: 0 },
            Element { id: 2, kind: ElementKind::Tri3, tag: 8, node_ids: vec![0, 2, 3], rank: 0 },
        ];
        let boundary_elements = vec![
            Element { id: 3, kind: ElementKind::Line2, tag: 10, node_ids: vec![0, 1], rank: 0 },
            Element { id: 4, kind: ElementKind::Line2, tag: 11, node_ids: vec![2, 3], rank: 0 },
        ];
        let mut boundary_tags = HashMap::new();
        boundary_tags.insert(10u32, BoundaryTag::WavePort { index: 1 });
        boundary_tags.insert(11u32, BoundaryTag::Ground);
        RemMesh {
            nodes,
            volume_elements,
            boundary_elements,
            domain_tags: Default::default(),
            boundary_tags,
            dim: 2,
            rank: 0,
            size: 1,
        }
    }

    #[test]
    fn te10_cutoff_wavenumber() {
        let w = 0.1; // 10 cm waveguide, k_c = π/0.1 ≈ 31.416 rad/m
        let mesh = rectangular_port_mesh(w, 20);
        let mode = compute_wave_port_mode(&mesh, 1).expect("mode not found");

        let kc_exact = std::f64::consts::PI / w;
        let rel_err = (mode.kc - kc_exact).abs() / kc_exact;
        assert!(
            rel_err < 0.01,
            "k_c = {:.4e}, expected {:.4e}, rel_err = {:.3e}",
            mode.kc, kc_exact, rel_err
        );
    }

    #[test]
    fn te10_mode_shape_sin() {
        let w = 1.0;
        let n_seg = 40;
        let mesh = rectangular_port_mesh(w, n_seg);
        let mode = compute_wave_port_mode(&mesh, 1).expect("mode not found");

        // Check shape ≈ sin(π x / W) at interior nodes
        let dx = w / n_seg as f64;
        let mut max_err = 0.0_f64;
        for i in 1..n_seg {
            let x = i as f64 * dx;
            let expected = (std::f64::consts::PI * x / w).sin();
            let got = mode.shape.get(&i).copied().unwrap_or(0.0);
            // May be negated — compare |got| vs expected (shape normalised to max=1 = sin(π/2))
            // expected peak at x=0.5 → sin(π/2)=1 → expected=1; got≈1
            let err = (got.abs() - expected.abs()).abs();
            if err > max_err { max_err = err; }
        }
        assert!(max_err < 0.02, "max mode shape error = {max_err:.4e}");
    }

    #[test]
    fn te_impedance_above_cutoff() {
        let w = 0.1_f64; // a-band waveguide
        let mesh = rectangular_port_mesh(w, 20);
        let mode = compute_wave_port_mode(&mesh, 1).unwrap();
        let fc = mode.kc * C0 / (2.0 * std::f64::consts::PI);
        let f_op = 2.0 * fc; // well above cutoff
        let z = mode.te_impedance(f_op);
        // Z_TE should be > 377 Ω at f=2*fc
        assert!(z > 377.0, "Z_TE = {z:.1} Ω");
        // should be finite
        assert!(z.is_finite(), "Z_TE not finite");
    }

    #[test]
    fn te_impedance_below_cutoff() {
        let w = 0.1_f64;
        let mesh = rectangular_port_mesh(w, 20);
        let mode = compute_wave_port_mode(&mesh, 1).unwrap();
        let fc = mode.kc * C0 / (2.0 * std::f64::consts::PI);
        let z = mode.te_impedance(0.5 * fc); // below cutoff → evanescent
        assert!(!z.is_finite(), "expected inf below cutoff, got {z}");
    }

    #[test]
    fn support_region_logging_helper_handles_non_tri3_mesh() {
        let mesh = rectangular_port_mesh(1.0, 4);
        let port_elems: Vec<_> = mesh.boundary_elements.iter().collect();
        log_port_support_region(&mesh, 1, &port_elems);
    }

    #[test]
    fn collect_port_support_region_summarizes_adjacent_triangles() {
        let mesh = tri3_support_mesh();
        let summary = collect_port_support_region(&mesh, 1).expect("support-region summary should exist");
        assert_eq!(summary.port_index, 1);
        assert_eq!(summary.n_volume_elements, 1);
        assert_eq!(summary.n_nodes, 3);
        assert_eq!(summary.domain_tags, vec![7]);
        assert!((summary.boundary_length - 1.0).abs() < 1e-12);
        assert_eq!(summary.x_min, 0.0);
        assert_eq!(summary.y_min, 0.0);
        assert_eq!(summary.x_max, 1.0);
        assert_eq!(summary.y_max, 0.0);
    }
}
