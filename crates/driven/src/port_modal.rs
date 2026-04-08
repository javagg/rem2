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

use rem_mesh::{RemMesh, BoundaryTag, ElementKind};
use rem_core::TripletMatrix;
use nalgebra::DMatrix;
use std::collections::HashMap;

const MU0: f64 = 1.256_637_061_4e-6; // H/m
const C0:  f64 = 2.997_924_58e8;     // m/s

/// Result of a wave-port cross-section eigenvalue solve.
#[derive(Debug, Clone)]
pub struct PortMode {
    /// Cutoff wavenumber k_c = √λ₁  [rad/m]
    pub kc: f64,
    /// Mode shape: global node index → normalised excitation value ∈ [-1, 1].
    /// Only WavePort nodes appear; all others default to 0.
    pub shape: HashMap<usize, f64>,
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

    /// `true` when the port is above cutoff at `freq_hz`.
    pub fn is_propagating(&self, freq_hz: f64) -> bool {
        let k = 2.0 * std::f64::consts::PI * freq_hz / C0;
        k > self.kc
    }
}

/// Compute the fundamental TE mode for `port_idx` on `mesh`.
///
/// Returns `None` when:
/// - No Line2 elements are found for this port tag.
/// - The cross-section has fewer than 3 free DOFs (degenerate geometry).
pub fn compute_wave_port_mode(mesh: &RemMesh, port_idx: u32) -> Option<PortMode> {
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

    let (lambda1, y1) = pairs.into_iter().next()?;
    let kc = lambda1.sqrt();

    // 9. Back-transform to get x = L^{-T} y
    let y1_vec = nalgebra::DVector::from_vec(y1);
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
        "WavePort {port_idx}: k_c = {kc:.4e} rad/m, f_cutoff = {:.4e} Hz, {} free DOFs",
        kc * C0 / (2.0 * std::f64::consts::PI),
        n_free
    );

    Some(PortMode { kc, shape: shape_map })
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
}
