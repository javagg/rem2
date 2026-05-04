//! MoM lumped-port model for S-parameter extraction.
//!
//! A lumped port is defined over a set of surface triangles (identified by
//! boundary attribute IDs).  The port produces:
//! - An excitation RHS vector (one solve per active port).
//! - V/I extraction from the solved surface current coefficients.

use crate::surface_mesh::SurfaceMesh;
use crate::basis::rwg::RwgBasis;
use nalgebra::{DMatrix, DVector};
use num_complex::Complex64;
use rem_core::RemResult;
use std::collections::{HashMap, HashSet};

/// Lumped port: a set of RWG indices + excitation direction + reference Z₀.
#[derive(Debug, Clone)]
pub struct MomLumpedPort {
    /// 1-based port index (matches Touchstone ordering).
    pub index: u32,
    /// Indices into the RWG basis array that are on this port's surface.
    pub rwg_indices: Vec<usize>,
    /// Dominant field direction unit vector [x, y, z].
    pub direction: [f64; 3],
    /// Optional modal profile sampled at global mesh nodes (WavePort mode >= 1).
    pub modal_profile: Option<HashMap<usize, f64>>,
    /// Reference impedance Z₀ [Ω].
    pub z0: f64,
}

impl MomLumpedPort {
    /// Find the RWG basis functions whose shared edge lies on a face that
    /// belongs to one of `port_attrs`.  A face "belongs to a port" if its
    /// attribute tag (stored in `SurfaceMesh::face_attrs`) matches.
    ///
    /// `direction_str` — "x", "y", or "z"; mapped to unit vector.
    pub fn from_surface(
        surf: &SurfaceMesh,
        bases: &[RwgBasis],
        port_attrs: &[u32],
        index: u32,
        direction_str: &str,
        port_type: &str,
        mode: u32,
        z0: f64,
    ) -> RemResult<Self> {
        let port_attr_set: HashSet<u32> = port_attrs.iter().copied().collect();

        let rwg_indices: Vec<usize> = bases.iter().enumerate()
            .filter(|(_, b)| {
                let plus_attr  = surf.face_attrs.get(b.plus_face).copied().unwrap_or(0);
                let minus_attr = surf.face_attrs.get(b.minus_face).copied().unwrap_or(0);
                port_attr_set.contains(&plus_attr) && port_attr_set.contains(&minus_attr)
            })
            .map(|(i, _)| i)
            .collect();

        let modal_profile = if port_type.eq_ignore_ascii_case("waveport") {
            build_waveport_mode_profile(surf, &port_attr_set, mode)
        } else {
            None
        };

        Ok(Self {
            index,
            rwg_indices,
            direction: direction_vec(direction_str),
            modal_profile,
            z0,
        })
    }

    /// Build the N-element excitation RHS for this port (one entry per RWG).
    ///
    /// For RWG basis, the port contribution at basis m is:
    ///   V_m = -∫_{port} f_m(r) · d̂ dS   (d̂ = direction unit vector)
    ///
    /// Only the `rwg_indices` of this port get non-zero entries; all others are 0.
    /// The excitation amplitude `v0` is nominally 1 V.
    pub fn excitation_rhs(
        &self,
        surf: &SurfaceMesh,
        bases: &[RwgBasis],
        n_total: usize,
        v0: Complex64,
    ) -> Vec<Complex64> {
        let mut rhs = vec![Complex64::ZERO; n_total];
        let d = self.direction;
        for &mi in &self.rwg_indices {
            let b = &bases[mi];
            let mut val = Complex64::ZERO;
            // Integrate f_m · d̂ over both support triangles using centroid quadrature.
            // For WavePort mode n, each face contribution is weighted by the
            // average nodal modal amplitude on that triangle.
            for &(face_idx, in_plus) in &[(b.plus_face, true), (b.minus_face, false)] {
                let face = &surf.faces[face_idx];
                let fm = b.eval(&face.centroid, surf, in_plus);
                let dot = d[0]*fm[0] + d[1]*fm[1] + d[2]*fm[2];
                let mode_scale = if let Some(profile) = &self.modal_profile {
                    let n0 = face.nodes[0];
                    let n1 = face.nodes[1];
                    let n2 = face.nodes[2];
                    let v0 = profile.get(&n0).copied().unwrap_or(0.0);
                    let v1 = profile.get(&n1).copied().unwrap_or(0.0);
                    let v2 = profile.get(&n2).copied().unwrap_or(0.0);
                    (v0 + v1 + v2) / 3.0
                } else {
                    1.0
                };
                val += Complex64::new(dot * face.area * mode_scale, 0.0);
            }
            rhs[mi] = -v0 * val;
        }
        rhs
    }

    /// Extract port current I from solved surface current coefficients.
    ///
    /// I ≈ Σ_{m ∈ port} a_m * (div_m_plus * A_plus + div_m_minus * A_minus)
    ///
    /// where a_m are the solved RWG coefficients.
    pub fn extract_current(
        &self,
        surf: &SurfaceMesh,
        bases: &[RwgBasis],
        coeffs: &[Complex64],
    ) -> Complex64 {
        let mut i_port = Complex64::ZERO;
        for &mi in &self.rwg_indices {
            if mi >= coeffs.len() { continue; }
            let b = &bases[mi];
            let div_p = b.divergence(surf, true);
            let div_m = b.divergence(surf, false);
            let area_p = surf.faces[b.plus_face].area;
            let area_m = surf.faces[b.minus_face].area;
            let contrib = div_p * area_p + div_m * area_m;
            i_port += coeffs[mi] * contrib;
        }
        i_port
    }
}

fn direction_vec(s: &str) -> [f64; 3] {
    match s.to_lowercase().as_str() {
        "y" => [0.0, 1.0, 0.0],
        "z" => [0.0, 0.0, 1.0],
        _   => [1.0, 0.0, 0.0],  // default "x"
    }
}

fn edge_key(i: usize, j: usize) -> (usize, usize) {
    if i < j { (i, j) } else { (j, i) }
}

fn build_waveport_mode_profile(
    surf: &SurfaceMesh,
    port_attr_set: &HashSet<u32>,
    mode: u32,
) -> Option<HashMap<usize, f64>> {
    let mode = mode.max(1) as usize;

    let port_faces: Vec<usize> = surf.face_attrs.iter().enumerate()
        .filter_map(|(fi, &attr)| if port_attr_set.contains(&attr) { Some(fi) } else { None })
        .collect();
    if port_faces.is_empty() {
        return None;
    }

    let mut port_nodes = HashSet::<usize>::new();
    let mut edge_count = HashMap::<(usize, usize), usize>::new();
    let mut edge_len = HashMap::<(usize, usize), f64>::new();

    for &fi in &port_faces {
        let tri = &surf.faces[fi];
        let ns = tri.nodes;
        for &n in &ns { port_nodes.insert(n); }
        for &(a, b) in &[(ns[0], ns[1]), (ns[1], ns[2]), (ns[2], ns[0])] {
            let key = edge_key(a, b);
            *edge_count.entry(key).or_insert(0) += 1;
            edge_len.entry(key).or_insert_with(|| {
                let pa = surf.nodes[a];
                let pb = surf.nodes[b];
                let dx = pb[0] - pa[0];
                let dy = pb[1] - pa[1];
                let dz = pb[2] - pa[2];
                (dx * dx + dy * dy + dz * dz).sqrt().max(1.0e-12)
            });
        }
    }

    let mut local_to_global: Vec<usize> = port_nodes.iter().copied().collect();
    local_to_global.sort_unstable();
    let global_to_local: HashMap<usize, usize> = local_to_global.iter().enumerate()
        .map(|(li, &gi)| (gi, li))
        .collect();

    let n_local = local_to_global.len();
    if n_local < 3 {
        return None;
    }

    let mut constrained = vec![false; n_local];
    for (edge, &cnt) in &edge_count {
        if cnt == 1 {
            let (a, b) = *edge;
            if let Some(&la) = global_to_local.get(&a) { constrained[la] = true; }
            if let Some(&lb) = global_to_local.get(&b) { constrained[lb] = true; }
        }
    }

    let free_dofs: Vec<usize> = (0..n_local).filter(|&i| !constrained[i]).collect();
    if free_dofs.len() < 2 {
        return None;
    }

    let free_map: HashMap<usize, usize> = free_dofs.iter().enumerate().map(|(i, &g)| (g, i)).collect();

    let mut k = DMatrix::<f64>::zeros(n_local, n_local);
    let mut m = DMatrix::<f64>::zeros(n_local, n_local);
    for (edge, &len) in &edge_len {
        let (a, b) = *edge;
        let la = global_to_local[&a];
        let lb = global_to_local[&b];

        let inv_l = 1.0 / len;
        k[(la, la)] += inv_l;
        k[(la, lb)] -= inv_l;
        k[(lb, la)] -= inv_l;
        k[(lb, lb)] += inv_l;

        let l3 = len / 3.0;
        let l6 = len / 6.0;
        m[(la, la)] += l3;
        m[(la, lb)] += l6;
        m[(lb, la)] += l6;
        m[(lb, lb)] += l3;
    }

    let n_free = free_dofs.len();
    let mut kf = DMatrix::<f64>::zeros(n_free, n_free);
    let mut mf = DMatrix::<f64>::zeros(n_free, n_free);
    for i in 0..n_local {
        let Some(&ri) = free_map.get(&i) else { continue };
        for j in 0..n_local {
            let Some(&rj) = free_map.get(&j) else { continue };
            kf[(ri, rj)] = k[(i, j)];
            mf[(ri, rj)] = m[(i, j)];
        }
    }

    let chol = mf.clone().cholesky()?;
    let l = chol.l();
    let linv_k = l.solve_lower_triangular(&kf)?;
    let lt = l.transpose();
    let inv_lt = lt.solve_upper_triangular(&DMatrix::<f64>::identity(n_free, n_free))?;
    let a_sym = &linv_k * inv_lt;

    let eig = a_sym.symmetric_eigen();
    let mut pairs: Vec<(f64, Vec<f64>)> = eig.eigenvalues.iter().copied()
        .zip(eig.eigenvectors.column_iter().map(|c| c.iter().copied().collect::<Vec<_>>()))
        .filter(|(lam, _)| *lam > 1.0e-12)
        .collect();
    pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    if pairs.len() < mode {
        return None;
    }

    let (_, y_n) = pairs.remove(mode - 1);
    let y = DVector::from_vec(y_n);
    let x_free = lt.solve_upper_triangular(&y)?;

    let mut x = vec![0.0_f64; n_local];
    for (fi, &li) in free_dofs.iter().enumerate() {
        x[li] = x_free[fi];
    }
    let max_abs = x.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
    if max_abs < 1.0e-12 {
        return None;
    }

    let mut shape = HashMap::<usize, f64>::new();
    for (li, &gi) in local_to_global.iter().enumerate() {
        let v = x[li] / max_abs;
        if v.abs() > 1.0e-8 {
            shape.insert(gi, v);
        }
    }
    Some(shape)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::surface_mesh::{SurfaceMesh, TriFace, SharedEdge, tri_geometry};
    use crate::basis::rwg::generate_rwg_bases;

    /// Build a two-triangle surface with attribute tags.
    fn two_tri_port_surf() -> SurfaceMesh {
        let nodes = vec![
            [0.0_f64, 0.0, 0.0],
            [1.0,     0.0, 0.0],
            [0.5,     1.0, 0.0],
            [-0.5,    1.0, 0.0],
        ];
        let (c0, n0, a0) = tri_geometry(&nodes[0], &nodes[1], &nodes[2]);
        let (c1, n1, a1) = tri_geometry(&nodes[0], &nodes[2], &nodes[3]);
        let faces = vec![
            TriFace { nodes: [0,1,2], centroid: c0, normal: n0, area: a0 },
            TriFace { nodes: [0,2,3], centroid: c1, normal: n1, area: a1 },
        ];
        let edges = vec![SharedEdge {
            nodes: [0, 2],
            plus_face: 0,
            minus_face: 1,
            length: (nodes[2][0].powi(2) + nodes[2][1].powi(2)).sqrt(),
        }];
        SurfaceMesh {
            nodes,
            faces,
            edges,
            boundary_edges: vec![[0,1],[1,2],[2,3],[3,0]],
            face_attrs: vec![1, 1],   // both faces tagged as attr 1
            global_node_ids: vec![],
        }
    }

    #[test]
    fn from_surface_finds_rwg_on_attr() {
        let surf = two_tri_port_surf();
        let bases = generate_rwg_bases(&surf);
        assert_eq!(bases.len(), 1, "should have 1 shared edge");
        let port = MomLumpedPort::from_surface(
            &surf, &bases, &[1], 1, "x", "Lumped", 1, 50.0
        ).expect("port construction failed");
        assert_eq!(port.rwg_indices.len(), 1, "all RWG on attr-1 surface should be found");
        assert_eq!(port.rwg_indices[0], 0);
        assert_eq!(port.z0, 50.0);
    }

    #[test]
    fn from_surface_empty_when_no_match() {
        let surf = two_tri_port_surf();
        let bases = generate_rwg_bases(&surf);
        let port = MomLumpedPort::from_surface(
            &surf, &bases, &[99], 1, "x", "Lumped", 1, 50.0
        ).expect("should succeed even with no match");
        assert!(port.rwg_indices.is_empty());
    }

    #[test]
    fn excitation_rhs_zero_outside_port() {
        let surf = two_tri_port_surf();
        let bases = generate_rwg_bases(&surf);
        let port = MomLumpedPort::from_surface(&surf, &bases, &[1], 1, "x", "Lumped", 1, 50.0).unwrap();
        let rhs = port.excitation_rhs(&surf, &bases, bases.len(), Complex64::new(1.0, 0.0));
        assert_eq!(rhs.len(), bases.len());
        // All entries on the port should be non-zero (there's 1 RWG on port attr 1)
        assert!(rhs[0].norm() > 0.0, "port RWG should have non-zero excitation");
    }

    #[test]
    fn extract_current_finite() {
        let surf = two_tri_port_surf();
        let bases = generate_rwg_bases(&surf);
        let port = MomLumpedPort::from_surface(&surf, &bases, &[1], 1, "x", "Lumped", 1, 50.0).unwrap();
        let coeffs = vec![Complex64::new(1.0, 0.5); bases.len()];
        let i = port.extract_current(&surf, &bases, &coeffs);
        assert!(i.re.is_finite() && i.im.is_finite());
    }
}
