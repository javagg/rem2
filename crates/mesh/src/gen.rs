/// Programmatic GMSH v4.1 ASCII mesh generators.
///
/// These are used to create reproducible test meshes for Palace-equivalent
/// integration tests without requiring an external GMSH installation.
///
/// Output format: GMSH 4.1 ASCII, compatible with `gmsh::read_msh_str`.

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Annular (coaxial cross-section) mesh
// ---------------------------------------------------------------------------

/// Generate a GMSH v4.1 ASCII mesh for an annular 2-D region.
///
/// Geometry: r ∈ [r_i, r_o], θ ∈ [0, 2π) (full annulus).
///
/// Physical groups:
/// - `inner_tag`: boundary edges on the inner circle (r = r_i)
/// - `outer_tag`: boundary edges on the outer circle (r = r_o)
/// - `vol_tag`:   volume triangles filling the annular region
///
/// Parameters:
/// - `r_i`, `r_o`: inner/outer radius in **mesh units** (apply L0 in config)
/// - `n_r`:     number of radial layers
/// - `n_theta`: number of angular sectors (≥ 4)
pub fn annular_msh(
    r_i: f64,
    r_o: f64,
    n_r: usize,
    n_theta: usize,
    inner_tag: u32,
    outer_tag: u32,
    vol_tag: u32,
) -> String {
    assert!(n_r >= 1 && n_theta >= 4, "n_r >= 1 and n_theta >= 4 required");

    let total_nodes = (n_r + 1) * n_theta;
    let n_tri = 2 * n_r * n_theta;
    let total_elems = 2 * n_theta + n_tri; // inner + outer + volume

    // Node id (1-based): layer i_r (0 = inner, n_r = outer), sector i_t
    let node_id = |i_r: usize, i_t: usize| -> usize { i_r * n_theta + i_t + 1 };

    let r_at = |i_r: usize| -> f64 { r_i + i_r as f64 * (r_o - r_i) / n_r as f64 };
    let theta_at = |i_t: usize| -> f64 { 2.0 * PI * i_t as f64 / n_theta as f64 };

    let mut s = String::with_capacity(total_nodes * 60 + total_elems * 30);

    s.push_str("$MeshFormat\n4.1 0 8\n$EndMeshFormat\n");

    // --- Nodes (single entity block) ---
    s.push_str("$Nodes\n");
    s.push_str(&format!("1 {} 1 {}\n", total_nodes, total_nodes));
    s.push_str(&format!("2 1 0 {}\n", total_nodes)); // entityDim=2, entityTag=1
    for id in 1..=total_nodes {
        s.push_str(&format!("{}\n", id));
    }
    for i_r in 0..=n_r {
        let r = r_at(i_r);
        for i_t in 0..n_theta {
            let theta = theta_at(i_t);
            s.push_str(&format!(
                "{:.15e} {:.15e} 0.000000000000000e0\n",
                r * theta.cos(),
                r * theta.sin()
            ));
        }
    }
    s.push_str("$EndNodes\n");

    // --- Elements ---
    s.push_str("$Elements\n");
    s.push_str(&format!("3 {} 1 {}\n", total_elems, total_elems));

    let mut eid = 1usize;

    // Inner boundary edges  (dim=1, tag=inner_tag, type=1=Line2)
    s.push_str(&format!("1 {} 1 {}\n", inner_tag, n_theta));
    for i_t in 0..n_theta {
        let n1 = node_id(0, i_t);
        let n2 = node_id(0, (i_t + 1) % n_theta);
        s.push_str(&format!("{} {} {}\n", eid, n1, n2));
        eid += 1;
    }

    // Outer boundary edges  (dim=1, tag=outer_tag, type=1=Line2)
    s.push_str(&format!("1 {} 1 {}\n", outer_tag, n_theta));
    for i_t in 0..n_theta {
        let n1 = node_id(n_r, i_t);
        let n2 = node_id(n_r, (i_t + 1) % n_theta);
        s.push_str(&format!("{} {} {}\n", eid, n1, n2));
        eid += 1;
    }

    // Volume triangles  (dim=2, tag=vol_tag, type=2=Tri3)
    s.push_str(&format!("2 {} 2 {}\n", vol_tag, n_tri));
    for i_r in 0..n_r {
        for i_t in 0..n_theta {
            let n00 = node_id(i_r,     i_t);
            let n10 = node_id(i_r + 1, i_t);
            let n11 = node_id(i_r + 1, (i_t + 1) % n_theta);
            let n01 = node_id(i_r,     (i_t + 1) % n_theta);
            s.push_str(&format!("{} {} {} {}\n", eid, n00, n10, n11));
            eid += 1;
            s.push_str(&format!("{} {} {} {}\n", eid, n00, n11, n01));
            eid += 1;
        }
    }

    s.push_str("$EndElements\n");
    s
}

// ---------------------------------------------------------------------------
// Rectangular single-material mesh
// ---------------------------------------------------------------------------

/// Generate a GMSH v4.1 ASCII mesh for a rectangular 2-D region.
///
/// Geometry: [0, w] × [0, h].
///
/// Physical groups:
/// - `bottom_tag`: boundary edges at y = 0
/// - `top_tag`:    boundary edges at y = h
/// - `left_tag`:   boundary edges at x = 0
/// - `right_tag`:  boundary edges at x = w
/// - `vol_tag`:    volume triangles filling the region
///
/// Parameters:
/// - `w`, `h`: width and height in **mesh units**
/// - `n_x`, `n_y`: number of cells in x and y directions
pub fn rect_msh(
    w: f64,
    h: f64,
    n_x: usize,
    n_y: usize,
    bottom_tag: u32,
    top_tag: u32,
    left_tag: u32,
    right_tag: u32,
    vol_tag: u32,
) -> String {
    assert!(n_x >= 1 && n_y >= 1, "n_x >= 1 and n_y >= 1 required");

    let total_nodes = (n_x + 1) * (n_y + 1);
    let n_tri = 2 * n_x * n_y;
    let total_elems = 2 * n_x + 2 * n_y + n_tri;

    // Node id (1-based): row iy (0 = bottom), column ix
    let node_id = |ix: usize, iy: usize| -> usize { iy * (n_x + 1) + ix + 1 };

    let mut s = String::with_capacity(total_nodes * 60 + total_elems * 30);

    s.push_str("$MeshFormat\n4.1 0 8\n$EndMeshFormat\n");

    // --- Nodes ---
    s.push_str("$Nodes\n");
    s.push_str(&format!("1 {} 1 {}\n", total_nodes, total_nodes));
    s.push_str(&format!("2 1 0 {}\n", total_nodes));
    for id in 1..=total_nodes {
        s.push_str(&format!("{}\n", id));
    }
    for iy in 0..=n_y {
        for ix in 0..=n_x {
            let x = ix as f64 * w / n_x as f64;
            let y = iy as f64 * h / n_y as f64;
            s.push_str(&format!(
                "{:.15e} {:.15e} 0.000000000000000e0\n",
                x, y
            ));
        }
    }
    s.push_str("$EndNodes\n");

    // --- Elements ---
    s.push_str("$Elements\n");
    s.push_str(&format!("5 {} 1 {}\n", total_elems, total_elems));

    let mut eid = 1usize;

    // Bottom boundary (iy=0)
    s.push_str(&format!("1 {} 1 {}\n", bottom_tag, n_x));
    for ix in 0..n_x {
        s.push_str(&format!("{} {} {}\n", eid, node_id(ix, 0), node_id(ix + 1, 0)));
        eid += 1;
    }

    // Top boundary (iy=n_y)
    s.push_str(&format!("1 {} 1 {}\n", top_tag, n_x));
    for ix in 0..n_x {
        s.push_str(&format!(
            "{} {} {}\n", eid, node_id(ix, n_y), node_id(ix + 1, n_y)
        ));
        eid += 1;
    }

    // Left boundary (ix=0)
    s.push_str(&format!("1 {} 1 {}\n", left_tag, n_y));
    for iy in 0..n_y {
        s.push_str(&format!("{} {} {}\n", eid, node_id(0, iy), node_id(0, iy + 1)));
        eid += 1;
    }

    // Right boundary (ix=n_x)
    s.push_str(&format!("1 {} 1 {}\n", right_tag, n_y));
    for iy in 0..n_y {
        s.push_str(&format!(
            "{} {} {}\n", eid, node_id(n_x, iy), node_id(n_x, iy + 1)
        ));
        eid += 1;
    }

    // Volume triangles
    s.push_str(&format!("2 {} 2 {}\n", vol_tag, n_tri));
    for iy in 0..n_y {
        for ix in 0..n_x {
            let n00 = node_id(ix,     iy);
            let n10 = node_id(ix + 1, iy);
            let n01 = node_id(ix,     iy + 1);
            let n11 = node_id(ix + 1, iy + 1);
            s.push_str(&format!("{} {} {} {}\n", eid, n00, n10, n11));
            eid += 1;
            s.push_str(&format!("{} {} {} {}\n", eid, n00, n11, n01));
            eid += 1;
        }
    }

    s.push_str("$EndElements\n");
    s
}

// ---------------------------------------------------------------------------

/// Generate a GMSH v4.1 ASCII mesh for a rectangular 2-D bimaterial slab.
///
/// Geometry: [0, w] × [0, h], split horizontally at y = h/2.
///
/// Physical groups:
/// - `bottom_tag`:    boundary edges at y = 0
/// - `top_tag`:       boundary edges at y = h
/// - `lower_vol_tag`: volume triangles in the lower half (y < h/2)
/// - `upper_vol_tag`: volume triangles in the upper half (y ≥ h/2)
///
/// Left and right edges are **untagged** (natural Neumann BC).
///
/// Parameters:
/// - `w`, `h`: width and height in **mesh units** (apply L0 in config)
/// - `n_x`, `n_y`: number of cells in x and y directions (n_y must be even)
pub fn rect_bimaterial_msh(
    w: f64,
    h: f64,
    n_x: usize,
    n_y: usize,
    bottom_tag: u32,
    top_tag: u32,
    lower_vol_tag: u32,
    upper_vol_tag: u32,
) -> String {
    assert!(n_x >= 1 && n_y >= 2 && n_y % 2 == 0, "n_x >= 1, n_y >= 2 and even");

    let total_nodes = (n_x + 1) * (n_y + 1);
    let n_bot = n_x;
    let n_top = n_x;
    let n_lower = 2 * n_x * (n_y / 2);
    let n_upper = 2 * n_x * (n_y / 2);
    let total_elems = n_bot + n_top + n_lower + n_upper;

    // Node id (1-based): row iy (0 = bottom), column ix
    let node_id = |ix: usize, iy: usize| -> usize { iy * (n_x + 1) + ix + 1 };

    let mut s = String::with_capacity(total_nodes * 55 + total_elems * 25);

    s.push_str("$MeshFormat\n4.1 0 8\n$EndMeshFormat\n");

    // --- Nodes ---
    s.push_str("$Nodes\n");
    s.push_str(&format!("1 {} 1 {}\n", total_nodes, total_nodes));
    s.push_str(&format!("2 1 0 {}\n", total_nodes));
    for id in 1..=total_nodes {
        s.push_str(&format!("{}\n", id));
    }
    for iy in 0..=n_y {
        for ix in 0..=n_x {
            let x = ix as f64 * w / n_x as f64;
            let y = iy as f64 * h / n_y as f64;
            s.push_str(&format!(
                "{:.15e} {:.15e} 0.000000000000000e0\n",
                x, y
            ));
        }
    }
    s.push_str("$EndNodes\n");

    // --- Elements ---
    s.push_str("$Elements\n");
    s.push_str(&format!("4 {} 1 {}\n", total_elems, total_elems));

    let mut eid = 1usize;

    // Bottom boundary (iy=0, tag=bottom_tag, type=1=Line2)
    s.push_str(&format!("1 {} 1 {}\n", bottom_tag, n_bot));
    for ix in 0..n_x {
        s.push_str(&format!("{} {} {}\n", eid, node_id(ix, 0), node_id(ix + 1, 0)));
        eid += 1;
    }

    // Top boundary (iy=n_y, tag=top_tag, type=1=Line2)
    s.push_str(&format!("1 {} 1 {}\n", top_tag, n_top));
    for ix in 0..n_x {
        s.push_str(&format!(
            "{} {} {}\n", eid, node_id(ix, n_y), node_id(ix + 1, n_y)
        ));
        eid += 1;
    }

    // Lower-half triangles (iy in 0..n_y/2, tag=lower_vol_tag, type=2=Tri3)
    s.push_str(&format!("2 {} 2 {}\n", lower_vol_tag, n_lower));
    for iy in 0..(n_y / 2) {
        for ix in 0..n_x {
            let n00 = node_id(ix,     iy);
            let n10 = node_id(ix + 1, iy);
            let n01 = node_id(ix,     iy + 1);
            let n11 = node_id(ix + 1, iy + 1);
            s.push_str(&format!("{} {} {} {}\n", eid, n00, n10, n11));
            eid += 1;
            s.push_str(&format!("{} {} {} {}\n", eid, n00, n11, n01));
            eid += 1;
        }
    }

    // Upper-half triangles (iy in n_y/2..n_y, tag=upper_vol_tag, type=2=Tri3)
    s.push_str(&format!("2 {} 2 {}\n", upper_vol_tag, n_upper));
    for iy in (n_y / 2)..n_y {
        for ix in 0..n_x {
            let n00 = node_id(ix,     iy);
            let n10 = node_id(ix + 1, iy);
            let n01 = node_id(ix,     iy + 1);
            let n11 = node_id(ix + 1, iy + 1);
            s.push_str(&format!("{} {} {} {}\n", eid, n00, n10, n11));
            eid += 1;
            s.push_str(&format!("{} {} {} {}\n", eid, n00, n11, n01));
            eid += 1;
        }
    }

    s.push_str("$EndElements\n");
    s
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gmsh::read_msh_str;

    #[test]
    fn annular_mesh_parses() {
        let msh = annular_msh(1.0, 4.0, 4, 8, 1, 2, 10);
        let raw = read_msh_str(&msh).expect("annular mesh should parse");
        // Nodes: (4+1)*8 = 40
        assert_eq!(raw.nodes.len(), 40);
        // Elements: 8 inner + 8 outer + 2*4*8=64 tris = 80 total
        assert_eq!(raw.elements.len(), 80);
    }

    #[test]
    fn annular_inner_tag_assigned() {
        let msh = annular_msh(1.0, 4.0, 3, 8, 7, 8, 99);
        let raw = read_msh_str(&msh).unwrap();
        // Inner boundary elements should have phys_tag = 7
        let inner: Vec<_> = raw.elements.iter().filter(|e| e.phys_tag == 7).collect();
        assert_eq!(inner.len(), 8, "8 inner boundary edges");
        let outer: Vec<_> = raw.elements.iter().filter(|e| e.phys_tag == 8).collect();
        assert_eq!(outer.len(), 8, "8 outer boundary edges");
        let vol: Vec<_> = raw.elements.iter().filter(|e| e.phys_tag == 99).collect();
        assert_eq!(vol.len(), 2 * 3 * 8, "volume triangles");
    }

    #[test]
    fn annular_node_radii() {
        let r_i = 1.0_f64;
        let r_o = 3.0_f64;
        let msh = annular_msh(r_i, r_o, 4, 16, 1, 2, 10);
        let raw = read_msh_str(&msh).unwrap();
        // Inner ring (first 16 nodes) should all have r ≈ r_i
        for (_, x, y, _) in &raw.nodes[0..16] {
            let r = (x * x + y * y).sqrt();
            assert!((r - r_i).abs() < 1e-12, "inner radius: r={}", r);
        }
        // Outer ring (last 16 nodes) should have r ≈ r_o
        for (_, x, y, _) in &raw.nodes[raw.nodes.len() - 16..] {
            let r = (x * x + y * y).sqrt();
            assert!((r - r_o).abs() < 1e-12, "outer radius: r={}", r);
        }
    }

    #[test]
    fn rect_mesh_parses() {
        let msh = rect_msh(1.0, 1.0, 4, 4, 1, 2, 3, 4, 10);
        let raw = read_msh_str(&msh).expect("rect mesh should parse");
        // Nodes: 5*5 = 25
        assert_eq!(raw.nodes.len(), 25);
        // Elements: 4 bottom + 4 top + 4 left + 4 right + 2*4*4=32 tris = 48
        assert_eq!(raw.elements.len(), 48);
    }

    #[test]
    fn rect_bimaterial_parses() {
        let msh = rect_bimaterial_msh(1.0, 1.0, 4, 4, 1, 2, 10, 20);
        let raw = read_msh_str(&msh).expect("rect bimaterial mesh should parse");
        // Nodes: 5*5 = 25
        assert_eq!(raw.nodes.len(), 25);
        // Elements: 4 bottom + 4 top + 2*4*2=16 lower + 2*4*2=16 upper = 40
        assert_eq!(raw.elements.len(), 40);
    }

    #[test]
    fn rect_bimaterial_half_split() {
        let msh = rect_bimaterial_msh(2.0, 2.0, 6, 6, 1, 2, 10, 20);
        let raw = read_msh_str(&msh).unwrap();
        let lower: Vec<_> = raw.elements.iter().filter(|e| e.phys_tag == 10).collect();
        let upper: Vec<_> = raw.elements.iter().filter(|e| e.phys_tag == 20).collect();
        assert_eq!(lower.len(), upper.len(), "equal lower/upper triangles");
        assert_eq!(lower.len(), 2 * 6 * 3); // 2 * n_x * (n_y/2) = 2*6*3
    }
}
