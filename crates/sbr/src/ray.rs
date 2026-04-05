//! Ray data structures for SBR+ solver.
//!
//! A `Ray` carries position, direction, and the electromagnetic field
//! it transports. Each bounce updates direction and field via Fresnel
//! reflection; the `weight` field drives early termination.

use num_complex::Complex64;

// ---------------------------------------------------------------------------
// Core types
// ---------------------------------------------------------------------------

/// A single propagating ray.
#[derive(Debug, Clone)]
pub struct Ray {
    /// Ray origin [m]
    pub origin: [f64; 3],
    /// Unit propagation direction
    pub dir: [f64; 3],
    /// Complex electric field amplitude [V/m] (3 Cartesian components)
    pub e_field: [Complex64; 3],
    /// Complex magnetic field amplitude [A/m] (3 Cartesian components)
    pub h_field: [Complex64; 3],
    /// Current bounce count (0 = incident ray)
    pub bounce: usize,
    /// Energy weight ∈ (0, 1]; used for early termination
    pub weight: f64,
}

impl Ray {
    /// Construct an incident ray with unit weight.
    pub fn new(origin: [f64; 3], dir: [f64; 3], e_field: [Complex64; 3], h_field: [Complex64; 3]) -> Self {
        Self { origin, dir, e_field, h_field, bounce: 0, weight: 1.0 }
    }
}

// ---------------------------------------------------------------------------
// Hit record
// ---------------------------------------------------------------------------

/// Result of a ray–triangle intersection test.
#[derive(Debug, Clone)]
pub struct RayHit {
    /// Parametric distance along ray: r_hit = origin + t * dir
    pub t: f64,
    /// Index into `SurfaceMesh::faces`
    pub face_idx: usize,
    /// Barycentric coordinates (u, v, 1−u−v)
    pub bary: [f64; 3],
    /// World-space hit point [m]
    pub point: [f64; 3],
    /// Outward unit surface normal at hit point
    pub normal: [f64; 3],
}

// ---------------------------------------------------------------------------
// Möller-Trumbore ray–triangle intersection
// ---------------------------------------------------------------------------

const MT_EPS: f64 = 1e-10;

/// Möller-Trumbore intersection test.
///
/// Returns `Some(t, u, v)` if the ray hits the triangle formed by `p0, p1, p2`
/// with `t > t_min`.  The ray hits the front face when `det > 0` (normal
/// points toward origin) — both sides are tested (double-sided).
pub fn ray_triangle(
    origin: &[f64; 3],
    dir: &[f64; 3],
    p0: &[f64; 3],
    p1: &[f64; 3],
    p2: &[f64; 3],
    t_min: f64,
) -> Option<(f64, f64, f64)> {
    let e1 = sub3(p1, p0);
    let e2 = sub3(p2, p0);
    let h = cross3(dir, &e2);
    let det = dot3(&e1, &h);

    if det.abs() < MT_EPS {
        return None; // parallel
    }

    let f = 1.0 / det;
    let s = sub3(origin, p0);
    let u = f * dot3(&s, &h);
    if !(0.0..=1.0).contains(&u) {
        return None;
    }

    let q = cross3(&s, &e1);
    let v = f * dot3(dir, &q);
    if v < 0.0 || u + v > 1.0 {
        return None;
    }

    let t = f * dot3(&e2, &q);
    if t > t_min { Some((t, u, v)) } else { None }
}

// ---------------------------------------------------------------------------
// Vector utilities (inline, no nalgebra dependency in this module)
// ---------------------------------------------------------------------------

#[inline]
pub fn dot3(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
}

#[inline]
pub fn cross3(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0],
    ]
}

#[inline]
pub fn sub3(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [a[0]-b[0], a[1]-b[1], a[2]-b[2]]
}

#[inline]
pub fn add3(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [a[0]+b[0], a[1]+b[1], a[2]+b[2]]
}

#[inline]
pub fn scale3(a: &[f64; 3], s: f64) -> [f64; 3] {
    [a[0]*s, a[1]*s, a[2]*s]
}

#[inline]
pub fn normalize3(v: [f64; 3]) -> [f64; 3] {
    let len = (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]).sqrt();
    if len < 1e-14 { [1.0, 0.0, 0.0] } else { [v[0]/len, v[1]/len, v[2]/len] }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hit_unit_triangle() {
        let o = [0.1, 0.1, 1.0];
        let d = [0.0, 0.0, -1.0];
        let p0 = [0.0, 0.0, 0.0];
        let p1 = [1.0, 0.0, 0.0];
        let p2 = [0.0, 1.0, 0.0];
        let hit = ray_triangle(&o, &d, &p0, &p1, &p2, 0.0);
        assert!(hit.is_some());
        let (t, u, v) = hit.unwrap();
        assert!((t - 1.0).abs() < 1e-12);
        assert!((u - 0.1).abs() < 1e-12);
        assert!((v - 0.1).abs() < 1e-12);
    }

    #[test]
    fn miss_parallel_ray() {
        let o = [0.5, 0.5, 1.0];
        let d = [1.0, 0.0, 0.0]; // parallel to triangle plane
        let p0 = [0.0, 0.0, 0.0];
        let p1 = [1.0, 0.0, 0.0];
        let p2 = [0.0, 1.0, 0.0];
        assert!(ray_triangle(&o, &d, &p0, &p1, &p2, 0.0).is_none());
    }

    #[test]
    fn miss_behind_origin() {
        let o = [0.1, 0.1, -1.0];
        let d = [0.0, 0.0, -1.0]; // pointing away
        let p0 = [0.0, 0.0, 0.0];
        let p1 = [1.0, 0.0, 0.0];
        let p2 = [0.0, 1.0, 0.0];
        assert!(ray_triangle(&o, &d, &p0, &p1, &p2, 0.0).is_none());
    }
}
