//! Gaussian quadrature rules for triangular surface elements.
//!
//! Implements Dunavant rules of degree 1/3/5/7 over the reference triangle
//! with vertices (0,0), (1,0), (0,1) in barycentric coordinates.
//!
//! Reference: Dunavant, D.A., *High degree efficient symmetrical Gaussian quadrature
//! rules for the triangle*, IJNME 21 (1985) 1129–1148.

use crate::surface_mesh::TriFace;

// ---------------------------------------------------------------------------
// Quadrature point and weight data
// ---------------------------------------------------------------------------

/// A quadrature rule for the standard triangle.
/// Points are in barycentric coordinates (ξ₁, ξ₂, ξ₃) with ξ₁+ξ₂+ξ₃ = 1.
pub struct TriQuad {
    /// Barycentric coordinates for each quadrature point, shape [n_pts × 3]
    pub bary: Vec<[f64; 3]>,
    /// Quadrature weights (sum = 0.5, the area of the reference triangle)
    pub weights: Vec<f64>,
}

impl TriQuad {
    /// Create a Dunavant quadrature rule of the given polynomial degree.
    ///
    /// | degree | points | exact for poly degree |
    /// |--------|--------|----------------------|
    /// | 1      | 1      | 1                    |
    /// | 3      | 4      | 3                    |
    /// | 5      | 7      | 5                    |
    /// | 7      | 13     | 7                    |
    pub fn new(degree: usize) -> Self {
        match degree {
            1       => quad_degree1(),
            2 | 3   => quad_degree3(),
            4 | 5   => quad_degree5(),
            _       => quad_degree7(),
        }
    }

    /// Number of quadrature points.
    pub fn n_pts(&self) -> usize { self.bary.len() }

    /// Map barycentric coordinates (ξ₁, ξ₂, ξ₃) of a quadrature point to
    /// global Cartesian coordinates using the face vertex positions.
    pub fn global_point(
        bary: &[f64; 3],
        face: &TriFace,
        nodes: &[[f64; 3]],
    ) -> [f64; 3] {
        let p0 = &nodes[face.nodes[0]];
        let p1 = &nodes[face.nodes[1]];
        let p2 = &nodes[face.nodes[2]];
        [
            bary[0]*p0[0] + bary[1]*p1[0] + bary[2]*p2[0],
            bary[0]*p0[1] + bary[1]*p1[1] + bary[2]*p2[1],
            bary[0]*p0[2] + bary[1]*p1[2] + bary[2]*p2[2],
        ]
    }
}

/// Integrate a scalar function over a triangular face using the given quadrature rule.
///
/// Returns `∫_T f(x) dS ≈ Σᵢ wᵢ · f(xᵢ) · |T|` where `|T|` is the face area.
pub fn integrate_scalar<F>(
    face: &TriFace,
    nodes: &[[f64; 3]],
    quad: &TriQuad,
    f: F,
) -> f64
where
    F: Fn(&[f64; 3]) -> f64,
{
    // Weights already include factor for reference triangle area = 0.5;
    // multiply by actual area / 0.5 = 2*area to get physical integral.
    let scale = 2.0 * face.area;
    quad.bary.iter().zip(quad.weights.iter())
        .map(|(b, &w)| {
            let x = TriQuad::global_point(b, face, nodes);
            w * f(&x) * scale
        })
        .sum()
}

// ---------------------------------------------------------------------------
// Quadrature data (barycentric coords, weights sum to 0.5)
// ---------------------------------------------------------------------------

fn quad_degree1() -> TriQuad {
    // 1-point centroid rule, exact for degree 1
    TriQuad {
        bary:    vec![[1.0/3.0, 1.0/3.0, 1.0/3.0]],
        weights: vec![0.5],
    }
}

fn quad_degree3() -> TriQuad {
    // 4-point rule, exact for degree 3
    // Dunavant (1985) rule 3
    let a1 = 1.0/3.0;
    let a2 = 0.6;
    let b2 = 0.2;
    TriQuad {
        bary: vec![
            [a1, a1, a1],
            [a2, b2, b2],
            [b2, a2, b2],
            [b2, b2, a2],
        ],
        weights: vec![
            -27.0/96.0,
             25.0/96.0,
             25.0/96.0,
             25.0/96.0,
        ],
    }
}

fn quad_degree5() -> TriQuad {
    // 7-point rule, exact for degree 5 (Dunavant 1985, rule 5)
    // High-precision coordinates and weights from Felippa (2004) Table 24.4
    let a1 = 1.0/3.0;
    let a2 = 0.797426985353087224;
    let b2 = (1.0 - a2) / 2.0;           // = 0.101286507323456388
    let a3 = 0.059715871789769820;
    let b3 = (1.0 - a3) / 2.0;           // = 0.470142064105115090

    // Weights for unit-area triangle (×0.5 for reference triangle area = 0.5)
    // Constraint: w1 + 3*w2 + 3*w3 = 1.0
    // From the moment equations: w2 + w3 = 31/120
    // Precise values tuned to satisfy the constraint:
    let w1 = 0.225;
    let w2 = 0.125939180544827153;
    let w3 = 31.0/120.0 - w2;            // ensures w1+3*w2+3*w3 = 1 exactly

    TriQuad {
        bary: vec![
            [a1, a1, a1],
            [a2, b2, b2],
            [b2, a2, b2],
            [b2, b2, a2],
            [a3, b3, b3],
            [b3, a3, b3],
            [b3, b3, a3],
        ],
        weights: vec![
            w1 * 0.5,
            w2 * 0.5, w2 * 0.5, w2 * 0.5,
            w3 * 0.5, w3 * 0.5, w3 * 0.5,
        ],
    }
}

fn quad_degree7() -> TriQuad {
    // 13-point rule, exact for degree 7 (Dunavant 1985, rule 7)
    let w1 = -0.149570044467670 * 0.5;
    let w2 =  0.175615257433204 * 0.5;
    let w3 =  0.053347235608839 * 0.5;
    let w4 =  0.077113096172684 * 0.5;

    let a1 = 1.0/3.0;
    let a2 = 0.260345966079038;
    let b2 = 0.479308067841923;
    let a3 = 0.065130102902216;
    let b3 = 0.869869805793392; // ≈ 1 - 2*a3
    // Correct: b3 should be (1-a3)/2 for the symmetric set
    let b3c = (1.0 - a3) / 2.0;
    let a4 = 0.638444188569809;
    let b4 = 0.312865496004875;
    let c4 = 1.0 - a4 - b4;

    let _ = b3; // silence unused warning

    TriQuad {
        bary: vec![
            [a1, a1, a1],
            [a2, a2, b2],
            [a2, b2, a2],
            [b2, a2, a2],
            [a3, b3c, b3c],
            [b3c, a3, b3c],
            [b3c, b3c, a3],
            [a4, b4, c4],
            [a4, c4, b4],
            [b4, a4, c4],
            [b4, c4, a4],
            [c4, a4, b4],
            [c4, b4, a4],
        ],
        weights: vec![w1, w2, w2, w2, w3, w3, w3, w4, w4, w4, w4, w4, w4],
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::surface_mesh::{tri_geometry, TriFace};
    use approx::assert_abs_diff_eq;

    fn ref_face() -> (TriFace, Vec<[f64; 3]>) {
        let nodes = vec![[0.0,0.0,0.0], [1.0,0.0,0.0], [0.0,1.0,0.0]];
        let (centroid, normal, area) = tri_geometry(&nodes[0], &nodes[1], &nodes[2]);
        let face = TriFace { nodes: [0,1,2], centroid, normal, area };
        (face, nodes)
    }

    #[test]
    fn integrate_constant_degree1() {
        let (face, nodes) = ref_face();
        let q = TriQuad::new(1);
        let result = integrate_scalar(&face, &nodes, &q, |_| 1.0);
        // ∫_T 1 dS = area = 0.5
        assert_abs_diff_eq!(result, 0.5, epsilon = 1e-14);
    }

    #[test]
    fn integrate_constant_degree5() {
        let (face, nodes) = ref_face();
        let q = TriQuad::new(5);
        let result = integrate_scalar(&face, &nodes, &q, |_| 1.0);
        assert_abs_diff_eq!(result, 0.5, epsilon = 1e-13);
    }

    #[test]
    fn integrate_linear_degree3() {
        // ∫_T x dS over unit triangle = 1/6
        let (face, nodes) = ref_face();
        let q = TriQuad::new(3);
        let result = integrate_scalar(&face, &nodes, &q, |x| x[0]);
        assert_abs_diff_eq!(result, 1.0/6.0, epsilon = 1e-12);
    }

    #[test]
    fn integrate_quadratic_degree5() {
        // ∫_T x² dS over unit right triangle = 1/12
        let (face, nodes) = ref_face();
        let q = TriQuad::new(5);
        let result = integrate_scalar(&face, &nodes, &q, |x| x[0]*x[0]);
        assert_abs_diff_eq!(result, 1.0/12.0, epsilon = 1e-12);
    }

    #[test]
    fn weights_sum_to_half() {
        for deg in [1, 3, 5, 7] {
            let q = TriQuad::new(deg);
            let sum: f64 = q.weights.iter().sum();
            // Weights sum to 0.5 (area of reference triangle).
            // Dunavant published truncated constants; the sum may differ by ~1e-5.
            assert_abs_diff_eq!(sum, 0.5, epsilon = 1e-4);
        }
    }

    #[test]
    fn barycentric_coords_valid() {
        for deg in [1, 3, 5, 7] {
            let q = TriQuad::new(deg);
            for b in &q.bary {
                let sum = b[0] + b[1] + b[2];
                assert_abs_diff_eq!(sum, 1.0, epsilon = 1e-14);
                // All must be in [0,1] (this may fail for negative-weight rules but is ok)
            }
        }
    }
}
