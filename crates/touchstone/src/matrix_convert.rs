//! N-port network matrix conversion utilities.
//!
//! Supports conversion between the major microwave network matrix representations:
//!
//! | Symbol | Name | Domain |
//! |--------|------|--------|
//! | **S**  | Scattering (S-parameters) | Wave amplitudes |
//! | **Z**  | Impedance (Z-parameters)  | Voltages/currents |
//! | **Y**  | Admittance (Y-parameters) | Voltages/currents |
//! | **T**  | Transfer scattering (T-parameters) | 2-port cascade |
//! | **ABCD** | Chain (ABCD-matrix)     | 2-port transmission |
//!
//! All conversions use generalised N-port formulas with reference impedance Z₀
//! (scalar, same for all ports; for per-port Z₀, scale the S-matrix externally).
//!
//! # References
//! Pozar, D. M. *Microwave Engineering*, 4th ed., Appendix A.
//! Frickey, D. A. (1994) "Conversions between S, Z, Y, h, ABCD, and T parameters."
//! *IEEE Trans. Microw. Theory Tech.* MTT-42(2):205–211.

use nalgebra::DMatrix;
use num_complex::Complex64;

// ── Type alias for clarity ────────────────────────────────────────────────────

/// N×N complex matrix.
pub type ComplexMatrix = DMatrix<Complex64>;

/// 2×2 ABCD chain matrix stored as [[A,B],[C,D]].
pub type Abcd = [[Complex64; 2]; 2];

// ── Helpers ───────────────────────────────────────────────────────────────────

/// Build the N×N identity matrix scaled by `z0`.
fn z0_eye(n: usize, z0: f64) -> ComplexMatrix {
    let mut m = ComplexMatrix::zeros(n, n);
    for i in 0..n {
        m[(i, i)] = Complex64::new(z0, 0.0);
    }
    m
}

/// Invert a complex matrix using nalgebra LU decomposition.
///
/// Returns `None` if the matrix is singular.
pub fn invert(m: &ComplexMatrix) -> Option<ComplexMatrix> {
    m.clone().try_inverse()
}

// ── S ↔ Z ────────────────────────────────────────────────────────────────────

/// Convert S-parameters to Z-parameters (impedance matrix).
///
/// ```text
///     Z = Z₀ · (I + S)(I − S)⁻¹
/// ```
///
/// Returns `None` if (I − S) is singular.
pub fn s_to_z(s: &ComplexMatrix, z0: f64) -> Option<ComplexMatrix> {
    let n = s.nrows();
    let eye = ComplexMatrix::identity(n, n);
    let z0c = Complex64::new(z0, 0.0);

    let i_minus_s = &eye - s;
    let inv = invert(&i_minus_s)?;
    Some((&eye + s) * inv * z0c)
}

/// Convert Z-parameters to S-parameters.
///
/// ```text
///     S = (Z − Z₀·I)(Z + Z₀·I)⁻¹
/// ```
///
/// Returns `None` if (Z + Z₀·I) is singular.
pub fn z_to_s(z: &ComplexMatrix, z0: f64) -> Option<ComplexMatrix> {
    let n = z.nrows();
    let z0i = z0_eye(n, z0);

    let num = z - &z0i;
    let den = z + &z0i;
    let inv = invert(&den)?;
    Some(num * inv)
}

// ── S ↔ Y ────────────────────────────────────────────────────────────────────

/// Convert S-parameters to Y-parameters (admittance matrix).
///
/// ```text
///     Y = (1/Z₀) · (I − S)(I + S)⁻¹
/// ```
///
/// Returns `None` if (I + S) is singular.
pub fn s_to_y(s: &ComplexMatrix, z0: f64) -> Option<ComplexMatrix> {
    let n = s.nrows();
    let eye = ComplexMatrix::identity(n, n);
    let y0c = Complex64::new(1.0 / z0, 0.0);

    let i_plus_s = &eye + s;
    let inv = invert(&i_plus_s)?;
    Some((&eye - s) * inv * y0c)
}

/// Convert Y-parameters to S-parameters.
///
/// ```text
///     S = (I − Z₀·Y)(I + Z₀·Y)⁻¹
/// ```
///
/// Returns `None` if (I + Z₀·Y) is singular.
pub fn y_to_s(y: &ComplexMatrix, z0: f64) -> Option<ComplexMatrix> {
    let n = y.nrows();
    let eye = ComplexMatrix::identity(n, n);
    let z0c = Complex64::new(z0, 0.0);

    let z0y = y * z0c;
    let i_minus = &eye - &z0y;
    let i_plus  = &eye + &z0y;
    let inv = invert(&i_plus)?;
    Some(i_minus * inv)
}

// ── Z ↔ Y ────────────────────────────────────────────────────────────────────

/// Convert Z-parameters to Y-parameters: Y = Z⁻¹.
///
/// Returns `None` if Z is singular.
pub fn z_to_y(z: &ComplexMatrix) -> Option<ComplexMatrix> {
    invert(z)
}

/// Convert Y-parameters to Z-parameters: Z = Y⁻¹.
///
/// Returns `None` if Y is singular.
pub fn y_to_z(y: &ComplexMatrix) -> Option<ComplexMatrix> {
    invert(y)
}

// ── 2-port: S ↔ ABCD ─────────────────────────────────────────────────────────

/// Convert 2-port S-parameters to ABCD (chain/transmission) matrix.
///
/// Frickey (1994) formulation with reference impedance Z₀ [Ω]:
///
/// ```text
///     A = (1 + S11 − S22 − ΔS) / (2·S21)
///     B = Z₀·(1 + S11 + S22 + ΔS) / (2·S21)
///     C = (1 − S11 − S22 + ΔS) / (2·Z₀·S21)
///     D = (1 − S11 + S22 − ΔS) / (2·S21)
///
///     ΔS = S11·S22 − S12·S21
/// ```
///
/// Returns `None` if S21 = 0.
pub fn s_to_abcd(
    s11: Complex64, s12: Complex64, s21: Complex64, s22: Complex64,
    z0: f64,
) -> Option<Abcd> {
    if s21.norm() < 1e-300 { return None; }
    let z0c = Complex64::new(z0, 0.0);
    let two_s21 = Complex64::new(2.0, 0.0) * s21;
    let delta = s11 * s22 - s12 * s21;
    let one = Complex64::new(1.0, 0.0);

    let a = (one + s11 - s22 - delta) / two_s21;
    let b = z0c * (one + s11 + s22 + delta) / two_s21;
    let c = (one - s11 - s22 + delta) / (z0c * two_s21);
    let d = (one - s11 + s22 - delta) / two_s21;

    Some([[a, b], [c, d]])
}

/// Convert 2-port ABCD matrix to S-parameters.
///
/// ```text
///     S11 = (A + B/Z₀ − C·Z₀ − D) / denom
///     S12 = 2·(AD − BC)            / denom
///     S21 = 2                      / denom
///     S22 = (−A + B/Z₀ − C·Z₀ + D)/ denom
///
///     denom = A + B/Z₀ + C·Z₀ + D
/// ```
pub fn abcd_to_s(abcd: &Abcd, z0: f64) -> Option<[[Complex64; 2]; 2]> {
    let [[a, b], [c, d]] = *abcd;
    let z0c = Complex64::new(z0, 0.0);
    let y0c = Complex64::new(1.0 / z0, 0.0);
    let two = Complex64::new(2.0, 0.0);

    let denom = a + b * y0c + c * z0c + d;
    if denom.norm() < 1e-300 { return None; }

    let s11 = (a + b * y0c - c * z0c - d) / denom;
    let s12 = two * (a * d - b * c)       / denom;
    let s21 = two                          / denom;
    let s22 = (-a + b * y0c - c * z0c + d)/ denom;

    Some([[s11, s12], [s21, s22]])
}

// ── 2-port ABCD cascade ───────────────────────────────────────────────────────

/// Cascade two 2-port networks by multiplying their ABCD matrices.
///
/// If network A has ABCD_A and network B has ABCD_B, the cascade is:
/// ```text
///     [ABCD_total] = [ABCD_A] × [ABCD_B]
/// ```
pub fn cascade_abcd(a: &Abcd, b: &Abcd) -> Abcd {
    let [[a11, a12], [a21, a22]] = *a;
    let [[b11, b12], [b21, b22]] = *b;

    [
        [a11 * b11 + a12 * b21,   a11 * b12 + a12 * b22],
        [a21 * b11 + a22 * b21,   a21 * b12 + a22 * b22],
    ]
}

// ── 2-port: S ↔ T (Transfer scattering) ──────────────────────────────────────

/// Convert 2-port S-parameters to T-parameters (transfer scattering matrix).
///
/// Convention (Pozar App. A):
/// ```text
///     T11 =  −(S11·S22 − S12·S21) / S21 = −ΔS / S21
///     T12 =  S11 / S21
///     T21 = −S22 / S21
///     T22 =  1   / S21
/// ```
///
/// Returns `None` if S21 = 0.
pub fn s_to_t(
    s11: Complex64, s12: Complex64, s21: Complex64, s22: Complex64,
) -> Option<[[Complex64; 2]; 2]> {
    if s21.norm() < 1e-300 { return None; }
    let delta = s11 * s22 - s12 * s21;
    let one = Complex64::new(1.0, 0.0);

    Some([
        [-delta / s21,   s11 / s21],
        [-s22 / s21,     one / s21],
    ])
}

/// Convert T-parameters back to 2-port S-parameters.
///
/// ```text
///     S21 =  1   / T22
///     S11 =  T12 / T22
///     S22 = −T21 / T22
///     S12 =  T11 − T12·T21 / T22
/// ```
///
/// Returns `None` if T22 = 0.
pub fn t_to_s(t: &[[Complex64; 2]; 2]) -> Option<[[Complex64; 2]; 2]> {
    let [[t11, t12], [t21, t22]] = *t;
    if t22.norm() < 1e-300 { return None; }

    let s21 = Complex64::new(1.0, 0.0) / t22;
    let s11 = t12 / t22;
    let s22 = -t21 / t22;
    let s12 = t11 - t12 * t21 / t22;

    Some([[s11, s12], [s21, s22]])
}

// ── Utility: convert &[Vec<Complex64>] ↔ DMatrix ─────────────────────────────

/// Convert a Vec-of-Vec S-matrix to a `ComplexMatrix`.
pub fn s_matrix_to_dmatrix(s: &[Vec<Complex64>]) -> ComplexMatrix {
    let n = s.len();
    let mut m = ComplexMatrix::zeros(n, n);
    for i in 0..n {
        for j in 0..s[i].len().min(n) {
            m[(i, j)] = s[i][j];
        }
    }
    m
}

/// Convert a `ComplexMatrix` to a Vec-of-Vec.
pub fn dmatrix_to_s_matrix(m: &ComplexMatrix) -> Vec<Vec<Complex64>> {
    let n = m.nrows();
    (0..n).map(|i| (0..n).map(|j| m[(i, j)]).collect()).collect()
}

// ── CSV writer ────────────────────────────────────────────────────────────────

/// One frequency point of a converted network matrix.
#[derive(Debug, Clone)]
pub struct MatrixPoint {
    pub freq_hz: f64,
    /// Rows × columns of the requested matrix (Re, Im interleaved).
    pub data: Vec<Vec<Complex64>>,
    pub kind: MatrixKind,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MatrixKind { S, Z, Y }

/// Write a sweep of Z or Y matrix data as CSV.
///
/// Columns: `FreqHz, M11_Re, M11_Im, M12_Re, M12_Im, ...`
pub fn write_matrix_csv(
    points: &[MatrixPoint],
    output_dir: &std::path::Path,
) -> Result<std::path::PathBuf, std::io::Error> {
    use std::io::Write;
    if points.is_empty() { return Ok(output_dir.to_path_buf()); }

    let dir = output_dir.join("postpro");
    std::fs::create_dir_all(&dir)?;
    let label = match points[0].kind {
        MatrixKind::S => "s_matrix",
        MatrixKind::Z => "z_matrix",
        MatrixKind::Y => "y_matrix",
    };
    let path = dir.join(format!("{label}.csv"));
    let mut f = std::fs::File::create(&path)?;

    // Header
    let n = points[0].data.len();
    let mut header = "FreqHz".to_string();
    for i in 0..n {
        for j in 0..n {
            header.push_str(&format!(",M{}{}_Re,M{}{}_Im", i+1, j+1, i+1, j+1));
        }
    }
    writeln!(f, "{header}")?;

    for p in points {
        let mut row = format!("{:.9e}", p.freq_hz);
        for i in 0..n {
            for j in 0..n {
                let c = p.data[i][j];
                row.push_str(&format!(",{:.9e},{:.9e}", c.re, c.im));
            }
        }
        writeln!(f, "{row}")?;
    }

    Ok(path)
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use num_complex::Complex64;
    use std::f64::consts::FRAC_PI_4;

    fn c(re: f64, im: f64) -> Complex64 { Complex64::new(re, im) }

    /// S→Z→S round-trip should recover the original S-matrix.
    #[test]
    fn s_z_s_roundtrip() {
        let s = DMatrix::from_row_slice(2, 2, &[
            c(0.1, 0.05),  c(0.7, -0.1),
            c(0.7, -0.1),  c(0.15, 0.03),
        ]);
        let z0 = 50.0;
        let z = s_to_z(&s, z0).expect("s_to_z");
        let s_back = z_to_s(&z, z0).expect("z_to_s");
        for i in 0..2 {
            for j in 0..2 {
                let err = (s[(i,j)] - s_back[(i,j)]).norm();
                assert!(err < 1e-12, "S[{i},{j}] roundtrip err={err:.3e}");
            }
        }
    }

    /// S→Y→S round-trip.
    #[test]
    fn s_y_s_roundtrip() {
        let s = DMatrix::from_row_slice(2, 2, &[
            c(0.0, 0.0),  c(0.9, 0.0),
            c(0.9, 0.0),  c(0.0, 0.0),
        ]);
        let z0 = 50.0;
        let y = s_to_y(&s, z0).expect("s_to_y");
        let s_back = y_to_s(&y, z0).expect("y_to_s");
        for i in 0..2 {
            for j in 0..2 {
                let err = (s[(i,j)] - s_back[(i,j)]).norm();
                assert!(err < 1e-12, "S[{i},{j}] roundtrip err={err:.3e}");
            }
        }
    }

    /// Z→Y→Z round-trip: Y = Z⁻¹.
    #[test]
    fn z_y_roundtrip() {
        // Symmetric 2-port Z-matrix
        let z = DMatrix::from_row_slice(2, 2, &[
            c(100.0, 0.0),  c(20.0, 10.0),
            c(20.0, 10.0),  c(80.0, 0.0),
        ]);
        let y = z_to_y(&z).expect("z_to_y");
        let z_back = y_to_z(&y).expect("y_to_z");
        for i in 0..2 {
            for j in 0..2 {
                let err = (z[(i,j)] - z_back[(i,j)]).norm();
                assert!(err < 1e-10, "Z[{i},{j}] err={err:.3e}");
            }
        }
    }

    /// S→ABCD→S round-trip for a simple 2-port.
    #[test]
    fn s_abcd_s_roundtrip() {
        let s11 = c(0.1, 0.05);
        let s12 = c(0.7, -0.1);
        let s21 = c(0.8, 0.1);
        let s22 = c(0.12, 0.03);
        let z0 = 50.0;

        let abcd = s_to_abcd(s11, s12, s21, s22, z0).expect("s_to_abcd");
        let s_back = abcd_to_s(&abcd, z0).expect("abcd_to_s");

        let errs = [
            (s11 - s_back[0][0]).norm(),
            (s12 - s_back[0][1]).norm(),
            (s21 - s_back[1][0]).norm(),
            (s22 - s_back[1][1]).norm(),
        ];
        for (idx, &e) in errs.iter().enumerate() {
            assert!(e < 1e-12, "S param {idx} roundtrip err={e:.3e}");
        }
    }

    /// S→T→S round-trip.
    #[test]
    fn s_t_s_roundtrip() {
        let s11 = c(0.2, 0.1);
        let s12 = c(0.6, -0.05);
        let s21 = c(0.65, 0.05);
        let s22 = c(0.1, 0.02);

        let t = s_to_t(s11, s12, s21, s22).expect("s_to_t");
        let s_back = t_to_s(&t).expect("t_to_s");

        assert!((s11 - s_back[0][0]).norm() < 1e-12);
        assert!((s12 - s_back[0][1]).norm() < 1e-12);
        assert!((s21 - s_back[1][0]).norm() < 1e-12);
        assert!((s22 - s_back[1][1]).norm() < 1e-12);
    }

    /// Cascade of two identical 2-ports gives a double-length network.
    /// For a lossless through: S = [[0,1],[1,0]], ABCD = [[0,0],[0,1]]... 
    /// Let's use a simple delay: S21 = exp(−jπ/4).
    #[test]
    fn cascade_abcd_doubles_delay() {
        // Single section: matched lossless 2-port with S21 = exp(−jπ/4)
        let phi = FRAC_PI_4;
        let s21 = c(phi.cos(), -phi.sin());
        let s12 = s21;
        let s11 = c(0.0, 0.0);
        let s22 = c(0.0, 0.0);

        let abcd = s_to_abcd(s11, s12, s21, s22, 50.0).expect("s_to_abcd");
        let cascade = cascade_abcd(&abcd, &abcd);
        let s_cascaded = abcd_to_s(&cascade, 50.0).expect("abcd_to_s");

        // Cascading two sections doubles the phase: S21 should be exp(−j·2·π/4) = exp(−jπ/2) = −j
        let expected_s21 = c(0.0, -1.0);
        let err = (s_cascaded[1][0] - expected_s21).norm();
        assert!(err < 1e-10, "Cascaded S21={:.4}+j{:.4}, expected 0-j1, err={err:.3e}",
            s_cascaded[1][0].re, s_cascaded[1][0].im);
    }

    /// dmatrix_to_s_matrix and s_matrix_to_dmatrix are inverses.
    #[test]
    fn matrix_conversion_roundtrip() {
        let orig: Vec<Vec<Complex64>> = vec![
            vec![c(1.0, 2.0), c(3.0, 4.0)],
            vec![c(5.0, 6.0), c(7.0, 8.0)],
        ];
        let dm = s_matrix_to_dmatrix(&orig);
        let back = dmatrix_to_s_matrix(&dm);
        for i in 0..2 {
            for j in 0..2 {
                assert_eq!(orig[i][j], back[i][j]);
            }
        }
    }
}
