//! Near-field data point and CSV I/O.
//!
//! A `NearFieldPoint` carries spatial coordinates and full E/H field vectors
//! (complex-valued) sampled at that location.  Different solvers populate
//! only the fields they can compute — unavailable fields are zero.
//!
//! CSV format header:
//! ```text
//! x,y,z,Ex_re,Ex_im,Ey_re,Ey_im,Ez_re,Ez_im,Hx_re,Hx_im,Hy_re,Hy_im,Hz_re,Hz_im
//! ```

use num_complex::Complex64;
use std::io::{BufRead, Write};
use std::path::Path;

use crate::{RemError, RemResult};

// ---------------------------------------------------------------------------
// Data point
// ---------------------------------------------------------------------------

/// One sampled point with E and H field vectors.
#[derive(Debug, Clone, Default)]
pub struct NearFieldPoint {
    /// Position [m]
    pub x: f64,
    pub y: f64,
    pub z: f64,
    /// Electric field E [V/m]
    pub ex: Complex64,
    pub ey: Complex64,
    pub ez: Complex64,
    /// Magnetic field H [A/m]
    pub hx: Complex64,
    pub hy: Complex64,
    pub hz: Complex64,
}

impl NearFieldPoint {
    /// Construct a point with real-only E field and zero H (typical for
    /// scalar-potential solvers like Driven with real phi).
    pub fn from_real_e(x: f64, y: f64, z: f64, e: [f64; 3]) -> Self {
        Self {
            x, y, z,
            ex: Complex64::new(e[0], 0.0),
            ey: Complex64::new(e[1], 0.0),
            ez: Complex64::new(e[2], 0.0),
            hx: Complex64::ZERO,
            hy: Complex64::ZERO,
            hz: Complex64::ZERO,
        }
    }

    /// Construct a point with full complex E and H fields.
    pub fn from_complex(
        x: f64, y: f64, z: f64,
        ex: Complex64, ey: Complex64, ez: Complex64,
        hx: Complex64, hy: Complex64, hz: Complex64,
    ) -> Self {
        Self { x, y, z, ex, ey, ez, hx, hy, hz }
    }
}

// ---------------------------------------------------------------------------
// CSV write
// ---------------------------------------------------------------------------

const HEADER: &str = "x,y,z,Ex_re,Ex_im,Ey_re,Ey_im,Ez_re,Ez_im,Hx_re,Hx_im,Hy_re,Hy_im,Hz_re,Hz_im";

/// Write near-field points to a CSV file.
///
/// The file is overwritten each call.  Creates parent directories as needed.
pub fn write_near_field_csv(path: &Path, points: &[NearFieldPoint]) -> RemResult<()> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent).map_err(RemError::Io)?;
        }
    }
    let mut f = std::fs::File::create(path).map_err(RemError::Io)?;
    writeln!(f, "{HEADER}").map_err(RemError::Io)?;
    for p in points {
        writeln!(
            f,
            "{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e}",
            p.x, p.y, p.z,
            p.ex.re, p.ex.im, p.ey.re, p.ey.im, p.ez.re, p.ez.im,
            p.hx.re, p.hx.im, p.hy.re, p.hy.im, p.hz.re, p.hz.im,
        ).map_err(RemError::Io)?;
    }
    log::info!("[REM] Near-field CSV written: {} points → {}", points.len(), path.display());
    Ok(())
}

// ---------------------------------------------------------------------------
// CSV read
// ---------------------------------------------------------------------------

/// Read near-field points from a CSV file.
///
/// The first non-comment line must be the standard header.  Lines starting
/// with `#` or `!` are treated as comments and skipped.
pub fn read_near_field_csv(path: &Path) -> RemResult<Vec<NearFieldPoint>> {
    let f = std::fs::File::open(path).map_err(RemError::Io)?;
    let reader = std::io::BufReader::new(f);
    let mut lines = reader.lines();
    let mut points = Vec::new();

    // Skip comments and find header
    loop {
        let line = lines.next().transpose().map_err(RemError::Io)?;
        match line {
            Some(s) if s.starts_with('#') || s.starts_with('!') => continue,
            Some(_header) => break, // header found and consumed
            None => {
                return Err(RemError::Config("near-field CSV is empty or has no header".to_string()));
            }
        }
    }

    for line_result in lines {
        let line = line_result.map_err(RemError::Io)?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') || line.starts_with('!') {
            continue;
        }
        let vals: Vec<&str> = line.split(',').collect();
        if vals.len() < 14 {
            return Err(RemError::NotImplemented(format!(
                "near-field CSV line has {} values, expected 14: '{}'",
                vals.len(), line
            )));
        }
        let parse_f64 = |s: &str| s.trim().parse::<f64>()
            .map_err(|e| RemError::Config(e.to_string()));
        let x = parse_f64(vals[0])?;
        let y = parse_f64(vals[1])?;
        let z = parse_f64(vals[2])?;
        let ex_re = parse_f64(vals[3])?;
        let ex_im = parse_f64(vals[4])?;
        let ey_re = parse_f64(vals[5])?;
        let ey_im = parse_f64(vals[6])?;
        let ez_re = parse_f64(vals[7])?;
        let ez_im = parse_f64(vals[8])?;
        let hx_re = parse_f64(vals[9])?;
        let hx_im = parse_f64(vals[10])?;
        let hy_re = parse_f64(vals[11])?;
        let hy_im = parse_f64(vals[12])?;
        let hz_re = parse_f64(vals[13])?;
        let hz_im = parse_f64(vals[14])?;

        points.push(NearFieldPoint {
            x, y, z,
            ex: Complex64::new(ex_re, ex_im),
            ey: Complex64::new(ey_re, ey_im),
            ez: Complex64::new(ez_re, ez_im),
            hx: Complex64::new(hx_re, hx_im),
            hy: Complex64::new(hy_re, hy_im),
            hz: Complex64::new(hz_re, hz_im),
        });
    }

    log::info!("[REM] Near-field CSV loaded: {} points from {}", points.len(), path.display());
    Ok(points)
}

// ---------------------------------------------------------------------------
// Spatial interpolation (IDW — inverse distance weighting)
// ---------------------------------------------------------------------------

/// Interpolate the electric field at `target` from the `k` nearest points
/// using inverse distance weighting.  Falls back to nearest-neighbour if
/// fewer than `k` points exist or if the target coincides with a data point.
pub fn interpolate_e_at(
    target: [f64; 3],
    points: &[NearFieldPoint],
    k: usize,
) -> Complex64 {
    if points.is_empty() {
        return Complex64::ZERO;
    }

    let mut dists: Vec<(f64, &NearFieldPoint)> = points.iter()
        .map(|p| {
            let dx = p.x - target[0];
            let dy = p.y - target[1];
            let dz = p.z - target[2];
            let d: f64 = (dx * dx + dy * dy + dz * dz).sqrt();
            (d, p)
        })
        .collect();

    // Partial sort to get k nearest
    let k = k.min(dists.len());
    dists.select_nth_unstable_by(k - 1, |a: &(f64, &NearFieldPoint), b: &(f64, &NearFieldPoint)| a.0.partial_cmp(&b.0).unwrap());
    dists.truncate(k);

    // Check for exact coincidence (distance ≈ 0)
    if dists[0].0 < 1e-300 {
        return dists[0].1.ex;
    }

    let mut sum_w = 0.0_f64;
    let mut ex_sum = Complex64::ZERO;
    let mut ey_sum = Complex64::ZERO;
    let mut ez_sum = Complex64::ZERO;

    for &(d, p) in &dists {
        let w = 1.0_f64 / d.max(1e-300);
        sum_w += w;
        ex_sum += w * p.ex;
        ey_sum += w * p.ey;
        ez_sum += w * p.ez;
    }

    // Return |E| magnitude as a scalar in the direction of nearest E
    let ex = ex_sum / sum_w;
    let ey = ey_sum / sum_w;
    let ez = ez_sum / sum_w;

    // Return the interpolated E vector norm as the effective field amplitude
    (ex * ex + ey * ey + ez * ez).sqrt()
}

/// Same as `interpolate_e_at` but returns the full interpolated E vector.
pub fn interpolate_e_vec_at(
    target: [f64; 3],
    points: &[NearFieldPoint],
    k: usize,
) -> [Complex64; 3] {
    if points.is_empty() {
        return [Complex64::ZERO; 3];
    }

    let mut dists: Vec<(f64, &NearFieldPoint)> = points.iter()
        .map(|p| {
            let dx = p.x - target[0];
            let dy = p.y - target[1];
            let dz = p.z - target[2];
            let d: f64 = (dx * dx + dy * dy + dz * dz).sqrt();
            (d, p)
        })
        .collect();

    let k = k.min(dists.len());
    dists.select_nth_unstable_by(k - 1, |a: &(f64, &NearFieldPoint), b: &(f64, &NearFieldPoint)| a.0.partial_cmp(&b.0).unwrap());
    dists.truncate(k);

    if dists[0].0 < 1e-300 {
        return [dists[0].1.ex, dists[0].1.ey, dists[0].1.ez];
    }

    let mut sum_w = 0.0_f64;
    let mut ex_sum = Complex64::ZERO;
    let mut ey_sum = Complex64::ZERO;
    let mut ez_sum = Complex64::ZERO;

    for &(d, p) in &dists {
        let w = 1.0_f64 / d.max(1e-300);
        sum_w += w;
        ex_sum += w * p.ex;
        ey_sum += w * p.ey;
        ez_sum += w * p.ez;
    }

    [ex_sum / sum_w, ey_sum / sum_w, ez_sum / sum_w]
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn tmp_path(name: &str) -> std::path::PathBuf {
        std::env::temp_dir().join(name)
    }

    #[test]
    fn write_and_read_roundtrip() {
        let pts = vec![
            NearFieldPoint::from_complex(
                0.0, 0.0, 0.0,
                Complex64::new(1.0, 0.0), Complex64::new(0.5, 0.1), Complex64::new(0.0, 0.0),
                Complex64::new(0.0, 0.0), Complex64::new(0.0, 1.0), Complex64::new(0.0, 0.0),
            ),
            NearFieldPoint::from_complex(
                1.0, 0.0, 0.0,
                Complex64::new(0.8, -0.2), Complex64::new(0.3, 0.05), Complex64::new(0.1, 0.0),
                Complex64::new(0.0, 0.0), Complex64::new(0.0, 0.0), Complex64::new(0.0, 0.0),
            ),
        ];

        let path = tmp_path("nf_roundtrip.csv");
        let _ = std::fs::remove_file(&path);
        write_near_field_csv(&path, &pts).expect("write failed");

        let loaded = read_near_field_csv(&path).expect("read failed");
        assert_eq!(loaded.len(), 2);

        // Check first point
        assert!((loaded[0].x - 0.0).abs() < 1e-12);
        assert!((loaded[0].ex.re - 1.0).abs() < 1e-12);
        assert!((loaded[0].ex.im - 0.0).abs() < 1e-12);
        assert!((loaded[0].hy.im - 1.0).abs() < 1e-12);

        // Check second point
        assert!((loaded[1].ex.re - 0.8).abs() < 1e-6);
        assert!((loaded[1].ex.im + 0.2).abs() < 1e-6);

        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn read_empty_file_returns_ok() {
        let path = tmp_path("nf_empty.csv");
        let _ = std::fs::remove_file(&path);
        let mut f = std::fs::File::create(&path).unwrap();
        writeln!(f, "{HEADER}").unwrap();
        let pts = read_near_field_csv(&path).expect("should succeed");
        assert!(pts.is_empty());
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn read_missing_file_returns_error() {
        let path = tmp_path("nf_nonexistent.csv");
        let _ = std::fs::remove_file(&path);
        assert!(read_near_field_csv(&path).is_err());
    }

    #[test]
    fn read_malformed_line_returns_error() {
        let path = tmp_path("nf_malformed.csv");
        let _ = std::fs::remove_file(&path);
        let mut f = std::fs::File::create(&path).unwrap();
        writeln!(f, "{HEADER}").unwrap();
        writeln!(f, "0.0,0.0,0.0,bad").unwrap(); // only 4 fields
        assert!(read_near_field_csv(&path).is_err());
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn comment_lines_are_skipped() {
        let path = tmp_path("nf_comments.csv");
        let _ = std::fs::remove_file(&path);
        let mut f = std::fs::File::create(&path).unwrap();
        writeln!(f, "# Generated by REM MoM solver").unwrap();
        writeln!(f, "! Frequency: 1 GHz").unwrap();
        writeln!(f, "{HEADER}").unwrap();
        writeln!(f, "0.0,0.0,0.0,1,0,0,0,0,0,0,0,0,0,0,0").unwrap();
        let pts = read_near_field_csv(&path).expect("should succeed");
        assert_eq!(pts.len(), 1);
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn from_real_e_zero_h() {
        let p = NearFieldPoint::from_real_e(1.0, 2.0, 3.0, [10.0, 0.0, 5.0]);
        assert_eq!(p.ex.re, 10.0);
        assert_eq!(p.ex.im, 0.0);
        assert_eq!(p.hx, Complex64::ZERO);
        assert_eq!(p.hz, Complex64::ZERO);
    }

    #[test]
    fn interpolate_exact_match() {
        let pts = vec![
            NearFieldPoint::from_complex(
                0.0, 0.0, 0.0,
                Complex64::new(1.0, 0.0), Complex64::ZERO, Complex64::ZERO,
                Complex64::ZERO, Complex64::ZERO, Complex64::ZERO,
            ),
        ];
        let e = interpolate_e_at([0.0, 0.0, 0.0], &pts, 3);
        assert!((e.re - 1.0).abs() < 1e-12);
    }

    #[test]
    fn interpolate_empty_points_returns_zero() {
        let pts: Vec<NearFieldPoint> = vec![];
        let e = interpolate_e_at([0.0, 0.0, 0.0], &pts, 3);
        assert_eq!(e, Complex64::ZERO);
    }

    #[test]
    fn interpolate_idw_symmetric() {
        // Two equidistant points with same E → result should be same E
        let pts = vec![
            NearFieldPoint::from_complex(
                -1.0, 0.0, 0.0,
                Complex64::new(1.0, 0.0), Complex64::ZERO, Complex64::ZERO,
                Complex64::ZERO, Complex64::ZERO, Complex64::ZERO,
            ),
            NearFieldPoint::from_complex(
                1.0, 0.0, 0.0,
                Complex64::new(1.0, 0.0), Complex64::ZERO, Complex64::ZERO,
                Complex64::ZERO, Complex64::ZERO, Complex64::ZERO,
            ),
        ];
        let e = interpolate_e_at([0.0, 0.0, 0.0], &pts, 3);
        assert!((e.re - 1.0).abs() < 1e-12, "IDW of equal E should give same E, got {}", e.re);
    }

    #[test]
    fn interpolate_e_vec_returns_components() {
        let pts = vec![
            NearFieldPoint::from_complex(
                0.0, 0.0, 0.0,
                Complex64::new(3.0, 0.0), Complex64::new(4.0, 0.0), Complex64::ZERO,
                Complex64::ZERO, Complex64::ZERO, Complex64::ZERO,
            ),
        ];
        let [ex, ey, ez] = interpolate_e_vec_at([0.0, 0.0, 0.0], &pts, 3);
        assert!((ex.re - 3.0).abs() < 1e-12);
        assert!((ey.re - 4.0).abs() < 1e-12);
        assert!(ez.re.abs() < 1e-12);
    }

    #[test]
    fn write_creates_parent_dirs() {
        let path = tmp_path("nf_subdir/nested/test.csv");
        let _ = std::fs::remove_dir_all(path.parent().unwrap());
        let pts = vec![NearFieldPoint::default()];
        write_near_field_csv(&path, &pts).expect("should create dirs and write");
        let loaded = read_near_field_csv(&path).expect("should read back");
        assert_eq!(loaded.len(), 1);
        let _ = std::fs::remove_dir_all(path.parent().unwrap().parent().unwrap());
    }

    #[test]
    fn interpolate_single_point_uses_that_point() {
        let pts = vec![
            NearFieldPoint::from_complex(
                0.0, 0.0, 0.0,
                Complex64::new(5.0, 2.0), Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0),
                Complex64::ZERO, Complex64::ZERO, Complex64::ZERO,
            ),
            NearFieldPoint::from_complex(
                10.0, 10.0, 10.0,
                Complex64::new(0.0, 0.0), Complex64::ZERO, Complex64::ZERO,
                Complex64::ZERO, Complex64::ZERO, Complex64::ZERO,
            ),
        ];
        // Very close to first point → should be dominated by it
        let e = interpolate_e_at([1e-15, 1e-15, 1e-15], &pts, 3);
        assert!((e.re - 5.0).abs() < 0.1, "should be dominated by nearest point");
    }

    #[test]
    fn interpolate_k_larger_than_points() {
        let pts = vec![
            NearFieldPoint::from_complex(
                0.0, 0.0, 0.0,
                Complex64::new(1.0, 0.0), Complex64::ZERO, Complex64::ZERO,
                Complex64::ZERO, Complex64::ZERO, Complex64::ZERO,
            ),
        ];
        // Ask for k=10 but only 1 point available → should still work
        let e = interpolate_e_at([1.0, 0.0, 0.0], &pts, 10);
        assert!(e.re.is_finite());
    }
}
