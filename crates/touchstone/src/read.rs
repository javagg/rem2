//! Touchstone (.s{N}p) file reader — RI, MA, DB formats supported.
use num_complex::Complex64;
use std::f64::consts::PI;

#[derive(Debug, Clone)]
pub struct TouchstoneFile {
    pub n_ports: usize,
    pub freqs_hz: Vec<f64>,
    /// S-matrix row-major per frequency: s_data[f][i*n+j] = S_{i+1,j+1}
    pub s_data: Vec<Vec<Complex64>>,
    pub z0: f64,
}

#[derive(Debug, thiserror::Error)]
pub enum TsReadError {
    #[error("parse error: {0}")]
    Parse(String),
    #[error("io error: {0}")]
    Io(#[from] std::io::Error),
}

// ─── Option-line format enum ─────────────────────────────────────────────────
#[derive(Debug, Clone, Copy, PartialEq)]
enum Fmt { Ri, Ma, Db }

// ─── Public API ──────────────────────────────────────────────────────────────

/// Parse a Touchstone file from a string.
///
/// Supports RI, MA, and DB data formats. Handles single-line (N≤2) and
/// multi-line (N≥3) frequency blocks per Touchstone 1.0 specification.
///
/// # Example
/// ```
/// let ts = rem_touchstone::read_snp("# GHz S RI R 50\n1.0e9  0.5 -0.3\n").unwrap();
/// assert_eq!(ts.n_ports, 1);
/// ```
pub fn read_snp(content: &str) -> Result<TouchstoneFile, TsReadError> {
    let mut freq_scale = 1.0e9_f64; // default GHz
    let mut fmt = Fmt::Ri;
    let mut z0 = 50.0_f64;
    let mut n_ports_from_option: Option<usize> = None;

    // Collect all numeric tokens from data lines, tracking freq-block structure
    let mut data_tokens: Vec<f64> = Vec::new();
    let mut option_parsed = false;

    for line in content.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('!') {
            continue;
        }
        if trimmed.starts_with('#') {
            if !option_parsed {
                parse_option_line(trimmed, &mut freq_scale, &mut fmt, &mut z0, &mut n_ports_from_option)?;
                option_parsed = true;
            }
            continue;
        }
        // Data line: tokenize all numeric values
        for tok in trimmed.split_whitespace() {
            let v: f64 = tok.parse().map_err(|_| {
                TsReadError::Parse(format!("non-numeric token in data: '{tok}'"))
            })?;
            data_tokens.push(v);
        }
    }

    if data_tokens.is_empty() {
        return Err(TsReadError::Parse("no data found in file".to_string()));
    }

    // Infer n_ports from the first data line (the first cluster of tokens
    // beginning with the frequency). We scan until we hit what looks like the
    // second frequency marker.
    let n_ports = infer_n_ports(&data_tokens)?;

    // n_ports² S-pairs per frequency = 2*n_ports² scalar values + 1 freq token
    let per_freq = 1 + 2 * n_ports * n_ports;
    if data_tokens.len() % per_freq != 0 {
        return Err(TsReadError::Parse(format!(
            "token count {} is not a multiple of {} (n_ports={n_ports})",
            data_tokens.len(), per_freq
        )));
    }

    let n_freqs = data_tokens.len() / per_freq;
    let mut freqs_hz = Vec::with_capacity(n_freqs);
    let mut s_data: Vec<Vec<Complex64>> = Vec::with_capacity(n_freqs);

    for fi in 0..n_freqs {
        let base = fi * per_freq;
        freqs_hz.push(data_tokens[base] * freq_scale);
        let mut row: Vec<Complex64> = Vec::with_capacity(n_ports * n_ports);
        for pair in 0..n_ports * n_ports {
            let a = data_tokens[base + 1 + 2 * pair];
            let b = data_tokens[base + 1 + 2 * pair + 1];
            row.push(to_ri(a, b, fmt));
        }
        s_data.push(row);
    }

    Ok(TouchstoneFile { n_ports, freqs_hz, s_data, z0 })
}

// ─── Helpers ─────────────────────────────────────────────────────────────────

fn parse_option_line(
    line: &str,
    freq_scale: &mut f64,
    fmt: &mut Fmt,
    z0: &mut f64,
    _n_ports: &mut Option<usize>,
) -> Result<(), TsReadError> {
    // Format: # <unit> S <fmt> R <z0>
    let upper = line.to_ascii_uppercase();
    let parts: Vec<&str> = upper.split_whitespace().collect();

    for (i, &tok) in parts.iter().enumerate() {
        match tok {
            "HZ"  => *freq_scale = 1.0,
            "KHZ" => *freq_scale = 1e3,
            "MHZ" => *freq_scale = 1e6,
            "GHZ" => *freq_scale = 1e9,
            "RI"  => *fmt = Fmt::Ri,
            "MA"  => *fmt = Fmt::Ma,
            "DB"  => *fmt = Fmt::Db,
            "R"   => {
                if let Some(&z_str) = parts.get(i + 1) {
                    *z0 = z_str.parse().unwrap_or(50.0);
                }
            }
            _ => {}
        }
    }
    Ok(())
}

/// Infer N (number of ports) from the flat token stream using frequency monotonicity.
///
/// For a valid Touchstone file, the frequency column (every `per = 1+2N²`-th token)
/// must be strictly positive and monotonically non-decreasing.  We try N=1,2,…
/// and return the first N for which the token count divides evenly AND the
/// extracted frequencies are monotonically non-decreasing.
fn infer_n_ports(tokens: &[f64]) -> Result<usize, TsReadError> {
    for n in 1..=32_usize {
        let per = 1 + 2 * n * n;
        if tokens.len() % per != 0 {
            continue;
        }
        let n_freqs = tokens.len() / per;
        // Check frequency monotonicity
        let freqs: Vec<f64> = (0..n_freqs).map(|i| tokens[i * per]).collect();
        let mono = freqs.windows(2).all(|w| w[1] >= w[0]);
        let positive = freqs.iter().all(|&f| f >= 0.0);
        if mono && positive {
            return Ok(n);
        }
    }
    Err(TsReadError::Parse(format!(
        "cannot infer n_ports from {} tokens (tried N=1..32 with monotonicity check)",
        tokens.len()
    )))
}

#[inline]
fn to_ri(a: f64, b: f64, fmt: Fmt) -> Complex64 {
    match fmt {
        Fmt::Ri => Complex64::new(a, b),
        Fmt::Ma => Complex64::from_polar(a, b * PI / 180.0),
        Fmt::Db => {
            let mag = 10.0_f64.powf(a / 20.0);
            Complex64::from_polar(mag, b * PI / 180.0)
        }
    }
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::write::{write_snp, TsFormat, TsFreqUnit};

    #[test]
    fn round_trip_s1p_ri() {
        let freqs = vec![1e9, 2e9, 3e9];
        let s_data = vec![
            vec![Complex64::new( 0.5, -0.3)],
            vec![Complex64::new( 0.4, -0.2)],
            vec![Complex64::new( 0.3, -0.1)],
        ];
        let text = write_snp(&freqs, &s_data, 1, 50.0, TsFormat::Ri, TsFreqUnit::Ghz);
        let ts = read_snp(&text).unwrap();
        assert_eq!(ts.n_ports, 1);
        assert_eq!(ts.freqs_hz.len(), 3);
        assert!((ts.freqs_hz[0] - 1e9).abs() < 1e3, "freq[0]={}", ts.freqs_hz[0]);
        assert!((ts.s_data[0][0].re - 0.5).abs() < 1e-6, "Re(S11)={}", ts.s_data[0][0].re);
        assert!((ts.s_data[0][0].im - (-0.3)).abs() < 1e-6, "Im(S11)={}", ts.s_data[0][0].im);
    }

    #[test]
    fn round_trip_s2p_ri() {
        let freqs = vec![1e9, 5e9];
        let row: Vec<Complex64> = vec![
            Complex64::new(0.1,  0.0),
            Complex64::new(0.9,  0.01),
            Complex64::new(0.9, -0.01),
            Complex64::new(0.1,  0.0),
        ];
        let s_data = vec![row.clone(), row];
        let text = write_snp(&freqs, &s_data, 2, 50.0, TsFormat::Ri, TsFreqUnit::Ghz);
        let ts = read_snp(&text).unwrap();
        assert_eq!(ts.n_ports, 2);
        assert_eq!(ts.freqs_hz.len(), 2);
        assert!((ts.s_data[0][0].re - 0.1).abs() < 1e-6);
        assert!((ts.s_data[0][1].re - 0.9).abs() < 1e-6);
    }

    #[test]
    fn round_trip_s3p_multiline() {
        let freqs = vec![1e9, 2e9];
        let row: Vec<Complex64> = (0..9)
            .map(|i| Complex64::new(i as f64 * 0.1, 0.0))
            .collect();
        let s_data = vec![row.clone(), row];
        let text = write_snp(&freqs, &s_data, 3, 50.0, TsFormat::Ri, TsFreqUnit::Ghz);
        let ts = read_snp(&text).unwrap();
        assert_eq!(ts.n_ports, 3);
        assert_eq!(ts.freqs_hz.len(), 2);
        assert!((ts.s_data[1][8].re - 0.8).abs() < 1e-6,
            "S33[f=2] wrong: {}", ts.s_data[1][8].re);
    }

    #[test]
    fn parse_ma_format() {
        // MA: magnitude (linear) + angle in degrees
        let content = "# GHz S MA R 50\n1.0  1.0  0.0\n";
        let ts = read_snp(content).unwrap();
        assert_eq!(ts.n_ports, 1);
        assert!((ts.s_data[0][0].re - 1.0).abs() < 1e-6);
        assert!((ts.s_data[0][0].im).abs() < 1e-6);
    }

    #[test]
    fn parse_db_format() {
        // DB: 0 dB = magnitude 1.0
        let content = "# GHz S DB R 50\n1.0  0.0  0.0\n";
        let ts = read_snp(content).unwrap();
        assert!((ts.s_data[0][0].re - 1.0).abs() < 1e-6,
            "DB parse wrong: re={}", ts.s_data[0][0].re);
    }

    #[test]
    fn skips_comment_lines() {
        let content = "! This is a comment\n# GHz S RI R 50\n! another comment\n1.0  0.5  -0.3\n";
        let ts = read_snp(content).unwrap();
        assert_eq!(ts.n_ports, 1);
        assert!((ts.s_data[0][0].re - 0.5).abs() < 1e-6);
    }
}
