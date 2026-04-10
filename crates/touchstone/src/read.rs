//! Touchstone reader — minimal stub for now.
use num_complex::Complex64;

#[derive(Debug)]
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

/// Parse a Touchstone file from a string.
/// Only RI format supported currently.
pub fn read_snp(content: &str) -> Result<TouchstoneFile, TsReadError> {
    let _ = content;
    Err(TsReadError::Parse("read_snp not yet implemented".to_string()))
}
