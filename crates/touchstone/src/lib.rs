//! Touchstone (.s{N}p) file read and write utilities.
//!
//! Supports N-port Touchstone 1.0 format, RI data only (MA/DB planned).
//! Used by `rem-mom` (S-parameter output) and `rem-driven` (VF circuit synthesis).

pub mod write;
pub mod read;

pub use write::{write_snp, TsFormat, TsFreqUnit};
pub use read::{read_snp, TouchstoneFile, TsReadError};
