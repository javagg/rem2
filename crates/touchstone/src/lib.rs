//! Touchstone (.s{N}p) file read and write utilities.
//!
//! Supports N-port Touchstone 1.0 format, RI data only (MA/DB planned).
//! Used by `rem-mom` (S-parameter output) and `rem-driven` (VF circuit synthesis).

pub mod write;
pub mod read;
pub mod matrix_convert;

pub use write::{write_snp, TsFormat, TsFreqUnit};
pub use read::{read_snp, TouchstoneFile, TsReadError};
pub use matrix_convert::{
    s_to_z, z_to_s, s_to_y, y_to_s, z_to_y, y_to_z,
    s_to_abcd, abcd_to_s, s_to_t, t_to_s, cascade_abcd,
    s_matrix_to_dmatrix, dmatrix_to_s_matrix,
    MatrixPoint, MatrixKind, write_matrix_csv,
};
