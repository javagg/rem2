pub mod constants;
pub mod error;
pub mod memory;
pub mod sparse;

pub use constants::{EPS0, MU0, C0, ETA0, NU0};
pub use error::{RemError, RemResult};
pub use memory::report_peak_memory;
pub use sparse::{CsrMatrix, TripletMatrix, SolveResult, solve_pcg};
