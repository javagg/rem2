pub mod constants;
pub mod error;
pub mod sparse;

pub use error::{RemError, RemResult};
pub use sparse::{CsrMatrix, TripletMatrix, SolveResult, solve_pcg};
