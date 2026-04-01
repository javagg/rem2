use thiserror::Error;

#[derive(Debug, Error)]
pub enum RemError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    #[error("Config parse error: {0}")]
    Config(String),

    #[error("Mesh error: {0}")]
    Mesh(String),

    #[error("Solver did not converge after {max_iter} iterations (residual={residual:.3e})")]
    SolverNotConverged { max_iter: usize, residual: f64 },

    #[error("Unknown problem type: '{0}'")]
    UnknownProblemType(String),

    #[error("Unknown or unsupported file format")]
    UnknownFormat,

    #[error("Feature not yet implemented: {0}")]
    NotImplemented(String),

    #[error("{0}")]
    Other(String),
}

pub type RemResult<T> = Result<T, RemError>;

