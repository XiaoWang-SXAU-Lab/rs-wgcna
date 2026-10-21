use thiserror::Error;

#[derive(Debug, Error)]
pub enum WgcnaError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("parse error: {0}")]
    Csv(#[from] csv::Error),

    #[error("invalid input: {0}")]
    InvalidInput(String),

    #[error("numeric error: {0}")]
    Numeric(String),
}

pub type Result<T> = std::result::Result<T, WgcnaError>;
