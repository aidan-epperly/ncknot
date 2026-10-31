//! The error type.

use std::fmt;

/// Why an operation failed.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum KnotError {
    /// The input does not describe a diagram, braid or presentation.
    Invalid(String),
    /// The input is valid but the operation does not handle it (for example
    /// the signature of a diagram whose projection is disconnected).
    Unsupported(String),
}

impl fmt::Display for KnotError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            KnotError::Invalid(s) => write!(f, "invalid input: {s}"),
            KnotError::Unsupported(s) => write!(f, "unsupported: {s}"),
        }
    }
}

impl std::error::Error for KnotError {}

/// `Result<T, KnotError>`.
pub type Result<T> = std::result::Result<T, KnotError>;
