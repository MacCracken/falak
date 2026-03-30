//! Error types for falak.

use std::borrow::Cow;

/// Errors that can occur during orbital computations.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum FalakError {
    /// An input parameter is invalid or out of range.
    #[error("invalid parameter: {0}")]
    InvalidParameter(Cow<'static, str>),

    /// A mathematical operation failed (e.g., division by zero, domain error).
    #[error("math error: {0}")]
    MathError(Cow<'static, str>),

    /// An iterative algorithm failed to converge within the allowed iterations.
    #[error("convergence error: {message} (after {iterations} iterations)")]
    ConvergenceError {
        /// Description of what failed to converge.
        message: Cow<'static, str>,
        /// Number of iterations attempted.
        iterations: u32,
    },

    /// An ephemeris lookup or computation failed.
    #[error("ephemeris error: {0}")]
    EphemerisError(Cow<'static, str>),

    /// An I/O error occurred.
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
}

/// Convenience type alias for falak results.
pub type Result<T> = std::result::Result<T, FalakError>;
