//! Error types and input validation for falak.
//!
//! Provides [`FalakError`] for all fallible operations and validation helpers
//! (`require_finite`, `ensure_finite`) for guarding against NaN/Infinity
//! at API boundaries.

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

    /// Non-finite input detected (NaN or ±Infinity).
    #[error("non-finite input in {context}: {value}")]
    NonFinite {
        /// Which function or parameter had the bad value.
        context: &'static str,
        /// The non-finite value that was detected.
        value: f64,
    },

    /// An I/O error occurred.
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
}

/// Convenience type alias for falak results.
pub type Result<T> = std::result::Result<T, FalakError>;

// ── Input validation helpers ─────────────────────────────────────────────

/// Validate that a single `f64` input is finite (not NaN or ±Infinity).
///
/// Use at API boundaries to catch bad inputs early.
///
/// # Errors
///
/// Returns [`FalakError::NonFinite`] if the value is NaN or infinite.
#[inline]
pub fn require_finite(value: f64, context: &'static str) -> Result<()> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(FalakError::NonFinite { context, value })
    }
}

/// Validate that all `f64` values in a slice are finite.
///
/// # Errors
///
/// Returns [`FalakError::NonFinite`] on the first non-finite value found.
#[inline]
pub fn require_all_finite(values: &[f64], context: &'static str) -> Result<()> {
    for &v in values {
        require_finite(v, context)?;
    }
    Ok(())
}

/// Validate that a computed `f64` result is finite (catches overflow, 0/0, etc.).
///
/// Use after computation to guard outputs before returning to callers.
///
/// # Errors
///
/// Returns [`FalakError::MathError`] if the result is NaN or infinite.
#[inline]
pub fn ensure_finite(value: f64, context: &'static str) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(FalakError::MathError(
            format!("{context}: result is {value}").into(),
        ))
    }
}

/// Validate that all components of a 3D vector are finite.
///
/// # Errors
///
/// Returns [`FalakError::NonFinite`] on the first non-finite component.
#[inline]
pub fn require_finite_vec3(v: [f64; 3], context: &'static str) -> Result<()> {
    for &c in &v {
        require_finite(c, context)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn require_finite_accepts_normal() {
        assert!(require_finite(42.0, "test").is_ok());
        assert!(require_finite(0.0, "test").is_ok());
        assert!(require_finite(-1e300, "test").is_ok());
    }

    #[test]
    fn require_finite_rejects_nan() {
        let err = require_finite(f64::NAN, "position").unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("position"));
        assert!(msg.contains("NaN"));
    }

    #[test]
    fn require_finite_rejects_infinity() {
        assert!(require_finite(f64::INFINITY, "velocity").is_err());
        assert!(require_finite(f64::NEG_INFINITY, "velocity").is_err());
    }

    #[test]
    fn require_all_finite_rejects_first_bad() {
        assert!(require_all_finite(&[1.0, 2.0, 3.0], "test").is_ok());
        assert!(require_all_finite(&[1.0, f64::NAN, 3.0], "test").is_err());
    }

    #[test]
    fn ensure_finite_passes_through() {
        assert_eq!(ensure_finite(42.0, "test").unwrap(), 42.0);
    }

    #[test]
    fn ensure_finite_catches_nan() {
        assert!(ensure_finite(f64::NAN, "result").is_err());
    }

    #[test]
    fn ensure_finite_catches_infinity() {
        assert!(ensure_finite(f64::INFINITY, "result").is_err());
    }

    #[test]
    fn require_finite_vec3_ok() {
        assert!(require_finite_vec3([1.0, 2.0, 3.0], "pos").is_ok());
    }

    #[test]
    fn require_finite_vec3_catches_nan() {
        assert!(require_finite_vec3([1.0, f64::NAN, 3.0], "pos").is_err());
    }
}
