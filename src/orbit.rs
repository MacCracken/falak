//! Orbital elements and state vectors.

use serde::{Deserialize, Serialize};
use tracing::instrument;

use crate::error::{FalakError, Result};

/// Classical (Keplerian) orbital elements defining an orbit.
///
/// All angles are in radians. Semi-major axis is in the same unit system
/// as the gravitational parameter used (typically meters or kilometers).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[non_exhaustive]
pub struct OrbitalElements {
    /// Semi-major axis (distance units). For parabolic orbits (e = 1), this
    /// field stores the semi-latus rectum *p* instead, since *a* is undefined.
    pub semi_major_axis: f64,

    /// Eccentricity (dimensionless, 0 = circular, 0..1 = elliptical,
    /// 1 = parabolic, >1 = hyperbolic).
    pub eccentricity: f64,

    /// Inclination (radians, 0..pi).
    pub inclination: f64,

    /// Right ascension of the ascending node (radians, 0..2*pi).
    pub raan: f64,

    /// Argument of periapsis (radians, 0..2*pi).
    pub argument_of_periapsis: f64,

    /// True anomaly (radians, 0..2*pi).
    pub true_anomaly: f64,
}

impl OrbitalElements {
    /// Create a new set of orbital elements with validation.
    ///
    /// # Errors
    ///
    /// Returns [`FalakError::InvalidParameter`] if:
    /// - `eccentricity` is negative
    /// - `semi_major_axis` is not positive for elliptical/circular orbits (e < 1)
    /// - `semi_major_axis` is not negative for hyperbolic orbits (e > 1)
    /// - `inclination` is not in `[0, pi]`
    #[must_use = "returns the validated orbital elements"]
    #[instrument(level = "trace")]
    pub fn new(
        semi_major_axis: f64,
        eccentricity: f64,
        inclination: f64,
        raan: f64,
        argument_of_periapsis: f64,
        true_anomaly: f64,
    ) -> Result<Self> {
        if eccentricity < 0.0 {
            return Err(FalakError::InvalidParameter(
                format!("eccentricity must be non-negative, got {eccentricity}").into(),
            ));
        }

        if eccentricity < 1.0 && semi_major_axis <= 0.0 {
            return Err(FalakError::InvalidParameter(
                format!(
                    "semi-major axis must be positive for elliptical orbits, got {semi_major_axis}"
                )
                .into(),
            ));
        }

        if eccentricity > 1.0 && semi_major_axis >= 0.0 {
            return Err(FalakError::InvalidParameter(
                format!(
                    "semi-major axis must be negative for hyperbolic orbits, got {semi_major_axis}"
                )
                .into(),
            ));
        }

        // Parabolic: e == 1, a is not meaningful (use semi-latus rectum p instead)
        // We accept any finite a for parabolic and store p as semi_major_axis by convention.
        if eccentricity == 1.0 && semi_major_axis <= 0.0 {
            return Err(FalakError::InvalidParameter(
                format!(
                    "semi-latus rectum (stored as semi_major_axis) must be positive for parabolic orbits, got {semi_major_axis}"
                )
                .into(),
            ));
        }

        if !(0.0..=std::f64::consts::PI).contains(&inclination) {
            return Err(FalakError::InvalidParameter(
                format!("inclination must be in [0, pi], got {inclination}").into(),
            ));
        }

        Ok(Self {
            semi_major_axis,
            eccentricity,
            inclination,
            raan,
            argument_of_periapsis,
            true_anomaly,
        })
    }

    /// Periapsis distance (closest approach).
    ///
    /// For parabolic orbits (e = 1), returns p/2 where p is stored as `semi_major_axis`.
    /// For hyperbolic orbits (a < 0, e > 1), returns |a|(e − 1).
    #[must_use]
    #[inline]
    pub fn periapsis(&self) -> f64 {
        if (self.eccentricity - 1.0).abs() < 1e-10 {
            // Parabolic: periapsis = p / 2
            self.semi_major_axis / 2.0
        } else {
            self.semi_major_axis.abs() * (1.0 - self.eccentricity).abs()
        }
    }

    /// Apoapsis distance (farthest point).
    ///
    /// Only meaningful for elliptical orbits (e < 1).
    /// Returns `f64::INFINITY` for parabolic and hyperbolic orbits.
    #[must_use]
    #[inline]
    pub fn apoapsis(&self) -> f64 {
        if self.eccentricity >= 1.0 {
            f64::INFINITY
        } else {
            self.semi_major_axis * (1.0 + self.eccentricity)
        }
    }

    /// Semi-latus rectum (p = |a| × |1 − e²|).
    #[must_use]
    #[inline]
    pub fn semi_latus_rectum(&self) -> f64 {
        self.semi_major_axis.abs() * (1.0 - self.eccentricity * self.eccentricity).abs()
    }

    /// Whether this orbit is elliptical (e < 1).
    #[must_use]
    #[inline]
    pub fn is_elliptical(&self) -> bool {
        self.eccentricity < 1.0
    }

    /// Whether this orbit is circular (e ≈ 0).
    #[must_use]
    #[inline]
    pub fn is_circular(&self) -> bool {
        self.eccentricity < 1e-10
    }

    /// Whether this orbit is parabolic (e = 1).
    #[must_use]
    #[inline]
    pub fn is_parabolic(&self) -> bool {
        (self.eccentricity - 1.0).abs() < 1e-10
    }

    /// Whether this orbit is hyperbolic (e > 1).
    #[must_use]
    #[inline]
    pub fn is_hyperbolic(&self) -> bool {
        self.eccentricity > 1.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn valid_elliptical() {
        let orbit = OrbitalElements::new(7000.0e3, 0.01, 0.9, 1.0, 0.5, 0.0);
        assert!(orbit.is_ok());
        assert!(orbit.unwrap().is_elliptical());
    }

    #[test]
    fn valid_hyperbolic() {
        let orbit = OrbitalElements::new(-10000.0e3, 1.5, 0.5, 0.0, 0.0, 0.5);
        assert!(orbit.is_ok());
        assert!(orbit.unwrap().is_hyperbolic());
    }

    #[test]
    fn valid_parabolic() {
        let orbit = OrbitalElements::new(7000.0e3, 1.0, 0.0, 0.0, 0.0, 0.0);
        assert!(orbit.is_ok());
        assert!(orbit.unwrap().is_parabolic());
    }

    #[test]
    fn negative_sma_rejected_elliptical() {
        let orbit = OrbitalElements::new(-1.0, 0.5, 0.0, 0.0, 0.0, 0.0);
        assert!(matches!(orbit, Err(FalakError::InvalidParameter(_))));
    }

    #[test]
    fn positive_sma_rejected_hyperbolic() {
        let orbit = OrbitalElements::new(1000.0, 1.5, 0.0, 0.0, 0.0, 0.0);
        assert!(matches!(orbit, Err(FalakError::InvalidParameter(_))));
    }

    #[test]
    fn negative_eccentricity_rejected() {
        let orbit = OrbitalElements::new(7000.0e3, -0.1, 0.0, 0.0, 0.0, 0.0);
        assert!(matches!(orbit, Err(FalakError::InvalidParameter(_))));
    }

    #[test]
    fn periapsis_apoapsis_elliptical() {
        let orbit = OrbitalElements::new(10000.0, 0.5, 0.0, 0.0, 0.0, 0.0).unwrap();
        assert!((orbit.periapsis() - 5000.0).abs() < 1e-10);
        assert!((orbit.apoapsis() - 15000.0).abs() < 1e-10);
    }

    #[test]
    fn apoapsis_hyperbolic_infinite() {
        let orbit = OrbitalElements::new(-10000.0, 1.5, 0.0, 0.0, 0.0, 0.0).unwrap();
        assert!(orbit.apoapsis().is_infinite());
        assert!(orbit.periapsis() > 0.0);
    }

    #[test]
    fn semi_latus_rectum() {
        let orbit = OrbitalElements::new(10000.0, 0.5, 0.0, 0.0, 0.0, 0.0).unwrap();
        assert!((orbit.semi_latus_rectum() - 7500.0).abs() < 1e-10);
    }

    #[test]
    fn periapsis_parabolic() {
        // Parabolic: semi_major_axis stores p, periapsis = p/2
        let orbit = OrbitalElements::new(14000.0, 1.0, 0.0, 0.0, 0.0, 0.0).unwrap();
        assert!(
            (orbit.periapsis() - 7000.0).abs() < 1e-10,
            "parabolic periapsis: {}",
            orbit.periapsis()
        );
    }

    #[test]
    fn periapsis_hyperbolic() {
        // Hyperbolic: a = -10000, e = 1.5 → periapsis = |a| × |1-e| = 10000 × 0.5 = 5000
        let orbit = OrbitalElements::new(-10000.0, 1.5, 0.0, 0.0, 0.0, 0.0).unwrap();
        assert!(
            (orbit.periapsis() - 5000.0).abs() < 1e-10,
            "hyperbolic periapsis: {}",
            orbit.periapsis()
        );
    }

    #[test]
    fn orbit_type_queries() {
        let circ = OrbitalElements::new(7000.0, 0.0, 0.0, 0.0, 0.0, 0.0).unwrap();
        assert!(circ.is_circular());
        assert!(circ.is_elliptical());
        assert!(!circ.is_hyperbolic());
        assert!(!circ.is_parabolic());
    }
}
