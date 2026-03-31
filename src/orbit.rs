//! Orbital elements and state vectors.

use serde::{Deserialize, Serialize};

use crate::error::{FalakError, Result};

/// Classical (Keplerian) orbital elements defining an orbit.
///
/// All angles are in radians. Semi-major axis is in the same unit system
/// as the gravitational parameter used (typically meters or kilometers).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[non_exhaustive]
pub struct OrbitalElements {
    /// Semi-major axis (distance units).
    pub semi_major_axis: f64,

    /// Eccentricity (dimensionless, 0 = circular, 0..1 = elliptical).
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
    /// - `semi_major_axis` is not positive
    /// - `eccentricity` is not in `[0, 1)` (elliptical orbits only)
    /// - `inclination` is not in `[0, pi]`
    #[must_use = "returns the validated orbital elements"]
    pub fn new(
        semi_major_axis: f64,
        eccentricity: f64,
        inclination: f64,
        raan: f64,
        argument_of_periapsis: f64,
        true_anomaly: f64,
    ) -> Result<Self> {
        if semi_major_axis <= 0.0 {
            return Err(FalakError::InvalidParameter(
                format!("semi-major axis must be positive, got {semi_major_axis}").into(),
            ));
        }

        if !(0.0..1.0).contains(&eccentricity) {
            return Err(FalakError::InvalidParameter(
                format!("eccentricity must be in [0, 1) for elliptical orbits, got {eccentricity}")
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
    #[must_use]
    #[inline]
    pub fn periapsis(&self) -> f64 {
        self.semi_major_axis * (1.0 - self.eccentricity)
    }

    /// Apoapsis distance (farthest point).
    #[must_use]
    #[inline]
    pub fn apoapsis(&self) -> f64 {
        self.semi_major_axis * (1.0 + self.eccentricity)
    }

    /// Semi-latus rectum (p = a * (1 - e^2)).
    #[must_use]
    #[inline]
    pub fn semi_latus_rectum(&self) -> f64 {
        self.semi_major_axis * (1.0 - self.eccentricity * self.eccentricity)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn valid_orbit() {
        let orbit = OrbitalElements::new(7000.0e3, 0.01, 0.9, 1.0, 0.5, 0.0);
        assert!(orbit.is_ok());
    }

    #[test]
    fn negative_sma_rejected() {
        let orbit = OrbitalElements::new(-1.0, 0.5, 0.0, 0.0, 0.0, 0.0);
        assert!(matches!(orbit, Err(FalakError::InvalidParameter(_))));
    }

    #[test]
    fn eccentricity_out_of_range() {
        let orbit = OrbitalElements::new(7000.0e3, 1.5, 0.0, 0.0, 0.0, 0.0);
        assert!(matches!(orbit, Err(FalakError::InvalidParameter(_))));
    }

    #[test]
    fn periapsis_apoapsis() {
        let orbit = OrbitalElements::new(10000.0, 0.5, 0.0, 0.0, 0.0, 0.0).unwrap();
        assert!((orbit.periapsis() - 5000.0).abs() < 1e-10);
        assert!((orbit.apoapsis() - 15000.0).abs() < 1e-10);
    }

    #[test]
    fn semi_latus_rectum() {
        let orbit = OrbitalElements::new(10000.0, 0.5, 0.0, 0.0, 0.0, 0.0).unwrap();
        assert!((orbit.semi_latus_rectum() - 7500.0).abs() < 1e-10);
    }
}
