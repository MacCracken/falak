//! Reference frames — ECI, ECEF, perifocal, rotating frames, coordinate transforms.
//!
//! Provides transformations between common orbital mechanics reference frames.
//! All positions are in metres, angles in radians.

use tracing::instrument;

use crate::error::{FalakError, Result};

/// A 3D position in a specific reference frame.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
#[non_exhaustive]
pub struct Position {
    /// Coordinates `[x, y, z]` (metres).
    pub coords: [f64; 3],
}

impl Position {
    /// Create a new position.
    #[must_use]
    #[inline]
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { coords: [x, y, z] }
    }

    /// Magnitude (distance from origin).
    #[must_use]
    #[inline]
    pub fn magnitude(&self) -> f64 {
        let [x, y, z] = self.coords;
        (x * x + y * y + z * z).sqrt()
    }
}

// ── Perifocal ↔ ECI ───────────────────────────────────────────────────────

/// Transform a position from perifocal (PQW) frame to ECI frame.
///
/// Uses the standard rotation matrix R = R₃(-Ω) R₁(-i) R₃(-ω).
///
/// # Arguments
///
/// * `pqw` — Position in perifocal frame `[P, Q, W]`
/// * `raan` — Right ascension of ascending node Ω (radians)
/// * `inclination` — Orbital inclination i (radians)
/// * `arg_periapsis` — Argument of periapsis ω (radians)
#[must_use]
#[inline]
pub fn perifocal_to_eci(
    pqw: [f64; 3],
    raan: f64,
    inclination: f64,
    arg_periapsis: f64,
) -> [f64; 3] {
    let (cos_o, sin_o) = (raan.cos(), raan.sin());
    let (cos_w, sin_w) = (arg_periapsis.cos(), arg_periapsis.sin());
    let (cos_i, sin_i) = (inclination.cos(), inclination.sin());

    let r11 = cos_o * cos_w - sin_o * sin_w * cos_i;
    let r12 = -(cos_o * sin_w + sin_o * cos_w * cos_i);
    let r21 = sin_o * cos_w + cos_o * sin_w * cos_i;
    let r22 = -(sin_o * sin_w - cos_o * cos_w * cos_i);
    let r31 = sin_w * sin_i;
    let r32 = cos_w * sin_i;

    [
        r11 * pqw[0] + r12 * pqw[1],
        r21 * pqw[0] + r22 * pqw[1],
        r31 * pqw[0] + r32 * pqw[1],
    ]
}

/// Transform a position from ECI to perifocal (PQW) frame.
///
/// Inverse of [`perifocal_to_eci`] — uses the transpose of the rotation matrix.
#[must_use]
#[inline]
pub fn eci_to_perifocal(
    eci: [f64; 3],
    raan: f64,
    inclination: f64,
    arg_periapsis: f64,
) -> [f64; 3] {
    let (cos_o, sin_o) = (raan.cos(), raan.sin());
    let (cos_w, sin_w) = (arg_periapsis.cos(), arg_periapsis.sin());
    let (cos_i, sin_i) = (inclination.cos(), inclination.sin());

    // Transpose of the perifocal→ECI matrix
    let r11 = cos_o * cos_w - sin_o * sin_w * cos_i;
    let r21 = -(cos_o * sin_w + sin_o * cos_w * cos_i);
    let r12 = sin_o * cos_w + cos_o * sin_w * cos_i;
    let r22 = -(sin_o * sin_w - cos_o * cos_w * cos_i);
    let r13 = sin_w * sin_i;
    let r23 = cos_w * sin_i;

    [
        r11 * eci[0] + r12 * eci[1] + r13 * eci[2],
        r21 * eci[0] + r22 * eci[1] + r23 * eci[2],
        0.0, // W component is always 0 for orbits in the orbital plane
    ]
}

// ── ECI ↔ ECEF ────────────────────────────────────────────────────────────

/// Transform a position from ECI to ECEF (Earth-Centered Earth-Fixed).
///
/// Applies a rotation about the Z-axis by the Greenwich Mean Sidereal Time (GMST).
///
/// # Arguments
///
/// * `eci` — Position in ECI frame `[x, y, z]` (metres)
/// * `gmst` — Greenwich Mean Sidereal Time (radians)
#[must_use]
#[inline]
pub fn eci_to_ecef(eci: [f64; 3], gmst: f64) -> [f64; 3] {
    let (cos_g, sin_g) = (gmst.cos(), gmst.sin());
    [
        cos_g * eci[0] + sin_g * eci[1],
        -sin_g * eci[0] + cos_g * eci[1],
        eci[2],
    ]
}

/// Transform a position from ECEF to ECI.
///
/// Inverse of [`eci_to_ecef`] — rotates by −GMST about Z.
#[must_use]
#[inline]
pub fn ecef_to_eci(ecef: [f64; 3], gmst: f64) -> [f64; 3] {
    let (cos_g, sin_g) = (gmst.cos(), gmst.sin());
    [
        cos_g * ecef[0] - sin_g * ecef[1],
        sin_g * ecef[0] + cos_g * ecef[1],
        ecef[2],
    ]
}

// ── ECEF ↔ Geodetic ──────────────────────────────────────────────────────

/// WGS-84 semi-major axis (equatorial radius, metres).
pub const WGS84_A: f64 = 6_378_137.0;

/// WGS-84 flattening.
pub const WGS84_F: f64 = 1.0 / 298.257_223_563;

/// WGS-84 eccentricity squared.
pub const WGS84_E2: f64 = 2.0 * WGS84_F - WGS84_F * WGS84_F;

/// Geodetic coordinates (latitude, longitude, altitude).
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
#[non_exhaustive]
pub struct Geodetic {
    /// Latitude (radians, −π/2 to π/2).
    pub latitude: f64,
    /// Longitude (radians, −π to π).
    pub longitude: f64,
    /// Altitude above the WGS-84 ellipsoid (metres).
    pub altitude: f64,
}

/// Convert ECEF coordinates to geodetic (lat, lon, alt) using Bowring's method.
///
/// Uses iterative Bowring's method for WGS-84 ellipsoid.
///
/// # Errors
///
/// Returns [`FalakError::InvalidParameter`] if position is at the origin.
#[must_use = "returns the computed geodetic coordinates"]
#[instrument(level = "trace")]
pub fn ecef_to_geodetic(ecef: [f64; 3]) -> Result<Geodetic> {
    let [x, y, z] = ecef;
    let p = (x * x + y * y).sqrt();

    if p < 1e-30 && z.abs() < 1e-30 {
        return Err(FalakError::InvalidParameter(
            "ECEF position is at the origin".into(),
        ));
    }

    let longitude = y.atan2(x);

    // Bowring's iterative method
    let b = WGS84_A * (1.0 - WGS84_F);
    let e2 = WGS84_E2;
    let ep2 = (WGS84_A * WGS84_A - b * b) / (b * b);

    // Initial estimate
    let theta = z.atan2(p * (1.0 - WGS84_F));
    let mut lat = (z + ep2 * b * theta.sin().powi(3)).atan2(p - e2 * WGS84_A * theta.cos().powi(3));

    // Iterate (typically converges in 2-3 iterations)
    for _ in 0..10 {
        let sin_lat = lat.sin();
        let n = WGS84_A / (1.0 - e2 * sin_lat * sin_lat).sqrt();
        let lat_new = (z + e2 * n * sin_lat).atan2(p);
        if (lat_new - lat).abs() < 1e-12 {
            lat = lat_new;
            break;
        }
        lat = lat_new;
    }

    let sin_lat = lat.sin();
    let n = WGS84_A / (1.0 - e2 * sin_lat * sin_lat).sqrt();
    let altitude = if lat.cos().abs() > 1e-10 {
        p / lat.cos() - n
    } else {
        z.abs() / sin_lat.abs() - n * (1.0 - e2)
    };

    Ok(Geodetic {
        latitude: lat,
        longitude,
        altitude,
    })
}

/// Convert geodetic coordinates to ECEF.
#[must_use]
#[inline]
pub fn geodetic_to_ecef(geodetic: &Geodetic) -> [f64; 3] {
    let sin_lat = geodetic.latitude.sin();
    let cos_lat = geodetic.latitude.cos();
    let sin_lon = geodetic.longitude.sin();
    let cos_lon = geodetic.longitude.cos();

    let n = WGS84_A / (1.0 - WGS84_E2 * sin_lat * sin_lat).sqrt();

    [
        (n + geodetic.altitude) * cos_lat * cos_lon,
        (n + geodetic.altitude) * cos_lat * sin_lon,
        (n * (1.0 - WGS84_E2) + geodetic.altitude) * sin_lat,
    ]
}

// ── Rotating (synodic) frame ──────────────────────────────────────────────

/// Transform from inertial to rotating (synodic) frame.
///
/// The rotating frame co-rotates with a reference body at angular rate ω.
/// Used for CR3BP (Circular Restricted Three-Body Problem) analysis.
///
/// # Arguments
///
/// * `inertial` — Position in inertial frame `[x, y, z]`
/// * `angle` — Rotation angle ωt (radians)
#[must_use]
#[inline]
pub fn inertial_to_rotating(inertial: [f64; 3], angle: f64) -> [f64; 3] {
    let (cos_a, sin_a) = (angle.cos(), angle.sin());
    [
        cos_a * inertial[0] + sin_a * inertial[1],
        -sin_a * inertial[0] + cos_a * inertial[1],
        inertial[2],
    ]
}

/// Transform from rotating (synodic) frame to inertial frame.
#[must_use]
#[inline]
pub fn rotating_to_inertial(rotating: [f64; 3], angle: f64) -> [f64; 3] {
    let (cos_a, sin_a) = (angle.cos(), angle.sin());
    [
        cos_a * rotating[0] - sin_a * rotating[1],
        sin_a * rotating[0] + cos_a * rotating[1],
        rotating[2],
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::{FRAC_PI_2, TAU};

    // ── Position ─────────────────────────────────────────────────────

    #[test]
    fn position_magnitude() {
        let p = Position::new(3.0, 4.0, 0.0);
        assert!((p.magnitude() - 5.0).abs() < 1e-15);
    }

    // ── Perifocal ↔ ECI ──────────────────────────────────────────────

    #[test]
    fn perifocal_eci_identity() {
        // Zero angles → perifocal = ECI
        let pqw = [7e6, 0.0, 0.0];
        let eci = perifocal_to_eci(pqw, 0.0, 0.0, 0.0);
        assert!((eci[0] - 7e6).abs() < 1e-6);
        assert!(eci[1].abs() < 1e-6);
        assert!(eci[2].abs() < 1e-6);
    }

    #[test]
    fn perifocal_eci_roundtrip() {
        let pqw = [7e6, 3e6, 0.0];
        let raan = 1.2;
        let inc = 0.8;
        let aop = 0.5;
        let eci = perifocal_to_eci(pqw, raan, inc, aop);
        let pqw2 = eci_to_perifocal(eci, raan, inc, aop);
        assert!(
            (pqw[0] - pqw2[0]).abs() < 1e-4,
            "P: {} vs {}",
            pqw[0],
            pqw2[0]
        );
        assert!(
            (pqw[1] - pqw2[1]).abs() < 1e-4,
            "Q: {} vs {}",
            pqw[1],
            pqw2[1]
        );
    }

    #[test]
    fn perifocal_eci_polar() {
        // i = π/2, Ω = 0, ω = 0 → perifocal x maps to ECI x, y maps to ECI z
        let pqw = [7e6, 0.0, 0.0];
        let eci = perifocal_to_eci(pqw, 0.0, FRAC_PI_2, 0.0);
        assert!((eci[0] - 7e6).abs() < 1e-4);
        assert!(eci[1].abs() < 1e-4);
        assert!(eci[2].abs() < 1e-4);

        let pqw2 = [0.0, 7e6, 0.0];
        let eci2 = perifocal_to_eci(pqw2, 0.0, FRAC_PI_2, 0.0);
        assert!(eci2[0].abs() < 1e-4);
        assert!(eci2[1].abs() < 1e-4);
        assert!((eci2[2] - 7e6).abs() < 1e-4, "polar Q→Z: {}", eci2[2]);
    }

    // ── ECI ↔ ECEF ───────────────────────────────────────────────────

    #[test]
    fn eci_ecef_zero_gmst() {
        let eci = [7e6, 0.0, 0.0];
        let ecef = eci_to_ecef(eci, 0.0);
        assert!((ecef[0] - 7e6).abs() < 1e-6);
        assert!(ecef[1].abs() < 1e-6);
    }

    #[test]
    fn eci_ecef_roundtrip() {
        let eci = [7e6, 3e6, 1e6];
        let gmst = 1.5;
        let ecef = eci_to_ecef(eci, gmst);
        let eci2 = ecef_to_eci(ecef, gmst);
        for j in 0..3 {
            assert!(
                (eci[j] - eci2[j]).abs() < 1e-4,
                "axis {j}: {} vs {}",
                eci[j],
                eci2[j]
            );
        }
    }

    #[test]
    fn eci_ecef_quarter_turn() {
        // GMST = π/2 → ECI x-axis maps to ECEF y-axis
        let eci = [7e6, 0.0, 0.0];
        let ecef = eci_to_ecef(eci, FRAC_PI_2);
        assert!(ecef[0].abs() < 1e-4);
        assert!((ecef[1] - (-7e6)).abs() < 1e-4);
    }

    // ── ECEF ↔ Geodetic ──────────────────────────────────────────────

    #[test]
    fn geodetic_roundtrip_equator() {
        let geo = Geodetic {
            latitude: 0.0,
            longitude: 0.0,
            altitude: 0.0,
        };
        let ecef = geodetic_to_ecef(&geo);
        assert!((ecef[0] - WGS84_A).abs() < 1.0, "equator x: {}", ecef[0]);
        assert!(ecef[1].abs() < 1e-6);
        assert!(ecef[2].abs() < 1e-6);

        let geo2 = ecef_to_geodetic(ecef).unwrap();
        assert!((geo2.latitude - geo.latitude).abs() < 1e-10);
        assert!((geo2.longitude - geo.longitude).abs() < 1e-10);
        assert!(geo2.altitude.abs() < 1.0, "alt: {}", geo2.altitude);
    }

    #[test]
    fn geodetic_roundtrip_pole() {
        let geo = Geodetic {
            latitude: FRAC_PI_2,
            longitude: 0.0,
            altitude: 0.0,
        };
        let ecef = geodetic_to_ecef(&geo);
        let geo2 = ecef_to_geodetic(ecef).unwrap();
        assert!(
            (geo2.latitude - FRAC_PI_2).abs() < 1e-10,
            "pole lat: {}",
            geo2.latitude
        );
        assert!(geo2.altitude.abs() < 1.0, "pole alt: {}", geo2.altitude);
    }

    #[test]
    fn geodetic_roundtrip_arbitrary() {
        let geo = Geodetic {
            latitude: 0.7, // ~40° N
            longitude: -1.3,
            altitude: 35786.0e3, // GEO altitude
        };
        let ecef = geodetic_to_ecef(&geo);
        let geo2 = ecef_to_geodetic(ecef).unwrap();
        assert!(
            (geo2.latitude - geo.latitude).abs() < 1e-8,
            "lat: {} vs {}",
            geo.latitude,
            geo2.latitude
        );
        assert!(
            (geo2.longitude - geo.longitude).abs() < 1e-8,
            "lon: {} vs {}",
            geo.longitude,
            geo2.longitude
        );
        assert!(
            (geo2.altitude - geo.altitude).abs() < 1.0,
            "alt: {} vs {}",
            geo.altitude,
            geo2.altitude
        );
    }

    #[test]
    fn geodetic_origin_error() {
        assert!(ecef_to_geodetic([0.0, 0.0, 0.0]).is_err());
    }

    // ── Rotating frame ───────────────────────────────────────────────

    #[test]
    fn rotating_roundtrip() {
        let inertial = [5.0, 3.0, 1.0];
        let angle = 1.2;
        let rot = inertial_to_rotating(inertial, angle);
        let back = rotating_to_inertial(rot, angle);
        for j in 0..3 {
            assert!(
                (inertial[j] - back[j]).abs() < 1e-12,
                "axis {j}: {} vs {}",
                inertial[j],
                back[j]
            );
        }
    }

    #[test]
    fn rotating_full_revolution() {
        let pos = [7e6, 0.0, 0.0];
        let rot = inertial_to_rotating(pos, TAU);
        assert!((rot[0] - 7e6).abs() < 1e-4);
        assert!(rot[1].abs() < 1e-4);
    }

    #[test]
    fn rotating_quarter_turn() {
        let pos = [1.0, 0.0, 0.0];
        let rot = inertial_to_rotating(pos, FRAC_PI_2);
        assert!(rot[0].abs() < 1e-15);
        assert!((rot[1] - (-1.0)).abs() < 1e-15);
    }
}
