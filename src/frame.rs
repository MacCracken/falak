//! Reference frames — ECI, ECEF, perifocal, rotating frames, coordinate transforms.
//!
//! Provides transformations between common orbital mechanics reference frames.
//! All positions are in metres, angles in radians.

use tracing::instrument;

use crate::error::{FalakError, Result};

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

    // Transpose of the perifocal→ECI matrix (Rᵀ = R⁻¹ for rotation)
    let r11 = cos_o * cos_w - sin_o * sin_w * cos_i;
    let r21 = -(cos_o * sin_w + sin_o * cos_w * cos_i);
    let r12 = sin_o * cos_w + cos_o * sin_w * cos_i;
    let r22 = -(sin_o * sin_w - cos_o * cos_w * cos_i);
    let r13 = sin_w * sin_i;
    let r23 = cos_w * sin_i;
    let r33 = cos_i;

    [
        r11 * eci[0] + r12 * eci[1] + r13 * eci[2],
        r21 * eci[0] + r22 * eci[1] + r23 * eci[2],
        -sin_i * sin_o * eci[0] + sin_i * cos_o * eci[1] + r33 * eci[2],
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

// ── Precession ──────────────────────────────────────────────────────────

/// Compute IAU 2006 precession angles (ζ, θ, z) at a given Julian date.
///
/// These three angles define the precession rotation from the mean equator
/// and equinox of J2000.0 to the mean equator and equinox of date.
///
/// Returns `(zeta, theta, z)` in radians.
///
/// # Arguments
///
/// * `jd` — Julian date
///
/// # Reference
///
/// Lieske et al. (1977), IAU 2006 values from Capitaine et al. (2003).
#[must_use]
pub fn precession_angles(jd: f64) -> (f64, f64, f64) {
    let t = (jd - 2_451_545.0) / 36525.0; // Julian centuries from J2000
    let t2 = t * t;
    let t3 = t2 * t;

    // IAU 2006 precession angles in arcseconds
    let zeta_as = 2.650_545 + 2_306.083_227 * t + 0.299_3 * t2 + 0.017_998 * t3;
    let theta_as = 2_004.191_903 * t - 0.429_4 * t2 - 0.041_82 * t3;
    let z_as = -2.650_545 + 2_306.077_181 * t + 1.094_68 * t2 + 0.018_203 * t3;

    let as_to_rad = std::f64::consts::PI / (180.0 * 3600.0);
    (zeta_as * as_to_rad, theta_as * as_to_rad, z_as * as_to_rad)
}

/// Apply the precession rotation to transform a position from J2000.0
/// mean equator/equinox to the mean equator/equinox of date.
///
/// P = R₃(−z) · R₂(θ) · R₃(−ζ)
#[must_use]
pub fn precess_j2000_to_date(position: [f64; 3], jd: f64) -> [f64; 3] {
    let (zeta, theta, z) = precession_angles(jd);

    // R₃(−ζ)
    let p1 = rotate_z(position, -zeta);
    // R₂(θ)
    let p2 = rotate_y(p1, theta);
    // R₃(−z)
    rotate_z(p2, -z)
}

/// Inverse: precess from mean equator/equinox of date back to J2000.0.
#[must_use]
pub fn precess_date_to_j2000(position: [f64; 3], jd: f64) -> [f64; 3] {
    let (zeta, theta, z) = precession_angles(jd);

    // Inverse: R₃(ζ) · R₂(−θ) · R₃(z)
    let p1 = rotate_z(position, z);
    let p2 = rotate_y(p1, -theta);
    rotate_z(p2, zeta)
}

// ── Nutation ────────────────────────────────────────────────────────────

/// Compute nutation in longitude (Δψ) and nutation in obliquity (Δε) using
/// a truncated IAU 1980 nutation series (9 dominant terms).
///
/// Returns `(delta_psi, delta_epsilon)` in radians.
///
/// Accuracy: ~0.5 arcsecond (sufficient for most orbital mechanics applications).
///
/// # Arguments
///
/// * `jd` — Julian date
///
/// # Reference
///
/// Seidelmann (1982), truncated from the full 106-term series.
#[must_use]
pub fn nutation(jd: f64) -> (f64, f64) {
    let t = (jd - 2_451_545.0) / 36525.0;

    // Fundamental arguments (Delaunay variables) in radians
    // Mean anomaly of the Moon
    let l = (134.963_41 + 477_198.867_63 * t).to_radians();
    // Mean anomaly of the Sun
    let lp = (357.529_11 + 35_999.050_29 * t).to_radians();
    // Mean argument of latitude of the Moon
    let f = (93.271_91 + 483_202.017_53 * t).to_radians();
    // Mean elongation of the Moon from the Sun
    let d = (297.850_36 + 445_267.111_48 * t).to_radians();
    // Longitude of ascending node of the Moon
    let omega = (125.044_52 - 1_934.136_61 * t).to_radians();

    // Truncated nutation series: 9 dominant terms
    // Each row: [l_mult, lp_mult, f_mult, d_mult, omega_mult, dpsi_sin_coeff_as, deps_cos_coeff_as]
    #[rustfmt::skip]
    let terms: [(f64, f64, f64, f64, f64, f64, f64); 9] = [
        (0.0, 0.0, 0.0, 0.0, 1.0,  -17.2064,  9.2052),  // Ω
        (0.0, 0.0, 2.0, -2.0, 2.0, -1.3171,  0.5736),   // 2F-2D+2Ω
        (0.0, 0.0, 2.0, 0.0, 2.0,  -0.2276,  0.0977),   // 2F+2Ω
        (0.0, 0.0, 0.0, 0.0, 2.0,   0.2075, -0.0895),   // 2Ω
        (0.0, 1.0, 0.0, 0.0, 0.0,   0.1476, -0.0001),   // l'
        (1.0, 0.0, 0.0, 0.0, 0.0,  -0.0517,  0.0224),   // l
        (0.0, 0.0, 2.0, 0.0, 1.0,  -0.0387,  0.0200),   // 2F+Ω
        (0.0, 0.0, 2.0, -2.0, 1.0, -0.0301,  0.0129),   // 2F-2D+Ω
        (1.0, 0.0, 0.0, 0.0, 1.0,   0.0217, -0.0095),   // l+Ω
    ];

    let mut dpsi = 0.0; // nutation in longitude (arcsec)
    let mut deps = 0.0; // nutation in obliquity (arcsec)

    for &(lm, lpm, fm, dm, om, sp, ce) in &terms {
        let arg = lm * l + lpm * lp + fm * f + dm * d + om * omega;
        dpsi += sp * arg.sin();
        deps += ce * arg.cos();
    }

    let as_to_rad = std::f64::consts::PI / (180.0 * 3600.0);
    (dpsi * as_to_rad, deps * as_to_rad)
}

/// Mean obliquity of the ecliptic (radians) at a given Julian date.
///
/// IAU 2006 formula (Hilton et al. 2006).
#[must_use]
pub fn mean_obliquity(jd: f64) -> f64 {
    let t = (jd - 2_451_545.0) / 36525.0;
    let t2 = t * t;
    let t3 = t2 * t;

    // In arcseconds
    let eps_as = 84_381.406 - 46.836_769 * t - 0.000_183_1 * t2 + 0.002_003_4 * t3;
    eps_as * std::f64::consts::PI / (180.0 * 3600.0)
}

/// True obliquity of the ecliptic = mean obliquity + Δε.
#[must_use]
#[inline]
pub fn true_obliquity(jd: f64) -> f64 {
    let (_, deps) = nutation(jd);
    mean_obliquity(jd) + deps
}

/// Equation of the equinoxes: Δψ · cos(ε).
///
/// This is the difference between apparent and mean sidereal time (radians).
#[must_use]
#[inline]
pub fn equation_of_equinoxes(jd: f64) -> f64 {
    let (dpsi, _) = nutation(jd);
    dpsi * mean_obliquity(jd).cos()
}

// ── Rotation helpers ────────────────────────────────────────────────────

/// Rotate a 3D vector about the Z-axis by angle α.
#[must_use]
#[inline]
fn rotate_z(v: [f64; 3], alpha: f64) -> [f64; 3] {
    let (s, c) = alpha.sin_cos();
    [c * v[0] - s * v[1], s * v[0] + c * v[1], v[2]]
}

/// Rotate a 3D vector about the Y-axis by angle α.
#[must_use]
#[inline]
fn rotate_y(v: [f64; 3], alpha: f64) -> [f64; 3] {
    let (s, c) = alpha.sin_cos();
    [c * v[0] + s * v[2], v[1], -s * v[0] + c * v[2]]
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::{FRAC_PI_2, TAU};

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

    #[test]
    fn eci_to_perifocal_out_of_plane() {
        // For a polar orbit (i=π/2, Ω=0, ω=0), an ECI z-component should
        // map to the Q axis (perifocal y), not be dropped.
        let eci = [0.0, 0.0, 7e6]; // pure z in ECI
        let pqw = eci_to_perifocal(eci, 0.0, FRAC_PI_2, 0.0);
        // For i=π/2: W = cos(i)*z = 0, Q = cos(ω)*sin(i)*z = z
        // Actually the mapping is: the perifocal Q-axis maps to ECI z for polar orbit
        // So ECI z should map back to PQW Q
        assert!(
            (pqw[1] - 7e6).abs() < 1e-4,
            "out-of-plane should map to Q: {:?}",
            pqw
        );
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

    // ── Precession ──────────────────────────────────────────────────────

    const J2000: f64 = 2_451_545.0;

    #[test]
    fn precession_at_j2000_is_zero() {
        let (zeta, theta, z) = precession_angles(J2000);
        // At J2000, T=0, so angles should be near zero (small constant terms)
        assert!(zeta.abs() < 1e-4, "zeta at J2000: {zeta}");
        assert!(theta.abs() < 1e-6, "theta at J2000: {theta}");
        assert!(z.abs() < 1e-4, "z at J2000: {z}");
    }

    #[test]
    fn precession_grows_with_time() {
        let (_, theta1, _) = precession_angles(J2000 + 365.25 * 50.0); // 50 years
        let (_, theta2, _) = precession_angles(J2000 + 365.25 * 100.0); // 100 years
        assert!(
            theta2.abs() > theta1.abs(),
            "precession should grow with time"
        );
    }

    #[test]
    fn precession_roundtrip() {
        let pos = [1e7, 2e6, 3e6];
        let jd = J2000 + 365.25 * 50.0; // 50 years from J2000
        let precessed = precess_j2000_to_date(pos, jd);
        let back = precess_date_to_j2000(precessed, jd);
        for i in 0..3 {
            assert!(
                (back[i] - pos[i]).abs() < 1e-2,
                "roundtrip [{i}]: {} vs {}",
                back[i],
                pos[i]
            );
        }
    }

    #[test]
    fn precession_magnitude_50_years() {
        // In 50 years, θ ≈ 50 × 20.04"/yr ≈ 1002" ≈ 0.278°
        // Total precession (ζ+z) ≈ 50 × 46.1"/yr ≈ 0.64°
        let (zeta, theta, z) = precession_angles(J2000 + 365.25 * 50.0);
        let theta_deg = theta.to_degrees().abs();
        assert!(
            theta_deg > 0.2 && theta_deg < 0.4,
            "50-year θ = {theta_deg}°, expected ~0.28°"
        );
        // ζ and z should each be ~0.32°
        let zeta_deg = zeta.to_degrees().abs();
        let z_deg = z.to_degrees().abs();
        assert!(zeta_deg > 0.2 && zeta_deg < 0.5, "50-year ζ = {zeta_deg}°");
        assert!(z_deg > 0.2 && z_deg < 0.5, "50-year z = {z_deg}°");
    }

    // ── Nutation ────────────────────────────────────────────────────────

    #[test]
    fn nutation_magnitude() {
        // Nutation in longitude is typically ±17" max, ±9" in obliquity
        let (dpsi, deps) = nutation(J2000);
        let dpsi_as = dpsi.to_degrees() * 3600.0;
        let deps_as = deps.to_degrees() * 3600.0;
        assert!(
            dpsi_as.abs() < 20.0,
            "nutation in longitude = {dpsi_as:.1}\", expected < 20\""
        );
        assert!(
            deps_as.abs() < 12.0,
            "nutation in obliquity = {deps_as:.1}\", expected < 12\""
        );
    }

    #[test]
    fn nutation_varies() {
        // Nutation should differ at different times (18.6 year period)
        let (dpsi1, _) = nutation(J2000);
        let (dpsi2, _) = nutation(J2000 + 365.25 * 9.3); // ~half nutation period
        assert!(
            (dpsi1 - dpsi2).abs() > 1e-6,
            "nutation should vary over time"
        );
    }

    #[test]
    fn mean_obliquity_j2000() {
        // Mean obliquity at J2000 ≈ 23.439°
        let eps = mean_obliquity(J2000);
        let deg = eps.to_degrees();
        assert!(
            (deg - 23.439).abs() < 0.01,
            "mean obliquity = {deg}°, expected ~23.439°"
        );
    }

    #[test]
    fn true_obliquity_differs_from_mean() {
        let eps_mean = mean_obliquity(J2000);
        let eps_true = true_obliquity(J2000);
        // Difference should be the nutation in obliquity (< 10")
        let diff_as = (eps_true - eps_mean).to_degrees() * 3600.0;
        assert!(
            diff_as.abs() < 12.0,
            "true-mean obliquity diff = {diff_as:.1}\", expected < 12\""
        );
    }

    #[test]
    fn equation_of_equinoxes_magnitude() {
        let eqeq = equation_of_equinoxes(J2000);
        let sec = eqeq.to_degrees() * 3600.0;
        // EqEq ≈ Δψ cos(ε) ≈ ±17" × cos(23.4°) ≈ ±15.6"
        assert!(
            sec.abs() < 20.0,
            "equation of equinoxes = {sec:.1}\", expected < 20\""
        );
    }
}
