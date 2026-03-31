//! Orbital perturbations — J2 oblateness, atmospheric drag, solar radiation pressure,
//! third-body effects, and general perturbation theory.
//!
//! All functions return acceleration vectors `[ax, ay, az]` in m/s² that can be
//! added to the two-body acceleration for orbit propagation.

use tracing::instrument;

use crate::error::{FalakError, Result};

/// Earth's J2 zonal harmonic coefficient.
pub const J2_EARTH: f64 = 1.082_63e-3;

/// Earth's J3 zonal harmonic coefficient.
pub const J3_EARTH: f64 = -2.532_6e-6;

/// Earth's equatorial radius (metres) — identical to [`crate::frame::WGS84_A`].
pub const R_EARTH: f64 = 6_378_137.0;

/// Standard gravitational parameter of the Sun (m³/s²).
pub const MU_SUN: f64 = 1.327_124_400_18e20;

/// Solar radiation pressure at 1 AU (N/m²).
pub const SOLAR_PRESSURE_1AU: f64 = 4.56e-6;

/// 1 AU in metres.
pub const AU_METRES: f64 = 1.495_978_707e11;

// ── J2 oblateness ─────────────────────────────────────────────────────────

/// Compute the J2 oblateness perturbation acceleration.
///
/// Models the dominant effect of Earth's equatorial bulge on satellite orbits.
/// Causes nodal regression (RAAN drift) and apsidal rotation.
///
/// # Arguments
///
/// * `position` — Satellite position `[x, y, z]` in ECI (metres)
/// * `mu` — Central body gravitational parameter (m³/s²)
/// * `j2` — J2 zonal harmonic coefficient
/// * `r_body` — Central body equatorial radius (metres)
///
/// # Returns
///
/// Perturbation acceleration `[ax, ay, az]` in m/s².
#[must_use]
#[inline]
pub fn j2_acceleration(position: [f64; 3], mu: f64, j2: f64, r_body: f64) -> [f64; 3] {
    let [x, y, z] = position;
    let r2 = x * x + y * y + z * z;
    let r = r2.sqrt();
    let r5 = r2 * r2 * r;

    if r5 < 1e-30 {
        return [0.0, 0.0, 0.0];
    }

    let coeff = -1.5 * j2 * mu * r_body * r_body / r5;
    let z2_r2 = z * z / r2;

    [
        coeff * x * (1.0 - 5.0 * z2_r2),
        coeff * y * (1.0 - 5.0 * z2_r2),
        coeff * z * (3.0 - 5.0 * z2_r2),
    ]
}

/// Compute secular J2 drift rates for RAAN and argument of periapsis.
///
/// # Returns
///
/// `(d_raan_dt, d_omega_dt)` in rad/s.
///
/// # Errors
///
/// Returns [`FalakError::InvalidParameter`] if inputs are non-positive.
#[must_use = "returns (RAAN rate, arg periapsis rate)"]
#[instrument(level = "trace")]
pub fn j2_secular_rates(
    semi_major_axis: f64,
    eccentricity: f64,
    inclination: f64,
    mu: f64,
    j2: f64,
    r_body: f64,
) -> Result<(f64, f64)> {
    if semi_major_axis <= 0.0 {
        return Err(FalakError::InvalidParameter(
            format!("semi-major axis must be positive, got {semi_major_axis}").into(),
        ));
    }
    if mu <= 0.0 {
        return Err(FalakError::InvalidParameter(
            format!("mu must be positive, got {mu}").into(),
        ));
    }

    let a = semi_major_axis;
    let e = eccentricity;
    let n = (mu / (a * a * a)).sqrt(); // mean motion
    let p = a * (1.0 - e * e);

    if p <= 0.0 {
        return Err(FalakError::InvalidParameter(
            "semi-latus rectum must be positive".into(),
        ));
    }

    let ratio = r_body / p;
    let common = -1.5 * n * j2 * ratio * ratio;

    // RAAN rate: dΩ/dt = -1.5 n J₂ (R/p)² cos(i)
    let d_raan_dt = common * inclination.cos();

    // Arg periapsis rate: dω/dt = -1.5 n J₂ (R/p)² (2.5 sin²(i) - 2)
    // = 1.5 n J₂ (R/p)² (2 - 2.5 sin²(i))
    let d_omega_dt = -common * (2.0 - 2.5 * inclination.sin().powi(2));

    Ok((d_raan_dt, d_omega_dt))
}

// ── Higher-order zonal harmonics ──────────────────────────────────────────

/// Compute J3 perturbation acceleration.
///
/// J3 is the "pear-shaped" asymmetry — smaller than J2 but causes
/// eccentricity perturbations in near-polar orbits.
#[must_use]
#[inline]
pub fn j3_acceleration(position: [f64; 3], mu: f64, j3: f64, r_body: f64) -> [f64; 3] {
    let [x, y, z] = position;
    let r2 = x * x + y * y + z * z;
    let r = r2.sqrt();
    let r7 = r2 * r2 * r2 * r;

    if r7 < 1e-30 {
        return [0.0, 0.0, 0.0];
    }

    let re3 = r_body * r_body * r_body;
    let coeff = -2.5 * j3 * mu * re3 / r7;
    let z2_r2 = z * z / r2;

    [
        coeff * x * z * (3.0 - 7.0 * z2_r2),
        coeff * y * z * (3.0 - 7.0 * z2_r2),
        coeff * (z * z * (6.0 - 7.0 * z2_r2) - 0.6 * r2),
    ]
}

// ── Atmospheric drag ──────────────────────────────────────────────────────

/// Parameters for an exponential atmosphere model.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
#[non_exhaustive]
pub struct AtmosphereParams {
    /// Drag coefficient (typically 2.0–2.5).
    pub cd: f64,
    /// Cross-sectional area / mass (m²/kg).
    pub area_mass_ratio: f64,
    /// Reference atmospheric density (kg/m³).
    pub rho0: f64,
    /// Reference altitude (metres).
    pub h0: f64,
    /// Scale height (metres).
    pub h_scale: f64,
    /// Central body radius (metres).
    pub r_body: f64,
}

/// Compute atmospheric drag acceleration.
///
/// Uses an exponential atmosphere model. Drag opposes the velocity vector.
///
/// a_drag = −0.5 × ρ × Cd × (A/m) × v² × v̂
#[must_use]
pub fn drag_acceleration(
    position: [f64; 3],
    velocity: [f64; 3],
    params: &AtmosphereParams,
) -> [f64; 3] {
    let r =
        (position[0] * position[0] + position[1] * position[1] + position[2] * position[2]).sqrt();
    let altitude = r - params.r_body;

    if altitude < 0.0 || params.h_scale <= 0.0 {
        return [0.0, 0.0, 0.0];
    }

    let rho = params.rho0 * (-((altitude - params.h0) / params.h_scale)).exp();

    let v2 = velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2];
    let v = v2.sqrt();

    if v < 1e-30 {
        return [0.0, 0.0, 0.0];
    }

    let factor = -0.5 * rho * params.cd * params.area_mass_ratio * v;

    [
        factor * velocity[0],
        factor * velocity[1],
        factor * velocity[2],
    ]
}

// ── Solar radiation pressure ──────────────────────────────────────────────

/// Compute solar radiation pressure (SRP) acceleration.
///
/// a_srp = −P × Cr × (A/m) × (AU/d)² × r̂_sun
///
/// # Arguments
///
/// * `position` — Satellite position (metres, ECI)
/// * `sun_position` — Sun position (metres, ECI)
/// * `cr` — Radiation pressure coefficient (1.0 = absorber, 2.0 = reflector)
/// * `area_mass_ratio` — Cross-sectional area / mass (m²/kg)
#[must_use]
pub fn srp_acceleration(
    position: [f64; 3],
    sun_position: [f64; 3],
    cr: f64,
    area_mass_ratio: f64,
) -> [f64; 3] {
    // Vector from satellite to sun
    let dx = sun_position[0] - position[0];
    let dy = sun_position[1] - position[1];
    let dz = sun_position[2] - position[2];
    let d2 = dx * dx + dy * dy + dz * dz;
    let d = d2.sqrt();

    if d < 1e-10 {
        return [0.0, 0.0, 0.0];
    }

    // Pressure scales as 1/d² relative to 1 AU
    let pressure = SOLAR_PRESSURE_1AU * (AU_METRES * AU_METRES) / d2;
    let factor = -pressure * cr * area_mass_ratio / d;

    // Acceleration pushes satellite away from the sun
    [factor * dx, factor * dy, factor * dz]
}

// ── Third-body perturbation ───────────────────────────────────────────────

/// Compute third-body gravitational perturbation acceleration.
///
/// Accounts for the differential gravity: the third body pulls on the satellite
/// differently than it pulls on the central body.
///
/// a = μ₃ × (r₃ₛ/|r₃ₛ|³ − r₃/|r₃|³)
///
/// # Arguments
///
/// * `position` — Satellite position relative to central body (metres)
/// * `third_body_position` — Third body position relative to central body (metres)
/// * `mu_third` — Gravitational parameter of the third body (m³/s²)
#[must_use]
pub fn third_body_acceleration(
    position: [f64; 3],
    third_body_position: [f64; 3],
    mu_third: f64,
) -> [f64; 3] {
    // Vector from satellite to third body
    let d = [
        third_body_position[0] - position[0],
        third_body_position[1] - position[1],
        third_body_position[2] - position[2],
    ];
    let d_mag = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();

    // Distance of third body from central body
    let r3_mag = (third_body_position[0] * third_body_position[0]
        + third_body_position[1] * third_body_position[1]
        + third_body_position[2] * third_body_position[2])
        .sqrt();

    if d_mag < 1e-10 || r3_mag < 1e-10 {
        return [0.0, 0.0, 0.0];
    }

    let d3 = d_mag * d_mag * d_mag;
    let r3_3 = r3_mag * r3_mag * r3_mag;

    [
        mu_third * (d[0] / d3 - third_body_position[0] / r3_3),
        mu_third * (d[1] / d3 - third_body_position[1] / r3_3),
        mu_third * (d[2] / d3 - third_body_position[2] / r3_3),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    const MU_EARTH: f64 = 3.986_004_418e14;

    // ── J2 ───────────────────────────────────────────────────────────

    #[test]
    fn j2_equatorial_symmetry() {
        // Symmetric about z=0: satellite at (r,0,0) should have zero z-acceleration
        let pos = [7e6, 0.0, 0.0];
        let acc = j2_acceleration(pos, MU_EARTH, J2_EARTH, R_EARTH);
        assert!(acc[2].abs() < 1e-10, "equatorial z-accel: {}", acc[2]);
        // x-acceleration should be negative (inward, deceleration)
        assert!(acc[0] < 0.0, "equatorial x-accel should be negative");
    }

    #[test]
    fn j2_polar_different() {
        // Satellite at same radius but on z-axis: different acceleration pattern
        let r = 7e6;
        let eq = j2_acceleration([r, 0.0, 0.0], MU_EARTH, J2_EARTH, R_EARTH);
        let pol = j2_acceleration([0.0, 0.0, r], MU_EARTH, J2_EARTH, R_EARTH);
        // Magnitudes should differ
        let mag_eq = (eq[0] * eq[0] + eq[1] * eq[1] + eq[2] * eq[2]).sqrt();
        let mag_pol = (pol[0] * pol[0] + pol[1] * pol[1] + pol[2] * pol[2]).sqrt();
        assert!(
            (mag_eq - mag_pol).abs() > 1e-10,
            "equatorial and polar J2 should differ"
        );
    }

    #[test]
    fn j2_magnitude_order() {
        // J2 perturbation at LEO should be ~1e-2 m/s² (much less than gravity ~8 m/s²)
        let pos = [R_EARTH + 400e3, 0.0, 0.0];
        let acc = j2_acceleration(pos, MU_EARTH, J2_EARTH, R_EARTH);
        let mag = (acc[0] * acc[0] + acc[1] * acc[1] + acc[2] * acc[2]).sqrt();
        assert!(mag > 1e-3 && mag < 1e-1, "J2 LEO magnitude: {} m/s²", mag);
    }

    #[test]
    fn j2_zero_position() {
        let acc = j2_acceleration([0.0, 0.0, 0.0], MU_EARTH, J2_EARTH, R_EARTH);
        assert_eq!(acc, [0.0, 0.0, 0.0]);
    }

    // ── J2 secular rates ─────────────────────────────────────────────

    #[test]
    fn j2_secular_sun_sync() {
        // Sun-synchronous orbit: i ≈ 98.6° at ~786 km altitude gives ~0.9856°/day
        // (exact values depend on precise altitude/inclination combination)
        let a = R_EARTH + 786e3;
        let i = 98.6_f64.to_radians();
        let (d_raan, _) = j2_secular_rates(a, 0.0, i, MU_EARTH, J2_EARTH, R_EARTH).unwrap();
        let deg_per_day = d_raan.to_degrees() * 86400.0;
        assert!(
            (deg_per_day - 0.9856).abs() < 0.1,
            "sun-sync RAAN rate: {deg_per_day}°/day"
        );
    }

    #[test]
    fn j2_secular_equatorial() {
        // Equatorial orbit (i=0): max RAAN regression
        let a = 7e6;
        let (d_raan, d_omega) = j2_secular_rates(a, 0.0, 0.0, MU_EARTH, J2_EARTH, R_EARTH).unwrap();
        assert!(d_raan < 0.0, "equatorial RAAN should regress");
        assert!(d_omega > 0.0, "equatorial ω should advance");
    }

    #[test]
    fn j2_secular_critical_inclination() {
        // At i = 63.4° (critical inclination), dω/dt = 0
        let a = 26_560e3; // Molniya SMA
        let i = 63.4_f64.to_radians();
        let (_, d_omega) = j2_secular_rates(a, 0.7, i, MU_EARTH, J2_EARTH, R_EARTH).unwrap();
        assert!(
            d_omega.abs() < 1e-10,
            "critical inclination dω/dt: {}",
            d_omega
        );
    }

    #[test]
    fn j2_secular_invalid() {
        assert!(j2_secular_rates(-1.0, 0.0, 0.0, MU_EARTH, J2_EARTH, R_EARTH).is_err());
        assert!(j2_secular_rates(7e6, 0.0, 0.0, -1.0, J2_EARTH, R_EARTH).is_err());
    }

    // ── J3 ───────────────────────────────────────────────────────────

    #[test]
    fn j3_much_smaller_than_j2() {
        let pos = [R_EARTH + 400e3, 0.0, R_EARTH + 400e3];
        let j2_acc = j2_acceleration(pos, MU_EARTH, J2_EARTH, R_EARTH);
        let j3_acc = j3_acceleration(pos, MU_EARTH, J3_EARTH, R_EARTH);
        let j2_mag = (j2_acc[0] * j2_acc[0] + j2_acc[1] * j2_acc[1] + j2_acc[2] * j2_acc[2]).sqrt();
        let j3_mag = (j3_acc[0] * j3_acc[0] + j3_acc[1] * j3_acc[1] + j3_acc[2] * j3_acc[2]).sqrt();
        assert!(
            j3_mag < j2_mag * 0.01,
            "J3 should be <<< J2: J3={j3_mag}, J2={j2_mag}"
        );
    }

    // ── Drag ─────────────────────────────────────────────────────────

    fn earth_atmo() -> AtmosphereParams {
        AtmosphereParams {
            cd: 2.2,
            area_mass_ratio: 0.01,
            rho0: 1.225,
            h0: 0.0,
            h_scale: 8500.0,
            r_body: R_EARTH,
        }
    }

    #[test]
    fn drag_opposes_velocity() {
        let pos = [R_EARTH + 400e3, 0.0, 0.0];
        let vel = [0.0, 7700.0, 0.0];
        let acc = drag_acceleration(pos, vel, &earth_atmo());
        assert!(acc[1] < 0.0, "drag should oppose velocity: {:?}", acc);
        assert!(acc[0].abs() < 1e-20, "no x-drag for y-velocity");
    }

    #[test]
    fn drag_decreases_with_altitude() {
        let vel = [0.0, 7700.0, 0.0];
        let atmo = earth_atmo();
        let low = drag_acceleration([R_EARTH + 200e3, 0.0, 0.0], vel, &atmo);
        let high = drag_acceleration([R_EARTH + 800e3, 0.0, 0.0], vel, &atmo);
        let mag_low = low[1].abs();
        let mag_high = high[1].abs();
        assert!(
            mag_low > mag_high * 10.0,
            "drag should decrease with altitude: low={mag_low}, high={mag_high}"
        );
    }

    #[test]
    fn drag_below_surface() {
        let acc = drag_acceleration([R_EARTH - 1.0, 0.0, 0.0], [0.0, 7700.0, 0.0], &earth_atmo());
        assert_eq!(acc, [0.0, 0.0, 0.0]);
    }

    // ── SRP ──────────────────────────────────────────────────────────

    #[test]
    fn srp_direction() {
        // Sun along +x, satellite at origin → acceleration should be along -x (pushed away)
        let sat = [0.0, 0.0, 0.0];
        let sun = [AU_METRES, 0.0, 0.0];
        let acc = srp_acceleration(sat, sun, 1.5, 0.01);
        assert!(acc[0] < 0.0, "SRP should push away from sun: {:?}", acc);
        assert!(acc[1].abs() < 1e-30);
    }

    #[test]
    fn srp_decreases_with_distance() {
        let sun = [AU_METRES, 0.0, 0.0];
        let near = srp_acceleration([0.0, 0.0, 0.0], sun, 1.5, 0.01);
        let far_sun = [2.0 * AU_METRES, 0.0, 0.0];
        let far = srp_acceleration([0.0, 0.0, 0.0], far_sun, 1.5, 0.01);
        assert!(
            near[0].abs() > far[0].abs() * 3.0,
            "SRP should decrease with distance"
        );
    }

    // ── Third body ───────────────────────────────────────────────────

    #[test]
    fn third_body_tidal() {
        // Moon at ~384400 km, satellite at LEO → small perturbation
        let sat = [R_EARTH + 400e3, 0.0, 0.0];
        let moon = [384_400e3, 0.0, 0.0];
        let mu_moon = 4.902_800e12;
        let acc = third_body_acceleration(sat, moon, mu_moon);
        let mag = (acc[0] * acc[0] + acc[1] * acc[1] + acc[2] * acc[2]).sqrt();
        // Lunar perturbation at LEO: ~1e-6 m/s²
        assert!(
            mag > 1e-8 && mag < 1e-4,
            "lunar perturbation magnitude: {mag}"
        );
    }

    #[test]
    fn third_body_zero() {
        let acc = third_body_acceleration([7e6, 0.0, 0.0], [0.0, 0.0, 0.0], 1e10);
        assert_eq!(acc, [0.0, 0.0, 0.0]);
    }
}
