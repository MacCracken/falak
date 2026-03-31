//! Orbital transfer maneuvers — Hohmann, bi-elliptic, and plane change.
//!
//! All functions take gravitational parameter μ (m³/s²) and orbit radii (m),
//! returning delta-v values in m/s and times in seconds.

use tracing::instrument;

use crate::error::{FalakError, Result};

/// Result of a Hohmann transfer between two circular orbits.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
#[non_exhaustive]
pub struct HohmannTransfer {
    /// Delta-v for the first burn (departure, m/s).
    pub delta_v1: f64,
    /// Delta-v for the second burn (arrival, m/s).
    pub delta_v2: f64,
    /// Total delta-v (m/s).
    pub total_delta_v: f64,
    /// Time of flight (seconds).
    pub time_of_flight: f64,
    /// Semi-major axis of the transfer ellipse (m).
    pub transfer_sma: f64,
}

/// Result of a bi-elliptic transfer.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
#[non_exhaustive]
pub struct BiEllipticTransfer {
    /// Delta-v for the first burn (m/s).
    pub delta_v1: f64,
    /// Delta-v at the intermediate radius (m/s).
    pub delta_v2: f64,
    /// Delta-v for the final circularisation (m/s).
    pub delta_v3: f64,
    /// Total delta-v (m/s).
    pub total_delta_v: f64,
    /// Total time of flight (seconds).
    pub time_of_flight: f64,
    /// Intermediate apoapsis radius (m).
    pub intermediate_radius: f64,
}

/// Result of a simple plane change maneuver.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
#[non_exhaustive]
pub struct PlaneChange {
    /// Delta-v required (m/s).
    pub delta_v: f64,
    /// Inclination change (radians).
    pub delta_inclination: f64,
}

// ── Hohmann transfer ──────────────────────────────────────────────────────

/// Compute a Hohmann transfer between two circular orbits.
///
/// The Hohmann transfer is the minimum-energy two-impulse transfer between
/// coplanar circular orbits. Optimal when r₂/r₁ < 11.94.
///
/// # Errors
///
/// Returns [`FalakError::InvalidParameter`] if radii or μ are not positive.
#[must_use = "returns the computed transfer"]
#[instrument(level = "trace")]
pub fn hohmann(r1: f64, r2: f64, mu: f64) -> Result<HohmannTransfer> {
    validate_positive("r1", r1)?;
    validate_positive("r2", r2)?;
    validate_positive("mu", mu)?;

    let transfer_sma = (r1 + r2) / 2.0;

    // Vis-viva at departure (on transfer orbit at r1)
    let v_circ1 = (mu / r1).sqrt();
    let v_transfer_dep = (mu * (2.0 / r1 - 1.0 / transfer_sma)).sqrt();
    let delta_v1 = (v_transfer_dep - v_circ1).abs();

    // Vis-viva at arrival (on transfer orbit at r2)
    let v_circ2 = (mu / r2).sqrt();
    let v_transfer_arr = (mu * (2.0 / r2 - 1.0 / transfer_sma)).sqrt();
    let delta_v2 = (v_circ2 - v_transfer_arr).abs();

    let total_delta_v = delta_v1 + delta_v2;

    // Time of flight = half the period of the transfer ellipse
    let time_of_flight =
        std::f64::consts::PI * (transfer_sma * transfer_sma * transfer_sma / mu).sqrt();

    Ok(HohmannTransfer {
        delta_v1,
        delta_v2,
        total_delta_v,
        time_of_flight,
        transfer_sma,
    })
}

// ── Bi-elliptic transfer ──────────────────────────────────────────────────

/// Compute a bi-elliptic transfer between two circular orbits.
///
/// Uses an intermediate radius `r_intermediate` > max(r1, r2).
/// More efficient than Hohmann when r₂/r₁ > 11.94.
///
/// # Errors
///
/// Returns [`FalakError::InvalidParameter`] if radii or μ are not positive,
/// or if `r_intermediate` is not larger than both `r1` and `r2`.
#[must_use = "returns the computed transfer"]
#[instrument(level = "trace")]
pub fn bi_elliptic(r1: f64, r2: f64, r_intermediate: f64, mu: f64) -> Result<BiEllipticTransfer> {
    validate_positive("r1", r1)?;
    validate_positive("r2", r2)?;
    validate_positive("r_intermediate", r_intermediate)?;
    validate_positive("mu", mu)?;

    let r_max = r1.max(r2);
    if r_intermediate <= r_max {
        return Err(FalakError::InvalidParameter(
            format!(
                "intermediate radius ({r_intermediate}) must be greater than max(r1, r2) = {r_max}"
            )
            .into(),
        ));
    }

    // First transfer ellipse: r1 → r_intermediate
    let a1 = (r1 + r_intermediate) / 2.0;
    let v_circ1 = (mu / r1).sqrt();
    let v_dep = (mu * (2.0 / r1 - 1.0 / a1)).sqrt();
    let delta_v1 = (v_dep - v_circ1).abs();

    // At r_intermediate: transition between the two transfer ellipses
    let a2 = (r2 + r_intermediate) / 2.0;
    let v_arr1 = (mu * (2.0 / r_intermediate - 1.0 / a1)).sqrt();
    let v_dep2 = (mu * (2.0 / r_intermediate - 1.0 / a2)).sqrt();
    let delta_v2 = (v_dep2 - v_arr1).abs();

    // At r2: circularise
    let v_arr2 = (mu * (2.0 / r2 - 1.0 / a2)).sqrt();
    let v_circ2 = (mu / r2).sqrt();
    let delta_v3 = (v_circ2 - v_arr2).abs();

    let total_delta_v = delta_v1 + delta_v2 + delta_v3;

    // Time = half-period of each transfer ellipse
    let tof1 = std::f64::consts::PI * (a1 * a1 * a1 / mu).sqrt();
    let tof2 = std::f64::consts::PI * (a2 * a2 * a2 / mu).sqrt();
    let time_of_flight = tof1 + tof2;

    Ok(BiEllipticTransfer {
        delta_v1,
        delta_v2,
        delta_v3,
        total_delta_v,
        time_of_flight,
        intermediate_radius: r_intermediate,
    })
}

// ── Plane change ──────────────────────────────────────────────────────────

/// Compute delta-v for a simple plane change at a given orbital velocity.
///
/// Δv = 2v sin(Δi / 2)
///
/// Most efficient when performed at the lowest velocity point (apoapsis
/// for elliptical orbits, or at higher radii).
///
/// # Errors
///
/// Returns [`FalakError::InvalidParameter`] if velocity is not positive
/// or inclination change is negative.
#[must_use = "returns the computed plane change"]
#[instrument(level = "trace")]
pub fn plane_change(velocity: f64, delta_inclination: f64) -> Result<PlaneChange> {
    validate_positive("velocity", velocity)?;
    if delta_inclination < 0.0 {
        return Err(FalakError::InvalidParameter(
            format!("inclination change must be non-negative, got {delta_inclination}").into(),
        ));
    }

    let delta_v = 2.0 * velocity * (delta_inclination / 2.0).sin();

    Ok(PlaneChange {
        delta_v,
        delta_inclination,
    })
}

// ── Phasing orbit ─────────────────────────────────────────────────────────

/// Compute delta-v for a phasing maneuver to adjust orbital phase.
///
/// Returns the delta-v for the initial burn (equal in magnitude to the
/// return burn) and the phasing orbit period.
///
/// `phase_angle` is the desired phase advance in radians (positive = ahead).
/// `n_orbits` is the number of phasing orbits (1 = single-revolution).
///
/// # Errors
///
/// Returns [`FalakError::InvalidParameter`] if radius, μ, or n_orbits are invalid.
#[must_use = "returns (delta_v, phasing_period)"]
#[instrument(level = "trace")]
pub fn phasing(radius: f64, mu: f64, phase_angle: f64, n_orbits: u32) -> Result<(f64, f64)> {
    validate_positive("radius", radius)?;
    validate_positive("mu", mu)?;
    if n_orbits == 0 {
        return Err(FalakError::InvalidParameter(
            "n_orbits must be at least 1".into(),
        ));
    }

    let t_orig = crate::kepler::orbital_period(radius, mu)?;
    // Phasing orbit period: complete n orbits while target advances by phase_angle
    let fraction = phase_angle / std::f64::consts::TAU;
    let t_phase = t_orig * (n_orbits as f64 - fraction) / n_orbits as f64;

    // Phasing orbit SMA from desired period
    // T = 2π√(a³/μ) → a = (μ(T/2π)²)^(1/3)
    let a_phase = (mu * (t_phase / std::f64::consts::TAU).powi(2)).cbrt();

    // Delta-v: difference between circular and phasing orbit velocity at radius
    let v_circ = (mu / radius).sqrt();
    let v_phase = (mu * (2.0 / radius - 1.0 / a_phase)).sqrt();
    let delta_v = (v_phase - v_circ).abs();

    Ok((delta_v, t_phase))
}

// ── Helpers ───────────────────────────────────────────────────────────────

fn validate_positive(name: &str, value: f64) -> Result<()> {
    if value <= 0.0 {
        return Err(FalakError::InvalidParameter(
            format!("{name} must be positive, got {value}").into(),
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    const MU_EARTH: f64 = 3.986_004_418e14;
    const R_LEO: f64 = 6.671e6; // ~300 km altitude
    const R_GEO: f64 = 42_164.0e3;

    // ── Hohmann ──────────────────────────────────────────────────────

    #[test]
    fn hohmann_leo_to_geo() {
        let h = hohmann(R_LEO, R_GEO, MU_EARTH).unwrap();
        // Standard textbook: ~3.9 km/s total
        assert!(
            (h.total_delta_v - 3893.0).abs() < 50.0,
            "LEO→GEO Hohmann Δv: {} m/s",
            h.total_delta_v
        );
        // Time of flight: ~5.26 hours
        let tof_hours = h.time_of_flight / 3600.0;
        assert!((tof_hours - 5.26).abs() < 0.2, "TOF: {tof_hours} hours");
        assert!(h.delta_v1 > 0.0);
        assert!(h.delta_v2 > 0.0);
        assert!(h.delta_v1 > h.delta_v2); // departure burn is larger
    }

    #[test]
    fn hohmann_symmetric() {
        // Transfer up then down should have same total delta-v
        let up = hohmann(R_LEO, R_GEO, MU_EARTH).unwrap();
        let down = hohmann(R_GEO, R_LEO, MU_EARTH).unwrap();
        assert!(
            (up.total_delta_v - down.total_delta_v).abs() < 1e-6,
            "up: {} down: {}",
            up.total_delta_v,
            down.total_delta_v
        );
    }

    #[test]
    fn hohmann_same_orbit() {
        let h = hohmann(R_LEO, R_LEO, MU_EARTH).unwrap();
        assert!(h.total_delta_v < 1e-6, "same orbit Δv: {}", h.total_delta_v);
    }

    #[test]
    fn hohmann_invalid() {
        assert!(hohmann(-1.0, R_GEO, MU_EARTH).is_err());
        assert!(hohmann(R_LEO, -1.0, MU_EARTH).is_err());
        assert!(hohmann(R_LEO, R_GEO, -1.0).is_err());
    }

    // ── Bi-elliptic ──────────────────────────────────────────────────

    #[test]
    fn bi_elliptic_basic() {
        let r_inter = R_GEO * 3.0;
        let b = bi_elliptic(R_LEO, R_GEO, r_inter, MU_EARTH).unwrap();
        assert!(b.total_delta_v > 0.0);
        assert!(b.time_of_flight > 0.0);
        assert_eq!(b.intermediate_radius, r_inter);
        // Bi-elliptic should be more expensive than Hohmann for low ratios
        let h = hohmann(R_LEO, R_GEO, MU_EARTH).unwrap();
        assert!(
            b.total_delta_v > h.total_delta_v,
            "bi-elliptic {} should cost more than Hohmann {} for r2/r1={:.1}",
            b.total_delta_v,
            h.total_delta_v,
            R_GEO / R_LEO
        );
    }

    #[test]
    fn bi_elliptic_high_ratio() {
        // For r2/r1 > 11.94, bi-elliptic beats Hohmann
        let r1 = 7e6;
        let r2 = r1 * 15.0; // ratio = 15
        let r_inter = r2 * 2.0;
        let b = bi_elliptic(r1, r2, r_inter, MU_EARTH).unwrap();
        let h = hohmann(r1, r2, MU_EARTH).unwrap();
        assert!(
            b.total_delta_v < h.total_delta_v,
            "bi-elliptic {} should beat Hohmann {} at ratio 15",
            b.total_delta_v,
            h.total_delta_v
        );
    }

    #[test]
    fn bi_elliptic_invalid() {
        // r_intermediate must exceed both radii
        assert!(bi_elliptic(R_LEO, R_GEO, R_LEO, MU_EARTH).is_err());
        assert!(bi_elliptic(R_LEO, R_GEO, -1.0, MU_EARTH).is_err());
    }

    // ── Plane change ─────────────────────────────────────────────────

    #[test]
    fn plane_change_basic() {
        let v = (MU_EARTH / R_LEO).sqrt(); // circular LEO speed
        let di = 28.5_f64.to_radians(); // Cape Canaveral latitude
        let pc = plane_change(v, di).unwrap();
        // Should be substantial — several km/s for LEO at 28.5°
        assert!(pc.delta_v > 3000.0, "plane change Δv: {}", pc.delta_v);
        assert!((pc.delta_inclination - di).abs() < 1e-15);
    }

    #[test]
    fn plane_change_zero() {
        let pc = plane_change(7700.0, 0.0).unwrap();
        assert!(pc.delta_v < 1e-10);
    }

    #[test]
    fn plane_change_invalid() {
        assert!(plane_change(-1.0, 0.5).is_err());
        assert!(plane_change(7700.0, -0.1).is_err());
    }

    // ── Phasing ──────────────────────────────────────────────────────

    #[test]
    fn phasing_basic() {
        let r = 7e6;
        let phase = std::f64::consts::FRAC_PI_4; // 45° phase advance
        let (dv, t_phase) = phasing(r, MU_EARTH, phase, 1).unwrap();
        assert!(dv > 0.0, "phasing Δv should be positive");
        let t_orig = crate::kepler::orbital_period(r, MU_EARTH).unwrap();
        assert!(
            t_phase < t_orig,
            "phasing period should be shorter for positive phase advance"
        );
    }

    #[test]
    fn phasing_invalid() {
        assert!(phasing(-1.0, MU_EARTH, 1.0, 1).is_err());
        assert!(phasing(7e6, MU_EARTH, 1.0, 0).is_err());
    }
}
