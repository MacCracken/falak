//! Orbital transfer maneuvers — Hohmann, bi-elliptic, plane change, and Lambert.
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

// ── Combined plane change + altitude maneuver ───────────────────────────

/// Result of a combined plane change and altitude change maneuver.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
#[non_exhaustive]
pub struct CombinedManeuver {
    /// Single-impulse delta-v (m/s).
    pub delta_v: f64,
    /// Delta-v that would be required for separate Hohmann + plane change (m/s).
    pub separate_delta_v: f64,
    /// Delta-v savings from combining the maneuvers (m/s).
    pub savings: f64,
    /// Initial circular velocity (m/s).
    pub v_initial: f64,
    /// Final circular velocity (m/s).
    pub v_final: f64,
}

/// Compute a combined plane change and altitude change maneuver.
///
/// Performs both the orbit raise/lower and inclination change in a single
/// impulse, which is more efficient than doing them separately (especially
/// when the plane change is large).
///
/// The optimal strategy performs the plane change at the highest altitude
/// (lowest velocity) to minimise the delta-v cost. This function computes
/// the single-impulse cost assuming the burn occurs at the initial altitude.
///
/// Δv = √(v₁² + v₂² − 2·v₁·v₂·cos(Δi))
///
/// # Arguments
///
/// * `r1` — Initial circular orbit radius (m)
/// * `r2` — Final circular orbit radius (m)
/// * `delta_inclination` — Inclination change (radians, non-negative)
/// * `mu` — Gravitational parameter (m³/s²)
///
/// # Errors
///
/// Returns [`FalakError::InvalidParameter`] if radii, μ, or inclination are invalid.
#[must_use = "returns the computed combined maneuver"]
#[instrument(level = "trace")]
pub fn combined_maneuver(
    r1: f64,
    r2: f64,
    delta_inclination: f64,
    mu: f64,
) -> Result<CombinedManeuver> {
    validate_positive("r1", r1)?;
    validate_positive("r2", r2)?;
    validate_positive("mu", mu)?;
    if delta_inclination < 0.0 {
        return Err(FalakError::InvalidParameter(
            format!("inclination change must be non-negative, got {delta_inclination}").into(),
        ));
    }

    let v1 = (mu / r1).sqrt();
    let v2 = (mu / r2).sqrt();

    // Combined single-impulse delta-v (law of cosines in velocity space)
    let combined_dv = (v1 * v1 + v2 * v2 - 2.0 * v1 * v2 * delta_inclination.cos()).sqrt();

    // Separate maneuvers for comparison:
    // Hohmann Δv + plane change Δv at the optimal point (apoapsis of transfer)
    let hohmann_dv = (v1 - (mu * (2.0 / r1 - 2.0 / (r1 + r2))).sqrt()).abs()
        + (v2 - (mu * (2.0 / r2 - 2.0 / (r1 + r2))).sqrt()).abs();
    let plane_dv = 2.0 * v2.min(v1) * (delta_inclination / 2.0).sin();
    let separate_dv = hohmann_dv + plane_dv;

    Ok(CombinedManeuver {
        delta_v: combined_dv,
        separate_delta_v: separate_dv,
        savings: separate_dv - combined_dv,
        v_initial: v1,
        v_final: v2,
    })
}

// ── Lambert problem ──────────────────────────────────────────────────────

/// Solution to Lambert's problem.
///
/// Given two position vectors and a time of flight, the Lambert solver finds
/// the velocity vectors at departure and arrival that connect them on a
/// Keplerian arc.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
#[non_exhaustive]
pub struct LambertSolution {
    /// Departure velocity `[vx, vy, vz]` (m/s, inertial frame).
    pub v1: [f64; 3],
    /// Arrival velocity `[vx, vy, vz]` (m/s, inertial frame).
    pub v2: [f64; 3],
}

/// Maximum iterations for Lambert solver.
const LAMBERT_MAX_ITER: u32 = 50;

/// Convergence tolerance for Lambert solver.
const LAMBERT_TOL: f64 = 1e-12;

/// Solve Lambert's problem: find the transfer orbit connecting two positions
/// in a given time of flight.
///
/// Uses the universal-variable method with Stumpff functions (Bate, Mueller &
/// White) and Newton iteration on the universal variable *z*. Handles
/// elliptic, parabolic, and hyperbolic transfers.
///
/// # Arguments
///
/// * `r1_vec` — Departure position `[x, y, z]` (metres, inertial frame)
/// * `r2_vec` — Arrival position `[x, y, z]` (metres, inertial frame)
/// * `tof` — Time of flight (seconds, must be positive)
/// * `mu` — Gravitational parameter (m³/s²)
/// * `prograde` — If `true`, use the short-way (prograde) arc; if `false`, use
///   the long-way (retrograde) arc.
///
/// # Errors
///
/// Returns [`FalakError::InvalidParameter`] if `tof` or `mu` are not positive,
/// or if positions are degenerate (zero-length).
///
/// Returns [`FalakError::ConvergenceError`] if the solver fails to converge.
#[must_use = "returns the Lambert solution"]
#[instrument(level = "debug")]
pub fn lambert(
    r1_vec: [f64; 3],
    r2_vec: [f64; 3],
    tof: f64,
    mu: f64,
    prograde: bool,
) -> Result<LambertSolution> {
    validate_positive("tof", tof)?;
    validate_positive("mu", mu)?;

    let r1 = (r1_vec[0] * r1_vec[0] + r1_vec[1] * r1_vec[1] + r1_vec[2] * r1_vec[2]).sqrt();
    let r2 = (r2_vec[0] * r2_vec[0] + r2_vec[1] * r2_vec[1] + r2_vec[2] * r2_vec[2]).sqrt();

    if r1 < 1e-10 || r2 < 1e-10 {
        return Err(FalakError::InvalidParameter(
            "position vectors must be non-zero".into(),
        ));
    }

    // Transfer angle
    let cos_dnu =
        (r1_vec[0] * r2_vec[0] + r1_vec[1] * r2_vec[1] + r1_vec[2] * r2_vec[2]) / (r1 * r2);
    let cos_dnu = cos_dnu.clamp(-1.0, 1.0);

    // Cross product z-component determines direction
    let cross_z = r1_vec[0] * r2_vec[1] - r1_vec[1] * r2_vec[0];

    let mut dnu = cos_dnu.acos();
    if prograde {
        if cross_z < 0.0 {
            dnu = std::f64::consts::TAU - dnu;
        }
    } else if cross_z >= 0.0 {
        dnu = std::f64::consts::TAU - dnu;
    }

    // A parameter (Bate, Mueller & White eq 5.35)
    let sin_dnu = dnu.sin();
    let a_param = sin_dnu * (r1 * r2 / (1.0 - cos_dnu)).sqrt();

    // Relative threshold: a_param should be significant relative to the orbit scale
    if a_param.abs() < 1e-8 * (r1 + r2) {
        return Err(FalakError::InvalidParameter(
            "degenerate transfer geometry (positions are collinear with 180° transfer)".into(),
        ));
    }

    // Newton iteration on universal variable z to match TOF
    // z > 0 → elliptic, z = 0 → parabolic, z < 0 → hyperbolic
    let mut z = lambert_initial_z(r1, r2, tof, mu, a_param);

    for iter in 0..LAMBERT_MAX_ITER {
        let (c2, c3) = stumpff(z);
        let y = r1 + r2 + a_param * (z * c3 - 1.0) / c2.sqrt();

        if y < 0.0 {
            // y must be positive; adjust z upward
            z = z.abs() + 0.1;
            continue;
        }

        let x = (y / c2).sqrt();
        let t_z = (x * x * x * c3 + a_param * y.sqrt()) / mu.sqrt();

        let dt = t_z - tof;
        if dt.abs() < LAMBERT_TOL * tof.max(1.0) {
            // Converged — compute velocities
            return Ok(lambert_solve_velocities(
                r1_vec, r2_vec, r1, r2, y, a_param, mu, tof,
            ));
        }

        // Derivative dT/dz via finite difference (robust for all z regimes)
        let h = z.abs() * 1e-7 + 1e-8;
        let (c2p, c3p) = stumpff(z + h);
        let yp = r1 + r2 + a_param * ((z + h) * c3p - 1.0) / c2p.sqrt();
        let dt_dz = if yp > 0.0 {
            let xp = (yp / c2p).sqrt();
            let t_zp = (xp * xp * xp * c3p + a_param * yp.sqrt()) / mu.sqrt();
            (t_zp - t_z) / h
        } else {
            // y went negative — step is too aggressive, use smaller h
            let h2 = h * 0.01;
            let (c2p2, c3p2) = stumpff(z + h2);
            let yp2 = r1 + r2 + a_param * ((z + h2) * c3p2 - 1.0) / c2p2.sqrt();
            if yp2 > 0.0 {
                let xp2 = (yp2 / c2p2).sqrt();
                let t_zp2 = (xp2 * xp2 * xp2 * c3p2 + a_param * yp2.sqrt()) / mu.sqrt();
                (t_zp2 - t_z) / h2
            } else {
                1.0 // fallback
            }
        };

        if dt_dz.abs() < 1e-30 {
            return Err(FalakError::ConvergenceError {
                message: "Lambert solver: zero derivative".into(),
                iterations: iter,
            });
        }

        z -= dt / dt_dz;
    }

    Err(FalakError::ConvergenceError {
        message: format!("Lambert solver did not converge (z={z})").into(),
        iterations: LAMBERT_MAX_ITER,
    })
}

/// Stumpff functions C₂(z) and C₃(z).
///
/// C₂(z) = (1 - cos√z)/z     for z > 0
/// C₃(z) = (√z - sin√z)/√z³  for z > 0
/// With Taylor series near z = 0 and hyperbolic equivalents for z < 0.
#[must_use]
#[inline]
fn stumpff(z: f64) -> (f64, f64) {
    if z > 1e-6 {
        let sz = z.sqrt();
        let c2 = (1.0 - sz.cos()) / z;
        let c3 = (sz - sz.sin()) / (z * sz);
        (c2, c3)
    } else if z < -1e-6 {
        let sz = (-z).sqrt();
        let c2 = (sz.cosh() - 1.0) / (-z);
        let c3 = (sz.sinh() - sz) / ((-z) * sz);
        (c2, c3)
    } else {
        // Taylor series near z = 0
        // C2 = 1/2 - z/24 + z²/720 ...
        // C3 = 1/6 - z/120 + z²/5040 ...
        let c2 = 1.0 / 2.0 - z / 24.0 + z * z / 720.0;
        let c3 = 1.0 / 6.0 - z / 120.0 + z * z / 5040.0;
        (c2, c3)
    }
}

/// Initial guess for z based on the transfer geometry.
fn lambert_initial_z(r1: f64, r2: f64, tof: f64, mu: f64, a_param: f64) -> f64 {
    // Compute TOF at z=0 (parabolic) to decide direction
    let (c2_0, c3_0) = stumpff(0.0);
    let y0 = r1 + r2 + a_param * (0.0 * c3_0 - 1.0) / c2_0.sqrt();
    if y0 > 0.0 {
        let x0 = (y0 / c2_0).sqrt();
        let t0 = (x0 * x0 * x0 * c3_0 + a_param * y0.sqrt()) / mu.sqrt();
        if t0 > tof {
            // TOF at z=0 is too long → need hyperbolic (z < 0)
            // Use bisection to bracket: find z where y > 0 and T < tof
            let mut z = -0.5;
            for _ in 0..20 {
                let (c2, c3) = stumpff(z);
                let y = r1 + r2 + a_param * (z * c3 - 1.0) / c2.sqrt();
                if y < 0.0 {
                    z *= 0.5; // too negative
                } else {
                    let xz = (y / c2).sqrt();
                    let tz = (xz * xz * xz * c3 + a_param * y.sqrt()) / mu.sqrt();
                    if tz < tof * 0.5 {
                        z *= 0.5; // overshot
                    } else {
                        break;
                    }
                }
                z *= 2.0;
            }
            z
        } else {
            // TOF at z=0 is too short → need elliptic (z > 0)
            // Estimate z from target period
            let period = std::f64::consts::TAU * ((r1 + r2) / 2.0).powi(3).sqrt() / mu.sqrt();
            let ratio = tof / period.max(1e-10);
            // z ≈ (2π)² for one orbit; scale by ratio
            (std::f64::consts::TAU * ratio)
                .powi(2)
                .min(4.0 * std::f64::consts::PI * std::f64::consts::PI)
        }
    } else {
        // y0 < 0 means a_param is negative (retrograde/long-way)
        // Need higher z to make y positive
        1.0
    }
}

/// Compute departure and arrival velocities from the converged Lambert solution.
#[allow(clippy::too_many_arguments)]
fn lambert_solve_velocities(
    r1_vec: [f64; 3],
    r2_vec: [f64; 3],
    r1: f64,
    r2: f64,
    y: f64,
    _a_param: f64,
    mu: f64,
    _tof: f64,
) -> LambertSolution {
    // Lagrange coefficients (BMW eq 5.46)
    let f = 1.0 - y / r1;
    let g_dot = 1.0 - y / r2;
    let g = _a_param * (y / mu).sqrt();

    if g.abs() < 1e-30 {
        return LambertSolution {
            v1: [0.0; 3],
            v2: [0.0; 3],
        };
    }

    let g_inv = 1.0 / g;
    let v1 = [
        (r2_vec[0] - f * r1_vec[0]) * g_inv,
        (r2_vec[1] - f * r1_vec[1]) * g_inv,
        (r2_vec[2] - f * r1_vec[2]) * g_inv,
    ];
    let v2 = [
        (g_dot * r2_vec[0] - r1_vec[0]) * g_inv,
        (g_dot * r2_vec[1] - r1_vec[1]) * g_inv,
        (g_dot * r2_vec[2] - r1_vec[2]) * g_inv,
    ];

    LambertSolution { v1, v2 }
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

    // ── Lambert ─────────────────────────────────────────────────────────

    #[test]
    fn lambert_recovers_known_orbit() {
        // Create a known orbit, propagate to get two states, then verify
        // Lambert recovers the velocities.
        let elem = crate::orbit::OrbitalElements::new(7e6, 0.1, 0.5, 1.0, 0.5, 0.0).unwrap();
        let state1 = crate::kepler::elements_to_state(&elem, MU_EARTH).unwrap();

        let dt = 1200.0; // 20 minutes
        let state2 = crate::propagate::kepler_to_state(&elem, MU_EARTH, dt).unwrap();

        let sol = lambert(state1.position, state2.position, dt, MU_EARTH, true).unwrap();

        // Departure velocity should match state1.velocity
        for i in 0..3 {
            let diff = (sol.v1[i] - state1.velocity[i]).abs();
            assert!(
                diff < 1.0, // within 1 m/s
                "v1[{i}]: Lambert={:.3} vs expected={:.3}, diff={diff:.3}",
                sol.v1[i],
                state1.velocity[i]
            );
        }

        // Arrival velocity should match state2.velocity
        for i in 0..3 {
            let diff = (sol.v2[i] - state2.velocity[i]).abs();
            assert!(
                diff < 1.0,
                "v2[{i}]: Lambert={:.3} vs expected={:.3}, diff={diff:.3}",
                sol.v2[i],
                state2.velocity[i]
            );
        }
    }

    #[test]
    fn lambert_circular_orbit() {
        // Circular equatorial orbit: 90° transfer
        let r = 7e6;
        let v_circ = (MU_EARTH / r).sqrt();
        let r1 = [r, 0.0, 0.0];
        let r2 = [0.0, r, 0.0]; // 90° ahead

        // TOF for 90° on a circular orbit
        let period = crate::kepler::orbital_period(r, MU_EARTH).unwrap();
        let tof = period / 4.0;

        let sol = lambert(r1, r2, tof, MU_EARTH, true).unwrap();

        // Departure velocity should be approximately [0, v_circ, 0]
        assert!(
            sol.v1[0].abs() < 100.0,
            "v1_x should be ~0, got {}",
            sol.v1[0]
        );
        assert!(
            (sol.v1[1] - v_circ).abs() < 100.0,
            "v1_y should be ~{v_circ}, got {}",
            sol.v1[1]
        );
    }

    #[test]
    fn lambert_near_hohmann() {
        // Verify Lambert gives approximately Hohmann-like delta-v for a
        // near-180° transfer (exact 180° is degenerate for Lambert).
        let r1_mag = R_LEO;
        let r2_mag = R_GEO;
        // Offset slightly from 180° to avoid collinear degeneracy
        let angle = 179.0_f64.to_radians();
        let r1 = [r1_mag, 0.0, 0.0];
        let r2 = [r2_mag * angle.cos(), r2_mag * angle.sin(), 0.0];

        let h = hohmann(r1_mag, r2_mag, MU_EARTH).unwrap();

        let sol = lambert(r1, r2, h.time_of_flight, MU_EARTH, true).unwrap();

        let v_dep = (sol.v1[0] * sol.v1[0] + sol.v1[1] * sol.v1[1] + sol.v1[2] * sol.v1[2]).sqrt();
        let v_circ = (MU_EARTH / r1_mag).sqrt();
        let lambert_dv1 = (v_dep - v_circ).abs();

        // Should be within ~100 m/s of Hohmann (not exact due to geometry offset)
        assert!(
            (lambert_dv1 - h.delta_v1).abs() < 100.0,
            "Lambert Δv1={lambert_dv1:.1} vs Hohmann Δv1={:.1}",
            h.delta_v1
        );
    }

    #[test]
    fn lambert_invalid_inputs() {
        let r1 = [7e6, 0.0, 0.0];
        let r2 = [0.0, 7e6, 0.0];
        assert!(lambert(r1, r2, -1.0, MU_EARTH, true).is_err());
        assert!(lambert(r1, r2, 1000.0, -1.0, true).is_err());
        assert!(lambert([0.0; 3], r2, 1000.0, MU_EARTH, true).is_err());
    }

    #[test]
    fn lambert_retrograde() {
        // Long-way transfer (retrograde)
        let r1 = [7e6, 0.0, 0.0];
        let r2 = [0.0, 7e6, 0.0];
        let period = crate::kepler::orbital_period(7e6, MU_EARTH).unwrap();
        let tof = period * 0.75; // 270° the long way

        let sol = lambert(r1, r2, tof, MU_EARTH, false).unwrap();

        // Should produce finite velocities
        for i in 0..3 {
            assert!(sol.v1[i].is_finite(), "v1[{i}] is not finite");
            assert!(sol.v2[i].is_finite(), "v2[{i}] is not finite");
        }

        // Speed should be orbital-magnitude
        let v1_mag = (sol.v1[0] * sol.v1[0] + sol.v1[1] * sol.v1[1] + sol.v1[2] * sol.v1[2]).sqrt();
        assert!(
            v1_mag > 1000.0 && v1_mag < 20000.0,
            "retrograde v1_mag={v1_mag} should be orbital"
        );
    }

    // ── Combined maneuver ─────────────────────────────────────────────

    #[test]
    fn combined_maneuver_pure_altitude() {
        // With zero inclination change, should equal Hohmann delta-v
        let cm = combined_maneuver(R_LEO, R_GEO, 0.0, MU_EARTH).unwrap();
        // Combined with di=0: dv = |v1 - v2| (single impulse)
        // This differs from Hohmann (two impulses), so combined < Hohmann is not guaranteed
        assert!(cm.delta_v > 0.0);
        assert!(
            (cm.savings).abs() < cm.separate_delta_v,
            "savings should be reasonable"
        );
    }

    #[test]
    fn combined_maneuver_pure_plane_change() {
        // Same altitude, only plane change: should match plane_change function
        let r = 7e6;
        let di = 28.5_f64.to_radians();
        let cm = combined_maneuver(r, r, di, MU_EARTH).unwrap();
        let pc = plane_change((MU_EARTH / r).sqrt(), di).unwrap();
        assert!(
            (cm.delta_v - pc.delta_v).abs() < 1.0,
            "combined={:.1} vs plane_change={:.1}",
            cm.delta_v,
            pc.delta_v
        );
    }

    #[test]
    fn combined_maneuver_saves_over_separate() {
        // Large plane change + altitude change: combined should save delta-v
        let cm = combined_maneuver(R_LEO, R_GEO, 28.5_f64.to_radians(), MU_EARTH).unwrap();
        assert!(
            cm.savings > 0.0,
            "combined should save over separate: savings={:.1}",
            cm.savings
        );
        assert!(
            cm.delta_v < cm.separate_delta_v,
            "combined {:.1} should < separate {:.1}",
            cm.delta_v,
            cm.separate_delta_v
        );
    }

    #[test]
    fn combined_maneuver_invalid() {
        assert!(combined_maneuver(-1.0, R_GEO, 0.5, MU_EARTH).is_err());
        assert!(combined_maneuver(R_LEO, R_GEO, -0.1, MU_EARTH).is_err());
    }

    // ── Lambert ─────────────────────────────────────────────────────────

    #[test]
    fn lambert_hyperbolic_short_tof() {
        // Very short TOF forces a hyperbolic transfer (high energy)
        let r1 = [7e6, 0.0, 0.0];
        let r2 = [0.0, 7e6, 0.0]; // 90°
        let tof = 500.0; // much shorter than quarter-period (~1400s)

        let sol = lambert(r1, r2, tof, MU_EARTH, true).unwrap();

        // Speed should be higher than circular (hyperbolic excess)
        let v_circ = (MU_EARTH / 7e6).sqrt();
        let v1_mag = (sol.v1[0] * sol.v1[0] + sol.v1[1] * sol.v1[1] + sol.v1[2] * sol.v1[2]).sqrt();
        assert!(
            v1_mag > v_circ * 1.1,
            "short TOF should give higher speed: v1={v1_mag:.0}, v_circ={v_circ:.0}"
        );
    }

    #[test]
    fn lambert_180_degrees_rejected() {
        // Exactly 180° (collinear, anti-parallel) is degenerate
        let r1 = [7e6, 0.0, 0.0];
        let r2 = [-7e6, 0.0, 0.0]; // exactly opposite
        let result = lambert(r1, r2, 3000.0, MU_EARTH, true);
        assert!(
            result.is_err(),
            "180° transfer should be rejected as degenerate"
        );
    }

    #[test]
    fn lambert_3d_inclined() {
        // Transfer with out-of-plane component
        let r1 = [7e6, 0.0, 0.0];
        let r2 = [0.0, 5e6, 5e6]; // inclined target
        let tof = 2000.0;

        let sol = lambert(r1, r2, tof, MU_EARTH, true).unwrap();

        // Velocity should have z-component (out of plane)
        assert!(
            sol.v1[2].abs() > 100.0,
            "3D transfer should have vz component: vz={}",
            sol.v1[2]
        );

        // Velocities should be finite
        for i in 0..3 {
            assert!(sol.v1[i].is_finite());
            assert!(sol.v2[i].is_finite());
        }
    }
}
