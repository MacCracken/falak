//! Kepler's laws — period, velocity, anomaly conversions (mean, eccentric, true).
//!
//! Covers elliptical (e < 1), parabolic (e = 1), and hyperbolic (e > 1) orbits.

use tracing::instrument;

use crate::error::{FalakError, Result};

/// Gravitational constant G (m³ kg⁻¹ s⁻²).
pub const G: f64 = 6.674_30e-11;

/// Default maximum iterations for Kepler's equation solver.
const MAX_ITERATIONS: u32 = 50;

/// Convergence tolerance for Kepler's equation solver.
const TOLERANCE: f64 = 1e-15;

// ── Orbital period and mean motion ────────────────────────────────────────

/// Orbital period for an elliptical orbit (seconds).
///
/// T = 2π √(a³ / μ)
///
/// # Errors
///
/// Returns [`FalakError::InvalidParameter`] if `semi_major_axis` or `mu` is not positive.
#[must_use = "returns the computed period"]
#[inline]
pub fn orbital_period(semi_major_axis: f64, mu: f64) -> Result<f64> {
    if semi_major_axis <= 0.0 {
        return Err(FalakError::InvalidParameter(
            format!("semi-major axis must be positive, got {semi_major_axis}").into(),
        ));
    }
    if mu <= 0.0 {
        return Err(FalakError::InvalidParameter(
            format!("gravitational parameter must be positive, got {mu}").into(),
        ));
    }
    let a = semi_major_axis;
    Ok(std::f64::consts::TAU * (a * a * a / mu).sqrt())
}

/// Mean motion for an elliptical orbit (rad/s).
///
/// n = √(μ / a³)
///
/// # Errors
///
/// Returns [`FalakError::InvalidParameter`] if `semi_major_axis` or `mu` is not positive.
#[must_use = "returns the computed mean motion"]
#[inline]
pub fn mean_motion(semi_major_axis: f64, mu: f64) -> Result<f64> {
    if semi_major_axis <= 0.0 {
        return Err(FalakError::InvalidParameter(
            format!("semi-major axis must be positive, got {semi_major_axis}").into(),
        ));
    }
    if mu <= 0.0 {
        return Err(FalakError::InvalidParameter(
            format!("gravitational parameter must be positive, got {mu}").into(),
        ));
    }
    let a = semi_major_axis;
    Ok((mu / (a * a * a)).sqrt())
}

// ── Vis-viva ──────────────────────────────────────────────────────────────

/// Orbital velocity at a given radius via the vis-viva equation (m/s).
///
/// v = √(μ (2/r − 1/a))
///
/// For hyperbolic orbits, `semi_major_axis` should be negative (convention).
///
/// # Errors
///
/// Returns [`FalakError::InvalidParameter`] if `radius` or `mu` is not positive.
#[must_use = "returns the computed velocity"]
#[inline]
pub fn vis_viva(radius: f64, semi_major_axis: f64, mu: f64) -> Result<f64> {
    if radius <= 0.0 {
        return Err(FalakError::InvalidParameter(
            format!("radius must be positive, got {radius}").into(),
        ));
    }
    if mu <= 0.0 {
        return Err(FalakError::InvalidParameter(
            format!("gravitational parameter must be positive, got {mu}").into(),
        ));
    }
    Ok((mu * (2.0 / radius - 1.0 / semi_major_axis)).sqrt())
}

// ── Kepler's equation (elliptical) ────────────────────────────────────────

/// Solve Kepler's equation for elliptical orbits: M = E − e sin(E).
///
/// Uses the Danby starter with Newton-Raphson iteration.
///
/// # Arguments
///
/// * `mean_anomaly` — Mean anomaly M (radians).
/// * `eccentricity` — Eccentricity e, must be in [0, 1).
///
/// # Errors
///
/// Returns [`FalakError::InvalidParameter`] if eccentricity is out of range.
/// Returns [`FalakError::ConvergenceError`] if the solver fails to converge.
#[must_use = "returns the computed eccentric anomaly"]
#[instrument(level = "trace")]
pub fn solve_kepler_elliptic(mean_anomaly: f64, eccentricity: f64) -> Result<f64> {
    if !(0.0..1.0).contains(&eccentricity) {
        return Err(FalakError::InvalidParameter(
            format!("eccentricity must be in [0, 1) for elliptical, got {eccentricity}").into(),
        ));
    }

    let e = eccentricity;
    // Normalise M to [0, 2π)
    let m = mean_anomaly.rem_euclid(std::f64::consts::TAU);

    // Danby starter: E₀ = M + e sin(M) + e² sin(2M)/2
    let mut ecc_anom = m + e * m.sin() + 0.5 * e * e * (2.0 * m).sin();

    for i in 0..MAX_ITERATIONS {
        let f = ecc_anom - e * ecc_anom.sin() - m;
        if f.abs() < TOLERANCE {
            return Ok(ecc_anom);
        }
        let f_prime = 1.0 - e * ecc_anom.cos();
        // Guard against near-zero derivative
        if f_prime.abs() < 1e-30 {
            return Err(FalakError::ConvergenceError {
                message: "Kepler elliptic: near-zero derivative".into(),
                iterations: i,
            });
        }
        ecc_anom -= f / f_prime;
    }

    Err(FalakError::ConvergenceError {
        message: format!("Kepler elliptic did not converge for M={mean_anomaly}, e={eccentricity}")
            .into(),
        iterations: MAX_ITERATIONS,
    })
}

/// Solve Kepler's equation for hyperbolic orbits: M = e sinh(H) − H.
///
/// Uses Newton-Raphson iteration.
///
/// # Arguments
///
/// * `mean_anomaly` — Hyperbolic mean anomaly M (radians).
/// * `eccentricity` — Eccentricity e, must be > 1.
///
/// # Errors
///
/// Returns [`FalakError::InvalidParameter`] if eccentricity is not > 1.
/// Returns [`FalakError::ConvergenceError`] if the solver fails to converge.
#[must_use = "returns the computed hyperbolic anomaly"]
#[instrument(level = "trace")]
pub fn solve_kepler_hyperbolic(mean_anomaly: f64, eccentricity: f64) -> Result<f64> {
    if eccentricity <= 1.0 {
        return Err(FalakError::InvalidParameter(
            format!("eccentricity must be > 1 for hyperbolic, got {eccentricity}").into(),
        ));
    }

    let e = eccentricity;
    let m = mean_anomaly;

    // Starter: H₀ = sign(M) × ln(2|M|/e + 1.8)
    let mut h = m.signum() * (2.0 * m.abs() / e + 1.8).ln();

    for i in 0..MAX_ITERATIONS {
        let f = e * h.sinh() - h - m;
        if f.abs() < TOLERANCE {
            return Ok(h);
        }
        let f_prime = e * h.cosh() - 1.0;
        if f_prime.abs() < 1e-30 {
            return Err(FalakError::ConvergenceError {
                message: "Kepler hyperbolic: near-zero derivative".into(),
                iterations: i,
            });
        }
        h -= f / f_prime;
    }

    Err(FalakError::ConvergenceError {
        message: format!(
            "Kepler hyperbolic did not converge for M={mean_anomaly}, e={eccentricity}"
        )
        .into(),
        iterations: MAX_ITERATIONS,
    })
}

// ── Anomaly conversions ───────────────────────────────────────────────────

/// Convert eccentric anomaly E to true anomaly ν (elliptical orbits).
///
/// tan(ν/2) = √((1+e)/(1−e)) × tan(E/2)
#[must_use]
#[inline]
pub fn eccentric_to_true_anomaly(eccentric_anomaly: f64, eccentricity: f64) -> f64 {
    let e = eccentricity;
    let half_e = eccentric_anomaly / 2.0;
    let factor = ((1.0 + e) / (1.0 - e)).sqrt();
    2.0 * (factor * half_e.tan()).atan()
}

/// Convert true anomaly ν to eccentric anomaly E (elliptical orbits).
///
/// tan(E/2) = √((1−e)/(1+e)) × tan(ν/2)
#[must_use]
#[inline]
pub fn true_to_eccentric_anomaly(true_anomaly: f64, eccentricity: f64) -> f64 {
    let e = eccentricity;
    let half_v = true_anomaly / 2.0;
    let factor = ((1.0 - e) / (1.0 + e)).sqrt();
    2.0 * (factor * half_v.tan()).atan()
}

/// Convert eccentric anomaly E to mean anomaly M (elliptical orbits).
///
/// M = E − e sin(E)
#[must_use]
#[inline]
pub fn eccentric_to_mean_anomaly(eccentric_anomaly: f64, eccentricity: f64) -> f64 {
    eccentric_anomaly - eccentricity * eccentric_anomaly.sin()
}

/// Convert mean anomaly M to true anomaly ν (elliptical orbits).
///
/// Solves Kepler's equation, then converts E → ν.
///
/// # Errors
///
/// Returns errors from [`solve_kepler_elliptic`].
#[must_use = "returns the computed true anomaly"]
pub fn mean_to_true_anomaly(mean_anomaly: f64, eccentricity: f64) -> Result<f64> {
    let e_anom = solve_kepler_elliptic(mean_anomaly, eccentricity)?;
    Ok(eccentric_to_true_anomaly(e_anom, eccentricity))
}

/// Convert true anomaly ν to mean anomaly M (elliptical orbits).
///
/// ν → E → M.
#[must_use]
#[inline]
pub fn true_to_mean_anomaly(true_anomaly: f64, eccentricity: f64) -> f64 {
    let e_anom = true_to_eccentric_anomaly(true_anomaly, eccentricity);
    eccentric_to_mean_anomaly(e_anom, eccentricity)
}

/// Convert hyperbolic anomaly H to true anomaly ν.
///
/// tan(ν/2) = √((e+1)/(e−1)) × tanh(H/2)
#[must_use]
#[inline]
pub fn hyperbolic_to_true_anomaly(hyperbolic_anomaly: f64, eccentricity: f64) -> f64 {
    let e = eccentricity;
    let factor = ((e + 1.0) / (e - 1.0)).sqrt();
    2.0 * (factor * (hyperbolic_anomaly / 2.0).tanh()).atan()
}

/// Convert true anomaly ν to hyperbolic anomaly H.
///
/// tanh(H/2) = √((e−1)/(e+1)) × tan(ν/2)
#[must_use]
#[inline]
pub fn true_to_hyperbolic_anomaly(true_anomaly: f64, eccentricity: f64) -> f64 {
    let e = eccentricity;
    let factor = ((e - 1.0) / (e + 1.0)).sqrt();
    2.0 * (factor * (true_anomaly / 2.0).tan()).atanh()
}

// ── Orbital radius ────────────────────────────────────────────────────────

/// Orbital radius at a given true anomaly.
///
/// r = a(1 − e²) / (1 + e cos(ν))
#[must_use]
#[inline]
pub fn radius_at_true_anomaly(semi_major_axis: f64, eccentricity: f64, true_anomaly: f64) -> f64 {
    let p = semi_major_axis * (1.0 - eccentricity * eccentricity);
    p / (1.0 + eccentricity * true_anomaly.cos())
}

// ── State vector <--> orbital elements ────────────────────────────────────

/// Cartesian state vector (position and velocity).
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
#[non_exhaustive]
pub struct StateVector {
    /// Position `[x, y, z]` (metres).
    pub position: [f64; 3],
    /// Velocity `[vx, vy, vz]` (m/s).
    pub velocity: [f64; 3],
}

/// Convert orbital elements to a Cartesian state vector.
///
/// Produces position and velocity in the inertial frame (ECI-like) from
/// classical Keplerian elements and the gravitational parameter μ.
///
/// # Errors
///
/// Returns [`FalakError::InvalidParameter`] if `mu` is not positive.
#[must_use = "returns the computed state vector"]
#[instrument(level = "trace", skip(elements))]
pub fn elements_to_state(elements: &crate::orbit::OrbitalElements, mu: f64) -> Result<StateVector> {
    if mu <= 0.0 {
        return Err(FalakError::InvalidParameter(
            format!("gravitational parameter must be positive, got {mu}").into(),
        ));
    }

    let a = elements.semi_major_axis;
    let e = elements.eccentricity;
    let i = elements.inclination;
    let omega_big = elements.raan;
    let omega = elements.argument_of_periapsis;
    let nu = elements.true_anomaly;

    // Semi-latus rectum
    let p = a * (1.0 - e * e);
    let r = p / (1.0 + e * nu.cos());

    // Position in perifocal frame
    let r_pf = [r * nu.cos(), r * nu.sin(), 0.0];

    // Velocity in perifocal frame
    let sqrt_mu_p = (mu / p).sqrt();
    let v_pf = [-sqrt_mu_p * nu.sin(), sqrt_mu_p * (e + nu.cos()), 0.0];

    // Rotation matrix: perifocal → inertial
    let cos_o = omega_big.cos();
    let sin_o = omega_big.sin();
    let cos_w = omega.cos();
    let sin_w = omega.sin();
    let cos_i = i.cos();
    let sin_i = i.sin();

    // Row 1 of rotation matrix
    let r11 = cos_o * cos_w - sin_o * sin_w * cos_i;
    let r12 = -(cos_o * sin_w + sin_o * cos_w * cos_i);
    // Row 2
    let r21 = sin_o * cos_w + cos_o * sin_w * cos_i;
    let r22 = -(sin_o * sin_w - cos_o * cos_w * cos_i);
    // Row 3
    let r31 = sin_w * sin_i;
    let r32 = cos_w * sin_i;

    let position = [
        r11 * r_pf[0] + r12 * r_pf[1],
        r21 * r_pf[0] + r22 * r_pf[1],
        r31 * r_pf[0] + r32 * r_pf[1],
    ];

    let velocity = [
        r11 * v_pf[0] + r12 * v_pf[1],
        r21 * v_pf[0] + r22 * v_pf[1],
        r31 * v_pf[0] + r32 * v_pf[1],
    ];

    Ok(StateVector { position, velocity })
}

/// Convert a Cartesian state vector to classical orbital elements.
///
/// # Errors
///
/// Returns [`FalakError::InvalidParameter`] if `mu` is not positive or position is zero.
#[must_use = "returns the computed orbital elements"]
#[instrument(level = "trace", skip(state))]
pub fn state_to_elements(state: &StateVector, mu: f64) -> Result<crate::orbit::OrbitalElements> {
    if mu <= 0.0 {
        return Err(FalakError::InvalidParameter(
            format!("gravitational parameter must be positive, got {mu}").into(),
        ));
    }

    let r = &state.position;
    let v = &state.velocity;

    let r_mag = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]).sqrt();
    if r_mag < 1e-30 {
        return Err(FalakError::InvalidParameter(
            "position vector has zero magnitude".into(),
        ));
    }
    let v_mag = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();

    // Angular momentum: h = r × v
    let h = [
        r[1] * v[2] - r[2] * v[1],
        r[2] * v[0] - r[0] * v[2],
        r[0] * v[1] - r[1] * v[0],
    ];
    let h_mag = (h[0] * h[0] + h[1] * h[1] + h[2] * h[2]).sqrt();

    // Node vector: n = k × h (k = [0,0,1])
    let n = [-h[1], h[0], 0.0];
    let n_mag = (n[0] * n[0] + n[1] * n[1]).sqrt();

    // Eccentricity vector: e_vec = (v × h)/μ − r̂
    let v_cross_h = [
        v[1] * h[2] - v[2] * h[1],
        v[2] * h[0] - v[0] * h[2],
        v[0] * h[1] - v[1] * h[0],
    ];
    let e_vec = [
        v_cross_h[0] / mu - r[0] / r_mag,
        v_cross_h[1] / mu - r[1] / r_mag,
        v_cross_h[2] / mu - r[2] / r_mag,
    ];
    let eccentricity = (e_vec[0] * e_vec[0] + e_vec[1] * e_vec[1] + e_vec[2] * e_vec[2]).sqrt();

    // Semi-major axis: a = -μ / (2ε)  where ε = v²/2 - μ/r
    let energy = 0.5 * v_mag * v_mag - mu / r_mag;
    let semi_major_axis = if energy.abs() < 1e-30 {
        // Parabolic — return semi-latus rectum as a
        h_mag * h_mag / mu
    } else {
        -mu / (2.0 * energy)
    };

    // Inclination
    let inclination = (h[2] / h_mag).clamp(-1.0, 1.0).acos();

    // RAAN
    let raan = if n_mag > 1e-30 {
        let val = (n[0] / n_mag).clamp(-1.0, 1.0).acos();
        if n[1] >= 0.0 {
            val
        } else {
            std::f64::consts::TAU - val
        }
    } else {
        0.0
    };

    // Argument of periapsis
    let argument_of_periapsis = if n_mag > 1e-30 && eccentricity > 1e-10 {
        let dot_n_e = n[0] * e_vec[0] + n[1] * e_vec[1] + n[2] * e_vec[2];
        let val = (dot_n_e / (n_mag * eccentricity)).clamp(-1.0, 1.0).acos();
        if e_vec[2] >= 0.0 {
            val
        } else {
            std::f64::consts::TAU - val
        }
    } else {
        0.0
    };

    // True anomaly
    let true_anomaly = if eccentricity > 1e-10 {
        let dot_e_r = e_vec[0] * r[0] + e_vec[1] * r[1] + e_vec[2] * r[2];
        let val = (dot_e_r / (eccentricity * r_mag)).clamp(-1.0, 1.0).acos();
        let r_dot_v = r[0] * v[0] + r[1] * v[1] + r[2] * v[2];
        if r_dot_v >= 0.0 {
            val
        } else {
            std::f64::consts::TAU - val
        }
    } else {
        // Circular orbit — use argument of latitude
        if n_mag > 1e-30 {
            let dot_n_r = n[0] * r[0] + n[1] * r[1] + n[2] * r[2];
            let val = (dot_n_r / (n_mag * r_mag)).clamp(-1.0, 1.0).acos();
            if r[2] >= 0.0 {
                val
            } else {
                std::f64::consts::TAU - val
            }
        } else {
            // Circular equatorial — use true longitude
            let val = (r[0] / r_mag).clamp(-1.0, 1.0).acos();
            if r[1] >= 0.0 {
                val
            } else {
                std::f64::consts::TAU - val
            }
        }
    };

    // Use internal constructor bypass for e >= 1 (hyperbolic/parabolic)
    // since OrbitalElements::new only accepts e < 1
    Ok(crate::orbit::OrbitalElements {
        semi_major_axis,
        eccentricity,
        inclination,
        raan,
        argument_of_periapsis,
        true_anomaly,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::{FRAC_PI_2, PI, TAU};

    const MU_EARTH: f64 = 3.986_004_418e14; // m³/s²

    // ── Period and mean motion ────────────────────────────────────────

    #[test]
    fn period_leo() {
        // ISS-like: a ≈ 6.771e6 m → T ≈ 5545s
        let t = orbital_period(6.771e6, MU_EARTH).unwrap();
        assert!((t - 5545.0).abs() < 5.0, "LEO period: {t}");
    }

    #[test]
    fn period_geo() {
        // GEO: a ≈ 42164 km → T ≈ 86164s (sidereal day)
        let t = orbital_period(42_164.0e3, MU_EARTH).unwrap();
        assert!((t - 86164.0).abs() < 10.0, "GEO period: {t}");
    }

    #[test]
    fn period_invalid() {
        assert!(orbital_period(-1.0, MU_EARTH).is_err());
        assert!(orbital_period(7e6, -1.0).is_err());
    }

    #[test]
    fn mean_motion_consistency() {
        let a = 7e6;
        let t = orbital_period(a, MU_EARTH).unwrap();
        let n = mean_motion(a, MU_EARTH).unwrap();
        assert!((n - TAU / t).abs() < 1e-12);
    }

    // ── Vis-viva ─────────────────────────────────────────────────────

    #[test]
    fn vis_viva_circular() {
        let r = 7e6;
        let v = vis_viva(r, r, MU_EARTH).unwrap();
        let expected = (MU_EARTH / r).sqrt();
        assert!(
            (v - expected).abs() < 0.01,
            "circular speed: {v} vs {expected}"
        );
    }

    #[test]
    fn vis_viva_invalid() {
        assert!(vis_viva(-1.0, 7e6, MU_EARTH).is_err());
        assert!(vis_viva(7e6, 7e6, -1.0).is_err());
    }

    // ── Kepler's equation (elliptic) ─────────────────────────────────

    #[test]
    fn kepler_circular() {
        // e = 0 → E = M
        let e_anom = solve_kepler_elliptic(1.0, 0.0).unwrap();
        assert!((e_anom - 1.0).abs() < TOLERANCE);
    }

    #[test]
    fn kepler_moderate_ecc() {
        let m = 1.5;
        let e = 0.5;
        let ea = solve_kepler_elliptic(m, e).unwrap();
        // Verify: M = E - e sin(E)
        let m_check = ea - e * ea.sin();
        assert!((m_check - m).abs() < 1e-12, "residual: {}", m_check - m);
    }

    #[test]
    fn kepler_high_ecc() {
        let m = 0.1;
        let e = 0.99;
        let ea = solve_kepler_elliptic(m, e).unwrap();
        let m_check = ea - e * ea.sin();
        assert!((m_check - m).abs() < 1e-12);
    }

    #[test]
    fn kepler_negative_m() {
        // Should handle negative mean anomaly
        let ea = solve_kepler_elliptic(-1.0, 0.3).unwrap();
        let m_check = ea - 0.3 * ea.sin();
        // rem_euclid normalises to [0, 2π)
        let m_norm = (-1.0_f64).rem_euclid(TAU);
        assert!((m_check - m_norm).abs() < 1e-12);
    }

    #[test]
    fn kepler_invalid_ecc() {
        assert!(solve_kepler_elliptic(1.0, -0.1).is_err());
        assert!(solve_kepler_elliptic(1.0, 1.0).is_err());
    }

    // ── Kepler's equation (hyperbolic) ───────────────────────────────

    #[test]
    fn kepler_hyperbolic_basic() {
        let m = 2.0;
        let e = 1.5;
        let h = solve_kepler_hyperbolic(m, e).unwrap();
        let m_check = e * h.sinh() - h;
        assert!((m_check - m).abs() < 1e-12, "residual: {}", m_check - m);
    }

    #[test]
    fn kepler_hyperbolic_negative() {
        let m = -3.0;
        let e = 2.0;
        let h = solve_kepler_hyperbolic(m, e).unwrap();
        let m_check = e * h.sinh() - h;
        assert!((m_check - m).abs() < 1e-12);
    }

    #[test]
    fn kepler_hyperbolic_invalid() {
        assert!(solve_kepler_hyperbolic(1.0, 0.5).is_err());
        assert!(solve_kepler_hyperbolic(1.0, 1.0).is_err());
    }

    // ── Anomaly conversions ──────────────────────────────────────────

    #[test]
    fn eccentric_true_roundtrip() {
        let e = 0.3;
        for deg in (0..360).step_by(15) {
            let nu = deg as f64 * PI / 180.0;
            let ea = true_to_eccentric_anomaly(nu, e);
            let nu2 = eccentric_to_true_anomaly(ea, e);
            // Normalise to [0, 2π) for comparison
            let diff = (nu - nu2).rem_euclid(TAU);
            let diff = diff.min(TAU - diff);
            assert!(diff < 1e-12, "roundtrip failed at {deg}°: diff={diff}");
        }
    }

    #[test]
    fn mean_true_roundtrip() {
        let e = 0.4;
        for deg in (0..360).step_by(30) {
            let nu = deg as f64 * PI / 180.0;
            let m = true_to_mean_anomaly(nu, e);
            let nu2 = mean_to_true_anomaly(m, e).unwrap();
            let diff = (nu - nu2).rem_euclid(TAU);
            let diff = diff.min(TAU - diff);
            assert!(diff < 1e-10, "roundtrip failed at {deg}°: diff={diff}");
        }
    }

    #[test]
    fn eccentric_to_mean_identity_circular() {
        // e = 0 → M = E
        let m = eccentric_to_mean_anomaly(2.0, 0.0);
        assert!((m - 2.0).abs() < 1e-15);
    }

    #[test]
    fn hyperbolic_true_roundtrip() {
        let e = 2.0;
        for &nu in &[0.1, 0.5, 1.0, -0.5] {
            let h = true_to_hyperbolic_anomaly(nu, e);
            let nu2 = hyperbolic_to_true_anomaly(h, e);
            assert!((nu - nu2).abs() < 1e-12, "roundtrip failed for ν={nu}");
        }
    }

    // ── Radius ───────────────────────────────────────────────────────

    #[test]
    fn radius_periapsis() {
        let r = radius_at_true_anomaly(10000.0, 0.5, 0.0);
        assert!((r - 5000.0).abs() < 1e-10);
    }

    #[test]
    fn radius_apoapsis() {
        let r = radius_at_true_anomaly(10000.0, 0.5, PI);
        assert!((r - 15000.0).abs() < 1e-10, "apoapsis: {r}");
    }

    // ── State vector conversions ─────────────────────────────────────

    #[test]
    fn elements_state_roundtrip() {
        let original = crate::orbit::OrbitalElements::new(7e6, 0.1, 0.5, 1.0, 0.5, 0.8).unwrap();

        let state = elements_to_state(&original, MU_EARTH).unwrap();
        let recovered = state_to_elements(&state, MU_EARTH).unwrap();

        assert!(
            (original.semi_major_axis - recovered.semi_major_axis).abs() < 1.0,
            "a: {} vs {}",
            original.semi_major_axis,
            recovered.semi_major_axis
        );
        assert!(
            (original.eccentricity - recovered.eccentricity).abs() < 1e-10,
            "e: {} vs {}",
            original.eccentricity,
            recovered.eccentricity
        );
        assert!(
            (original.inclination - recovered.inclination).abs() < 1e-10,
            "i: {} vs {}",
            original.inclination,
            recovered.inclination
        );
        assert!(
            (original.raan - recovered.raan).abs() < 1e-10,
            "Ω: {} vs {}",
            original.raan,
            recovered.raan
        );
        assert!(
            (original.argument_of_periapsis - recovered.argument_of_periapsis).abs() < 1e-10,
            "ω: {} vs {}",
            original.argument_of_periapsis,
            recovered.argument_of_periapsis
        );
        let diff = (original.true_anomaly - recovered.true_anomaly).rem_euclid(TAU);
        let diff = diff.min(TAU - diff);
        assert!(
            diff < 1e-10,
            "ν: {} vs {}",
            original.true_anomaly,
            recovered.true_anomaly
        );
    }

    #[test]
    fn elements_state_circular_equatorial() {
        let elements = crate::orbit::OrbitalElements::new(7e6, 0.0, 0.0, 0.0, 0.0, 0.0).unwrap();

        let state = elements_to_state(&elements, MU_EARTH).unwrap();
        // Position should be along +x
        assert!((state.position[0] - 7e6).abs() < 1.0);
        assert!(state.position[1].abs() < 1e-6);
        assert!(state.position[2].abs() < 1e-6);
        // Velocity should be along +y
        let v_circ = (MU_EARTH / 7e6).sqrt();
        assert!(state.velocity[0].abs() < 1e-6);
        assert!((state.velocity[1] - v_circ).abs() < 0.01);
    }

    #[test]
    fn elements_state_polar() {
        let elements =
            crate::orbit::OrbitalElements::new(7e6, 0.0, FRAC_PI_2, 0.0, 0.0, 0.0).unwrap();

        let state = elements_to_state(&elements, MU_EARTH).unwrap();
        // Polar orbit at ν=0 → position along +x, velocity has z component
        assert!((state.position[0] - 7e6).abs() < 1.0);
        assert!(
            state.velocity[2].abs() > 1000.0,
            "polar orbit should have z-velocity"
        );
    }

    #[test]
    fn state_to_elements_invalid() {
        let state = StateVector {
            position: [7e6, 0.0, 0.0],
            velocity: [0.0, 7500.0, 0.0],
        };
        assert!(state_to_elements(&state, -1.0).is_err());
    }
}
