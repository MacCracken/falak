//! Orbit propagation — two-body (Kepler) and perturbed (Cowell's method).
//!
//! Propagates a spacecraft state forward in time using either analytic
//! two-body mechanics or numerical integration with perturbation forces.

use tracing::instrument;

use crate::error::{FalakError, Result};
use crate::kepler::StateVector;

// ── Two-body (Kepler) propagation ─────────────────────────────────────────

/// Propagate an orbit analytically using two-body mechanics (Kepler problem).
///
/// Advances the orbital elements by a time `dt` using mean anomaly propagation:
/// M(t) = M₀ + n × dt, then converts back to true anomaly.
///
/// This is exact for unperturbed orbits — no numerical error accumulates.
///
/// # Arguments
///
/// * `elements` — Initial orbital elements (elliptical only, e < 1).
/// * `mu` — Gravitational parameter (m³/s²).
/// * `dt` — Time to propagate (seconds, can be negative for backward propagation).
///
/// # Errors
///
/// Returns errors if elements are not elliptical or if Kepler's equation fails.
#[must_use = "returns the propagated orbital elements"]
#[instrument(level = "trace", skip(elements))]
pub fn kepler(
    elements: &crate::orbit::OrbitalElements,
    mu: f64,
    dt: f64,
) -> Result<crate::orbit::OrbitalElements> {
    if !elements.is_elliptical() {
        return Err(FalakError::InvalidParameter(
            format!(
                "two-body propagation requires elliptical orbit, got e={}",
                elements.eccentricity
            )
            .into(),
        ));
    }

    let n = crate::kepler::mean_motion(elements.semi_major_axis, mu)?;

    // Current true anomaly → mean anomaly
    let m0 = crate::kepler::true_to_mean_anomaly(elements.true_anomaly, elements.eccentricity);

    // Advance mean anomaly
    let m_new = (m0 + n * dt).rem_euclid(std::f64::consts::TAU);

    // Mean anomaly → true anomaly
    let nu_new = crate::kepler::mean_to_true_anomaly(m_new, elements.eccentricity)?;

    crate::orbit::OrbitalElements::new(
        elements.semi_major_axis,
        elements.eccentricity,
        elements.inclination,
        elements.raan,
        elements.argument_of_periapsis,
        nu_new,
    )
}

/// Propagate an orbit analytically and return the state vector at time `dt`.
///
/// Convenience wrapper: propagates elements, then converts to state vector.
#[must_use = "returns the propagated state vector"]
#[instrument(level = "trace", skip(elements))]
pub fn kepler_to_state(
    elements: &crate::orbit::OrbitalElements,
    mu: f64,
    dt: f64,
) -> Result<StateVector> {
    let new_elements = kepler(elements, mu, dt)?;
    crate::kepler::elements_to_state(&new_elements, mu)
}

// ── Cowell's method (perturbed propagation) ───────────────────────────────

/// Perturbation force function signature.
///
/// Takes `(position, velocity, time)` and returns acceleration `[ax, ay, az]`.
pub type PerturbationFn = dyn Fn([f64; 3], [f64; 3], f64) -> [f64; 3];

/// Propagate a state vector with perturbations using Cowell's method (RK4).
///
/// Directly integrates the equations of motion:
/// - dx/dt = v
/// - dv/dt = −μr/|r|³ + a_perturb(r, v, t)
///
/// # Arguments
///
/// * `state` — Initial state vector (position and velocity).
/// * `mu` — Central body gravitational parameter (m³/s²).
/// * `total_time` — Total propagation time (seconds).
/// * `dt` — Integration time step (seconds).
/// * `perturbation` — Optional perturbation acceleration function.
///
/// # Errors
///
/// Returns [`FalakError::InvalidParameter`] if `dt` is not positive or `mu` is not positive.
#[instrument(level = "debug", skip(state, perturbation))]
pub fn cowell(
    state: &StateVector,
    mu: f64,
    total_time: f64,
    dt: f64,
    perturbation: Option<&PerturbationFn>,
) -> Result<StateVector> {
    if mu <= 0.0 {
        return Err(FalakError::InvalidParameter(
            format!("mu must be positive, got {mu}").into(),
        ));
    }
    if dt <= 0.0 {
        return Err(FalakError::InvalidParameter(
            format!("time step must be positive, got {dt}").into(),
        ));
    }

    let mut pos = state.position;
    let mut vel = state.velocity;
    let mut t = 0.0;

    while t < total_time {
        let step_dt = dt.min(total_time - t);
        rk4_step(&mut pos, &mut vel, mu, t, step_dt, perturbation);
        t += step_dt;
    }

    Ok(StateVector {
        position: pos,
        velocity: vel,
    })
}

/// Single RK4 step for the two-body + perturbation ODE.
fn rk4_step(
    pos: &mut [f64; 3],
    vel: &mut [f64; 3],
    mu: f64,
    t: f64,
    dt: f64,
    perturbation: Option<&PerturbationFn>,
) {
    let accel = |p: [f64; 3], v: [f64; 3], time: f64| -> [f64; 3] {
        let r2 = p[0] * p[0] + p[1] * p[1] + p[2] * p[2];
        let r = r2.sqrt();
        let r3 = r2 * r;
        let mu_r3 = if r3 > 1e-30 { mu / r3 } else { 0.0 };

        let mut a = [-mu_r3 * p[0], -mu_r3 * p[1], -mu_r3 * p[2]];

        if let Some(perturb) = perturbation {
            let ap = perturb(p, v, time);
            a[0] += ap[0];
            a[1] += ap[1];
            a[2] += ap[2];
        }
        a
    };

    // k1
    let k1v = accel(*pos, *vel, t);
    let k1x = *vel;

    // k2
    let p2 = [
        pos[0] + 0.5 * dt * k1x[0],
        pos[1] + 0.5 * dt * k1x[1],
        pos[2] + 0.5 * dt * k1x[2],
    ];
    let v2 = [
        vel[0] + 0.5 * dt * k1v[0],
        vel[1] + 0.5 * dt * k1v[1],
        vel[2] + 0.5 * dt * k1v[2],
    ];
    let k2v = accel(p2, v2, t + 0.5 * dt);
    let k2x = v2;

    // k3
    let p3 = [
        pos[0] + 0.5 * dt * k2x[0],
        pos[1] + 0.5 * dt * k2x[1],
        pos[2] + 0.5 * dt * k2x[2],
    ];
    let v3 = [
        vel[0] + 0.5 * dt * k2v[0],
        vel[1] + 0.5 * dt * k2v[1],
        vel[2] + 0.5 * dt * k2v[2],
    ];
    let k3v = accel(p3, v3, t + 0.5 * dt);
    let k3x = v3;

    // k4
    let p4 = [
        pos[0] + dt * k3x[0],
        pos[1] + dt * k3x[1],
        pos[2] + dt * k3x[2],
    ];
    let v4 = [
        vel[0] + dt * k3v[0],
        vel[1] + dt * k3v[1],
        vel[2] + dt * k3v[2],
    ];
    let k4v = accel(p4, v4, t + dt);
    let k4x = v4;

    // Combine
    let dt6 = dt / 6.0;
    for i in 0..3 {
        pos[i] += dt6 * (k1x[i] + 2.0 * k2x[i] + 2.0 * k3x[i] + k4x[i]);
        vel[i] += dt6 * (k1v[i] + 2.0 * k2v[i] + 2.0 * k3v[i] + k4v[i]);
    }
}

/// Propagate a state vector with two-body gravity only (no perturbations).
///
/// Numerically integrates the two-body problem. For exact results without
/// accumulated error, prefer [`kepler`] with analytic propagation.
#[must_use = "returns the propagated state vector"]
#[instrument(level = "debug", skip(state))]
pub fn two_body(state: &StateVector, mu: f64, total_time: f64, dt: f64) -> Result<StateVector> {
    cowell(state, mu, total_time, dt, None)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::TAU;

    const MU_EARTH: f64 = 3.986_004_418e14;

    fn leo_elements() -> crate::orbit::OrbitalElements {
        crate::orbit::OrbitalElements::new(7e6, 0.01, 0.5, 1.0, 0.5, 0.0).unwrap()
    }

    fn leo_state() -> StateVector {
        crate::kepler::elements_to_state(&leo_elements(), MU_EARTH).unwrap()
    }

    // ── Kepler propagation ───────────────────────────────────────────

    #[test]
    fn kepler_full_period() {
        let elem = leo_elements();
        let period = crate::kepler::orbital_period(elem.semi_major_axis, MU_EARTH).unwrap();
        let prop = kepler(&elem, MU_EARTH, period).unwrap();

        // After one full period, true anomaly should return to start
        let diff = (prop.true_anomaly - elem.true_anomaly).rem_euclid(TAU);
        let diff = diff.min(TAU - diff);
        assert!(diff < 1e-8, "true anomaly drift: {diff}");

        // All other elements should be unchanged
        assert!((prop.semi_major_axis - elem.semi_major_axis).abs() < 1e-6);
        assert!((prop.eccentricity - elem.eccentricity).abs() < 1e-12);
    }

    #[test]
    fn kepler_half_period() {
        let elem = crate::orbit::OrbitalElements::new(7e6, 0.0, 0.0, 0.0, 0.0, 0.0).unwrap();
        let period = crate::kepler::orbital_period(7e6, MU_EARTH).unwrap();
        let prop = kepler(&elem, MU_EARTH, period / 2.0).unwrap();

        // Circular orbit: half period → ν ≈ π
        assert!(
            (prop.true_anomaly - std::f64::consts::PI).abs() < 1e-6,
            "half period ν: {}",
            prop.true_anomaly
        );
    }

    #[test]
    fn kepler_backward() {
        let elem = leo_elements();
        let period = crate::kepler::orbital_period(elem.semi_major_axis, MU_EARTH).unwrap();

        // Forward then backward should return to start
        let fwd = kepler(&elem, MU_EARTH, period * 0.3).unwrap();
        let back = kepler(&fwd, MU_EARTH, -period * 0.3).unwrap();

        let diff = (back.true_anomaly - elem.true_anomaly).rem_euclid(TAU);
        let diff = diff.min(TAU - diff);
        assert!(diff < 1e-8, "forward+backward drift: {diff}");
    }

    #[test]
    fn kepler_to_state_roundtrip() {
        let elem = leo_elements();
        let state = crate::kepler::elements_to_state(&elem, MU_EARTH).unwrap();
        let state2 = kepler_to_state(&elem, MU_EARTH, 0.0).unwrap();

        for i in 0..3 {
            assert!(
                (state.position[i] - state2.position[i]).abs() < 1e-4,
                "position[{i}]: {} vs {}",
                state.position[i],
                state2.position[i]
            );
        }
    }

    #[test]
    fn kepler_hyperbolic_rejected() {
        let hyp = crate::orbit::OrbitalElements::new(-1e7, 1.5, 0.0, 0.0, 0.0, 0.5).unwrap();
        assert!(kepler(&hyp, MU_EARTH, 100.0).is_err());
    }

    // ── Cowell propagation ───────────────────────────────────────────

    #[test]
    fn cowell_two_body_orbit() {
        let state = leo_state();
        let period = crate::kepler::orbital_period(7e6, MU_EARTH).unwrap();

        let result = cowell(&state, MU_EARTH, period, 10.0, None).unwrap();

        // Should return near start position
        let dr = [
            result.position[0] - state.position[0],
            result.position[1] - state.position[1],
            result.position[2] - state.position[2],
        ];
        let dr_mag = (dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]).sqrt();
        let frac = dr_mag / 7e6;
        assert!(frac < 0.001, "Cowell orbit closure: {frac:.6}");
    }

    #[test]
    fn cowell_with_j2() {
        let state = leo_state();

        // Propagate with J2 perturbation
        let j2_perturb = |pos: [f64; 3], _vel: [f64; 3], _t: f64| -> [f64; 3] {
            crate::perturbation::j2_acceleration(
                pos,
                MU_EARTH,
                crate::perturbation::J2_EARTH,
                crate::perturbation::R_EARTH,
            )
        };

        let result = cowell(&state, MU_EARTH, 1000.0, 10.0, Some(&j2_perturb)).unwrap();

        // Orbit should still be roughly the same radius (J2 doesn't change energy much)
        let r =
            (result.position[0].powi(2) + result.position[1].powi(2) + result.position[2].powi(2))
                .sqrt();
        assert!(
            (r - 7e6).abs() / 7e6 < 0.01,
            "J2 perturbed orbit radius: {r}"
        );
    }

    #[test]
    fn cowell_matches_kepler_unperturbed() {
        let elem = leo_elements();
        let state = crate::kepler::elements_to_state(&elem, MU_EARTH).unwrap();
        let dt_prop = 500.0;

        // Analytic Kepler propagation
        let kepler_state = kepler_to_state(&elem, MU_EARTH, dt_prop).unwrap();

        // Numerical Cowell propagation (no perturbations)
        let cowell_state = cowell(&state, MU_EARTH, dt_prop, 1.0, None).unwrap();

        // Should agree closely
        for i in 0..3 {
            let diff = (kepler_state.position[i] - cowell_state.position[i]).abs();
            assert!(
                diff < 100.0, // within 100m over 500s
                "position[{i}] diff: {diff}"
            );
        }
    }

    #[test]
    fn two_body_convenience() {
        let state = leo_state();
        let result = two_body(&state, MU_EARTH, 100.0, 10.0).unwrap();
        let r =
            (result.position[0].powi(2) + result.position[1].powi(2) + result.position[2].powi(2))
                .sqrt();
        // e=0.01 → radius varies by ±1% from SMA, so allow 2%
        assert!((r - 7e6).abs() / 7e6 < 0.02, "orbit radius: {r}");
    }

    #[test]
    fn cowell_invalid() {
        let state = leo_state();
        assert!(cowell(&state, -1.0, 100.0, 10.0, None).is_err());
        assert!(cowell(&state, MU_EARTH, 100.0, -1.0, None).is_err());
    }
}
