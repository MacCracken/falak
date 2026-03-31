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

/// Perturbation acceleration function signature.
///
/// # Parameters
///
/// * `position` — Satellite position `[x, y, z]` in ECI (metres)
/// * `velocity` — Satellite velocity `[vx, vy, vz]` in ECI (m/s)
/// * `time` — Elapsed time since propagation start (seconds)
///
/// # Returns
///
/// Perturbation acceleration `[ax, ay, az]` in m/s². This is added to the
/// two-body gravitational acceleration during numerical integration.
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

// ── Encke's method ────────────────────────────────────────────────────────

/// Propagate using Encke's method: integrate the deviation from a reference
/// two-body orbit.
///
/// More accurate than Cowell for nearly-Keplerian orbits because the integrator
/// handles only the small perturbation deviation, not the full gravity field.
///
/// # Arguments
///
/// * `elements` — Initial orbital elements (elliptical only).
/// * `mu` — Central body gravitational parameter (m³/s²).
/// * `total_time` — Propagation time (seconds).
/// * `dt` — Integration time step (seconds).
/// * `perturbation` — Perturbation acceleration function.
///
/// # Errors
///
/// Returns errors if elements are not elliptical or parameters are invalid.
#[must_use = "returns the propagated state vector"]
#[instrument(level = "debug", skip(elements, perturbation))]
pub fn encke(
    elements: &crate::orbit::OrbitalElements,
    mu: f64,
    total_time: f64,
    dt: f64,
    perturbation: &PerturbationFn,
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
    if !elements.is_elliptical() {
        return Err(FalakError::InvalidParameter(
            "Encke's method requires elliptical reference orbit".into(),
        ));
    }

    // Deviation from reference: δr, δv
    let mut dr = [0.0; 3];
    let mut dv = [0.0; 3];
    let mut t = 0.0;

    while t < total_time {
        let step_dt = dt.min(total_time - t);

        // Reference state at current time
        let ref_state = kepler_to_state(elements, mu, t)?;

        // Actual position = reference + deviation
        let pos = [
            ref_state.position[0] + dr[0],
            ref_state.position[1] + dr[1],
            ref_state.position[2] + dr[2],
        ];
        let vel = [
            ref_state.velocity[0] + dv[0],
            ref_state.velocity[1] + dv[1],
            ref_state.velocity[2] + dv[2],
        ];

        // Encke acceleration: a_perturb + μ(F(q)·r_ref/r_ref³ - δr/r³)
        // where F(q) corrects for the difference in inverse-cube terms
        let r_ref_sq = ref_state.position[0].powi(2)
            + ref_state.position[1].powi(2)
            + ref_state.position[2].powi(2);
        let r_sq = pos[0].powi(2) + pos[1].powi(2) + pos[2].powi(2);

        // Battin's F(q) function for numerical stability
        let q = (dr[0] * (dr[0] - 2.0 * pos[0])
            + dr[1] * (dr[1] - 2.0 * pos[1])
            + dr[2] * (dr[2] - 2.0 * pos[2]))
            / r_sq;
        let fq = q * (3.0 + 3.0 * q + q * q) / (1.0 + (1.0 + q).sqrt().powi(3));

        let r_ref_3 = r_ref_sq * r_ref_sq.sqrt();
        let a_perturb = perturbation(pos, vel, t);

        // Deviation acceleration
        let mut da = [0.0; 3];
        for i in 0..3 {
            da[i] = a_perturb[i] - mu / r_ref_3 * (dr[i] + fq * ref_state.position[i]);
        }

        // RK4 on the deviation with proper time offsets for each stage
        let encke_accel_at = |dr_t: [f64; 3], dv_t: [f64; 3], time_offset: f64| -> [f64; 3] {
            let ref_t = kepler_to_state(elements, mu, t + time_offset).unwrap_or(ref_state.clone());
            let pos_t = [
                ref_t.position[0] + dr_t[0],
                ref_t.position[1] + dr_t[1],
                ref_t.position[2] + dr_t[2],
            ];
            let vel_t = [
                ref_t.velocity[0] + dv_t[0],
                ref_t.velocity[1] + dv_t[1],
                ref_t.velocity[2] + dv_t[2],
            ];
            let r_sq_t = pos_t[0].powi(2) + pos_t[1].powi(2) + pos_t[2].powi(2);
            let q_t = (dr_t[0] * (dr_t[0] - 2.0 * pos_t[0])
                + dr_t[1] * (dr_t[1] - 2.0 * pos_t[1])
                + dr_t[2] * (dr_t[2] - 2.0 * pos_t[2]))
                / r_sq_t;
            let fq_t = q_t * (3.0 + 3.0 * q_t + q_t * q_t) / (1.0 + (1.0 + q_t).sqrt().powi(3));
            let r_ref_sq_t =
                ref_t.position[0].powi(2) + ref_t.position[1].powi(2) + ref_t.position[2].powi(2);
            let r_ref_3_t = r_ref_sq_t * r_ref_sq_t.sqrt();
            let ap = perturbation(pos_t, vel_t, t + time_offset);
            let mut a = [0.0; 3];
            for i in 0..3 {
                a[i] = ap[i] - mu / r_ref_3_t * (dr_t[i] + fq_t * ref_t.position[i]);
            }
            a
        };

        // RK4 stages with correct time offsets
        let k1v = da;
        let k1x = dv;

        let half_dt = step_dt * 0.5;
        let dr2 = [
            dr[0] + half_dt * k1x[0],
            dr[1] + half_dt * k1x[1],
            dr[2] + half_dt * k1x[2],
        ];
        let dv2 = [
            dv[0] + half_dt * k1v[0],
            dv[1] + half_dt * k1v[1],
            dv[2] + half_dt * k1v[2],
        ];
        let k2v = encke_accel_at(dr2, dv2, half_dt);
        let k2x = dv2;

        let dr3 = [
            dr[0] + half_dt * k2x[0],
            dr[1] + half_dt * k2x[1],
            dr[2] + half_dt * k2x[2],
        ];
        let dv3 = [
            dv[0] + half_dt * k2v[0],
            dv[1] + half_dt * k2v[1],
            dv[2] + half_dt * k2v[2],
        ];
        let k3v = encke_accel_at(dr3, dv3, half_dt);
        let k3x = dv3;

        let dr4 = [
            dr[0] + step_dt * k3x[0],
            dr[1] + step_dt * k3x[1],
            dr[2] + step_dt * k3x[2],
        ];
        let dv4 = [
            dv[0] + step_dt * k3v[0],
            dv[1] + step_dt * k3v[1],
            dv[2] + step_dt * k3v[2],
        ];
        let k4v = encke_accel_at(dr4, dv4, step_dt); // k4 at full step
        let k4x = dv4;

        let dt6 = step_dt / 6.0;
        for i in 0..3 {
            dr[i] += dt6 * (k1x[i] + 2.0 * k2x[i] + 2.0 * k3x[i] + k4x[i]);
            dv[i] += dt6 * (k1v[i] + 2.0 * k2v[i] + 2.0 * k3v[i] + k4v[i]);
        }

        t += step_dt;
    }

    // Final state = reference at total_time + deviation
    let ref_final = kepler_to_state(elements, mu, total_time)?;
    Ok(StateVector {
        position: [
            ref_final.position[0] + dr[0],
            ref_final.position[1] + dr[1],
            ref_final.position[2] + dr[2],
        ],
        velocity: [
            ref_final.velocity[0] + dv[0],
            ref_final.velocity[1] + dv[1],
            ref_final.velocity[2] + dv[2],
        ],
    })
}

// ── Mean elements (Brouwer J2 theory) ────────────────────────────────────

/// Mean orbital elements with J2 short-period variations removed.
///
/// Mean elements change only due to secular and long-period effects,
/// making them suitable for long-term orbit prediction without numerical
/// integration.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
#[non_exhaustive]
pub struct MeanElements {
    /// Mean semi-major axis (m).
    pub semi_major_axis: f64,
    /// Mean eccentricity.
    pub eccentricity: f64,
    /// Mean inclination (radians).
    pub inclination: f64,
    /// Mean RAAN (radians).
    pub raan: f64,
    /// Mean argument of periapsis (radians).
    pub argument_of_periapsis: f64,
    /// Mean anomaly (radians).
    pub mean_anomaly: f64,
}

/// Convert osculating orbital elements to mean elements by removing
/// J2 short-period variations (first-order Brouwer theory).
///
/// The short-period terms oscillate at the orbital period and average to
/// zero over one revolution. Removing them gives the slowly-varying mean
/// elements.
///
/// # Arguments
///
/// * `osc` — Osculating orbital elements
/// * `mu` — Gravitational parameter (m³/s²)
/// * `j2` — J2 zonal harmonic coefficient
/// * `r_body` — Central body equatorial radius (m)
///
/// # Errors
///
/// Returns an error if the orbit is not elliptical or parameters are invalid.
#[must_use = "returns the mean elements"]
#[instrument(level = "trace", skip(osc))]
pub fn osculating_to_mean(
    osc: &crate::orbit::OrbitalElements,
    mu: f64,
    j2: f64,
    r_body: f64,
) -> Result<MeanElements> {
    if !osc.is_elliptical() {
        return Err(FalakError::InvalidParameter(
            "osculating-to-mean conversion requires elliptical orbit".into(),
        ));
    }
    if mu <= 0.0 || r_body <= 0.0 {
        return Err(FalakError::InvalidParameter(
            "mu and r_body must be positive".into(),
        ));
    }

    let a = osc.semi_major_axis;
    let e = osc.eccentricity;
    let i = osc.inclination;
    let omega = osc.argument_of_periapsis;
    let nu = osc.true_anomaly;

    let p = a * (1.0 - e * e);
    let eta = (1.0 - e * e).sqrt();
    let gamma2 = j2 * (r_body / p) * (r_body / p); // J2 × (R/p)²

    let sin_i = i.sin();
    let cos_i = i.cos();
    let sin2_i = sin_i * sin_i;

    // Eccentric anomaly from true anomaly
    let ea = 2.0 * ((nu / 2.0).tan() / ((1.0 + e) / (1.0 - e)).sqrt()).atan();
    let m_osc = ea - e * ea.sin(); // osculating mean anomaly

    // J2 short-period corrections (first-order Brouwer)
    let r = p / (1.0 + e * nu.cos());
    let a_r = a / r;

    // Short-period variation in semi-major axis
    let da_sp = a * gamma2 * ((a_r).powi(3) - 1.0 / eta.powi(3)) * (1.0 - 1.5 * sin2_i)
        + a * gamma2 * 1.5 * sin2_i * (a_r).powi(3) * (2.0 * (omega + nu)).cos();

    // Short-period variation in eccentricity
    let de_sp = gamma2 / (4.0 * e)
        * eta
        * ((1.0 - 1.5 * sin2_i) * (1.0 - (a_r).powi(2) + 2.0 * e * (a_r - 1.0 / (1.0 + eta)))
            + 1.5 * sin2_i * ((a_r).powi(2) - (a_r)) * (2.0 * (omega + nu)).cos());

    // Short-period variation in inclination
    let di_sp = -gamma2 * 0.75 * sin_i * cos_i * (a_r).powi(2) * (2.0 * (omega + nu)).cos();

    // Short-period variation in RAAN
    let draan_sp =
        -gamma2 * 1.5 * cos_i * (a_r).powi(2) * (omega + nu - m_osc + e * nu.sin()).sin();

    // Short-period variation in argument of periapsis
    let domega_sp = gamma2
        * (0.75 * (5.0 * cos_i * cos_i - 1.0) / e * ((a_r).powi(2) - (a_r)) * nu.sin()
            + 1.5 * sin2_i / (2.0 * e) * (a_r).powi(2) * (2.0 * (omega + nu)).sin());

    // Short-period variation in mean anomaly (includes secular removed)
    let dm_sp = gamma2 * eta / (2.0 * e)
        * ((1.0 - 1.5 * sin2_i) * ((a_r).powi(2) + a_r + 1.0 / (1.0 + eta)) * e * nu.sin()
            - 1.5 * sin2_i * (a_r).powi(2) * (2.0 * (omega + nu)).sin());

    // Mean = osculating - short-period
    Ok(MeanElements {
        semi_major_axis: a - da_sp,
        eccentricity: e - de_sp,
        inclination: i - di_sp,
        raan: osc.raan - draan_sp,
        argument_of_periapsis: omega - domega_sp,
        mean_anomaly: (m_osc - dm_sp).rem_euclid(std::f64::consts::TAU),
    })
}

/// Propagate mean elements forward in time using secular J2 rates.
///
/// This is much faster than numerical integration and accurate for
/// long-term prediction (days to months) when only J2 effects matter.
///
/// # Arguments
///
/// * `mean` — Mean orbital elements
/// * `mu` — Gravitational parameter (m³/s²)
/// * `j2` — J2 coefficient
/// * `r_body` — Body equatorial radius (m)
/// * `dt` — Time to propagate (seconds)
#[must_use = "returns the propagated mean elements"]
#[instrument(level = "trace", skip(mean))]
pub fn propagate_mean_elements(
    mean: &MeanElements,
    mu: f64,
    j2: f64,
    r_body: f64,
    dt: f64,
) -> Result<MeanElements> {
    let a = mean.semi_major_axis;
    let e = mean.eccentricity;
    let i = mean.inclination;

    if a <= 0.0 || mu <= 0.0 {
        return Err(FalakError::InvalidParameter(
            "semi-major axis and mu must be positive".into(),
        ));
    }

    let n = (mu / (a * a * a)).sqrt(); // mean motion
    let p = a * (1.0 - e * e);
    let ratio = r_body / p;
    let ratio2 = ratio * ratio;
    let sin2_i = i.sin().powi(2);

    // Secular rates (Brouwer first-order)
    let d_raan = -1.5 * n * j2 * ratio2 * i.cos();
    let d_omega = 1.5 * n * j2 * ratio2 * (2.0 - 2.5 * sin2_i);
    let d_m = n * (1.0 + 1.5 * j2 * ratio2 * (1.0 - 1.5 * sin2_i) / (1.0 - e * e).sqrt());

    Ok(MeanElements {
        semi_major_axis: a, // a is constant under J2 secular
        eccentricity: e,    // e is constant under J2 secular (first order)
        inclination: i,     // i is constant under J2 secular (first order)
        raan: (mean.raan + d_raan * dt).rem_euclid(std::f64::consts::TAU),
        argument_of_periapsis: (mean.argument_of_periapsis + d_omega * dt)
            .rem_euclid(std::f64::consts::TAU),
        mean_anomaly: (mean.mean_anomaly + d_m * dt).rem_euclid(std::f64::consts::TAU),
    })
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

    // ── Encke propagation ────────────────────────────────────────────

    #[test]
    fn encke_matches_cowell_j2() {
        let elem = leo_elements();
        let state = leo_state();
        let dt_prop = 500.0;

        let j2_perturb = |pos: [f64; 3], _vel: [f64; 3], _t: f64| -> [f64; 3] {
            crate::perturbation::j2_acceleration(
                pos,
                MU_EARTH,
                crate::perturbation::J2_EARTH,
                crate::perturbation::R_EARTH,
            )
        };

        let cowell_result = cowell(&state, MU_EARTH, dt_prop, 1.0, Some(&j2_perturb)).unwrap();
        let encke_result = encke(&elem, MU_EARTH, dt_prop, 1.0, &j2_perturb).unwrap();

        // Should agree within ~1 km over 500s with J2
        for i in 0..3 {
            let diff = (cowell_result.position[i] - encke_result.position[i]).abs();
            assert!(diff < 1000.0, "position[{i}] diff: {diff} m");
        }
    }

    #[test]
    fn encke_high_eccentricity() {
        // Molniya-like orbit: high eccentricity (e=0.74), verify Encke stays stable
        let elem = crate::orbit::OrbitalElements::new(26_560e3, 0.74, 1.1, 0.0, 4.7, 0.0).unwrap();
        let state = crate::kepler::elements_to_state(&elem, MU_EARTH).unwrap();

        let j2_perturb = |pos: [f64; 3], _vel: [f64; 3], _t: f64| -> [f64; 3] {
            crate::perturbation::j2_acceleration(
                pos,
                MU_EARTH,
                crate::perturbation::J2_EARTH,
                crate::perturbation::R_EARTH,
            )
        };

        let dt_prop = 1000.0;

        let cowell_result = cowell(&state, MU_EARTH, dt_prop, 1.0, Some(&j2_perturb)).unwrap();
        let encke_result = encke(&elem, MU_EARTH, dt_prop, 1.0, &j2_perturb).unwrap();

        // Both should produce finite, physically reasonable results
        for i in 0..3 {
            assert!(
                cowell_result.position[i].is_finite(),
                "Cowell pos[{i}] is not finite"
            );
            assert!(
                encke_result.position[i].is_finite(),
                "Encke pos[{i}] is not finite"
            );
        }

        // Encke and Cowell should agree within ~5 km for high-e orbit over 1000s
        for i in 0..3 {
            let diff = (cowell_result.position[i] - encke_result.position[i]).abs();
            assert!(
                diff < 5000.0,
                "high-e position[{i}] Cowell vs Encke diff: {diff:.1} m"
            );
        }
    }

    #[test]
    fn encke_very_high_eccentricity() {
        // e=0.95 — near-parabolic, stress test for deviation tracking
        let elem = crate::orbit::OrbitalElements::new(40_000e3, 0.95, 0.5, 0.0, 0.0, 0.0).unwrap();

        let j2_perturb = |pos: [f64; 3], _vel: [f64; 3], _t: f64| -> [f64; 3] {
            crate::perturbation::j2_acceleration(
                pos,
                MU_EARTH,
                crate::perturbation::J2_EARTH,
                crate::perturbation::R_EARTH,
            )
        };

        // Short propagation to avoid periapsis passage instability
        let dt_prop = 500.0;
        let encke_result = encke(&elem, MU_EARTH, dt_prop, 1.0, &j2_perturb).unwrap();

        // Result should be finite and at a reasonable distance
        let r = (encke_result.position[0].powi(2)
            + encke_result.position[1].powi(2)
            + encke_result.position[2].powi(2))
        .sqrt();
        assert!(r.is_finite(), "result radius is not finite");
        assert!(
            r > 1e6 && r < 1e9,
            "high-e orbit radius should be reasonable, got {r:.0} m"
        );
    }

    #[test]
    fn encke_invalid() {
        let elem = leo_elements();
        let no_perturb = |_p: [f64; 3], _v: [f64; 3], _t: f64| -> [f64; 3] { [0.0; 3] };
        assert!(encke(&elem, -1.0, 100.0, 10.0, &no_perturb).is_err());
        assert!(encke(&elem, MU_EARTH, 100.0, -1.0, &no_perturb).is_err());
    }

    // ── Mean elements ───────────────────────────────────────────────────

    #[test]
    fn osculating_to_mean_preserves_scale() {
        let elem = leo_elements();
        let mean = osculating_to_mean(
            &elem,
            MU_EARTH,
            crate::perturbation::J2_EARTH,
            crate::perturbation::R_EARTH,
        )
        .unwrap();

        // Mean and osculating should differ by small J2 corrections
        let da = (mean.semi_major_axis - elem.semi_major_axis).abs();
        assert!(
            da < 5000.0, // < 5 km difference for LEO (J2 short-period is ~1-2 km)
            "mean-osc SMA difference = {da:.1} m"
        );

        let de = (mean.eccentricity - elem.eccentricity).abs();
        assert!(de < 0.01, "mean-osc ecc difference = {de:.6}");

        let di = (mean.inclination - elem.inclination).abs();
        assert!(di < 0.001, "mean-osc inc difference = {di:.6} rad");
    }

    #[test]
    fn mean_elements_secular_raan_drift() {
        let elem =
            crate::orbit::OrbitalElements::new(7e6, 0.001, 51.6_f64.to_radians(), 0.0, 0.0, 0.0)
                .unwrap();
        let mean = osculating_to_mean(
            &elem,
            MU_EARTH,
            crate::perturbation::J2_EARTH,
            crate::perturbation::R_EARTH,
        )
        .unwrap();

        // Propagate mean elements for 1 day
        let mean_1d = propagate_mean_elements(
            &mean,
            MU_EARTH,
            crate::perturbation::J2_EARTH,
            crate::perturbation::R_EARTH,
            86400.0,
        )
        .unwrap();

        // RAAN should drift (ISS-like orbit drifts ~5°/day westward)
        let mut raan_drift = mean_1d.raan - mean.raan;
        if raan_drift > std::f64::consts::PI {
            raan_drift -= TAU;
        }
        if raan_drift < -std::f64::consts::PI {
            raan_drift += TAU;
        }
        let raan_drift_deg = raan_drift.to_degrees();
        assert!(
            raan_drift_deg < -3.0 && raan_drift_deg > -7.0,
            "RAAN drift = {raan_drift_deg:.2}°/day, expected ~-5°"
        );

        // SMA, ecc, inc should be unchanged (secular J2)
        assert!(
            (mean_1d.semi_major_axis - mean.semi_major_axis).abs() < 1e-6,
            "SMA should be constant under secular J2"
        );
    }

    #[test]
    fn mean_elements_invalid() {
        let hyp = crate::orbit::OrbitalElements::new(-1e7, 1.5, 0.0, 0.0, 0.0, 0.5).unwrap();
        assert!(
            osculating_to_mean(
                &hyp,
                MU_EARTH,
                crate::perturbation::J2_EARTH,
                crate::perturbation::R_EARTH,
            )
            .is_err()
        );
    }
}
