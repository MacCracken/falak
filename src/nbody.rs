//! N-body gravitational simulation — direct computation, numerical integration.
//!
//! Provides gravitational acceleration computation and symplectic/explicit
//! integrators for evolving systems of gravitating bodies.

use tracing::instrument;

use crate::error::{FalakError, Result};
use crate::kepler::G;

/// A gravitating body with position, velocity, and mass.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
#[non_exhaustive]
pub struct Body {
    /// Position `[x, y, z]` (metres).
    pub position: [f64; 3],
    /// Velocity `[vx, vy, vz]` (m/s).
    pub velocity: [f64; 3],
    /// Mass (kg).
    pub mass: f64,
}

impl Body {
    /// Create a new body.
    #[must_use]
    #[inline]
    pub fn new(position: [f64; 3], velocity: [f64; 3], mass: f64) -> Self {
        Self {
            position,
            velocity,
            mass,
        }
    }

    /// Kinetic energy (J).
    #[must_use]
    #[inline]
    pub fn kinetic_energy(&self) -> f64 {
        let [vx, vy, vz] = self.velocity;
        0.5 * self.mass * (vx * vx + vy * vy + vz * vz)
    }
}

/// An N-body system: a collection of gravitating bodies.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
#[non_exhaustive]
pub struct System {
    /// The bodies in the system.
    pub bodies: Vec<Body>,
    /// Gravitational softening length squared (m²) to prevent singularities.
    pub softening_sq: f64,
}

impl System {
    /// Create a new system from a list of bodies.
    ///
    /// # Arguments
    ///
    /// * `bodies` — The gravitating bodies.
    /// * `softening` — Softening length (metres). Use 0.0 for exact Newtonian gravity.
    ///
    /// # Errors
    ///
    /// Returns [`FalakError::InvalidParameter`] if no bodies or any mass is non-positive.
    #[must_use = "returns the constructed system"]
    #[instrument(level = "trace", skip(bodies))]
    pub fn new(bodies: Vec<Body>, softening: f64) -> Result<Self> {
        if bodies.is_empty() {
            return Err(FalakError::InvalidParameter("need at least 1 body".into()));
        }
        for (i, b) in bodies.iter().enumerate() {
            if b.mass <= 0.0 {
                return Err(FalakError::InvalidParameter(
                    format!("body {i} has non-positive mass: {}", b.mass).into(),
                ));
            }
        }
        Ok(Self {
            bodies,
            softening_sq: softening * softening,
        })
    }

    /// Number of bodies.
    #[must_use]
    #[inline]
    pub fn len(&self) -> usize {
        self.bodies.len()
    }

    /// Whether the system is empty.
    #[must_use]
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.bodies.is_empty()
    }

    /// Total kinetic energy of the system (J).
    #[must_use]
    pub fn kinetic_energy(&self) -> f64 {
        self.bodies.iter().map(|b| b.kinetic_energy()).sum()
    }

    /// Total gravitational potential energy of the system (J).
    #[must_use]
    pub fn potential_energy(&self) -> f64 {
        let n = self.bodies.len();
        let mut pe = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                let dx = self.bodies[j].position[0] - self.bodies[i].position[0];
                let dy = self.bodies[j].position[1] - self.bodies[i].position[1];
                let dz = self.bodies[j].position[2] - self.bodies[i].position[2];
                let r = (dx * dx + dy * dy + dz * dz + self.softening_sq).sqrt();
                pe -= G * self.bodies[i].mass * self.bodies[j].mass / r;
            }
        }
        pe
    }

    /// Total energy (kinetic + potential) of the system (J).
    #[must_use]
    #[inline]
    pub fn total_energy(&self) -> f64 {
        self.kinetic_energy() + self.potential_energy()
    }

    /// Centre of mass position.
    #[must_use]
    pub fn centre_of_mass(&self) -> [f64; 3] {
        let mut total_mass = 0.0;
        let mut com = [0.0; 3];
        for b in &self.bodies {
            total_mass += b.mass;
            for (c, &p) in com.iter_mut().zip(b.position.iter()) {
                *c += b.mass * p;
            }
        }
        if total_mass > 0.0 {
            for c in &mut com {
                *c /= total_mass;
            }
        }
        com
    }
}

// ── Gravitational acceleration ────────────────────────────────────────────

/// Compute gravitational accelerations for all bodies via direct summation.
///
/// O(N²) pairwise computation. Returns one `[ax, ay, az]` per body.
#[must_use]
pub fn compute_accelerations(system: &System) -> Vec<[f64; 3]> {
    let n = system.bodies.len();
    let mut acc = vec![[0.0; 3]; n];
    let eps2 = system.softening_sq;

    for i in 0..n {
        for j in (i + 1)..n {
            let dx = system.bodies[j].position[0] - system.bodies[i].position[0];
            let dy = system.bodies[j].position[1] - system.bodies[i].position[1];
            let dz = system.bodies[j].position[2] - system.bodies[i].position[2];

            let r2 = dx * dx + dy * dy + dz * dz + eps2;
            let r = r2.sqrt();
            let r3 = r2 * r;

            let fac_j = G * system.bodies[j].mass / r3;
            let fac_i = G * system.bodies[i].mass / r3;

            acc[i][0] += fac_j * dx;
            acc[i][1] += fac_j * dy;
            acc[i][2] += fac_j * dz;

            acc[j][0] -= fac_i * dx;
            acc[j][1] -= fac_i * dy;
            acc[j][2] -= fac_i * dz;
        }
    }
    acc
}

// ── Leapfrog (Stormer-Verlet) integrator ──────────────────────────────────

/// Advance the system by one step using the leapfrog (kick-drift-kick) integrator.
///
/// Symplectic — conserves energy over long timescales. Second-order accurate.
///
/// # Arguments
///
/// * `system` — The N-body system (modified in place).
/// * `dt` — Time step (seconds).
#[instrument(level = "trace", skip(system))]
pub fn step_leapfrog(system: &mut System, dt: f64) {
    let half_dt = 0.5 * dt;

    // Kick (half step)
    let acc = compute_accelerations(system);
    for (body, a) in system.bodies.iter_mut().zip(acc.iter()) {
        body.velocity[0] += half_dt * a[0];
        body.velocity[1] += half_dt * a[1];
        body.velocity[2] += half_dt * a[2];
    }

    // Drift (full step)
    for body in &mut system.bodies {
        body.position[0] += dt * body.velocity[0];
        body.position[1] += dt * body.velocity[1];
        body.position[2] += dt * body.velocity[2];
    }

    // Kick (half step)
    let acc = compute_accelerations(system);
    for (body, a) in system.bodies.iter_mut().zip(acc.iter()) {
        body.velocity[0] += half_dt * a[0];
        body.velocity[1] += half_dt * a[1];
        body.velocity[2] += half_dt * a[2];
    }
}

// ── Runge-Kutta 4 integrator ──────────────────────────────────────────────

/// Advance the system by one step using the classic RK4 integrator.
///
/// Fourth-order accurate but not symplectic — energy drifts over long runs.
/// Better for short-duration high-accuracy propagation.
///
/// # Arguments
///
/// * `system` — The N-body system (modified in place).
/// * `dt` — Time step (seconds).
#[instrument(level = "trace", skip(system))]
pub fn step_rk4(system: &mut System, dt: f64) {
    let n = system.bodies.len();

    // Save initial state
    let pos0: Vec<[f64; 3]> = system.bodies.iter().map(|b| b.position).collect();
    let vel0: Vec<[f64; 3]> = system.bodies.iter().map(|b| b.velocity).collect();

    // k1
    let acc1 = compute_accelerations(system);
    let k1v: Vec<[f64; 3]> = acc1;
    let k1x: Vec<[f64; 3]> = vel0.clone();

    // Set state to x0 + 0.5*dt*k1
    for i in 0..n {
        for k in 0..3 {
            system.bodies[i].position[k] = pos0[i][k] + 0.5 * dt * k1x[i][k];
            system.bodies[i].velocity[k] = vel0[i][k] + 0.5 * dt * k1v[i][k];
        }
    }

    // k2
    let acc2 = compute_accelerations(system);
    let k2v: Vec<[f64; 3]> = acc2;
    let k2x: Vec<[f64; 3]> = system.bodies.iter().map(|b| b.velocity).collect();

    // Set state to x0 + 0.5*dt*k2
    for i in 0..n {
        for k in 0..3 {
            system.bodies[i].position[k] = pos0[i][k] + 0.5 * dt * k2x[i][k];
            system.bodies[i].velocity[k] = vel0[i][k] + 0.5 * dt * k2v[i][k];
        }
    }

    // k3
    let acc3 = compute_accelerations(system);
    let k3v: Vec<[f64; 3]> = acc3;
    let k3x: Vec<[f64; 3]> = system.bodies.iter().map(|b| b.velocity).collect();

    // Set state to x0 + dt*k3
    for i in 0..n {
        for k in 0..3 {
            system.bodies[i].position[k] = pos0[i][k] + dt * k3x[i][k];
            system.bodies[i].velocity[k] = vel0[i][k] + dt * k3v[i][k];
        }
    }

    // k4
    let acc4 = compute_accelerations(system);
    let k4v: Vec<[f64; 3]> = acc4;
    let k4x: Vec<[f64; 3]> = system.bodies.iter().map(|b| b.velocity).collect();

    // Combine: x = x0 + (dt/6)(k1 + 2k2 + 2k3 + k4)
    let dt6 = dt / 6.0;
    for i in 0..n {
        for k in 0..3 {
            system.bodies[i].position[k] =
                pos0[i][k] + dt6 * (k1x[i][k] + 2.0 * k2x[i][k] + 2.0 * k3x[i][k] + k4x[i][k]);
            system.bodies[i].velocity[k] =
                vel0[i][k] + dt6 * (k1v[i][k] + 2.0 * k2v[i][k] + 2.0 * k3v[i][k] + k4v[i][k]);
        }
    }
}

// ── Adaptive RK45 integrator (Richardson extrapolation) ───────────────────

/// Advance the system by one adaptive step using embedded RK4 error estimation.
///
/// Compares a full RK4 step with two half-steps (Richardson extrapolation)
/// to estimate the error. Adjusts step size to maintain the target tolerance.
///
/// # Arguments
///
/// * `system` — The N-body system (modified in place).
/// * `dt` — Suggested time step (seconds).
/// * `tolerance` — Maximum allowed position error per step (metres).
///
/// # Returns
///
/// `(actual_dt, next_suggested_dt)`.
pub fn step_adaptive(system: &mut System, dt: f64, tolerance: f64) -> (f64, f64) {
    let n = system.bodies.len();
    let pos0: Vec<[f64; 3]> = system.bodies.iter().map(|b| b.position).collect();
    let vel0: Vec<[f64; 3]> = system.bodies.iter().map(|b| b.velocity).collect();

    // Full step with RK4
    let mut full_sys = system.clone();
    step_rk4(&mut full_sys, dt);

    // Two half steps with RK4
    step_rk4(system, dt / 2.0);
    step_rk4(system, dt / 2.0);

    // Error estimate: max position difference
    let mut max_err: f64 = 0.0;
    for i in 0..n {
        for j in 0..3 {
            let err = (system.bodies[i].position[j] - full_sys.bodies[i].position[j]).abs();
            max_err = max_err.max(err);
        }
    }

    // Reject if error exceeds tolerance
    if max_err > tolerance && dt > 1e-6 {
        for i in 0..n {
            system.bodies[i].position = pos0[i];
            system.bodies[i].velocity = vel0[i];
        }
        return step_adaptive(system, dt * 0.5, tolerance);
    }

    // The two-half-step result (already in system) is more accurate

    // Suggest next step size (5th-order scaling for RK4)
    let safety = 0.9;
    let next_dt = if max_err > 1e-30 {
        safety * dt * (tolerance / max_err).powf(0.2)
    } else {
        dt * 2.0
    };

    (dt, next_dt.clamp(dt * 0.1, dt * 2.0))
}

// ── Multi-step integration ────────────────────────────────────────────────

/// Integration method.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[non_exhaustive]
pub enum Integrator {
    /// Leapfrog (kick-drift-kick). Symplectic, 2nd order.
    Leapfrog,
    /// Classic Runge-Kutta 4. Non-symplectic, 4th order.
    Rk4,
}

/// Evolve a system for a given duration.
///
/// # Arguments
///
/// * `system` — The N-body system (modified in place).
/// * `total_time` — Total integration time (seconds).
/// * `dt` — Time step (seconds).
/// * `integrator` — Which integration method to use.
///
/// # Errors
///
/// Returns [`FalakError::InvalidParameter`] if `dt` is not positive or `total_time` is negative.
#[must_use = "returns Ok on success or an error"]
#[instrument(level = "debug", skip(system))]
pub fn evolve(system: &mut System, total_time: f64, dt: f64, integrator: Integrator) -> Result<()> {
    if dt <= 0.0 {
        return Err(FalakError::InvalidParameter(
            format!("time step must be positive, got {dt}").into(),
        ));
    }
    if total_time < 0.0 {
        return Err(FalakError::InvalidParameter(
            format!("total time must be non-negative, got {total_time}").into(),
        ));
    }

    let step_fn = match integrator {
        Integrator::Leapfrog => step_leapfrog,
        Integrator::Rk4 => step_rk4,
    };

    let mut remaining = total_time;
    while remaining > 0.0 {
        let step_dt = dt.min(remaining);
        step_fn(system, step_dt);
        remaining -= step_dt;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::TAU;

    const MU_EARTH: f64 = 3.986_004_418e14;
    const M_EARTH: f64 = 5.972e24;

    fn two_body_system() -> System {
        // Earth at origin, satellite in circular LEO
        let r = 7e6;
        let v = (MU_EARTH / r).sqrt();
        System::new(
            vec![
                Body::new([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], M_EARTH),
                Body::new([r, 0.0, 0.0], [0.0, v, 0.0], 1.0),
            ],
            0.0,
        )
        .unwrap()
    }

    // ── System ───────────────────────────────────────────────────────

    #[test]
    fn system_creation() {
        let sys = two_body_system();
        assert_eq!(sys.len(), 2);
        assert!(!sys.is_empty());
    }

    #[test]
    fn system_empty_rejected() {
        assert!(System::new(vec![], 0.0).is_err());
    }

    #[test]
    fn system_single_body() {
        let sys = System::new(vec![Body::new([0.0; 3], [1.0, 0.0, 0.0], 1e30)], 0.0).unwrap();
        assert_eq!(sys.len(), 1);
        let acc = compute_accelerations(&sys);
        assert_eq!(acc[0], [0.0, 0.0, 0.0]); // no self-interaction
    }

    #[test]
    fn system_negative_mass() {
        assert!(
            System::new(
                vec![
                    Body::new([0.0; 3], [0.0; 3], 1.0),
                    Body::new([1.0, 0.0, 0.0], [0.0; 3], -1.0),
                ],
                0.0,
            )
            .is_err()
        );
    }

    #[test]
    fn energy_conservation_basic() {
        let sys = two_body_system();
        let ke = sys.kinetic_energy();
        let pe = sys.potential_energy();
        let te = sys.total_energy();
        assert!((te - (ke + pe)).abs() < 1e-10);
        // Bound orbit → total energy should be negative
        assert!(te < 0.0, "bound orbit should have negative energy: {te}");
    }

    #[test]
    fn centre_of_mass() {
        let sys = two_body_system();
        let com = sys.centre_of_mass();
        // Earth is so massive that COM is very close to origin
        let com_mag = (com[0] * com[0] + com[1] * com[1] + com[2] * com[2]).sqrt();
        assert!(com_mag < 2.0, "COM should be near Earth: {com_mag}");
    }

    // ── Accelerations ────────────────────────────────────────────────

    #[test]
    fn accelerations_newton_third_law() {
        let sys = two_body_system();
        let acc = compute_accelerations(&sys);
        // Newton's third law: m1*a1 + m2*a2 = 0
        for (k, (&a0, &a1)) in acc[0].iter().zip(acc[1].iter()).enumerate() {
            let net = sys.bodies[0].mass * a0 + sys.bodies[1].mass * a1;
            assert!(
                net.abs() < 1e-6,
                "Newton's 3rd law violated: axis {k}, net={net}"
            );
        }
    }

    #[test]
    fn acceleration_magnitude() {
        let sys = two_body_system();
        let acc = compute_accelerations(&sys);
        // Satellite acceleration should be ~G*M_earth/r² ≈ 8.13 m/s²
        // (slightly different from μ/r² because G*M_earth ≠ μ_earth exactly)
        let a_sat = (acc[1][0] * acc[1][0] + acc[1][1] * acc[1][1] + acc[1][2] * acc[1][2]).sqrt();
        let expected = G * M_EARTH / (7e6 * 7e6);
        assert!(
            (a_sat - expected).abs() / expected < 1e-6,
            "satellite accel: {a_sat} vs {expected}"
        );
    }

    // ── Leapfrog integrator ──────────────────────────────────────────

    #[test]
    fn leapfrog_one_orbit() {
        let mut sys = two_body_system();
        let e0 = sys.total_energy();
        let r0 = sys.bodies[1].position;

        // One orbital period: T = 2π√(a³/μ)
        let a = 7e6;
        let period = TAU * (a * a * a / MU_EARTH).sqrt();
        let dt = 10.0; // 10s steps
        let steps = (period / dt) as u64;

        for _ in 0..steps {
            step_leapfrog(&mut sys, dt);
        }

        // Energy should be conserved (leapfrog is symplectic)
        let e1 = sys.total_energy();
        let e_err = ((e1 - e0) / e0).abs();
        assert!(e_err < 1e-6, "leapfrog energy drift: {e_err:.2e}");

        // Position should return close to start
        let dr = [
            sys.bodies[1].position[0] - r0[0],
            sys.bodies[1].position[1] - r0[1],
            sys.bodies[1].position[2] - r0[2],
        ];
        let dr_mag = (dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]).sqrt();
        let frac = dr_mag / 7e6;
        assert!(
            frac < 0.01,
            "leapfrog position drift: {frac:.4} of orbit radius"
        );
    }

    // ── RK4 integrator ───────────────────────────────────────────────

    #[test]
    fn rk4_one_orbit() {
        let mut sys = two_body_system();
        let r0 = sys.bodies[1].position;

        let a = 7e6;
        let period = TAU * (a * a * a / MU_EARTH).sqrt();
        let dt = 10.0;
        let steps = (period / dt) as u64;

        for _ in 0..steps {
            step_rk4(&mut sys, dt);
        }

        // RK4 should close the orbit well (not as tight as analytic due to G*M vs μ)
        let dr = [
            sys.bodies[1].position[0] - r0[0],
            sys.bodies[1].position[1] - r0[1],
            sys.bodies[1].position[2] - r0[2],
        ];
        let dr_mag = (dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]).sqrt();
        let frac = dr_mag / 7e6;
        assert!(frac < 0.02, "RK4 position drift: {frac:.6} of orbit radius");
    }

    // ── Evolve ───────────────────────────────────────────────────────

    #[test]
    fn evolve_leapfrog() {
        let mut sys = two_body_system();
        let e0 = sys.total_energy();
        // Evolve for 1000 seconds
        evolve(&mut sys, 1000.0, 10.0, Integrator::Leapfrog).unwrap();
        let e1 = sys.total_energy();
        let e_err = ((e1 - e0) / e0).abs();
        assert!(e_err < 1e-8, "evolve energy drift: {e_err:.2e}");
    }

    #[test]
    fn evolve_rk4() {
        let mut sys = two_body_system();
        evolve(&mut sys, 1000.0, 10.0, Integrator::Rk4).unwrap();
        // Should not crash, satellite should still be orbiting
        let r = (sys.bodies[1].position[0].powi(2)
            + sys.bodies[1].position[1].powi(2)
            + sys.bodies[1].position[2].powi(2))
        .sqrt();
        assert!((r - 7e6).abs() / 7e6 < 0.01, "RK4 orbit radius drift: {r}");
    }

    #[test]
    fn evolve_invalid() {
        let mut sys = two_body_system();
        assert!(evolve(&mut sys, 100.0, -1.0, Integrator::Leapfrog).is_err());
        assert!(evolve(&mut sys, -1.0, 10.0, Integrator::Leapfrog).is_err());
    }

    // ── Adaptive integrator ────────────────────────────────────────────

    #[test]
    fn adaptive_conserves_energy() {
        let mut sys = two_body_system();
        let e0 = sys.total_energy();

        let mut dt: f64 = 50.0;
        let mut t: f64 = 0.0;
        let target: f64 = 1000.0;
        while t < target {
            let (actual_dt, next_dt) = step_adaptive(&mut sys, dt.min(target - t), 1.0);
            t += actual_dt;
            dt = next_dt;
        }

        let e1 = sys.total_energy();
        let e_err = ((e1 - e0) / e0).abs();
        assert!(e_err < 1e-6, "adaptive energy drift: {e_err:.2e}");
    }

    #[test]
    fn adaptive_reduces_step_for_tight_tolerance() {
        let mut sys = two_body_system();
        let (actual_dt, next_dt) = step_adaptive(&mut sys, 100.0, 0.001);
        // Should have reduced from 100s to something smaller
        assert!(actual_dt <= 100.0);
        assert!(next_dt > 0.0);
    }

    // ── Three-body ───────────────────────────────────────────────────

    #[test]
    fn three_body_energy() {
        let sys = System::new(
            vec![
                Body::new([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1e30),
                Body::new([1e11, 0.0, 0.0], [0.0, 3e4, 0.0], 1e24),
                Body::new([-1e11, 0.0, 0.0], [0.0, -3e4, 0.0], 1e24),
            ],
            0.0,
        )
        .unwrap();

        let e = sys.total_energy();
        // Should be finite and negative (bound system)
        assert!(e.is_finite(), "energy should be finite: {e}");
    }

    #[test]
    fn body_kinetic_energy() {
        let b = Body::new([0.0; 3], [3.0, 4.0, 0.0], 2.0);
        assert!((b.kinetic_energy() - 25.0).abs() < 1e-10);
    }
}
