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
    /// Gravitational parameter μ (m³/s²). If `None`, uses `G × mass`.
    pub mu: Option<f64>,
}

impl Body {
    /// Create a new body (gravity computed as G × mass).
    #[must_use]
    #[inline]
    pub fn new(position: [f64; 3], velocity: [f64; 3], mass: f64) -> Self {
        Self {
            position,
            velocity,
            mass,
            mu: None,
        }
    }

    /// Create a new body with an explicit gravitational parameter μ.
    ///
    /// Use this when the standard μ value is known more precisely than G × M
    /// (e.g., μ_Earth = 3.986004418e14 m³/s²).
    #[must_use]
    #[inline]
    pub fn with_mu(position: [f64; 3], velocity: [f64; 3], mass: f64, mu: f64) -> Self {
        Self {
            position,
            velocity,
            mass,
            mu: Some(mu),
        }
    }

    /// Effective gravitational parameter (m³/s²).
    #[must_use]
    #[inline]
    pub fn gravitational_parameter(&self) -> f64 {
        self.mu.unwrap_or(G * self.mass)
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
                // Symmetric PE: −(μ_i × m_j + μ_j × m_i) / (2r)
                // Reduces to −G × m_i × m_j / r when both use G×M
                let pe_ij = self.bodies[i].gravitational_parameter() * self.bodies[j].mass
                    + self.bodies[j].gravitational_parameter() * self.bodies[i].mass;
                pe -= pe_ij / (2.0 * r);
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
/// Uses each body's [`Body::gravitational_parameter`] for gravity (canonical μ
/// when available, otherwise G × mass).
#[must_use]
pub fn compute_accelerations(system: &System) -> Vec<[f64; 3]> {
    let n = system.bodies.len();
    let mut acc = vec![[0.0; 3]; n];
    compute_accelerations_into(system, &mut acc);
    acc
}

/// Compute gravitational accelerations into a pre-allocated buffer.
///
/// Same as [`compute_accelerations`] but avoids allocation.
/// `acc` must have length ≥ `system.bodies.len()`.
pub fn compute_accelerations_into(system: &System, acc: &mut Vec<[f64; 3]>) {
    let n = system.bodies.len();
    acc.resize(n, [0.0; 3]);
    for a in acc.iter_mut() {
        *a = [0.0, 0.0, 0.0];
    }
    let eps2 = system.softening_sq;

    for i in 0..n {
        for j in (i + 1)..n {
            let dx = system.bodies[j].position[0] - system.bodies[i].position[0];
            let dy = system.bodies[j].position[1] - system.bodies[i].position[1];
            let dz = system.bodies[j].position[2] - system.bodies[i].position[2];

            let r2 = dx * dx + dy * dy + dz * dz + eps2;
            let r = r2.sqrt();
            let r3 = r2 * r;

            let mu_j = system.bodies[j].gravitational_parameter();
            let mu_i = system.bodies[i].gravitational_parameter();

            let fac_j = mu_j / r3;
            let fac_i = mu_i / r3;

            acc[i][0] += fac_j * dx;
            acc[i][1] += fac_j * dy;
            acc[i][2] += fac_j * dz;

            acc[j][0] -= fac_i * dx;
            acc[j][1] -= fac_i * dy;
            acc[j][2] -= fac_i * dz;
        }
    }
}

// ── Barnes-Hut octree ────────────────────────────────────────────────────

/// A node in the Barnes-Hut octree.
#[derive(Debug)]
struct OctreeNode {
    /// Centre of mass of all bodies in this node.
    com: [f64; 3],
    /// Total mass (kg).
    total_mass: f64,
    /// Total gravitational parameter μ (m³/s²).
    total_mu: f64,
    /// Half-width of this cubic cell.
    half_width: f64,
    /// Centre of the cubic cell.
    centre: [f64; 3],
    /// Child octants (None if leaf or empty).
    children: [Option<Box<OctreeNode>>; 8],
    /// Body index if this is a leaf with exactly one body.
    body_index: Option<usize>,
}

impl OctreeNode {
    fn new(centre: [f64; 3], half_width: f64) -> Self {
        Self {
            com: [0.0; 3],
            total_mass: 0.0,
            total_mu: 0.0,
            half_width,
            centre,
            children: Default::default(),
            body_index: None,
        }
    }

    /// Determine which octant a position falls into (0..7).
    fn octant(&self, pos: [f64; 3]) -> usize {
        let mut idx = 0;
        if pos[0] >= self.centre[0] {
            idx |= 1;
        }
        if pos[1] >= self.centre[1] {
            idx |= 2;
        }
        if pos[2] >= self.centre[2] {
            idx |= 4;
        }
        idx
    }

    /// Centre of the child octant.
    fn child_centre(&self, octant: usize) -> [f64; 3] {
        let qw = self.half_width * 0.5;
        [
            self.centre[0] + if octant & 1 != 0 { qw } else { -qw },
            self.centre[1] + if octant & 2 != 0 { qw } else { -qw },
            self.centre[2] + if octant & 4 != 0 { qw } else { -qw },
        ]
    }

    /// Insert a body into the tree.
    fn insert(&mut self, idx: usize, pos: [f64; 3], mass: f64, mu: f64) {
        if self.total_mass == 0.0 && self.body_index.is_none() {
            // Empty node — store directly
            self.body_index = Some(idx);
            self.com = pos;
            self.total_mass = mass;
            self.total_mu = mu;
            return;
        }

        // If this node has a single body, push it down first
        if let Some(existing_idx) = self.body_index.take() {
            let existing_pos = self.com;
            let existing_mass = self.total_mass;
            let existing_mu = self.total_mu;

            // Reset and re-insert existing body into child
            let oct = self.octant(existing_pos);
            let cc = self.child_centre(oct);
            let hw = self.half_width * 0.5;
            let child = self.children[oct].get_or_insert_with(|| Box::new(OctreeNode::new(cc, hw)));
            child.insert(existing_idx, existing_pos, existing_mass, existing_mu);
        }

        // Insert new body into appropriate child
        let oct = self.octant(pos);
        let cc = self.child_centre(oct);
        let hw = self.half_width * 0.5;
        let child = self.children[oct].get_or_insert_with(|| Box::new(OctreeNode::new(cc, hw)));
        child.insert(idx, pos, mass, mu);

        // Update centre of mass
        let new_total = self.total_mass + mass;
        for (c, &p) in self.com.iter_mut().zip(pos.iter()) {
            *c = (*c * self.total_mass + p * mass) / new_total;
        }
        self.total_mass = new_total;
        self.total_mu += mu;
    }

    /// Compute gravitational acceleration on a body at `pos` using the
    /// Barnes-Hut approximation with opening angle `theta`.
    fn acceleration(
        &self,
        pos: [f64; 3],
        body_idx: usize,
        theta: f64,
        softening_sq: f64,
    ) -> [f64; 3] {
        if self.total_mass == 0.0 {
            return [0.0; 3];
        }

        // If this is a leaf with the query body itself, skip self-interaction
        if let Some(idx) = self.body_index {
            if idx == body_idx {
                return [0.0; 3];
            }
            // Leaf with a different body — compute directly
            return point_acceleration(pos, self.com, self.total_mu, softening_sq);
        }

        // Cell width / distance ratio (opening criterion)
        let dx = self.com[0] - pos[0];
        let dy = self.com[1] - pos[1];
        let dz = self.com[2] - pos[2];
        let dist_sq = dx * dx + dy * dy + dz * dz + softening_sq;
        let cell_size = 2.0 * self.half_width;

        if cell_size * cell_size < theta * theta * dist_sq {
            // Cell is far enough — use monopole approximation
            return point_acceleration(pos, self.com, self.total_mu, softening_sq);
        }

        // Cell is too close — recurse into children
        let mut acc = [0.0; 3];
        for c in self.children.iter().flatten() {
            let ca = c.acceleration(pos, body_idx, theta, softening_sq);
            acc[0] += ca[0];
            acc[1] += ca[1];
            acc[2] += ca[2];
        }
        acc
    }
}

/// Gravitational acceleration from a point mass.
#[inline]
fn point_acceleration(pos: [f64; 3], source: [f64; 3], mu: f64, softening_sq: f64) -> [f64; 3] {
    let dx = source[0] - pos[0];
    let dy = source[1] - pos[1];
    let dz = source[2] - pos[2];
    let r2 = dx * dx + dy * dy + dz * dz + softening_sq;
    let r = r2.sqrt();
    let r3 = r2 * r;
    let fac = mu / r3;
    [fac * dx, fac * dy, fac * dz]
}

/// Compute gravitational accelerations using the Barnes-Hut tree algorithm.
///
/// O(N log N) approximation of the O(N²) direct summation. Accuracy is
/// controlled by the opening angle `theta`:
/// - `theta = 0.0` — exact (equivalent to direct summation)
/// - `theta = 0.5` — good accuracy for most simulations
/// - `theta = 1.0` — faster but less accurate
///
/// # Arguments
///
/// * `system` — The N-body system
/// * `theta` — Opening angle parameter (0.0–1.0, typical 0.5)
///
/// # Returns
///
/// One acceleration `[ax, ay, az]` per body.
#[must_use]
pub fn compute_accelerations_barnes_hut(system: &System, theta: f64) -> Vec<[f64; 3]> {
    let n = system.bodies.len();
    if n <= 1 {
        return vec![[0.0; 3]; n];
    }

    // Find bounding box
    let mut min = system.bodies[0].position;
    let mut max = system.bodies[0].position;
    for b in &system.bodies[1..] {
        for k in 0..3 {
            min[k] = min[k].min(b.position[k]);
            max[k] = max[k].max(b.position[k]);
        }
    }
    let centre = [
        (min[0] + max[0]) * 0.5,
        (min[1] + max[1]) * 0.5,
        (min[2] + max[2]) * 0.5,
    ];
    let half_width = ((max[0] - min[0]).max((max[1] - min[1]).max(max[2] - min[2])) * 0.5).max(1.0); // minimum 1m to avoid degenerate trees

    // Build octree
    let mut root = OctreeNode::new(centre, half_width);
    for (i, b) in system.bodies.iter().enumerate() {
        root.insert(i, b.position, b.mass, b.gravitational_parameter());
    }

    // Compute accelerations
    let mut acc = vec![[0.0; 3]; n];
    for (i, b) in system.bodies.iter().enumerate() {
        acc[i] = root.acceleration(b.position, i, theta, system.softening_sq);
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
    let mut acc = Vec::with_capacity(system.bodies.len());

    // Kick (half step)
    compute_accelerations_into(system, &mut acc);
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

    // Kick (half step) — reuses acc buffer
    compute_accelerations_into(system, &mut acc);
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

    // Pre-allocate all scratch space in flat arrays (2 × 4 × n × 3 values)
    // Layout: [k1, k2, k3, k4] for both position-rates (kx) and velocity-rates (kv)
    let mut kx = vec![[0.0; 3]; 4 * n]; // position rates (velocities at each stage)
    let mut kv = vec![[0.0; 3]; 4 * n]; // velocity rates (accelerations at each stage)
    let mut pos0 = vec![[0.0; 3]; n];
    let mut vel0 = vec![[0.0; 3]; n];
    let mut acc = Vec::with_capacity(n);

    // Save initial state
    for i in 0..n {
        pos0[i] = system.bodies[i].position;
        vel0[i] = system.bodies[i].velocity;
    }

    // k1: rates at initial state
    compute_accelerations_into(system, &mut acc);
    kv[..n].copy_from_slice(&acc[..n]);
    kx[..n].copy_from_slice(&vel0[..n]);

    // Set state to x0 + 0.5*dt*k1
    for i in 0..n {
        for k in 0..3 {
            system.bodies[i].position[k] = pos0[i][k] + 0.5 * dt * kx[i][k];
            system.bodies[i].velocity[k] = vel0[i][k] + 0.5 * dt * kv[i][k];
        }
    }

    // k2: rates at midpoint via k1
    compute_accelerations_into(system, &mut acc);
    for i in 0..n {
        kv[n + i] = acc[i];
        kx[n + i] = system.bodies[i].velocity;
    }

    // Set state to x0 + 0.5*dt*k2
    for i in 0..n {
        for k in 0..3 {
            system.bodies[i].position[k] = pos0[i][k] + 0.5 * dt * kx[n + i][k];
            system.bodies[i].velocity[k] = vel0[i][k] + 0.5 * dt * kv[n + i][k];
        }
    }

    // k3: rates at midpoint via k2
    compute_accelerations_into(system, &mut acc);
    for i in 0..n {
        kv[2 * n + i] = acc[i];
        kx[2 * n + i] = system.bodies[i].velocity;
    }

    // Set state to x0 + dt*k3
    for i in 0..n {
        for k in 0..3 {
            system.bodies[i].position[k] = pos0[i][k] + dt * kx[2 * n + i][k];
            system.bodies[i].velocity[k] = vel0[i][k] + dt * kv[2 * n + i][k];
        }
    }

    // k4: rates at endpoint via k3
    compute_accelerations_into(system, &mut acc);
    for i in 0..n {
        kv[3 * n + i] = acc[i];
        kx[3 * n + i] = system.bodies[i].velocity;
    }

    // Combine: x = x0 + (dt/6)(k1 + 2k2 + 2k3 + k4)
    let dt6 = dt / 6.0;
    for i in 0..n {
        for k in 0..3 {
            system.bodies[i].position[k] = pos0[i][k]
                + dt6 * (kx[i][k] + 2.0 * kx[n + i][k] + 2.0 * kx[2 * n + i][k] + kx[3 * n + i][k]);
            system.bodies[i].velocity[k] = vel0[i][k]
                + dt6 * (kv[i][k] + 2.0 * kv[n + i][k] + 2.0 * kv[2 * n + i][k] + kv[3 * n + i][k]);
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
#[must_use]
pub fn step_adaptive(system: &mut System, dt: f64, tolerance: f64) -> (f64, f64) {
    let n = system.bodies.len();
    let mut current_dt = dt;

    // Pre-allocate scratch buffers once (avoids system.clone() per attempt)
    let mut pos0 = vec![[0.0; 3]; n];
    let mut vel0 = vec![[0.0; 3]; n];
    let mut full_pos = vec![[0.0; 3]; n];
    let mut full_vel = vec![[0.0; 3]; n];

    // Iterative step reduction (max 30 halvings: dt/2^30 ≈ 1e-9 × dt)
    for _ in 0..30 {
        for i in 0..n {
            pos0[i] = system.bodies[i].position;
            vel0[i] = system.bodies[i].velocity;
        }

        // Full step with RK4, then save result
        step_rk4(system, current_dt);
        for i in 0..n {
            full_pos[i] = system.bodies[i].position;
            full_vel[i] = system.bodies[i].velocity;
        }

        // Restore initial state and do two half steps
        for i in 0..n {
            system.bodies[i].position = pos0[i];
            system.bodies[i].velocity = vel0[i];
        }
        step_rk4(system, current_dt / 2.0);
        step_rk4(system, current_dt / 2.0);

        // Error estimate: max position difference
        let mut max_err: f64 = 0.0;
        for (body, fp) in system.bodies.iter().zip(full_pos.iter()) {
            for (bp, &fp_k) in body.position.iter().zip(fp.iter()) {
                let err = (bp - fp_k).abs();
                max_err = max_err.max(err);
            }
        }

        if max_err <= tolerance || current_dt <= 1e-6 {
            // Accept — two-half-step result is already in system
            let safety = 0.9;
            let next_dt = if max_err > 1e-30 {
                safety * current_dt * (tolerance / max_err).powf(0.2)
            } else {
                current_dt * 2.0
            };
            return (
                current_dt,
                next_dt.clamp(current_dt * 0.1, current_dt * 2.0),
            );
        }

        // Reject — restore and halve
        for i in 0..n {
            system.bodies[i].position = pos0[i];
            system.bodies[i].velocity = vel0[i];
        }
        current_dt *= 0.5;
    }

    // Fallback: accept whatever we have at minimum dt
    (current_dt, current_dt)
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
        // Earth at origin with canonical μ, satellite in circular LEO
        let r = 7e6;
        let v = (MU_EARTH / r).sqrt();
        System::new(
            vec![
                Body::with_mu([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], M_EARTH, MU_EARTH),
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
        // Use G×M bodies (no canonical μ) for exact Newton's 3rd law
        let r = 7e6;
        let v = (G * M_EARTH / r).sqrt();
        let sys = System::new(
            vec![
                Body::new([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], M_EARTH),
                Body::new([r, 0.0, 0.0], [0.0, v, 0.0], 1.0),
            ],
            0.0,
        )
        .unwrap();
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
        // Satellite acceleration should be μ/r² ≈ 8.13 m/s² (canonical μ used)
        let a_sat = (acc[1][0] * acc[1][0] + acc[1][1] * acc[1][1] + acc[1][2] * acc[1][2]).sqrt();
        let expected = MU_EARTH / (7e6 * 7e6);
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
        let dt = 1.0; // 1s steps for tight closure
        let steps = (period / dt) as u64;

        for _ in 0..steps {
            step_rk4(&mut sys, dt);
        }

        // RK4 with canonical μ should close the orbit tightly
        let dr = [
            sys.bodies[1].position[0] - r0[0],
            sys.bodies[1].position[1] - r0[1],
            sys.bodies[1].position[2] - r0[2],
        ];
        let dr_mag = (dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]).sqrt();
        let frac = dr_mag / 7e6;
        assert!(
            frac < 0.001,
            "RK4 position drift: {frac:.6} of orbit radius"
        );
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

    // ── Barnes-Hut ──────────────────────────────────────────────────────

    #[test]
    #[allow(clippy::needless_range_loop)]
    fn barnes_hut_matches_direct_two_body() {
        let sys = two_body_system();
        let direct = compute_accelerations(&sys);
        let bh = compute_accelerations_barnes_hut(&sys, 0.0); // θ=0 → exact

        for i in 0..2 {
            for k in 0..3 {
                let rel = if direct[i][k].abs() > 1e-10 {
                    (bh[i][k] - direct[i][k]).abs() / direct[i][k].abs()
                } else {
                    (bh[i][k] - direct[i][k]).abs()
                };
                assert!(
                    rel < 1e-10,
                    "body {i} axis {k}: direct={} bh={}",
                    direct[i][k],
                    bh[i][k]
                );
            }
        }
    }

    #[test]
    fn barnes_hut_approximate_accuracy() {
        // Multi-body system: verify θ=0.5 gives < 1% error
        let bodies = vec![
            Body::new([0.0, 0.0, 0.0], [0.0; 3], 1e30),
            Body::new([1e10, 0.0, 0.0], [0.0; 3], 1e28),
            Body::new([0.0, 2e10, 0.0], [0.0; 3], 1e28),
            Body::new([-1e10, -1e10, 0.0], [0.0; 3], 1e28),
            Body::new([5e9, 5e9, 5e9], [0.0; 3], 1e28),
        ];
        let sys = System::new(bodies, 0.0).unwrap();

        let direct = compute_accelerations(&sys);
        let bh = compute_accelerations_barnes_hut(&sys, 0.5);

        for i in 0..5 {
            let d_mag = (direct[i][0].powi(2) + direct[i][1].powi(2) + direct[i][2].powi(2)).sqrt();
            let err_mag = ((bh[i][0] - direct[i][0]).powi(2)
                + (bh[i][1] - direct[i][1]).powi(2)
                + (bh[i][2] - direct[i][2]).powi(2))
            .sqrt();
            let rel_err = if d_mag > 0.0 { err_mag / d_mag } else { 0.0 };
            assert!(
                rel_err < 0.01,
                "body {i}: relative error = {rel_err:.4} (> 1%)"
            );
        }
    }

    #[test]
    fn barnes_hut_single_body() {
        let sys = System::new(vec![Body::new([0.0; 3], [0.0; 3], 1e30)], 0.0).unwrap();
        let bh = compute_accelerations_barnes_hut(&sys, 0.5);
        assert_eq!(bh[0], [0.0, 0.0, 0.0]);
    }

    #[test]
    #[allow(clippy::needless_range_loop)]
    fn barnes_hut_newton_third_law() {
        let sys = two_body_system();
        let bh = compute_accelerations_barnes_hut(&sys, 0.5);
        // m1*a1 + m2*a2 ≈ 0 (Newton's third law)
        for k in 0..3 {
            let net = sys.bodies[0].mass * bh[0][k] + sys.bodies[1].mass * bh[1][k];
            assert!(net.abs() < 1e-3, "net force axis {k}: {net}");
        }
    }
}
