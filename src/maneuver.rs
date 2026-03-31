//! Spacecraft maneuvers — delta-v budgets, impulsive burns, continuous thrust.
//!
//! Provides types and functions for planning orbital maneuvers.

use tracing::instrument;

use crate::error::{FalakError, Result};

/// An impulsive burn (instantaneous velocity change).
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
#[non_exhaustive]
pub struct ImpulsiveBurn {
    /// Delta-v vector in the burn frame `[prograde, normal, radial]` (m/s).
    pub delta_v: [f64; 3],
    /// Total delta-v magnitude (m/s).
    pub magnitude: f64,
}

impl ImpulsiveBurn {
    /// Create a new impulsive burn from a delta-v vector.
    ///
    /// The vector is in the local orbital frame: `[prograde, normal, radial]`.
    #[must_use]
    #[inline]
    pub fn new(prograde: f64, normal: f64, radial: f64) -> Self {
        let magnitude = (prograde * prograde + normal * normal + radial * radial).sqrt();
        Self {
            delta_v: [prograde, normal, radial],
            magnitude,
        }
    }

    /// Create a prograde-only burn.
    #[must_use]
    #[inline]
    pub fn prograde(delta_v: f64) -> Self {
        Self::new(delta_v, 0.0, 0.0)
    }

    /// Create a retrograde burn (negative prograde).
    #[must_use]
    #[inline]
    pub fn retrograde(delta_v: f64) -> Self {
        Self::new(-delta_v.abs(), 0.0, 0.0)
    }

    /// Create a normal (out-of-plane) burn.
    #[must_use]
    #[inline]
    pub fn normal(delta_v: f64) -> Self {
        Self::new(0.0, delta_v, 0.0)
    }

    /// Create a radial burn.
    #[must_use]
    #[inline]
    pub fn radial(delta_v: f64) -> Self {
        Self::new(0.0, 0.0, delta_v)
    }
}

/// A burn scheduled at a specific time.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
#[non_exhaustive]
pub struct ScheduledBurn {
    /// Time offset from plan start (seconds).
    pub time: f64,
    /// The burn itself.
    pub burn: ImpulsiveBurn,
}

/// A maneuver plan: an ordered sequence of timed burns.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
#[non_exhaustive]
pub struct ManeuverPlan {
    /// Ordered sequence of burns (sorted by time).
    pub burns: Vec<ScheduledBurn>,
}

impl ManeuverPlan {
    /// Create a new empty maneuver plan.
    #[must_use]
    pub fn new() -> Self {
        Self { burns: Vec::new() }
    }

    /// Add a burn at a given time offset. Burns are kept sorted by time.
    pub fn add_burn(&mut self, time: f64, burn: ImpulsiveBurn) {
        let idx = self
            .burns
            .binary_search_by(|b| {
                b.time
                    .partial_cmp(&time)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .unwrap_or_else(|i| i);
        self.burns.insert(idx, ScheduledBurn { time, burn });
    }

    /// Remove a burn by index. Returns the removed burn.
    pub fn remove_burn(&mut self, index: usize) -> Option<ScheduledBurn> {
        if index < self.burns.len() {
            Some(self.burns.remove(index))
        } else {
            None
        }
    }

    /// Total delta-v budget (m/s).
    #[must_use]
    pub fn total_delta_v(&self) -> f64 {
        self.burns.iter().map(|b| b.burn.magnitude).sum()
    }

    /// Total time span from first to last burn (seconds).
    #[must_use]
    pub fn total_time(&self) -> f64 {
        match (self.burns.first(), self.burns.last()) {
            (Some(first), Some(last)) => last.time - first.time,
            _ => 0.0,
        }
    }

    /// Number of burns in the plan.
    #[must_use]
    #[inline]
    pub fn len(&self) -> usize {
        self.burns.len()
    }

    /// Whether the plan has no burns.
    #[must_use]
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.burns.is_empty()
    }
}

impl Default for ManeuverPlan {
    fn default() -> Self {
        Self::new()
    }
}

// ── Delta-v budget ────────────────────────────────────────────────────────

/// Compute the total delta-v required for a sequence of burns.
#[must_use]
#[inline]
pub fn total_delta_v(burns: &[ImpulsiveBurn]) -> f64 {
    burns.iter().map(|b| b.magnitude).sum()
}

// ── Rocket equation ───────────────────────────────────────────────────────

/// Compute propellant mass fraction via the Tsiolkovsky rocket equation.
///
/// m_propellant / m_initial = 1 − exp(−Δv / (Isp × g₀))
///
/// # Arguments
///
/// * `delta_v` — Required delta-v (m/s)
/// * `isp` — Specific impulse (seconds)
///
/// # Errors
///
/// Returns [`FalakError::InvalidParameter`] if delta_v is negative or isp is not positive.
#[must_use = "returns the propellant mass fraction"]
#[instrument(level = "trace")]
pub fn propellant_mass_fraction(delta_v: f64, isp: f64) -> Result<f64> {
    if delta_v < 0.0 {
        return Err(FalakError::InvalidParameter(
            format!("delta_v must be non-negative, got {delta_v}").into(),
        ));
    }
    if isp <= 0.0 {
        return Err(FalakError::InvalidParameter(
            format!("specific impulse must be positive, got {isp}").into(),
        ));
    }

    const G0: f64 = 9.806_65;
    let ve = isp * G0; // exhaust velocity
    Ok(1.0 - (-delta_v / ve).exp())
}

/// Compute maximum achievable delta-v from propellant mass fraction.
///
/// Δv = Isp × g₀ × ln(1 / (1 − fraction))
///
/// # Errors
///
/// Returns [`FalakError::InvalidParameter`] if fraction is not in [0, 1) or isp is not positive.
#[must_use = "returns the maximum delta-v"]
#[instrument(level = "trace")]
pub fn max_delta_v(propellant_fraction: f64, isp: f64) -> Result<f64> {
    if !(0.0..1.0).contains(&propellant_fraction) {
        return Err(FalakError::InvalidParameter(
            format!("propellant fraction must be in [0, 1), got {propellant_fraction}").into(),
        ));
    }
    if isp <= 0.0 {
        return Err(FalakError::InvalidParameter(
            format!("specific impulse must be positive, got {isp}").into(),
        ));
    }

    const G0: f64 = 9.806_65;
    let ve = isp * G0;
    Ok(ve * (1.0 / (1.0 - propellant_fraction)).ln())
}

// ── Escape and capture ────────────────────────────────────────────────────

/// Delta-v required to escape from a circular orbit to hyperbolic excess velocity v_inf.
///
/// Δv = √(v_inf² + v_esc²) − v_circ
///
/// # Errors
///
/// Returns [`FalakError::InvalidParameter`] if inputs are not positive.
#[must_use = "returns the escape delta-v"]
#[instrument(level = "trace")]
pub fn escape_delta_v(radius: f64, mu: f64, v_infinity: f64) -> Result<f64> {
    if radius <= 0.0 {
        return Err(FalakError::InvalidParameter(
            format!("radius must be positive, got {radius}").into(),
        ));
    }
    if mu <= 0.0 {
        return Err(FalakError::InvalidParameter(
            format!("mu must be positive, got {mu}").into(),
        ));
    }
    if v_infinity < 0.0 {
        return Err(FalakError::InvalidParameter(
            format!("v_infinity must be non-negative, got {v_infinity}").into(),
        ));
    }

    let v_circ = (mu / radius).sqrt();
    let v_esc = (2.0 * mu / radius).sqrt();
    let v_departure = (v_infinity * v_infinity + v_esc * v_esc).sqrt();

    Ok(v_departure - v_circ)
}

/// Delta-v required to capture into a circular orbit from a hyperbolic approach.
///
/// Same magnitude as escape (by symmetry).
#[must_use = "returns the capture delta-v"]
#[instrument(level = "trace")]
pub fn capture_delta_v(radius: f64, mu: f64, v_infinity: f64) -> Result<f64> {
    escape_delta_v(radius, mu, v_infinity)
}

// ── Oberth effect ─────────────────────────────────────────────────────────

/// Compute the effective delta-v gained from the Oberth effect.
///
/// Burning at higher velocity is more efficient. This returns the ratio
/// of kinetic energy gained vs. a burn at infinity.
///
/// ratio = (v + Δv)² − v² = 2vΔv + Δv²
///
/// Compared to burning from rest: Δv²
///
/// Oberth factor = (2v/Δv + 1)
#[must_use]
#[inline]
pub fn oberth_factor(orbital_velocity: f64, delta_v: f64) -> f64 {
    if delta_v.abs() < 1e-30 {
        return 1.0;
    }
    2.0 * orbital_velocity / delta_v.abs() + 1.0
}

// ── Low-thrust trajectory modelling ──────────────────────────────────────

/// Result of a low-thrust spiral transfer computation.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
#[non_exhaustive]
pub struct LowThrustTransfer {
    /// Total delta-v for the spiral (m/s).
    pub delta_v: f64,
    /// Transfer time (seconds).
    pub transfer_time: f64,
    /// Initial circular velocity (m/s).
    pub v_initial: f64,
    /// Final circular velocity (m/s).
    pub v_final: f64,
    /// Number of revolutions (approximate).
    pub revolutions: f64,
}

/// Compute a continuous low-thrust spiral transfer between circular orbits.
///
/// Uses the constant-tangential-thrust approximation: the spacecraft
/// spirals outward (or inward) with thrust always aligned along the
/// velocity vector. The delta-v for a low-thrust spiral is simply:
///
/// Δv = |v₁ − v₂| = |√(μ/r₁) − √(μ/r₂)|
///
/// This is less efficient than a Hohmann transfer but representative of
/// electric propulsion trajectories.
///
/// # Arguments
///
/// * `r1` — Initial circular orbit radius (m)
/// * `r2` — Final circular orbit radius (m)
/// * `mu` — Gravitational parameter (m³/s²)
/// * `thrust_acceleration` — Continuous thrust acceleration (m/s²)
///
/// # Errors
///
/// Returns [`FalakError::InvalidParameter`] if inputs are invalid.
#[must_use = "returns the low-thrust transfer parameters"]
#[instrument(level = "trace")]
pub fn low_thrust_spiral(
    r1: f64,
    r2: f64,
    mu: f64,
    thrust_acceleration: f64,
) -> Result<LowThrustTransfer> {
    if r1 <= 0.0 || r2 <= 0.0 || mu <= 0.0 {
        return Err(FalakError::InvalidParameter(
            "r1, r2, and mu must be positive".into(),
        ));
    }
    if thrust_acceleration <= 0.0 {
        return Err(FalakError::InvalidParameter(
            format!("thrust acceleration must be positive, got {thrust_acceleration}").into(),
        ));
    }

    let v1 = (mu / r1).sqrt();
    let v2 = (mu / r2).sqrt();
    let delta_v = (v1 - v2).abs();
    let transfer_time = delta_v / thrust_acceleration;

    // Approximate number of revolutions
    let mean_period = std::f64::consts::PI * ((r1 + r2) / 2.0).powi(3).sqrt() * 2.0 / mu.sqrt();
    let revolutions = transfer_time / mean_period;

    Ok(LowThrustTransfer {
        delta_v,
        transfer_time,
        v_initial: v1,
        v_final: v2,
        revolutions,
    })
}

/// Edelbaum low-thrust delta-v for combined altitude and plane change.
///
/// For continuous low-thrust trajectories that must change both altitude
/// and inclination, the Edelbaum approximation gives:
///
/// Δv = √(v₁² + v₂² − 2·v₁·v₂·cos(π/2 · Δi))
///
/// Note the factor of π/2 on the inclination change — this accounts for
/// the continuous steering optimisation and is specific to low-thrust.
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
/// Returns [`FalakError::InvalidParameter`] if inputs are invalid.
#[must_use = "returns the Edelbaum delta-v (m/s)"]
#[instrument(level = "trace")]
pub fn edelbaum_delta_v(r1: f64, r2: f64, delta_inclination: f64, mu: f64) -> Result<f64> {
    if r1 <= 0.0 || r2 <= 0.0 || mu <= 0.0 {
        return Err(FalakError::InvalidParameter(
            "r1, r2, and mu must be positive".into(),
        ));
    }
    if delta_inclination < 0.0 {
        return Err(FalakError::InvalidParameter(
            format!("inclination change must be non-negative, got {delta_inclination}").into(),
        ));
    }

    let v1 = (mu / r1).sqrt();
    let v2 = (mu / r2).sqrt();

    // Edelbaum formula: Δv = √(v₁² + v₂² − 2v₁v₂cos(π/2 · Δi))
    let dv = (v1 * v1 + v2 * v2
        - 2.0 * v1 * v2 * (std::f64::consts::FRAC_PI_2 * delta_inclination).cos())
    .sqrt();

    Ok(dv)
}

#[cfg(test)]
mod tests {
    use super::*;

    const MU_EARTH: f64 = 3.986_004_418e14;
    const R_LEO: f64 = 6.671e6;

    // ── ImpulsiveBurn ────────────────────────────────────────────────

    #[test]
    fn burn_prograde() {
        let b = ImpulsiveBurn::prograde(100.0);
        assert!((b.magnitude - 100.0).abs() < 1e-10);
        assert!((b.delta_v[0] - 100.0).abs() < 1e-10);
    }

    #[test]
    fn burn_retrograde() {
        let b = ImpulsiveBurn::retrograde(100.0);
        assert!((b.magnitude - 100.0).abs() < 1e-10);
        assert!((b.delta_v[0] - (-100.0)).abs() < 1e-10);
    }

    #[test]
    fn burn_combined() {
        let b = ImpulsiveBurn::new(3.0, 4.0, 0.0);
        assert!((b.magnitude - 5.0).abs() < 1e-10);
    }

    // ── ManeuverPlan ──────────────────────────────────────────────────

    #[test]
    fn plan_empty() {
        let plan = ManeuverPlan::new();
        assert!(plan.is_empty());
        assert_eq!(plan.len(), 0);
        assert!(plan.total_delta_v() < 1e-15);
        assert!(plan.total_time() < 1e-15);
    }

    #[test]
    fn plan_add_burns() {
        let mut plan = ManeuverPlan::new();
        plan.add_burn(0.0, ImpulsiveBurn::prograde(100.0));
        plan.add_burn(3600.0, ImpulsiveBurn::prograde(200.0));
        assert_eq!(plan.len(), 2);
        assert!((plan.total_delta_v() - 300.0).abs() < 1e-10);
        assert!((plan.total_time() - 3600.0).abs() < 1e-10);
    }

    #[test]
    fn plan_sorts_by_time() {
        let mut plan = ManeuverPlan::new();
        plan.add_burn(200.0, ImpulsiveBurn::prograde(50.0));
        plan.add_burn(100.0, ImpulsiveBurn::prograde(100.0));
        plan.add_burn(300.0, ImpulsiveBurn::retrograde(75.0));
        assert!((plan.burns[0].time - 100.0).abs() < 1e-10);
        assert!((plan.burns[1].time - 200.0).abs() < 1e-10);
        assert!((plan.burns[2].time - 300.0).abs() < 1e-10);
    }

    #[test]
    fn plan_remove_burn() {
        let mut plan = ManeuverPlan::new();
        plan.add_burn(0.0, ImpulsiveBurn::prograde(100.0));
        plan.add_burn(100.0, ImpulsiveBurn::prograde(200.0));
        assert!((plan.total_delta_v() - 300.0).abs() < 1e-10);

        plan.remove_burn(0);
        assert_eq!(plan.len(), 1);
        assert!((plan.total_delta_v() - 200.0).abs() < 1e-10);
    }

    #[test]
    fn plan_default() {
        let plan = ManeuverPlan::default();
        assert!(plan.is_empty());
    }

    #[test]
    fn total_dv() {
        let burns = vec![
            ImpulsiveBurn::prograde(100.0),
            ImpulsiveBurn::prograde(200.0),
            ImpulsiveBurn::normal(50.0),
        ];
        assert!((total_delta_v(&burns) - 350.0).abs() < 1e-10);
    }

    // ── Rocket equation ──────────────────────────────────────────────

    #[test]
    fn propellant_fraction_basic() {
        // Isp=300s, Δv=3000 m/s → fraction ≈ 0.641
        let frac = propellant_mass_fraction(3000.0, 300.0).unwrap();
        assert!(frac > 0.6 && frac < 0.7, "fraction: {frac}");
    }

    #[test]
    fn propellant_fraction_zero_dv() {
        let frac = propellant_mass_fraction(0.0, 300.0).unwrap();
        assert!(frac.abs() < 1e-15);
    }

    #[test]
    fn propellant_fraction_invalid() {
        assert!(propellant_mass_fraction(-1.0, 300.0).is_err());
        assert!(propellant_mass_fraction(1000.0, -1.0).is_err());
    }

    #[test]
    fn max_dv_roundtrip() {
        let dv = 5000.0;
        let isp = 350.0;
        let frac = propellant_mass_fraction(dv, isp).unwrap();
        let dv2 = max_delta_v(frac, isp).unwrap();
        assert!((dv - dv2).abs() < 1e-6, "roundtrip: {dv} vs {dv2}");
    }

    #[test]
    fn max_dv_invalid() {
        assert!(max_delta_v(-0.1, 300.0).is_err());
        assert!(max_delta_v(1.0, 300.0).is_err()); // fraction = 1.0 → infinite Δv
        assert!(max_delta_v(0.5, -1.0).is_err());
    }

    // ── Escape / capture ─────────────────────────────────────────────

    #[test]
    fn escape_from_leo() {
        // Escape from LEO: Δv ≈ 3.2 km/s (for v_inf = 0)
        let dv = escape_delta_v(R_LEO, MU_EARTH, 0.0).unwrap();
        let v_circ = (MU_EARTH / R_LEO).sqrt();
        let v_esc = (2.0 * MU_EARTH / R_LEO).sqrt();
        let expected = v_esc - v_circ;
        assert!((dv - expected).abs() < 1.0, "escape Δv: {dv} vs {expected}");
    }

    #[test]
    fn escape_with_vinf() {
        // Escape with v_inf = 3 km/s for interplanetary transfer
        let dv = escape_delta_v(R_LEO, MU_EARTH, 3000.0).unwrap();
        let dv_bare = escape_delta_v(R_LEO, MU_EARTH, 0.0).unwrap();
        assert!(dv > dv_bare, "v_inf should increase Δv");
    }

    #[test]
    fn capture_symmetric() {
        let esc = escape_delta_v(R_LEO, MU_EARTH, 2000.0).unwrap();
        let cap = capture_delta_v(R_LEO, MU_EARTH, 2000.0).unwrap();
        assert!((esc - cap).abs() < 1e-10);
    }

    #[test]
    fn escape_invalid() {
        assert!(escape_delta_v(-1.0, MU_EARTH, 0.0).is_err());
        assert!(escape_delta_v(R_LEO, -1.0, 0.0).is_err());
        assert!(escape_delta_v(R_LEO, MU_EARTH, -1.0).is_err());
    }

    // ── Oberth effect ────────────────────────────────────────────────

    #[test]
    fn oberth_at_periapsis() {
        // Burning at 7700 m/s with 1000 m/s Δv
        let factor = oberth_factor(7700.0, 1000.0);
        // = 2*7700/1000 + 1 = 16.4
        assert!((factor - 16.4).abs() < 0.01, "Oberth factor: {factor}");
    }

    #[test]
    fn oberth_zero_velocity() {
        // At rest: factor = 1 (no Oberth benefit)
        let factor = oberth_factor(0.0, 1000.0);
        assert!((factor - 1.0).abs() < 1e-10);
    }

    #[test]
    fn oberth_zero_dv() {
        let factor = oberth_factor(7700.0, 0.0);
        assert!((factor - 1.0).abs() < 1e-10);
    }

    // ── Low-thrust ──────────────────────────────────────────────────────

    #[test]
    fn low_thrust_spiral_leo_to_geo() {
        let r_leo = 6.778e6;
        let r_geo = 42_164e3;
        let accel = 1e-4; // ~10 mN on 100 kg spacecraft

        let lt = low_thrust_spiral(r_leo, r_geo, MU_EARTH, accel).unwrap();

        // Spiral Δv = |v1 - v2| ≈ |7.67 - 3.07| ≈ 4.6 km/s
        assert!(
            (lt.delta_v - 4600.0).abs() < 200.0,
            "spiral Δv = {:.0} m/s, expected ~4600",
            lt.delta_v
        );

        // Transfer time = Δv / accel
        let expected_time = lt.delta_v / accel;
        assert!((lt.transfer_time - expected_time).abs() < 1.0);

        // Many revolutions for low-thrust
        assert!(
            lt.revolutions > 10.0,
            "low-thrust should take many orbits: {}",
            lt.revolutions
        );
    }

    #[test]
    fn low_thrust_spiral_vs_hohmann() {
        // Low-thrust spiral Δv should be higher than Hohmann
        let r1 = 6.778e6;
        let r2 = 42_164e3;
        let lt = low_thrust_spiral(r1, r2, MU_EARTH, 1e-3).unwrap();
        let h = crate::transfer::hohmann(r1, r2, MU_EARTH).unwrap();
        assert!(
            lt.delta_v > h.total_delta_v,
            "spiral Δv={:.0} should exceed Hohmann Δv={:.0}",
            lt.delta_v,
            h.total_delta_v
        );
    }

    #[test]
    fn low_thrust_invalid() {
        assert!(low_thrust_spiral(-1.0, 7e6, MU_EARTH, 1e-4).is_err());
        assert!(low_thrust_spiral(7e6, 42e6, MU_EARTH, -1e-4).is_err());
    }

    #[test]
    fn edelbaum_pure_altitude() {
        // With zero inclination change, Edelbaum reduces to |v1 - v2|
        let dv = edelbaum_delta_v(6.778e6, 42_164e3, 0.0, MU_EARTH).unwrap();
        let v1 = (MU_EARTH / 6.778e6).sqrt();
        let v2 = (MU_EARTH / 42_164e3).sqrt();
        assert!(
            (dv - (v1 - v2).abs()).abs() < 1.0,
            "Edelbaum with Δi=0 should equal |v1-v2|: dv={dv:.0}"
        );
    }

    #[test]
    fn edelbaum_pure_plane_change() {
        // Same altitude: Edelbaum Δv = 2v sin(π/4 · Δi)
        let r = 7e6;
        let di = 28.5_f64.to_radians();
        let dv = edelbaum_delta_v(r, r, di, MU_EARTH).unwrap();
        let v = (MU_EARTH / r).sqrt();
        let expected = 2.0 * v * (std::f64::consts::FRAC_PI_4 * di).sin();
        assert!(
            (dv - expected).abs() < 1.0,
            "Edelbaum plane change: dv={dv:.0} vs expected={expected:.0}"
        );
    }

    #[test]
    fn edelbaum_invalid() {
        assert!(edelbaum_delta_v(-1.0, 7e6, 0.0, MU_EARTH).is_err());
        assert!(edelbaum_delta_v(7e6, 42e6, -0.1, MU_EARTH).is_err());
    }
}
