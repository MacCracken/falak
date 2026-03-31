//! Circular Restricted Three-Body Problem (CR3BP).
//!
//! Models the motion of a massless body under the gravitational influence of
//! two massive bodies (primaries) in circular orbit about their barycentre.
//! All quantities are in the **synodic (rotating) frame** normalised so that
//! the total mass is 1, the distance between primaries is 1, and the angular
//! velocity is 1.
//!
//! Key results:
//! - Five Lagrange points (L1–L5)
//! - Jacobi constant (energy-like integral of motion)
//! - Equations of motion in the rotating frame
//! - Zero-velocity curves

use tracing::instrument;

use crate::error::{FalakError, Result};

/// Mass ratio μ = m₂ / (m₁ + m₂), where m₂ is the smaller primary.
///
/// Common values:
/// - Sun–Earth: ≈ 3.003e-6
/// - Earth–Moon: ≈ 0.01215
/// - Sun–Jupiter: ≈ 9.537e-4
pub type MassRatio = f64;

/// State in the synodic (rotating) frame: `[x, y, z, vx, vy, vz]`.
///
/// The x-axis points from the barycentre toward the smaller primary (m₂).
/// The z-axis is along the rotation axis. The y-axis completes the
/// right-handed triad.
pub type SynodicState = [f64; 6];

/// The five Lagrange (libration) points in the CR3BP.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
#[non_exhaustive]
pub struct LagrangePoints {
    /// L1: between the two primaries.
    pub l1: [f64; 3],
    /// L2: beyond the smaller primary.
    pub l2: [f64; 3],
    /// L3: beyond the larger primary (opposite side).
    pub l3: [f64; 3],
    /// L4: leading equilateral point (60° ahead of m₂).
    pub l4: [f64; 3],
    /// L5: trailing equilateral point (60° behind m₂).
    pub l5: [f64; 3],
}

/// Compute all five Lagrange points for a given mass ratio.
///
/// Uses Newton's method to find the collinear points (L1, L2, L3) and
/// the analytic solution for the triangular points (L4, L5).
///
/// # Arguments
///
/// * `mu` — Mass ratio m₂/(m₁+m₂) where m₂ is the smaller body (0 < μ < 0.5).
///
/// # Errors
///
/// Returns [`FalakError::InvalidParameter`] if μ is not in (0, 0.5].
#[must_use = "returns the five Lagrange points"]
#[instrument(level = "trace")]
pub fn lagrange_points(mu: MassRatio) -> Result<LagrangePoints> {
    if mu <= 0.0 || mu > 0.5 {
        return Err(FalakError::InvalidParameter(
            format!("mass ratio must be in (0, 0.5], got {mu}").into(),
        ));
    }

    let l1_x = find_collinear_l1(mu);
    let l2_x = find_collinear_l2(mu);
    let l3_x = find_collinear_l3(mu);

    Ok(LagrangePoints {
        l1: [l1_x, 0.0, 0.0],
        l2: [l2_x, 0.0, 0.0],
        l3: [l3_x, 0.0, 0.0],
        // L4: equilateral, 60° ahead of m₂
        l4: [0.5 - mu, 3.0_f64.sqrt() / 2.0, 0.0],
        // L5: equilateral, 60° behind m₂
        l5: [0.5 - mu, -3.0_f64.sqrt() / 2.0, 0.0],
    })
}

/// Compute the Jacobi constant (integral of motion) for a state in the CR3BP.
///
/// C_J = 2Ω(x,y,z) - (vx² + vy² + vz²)
///
/// where Ω is the pseudo-potential including centrifugal and gravitational terms.
///
/// Higher C_J means lower energy; the zero-velocity surface C_J = 2Ω bounds
/// the regions accessible to the particle.
///
/// # Arguments
///
/// * `state` — Position and velocity in the synodic frame `[x, y, z, vx, vy, vz]`.
/// * `mu` — Mass ratio.
#[must_use]
#[inline]
pub fn jacobi_constant(state: &SynodicState, mu: MassRatio) -> f64 {
    let [x, y, z, vx, vy, vz] = *state;
    let omega = pseudo_potential(x, y, z, mu);
    2.0 * omega - (vx * vx + vy * vy + vz * vz)
}

/// Pseudo-potential Ω(x, y, z) in the synodic frame.
///
/// Ω = ½(x² + y²) + (1-μ)/r₁ + μ/r₂
///
/// where r₁, r₂ are distances to the two primaries.
#[must_use]
#[inline]
pub fn pseudo_potential(x: f64, y: f64, z: f64, mu: MassRatio) -> f64 {
    let r1 = distance_to_m1(x, y, z, mu);
    let r2 = distance_to_m2(x, y, z, mu);
    0.5 * (x * x + y * y) + (1.0 - mu) / r1 + mu / r2
}

/// Evaluate the zero-velocity surface value: C_J = 2Ω(x, y, z).
///
/// Points where C_J ≤ 2Ω are accessible; points where C_J > 2Ω are forbidden.
#[must_use]
#[inline]
pub fn zero_velocity_value(x: f64, y: f64, z: f64, mu: MassRatio) -> f64 {
    2.0 * pseudo_potential(x, y, z, mu)
}

/// Equations of motion in the CR3BP rotating frame.
///
/// Given state `[x, y, z, vx, vy, vz]`, returns the derivative
/// `[vx, vy, vz, ax, ay, az]` where the accelerations include
/// Coriolis, centrifugal, and gravitational terms.
///
/// # Arguments
///
/// * `state` — Current state in the synodic frame.
/// * `mu` — Mass ratio.
#[must_use]
pub fn equations_of_motion(state: &SynodicState, mu: MassRatio) -> SynodicState {
    let [x, y, z, vx, vy, _vz] = *state;

    let r1 = distance_to_m1(x, y, z, mu);
    let r2 = distance_to_m2(x, y, z, mu);
    let r1_3 = r1 * r1 * r1;
    let r2_3 = r2 * r2 * r2;

    // Positions of primaries in synodic frame
    // m₁ at (-μ, 0, 0), m₂ at (1-μ, 0, 0)
    let ax = 2.0 * vy + x - (1.0 - mu) * (x + mu) / r1_3 - mu * (x - 1.0 + mu) / r2_3;
    let ay = -2.0 * vx + y - (1.0 - mu) * y / r1_3 - mu * y / r2_3;
    let az = -(1.0 - mu) * z / r1_3 - mu * z / r2_3;

    [vx, vy, _vz, ax, ay, az]
}

/// Distance from point (x, y, z) to the larger primary m₁ at (−μ, 0, 0).
#[must_use]
#[inline]
fn distance_to_m1(x: f64, y: f64, z: f64, mu: MassRatio) -> f64 {
    ((x + mu) * (x + mu) + y * y + z * z).sqrt()
}

/// Distance from point (x, y, z) to the smaller primary m₂ at (1−μ, 0, 0).
#[must_use]
#[inline]
fn distance_to_m2(x: f64, y: f64, z: f64, mu: MassRatio) -> f64 {
    ((x - 1.0 + mu) * (x - 1.0 + mu) + y * y + z * z).sqrt()
}

// ── Collinear Lagrange point solvers ─────────────────────────────────────

/// Find L1: between m₁ and m₂.
///
/// Solves the quintic equation for the collinear equilibrium between the
/// two primaries using Newton's method.
fn find_collinear_l1(mu: MassRatio) -> f64 {
    // L1 is between m₁ and m₂, at x = 1 - μ - r₂
    // where r₂ is the distance from m₂ to L1.
    // Initial guess: Hill sphere approximation
    let mut gamma = (mu / 3.0).cbrt();

    for _ in 0..50 {
        let x = 1.0 - mu - gamma;
        let f = collinear_force(x, mu);
        let df = collinear_force_derivative(x, mu);

        if df.abs() < 1e-30 {
            break;
        }

        let dx = f / df;
        gamma += dx; // x decreases → γ increases

        if dx.abs() < 1e-14 {
            break;
        }
    }

    1.0 - mu - gamma
}

/// Find L2: beyond m₂ (away from m₁).
fn find_collinear_l2(mu: MassRatio) -> f64 {
    // L2 is beyond m₂, at x = 1 - μ + γ
    let mut gamma = (mu / 3.0).cbrt();

    for _ in 0..50 {
        let x = 1.0 - mu + gamma;
        let f = collinear_force(x, mu);
        let df = collinear_force_derivative(x, mu);

        if df.abs() < 1e-30 {
            break;
        }

        let dx = f / df;
        gamma -= dx;

        if dx.abs() < 1e-14 {
            break;
        }
    }

    1.0 - mu + gamma
}

/// Find L3: beyond m₁ (opposite side from m₂).
fn find_collinear_l3(mu: MassRatio) -> f64 {
    // L3 is on the opposite side of m₁ from m₂, near x ≈ -1
    let mut x = -1.0 + 5.0 * mu / 12.0; // initial guess

    for _ in 0..50 {
        let f = collinear_force(x, mu);
        let df = collinear_force_derivative(x, mu);

        if df.abs() < 1e-30 {
            break;
        }

        let dx = f / df;
        x -= dx;

        if dx.abs() < 1e-14 {
            break;
        }
    }

    x
}

/// Force balance on the x-axis: f(x) = 0 at equilibrium points.
///
/// f(x) = x - (1-μ)(x+μ)/|x+μ|³ - μ(x-1+μ)/|x-1+μ|³
#[inline]
fn collinear_force(x: f64, mu: MassRatio) -> f64 {
    let d1 = x + mu;
    let d2 = x - 1.0 + mu;
    let r1 = d1.abs();
    let r2 = d2.abs();

    if r1 < 1e-15 || r2 < 1e-15 {
        return 0.0;
    }

    x - (1.0 - mu) * d1 / (r1 * r1 * r1) - mu * d2 / (r2 * r2 * r2)
}

/// Derivative of the collinear force balance.
#[inline]
fn collinear_force_derivative(x: f64, mu: MassRatio) -> f64 {
    let d1 = x + mu;
    let d2 = x - 1.0 + mu;
    let r1 = d1.abs();
    let r2 = d2.abs();

    if r1 < 1e-15 || r2 < 1e-15 {
        return 1.0;
    }

    let r1_3 = r1 * r1 * r1;
    let r2_3 = r2 * r2 * r2;

    1.0 - (1.0 - mu) * (1.0 / r1_3 - 3.0 * d1 * d1 / (r1_3 * r1 * r1))
        - mu * (1.0 / r2_3 - 3.0 * d2 * d2 / (r2_3 * r2 * r2))
}

#[cfg(test)]
mod tests {
    use super::*;

    // Earth-Moon mass ratio
    const MU_EM: f64 = 0.012_150_585;

    // Sun-Earth mass ratio
    const MU_SE: f64 = 3.003_467e-6;

    #[test]
    fn lagrange_points_earth_moon() {
        let lp = lagrange_points(MU_EM).unwrap();

        // L1: between Earth and Moon, ~84% of the way from Earth to Moon
        assert!(
            lp.l1[0] > 0.8 && lp.l1[0] < 0.9,
            "L1 x={} should be ~0.837",
            lp.l1[0]
        );
        assert!(lp.l1[1].abs() < 1e-10, "L1 should be on x-axis");

        // L2: beyond Moon, ~1.16
        assert!(
            lp.l2[0] > 1.1 && lp.l2[0] < 1.2,
            "L2 x={} should be ~1.156",
            lp.l2[0]
        );

        // L3: opposite side, ~−1
        assert!(
            lp.l3[0] > -1.01 && lp.l3[0] < -0.99,
            "L3 x={} should be ~-1.005",
            lp.l3[0]
        );

        // L4, L5: equilateral points at (0.5-μ, ±√3/2)
        let expected_x = 0.5 - MU_EM;
        let expected_y = 3.0_f64.sqrt() / 2.0;
        assert!((lp.l4[0] - expected_x).abs() < 1e-10, "L4 x");
        assert!((lp.l4[1] - expected_y).abs() < 1e-10, "L4 y");
        assert!((lp.l5[0] - expected_x).abs() < 1e-10, "L5 x");
        assert!((lp.l5[1] + expected_y).abs() < 1e-10, "L5 y");
    }

    #[test]
    fn lagrange_points_sun_earth() {
        let lp = lagrange_points(MU_SE).unwrap();

        // L1: ~0.99 AU (Earth-Sun L1, approximately 1.5M km sunward)
        assert!(
            lp.l1[0] > 0.98 && lp.l1[0] < 1.0,
            "SE L1 x={} should be ~0.990",
            lp.l1[0]
        );

        // L2: ~1.01 AU
        assert!(
            lp.l2[0] > 1.0 && lp.l2[0] < 1.02,
            "SE L2 x={} should be ~1.010",
            lp.l2[0]
        );
    }

    #[test]
    fn lagrange_points_are_equilibria() {
        let lp = lagrange_points(MU_EM).unwrap();

        // At each Lagrange point, the acceleration should be zero (equilibrium)
        for (name, point) in [
            ("L1", lp.l1),
            ("L2", lp.l2),
            ("L3", lp.l3),
            ("L4", lp.l4),
            ("L5", lp.l5),
        ] {
            let state = [point[0], point[1], point[2], 0.0, 0.0, 0.0];
            let deriv = equations_of_motion(&state, MU_EM);
            // vx, vy, vz should be zero (they are)
            // ax, ay, az should be near zero at equilibrium
            let accel_mag =
                (deriv[3] * deriv[3] + deriv[4] * deriv[4] + deriv[5] * deriv[5]).sqrt();
            assert!(
                accel_mag < 1e-10,
                "{name}: acceleration magnitude {accel_mag:.2e} should be ~0"
            );
        }
    }

    #[test]
    fn jacobi_constant_conserved() {
        // At L1, the Jacobi constant should match the zero-velocity value
        let lp = lagrange_points(MU_EM).unwrap();
        let state_l1 = [lp.l1[0], lp.l1[1], lp.l1[2], 0.0, 0.0, 0.0];
        let cj = jacobi_constant(&state_l1, MU_EM);
        let zv = zero_velocity_value(lp.l1[0], lp.l1[1], lp.l1[2], MU_EM);
        assert!(
            (cj - zv).abs() < 1e-10,
            "At L1 with zero velocity, C_J={cj} should equal 2Ω={zv}"
        );
    }

    #[test]
    fn jacobi_ordering() {
        // C_J(L1) < C_J(L2) < C_J(L3) < C_J(L4) = C_J(L5) for μ < ~0.038
        let lp = lagrange_points(MU_EM).unwrap();

        let cj = |p: [f64; 3]| jacobi_constant(&[p[0], p[1], p[2], 0.0, 0.0, 0.0], MU_EM);

        let cj1 = cj(lp.l1);
        let cj2 = cj(lp.l2);
        let cj3 = cj(lp.l3);
        let cj4 = cj(lp.l4);
        let cj5 = cj(lp.l5);

        // Standard ordering for small μ: L1 > L2 > L3 > L4 = L5
        // (higher C_J = more restricted = lower energy)
        assert!(cj1 > cj2, "C_J(L1)={cj1} should > C_J(L2)={cj2}");
        assert!(cj2 > cj3, "C_J(L2)={cj2} should > C_J(L3)={cj3}");
        assert!(cj3 > cj4, "C_J(L3)={cj3} should > C_J(L4)={cj4}");
        assert!(
            (cj4 - cj5).abs() < 1e-10,
            "C_J(L4)={cj4} should ≈ C_J(L5)={cj5}"
        );
    }

    #[test]
    fn equations_of_motion_symmetry() {
        // CR3BP symmetry: (x, y, vx, vy) → (x, -y, -vx, vy) under time reversal
        // More useful: gravitational + centrifugal terms are symmetric in y.
        // Test: at the same x with y=0, ay should be zero (on the x-axis)
        let state = [0.5, 0.0, 0.0, 0.0, 0.3, 0.0];
        let d = equations_of_motion(&state, MU_EM);
        assert!(
            d[4].abs() < 1e-12,
            "ay on x-axis should be zero (centrifugal + gravity cancel): {}",
            d[4]
        );

        // Verify gravitational symmetry: positions mirrored in y give opposite ay
        // (with zero velocities to remove Coriolis)
        let s1 = [0.5, 0.3, 0.0, 0.0, 0.0, 0.0];
        let s2 = [0.5, -0.3, 0.0, 0.0, 0.0, 0.0];
        let d1 = equations_of_motion(&s1, MU_EM);
        let d2 = equations_of_motion(&s2, MU_EM);
        assert!(
            (d1[3] - d2[3]).abs() < 1e-12,
            "ax should be same: {} vs {}",
            d1[3],
            d2[3]
        );
        assert!(
            (d1[4] + d2[4]).abs() < 1e-12,
            "ay should be negated: {} vs {}",
            d1[4],
            d2[4]
        );
    }

    #[test]
    fn pseudo_potential_at_primaries() {
        // At m₂ (x=1-μ), potential should be very large (μ/r₂ → ∞)
        // At distance 0.001, μ/r₂ ≈ 0.012/0.001 ≈ 12
        let val = pseudo_potential(1.0 - MU_EM + 0.001, 0.0, 0.0, MU_EM);
        assert!(val > 10.0, "potential near m₂ should be large: {val}");

        // Closer should be even larger
        let val_closer = pseudo_potential(1.0 - MU_EM + 0.0001, 0.0, 0.0, MU_EM);
        assert!(
            val_closer > val,
            "closer to m₂ should increase potential: {val_closer} vs {val}"
        );
    }

    #[test]
    fn lagrange_invalid_mu() {
        assert!(lagrange_points(0.0).is_err());
        assert!(lagrange_points(-0.1).is_err());
        assert!(lagrange_points(0.6).is_err());
    }
}
