//! Soorat integration — visualization data structures for orbital mechanics.
//!
//! Provides structured types that soorat can render: orbit paths,
//! planetary positions, transfer trajectories, and ground tracks.

use serde::{Deserialize, Serialize};

// ── Orbit path ─────────────────────────────────────────────────────────────

/// Orbital trajectory points for line rendering.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[non_exhaustive]
pub struct OrbitPath {
    /// Points along the orbit `[x, y, z]` in metres (or km, consumer-defined).
    pub points: Vec<[f64; 3]>,
    /// Velocity magnitude at each point (for color mapping).
    pub speeds: Vec<f64>,
    /// Orbital period (seconds).
    pub period: f64,
    /// Semi-major axis (same unit as points).
    pub semi_major_axis: f64,
    /// Eccentricity.
    pub eccentricity: f64,
}

impl OrbitPath {
    /// Generate orbit path from orbital elements and gravitational parameter.
    ///
    /// `mu`: gravitational parameter (m³/s²).
    /// `num_points`: number of points around the orbit.
    ///
    /// # Errors
    ///
    /// Returns [`crate::FalakError::InvalidParameter`] if orbital parameters or μ are invalid.
    #[must_use = "returns the computed orbit path"]
    pub fn from_elements(
        elements: &crate::orbit::OrbitalElements,
        mu: f64,
        num_points: usize,
    ) -> crate::error::Result<Self> {
        let a = elements.semi_major_axis;
        let e = elements.eccentricity;
        let period = crate::kepler::orbital_period(a, mu)?;

        let mut points = Vec::with_capacity(num_points);
        let mut speeds = Vec::with_capacity(num_points);

        for i in 0..num_points {
            let true_anom = std::f64::consts::TAU * i as f64 / num_points as f64;
            let r = crate::kepler::radius_at_true_anomaly(a, e, true_anom);

            // Perifocal coordinates
            let px = r * true_anom.cos();
            let py = r * true_anom.sin();

            // Simplified: perifocal frame (skip full rotation to ECI)
            points.push([px, py, 0.0]);

            // Vis-viva speed
            let v = crate::kepler::vis_viva(r, a, mu)?;
            speeds.push(v);
        }

        Ok(Self {
            points,
            speeds,
            period,
            semi_major_axis: a,
            eccentricity: e,
        })
    }

    /// Generate orbit path in the ECI frame from orbital elements.
    ///
    /// Unlike [`from_elements`](Self::from_elements), this applies the full
    /// perifocal → ECI rotation using RAAN, inclination, and argument of periapsis.
    ///
    /// # Errors
    ///
    /// Returns [`crate::FalakError::InvalidParameter`] if orbital parameters or μ are invalid.
    #[must_use = "returns the computed orbit path"]
    pub fn from_elements_eci(
        elements: &crate::orbit::OrbitalElements,
        mu: f64,
        num_points: usize,
    ) -> crate::error::Result<Self> {
        let a = elements.semi_major_axis;
        let e = elements.eccentricity;
        let period = crate::kepler::orbital_period(a, mu)?;

        let mut points = Vec::with_capacity(num_points);
        let mut speeds = Vec::with_capacity(num_points);

        for i in 0..num_points {
            let true_anom = std::f64::consts::TAU * i as f64 / num_points as f64;
            let r = crate::kepler::radius_at_true_anomaly(a, e, true_anom);

            // Perifocal coordinates
            let pqw = [r * true_anom.cos(), r * true_anom.sin(), 0.0];

            // Rotate to ECI
            let eci = crate::frame::perifocal_to_eci(
                pqw,
                elements.raan,
                elements.inclination,
                elements.argument_of_periapsis,
            );
            points.push(eci);

            let v = crate::kepler::vis_viva(r, a, mu)?;
            speeds.push(v);
        }

        Ok(Self {
            points,
            speeds,
            period,
            semi_major_axis: a,
            eccentricity: e,
        })
    }
}

// ── Planetary positions ────────────────────────────────────────────────────

/// Body positions for instanced sphere rendering.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[non_exhaustive]
pub struct PlanetaryPositions {
    /// Bodies with position and properties.
    pub bodies: Vec<CelestialBody>,
}

/// A celestial body for rendering.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[non_exhaustive]
pub struct CelestialBody {
    /// Name.
    pub name: String,
    /// Position `[x, y, z]` in AU or metres.
    pub position: [f64; 3],
    /// Radius for display scaling (in same units as position).
    pub radius: f64,
    /// Mass (kg) for label/info.
    pub mass: f64,
}

// ── Transfer trajectory ────────────────────────────────────────────────────

/// Transfer trajectory for colored line rendering.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[non_exhaustive]
pub struct TransferTrajectory {
    /// Points along the transfer `[x, y, z]`.
    pub points: Vec<[f64; 3]>,
    /// Delta-v markers: `(point_index, delta_v_ms)`.
    pub burns: Vec<(usize, f64)>,
    /// Total delta-v (m/s).
    pub total_delta_v: f64,
    /// Transfer time (seconds).
    pub transfer_time: f64,
}

impl TransferTrajectory {
    /// Generate a transfer trajectory from a Lambert solution.
    ///
    /// Propagates the departure state (r1, v1) numerically using Cowell's
    /// method to produce renderable trajectory points.
    ///
    /// # Arguments
    ///
    /// * `r1` — Departure position (metres, ECI)
    /// * `r2` — Arrival position (metres, ECI)
    /// * `solution` — Lambert solution containing departure/arrival velocities
    /// * `tof` — Time of flight (seconds)
    /// * `mu` — Gravitational parameter (m³/s²)
    /// * `num_points` — Number of trajectory points to generate
    ///
    /// # Errors
    ///
    /// Returns an error if propagation fails.
    #[must_use = "returns the computed transfer trajectory"]
    pub fn from_lambert(
        r1: [f64; 3],
        r2: [f64; 3],
        solution: &crate::transfer::LambertSolution,
        tof: f64,
        mu: f64,
        num_points: usize,
    ) -> crate::error::Result<Self> {
        let num_points = num_points.max(2);
        let dt_seg = tof / (num_points - 1) as f64;
        let dt_step = dt_seg.min(10.0); // integration step ≤ 10s

        let mut points = Vec::with_capacity(num_points);
        let mut state = crate::kepler::StateVector {
            position: r1,
            velocity: solution.v1,
        };

        points.push(state.position);

        for i in 1..num_points {
            let seg_time = dt_seg.min(tof - (i - 1) as f64 * dt_seg);
            state = crate::propagate::two_body(&state, mu, seg_time, dt_step)?;
            points.push(state.position);
        }

        // Compute departure and arrival delta-v magnitudes
        // (departure relative to circular orbit at r1, arrival relative to circular at r2)
        let r1_mag = (r1[0] * r1[0] + r1[1] * r1[1] + r1[2] * r1[2]).sqrt();
        let r2_mag = (r2[0] * r2[0] + r2[1] * r2[1] + r2[2] * r2[2]).sqrt();
        let v_circ1 = (mu / r1_mag).sqrt();
        let v_circ2 = (mu / r2_mag).sqrt();
        let v1_mag =
            (solution.v1[0].powi(2) + solution.v1[1].powi(2) + solution.v1[2].powi(2)).sqrt();
        let v2_mag =
            (solution.v2[0].powi(2) + solution.v2[1].powi(2) + solution.v2[2].powi(2)).sqrt();
        let dv1 = (v1_mag - v_circ1).abs();
        let dv2 = (v2_mag - v_circ2).abs();

        Ok(Self {
            points,
            burns: vec![(0, dv1), (num_points - 1, dv2)],
            total_delta_v: dv1 + dv2,
            transfer_time: tof,
        })
    }
}

// ── Ground track ───────────────────────────────────────────────────────────

/// Ground track for map overlay rendering.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[non_exhaustive]
pub struct GroundTrack {
    /// Sub-satellite points `[latitude_deg, longitude_deg]`.
    pub points: Vec<[f64; 2]>,
    /// Altitude at each point (km).
    pub altitudes: Vec<f64>,
    /// Orbit number at each point (for color coding successive orbits).
    pub orbit_numbers: Vec<u32>,
}

impl GroundTrack {
    /// Compute a ground track from orbital elements over multiple orbits.
    ///
    /// Propagates the orbit analytically, converts each position to ECEF
    /// (using GMST at the epoch + elapsed time), then to geodetic coordinates.
    ///
    /// # Arguments
    ///
    /// * `elements` — Initial orbital elements (elliptical).
    /// * `mu` — Gravitational parameter (m³/s²).
    /// * `epoch_jd` — Julian Date of the epoch (for GMST computation).
    /// * `num_orbits` — Number of orbits to trace.
    /// * `points_per_orbit` — Number of sample points per orbit.
    ///
    /// # Errors
    ///
    /// Returns errors if elements are not elliptical or parameters are invalid.
    #[must_use = "returns the computed ground track"]
    pub fn from_elements(
        elements: &crate::orbit::OrbitalElements,
        mu: f64,
        epoch_jd: f64,
        num_orbits: u32,
        points_per_orbit: usize,
    ) -> crate::error::Result<Self> {
        let period = crate::kepler::orbital_period(elements.semi_major_axis, mu)?;
        let total_points = num_orbits as usize * points_per_orbit;
        let total_time = period * num_orbits as f64;
        let dt = total_time / total_points as f64;

        let gmst0 = crate::ephemeris::gmst(epoch_jd);
        let earth_rate = crate::ephemeris::EARTH_ROTATION_RATE;

        let mut points = Vec::with_capacity(total_points);
        let mut altitudes = Vec::with_capacity(total_points);
        let mut orbit_numbers = Vec::with_capacity(total_points);

        for i in 0..total_points {
            let t = dt * i as f64;
            let orbit_num = (t / period) as u32;

            // Propagate orbit to time t
            let prop = crate::propagate::kepler(elements, mu, t)?;
            let state = crate::kepler::elements_to_state(&prop, mu)?;

            // GMST at this time
            let gmst_t = gmst0 + earth_rate * t;

            // ECI → ECEF → geodetic
            let ecef = crate::frame::eci_to_ecef(state.position, gmst_t);
            let geo = crate::frame::ecef_to_geodetic(ecef)?;

            points.push([geo.latitude.to_degrees(), geo.longitude.to_degrees()]);
            altitudes.push(geo.altitude / 1000.0); // metres → km
            orbit_numbers.push(orbit_num);
        }

        Ok(Self {
            points,
            altitudes,
            orbit_numbers,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn orbit_path_circular() {
        let elements = crate::orbit::OrbitalElements::new(7e6, 0.0, 0.0, 0.0, 0.0, 0.0).unwrap();
        let path = OrbitPath::from_elements(&elements, 3.986e14, 36).unwrap();
        assert_eq!(path.points.len(), 36);
        assert_eq!(path.speeds.len(), 36);
        assert!(path.period > 0.0);
        // Circular orbit → all speeds roughly equal
        let v_avg = path.speeds.iter().sum::<f64>() / path.speeds.len() as f64;
        for v in &path.speeds {
            assert!((v - v_avg).abs() / v_avg < 0.01, "speed variance too high");
        }
    }

    #[test]
    fn orbit_path_elliptical() {
        let elements = crate::orbit::OrbitalElements::new(10e6, 0.5, 0.0, 0.0, 0.0, 0.0).unwrap();
        let path = OrbitPath::from_elements(&elements, 3.986e14, 72).unwrap();
        assert_eq!(path.points.len(), 72);
        assert!((path.eccentricity - 0.5).abs() < 0.001);
        // Elliptical → speed varies
        let v_min = path.speeds.iter().cloned().fold(f64::MAX, f64::min);
        let v_max = path.speeds.iter().cloned().fold(f64::MIN, f64::max);
        assert!(
            v_max > v_min * 1.5,
            "elliptical should have speed variation"
        );
    }

    #[test]
    fn orbit_path_eci_inclined() {
        let elements = crate::orbit::OrbitalElements::new(7e6, 0.0, 0.5, 1.0, 0.5, 0.0).unwrap();
        let path = OrbitPath::from_elements_eci(&elements, 3.986e14, 36).unwrap();
        assert_eq!(path.points.len(), 36);
        // Inclined orbit in ECI should have non-zero z components
        let has_z = path.points.iter().any(|p| p[2].abs() > 1e3);
        assert!(has_z, "ECI inclined orbit should have z-components");
    }

    #[test]
    fn orbit_path_eci_equatorial_matches_perifocal() {
        // For equatorial, zero-arg orbit: ECI = perifocal
        let elements = crate::orbit::OrbitalElements::new(7e6, 0.0, 0.0, 0.0, 0.0, 0.0).unwrap();
        let pf = OrbitPath::from_elements(&elements, 3.986e14, 12).unwrap();
        let eci = OrbitPath::from_elements_eci(&elements, 3.986e14, 12).unwrap();
        for (p, e) in pf.points.iter().zip(eci.points.iter()) {
            for k in 0..3 {
                assert!(
                    (p[k] - e[k]).abs() < 1.0,
                    "equatorial perifocal vs ECI mismatch"
                );
            }
        }
    }

    #[test]
    fn ground_track_iss_like() {
        // ISS-like orbit: ~400km altitude, 51.6° inclination
        let elements =
            crate::orbit::OrbitalElements::new(6.771e6, 0.0005, 0.9, 0.0, 0.0, 0.0).unwrap();
        let jd = crate::ephemeris::J2000_JD;
        let gt = GroundTrack::from_elements(&elements, 3.986e14, jd, 1, 36).unwrap();

        assert_eq!(gt.points.len(), 36);
        assert_eq!(gt.altitudes.len(), 36);
        assert_eq!(gt.orbit_numbers.len(), 36);

        // Latitude should be within ±inclination (in degrees, ~51.6°)
        for pt in &gt.points {
            assert!(pt[0].abs() < 60.0, "latitude out of range: {}°", pt[0]);
        }

        // Altitude should be ~400 km
        for alt in &gt.altitudes {
            assert!((*alt - 400.0).abs() < 50.0, "altitude: {} km", alt);
        }

        // All orbit 0
        assert!(gt.orbit_numbers.iter().all(|&n| n == 0));
    }

    #[test]
    fn ground_track_two_orbits() {
        let elements = crate::orbit::OrbitalElements::new(7e6, 0.0, 0.5, 0.0, 0.0, 0.0).unwrap();
        let gt = GroundTrack::from_elements(&elements, 3.986e14, crate::ephemeris::J2000_JD, 2, 18)
            .unwrap();
        assert_eq!(gt.points.len(), 36); // 2 × 18
        // Should have both orbit 0 and orbit 1
        assert!(gt.orbit_numbers.contains(&0));
        assert!(gt.orbit_numbers.contains(&1));
    }

    #[test]
    fn planetary_positions_serializes() {
        let pp = PlanetaryPositions {
            bodies: vec![CelestialBody {
                name: "Earth".into(),
                position: [1.0, 0.0, 0.0],
                radius: 6.371e6,
                mass: 5.972e24,
            }],
        };
        let json = serde_json::to_string(&pp);
        assert!(json.is_ok());
    }

    #[test]
    fn transfer_trajectory_manual() {
        let tt = TransferTrajectory {
            points: vec![[7e6, 0.0, 0.0], [0.0, 4.2e7, 0.0]],
            burns: vec![(0, 2450.0), (1, 1680.0)],
            total_delta_v: 4130.0,
            transfer_time: 15768000.0,
        };
        assert_eq!(tt.burns.len(), 2);
    }

    #[test]
    fn ground_track_manual() {
        let gt = GroundTrack {
            points: vec![[28.5, -80.6], [30.0, -75.0]],
            altitudes: vec![400.0, 400.0],
            orbit_numbers: vec![1, 1],
        };
        assert_eq!(gt.points.len(), 2);
    }

    #[test]
    fn transfer_trajectory_from_lambert() {
        let mu = 3.986_004_418e14;
        let r1 = [7e6, 0.0, 0.0];
        let r2 = [0.0, 7e6, 0.0]; // 90° transfer
        let period = crate::kepler::orbital_period(7e6, mu).unwrap();
        let tof = period / 4.0;

        let sol = crate::transfer::lambert(r1, r2, tof, mu, true).unwrap();
        let tt = TransferTrajectory::from_lambert(r1, r2, &sol, tof, mu, 50).unwrap();

        assert_eq!(tt.points.len(), 50);
        assert_eq!(tt.burns.len(), 2);
        assert!(tt.total_delta_v > 0.0);
        assert!((tt.transfer_time - tof).abs() < 1e-6);

        // First point should be near r1
        let p0 = tt.points[0];
        assert!((p0[0] - r1[0]).abs() < 1.0);

        // Last point should be near r2
        let pn = tt.points[49];
        let dist_to_r2 =
            ((pn[0] - r2[0]).powi(2) + (pn[1] - r2[1]).powi(2) + (pn[2] - r2[2]).powi(2)).sqrt();
        assert!(
            dist_to_r2 < 10_000.0, // within 10 km
            "last point should be near r2, dist={dist_to_r2:.0} m"
        );
    }
}
