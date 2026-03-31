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
    #[must_use]
    pub fn from_elements(
        elements: &crate::orbit::OrbitalElements,
        mu: f64,
        num_points: usize,
    ) -> Self {
        let a = elements.semi_major_axis;
        let e = elements.eccentricity;
        // Kepler's third law via kepler module
        let period = crate::kepler::orbital_period(a, mu).unwrap_or(0.0);

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
            let v = crate::kepler::vis_viva(r, a, mu).unwrap_or(0.0);
            speeds.push(v);
        }

        Self {
            points,
            speeds,
            period,
            semi_major_axis: a,
            eccentricity: e,
        }
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn orbit_path_circular() {
        let elements = crate::orbit::OrbitalElements::new(7e6, 0.0, 0.0, 0.0, 0.0, 0.0).unwrap();
        let path = OrbitPath::from_elements(&elements, 3.986e14, 36);
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
        let path = OrbitPath::from_elements(&elements, 3.986e14, 72);
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
}
