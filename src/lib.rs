#![warn(missing_docs)]
#![forbid(unsafe_code)]

//! Falak — Orbital Mechanics & Celestial Dynamics
//!
//! Arabic/Persian: فلک (falak) — sky, celestial sphere
//!
//! # Architecture
//!
//! Falak provides the mechanics of orbits for the AGNOS science stack.
//! Built on [hisab](https://crates.io/crates/hisab) for math foundations,
//! with optional [impetus](https://crates.io/crates/impetus) integration for
//! physics coupling.
//!
//! ```text
//! hisab (math) ──┐
//!                ├── falak (orbits) ──┬── kiran (game engine)
//! impetus (phys)─┘                   ├── joshua (simulation)
//!                                    └── soorat (rendering)
//! ```
//!
//! ## Domain Modules
//!
//! | Module | Feature | Description |
//! |--------|---------|-------------|
//! | [`orbit`] | always | Classical Keplerian orbital elements |
//! | [`kepler`] | always | Kepler's equation, anomaly conversions, state vectors |
//! | [`transfer`] | always | Hohmann, bi-elliptic, Lambert, combined maneuvers |
//! | [`perturbation`] | always | J2/J3, drag, SRP, third-body accelerations |
//! | [`nbody`] | always | N-body simulation: direct O(N²), Barnes-Hut O(N log N) |
//! | [`ephemeris`] | always | Julian dates, GMST, planets, Moon, eclipses, rise/set |
//! | [`frame`] | always | ECI/ECEF/perifocal/geodetic, precession, nutation |
//! | [`maneuver`] | always | Burns, rocket equation, escape/capture, low-thrust |
//! | [`propagate`] | always | Kepler/Cowell/Encke propagation, mean elements |
//! | [`cr3bp`] | always | Lagrange points, Jacobi constant, CR3BP dynamics |
//! | [`bridge`] | always | Cross-crate conversions (tara/impetus/badal) |
//! | [`integration`] | always | Consumer APIs (soorat orbit paths, ground tracks) |
//!
//! # Examples
//!
//! ## Compute an orbital period
//!
//! ```
//! let mu_earth = 3.986_004_418e14; // m³/s²
//! let sma = 6_778e3; // ISS altitude
//! let period = falak::kepler::orbital_period(sma, mu_earth).unwrap();
//! assert!((period / 60.0 - 92.6).abs() < 0.5); // ~92.6 minutes
//! ```
//!
//! ## Hohmann transfer LEO → GEO
//!
//! ```
//! let mu = 3.986_004_418e14;
//! let h = falak::transfer::hohmann(6_778e3, 42_164e3, mu).unwrap();
//! assert!((h.total_delta_v - 3854.0).abs() < 100.0); // ~3.85 km/s
//! ```
//!
//! ## Earth-Moon Lagrange points
//!
//! ```
//! let lp = falak::cr3bp::lagrange_points(0.012_150_585).unwrap();
//! assert!(lp.l1[0] > 0.8 && lp.l1[0] < 0.9); // L1 between Earth and Moon
//! ```

/// Cross-crate bridges — primitive-value conversions from other AGNOS science crates.
pub mod bridge;
/// Error types and input validation for orbital computations.
pub mod error;
/// Integration APIs for downstream consumers (soorat rendering).
pub mod integration;

/// Classical Keplerian orbital elements.
pub mod orbit;

/// Kepler's laws — period, velocity, anomaly conversions (mean, eccentric, true),
/// state vector ↔ orbital elements conversion.
pub mod kepler;

/// Orbital transfer maneuvers — Hohmann, bi-elliptic, Lambert, plane change,
/// phasing, combined altitude + plane change.
pub mod transfer;

/// Orbital perturbations — J2/J3 oblateness, atmospheric drag, solar radiation
/// pressure, third-body gravitational effects.
pub mod perturbation;

/// N-body gravitational simulation — direct O(N²) summation, Barnes-Hut O(N log N),
/// leapfrog (symplectic), RK4, adaptive RK45 integrators.
pub mod nbody;

/// Ephemeris — Julian dates, sidereal time, planetary/lunar positions,
/// eclipse prediction, rise/set/transit times.
pub mod ephemeris;

/// Reference frames — ECI, ECEF, perifocal, rotating, geodetic transforms,
/// IAU 2006 precession, IAU 1980 nutation.
pub mod frame;

/// Spacecraft maneuvers — impulsive burns, rocket equation, escape/capture,
/// low-thrust spirals, Edelbaum approximation.
pub mod maneuver;

/// Orbit propagation — analytic two-body (Kepler), perturbed (Cowell/Encke),
/// Brouwer mean elements with secular J2 rates.
pub mod propagate;

/// Circular Restricted Three-Body Problem — Lagrange points, Jacobi constant,
/// zero-velocity curves, equations of motion in the synodic frame.
pub mod cr3bp;

/// Structured logging via `FALAK_LOG` env var (feature-gated).
#[cfg(feature = "logging")]
pub mod logging;

pub use error::FalakError;
