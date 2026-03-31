#![warn(missing_docs)]
#![forbid(unsafe_code)]

//! Falak — Orbital mechanics and celestial dynamics for AGNOS
//!
//! Arabic/Persian: فلک (falak) — sky, celestial sphere
//!
//! Provides Keplerian orbit mechanics, transfer maneuvers, perturbation
//! models, N-body simulation, ephemeris computation, reference frames,
//! and spacecraft maneuver planning. Built on
//! [hisab](https://crates.io/crates/hisab) for math foundations.
//!
//! # Modules
//!
//! - [`error`] — Error types for orbital computations
//! - [`orbit`] — Orbital elements and state vectors
//! - [`kepler`] — Kepler's laws, anomaly conversions, orbital period and velocity
//! - [`transfer`] — Hohmann, bi-elliptic, and Lambert transfer maneuvers
//! - [`perturbation`] — J2 oblateness, drag, solar radiation pressure, third-body effects
//! - [`nbody`] — N-body gravitational simulation with numerical integration
//! - [`ephemeris`] — Planetary positions, Julian date, coordinate transforms
//! - [`frame`] — Reference frames (ECI, ECEF, perifocal, rotating)
//! - [`maneuver`] — Delta-v budgets, impulsive burns, continuous thrust

/// Cross-crate bridges — primitive-value conversions from other AGNOS science crates.
pub mod bridge;
pub mod error;
/// Integration APIs for downstream consumers (soorat rendering).
pub mod integration;

/// Orbital elements and state vectors.
pub mod orbit;

/// Kepler's laws — period, velocity, anomaly conversions (mean, eccentric, true).
pub mod kepler;

/// Orbital transfer maneuvers — Hohmann, bi-elliptic, Lambert problem.
pub mod transfer;

/// Orbital perturbations — J2 oblateness, drag, solar radiation pressure, third body.
pub mod perturbation;

/// N-body gravitational simulation — numerical integration, symplectic methods.
pub mod nbody;

/// Ephemeris computation — planetary positions, Julian date, coordinate transforms.
pub mod ephemeris;

/// Reference frames — ECI, ECEF, perifocal, rotating frames, coordinate transforms.
pub mod frame;

/// Spacecraft maneuvers — delta-v budgets, impulsive burns, continuous thrust.
pub mod maneuver;

#[cfg(feature = "logging")]
/// Structured logging for falak via `FALAK_LOG` env var.
pub mod logging;

pub use error::FalakError;
