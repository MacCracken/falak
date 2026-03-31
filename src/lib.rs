#![warn(missing_docs)]
#![forbid(unsafe_code)]

//! Falak — Orbital mechanics and celestial dynamics for AGNOS
//!
//! Arabic/Persian: فلک (falak) — sky, celestial sphere
//!
//! Provides Keplerian orbit mechanics, transfer maneuvers, perturbation models,
//! N-body simulation, ephemeris computation, reference frames, and spacecraft
//! maneuver planning. Built on [hisab](https://crates.io/crates/hisab) for math
//! foundations.
//!
//! # Modules
//!
//! - [`error`] — Error types for orbital computations
//! - [`orbit`] — Orbital elements and state vectors
//! - [`kepler`] — Kepler's equation solver, anomaly conversions, period, velocity, state vectors
//! - [`bridge`] — Cross-crate conversions (tara/impetus/badal)
//! - [`integration`] — Downstream consumer APIs (soorat rendering)
//!
//! - [`transfer`] — Hohmann, bi-elliptic, plane change, phasing maneuvers
//! - [`frame`] — Reference frames (ECI, ECEF, perifocal, rotating, geodetic)
//! - [`ephemeris`] — Julian date, sidereal time, calendar conversions
//!
//! - [`perturbation`] — J2/J3, drag, SRP, third-body perturbation accelerations
//! - [`maneuver`] — Impulsive burns, rocket equation, escape/capture delta-v
//!
//! - [`nbody`] — N-body gravitational simulation with leapfrog and RK4 integrators
//! - [`propagate`] — Orbit propagation: analytic two-body and perturbed Cowell's method

/// Cross-crate bridges — primitive-value conversions from other AGNOS science crates.
pub mod bridge;
/// Error types for orbital computations.
pub mod error;
/// Integration APIs for downstream consumers (soorat rendering).
pub mod integration;

/// Orbital elements and state vectors.
pub mod orbit;

/// Kepler's laws — period, velocity, anomaly conversions (mean, eccentric, true),
/// state vector ↔ orbital elements conversion.
pub mod kepler;

/// Orbital transfer maneuvers — Hohmann, bi-elliptic, plane change, phasing.
pub mod transfer;

/// Orbital perturbations — J2/J3 oblateness, drag, SRP, third-body.
pub mod perturbation;

/// N-body gravitational simulation — direct summation, leapfrog, RK4.
pub mod nbody;

/// Ephemeris computation — Julian date, sidereal time, calendar conversions.
pub mod ephemeris;

/// Reference frames — ECI, ECEF, perifocal, rotating, geodetic.
pub mod frame;

/// Spacecraft maneuvers — impulsive burns, rocket equation, escape/capture.
pub mod maneuver;

/// Orbit propagation — two-body (Kepler) and perturbed (Cowell's method).
pub mod propagate;

/// Structured logging for falak via `FALAK_LOG` env var.
#[cfg(feature = "logging")]
pub mod logging;

pub use error::FalakError;
