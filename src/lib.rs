#![warn(missing_docs)]
#![forbid(unsafe_code)]

//! Falak — Orbital mechanics and celestial dynamics for AGNOS
//!
//! Arabic/Persian: فلک (falak) — sky, celestial sphere
//!
//! Provides Keplerian orbit mechanics with anomaly conversions, state vector
//! transforms, and cross-crate bridges for the AGNOS ecosystem. Built on
//! [hisab](https://crates.io/crates/hisab) for math foundations.
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
//!
//! ### Stubs (not yet implemented)
//! - [`perturbation`], [`nbody`], [`ephemeris`], [`frame`], [`maneuver`]

/// Cross-crate bridges — primitive-value conversions from other AGNOS science crates.
pub mod bridge;
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

/// Orbital perturbations (stub — not yet implemented).
pub mod perturbation;

/// N-body gravitational simulation (stub — not yet implemented).
pub mod nbody;

/// Ephemeris computation (stub — not yet implemented).
pub mod ephemeris;

/// Reference frames (stub — not yet implemented).
pub mod frame;

/// Spacecraft maneuvers (stub — not yet implemented).
pub mod maneuver;

#[cfg(feature = "logging")]
/// Structured logging for falak via `FALAK_LOG` env var.
pub mod logging;

pub use error::FalakError;
