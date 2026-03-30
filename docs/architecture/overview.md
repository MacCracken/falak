# Architecture

## Module Map

| Module | Files | Key Types | Purpose |
|--------|-------|-----------|---------|
| `error` | `error.rs` | `FalakError`, `Result<T>` | Shared error types with `Cow<'static, str>` messages |
| `orbit` | `orbit.rs` | `OrbitalElements` | Classical Keplerian orbital elements, state vectors |
| `kepler` | `kepler.rs` | (planned) | Kepler's equation, anomaly conversions, orbital period/velocity |
| `transfer` | `transfer.rs` | (planned) | Hohmann, bi-elliptic, Lambert transfer maneuvers |
| `perturbation` | `perturbation.rs` | (planned) | J2, drag, SRP, third-body perturbation models |
| `nbody` | `nbody.rs` | (planned) | N-body gravitational simulation, symplectic integrators |
| `ephemeris` | `ephemeris.rs` | (planned) | Planetary positions, Julian date, VSOP87 |
| `frame` | `frame.rs` | (planned) | ECI, ECEF, perifocal, rotating frame transforms |
| `maneuver` | `maneuver.rs` | (planned) | Delta-v budgets, impulsive/continuous thrust |
| `logging` | `logging.rs` | (init functions) | tracing-subscriber setup |

## Design Principles

- Flat library crate -- no internal binaries
- Feature-gated optional modules -- consumers pull only what they need
- f64 precision throughout -- astrodynamics demands double precision
- `#[must_use]` on all pure functions, `#[non_exhaustive]` on all public enums
- `#![warn(missing_docs)]` and `#![forbid(unsafe_code)]`
- `Cow<'static, str>` in error variants -- zero alloc on static messages
- `#[inline]` on hot-path functions
- tracing instrumentation on multi-step computations

## Data Flow

```
orbit ──> orbital elements (Keplerian state definition)
  │
  v
kepler ──> anomaly conversions, period, velocity (Kepler's equation solver)
  │
  v
transfer ──> maneuver planning (Hohmann, Lambert, bi-elliptic)
  │
  v
maneuver ──> delta-v budgets, burn execution
  │
  v
perturbation ──> real-world corrections (J2, drag, SRP)
  │
  v
nbody ──> multi-body simulation (numerical integration)

ephemeris ──> planetary positions, time systems
frame ──> coordinate transforms between reference frames
```

## Dependency Stack

```
falak
  └── hisab (linear algebra, geometry, calculus, spatial structures)
  └── impetus (optional: physics engine for force/torque integration)
```

## Consumers

- **kiran** — game engine: orbital mechanics for space simulations
- **joshua** — simulation: celestial body dynamics, mission planning
