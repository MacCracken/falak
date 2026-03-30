# Falak

> **Falak** (Arabic/Persian: فلک -- sky/celestial sphere) -- orbital mechanics engine for the AGNOS science stack

Keplerian orbits, transfer maneuvers, perturbation models, N-body simulation, ephemeris computation, reference frames, and spacecraft maneuver planning. Built on [hisab](https://crates.io/crates/hisab) for math foundations.

## Modules

| Module | Description |
|--------|-------------|
| `orbit` | Classical orbital elements, state vectors, periapsis/apoapsis |
| `kepler` | Kepler's laws, orbital period, velocity, anomaly conversions (mean, eccentric, true) |
| `transfer` | Hohmann transfers, bi-elliptic transfers, Lambert problem |
| `perturbation` | J2 oblateness, atmospheric drag, solar radiation pressure, third-body effects |
| `nbody` | N-body gravitational simulation, numerical integration, symplectic methods |
| `ephemeris` | Planetary positions, Julian date conversions, coordinate transforms |
| `frame` | Reference frames (ECI, ECEF, perifocal, rotating), coordinate transforms |
| `maneuver` | Delta-v budgets, impulsive burns, continuous thrust profiles |
| `error` | Unified error types for all modules |

## Feature Flags

| Feature | Default | Description |
|---------|---------|-------------|
| `logging` | no | Structured logging via `FALAK_LOG` env var |
| `physics` | no | Physics engine integration via impetus |

## Quick Start

```rust
use falak::orbit::OrbitalElements;

// Low Earth Orbit: 400km altitude, nearly circular, 51.6 degree inclination (ISS-like)
let iss_orbit = OrbitalElements::new(
    6_778_000.0,                  // semi-major axis (m)
    0.0001,                       // eccentricity (nearly circular)
    0.9006,                       // inclination (51.6 degrees in radians)
    1.2,                          // RAAN
    0.5,                          // argument of periapsis
    0.0,                          // true anomaly
).unwrap();

assert!(iss_orbit.periapsis() < iss_orbit.apoapsis());
```

## Dependency Stack

```
falak
  └── hisab (linear algebra, geometry, calculus)
  └── impetus (optional: physics engine integration)
```

## Consumers

- **kiran** — game engine (orbital mechanics for space simulations)
- **joshua** — simulation manager (celestial body dynamics)

## License

GPL-3.0-only. See [LICENSE](LICENSE).
