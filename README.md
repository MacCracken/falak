# Falak (فلک) — Orbital Mechanics & Celestial Dynamics

> **Falak** (Arabic/Persian: فلک — sky, celestial sphere) — orbital mechanics and celestial dynamics for the AGNOS science stack

Keplerian orbits, transfer maneuvers, perturbation models, N-body simulation, ephemeris computation, reference frames, spacecraft maneuver planning, orbit propagation, and the circular restricted three-body problem. Built on [hisab](https://crates.io/crates/hisab) for math foundations.

[![License: GPL-3.0-only](https://img.shields.io/badge/license-GPL--3.0--only-blue.svg)](LICENSE)

## Modules

| Module | Description |
|--------|-------------|
| `orbit` | Classical orbital elements, state vectors, periapsis/apoapsis, Keplerian parameters |
| `kepler` | Kepler's equation solver, anomaly conversions (mean, eccentric, true), period, velocity, state vectors |
| `transfer` | Hohmann transfers, bi-elliptic transfers, plane changes, phasing orbits, Lambert problem |
| `perturbation` | J2/J3 oblateness, atmospheric drag, solar radiation pressure, third-body effects |
| `nbody` | N-body gravitational simulation, direct summation, leapfrog and RK4 integrators |
| `ephemeris` | Julian date conversions, sidereal time, planetary and lunar positions |
| `frame` | Reference frames (ECI, ECEF, perifocal, rotating, geodetic), coordinate transforms |
| `maneuver` | Delta-v budgets, impulsive burns, rocket equation, escape/capture maneuvers |
| `propagate` | Orbit propagation: analytic two-body (Kepler) and perturbed Cowell's method |
| `cr3bp` | Circular Restricted Three-Body Problem: Lagrange points, Jacobi constant, zero-velocity curves |
| `bridge` | Cross-crate conversions — primitive-value bridges to other AGNOS science crates |
| `integration::soorat` | Downstream consumer API for soorat rendering integration |
| `error` | Unified error types (`FalakError`) for all orbital computations |
| `logging` | Structured logging via `FALAK_LOG` env var (feature-gated) |

## Feature Flags

| Feature | Default | Description |
|---------|---------|-------------|
| `default` | — | No optional features enabled |
| `soorat-compat` | no | Compatibility layer for soorat rendering integration |
| `logging` | no | Structured tracing via `FALAK_LOG` env var |
| `physics` | no | Physics engine integration via impetus |

```toml
[dependencies]
falak = { version = "0.2", features = ["logging"] }
```

## Quick Start

### LEO Orbit — Period Calculation

```rust
use falak::orbit::OrbitalElements;
use falak::kepler;

// ISS-like orbit: 400 km altitude, nearly circular, 51.6 deg inclination
let iss = OrbitalElements::new(
    6_778_000.0,  // semi-major axis (m)
    0.0001,       // eccentricity (nearly circular)
    0.9006,       // inclination (51.6 deg in radians)
    1.2,          // RAAN
    0.5,          // argument of periapsis
    0.0,          // true anomaly
).unwrap();

let period = kepler::orbital_period(iss.semi_major_axis(), kepler::MU_EARTH);
let period_min = period / 60.0;
println!("ISS period: {period_min:.1} min"); // ~92.6 min
```

### Hohmann Transfer — LEO to GEO

```rust
use falak::transfer::hohmann;

let r_leo = 6_778_000.0;   // 400 km altitude (m)
let r_geo = 42_164_000.0;  // GEO radius (m)
let mu = 3.986004418e14;   // Earth gravitational parameter (m^3/s^2)

let result = hohmann::hohmann_transfer(r_leo, r_geo, mu).unwrap();
println!("Total delta-v: {:.0} m/s", result.total_delta_v()); // ~3854 m/s
```

### Lambert Problem — Interplanetary Transfer

```rust
use falak::transfer::lambert;

let r1 = [1.0e11, 0.0, 0.0];       // departure position (m)
let r2 = [0.0, 1.5e11, 0.0];       // arrival position (m)
let tof = 200.0 * 86400.0;          // time of flight (s)
let mu = 1.327e20;                   // Sun gravitational parameter

let solution = lambert::solve(r1, r2, tof, mu, false).unwrap();
println!("Departure v: {:.0} m/s", solution.v1_magnitude());
```

### N-Body Simulation

```rust
use falak::nbody::{Body, NBodySimulation};

let bodies = vec![
    Body::new("Sun",   1.989e30, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]),
    Body::new("Earth", 5.972e24, [1.496e11, 0.0, 0.0], [0.0, 29_780.0, 0.0]),
    Body::new("Mars",  6.417e23, [2.279e11, 0.0, 0.0], [0.0, 24_077.0, 0.0]),
];

let mut sim = NBodySimulation::new(bodies);
sim.step_leapfrog(3600.0); // advance 1 hour
```

## Validated Results

All implementations validated against known astrodynamics results:

| Test | Expected | Verified |
|------|----------|----------|
| ISS orbital period (400 km) | ~92.6 min | 92.6 min |
| Hohmann LEO to GEO total delta-v | ~3854 m/s | 3854 m/s |
| GEO orbital period | ~23.93 hr | 23.93 hr |
| J2 RAAN drift (sun-synchronous) | ~0.9856 deg/day | confirmed |
| Lagrange point L1 (Earth-Moon) | ~0.8369 (normalized) | confirmed |
| Jacobi constant conservation | invariant under CR3BP | confirmed |
| Kepler's equation (Newton-Raphson) | converges for all e < 1 | confirmed |
| Vis-viva equation | v = sqrt(mu(2/r - 1/a)) | confirmed |
| Hohmann transfer orbit (energy) | minimum two-impulse transfer | confirmed |
| N-body energy conservation | leapfrog symplectic | confirmed |

## Relationship to AGNOS Science Stack

```
hisab (math) ──┐
               ├── falak (orbits) ──┬── kiran (game engine)
impetus (phys)─┘                   ├── joshua (simulation)
                                   └── soorat (rendering)
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
- **impetus** — physics engine (integration bridge)

## Documentation

- [Architecture Overview](docs/architecture/overview.md) — module map, data flow, dependency stack
- [Development Roadmap](docs/development/roadmap.md) — milestones and planned work

## License

GPL-3.0-only. See [LICENSE](LICENSE).
