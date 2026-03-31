# Architecture Overview

> Falak — Orbital mechanics and celestial dynamics for AGNOS

## System Diagram

```
                    ┌─────────────────────────────────────────┐
                    │           hisab (math foundation)        │
                    │  Linear algebra, geometry, calculus,      │
                    │  spatial structures                       │
                    └──────────────────┬──────────────────────┘
                                       │
                    ┌──────────────────┴──────────────────────┐
                    │              falak                        │
                    │    orbital mechanics & celestial dynamics │
                    ├─────────────────────────────────────────┤
                    │                                         │
                    │  ┌───────────┐   ┌──────────────┐      │
                    │  │  orbit    │   │  ephemeris    │      │
                    │  │  elements │   │  Julian date  │      │
                    │  │  kepler   │   │  planets/moon │      │
                    │  │  transfer │   │  eclipse      │      │
                    │  └─────┬─────┘   └──────┬───────┘      │
                    │        │                │               │
                    │  ┌─────┴─────┐   ┌──────┴───────┐      │
                    │  │maneuver   │   │  frame       │      │
                    │  │burns, dv  │   │  ECI/ECEF    │      │
                    │  │low-thrust │   │  geodetic    │      │
                    │  │Edelbaum   │   │  precession  │      │
                    │  └─────┬─────┘   └──────┬───────┘      │
                    │        │                │               │
                    │  ┌─────┴────────────────┴───────┐      │
                    │  │     perturbation + propagate   │      │
                    │  │  J2/J3, drag, SRP, third-body │      │
                    │  │  Kepler/Cowell/Encke/Brouwer  │      │
                    │  └──────────────┬────────────────┘      │
                    │                 │                        │
                    │  ┌──────────────┴────────────────┐      │
                    │  │  nbody          cr3bp          │      │
                    │  │  N-body sim     Lagrange pts   │      │
                    │  │  Barnes-Hut     Jacobi const   │      │
                    │  └──────────────────────────────┘       │
                    │                                         │
                    └─────────────────┬───────────────────────┘
                                      │
              ┌───────────────────────┼───────────────────────┐
              │                       │                       │
       ┌──────┴──────┐      ┌────────┴────────┐    ┌────────┴───────┐
       │   kiran     │      │    joshua       │    │    soorat      │
       │ (game       │      │ (celestial     │    │ (visual world) │
       │  engine)    │      │  simulation)   │    │                │
       └─────────────┘      └─────────────────┘    └────────────────┘
```

## Module Structure

### `error` (always available)

Error types and input validation for all fallible operations.

| File | Key Types | Purpose |
|------|-----------|---------|
| `error.rs` | `FalakError`, `Result<T>` | Six error variants with `Cow<'static, str>` messages |
| | `require_finite()`, `ensure_finite()` | Input/output NaN/Inf guards at API boundaries |

### `orbit` (always available)

Classical Keplerian orbital elements.

| File | Key Types | Purpose |
|------|-----------|---------|
| `orbit.rs` | `OrbitalElements` | Six classical elements (a, e, i, RAAN, omega, nu) with validation |

### `kepler` (always available)

Kepler's equation solver, anomaly conversions, state vector transforms.

| File | Key Types | Purpose |
|------|-----------|---------|
| `kepler.rs` | `StateVector` | Position + velocity state `[x, y, z]` pairs |
| | `solve_kepler_elliptic()` | Newton-Raphson Kepler solver for e < 1 |
| | `solve_kepler_hyperbolic()` | Kepler solver for e > 1 |
| | `elements_to_state()`, `state_to_elements()` | Orbital elements <-> state vector conversion |
| | Anomaly conversions | Mean, eccentric, true, hyperbolic anomaly transforms |

### `transfer` (always available)

Orbital transfer maneuvers -- Hohmann, bi-elliptic, plane change, Lambert.

| File | Key Types | Purpose |
|------|-----------|---------|
| `transfer.rs` | `HohmannTransfer` | Two-impulse minimum-energy coplanar transfer |
| | `BiEllipticTransfer` | Three-impulse transfer via intermediate apoapsis |
| | `PlaneChange` | Simple inclination change maneuver |
| | `LambertSolution` | Two-point boundary value problem solver |
| | `CombinedManeuver` | Plane change + orbit raise in single burn |

### `perturbation` (always available)

Perturbation accelerations for orbit propagation.

| File | Key Types | Purpose |
|------|-----------|---------|
| `perturbation.rs` | `j2_acceleration()`, `j3_acceleration()` | Oblateness perturbation (zonal harmonics) |
| | `drag_acceleration()` | Atmospheric drag with `AtmosphereParams` |
| | `srp_acceleration()` | Solar radiation pressure at arbitrary AU |
| | `third_body_acceleration()` | Third-body gravitational perturbation |
| | `j2_secular_rates()` | Secular RAAN/omega drift rates |

### `nbody` (always available)

N-body gravitational simulation with multiple integrators.

| File | Key Types | Purpose |
|------|-----------|---------|
| `nbody.rs` | `Body` | Gravitating body with position, velocity, mass, optional mu |
| | `System` | Collection of bodies with softening parameter |
| | `Integrator` | Enum: `Leapfrog`, `Rk4` |
| | `step_leapfrog()`, `step_rk4()` | Symplectic and explicit integrators |
| | `step_adaptive()` | Adaptive-step RK4 with error control |
| | `compute_accelerations_barnes_hut()` | O(N log N) tree-based gravity |

### `ephemeris` (always available)

Julian dates, sidereal time, planetary/lunar positions, eclipse detection.

| File | Key Types | Purpose |
|------|-----------|---------|
| `ephemeris.rs` | `PlanetaryPosition`, `Planet` | Simplified planetary ephemeris (ecliptic lon/lat/distance) |
| | `LunarPosition` | Lunar position model |
| | `RiseTransitSet` | Observer-centric rise/transit/set times |
| | `EclipseInfo`, `EclipseState` | Cylindrical and conical eclipse shadow models |
| | `calendar_to_jd()`, `gmst()` | Julian date conversions, Greenwich sidereal time |

### `frame` (always available)

Reference frame transformations.

| File | Key Types | Purpose |
|------|-----------|---------|
| `frame.rs` | `Geodetic` | Latitude, longitude, altitude (WGS84) |
| | `perifocal_to_eci()`, `eci_to_perifocal()` | Perifocal <-> ECI rotation |
| | `eci_to_ecef()`, `ecef_to_eci()` | Inertial <-> rotating Earth |
| | `ecef_to_geodetic()`, `geodetic_to_ecef()` | Cartesian <-> WGS84 geodetic |
| | `inertial_to_rotating()`, `rotating_to_inertial()` | Generic rotating frame transform |
| | `precession_angles()`, `nutation()` | IAU precession/nutation models |

### `maneuver` (always available)

Spacecraft maneuver planning -- burns, rocket equation, low-thrust.

| File | Key Types | Purpose |
|------|-----------|---------|
| `maneuver.rs` | `ImpulsiveBurn` | Delta-v vector in local orbital frame [prograde, normal, radial] |
| | `ManeuverPlan`, `ScheduledBurn` | Ordered sequence of timed burns |
| | `propellant_mass_fraction()`, `max_delta_v()` | Tsiolkovsky rocket equation |
| | `escape_delta_v()`, `capture_delta_v()` | Hyperbolic excess velocity burns |
| | `LowThrustTransfer` | Continuous low-thrust spiral transfer |
| | `edelbaum_delta_v()` | Edelbaum combined plane-change/orbit-raise |

### `propagate` (always available)

Orbit propagation -- analytic two-body, numerical Cowell/Encke, Brouwer mean elements.

| File | Key Types | Purpose |
|------|-----------|---------|
| `propagate.rs` | `kepler()` | Analytic two-body propagation via mean anomaly |
| | `cowell()` | Numerical propagation with arbitrary perturbation forces |
| | `encke()` | Encke's method -- propagate deviation from reference orbit |
| | `two_body()` | Multi-step two-body numerical propagation |
| | `MeanElements` | Brouwer mean orbital elements |
| | `osculating_to_mean()`, `propagate_mean_elements()` | Brouwer theory conversions |
| | `PerturbationFn` | Trait alias for perturbation acceleration callbacks |

### `cr3bp` (always available)

Circular Restricted Three-Body Problem.

| File | Key Types | Purpose |
|------|-----------|---------|
| `cr3bp.rs` | `LagrangePoints` | All five libration points (L1-L5) via Newton's method |
| | `MassRatio`, `SynodicState` | Type aliases for CR3BP normalized coordinates |
| | `jacobi_constant()` | Energy-like integral of motion |
| | `equations_of_motion()` | State derivative in the synodic frame |
| | `pseudo_potential()`, `zero_velocity_value()` | Effective potential and zero-velocity curves |

### `bridge` (always available)

Cross-crate primitive-value bridges for other AGNOS science crates.

| File | Key Types | Purpose |
|------|-----------|---------|
| `bridge.rs` | `stellar_mass_to_mu()` | Tara bridge: stellar mass -> gravitational parameter |
| | `luminosity_to_habitable_zone_au()` | Tara bridge: luminosity -> habitable zone distance |
| | `orbital_to_gravity_force()` | Impetus bridge: position -> gravitational force |
| | `escape_energy_deficit()` | Impetus bridge: bound/escaping orbit check |
| | `solar_distance_to_insolation()` | Badal bridge: distance -> insolation |
| | `eccentricity_to_climate_variation()` | Badal bridge: orbital -> climate forcing |

### `integration::soorat` (feature: `soorat-compat`)

Visualization data structures for soorat rendering.

| File | Key Types | Purpose |
|------|-----------|---------|
| `integration/soorat.rs` | `OrbitPath` | Trajectory points with velocity color mapping |
| | `GroundTrack` | Sub-satellite ground track (lat/lon) |
| | `TransferTrajectory` | Transfer orbit visualization with departure/arrival |
| | `PlanetaryPositions`, `CelestialBody` | Solar system snapshot for scene rendering |

### `logging` (feature: `logging`)

Structured tracing subscriber setup.

| File | Key Types | Purpose |
|------|-----------|---------|
| `logging.rs` | `init()`, `init_with_level()` | tracing-subscriber with `FALAK_LOG` env var filter |

## Data Flow

```
OrbitalElements (orbit)
    │
    ├──► kepler ──► StateVector, anomaly conversions, period, velocity
    │                    │
    │                    ├──► transfer ──► HohmannTransfer, LambertSolution
    │                    │                    │
    │                    │                    ▼
    │                    │               maneuver ──► ImpulsiveBurn, ManeuverPlan
    │                    │                    │
    │                    │                    ▼
    │                    │               propagate ──► kepler/cowell/encke propagation
    │                    │                    │
    │                    │                    ▼
    │                    │            perturbation ──► J2/J3, drag, SRP accelerations
    │                    │
    │                    └──► nbody ──► System evolution (leapfrog/RK4/Barnes-Hut)
    │
    └──► cr3bp ──► LagrangePoints, Jacobi constant

ephemeris ──► Julian dates, GMST, planetary/lunar positions, eclipse
    │
    ▼
frame ──► ECI/ECEF/perifocal/geodetic coordinate transforms

bridge ──► tara/impetus/badal primitive-value conversions (always available)

integration::soorat ──► OrbitPath, GroundTrack, TransferTrajectory (soorat-compat)
```

**Bridge pattern**: All bridge functions are dependency-free -- they take f64 primitives, not external crate types. The game engine orchestrates, calling science crates and bridge functions separately.

## Error Handling

Every public fallible function returns `Result<T, FalakError>`. Three-layer validation:

1. **Input**: `require_finite()` -- rejects NaN/Inf at API boundaries before any computation
2. **Domain**: Physics boundary checks with `tracing::warn!` (negative eccentricity, non-positive semi-major axis, invalid mass ratio, convergence failure)
3. **Output**: `ensure_finite()` -- catches overflow and 0/0 in computed results

`FalakError` variants:
- `InvalidParameter` -- out-of-range inputs
- `MathError` -- division by zero, domain errors
- `ConvergenceError` -- iterative solver failure (with iteration count)
- `EphemerisError` -- ephemeris lookup failure
- `NonFinite` -- NaN/Inf detected (with context and value)
- `Io` -- I/O errors (wraps `std::io::Error`)

All messages use `Cow<'static, str>` -- zero allocation on static messages, heap only when formatting dynamic values.

## Observability

Structured tracing via the `tracing` crate at three levels:

| Level | Usage | Example |
|-------|-------|---------|
| `error!` | Computational failures (convergence, invalid state) | `error!(%e, "Kepler solver diverged")` |
| `warn!` | Domain boundary violations (negative radius, bad eccentricity) | `warn!(e, "eccentricity out of range")` |
| `#[instrument]` | Span instrumentation on public functions | Entry/exit tracing with arguments |

Span levels: `trace` for individual computations (anomaly conversions, frame transforms), `debug` for higher-level aggregators (propagation, transfer planning).

Enable with the `logging` feature and `FALAK_LOG` env var (e.g., `FALAK_LOG=debug`).

## Constants

Distributed across modules, SI units throughout. Key sources:

| Module | Constants | Source |
|--------|-----------|--------|
| `kepler` | `G` (6.674e-11 m^3 kg^-1 s^-2) | CODATA 2018 |
| `perturbation` | `J2_EARTH`, `J3_EARTH`, `R_EARTH`, `MU_SUN`, `AU_METRES` | WGS84, IAU |
| `ephemeris` | `J2000_JD`, `UNIX_EPOCH_JD`, `SECONDS_PER_DAY`, `EARTH_ROTATION_RATE` | IAU/IERS 2010 |
| `frame` | `WGS84_A`, `WGS84_F`, `WGS84_E2` | WGS84 ellipsoid |

All angular values in radians, distances in metres, times in seconds.

## Feature Flags

| Feature | Default | Depends On | Description |
|---------|---------|------------|-------------|
| `soorat-compat` | no | -- | Soorat visualization types (OrbitPath, GroundTrack, etc.) |
| `logging` | no | `dep:tracing-subscriber` | Structured tracing subscriber with `FALAK_LOG` env filter |
| `physics` | no | `dep:impetus` | Impetus physics engine integration |

## Consumers

- **kiran** -- game engine: orbital mechanics for space simulations, mission planning UI
- **joshua** -- simulation: celestial body dynamics, multi-body evolution, mission analysis
- **soorat** -- visual engine: consumes `integration::soorat` types for orbit/ground-track rendering
- **impetus** -- physics engine: force/torque integration via bridge functions (optional `physics` feature)
