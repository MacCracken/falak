# Falak Roadmap

> **Falak** is the orbital mechanics engine. Math foundations come from [hisab](https://github.com/MacCracken/hisab). Physics integration is via [impetus](https://github.com/MacCracken/impetus).

## Scope

Falak owns the **mechanics of orbits**: how bodies move under gravity, how spacecraft transfer between orbits, how perturbations alter trajectories, and how to compute ephemerides. It provides the math; consumers decide what to do with it (simulate missions, render planets, plan maneuvers).

Falak does NOT own:
- **Physics engine** --> impetus (rigid bodies, forces, collisions)
- **Math primitives** --> hisab (vectors, geometry, calculus)
- **Rendering** --> soorat/kiran (visual display of orbits)
- **Game logic** --> joshua/kiran (mission scenarios, NPCs)

## Completed

| Phase | Release | Summary |
|-------|---------|---------|
| V0.1 | Foundation | Error types, OrbitalElements struct with validation, module stubs for all domains |

## Backlog

### Keplerian Orbits
- [ ] Kepler's equation solver (Newton-Raphson with Danby starter)
- [ ] Mean, eccentric, and true anomaly conversions
- [ ] Orbital period and mean motion
- [ ] Vis-viva equation (orbital velocity at any point)
- [ ] Hyperbolic and parabolic orbit support (e >= 1)
- [ ] State vector <--> orbital elements conversion

### Transfer Maneuvers
- [ ] Hohmann transfer (delta-v, time of flight)
- [ ] Bi-elliptic transfer (when more efficient than Hohmann)
- [ ] Lambert problem solver (given two positions and time)
- [ ] Plane change maneuvers
- [ ] Phasing orbits

### Perturbation Models
- [ ] J2 oblateness (secular and periodic)
- [ ] Higher-order zonal harmonics (J3, J4)
- [ ] Atmospheric drag (exponential atmosphere model)
- [ ] Solar radiation pressure
- [ ] Third-body perturbations (Sun, Moon)
- [ ] General perturbation theory (osculating elements)

### N-Body Simulation
- [ ] Direct N-body gravitational computation
- [ ] Leapfrog / Stormer-Verlet integrator
- [ ] Runge-Kutta 4/5 (Dormand-Prince) integrator
- [ ] Symplectic integrators for long-term stability
- [ ] Barnes-Hut tree approximation for large N
- [ ] Restricted three-body problem (Lagrange points)

### Ephemeris Computation
- [ ] Julian date / Modified Julian date conversions
- [ ] Sidereal time computation
- [ ] VSOP87 planetary positions (truncated series)
- [ ] Simple lunar ephemeris
- [ ] Rise/set/transit times
- [ ] Eclipse prediction

### Reference Frames
- [ ] ECI (Earth-Centered Inertial) frame
- [ ] ECEF (Earth-Centered Earth-Fixed) frame
- [ ] Perifocal frame
- [ ] Rotating (synodic) frame for CR3BP
- [ ] Frame transformation utilities
- [ ] Precession and nutation corrections

## Cross-Crate Bridges

- [ ] `bridge.rs` module — primitive-value conversions for cross-crate orbital mechanics
- [ ] **tara bridge**: stellar mass (kg) → gravitational parameter μ; luminosity (W) → habitable zone distance (AU)
- [ ] **impetus bridge**: orbital velocity [f32; 3] → gravitational force vector; escape velocity → kinetic energy threshold
- [ ] **badal bridge**: solar distance (AU), axial tilt (°) → seasonal insolation; orbital eccentricity → climate variation amplitude

## Soorat Integration

- [ ] `integration/soorat.rs` module — feature-gated `soorat-compat`
- [ ] **Orbit path**: elliptical/hyperbolic trajectory points for line rendering
- [ ] **Planetary positions**: body positions at epoch for instanced sphere rendering
- [ ] **Transfer trajectory**: Hohmann/Lambert arc points with delta-v markers for colored line rendering
- [ ] **Ground track**: sub-satellite point trace on planetary surface for map overlay rendering

## Future (demand-gated)

- Ground track computation
- Orbit determination from observations
- Conjunction analysis / collision avoidance
- Interplanetary trajectory design (patched conics)
- Low-thrust trajectory optimization
- Gravity assist (flyby) planning
- Satellite constellation design (Walker delta/star)
- TLE / SGP4 propagation
- Relativistic corrections
