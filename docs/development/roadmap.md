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

| Phase | Summary |
|-------|---------|
| Foundation | Error types, OrbitalElements struct with validation (elliptical/parabolic/hyperbolic) |
| Keplerian Orbits | Kepler's equation solver (elliptic + hyperbolic), all anomaly conversions, orbital period/mean motion, vis-viva, orbital radius, state vector ↔ elements |
| Transfer Maneuvers | Hohmann transfer, bi-elliptic transfer, plane change, phasing orbits |
| Perturbation Models | J2/J3 zonal harmonics, secular J2 rates, atmospheric drag (exponential model), solar radiation pressure, third-body perturbations |
| Reference Frames | Perifocal ↔ ECI, ECI ↔ ECEF, ECEF ↔ geodetic (WGS-84 Bowring), inertial ↔ rotating (synodic) |
| Ephemeris | Calendar ↔ Julian Date, MJD, Unix timestamp, GMST (IAU 1982), Julian centuries, day of year |
| Maneuvers | Impulsive burn types, delta-v budgets, Tsiolkovsky rocket equation, escape/capture delta-v, Oberth effect |
| Cross-Crate Bridges | Tara (stellar mass → μ, luminosity → HZ), impetus (gravity force, escape energy), badal (insolation, climate variation) |
| Soorat Integration | OrbitPath from elements, planetary position/transfer/ground track data types |
| Logging | Feature-gated structured tracing via FALAK_LOG |

## Backlog

### N-Body Simulation
- [ ] Direct N-body gravitational computation
- [ ] Leapfrog / Stormer-Verlet integrator
- [ ] Runge-Kutta 4/5 (Dormand-Prince) integrator
- [ ] Symplectic integrators for long-term stability
- [ ] Barnes-Hut tree approximation for large N
- [ ] Restricted three-body problem (Lagrange points)

### Orbit Propagation
- [ ] Two-body propagation (Kepler problem forward in time)
- [ ] Perturbed orbit propagation (J2 + drag + SRP + third-body)
- [ ] Cowell's method (direct integration of equations of motion)
- [ ] Encke's method (deviation from reference orbit)
- [ ] General perturbation theory (osculating → mean elements)

### Ephemeris — Extended
- [ ] VSOP87 planetary positions (truncated series)
- [ ] Simple lunar ephemeris
- [ ] Rise/set/transit times
- [ ] Eclipse prediction
- [ ] Precession and nutation corrections

### Transfer Maneuvers — Extended
- [ ] Lambert problem solver (given two positions and time)
- [ ] Combined plane change + altitude maneuvers

### Maneuver Planning — Extended
- [ ] Maneuver plan sequencing (multiple burns with timing)
- [ ] Continuous low-thrust trajectory modelling

### Cross-Crate Bridges — Extended
- [ ] **tara bridge**: additional stellar property conversions
- [ ] **impetus bridge**: orbital velocity → gravitational force vector
- [ ] **badal bridge**: additional climate coupling parameters

### Soorat Integration — Extended
- [ ] Orbit path in ECI frame (currently perifocal only)
- [ ] Transfer trajectory generation from Hohmann/Lambert results
- [ ] Ground track computation from orbital elements + GMST

## Future (demand-gated)

- Orbit determination from observations
- Conjunction analysis / collision avoidance
- Interplanetary trajectory design (patched conics)
- Gravity assist (flyby) planning
- Satellite constellation design (Walker delta/star)
- TLE / SGP4 propagation
- Relativistic corrections
- Higher-order zonal harmonics (J4+)
