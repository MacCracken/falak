# Falak Roadmap

> **Falak** is the orbital mechanics engine. Math foundations come from [hisab](https://github.com/MacCracken/hisab). Physics integration is via [impetus](https://github.com/MacCracken/impetus).

## Scope

Falak owns the **mechanics of orbits**: how bodies move under gravity, how spacecraft transfer between orbits, how perturbations alter trajectories, and how to compute ephemerides. It provides the math; consumers decide what to do with it (simulate missions, render planets, plan maneuvers).

Falak does NOT own:
- **Physics engine** --> impetus (rigid bodies, forces, collisions)
- **Math primitives** --> hisab (vectors, geometry, calculus)
- **Rendering** --> soorat/kiran (visual display of orbits)
- **Game logic** --> joshua/kiran (mission scenarios, NPCs)

## Backlog

### Orbit Propagation — Extended
- [ ] Encke's method (deviation from reference orbit)
- [ ] General perturbation theory (osculating → mean elements)

### N-Body Simulation — Extended
- [ ] Pre-allocated integration buffers (RK4 allocates ~10 Vecs per step — dominates runtime for large N)
- [ ] Canonical μ-based gravity option (current uses G×M which diverges from standard μ values)
- [ ] Dormand-Prince (RK45) adaptive step integrator
- [ ] Barnes-Hut tree approximation for large N
- [ ] Restricted three-body problem (Lagrange points)

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
