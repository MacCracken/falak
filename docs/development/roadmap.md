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
- [ ] General perturbation theory (osculating → mean elements)

### N-Body Simulation — Extended
- [ ] Barnes-Hut tree approximation for large N
- [ ] Restricted three-body problem (Lagrange points)

### Ephemeris — Extended
- [ ] Rise/set/transit times
- [ ] Eclipse prediction
- [ ] Precession and nutation corrections

### Transfer Maneuvers — Extended
- [ ] Lambert problem solver (given two positions and time)
- [ ] Combined plane change + altitude maneuvers

### Maneuver Planning — Extended
- [ ] Continuous low-thrust trajectory modelling

### Cross-Crate Bridges — Extended
- [ ] **tara bridge**: additional stellar property conversions
- [ ] **impetus bridge**: orbital velocity → gravitational force vector
- [ ] **badal bridge**: additional climate coupling parameters

### Soorat Integration — Extended
- [ ] Transfer trajectory generation from Hohmann/Lambert results

## Future (demand-gated)

- Orbit determination from observations
- Conjunction analysis / collision avoidance
- Interplanetary trajectory design (patched conics)
- Gravity assist (flyby) planning
- Satellite constellation design (Walker delta/star)
- TLE / SGP4 propagation
- Relativistic corrections
- Higher-order zonal harmonics (J4+)
