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
- [x] General perturbation theory (osculating → mean elements)

### N-Body Simulation — Extended
- [x] Barnes-Hut tree approximation for large N
- [x] Restricted three-body problem (Lagrange points)

### Ephemeris — Extended
- [x] Rise/set/transit times
- [x] Eclipse prediction
- [x] Precession and nutation corrections

### Transfer Maneuvers — Extended
- [x] Lambert problem solver (given two positions and time)
- [x] Combined plane change + altitude maneuvers

### Maneuver Planning — Extended
- [x] Continuous low-thrust trajectory modelling

### Cross-Crate Bridges — Extended
- [x] **tara bridge**: additional stellar property conversions
- [x] **impetus bridge**: orbital velocity → gravitational force vector
- [x] **badal bridge**: additional climate coupling parameters

### Soorat Integration — Extended
- [x] Transfer trajectory generation from Hohmann/Lambert results

## Future (demand-gated)

- Orbit determination from observations
- Conjunction analysis / collision avoidance
- Interplanetary trajectory design (patched conics)
- Gravity assist (flyby) planning
- Satellite constellation design (Walker delta/star)
- TLE / SGP4 propagation
- Relativistic corrections
- Higher-order zonal harmonics (J4+)
