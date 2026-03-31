# Falak Roadmap

> **Falak** is the orbital mechanics engine. Math foundations come from [hisab](https://github.com/MacCracken/hisab). Physics integration is via [impetus](https://github.com/MacCracken/impetus).

## Scope

Falak owns the **mechanics of orbits**: how bodies move under gravity, how spacecraft transfer between orbits, how perturbations alter trajectories, and how to compute ephemerides. It provides the math; consumers decide what to do with it (simulate missions, render planets, plan maneuvers).

Falak does NOT own:
- **Physics engine** → impetus (rigid bodies, forces, collisions)
- **Math primitives** → hisab (vectors, geometry, calculus)
- **Rendering** → soorat/kiran (visual display of orbits)
- **Game logic** → joshua/kiran (mission scenarios, NPCs)

## Dependencies

| Crate | Version | Role | Bridge |
|-------|---------|------|--------|
| [hisab](https://github.com/MacCracken/hisab) | 1.x | Math foundation (vectors, trig, calculus) | Direct dependency |
| [impetus](https://github.com/MacCracken/impetus) | 1.x | Physics integration (forces, energy) | Feature-gated (`physics`) |
| tara | — | Stellar astrophysics | Primitive-value bridge (`bridge.rs`) |
| badal | — | Weather/climate modelling | Primitive-value bridge (`bridge.rs`) |

## Completed

### v0.1.0 (2026-03-25)
- Core orbital elements, error types, module scaffolding

### v1.0.0 (2026-03-31)
- Kepler's equation solver (elliptic + hyperbolic), all anomaly conversions
- Hohmann, bi-elliptic, plane change, phasing transfers
- Lambert problem solver (universal variables + Stumpff functions)
- Combined plane change + altitude maneuver
- J2/J3 oblateness, atmospheric drag, SRP, third-body perturbations
- N-body simulation (direct O(N²), Barnes-Hut O(N log N))
- Leapfrog, RK4, adaptive RK45 integrators
- Restricted three-body problem (Lagrange points, Jacobi constant)
- Ephemeris: Julian dates, GMST, simplified planetary/lunar positions
- Eclipse prediction (cylindrical + conical shadow models)
- Rise/set/transit times (Meeus algorithm)
- Reference frame transforms (ECI, ECEF, perifocal, rotating, geodetic)
- IAU 2006 precession, IAU 1980 nutation, equation of equinoxes
- Impulsive burns, rocket equation, escape/capture, Oberth effect
- Low-thrust spiral transfers, Edelbaum approximation
- Analytic (Kepler), Cowell, Encke propagation methods
- General perturbation theory (Brouwer mean elements, secular J2 rates)
- Soorat integration (orbit paths, ground tracks, transfer trajectories)
- Cross-crate bridges (tara, impetus, badal)

## Future (demand-gated)

These features will be implemented when downstream consumers request them:

- Orbit determination from observations
- Conjunction analysis / collision avoidance
- Interplanetary trajectory design (patched conics)
- Gravity assist (flyby) planning
- Satellite constellation design (Walker delta/star)
- TLE / SGP4 propagation
- Relativistic corrections
- Higher-order zonal harmonics (J4+)
