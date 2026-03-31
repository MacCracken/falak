# Falak Roadmap

> **Falak** — the world's leading orbital mechanics engine in Rust. Math foundations from [hisab](https://github.com/MacCracken/hisab). Physics integration via [impetus](https://github.com/MacCracken/impetus).

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

---

## v1.1.0 — Precision & Interoperability

> **Theme**: Production-ready precision. Complete time systems, high-fidelity models, industry-standard formats. After 1.1.0, falak can match real-world tracking data.

### Time Systems & Precision
- [ ] Full time scale conversions: UTC ↔ TAI ↔ TT ↔ TDB (leap second table)
- [ ] UT1−UTC corrections (delta-UT1 from IERS data)
- [ ] IERS Earth Orientation Parameters (polar motion xp, yp)
- [ ] GAST (Greenwich Apparent Sidereal Time) from GMST + equation of equinoxes

### High-Fidelity Perturbations
- [ ] Higher-order zonal harmonics (J4–J6, configurable degree/order)
- [ ] NRLMSISE-00 atmosphere model (solar flux F10.7, geomagnetic Kp)
- [ ] Solid Earth tides (IERS 2010 conventions)
- [ ] Albedo + Earth infrared radiation pressure
- [ ] Relativistic corrections (Schwarzschild + Lense-Thirring)

### TLE & SGP4
- [ ] TLE parser (Two-Line Element set, standard + 3-line formats)
- [ ] SGP4/SDP4 propagator (Vallado's improved implementation)
- [ ] TLE ↔ osculating element conversion
- [ ] Catalog matching utilities

### Interoperability
- [ ] CCSDS OEM (Orbit Ephemeris Message) read/write
- [ ] CCSDS OPM (Orbit Parameter Message) read/write
- [ ] SP3 format read (GNSS precise ephemeris)

### Orbit Determination — Foundation
- [ ] Measurement models: range, range-rate (Doppler), angles (az/el)
- [ ] Batch least-squares orbit determination
- [ ] State transition matrix (STM) propagation
- [ ] Covariance propagation (linear)

### Propagation
- [ ] Event detection framework (apoapsis, periapsis, eclipse entry/exit, altitude crossing)
- [ ] Bulirsch-Stoer integrator (high-order variable-step)

### Testing & Validation
- [ ] Vallado appendix test cases (regression suite, ≥50 cases)
- [ ] IAU SOFA cross-validation (precession, nutation, sidereal time)
- [ ] Adversarial fuzzing with `proptest` (NaN, ±Inf, edge-case eccentricities)
- [ ] Benchmark suite expansion (all modules, 30+ benchmarks)

---

## v1.2.0 — Mission Design & Operations

> **Theme**: World-class mission design. Interplanetary trajectories, constellation planning, real-time operations. After 1.2.0, falak is the definitive Rust astrodynamics library.

### Interplanetary Trajectory Design
- [ ] Patched conics (sphere of influence transitions)
- [ ] Gravity assist (flyby) planning with B-plane targeting
- [ ] Powered flyby modelling
- [ ] Pork-chop plots (launch window analysis)
- [ ] Interplanetary Lambert arcs (heliocentric)

### Orbit Determination — Advanced
- [ ] Extended Kalman Filter (EKF)
- [ ] Unscented Kalman Filter (UKF)
- [ ] Consider parameters (drag coefficient, SRP coefficient estimation)
- [ ] Process noise tuning utilities

### Conjunction & Collision Avoidance
- [ ] Conjunction screening (miss distance, time of closest approach)
- [ ] Probability of collision (Alfano / Chan methods)
- [ ] Conjunction Data Message (CDM) support
- [ ] Maneuver planning for collision avoidance

### Constellation Design
- [ ] Walker delta/star pattern generation
- [ ] Ground coverage analysis (revisit time, coverage percentage)
- [ ] Station visibility windows (access intervals with mask angles)
- [ ] Link budget estimation (free-space path loss, SNR)

### Ground Operations
- [ ] Antenna pointing (azimuth/elevation from geodetic station)
- [ ] Doppler prediction (range-rate to frequency shift)
- [ ] Pass scheduling (multi-station, priority-based)

### Advanced Dynamics
- [ ] Finite burn modelling (continuous thrust with mass depletion)
- [ ] Manifold computation (stable/unstable manifolds at Lagrange points)
- [ ] Halo orbit families (differential correction)
- [ ] Lyapunov orbit computation

### Ephemeris — Full Precision
- [ ] JPL DE440+ reader (binary SPK kernel support)
- [ ] Full VSOP87 planetary theory
- [ ] High-precision lunar ephemeris (ELP 2000-82 or DE-based)

### Testing & Validation
- [ ] GMAT cross-validation suite (propagation, OD, maneuvers)
- [ ] 500+ validated test cases
- [ ] Performance benchmarks vs. C/Fortran reference implementations
- [ ] Formal verification of critical numerical paths

---

## Completed

### v0.1.0 (2026-03-25)
- Core orbital elements, error types, module scaffolding

### v1.0.0 (2026-03-31)
- Kepler's equation solver (elliptic + hyperbolic), all anomaly conversions
- Hohmann, bi-elliptic, plane change, phasing, combined maneuver transfers
- Lambert problem solver (universal variables + Stumpff functions)
- J2/J3 oblateness, atmospheric drag, SRP, third-body perturbations
- N-body simulation (direct O(N²), Barnes-Hut O(N log N))
- Leapfrog, RK4, adaptive RK45 integrators
- Restricted three-body problem (Lagrange points, Jacobi constant, CR3BP EOM)
- Ephemeris: Julian dates, GMST, simplified planetary/lunar positions
- Eclipse prediction (cylindrical + conical shadow models)
- Rise/set/transit times (Meeus algorithm)
- Reference frame transforms (ECI, ECEF, perifocal, rotating, geodetic)
- IAU 2006 precession, IAU 1980 nutation (9 terms), equation of equinoxes
- Impulsive burns, rocket equation, escape/capture, Oberth effect
- Low-thrust spiral transfers, Edelbaum approximation
- Analytic (Kepler), Cowell, Encke propagation methods
- General perturbation theory (Brouwer mean elements, secular J2 rates)
- Soorat integration (orbit paths, ground tracks, transfer trajectories)
- Cross-crate bridges (tara, impetus, badal)
- Input validation helpers (require_finite, ensure_finite)
- 5 runnable examples, 5 ADRs, comprehensive architecture docs
- 257 tests (249 unit + 5 integration + 3 doc-tests)
