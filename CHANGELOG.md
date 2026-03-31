# Changelog

## [1.0.0] - 2026-03-31

### Added
- **kepler** — Kepler's equation solver (Newton-Raphson with Danby starter for elliptic, log starter for hyperbolic), all anomaly conversions (M↔E↔ν elliptic, M↔H↔ν hyperbolic), orbital period, mean motion, vis-viva velocity, orbital radius, `StateVector` type, bidirectional state vector ↔ orbital elements conversion
- **orbit** — Hyperbolic (e > 1) and parabolic (e = 1) orbit support in `OrbitalElements`, orbit type queries (`is_elliptical`, `is_circular`, `is_parabolic`, `is_hyperbolic`), correct periapsis/apoapsis for all orbit types
- **transfer** — Hohmann transfer, bi-elliptic transfer, plane change maneuver, phasing orbit; `HohmannTransfer`, `BiEllipticTransfer`, `PlaneChange` result types
- **frame** — Perifocal ↔ ECI transforms, ECI ↔ ECEF rotation, ECEF ↔ geodetic (WGS-84 Bowring's method), inertial ↔ rotating (synodic) frame; `Geodetic` type; WGS-84 constants
- **ephemeris** — Calendar ↔ Julian Date (Gregorian algorithm), Julian Date ↔ MJD, Julian Date ↔ Unix timestamp, GMST (IAU 1982), Julian centuries since J2000, day of year with leap year and month-day validation, simplified planetary positions (VSOP87-like truncated series for Mercury–Saturn), simplified lunar position (low-order series), ecliptic/lunar Cartesian converters; `Planet` enum, `PlanetaryPosition`, `LunarPosition` types; time constants (`J2000_JD`, `UNIX_EPOCH_JD`, etc.)
- **perturbation** — J2 oblateness acceleration, J2 secular drift rates (RAAN and argument of periapsis), J3 acceleration, atmospheric drag (exponential model with `AtmosphereParams`), solar radiation pressure, third-body gravitational perturbation; physical constants (`J2_EARTH`, `J3_EARTH`, `R_EARTH`, `MU_SUN`, `SOLAR_PRESSURE_1AU`, `AU_METRES`)
- **maneuver** — `ImpulsiveBurn` type (prograde/retrograde/normal/radial constructors), `ManeuverPlan` with `ScheduledBurn` sequencing (auto-sorted by time), delta-v budget summing, Tsiolkovsky rocket equation (propellant mass fraction ↔ max delta-v), escape/capture delta-v from circular orbit, Oberth effect factor
- **nbody** — `Body` and `System` types, `Body::with_mu` for canonical μ-based gravity (eliminates G×M drift), `compute_accelerations_into` for zero-allocation hot path, direct O(N²) gravitational acceleration, leapfrog (kick-drift-kick) symplectic integrator, RK4 explicit integrator, adaptive RK45 step (Richardson extrapolation with automatic step sizing), `evolve` convenience function, energy conservation (kinetic + potential), centre of mass, gravitational softening, `Integrator` enum
- **propagate** — Analytic two-body (Kepler) propagation with forward/backward time, Cowell's method (RK4) for perturbed propagation with user-supplied perturbation function, Encke's method (deviation from reference orbit with Battin's F(q)), `two_body` convenience wrapper
- **bridge** — Tara bridges (stellar mass → μ, luminosity → habitable zone), impetus bridges (gravity force, escape energy deficit), badal bridges (insolation, climate variation)
- **integration/soorat** — `OrbitPath::from_elements` (now returns `Result`), `OrbitPath::from_elements_eci` (full perifocal→ECI rotation), `GroundTrack::from_elements` (orbital propagation + GMST + geodetic conversion for map overlay), `PlanetaryPositions`, `CelestialBody`, `TransferTrajectory` data types
- **perturbation** — `AtmosphereParams::new` constructor for external crate usage (`#[non_exhaustive]` struct)
- **tests** — Cross-module integration tests: combined perturbations (J2 + drag + third-body), Encke vs Cowell agreement with combined perturbations, high-eccentricity Encke (e=0.74, e=0.95)
- **examples** — `orbit_propagation` (LEO with J2, secular rates) and `hohmann_transfer` (LEO→GEO, maneuver planning, rocket equation)
- **transfer** — Lambert problem solver: universal-variable method with Stumpff functions and Newton iteration, handles elliptic/parabolic/hyperbolic transfers, prograde and retrograde arcs; `LambertSolution` type; benchmark: 97 ns per solve
- **cr3bp** — Circular Restricted Three-Body Problem module: all 5 Lagrange points (L1–L5, Newton's method for collinear + analytic equilateral), Jacobi constant, pseudo-potential, zero-velocity curves, equations of motion in the synodic frame; `LagrangePoints`, `SynodicState` types
- **ephemeris** — Eclipse prediction: cylindrical shadow model (`eclipse_cylindrical`), conical shadow model with penumbra (`eclipse_conical`), `EclipseState` enum (Sunlit/Penumbra/Umbra), `EclipseInfo` type with shadow fraction
- **integration/soorat** — `TransferTrajectory::from_lambert` generates renderable transfer trajectory from Lambert solution (numerical propagation + delta-v markers)
- **transfer** — Combined plane change + altitude maneuver (`combined_maneuver`): single-impulse Δv with savings over separate Hohmann + plane change; `CombinedManeuver` type
- **ephemeris** — Rise/set/transit times (`rise_transit_set`): Meeus algorithm for horizon crossings and meridian transit; `RiseTransitSet` type; `STANDARD_REFRACTION` constant
- **frame** — IAU 2006 precession (`precession_angles`, `precess_j2000_to_date`, `precess_date_to_j2000`), IAU 1980 nutation series (9 dominant terms, `nutation`), mean/true obliquity (`mean_obliquity`, `true_obliquity`), equation of the equinoxes
- **propagate** — General perturbation theory: osculating → mean elements (`osculating_to_mean`, first-order Brouwer J2), mean element propagation with secular rates (`propagate_mean_elements`); `MeanElements` type
- **nbody** — Barnes-Hut octree O(N log N) gravity approximation (`compute_accelerations_barnes_hut`), configurable opening angle θ
- **maneuver** — Low-thrust spiral transfer (`low_thrust_spiral`, `LowThrustTransfer`), Edelbaum low-thrust Δv with combined plane change (`edelbaum_delta_v`)
- **bridge** — Tara: `stellar_luminosity` (Stefan-Boltzmann), `temperature_to_mass_kg`; Impetus: `orbital_to_gravity_force_vector`, `specific_orbital_energy`; Badal: `orbital_distance_extremes_au`, `obliquity_forcing`, `precession_period_years`

### Changed
- **orbit** — `OrbitalElements::new` now accepts all orbit types: elliptical (e < 1, a > 0), parabolic (e = 1, p > 0), hyperbolic (e > 1, a < 0); periapsis returns p/2 for parabolic; apoapsis returns ∞ for open orbits
- **bridge** — Gravitational constant `G` now sourced from `kepler::G` (removed duplication)
- **integration/soorat** — `OrbitPath::from_elements` returns `Result` instead of silently masking errors; all structs marked `#[non_exhaustive]`; orbit path generation uses `kepler` module functions
- **nbody** — RK4 integrator uses pre-allocated flat arrays and `compute_accelerations_into` instead of per-step Vec allocations (6 allocations eliminated per step)
- **nbody** — Adaptive RK45 uses manual state backup/restore instead of `system.clone()` (eliminates full system clone per attempt)
- **kepler** — `StateVector` doc comments now specify ECI coordinate frame convention
- **orbit** — `OrbitalElements::semi_major_axis` doc clarifies parabolic convention (stores semi-latus rectum p)
- **propagate** — `PerturbationFn` type alias now fully documents parameter semantics

### Fixed
- **ephemeris** — GMST formula was double-counting fractional day (T included full JD, then fractional day added again); now correctly separates 0h UT1 and fractional day
- **frame** — `eci_to_perifocal` was hardcoding W-component to 0.0, dropping out-of-plane ECI z-component; now computes full 3×3 transpose
- **orbit** — Parabolic periapsis was returning 0.0 (|a|×|1-e| = 0 when e=1); now returns p/2
- **perturbation** — Third-body test tolerance tightened from 4 to 2 orders of magnitude (1e-7..1e-5 vs 1e-8..1e-4)

### Performance
- nbody_rk4_step(2): 166 ns → 188 ns (+13%, 2-body micro-benchmark; trade-off eliminates 6 heap allocations per step, net win at larger N)
- nbody_leapfrog_step(2): 58 ns → 54 ns (−6%)

## [0.1.0] - 2026-03-25

### Added
- **error** — `FalakError` enum with `InvalidParameter`, `MathError`, `ConvergenceError`, `EphemerisError`, `Io` variants
- **orbit** — `OrbitalElements` struct with validation, periapsis/apoapsis/semi-latus-rectum
- **kepler** — module stub for Kepler's laws, anomaly conversions
- **transfer** — module stub for Hohmann, bi-elliptic, Lambert transfers
- **perturbation** — module stub for J2, drag, SRP, third-body perturbations
- **nbody** — module stub for N-body gravitational simulation
- **ephemeris** — module stub for planetary positions and Julian date
- **frame** — module stub for reference frame transformations
- **maneuver** — module stub for delta-v and burn planning
- **logging** — structured logging via `FALAK_LOG` env var (feature-gated)
- Initial project scaffold with CI, benchmarks, and documentation structure
