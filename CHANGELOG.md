# Changelog

## [Unreleased]

### Added
- **kepler** — Kepler's equation solver (Newton-Raphson with Danby starter for elliptic, log starter for hyperbolic), all anomaly conversions (M↔E↔ν elliptic, M↔H↔ν hyperbolic), orbital period, mean motion, vis-viva velocity, orbital radius, `StateVector` type, bidirectional state vector ↔ orbital elements conversion
- **orbit** — Hyperbolic (e > 1) and parabolic (e = 1) orbit support in `OrbitalElements`, orbit type queries (`is_elliptical`, `is_circular`, `is_parabolic`, `is_hyperbolic`), correct periapsis/apoapsis for all orbit types
- **transfer** — Hohmann transfer, bi-elliptic transfer, plane change maneuver, phasing orbit; `HohmannTransfer`, `BiEllipticTransfer`, `PlaneChange` result types
- **frame** — Perifocal ↔ ECI transforms, ECI ↔ ECEF rotation, ECEF ↔ geodetic (WGS-84 Bowring's method), inertial ↔ rotating (synodic) frame; `Geodetic` type; WGS-84 constants
- **ephemeris** — Calendar ↔ Julian Date (Gregorian algorithm), Julian Date ↔ MJD, Julian Date ↔ Unix timestamp, GMST (IAU 1982), Julian centuries since J2000, day of year with leap year and month-day validation; time constants (`J2000_JD`, `UNIX_EPOCH_JD`, etc.)
- **perturbation** — J2 oblateness acceleration, J2 secular drift rates (RAAN and argument of periapsis), J3 acceleration, atmospheric drag (exponential model with `AtmosphereParams`), solar radiation pressure, third-body gravitational perturbation; physical constants (`J2_EARTH`, `J3_EARTH`, `R_EARTH`, `MU_SUN`, `SOLAR_PRESSURE_1AU`, `AU_METRES`)
- **maneuver** — `ImpulsiveBurn` type (prograde/retrograde/normal/radial constructors), delta-v budget summing, Tsiolkovsky rocket equation (propellant mass fraction ↔ max delta-v), escape/capture delta-v from circular orbit, Oberth effect factor
- **bridge** — Tara bridges (stellar mass → μ, luminosity → habitable zone), impetus bridges (gravity force, escape energy deficit), badal bridges (insolation, climate variation)
- **integration/soorat** — `OrbitPath::from_elements` (now returns `Result`), `PlanetaryPositions`, `CelestialBody`, `TransferTrajectory`, `GroundTrack` data types

### Changed
- **orbit** — `OrbitalElements::new` now accepts all orbit types: elliptical (e < 1, a > 0), parabolic (e = 1, p > 0), hyperbolic (e > 1, a < 0); periapsis returns p/2 for parabolic; apoapsis returns ∞ for open orbits
- **bridge** — Gravitational constant `G` now sourced from `kepler::G` (removed duplication)
- **integration/soorat** — `OrbitPath::from_elements` returns `Result` instead of silently masking errors; all structs marked `#[non_exhaustive]`; orbit path generation uses `kepler` module functions

### Fixed
- **ephemeris** — GMST formula was double-counting fractional day (T included full JD, then fractional day added again); now correctly separates 0h UT1 and fractional day
- **frame** — `eci_to_perifocal` was hardcoding W-component to 0.0, dropping out-of-plane ECI z-component; now computes full 3×3 transpose
- **orbit** — Parabolic periapsis was returning 0.0 (|a|×|1-e| = 0 when e=1); now returns p/2

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
