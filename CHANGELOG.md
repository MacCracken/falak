# Changelog

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
