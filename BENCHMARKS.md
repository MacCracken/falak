# Benchmarks

> Last run: 2026-03-31T23:06:07Z | Commit: `359cb3d` | Branch: `main` | Samples: 100

| Benchmark | Time | Unit |
|-----------|------|------|
| `OrbitalElements::new` | 4.5 ns | ns |
| `solve_kepler_elliptic(e=0.5)` | 109.3 ns | ns |
| `solve_kepler_elliptic(e=0.99)` | 244.6 ns | ns |
| `solve_kepler_hyperbolic(e=1.5)` | 186.8 ns | ns |
| `mean_to_true_anomaly(e=0.3)` | 148.7 ns | ns |
| `elements_to_state` | 41.1 ns | ns |
| `state_to_elements` | 57.1 ns | ns |
| `hohmann(LEO→GEO)` | 13.2 ns | ns |
| `bi_elliptic` | 25.3 ns | ns |
| `plane_change` | 15.2 ns | ns |
| `lambert(90°)` | 123.4 ns | ns |
| `perifocal_to_eci` | 0.6 ns | ps |
| `ecef_to_geodetic` | 220.9 ns | ns |
| `calendar_to_jd` | 22.0 ns | ns |
| `gmst` | 0.6 ns | ps |
| `propagate::kepler(1_period)` | 105.6 ns | ns |
