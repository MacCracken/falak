# ADR-001: ECI as Universal Coordinate Convention

**Status**: Accepted
**Date**: 2026-03-31

## Context

Orbital mechanics requires a consistent coordinate system across all modules. Position and velocity vectors flow between Keplerian conversion, propagation, Lambert solving, and perturbation analysis. Without a single agreed-upon frame, every module boundary becomes a potential source of silent sign-flip or axis-swap bugs.

## Decision

All position/velocity vectors in falak use an ECI-like inertial reference frame:

- **X-axis**: points toward the vernal equinox
- **Z-axis**: points along the celestial pole (north)
- **Y-axis**: completes the right-handed coordinate system

`StateVector` and all propagation outputs are defined in this frame. The frame convention is explicitly documented on the `StateVector` struct and in module-level docs.

## Consequences

- Frame transforms are required when converting to ECEF, geodetic, or body-fixed coordinates
- `StateVector` documentation explicitly states the assumed reference frame, preventing misuse
- Consumers (kiran, joshua) that need ground-track or sensor-relative coordinates must apply their own rotation
- All internal math (Kepler, Lambert, Brouwer, n-body) operates in a single consistent frame with no hidden transforms
