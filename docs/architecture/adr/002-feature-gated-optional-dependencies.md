# ADR-002: Feature-Gated Optional Dependencies

**Status**: Accepted
**Date**: 2026-03-25

## Context

Consumers may not need physics integration or logging overhead. A space sim that only converts orbital elements should not pay compile-time or binary-size costs for a full integration backend or structured logging infrastructure.

## Decision

Optional dependencies are placed behind Cargo feature gates:

- **`impetus`**: physics integration backend — gated behind the `impetus` feature
- **`tracing-subscriber`**: structured logging output — gated behind the `tracing-subscriber` feature
- **Core functionality** (Keplerian conversion, Lambert solver, propagation, perturbations) is always available with zero optional deps

The default feature set is empty. Consumers opt in to exactly what they need.

## Consequences

- Minimal consumers get a smaller, faster compile with no transitive dependency bloat
- CI must test all feature combinations (`--no-default-features`, `--all-features`, and individual gates) to prevent feature-gated code from silently breaking
- Adding a new optional dependency requires a new feature gate and a CI matrix entry
- Consumers (kiran, joshua, impetus) declare only the features they actually use in their `Cargo.toml`
