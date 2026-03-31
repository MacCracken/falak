# ADR-003: Dependency-Free Bridges

**Status**: Accepted
**Date**: 2026-03-25

## Context

falak must not depend on downstream crates (tara, badal) or create circular dependencies. The AGNOS ecosystem pattern (established by garjan and bodh) requires that bridge functions remain dependency-free so that science crates can be composed without locking the dependency graph.

## Decision

Bridge functions in falak:

- Accept **only f64 primitives and fixed-size arrays** — never import external types from upstream or downstream crates
- Return library-owned types (`StateVector`, `OrbitalElements`, `f64`) or primitives
- The **orchestration layer** (kiran, joshua) is responsible for calling other crates and passing primitive values into falak

## Consequences

- falak has zero compile-time coupling to tara, badal, or any other AGNOS science crate
- Bridge functions can be tested entirely with synthetic numerical data — no astronomy or atmosphere crate needed in the test harness
- Adding new upstream data sources (e.g., atmospheric drag from badal) requires no API changes in falak — just new bridge functions accepting f64
- The orchestration layer (kiran, joshua) is the only component that imports all crates simultaneously
