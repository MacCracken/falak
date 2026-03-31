# ADR-004: Universal Variables for Lambert Solver

**Status**: Accepted
**Date**: 2026-03-31

## Context

Lambert's problem (find the orbit connecting two position vectors in a given time of flight) can be solved via several methods: Battin's method, Izzo's algorithm, or universal variables with Stumpff functions. The choice affects code complexity, robustness across conic types, and performance.

## Decision

Use Stumpff functions with Newton iteration, following Bate/Mueller/White Chapter 5:

- **Stumpff C(z) and S(z)** provide a single formulation that handles elliptic (z > 0), parabolic (z = 0), and hyperbolic (z < 0) orbits without branching
- **Newton-Raphson iteration** converges to the universal variable for the given time of flight
- **Derivatives** are computed via finite differences rather than analytical expressions

## Consequences

- A single code path handles all conic sections — no separate elliptic/hyperbolic solvers to maintain
- Performance is 97 ns per solve, sufficient for real-time trajectory planning
- Finite-difference derivatives are simpler to implement and verify than analytical Stumpff derivatives, at the cost of one extra function evaluation per iteration
- The solver is well-documented in standard references (BMW Ch. 5, Vallado Ch. 5), making it auditable
- Multi-revolution solutions require additional logic on top of the base universal-variable solver
