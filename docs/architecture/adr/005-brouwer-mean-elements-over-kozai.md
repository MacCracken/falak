# ADR-005: Brouwer Mean Elements over Kozai

**Status**: Accepted
**Date**: 2026-03-31

## Context

Converting osculating orbital elements to mean elements is essential for long-term orbit prediction, station-keeping analysis, and filtering out short-period oscillations. Two main perturbation theories are widely used: Brouwer's theory and Kozai's theory. Both account for J2 oblateness, but differ in formulation and applicability.

## Decision

Use first-order Brouwer J2 theory with short-period removal:

- **Brouwer's method** removes short-period terms caused by J2, yielding mean elements that evolve smoothly under secular and long-period effects
- **First-order** in J2 is sufficient for the current accuracy requirements
- Short-period terms are computed and subtracted from osculating elements via closed-form expressions

## Consequences

- Valid for near-circular to moderate eccentricity orbits; breaks down for very high eccentricity or critical inclination (~63.4 deg) without additional corrections
- Secular rates of mean elements (RAAN drift, argument of perigee drift) match standard J2 drift formulas used throughout the astrodynamics community
- No atmospheric drag or higher-order gravitational harmonics (J3, J4, ...) are included — these can be added as future ADRs if consumers require them
- Kozai's method remains a viable alternative if critical-inclination orbits become important; switching would be a localized change in the mean-element module
