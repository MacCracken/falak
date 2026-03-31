//! Cross-crate bridges — convert primitive values from other AGNOS science crates
//! into falak orbital mechanics parameters and vice versa.
//!
//! Always available — takes primitive values (f64), no science crate deps.

// ── Tara bridges (stellar astrophysics) ────────────────────────────────────

/// Convert stellar mass (kg) to gravitational parameter μ (m³/s²).
///
/// μ = G × M, where G = 6.674e-11 m³/(kg·s²).
#[must_use]
#[inline]
pub fn stellar_mass_to_mu(mass_kg: f64) -> f64 {
    const G: f64 = 6.674_30e-11;
    G * mass_kg
}

/// Convert luminosity (W) to habitable zone distance (AU).
///
/// Approximate: d_HZ ≈ sqrt(L / L_sun) AU.
/// L_sun = 3.828e26 W.
#[must_use]
#[inline]
pub fn luminosity_to_habitable_zone_au(luminosity_w: f64) -> f64 {
    const L_SUN: f64 = 3.828e26;
    if luminosity_w <= 0.0 {
        return 0.0;
    }
    (luminosity_w / L_SUN).sqrt()
}

// ── Impetus bridges (physics) ──────────────────────────────────────────────

/// Convert orbital velocity `[vx, vy, vz]` (m/s) and position `[x, y, z]` (m)
/// to gravitational force magnitude (N) for a given body mass and central mass.
///
/// F = G × M × m / r²
#[must_use]
pub fn orbital_to_gravity_force(
    position: [f64; 3],
    central_mass_kg: f64,
    body_mass_kg: f64,
) -> f64 {
    const G: f64 = 6.674_30e-11;
    let r2 = position[0] * position[0] + position[1] * position[1] + position[2] * position[2];
    if r2 < 1e-20 {
        return 0.0;
    }
    G * central_mass_kg * body_mass_kg / r2
}

/// Convert orbital velocity magnitude (m/s) and gravitational parameter μ (m³/s²)
/// to escape velocity threshold.
///
/// v_escape = sqrt(2μ/r). Returns the kinetic energy deficit (positive = bound).
#[must_use]
#[inline]
pub fn escape_energy_deficit(velocity_ms: f64, mu: f64, radius_m: f64) -> f64 {
    if radius_m <= 0.0 {
        return 0.0;
    }
    let v_escape = (2.0 * mu / radius_m).sqrt();
    0.5 * (v_escape * v_escape - velocity_ms * velocity_ms)
}

// ── Badal bridges (weather/climate) ────────────────────────────────────────

/// Convert solar distance (AU) and axial tilt (degrees) to seasonal
/// insolation factor (0.0–2.0, relative to Earth average).
///
/// Simplified: S ∝ 1/d² × (1 + tilt_effect × cos(season_angle)).
#[must_use]
#[inline]
pub fn solar_distance_to_insolation(
    distance_au: f64,
    axial_tilt_deg: f64,
    season_angle_rad: f64,
) -> f64 {
    if distance_au <= 0.0 {
        return 0.0;
    }
    let base = 1.0 / (distance_au * distance_au);
    let tilt_effect = (axial_tilt_deg.to_radians()).sin() * season_angle_rad.cos();
    base * (1.0 + 0.5 * tilt_effect)
}

/// Convert orbital eccentricity to climate variation amplitude.
///
/// Higher eccentricity → larger seasonal temperature swings.
/// Returns a scaling factor (1.0 at e=0, increasing with e).
#[must_use]
#[inline]
pub fn eccentricity_to_climate_variation(eccentricity: f64) -> f64 {
    let e = eccentricity.clamp(0.0, 0.99);
    // Flux ratio: (1+e)²/(1-e)² at perihelion vs aphelion
    let peri = (1.0 + e) * (1.0 + e);
    let aph = (1.0 - e) * (1.0 - e);
    peri / aph
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sun_mu() {
        // Sun: ~1.989e30 kg → μ ≈ 1.327e20
        let mu = stellar_mass_to_mu(1.989e30);
        assert!((mu - 1.327e20).abs() < 1e18);
    }

    #[test]
    fn habitable_zone_sun() {
        let d = luminosity_to_habitable_zone_au(3.828e26);
        assert!((d - 1.0).abs() < 0.01);
    }

    #[test]
    fn habitable_zone_zero() {
        assert!((luminosity_to_habitable_zone_au(0.0)).abs() < f64::EPSILON);
    }

    #[test]
    fn gravity_force_earth_surface() {
        // Earth: M=5.972e24, m=1kg, r=6.371e6 → F ≈ 9.82N
        let f = orbital_to_gravity_force([6.371e6, 0.0, 0.0], 5.972e24, 1.0);
        assert!((f - 9.82).abs() < 0.1);
    }

    #[test]
    fn escape_energy_bound() {
        // LEO: v=7.7km/s, μ_earth=3.986e14, r=6.671e6
        let deficit = escape_energy_deficit(7700.0, 3.986e14, 6.671e6);
        assert!(deficit > 0.0, "LEO should be bound: {deficit}");
    }

    #[test]
    fn insolation_earth() {
        let s = solar_distance_to_insolation(1.0, 23.44, 0.0);
        assert!(s > 0.9 && s < 1.2, "Earth insolation: {s}");
    }

    #[test]
    fn climate_variation_circular() {
        let v = eccentricity_to_climate_variation(0.0);
        assert!((v - 1.0).abs() < 0.001);
    }

    #[test]
    fn climate_variation_eccentric() {
        let v = eccentricity_to_climate_variation(0.5);
        assert!(v > 5.0, "eccentric orbit should have large variation: {v}");
    }
}
