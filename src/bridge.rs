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
    crate::kepler::G * mass_kg
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

/// Convert position `[x, y, z]` (m) to gravitational force magnitude (N)
/// for a given body mass and central mass.
///
/// F = G × M × m / r²
#[must_use]
pub fn orbital_to_gravity_force(
    position: [f64; 3],
    central_mass_kg: f64,
    body_mass_kg: f64,
) -> f64 {
    let r2 = position[0] * position[0] + position[1] * position[1] + position[2] * position[2];
    if r2 < 1e-20 {
        return 0.0;
    }
    crate::kepler::G * central_mass_kg * body_mass_kg / r2
}

/// Compute the specific energy deficit relative to escape velocity (J/kg).
///
/// v_escape = sqrt(2μ/r). Returns 0.5 × (v_escape² − v²).
/// Positive → bound orbit, negative → escaping.
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

/// Convert stellar radius (m) and effective temperature (K) to luminosity (W).
///
/// L = 4π R² σ T⁴ (Stefan-Boltzmann law).
#[must_use]
#[inline]
pub fn stellar_luminosity(radius_m: f64, temperature_k: f64) -> f64 {
    const STEFAN_BOLTZMANN: f64 = 5.670_374e-8;
    if radius_m <= 0.0 || temperature_k <= 0.0 {
        return 0.0;
    }
    4.0 * std::f64::consts::PI * radius_m * radius_m * STEFAN_BOLTZMANN * temperature_k.powi(4)
}

/// Convert stellar spectral class temperature (K) to approximate mass (solar masses).
///
/// Empirical mass-luminosity-temperature relation for main sequence stars.
/// Returns mass in kg.
#[must_use]
#[inline]
pub fn temperature_to_mass_kg(temperature_k: f64) -> f64 {
    const M_SUN: f64 = 1.989e30;
    const T_SUN: f64 = 5778.0;
    if temperature_k <= 0.0 {
        return 0.0;
    }
    // M/M☉ ≈ (T/T☉)^(4/2.5) for main sequence (approximate)
    M_SUN * (temperature_k / T_SUN).powf(1.6)
}

// ── Impetus bridges — extended ───────────────────────────────────────────

/// Convert orbital state to gravitational force *vector* (N).
///
/// Unlike [`orbital_to_gravity_force`] which returns only the magnitude,
/// this returns the full `[Fx, Fy, Fz]` vector pointing toward the central body.
///
/// F⃗ = −G × M × m / r³ × r⃗
#[must_use]
pub fn orbital_to_gravity_force_vector(
    position: [f64; 3],
    central_mass_kg: f64,
    body_mass_kg: f64,
) -> [f64; 3] {
    let r2 = position[0] * position[0] + position[1] * position[1] + position[2] * position[2];
    if r2 < 1e-20 {
        return [0.0; 3];
    }
    let r = r2.sqrt();
    let factor = -crate::kepler::G * central_mass_kg * body_mass_kg / (r2 * r);
    [
        factor * position[0],
        factor * position[1],
        factor * position[2],
    ]
}

/// Convert orbital velocity and position to specific orbital energy (J/kg).
///
/// ε = v²/2 − μ/r (vis-viva). Negative → bound, positive → escaping.
#[must_use]
#[inline]
pub fn specific_orbital_energy(velocity_ms: f64, mu: f64, radius_m: f64) -> f64 {
    if radius_m <= 0.0 {
        return 0.0;
    }
    0.5 * velocity_ms * velocity_ms - mu / radius_m
}

// ── Badal bridges — extended ─────────────────────────────────────────────

/// Compute perihelion and aphelion distances from semi-major axis and eccentricity.
///
/// Returns `(perihelion_au, aphelion_au)` in AU. Useful for climate modelling
/// where the distance extremes drive temperature variation.
#[must_use]
#[inline]
pub fn orbital_distance_extremes_au(semi_major_axis_au: f64, eccentricity: f64) -> (f64, f64) {
    let perihelion = semi_major_axis_au * (1.0 - eccentricity);
    let aphelion = semi_major_axis_au * (1.0 + eccentricity);
    (perihelion, aphelion)
}

/// Compute Milankovitch obliquity forcing factor.
///
/// Maps axial tilt to a normalised forcing factor for climate models.
/// Higher obliquity → stronger seasonal contrast.
/// Returns a factor relative to Earth's current obliquity (23.44°).
#[must_use]
#[inline]
pub fn obliquity_forcing(axial_tilt_deg: f64) -> f64 {
    const EARTH_TILT: f64 = 23.44;
    if axial_tilt_deg <= 0.0 {
        return 0.0;
    }
    axial_tilt_deg.to_radians().sin() / EARTH_TILT.to_radians().sin()
}

/// Compute precession-driven insolation variation period (years).
///
/// The climatic precession period depends on the apsidal precession rate
/// and the axial precession rate. For Earth: ~21,000 years.
///
/// Simplified: T_prec ≈ 1 / |1/T_apsidal − 1/T_axial|
#[must_use]
#[inline]
pub fn precession_period_years(apsidal_period_yr: f64, axial_period_yr: f64) -> f64 {
    if apsidal_period_yr <= 0.0 || axial_period_yr <= 0.0 {
        return 0.0;
    }
    let rate_diff = (1.0 / apsidal_period_yr - 1.0 / axial_period_yr).abs();
    if rate_diff < 1e-30 {
        return f64::INFINITY;
    }
    1.0 / rate_diff
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

    // ── Extended bridges ─────────────────────────────────────────────

    #[test]
    fn stellar_luminosity_sun() {
        const R_SUN: f64 = 6.957e8;
        const T_SUN: f64 = 5778.0;
        let l = stellar_luminosity(R_SUN, T_SUN);
        assert!(
            (l - 3.828e26).abs() / 3.828e26 < 0.02,
            "Sun luminosity: {l:.3e}"
        );
    }

    #[test]
    fn temperature_to_mass_sun() {
        let m = temperature_to_mass_kg(5778.0);
        assert!((m - 1.989e30).abs() / 1.989e30 < 0.1, "Sun mass: {m:.3e}");
    }

    #[test]
    fn gravity_force_vector_direction() {
        let f = orbital_to_gravity_force_vector([7e6, 0.0, 0.0], 5.972e24, 1.0);
        assert!(f[0] < 0.0, "force should point inward: {}", f[0]);
        assert!(f[1].abs() < 1e-20);
        assert!(f[2].abs() < 1e-20);
        // Magnitude should match scalar version
        let mag = (f[0] * f[0] + f[1] * f[1] + f[2] * f[2]).sqrt();
        let scalar = orbital_to_gravity_force([7e6, 0.0, 0.0], 5.972e24, 1.0);
        assert!((mag - scalar).abs() < 1e-10);
    }

    #[test]
    fn specific_energy_bound() {
        // LEO is bound (negative energy)
        let e = specific_orbital_energy(7700.0, 3.986e14, 6.778e6);
        assert!(e < 0.0, "LEO should be bound: {e}");
    }

    #[test]
    fn orbital_distance_extremes() {
        let (peri, aph) = orbital_distance_extremes_au(1.0, 0.0167);
        assert!((peri - 0.983).abs() < 0.01, "Earth perihelion: {peri}");
        assert!((aph - 1.017).abs() < 0.01, "Earth aphelion: {aph}");
    }

    #[test]
    fn obliquity_forcing_earth() {
        let f = obliquity_forcing(23.44);
        assert!(
            (f - 1.0).abs() < 0.001,
            "Earth obliquity factor should be ~1.0: {f}"
        );
    }

    #[test]
    fn precession_period_earth() {
        // Earth: apsidal ~112,000 yr, axial ~25,770 yr
        // T = 1/|1/25770 - 1/112000| ≈ 33,500 yr
        let p = precession_period_years(112_000.0, 25_770.0);
        assert!(
            (p - 33_500.0).abs() < 1000.0,
            "precession period: {p:.0} yr, expected ~33500"
        );
    }
}
