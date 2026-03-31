//! Orbit propagation example — LEO satellite with J2 perturbation.
//!
//! Demonstrates:
//! - Creating orbital elements and converting to state vectors
//! - Analytic (Kepler) propagation
//! - Numerical (Cowell) propagation with J2 perturbation
//! - Computing orbital period and velocity

use falak::{kepler, orbit, perturbation, propagate};

const MU_EARTH: f64 = 3.986_004_418e14; // m³/s²

fn main() -> Result<(), falak::FalakError> {
    // ISS-like orbit: 400 km altitude, near-circular, 51.6° inclination
    let altitude = 400e3; // metres
    let sma = perturbation::R_EARTH + altitude;
    let elements = orbit::OrbitalElements::new(
        sma,
        0.001,                 // near-circular
        51.6_f64.to_radians(), // ISS inclination
        0.0,                   // RAAN
        0.0,                   // argument of periapsis
        0.0,                   // true anomaly (start at periapsis)
    )?;

    // Orbital period and velocity
    let period = kepler::orbital_period(sma, MU_EARTH)?;
    let velocity = kepler::vis_viva(sma, sma, MU_EARTH)?;
    println!("Orbital period: {:.1} min", period / 60.0);
    println!("Circular velocity: {:.1} m/s", velocity);

    // Convert to state vector (ECI)
    let state = kepler::elements_to_state(&elements, MU_EARTH)?;
    println!(
        "Initial position: [{:.0}, {:.0}, {:.0}] m",
        state.position[0], state.position[1], state.position[2]
    );

    // Analytic propagation (exact, no perturbations)
    let half_period = propagate::kepler(&elements, MU_EARTH, period / 2.0)?;
    println!(
        "After half orbit, true anomaly: {:.2}°",
        half_period.true_anomaly.to_degrees()
    );

    // Numerical propagation with J2 (shows orbital plane precession)
    let j2_perturb = |pos: [f64; 3], _vel: [f64; 3], _t: f64| -> [f64; 3] {
        perturbation::j2_acceleration(pos, MU_EARTH, perturbation::J2_EARTH, perturbation::R_EARTH)
    };

    let perturbed = propagate::cowell(&state, MU_EARTH, period, 10.0, Some(&j2_perturb))?;
    let r_final = (perturbed.position[0].powi(2)
        + perturbed.position[1].powi(2)
        + perturbed.position[2].powi(2))
    .sqrt();
    println!(
        "After 1 orbit (J2): altitude = {:.1} km",
        (r_final - perturbation::R_EARTH) / 1e3
    );

    // J2 secular rates (RAAN and arg-periapsis drift)
    let (d_raan, d_omega) = perturbation::j2_secular_rates(
        sma,
        0.001,
        51.6_f64.to_radians(),
        MU_EARTH,
        perturbation::J2_EARTH,
        perturbation::R_EARTH,
    )?;
    println!("RAAN drift: {:.4}°/day", d_raan.to_degrees() * 86400.0);
    println!(
        "Arg periapsis drift: {:.4}°/day",
        d_omega.to_degrees() * 86400.0
    );

    Ok(())
}
