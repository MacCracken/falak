//! Sun-Earth-Moon N-body simulation with energy conservation check.
//!
//! Demonstrates:
//! - Setting up a three-body system with canonical gravitational parameters
//! - Evolving with the symplectic leapfrog integrator
//! - Verifying energy conservation over one orbital period

use falak::nbody::{Body, Integrator, System, evolve};

const MU_SUN: f64 = 1.327_124_4e20; // m^3/s^2
const MU_EARTH: f64 = 3.986_004_418e14;
const MU_MOON: f64 = 4.902_800e12;
const AU: f64 = 1.496e11; // metres
const YEAR_S: f64 = 365.25 * 86_400.0; // seconds

fn main() -> Result<(), falak::FalakError> {
    let sun = Body::with_mu([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.989e30, MU_SUN);
    let v_earth = (MU_SUN / AU).sqrt(); // circular orbit speed
    let earth = Body::with_mu([AU, 0.0, 0.0], [0.0, v_earth, 0.0], 5.972e24, MU_EARTH);
    let r_moon = AU + 3.844e8;
    let v_moon = v_earth + (MU_EARTH / 3.844e8).sqrt();
    let moon = Body::with_mu([r_moon, 0.0, 0.0], [0.0, v_moon, 0.0], 7.342e22, MU_MOON);

    let mut sys = System::new(vec![sun, earth, moon], 0.0)?;
    let e0 = sys.total_energy();
    println!("=== Sun-Earth-Moon N-Body (1 year, leapfrog) ===");
    println!("Initial energy: {e0:.6e} J");

    evolve(&mut sys, YEAR_S, 3600.0, Integrator::Leapfrog)?;

    let e1 = sys.total_energy();
    let drift = ((e1 - e0) / e0).abs();
    println!("Final energy:   {e1:.6e} J");
    println!("Relative drift: {drift:.2e}");
    println!(
        "Earth final pos: ({:.4e}, {:.4e}) m",
        sys.bodies[1].position[0], sys.bodies[1].position[1]
    );

    Ok(())
}
