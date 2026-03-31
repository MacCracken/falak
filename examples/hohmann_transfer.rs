//! Hohmann transfer example — LEO to GEO.
//!
//! Demonstrates:
//! - Computing a Hohmann transfer between circular orbits
//! - Using the maneuver module for delta-v budgeting
//! - Tsiolkovsky rocket equation for propellant mass

use falak::{maneuver, perturbation, transfer};

const MU_EARTH: f64 = 3.986_004_418e14; // m³/s²

fn main() -> Result<(), falak::FalakError> {
    let r_leo = perturbation::R_EARTH + 400e3; // 400 km LEO
    let r_geo = 42_164e3; // GEO radius

    // Compute Hohmann transfer
    let hohmann = transfer::hohmann(r_leo, r_geo, MU_EARTH)?;
    println!("=== LEO → GEO Hohmann Transfer ===");
    println!("Δv₁ (departure burn):  {:.1} m/s", hohmann.delta_v1);
    println!("Δv₂ (arrival burn):    {:.1} m/s", hohmann.delta_v2);
    println!("Total Δv:              {:.1} m/s", hohmann.total_delta_v);
    println!(
        "Transfer time:         {:.1} hours",
        hohmann.time_of_flight / 3600.0
    );

    // Build a maneuver plan
    let mut plan = maneuver::ManeuverPlan::new();
    plan.add_burn(0.0, maneuver::ImpulsiveBurn::prograde(hohmann.delta_v1));
    plan.add_burn(
        hohmann.time_of_flight,
        maneuver::ImpulsiveBurn::prograde(hohmann.delta_v2),
    );
    println!(
        "\nManeuver plan: {} burns, total Δv = {:.1} m/s",
        plan.len(),
        plan.total_delta_v()
    );

    // Propellant budget (Tsiolkovsky)
    let isp = 300.0; // bipropellant engine
    let exhaust_velocity = isp * 9.80665; // m/s
    let prop_fraction =
        maneuver::propellant_mass_fraction(hohmann.total_delta_v, exhaust_velocity)?;
    println!(
        "\nWith Isp = {isp} s: propellant = {:.1}% of initial mass",
        prop_fraction * 100.0
    );

    // Compare: bi-elliptic transfer (only better for r_final/r_initial > 11.94)
    let r_intermediate = r_geo * 3.0; // swing out to 3× GEO
    let bielliptic = transfer::bi_elliptic(r_leo, r_geo, r_intermediate, MU_EARTH)?;
    println!(
        "\n=== Bi-Elliptic (r_int = {:.0} km) ===",
        r_intermediate / 1e3
    );
    println!("Total Δv:      {:.1} m/s", bielliptic.total_delta_v);
    println!(
        "Transfer time: {:.1} hours",
        bielliptic.time_of_flight / 3600.0
    );
    println!(
        "Bi-elliptic is {} than Hohmann by {:.1} m/s",
        if bielliptic.total_delta_v < hohmann.total_delta_v {
            "cheaper"
        } else {
            "more expensive"
        },
        (bielliptic.total_delta_v - hohmann.total_delta_v).abs()
    );

    // Plane change cost (e.g., changing inclination at GEO)
    let v_geo = (MU_EARTH / r_geo).sqrt();
    let plane_change = transfer::plane_change(v_geo, 28.5_f64.to_radians())?;
    println!(
        "\n=== GEO Plane Change (28.5°) ===\nΔv = {:.1} m/s",
        plane_change.delta_v
    );

    Ok(())
}
