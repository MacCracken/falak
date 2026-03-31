//! Earth-Moon Lagrange points and Jacobi constants.
//!
//! Demonstrates:
//! - Computing all five Lagrange points in the CR3BP
//! - Jacobi constant (energy integral) at each point
//! - Zero-velocity surface value at L1

use falak::cr3bp;

/// Earth-Moon mass ratio: μ = m_Moon / (m_Earth + m_Moon).
const MU_EM: f64 = 0.012_150_585;

fn main() -> Result<(), falak::FalakError> {
    let lp = cr3bp::lagrange_points(MU_EM)?;

    println!("=== Earth-Moon Lagrange Points (synodic frame) ===");
    for (name, pt) in [
        ("L1", lp.l1),
        ("L2", lp.l2),
        ("L3", lp.l3),
        ("L4", lp.l4),
        ("L5", lp.l5),
    ] {
        let state = [pt[0], pt[1], pt[2], 0.0, 0.0, 0.0];
        let cj = cr3bp::jacobi_constant(&state, MU_EM);
        println!("{name}: ({:+.6}, {:+.6})  C_J = {cj:.6}", pt[0], pt[1]);
    }

    // Zero-velocity surface at L1 (accessible region boundary)
    let zvs = cr3bp::zero_velocity_value(lp.l1[0], lp.l1[1], lp.l1[2], MU_EM);
    println!("\nZero-velocity surface at L1: 2*Omega = {zvs:.6}");
    println!("Particles with C_J > {zvs:.4} cannot cross between Earth and Moon.");

    Ok(())
}
