#[test]
fn error_type_exists() {
    let err = falak::error::FalakError::InvalidParameter("test".into());
    let msg = format!("{err}");
    assert!(msg.contains("test"));
}

#[test]
fn valid_orbit_creation() {
    let orbit = falak::orbit::OrbitalElements::new(7000.0e3, 0.01, 0.9, 1.0, 0.5, 0.0);
    assert!(orbit.is_ok());
}

#[test]
fn invalid_orbit_rejected() {
    let orbit = falak::orbit::OrbitalElements::new(-1.0, 0.5, 0.0, 0.0, 0.0, 0.0);
    assert!(orbit.is_err());
}

// ── Cross-module integration: combined perturbations ─────────────────

const MU_EARTH: f64 = 3.986_004_418e14;

/// Propagate LEO orbit with J2 + drag + lunar third-body combined.
///
/// Verifies that Cowell's method handles multiple simultaneous perturbations
/// and that the orbit remains physically reasonable over ~1 orbit.
#[test]
fn cowell_j2_drag_third_body_combined() {
    use falak::kepler;
    use falak::perturbation;
    use falak::propagate;

    // ISS-like orbit: 400 km, low eccentricity, 51.6° inclination
    let elements =
        falak::orbit::OrbitalElements::new(6_778e3, 0.001, 51.6_f64.to_radians(), 0.0, 0.0, 0.0)
            .unwrap();
    let state = kepler::elements_to_state(&elements, MU_EARTH).unwrap();
    let period = kepler::orbital_period(elements.semi_major_axis, MU_EARTH).unwrap();

    // Moon position (simplified: along +x at mean distance)
    let moon_pos = [384_400e3, 0.0, 0.0];
    let mu_moon = 4.902_800e12;

    // Drag parameters for ISS-like spacecraft
    let atmo = perturbation::AtmosphereParams::new(
        2.2,
        0.01, // ~100 kg/m²
        1.225,
        0.0,
        8500.0,
        perturbation::R_EARTH,
    );

    let combined_perturb = move |pos: [f64; 3], vel: [f64; 3], _t: f64| -> [f64; 3] {
        let j2 = perturbation::j2_acceleration(
            pos,
            MU_EARTH,
            perturbation::J2_EARTH,
            perturbation::R_EARTH,
        );
        let j3 = perturbation::j3_acceleration(
            pos,
            MU_EARTH,
            perturbation::J3_EARTH,
            perturbation::R_EARTH,
        );
        let drag = perturbation::drag_acceleration(pos, vel, &atmo);
        let lunar = perturbation::third_body_acceleration(pos, moon_pos, mu_moon);
        [
            j2[0] + j3[0] + drag[0] + lunar[0],
            j2[1] + j3[1] + drag[1] + lunar[1],
            j2[2] + j3[2] + drag[2] + lunar[2],
        ]
    };

    // Propagate for one full orbit
    let result =
        propagate::cowell(&state, MU_EARTH, period, 10.0, Some(&combined_perturb)).unwrap();

    // Orbit should remain roughly at same altitude (not crashed or escaped)
    let r_final =
        (result.position[0].powi(2) + result.position[1].powi(2) + result.position[2].powi(2))
            .sqrt();
    let alt_final = r_final - perturbation::R_EARTH;
    assert!(
        alt_final > 300e3 && alt_final < 500e3,
        "altitude should stay near 400 km, got {:.1} km",
        alt_final / 1e3
    );

    // Velocity should remain orbital (not crashed or escaped)
    let v_final =
        (result.velocity[0].powi(2) + result.velocity[1].powi(2) + result.velocity[2].powi(2))
            .sqrt();
    assert!(
        v_final > 7000.0 && v_final < 8000.0,
        "velocity should stay near 7.7 km/s, got {:.1} m/s",
        v_final
    );

    // Perturbed orbit should NOT close perfectly (unlike pure Kepler)
    let dr = [
        result.position[0] - state.position[0],
        result.position[1] - state.position[1],
        result.position[2] - state.position[2],
    ];
    let dr_mag = (dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]).sqrt();
    // J2 causes ~km-scale drift per orbit; pure Kepler would close within metres
    assert!(
        dr_mag > 100.0,
        "perturbed orbit should not close perfectly, drift = {dr_mag:.1} m"
    );
}

/// Verify Encke and Cowell agree when using combined perturbations.
#[test]
fn encke_vs_cowell_combined_perturbations() {
    use falak::kepler;
    use falak::perturbation;
    use falak::propagate;

    let elements = falak::orbit::OrbitalElements::new(7_000e3, 0.01, 0.9, 1.0, 0.5, 0.0).unwrap();
    let state = kepler::elements_to_state(&elements, MU_EARTH).unwrap();

    let moon_pos = [384_400e3, 0.0, 0.0];
    let mu_moon = 4.902_800e12;

    let combined = move |pos: [f64; 3], _vel: [f64; 3], _t: f64| -> [f64; 3] {
        let j2 = perturbation::j2_acceleration(
            pos,
            MU_EARTH,
            perturbation::J2_EARTH,
            perturbation::R_EARTH,
        );
        let lunar = perturbation::third_body_acceleration(pos, moon_pos, mu_moon);
        [j2[0] + lunar[0], j2[1] + lunar[1], j2[2] + lunar[2]]
    };

    let dt_prop = 500.0;

    let cowell_result = propagate::cowell(&state, MU_EARTH, dt_prop, 1.0, Some(&combined)).unwrap();
    let encke_result = propagate::encke(&elements, MU_EARTH, dt_prop, 1.0, &combined).unwrap();

    // Should agree within ~1 km over 500s
    for i in 0..3 {
        let diff = (cowell_result.position[i] - encke_result.position[i]).abs();
        assert!(
            diff < 1000.0,
            "position[{i}] Cowell vs Encke diff: {diff:.1} m"
        );
    }
}
