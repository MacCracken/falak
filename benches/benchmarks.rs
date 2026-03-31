use criterion::{Criterion, criterion_group, criterion_main};

fn orbit_creation(c: &mut Criterion) {
    c.bench_function("OrbitalElements::new", |b| {
        b.iter(|| {
            let _ = std::hint::black_box(falak::orbit::OrbitalElements::new(
                7000.0e3, 0.01, 0.9, 1.0, 0.5, 0.0,
            ));
        });
    });
}

fn kepler_elliptic(c: &mut Criterion) {
    c.bench_function("solve_kepler_elliptic(e=0.5)", |b| {
        b.iter(|| {
            let _ = std::hint::black_box(falak::kepler::solve_kepler_elliptic(1.5, 0.5));
        });
    });

    c.bench_function("solve_kepler_elliptic(e=0.99)", |b| {
        b.iter(|| {
            let _ = std::hint::black_box(falak::kepler::solve_kepler_elliptic(0.1, 0.99));
        });
    });
}

fn kepler_hyperbolic(c: &mut Criterion) {
    c.bench_function("solve_kepler_hyperbolic(e=1.5)", |b| {
        b.iter(|| {
            let _ = std::hint::black_box(falak::kepler::solve_kepler_hyperbolic(2.0, 1.5));
        });
    });
}

fn anomaly_conversions(c: &mut Criterion) {
    c.bench_function("mean_to_true_anomaly(e=0.3)", |b| {
        b.iter(|| {
            let _ = std::hint::black_box(falak::kepler::mean_to_true_anomaly(1.0, 0.3));
        });
    });
}

fn state_vector_roundtrip(c: &mut Criterion) {
    let elements = falak::orbit::OrbitalElements::new(7e6, 0.1, 0.5, 1.0, 0.5, 0.8).unwrap();
    let mu = 3.986_004_418e14;

    c.bench_function("elements_to_state", |b| {
        b.iter(|| {
            let _ = std::hint::black_box(falak::kepler::elements_to_state(&elements, mu));
        });
    });

    let state = falak::kepler::elements_to_state(&elements, mu).unwrap();
    c.bench_function("state_to_elements", |b| {
        b.iter(|| {
            let _ = std::hint::black_box(falak::kepler::state_to_elements(&state, mu));
        });
    });
}

fn transfer_functions(c: &mut Criterion) {
    let mu = 3.986_004_418e14;
    let r_leo = 6.671e6;
    let r_geo = 42_164.0e3;

    c.bench_function("hohmann(LEO→GEO)", |b| {
        b.iter(|| {
            let _ = std::hint::black_box(falak::transfer::hohmann(r_leo, r_geo, mu));
        });
    });

    c.bench_function("bi_elliptic", |b| {
        b.iter(|| {
            let _ =
                std::hint::black_box(falak::transfer::bi_elliptic(r_leo, r_geo, r_geo * 3.0, mu));
        });
    });

    c.bench_function("plane_change", |b| {
        b.iter(|| {
            let _ =
                std::hint::black_box(falak::transfer::plane_change(7700.0, 0.5_f64.to_radians()));
        });
    });
}

fn frame_transforms(c: &mut Criterion) {
    c.bench_function("perifocal_to_eci", |b| {
        b.iter(|| {
            let _ = std::hint::black_box(falak::frame::perifocal_to_eci(
                [7e6, 3e6, 0.0],
                1.2,
                0.8,
                0.5,
            ));
        });
    });

    c.bench_function("ecef_to_geodetic", |b| {
        b.iter(|| {
            let _ = std::hint::black_box(falak::frame::ecef_to_geodetic([
                4_000_000.0,
                3_000_000.0,
                4_500_000.0,
            ]));
        });
    });
}

fn ephemeris_functions(c: &mut Criterion) {
    c.bench_function("calendar_to_jd", |b| {
        b.iter(|| {
            let _ = std::hint::black_box(falak::ephemeris::calendar_to_jd(2024, 7, 15.75));
        });
    });

    c.bench_function("gmst", |b| {
        b.iter(|| {
            let _ = std::hint::black_box(falak::ephemeris::gmst(2_451_545.0));
        });
    });
}

fn perturbation_functions(c: &mut Criterion) {
    let pos = [falak::perturbation::R_EARTH + 400e3, 0.0, 0.0];

    c.bench_function("j2_acceleration", |b| {
        b.iter(|| {
            let _ = std::hint::black_box(falak::perturbation::j2_acceleration(
                pos,
                3.986e14,
                falak::perturbation::J2_EARTH,
                falak::perturbation::R_EARTH,
            ));
        });
    });
}

fn maneuver_functions(c: &mut Criterion) {
    c.bench_function("propellant_mass_fraction", |b| {
        b.iter(|| {
            let _ = std::hint::black_box(falak::maneuver::propellant_mass_fraction(3000.0, 300.0));
        });
    });

    c.bench_function("escape_delta_v", |b| {
        b.iter(|| {
            let _ =
                std::hint::black_box(falak::maneuver::escape_delta_v(6.671e6, 3.986e14, 3000.0));
        });
    });
}

fn bridge_functions(c: &mut Criterion) {
    c.bench_function("stellar_mass_to_mu", |b| {
        b.iter(|| {
            let _ = std::hint::black_box(falak::bridge::stellar_mass_to_mu(1.989e30));
        });
    });

    c.bench_function("orbital_to_gravity_force", |b| {
        b.iter(|| {
            let _ = std::hint::black_box(falak::bridge::orbital_to_gravity_force(
                [6.371e6, 0.0, 0.0],
                5.972e24,
                1.0,
            ));
        });
    });
}

criterion_group!(
    benches,
    orbit_creation,
    kepler_elliptic,
    kepler_hyperbolic,
    anomaly_conversions,
    state_vector_roundtrip,
    transfer_functions,
    frame_transforms,
    ephemeris_functions,
    perturbation_functions,
    maneuver_functions,
    bridge_functions,
);
criterion_main!(benches);
