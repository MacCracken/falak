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

criterion_group!(benches, orbit_creation);
criterion_main!(benches);
