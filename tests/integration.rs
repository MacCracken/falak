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
