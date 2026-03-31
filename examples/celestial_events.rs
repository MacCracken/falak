//! Celestial events — planetary positions, eclipses, sunrise/sunset.
//!
//! Demonstrates:
//! - Julian date conversion and planetary ephemeris
//! - Eclipse check for a satellite in LEO
//! - Sun rise/transit/set for an observer on Earth

use falak::ephemeris::{
    self, EclipseState, Planet, calendar_to_jd, eclipse_cylindrical, planetary_position,
    rise_transit_set,
};

fn main() -> Result<(), falak::FalakError> {
    // J2000.0 epoch: 2000-01-01 12:00 TT
    let jd = calendar_to_jd(2000, 1, 1.5)?;
    println!("=== Planetary Positions at J2000 ===");
    for (name, planet) in [
        ("Mercury", Planet::Mercury),
        ("Venus", Planet::Venus),
        ("Mars", Planet::Mars),
        ("Jupiter", Planet::Jupiter),
    ] {
        let pos = planetary_position(planet, jd);
        println!(
            "{name:8}: lon = {:6.2} deg, dist = {:.3} AU",
            pos.longitude.to_degrees(),
            pos.distance,
        );
    }

    // Eclipse check: ISS-like orbit on Earth's night side
    let r_earth = 6_378_137.0_f64;
    let sat_pos = [-(r_earth + 400e3), 0.0, 0.0]; // behind Earth
    let sun_pos = [1.496e11, 0.0, 0.0]; // Sun at +x
    let ecl = eclipse_cylindrical(sat_pos, sun_pos, r_earth);
    let state_str = match ecl.state {
        EclipseState::Sunlit => "sunlit",
        EclipseState::Penumbra => "penumbra",
        EclipseState::Umbra => "umbra",
        _ => "unknown",
    };
    println!("\n=== Eclipse Check (ISS night side) ===");
    println!(
        "State: {state_str}, shadow fraction: {:.2}",
        ecl.shadow_fraction
    );

    // Sunrise/sunset in London (51.5 N, 0.0 W) on 2025-06-21
    let jd_day = calendar_to_jd(2025, 6, 21.0)?;
    let sun = planetary_position(Planet::Earth, jd_day);
    let sun_lon = sun.longitude + std::f64::consts::PI; // geocentric
    // Approximate Sun RA/Dec (ecliptic lon -> equatorial, obliquity 23.44 deg)
    let obliq = 23.44_f64.to_radians();
    let ra = (sun_lon.sin() * obliq.cos()).atan2(sun_lon.cos());
    let dec = (obliq.sin() * sun_lon.sin()).asin();
    let lat = 51.5_f64.to_radians();
    let lon = 0.0_f64.to_radians();
    let h0 = (-0.8333_f64).to_radians(); // standard refraction
    let rts = rise_transit_set(jd_day, ra, dec, lat, lon, h0)?;

    println!("\n=== Sun Rise/Set — London, 2025-06-21 ===");
    if let Some(rise_jd) = rts.rise {
        let (_, _, d) = ephemeris::jd_to_calendar(rise_jd);
        println!(
            "Sunrise:  ~{:02.0}:{:02.0} UT",
            (d.fract() * 24.0).floor(),
            ((d.fract() * 24.0).fract() * 60.0).floor()
        );
    }
    if let Some(set_jd) = rts.set {
        let (_, _, d) = ephemeris::jd_to_calendar(set_jd);
        println!(
            "Sunset:   ~{:02.0}:{:02.0} UT",
            (d.fract() * 24.0).floor(),
            ((d.fract() * 24.0).fract() * 60.0).floor()
        );
    }

    Ok(())
}
