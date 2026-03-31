//! Ephemeris computation — Julian date, sidereal time, and time conversions.
//!
//! Provides the time foundation for orbital mechanics: Julian dates,
//! Modified Julian dates, and Greenwich Mean Sidereal Time.

use tracing::instrument;

use crate::error::{FalakError, Result};

/// Julian date of the J2000.0 epoch (2000 January 1, 12:00 TT).
pub const J2000_JD: f64 = 2_451_545.0;

/// Julian date of the Unix epoch (1970 January 1, 00:00 UTC).
pub const UNIX_EPOCH_JD: f64 = 2_440_587.5;

/// Modified Julian Date offset (MJD = JD − 2400000.5).
pub const MJD_OFFSET: f64 = 2_400_000.5;

/// Seconds per Julian day.
pub const SECONDS_PER_DAY: f64 = 86_400.0;

/// Julian century in days.
pub const DAYS_PER_JULIAN_CENTURY: f64 = 36_525.0;

// ── Julian Date conversions ───────────────────────────────────────────────

/// Convert a calendar date to Julian Date.
///
/// Uses the algorithm valid for dates after the Gregorian calendar reform (1582).
///
/// # Arguments
///
/// * `year` — Calendar year (e.g. 2024)
/// * `month` — Month (1–12)
/// * `day` — Day of month (1–31, may include fractional part for time)
///
/// # Errors
///
/// Returns [`FalakError::InvalidParameter`] if month is out of range.
#[must_use = "returns the computed Julian Date"]
#[instrument(level = "trace")]
pub fn calendar_to_jd(year: i32, month: u32, day: f64) -> Result<f64> {
    if !(1..=12).contains(&month) {
        return Err(FalakError::InvalidParameter(
            format!("month must be 1–12, got {month}").into(),
        ));
    }

    let (y, m) = if month <= 2 {
        (year as f64 - 1.0, month as f64 + 12.0)
    } else {
        (year as f64, month as f64)
    };

    let a = (y / 100.0).floor();
    let b = 2.0 - a + (a / 4.0).floor();

    Ok((365.25 * (y + 4716.0)).floor() + (30.6001 * (m + 1.0)).floor() + day + b - 1524.5)
}

/// Convert Julian Date to calendar date.
///
/// Returns `(year, month, day)` where day may include a fractional part.
#[must_use]
pub fn jd_to_calendar(jd: f64) -> (i32, u32, f64) {
    let jd = jd + 0.5;
    let z = jd.floor();
    let f = jd - z;

    let a = if z < 2_299_161.0 {
        z
    } else {
        let alpha = ((z - 1_867_216.25) / 36_524.25).floor();
        z + 1.0 + alpha - (alpha / 4.0).floor()
    };

    let b = a + 1524.0;
    let c = ((b - 122.1) / 365.25).floor();
    let d = (365.25 * c).floor();
    let e = ((b - d) / 30.6001).floor();

    let day = b - d - (30.6001 * e).floor() + f;
    let month = if e < 14.0 { e - 1.0 } else { e - 13.0 };
    let year = if month > 2.0 { c - 4716.0 } else { c - 4715.0 };

    (year as i32, month as u32, day)
}

/// Convert Julian Date to Modified Julian Date.
#[must_use]
#[inline]
pub fn jd_to_mjd(jd: f64) -> f64 {
    jd - MJD_OFFSET
}

/// Convert Modified Julian Date to Julian Date.
#[must_use]
#[inline]
pub fn mjd_to_jd(mjd: f64) -> f64 {
    mjd + MJD_OFFSET
}

/// Convert Unix timestamp (seconds since 1970-01-01 00:00 UTC) to Julian Date.
#[must_use]
#[inline]
pub fn unix_to_jd(unix_seconds: f64) -> f64 {
    UNIX_EPOCH_JD + unix_seconds / SECONDS_PER_DAY
}

/// Convert Julian Date to Unix timestamp.
#[must_use]
#[inline]
pub fn jd_to_unix(jd: f64) -> f64 {
    (jd - UNIX_EPOCH_JD) * SECONDS_PER_DAY
}

/// Julian centuries since J2000.0.
#[must_use]
#[inline]
pub fn julian_centuries_since_j2000(jd: f64) -> f64 {
    (jd - J2000_JD) / DAYS_PER_JULIAN_CENTURY
}

// ── Sidereal Time ─────────────────────────────────────────────────────────

/// Compute Greenwich Mean Sidereal Time (GMST) in radians.
///
/// Uses the IAU 1982 expression for GMST at 0h UT1, then adds the
/// fractional UT1 day scaled by the sidereal/solar ratio.
///
/// # Arguments
///
/// * `jd_ut1` — Julian Date in UT1 time scale
#[must_use]
#[inline]
pub fn gmst(jd_ut1: f64) -> f64 {
    // Separate into 0h UT1 and fractional day
    let jd_0h = (jd_ut1 + 0.5).floor() - 0.5;
    let frac_day = jd_ut1 - jd_0h;

    // Julian centuries from J2000.0 to 0h UT1 (NOT the full JD)
    let t0 = (jd_0h - J2000_JD) / DAYS_PER_JULIAN_CENTURY;

    // GMST at 0h UT1 in seconds (IAU 1982)
    let gmst_0h_sec =
        24_110.548_41 + 8_640_184.812_866 * t0 + 0.093_104 * t0 * t0 - 6.2e-6 * t0 * t0 * t0;

    // Add fractional day scaled by sidereal/solar ratio
    let gmst_total_sec = gmst_0h_sec + frac_day * SECONDS_PER_DAY * 1.002_737_909_350_795;

    // Convert to radians and normalise to [0, 2π)
    let gmst_rad = (gmst_total_sec / SECONDS_PER_DAY) * std::f64::consts::TAU;
    gmst_rad.rem_euclid(std::f64::consts::TAU)
}

// ── Day of year ───────────────────────────────────────────────────────────

/// Compute the day of year (1–366) from calendar date.
///
/// # Errors
///
/// Returns [`FalakError::InvalidParameter`] if month or day is out of range.
#[must_use = "returns the computed day of year"]
#[instrument(level = "trace")]
pub fn day_of_year(year: i32, month: u32, day: u32) -> Result<u32> {
    if !(1..=12).contains(&month) {
        return Err(FalakError::InvalidParameter(
            format!("month must be 1–12, got {month}").into(),
        ));
    }
    let is_leap = (year % 4 == 0 && year % 100 != 0) || year % 400 == 0;
    let days_in_month: [u32; 12] = if is_leap {
        [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    } else {
        [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    };

    let max_day = days_in_month[month as usize - 1];
    if day == 0 || day > max_day {
        return Err(FalakError::InvalidParameter(
            format!("day must be 1–{max_day} for month {month}, got {day}").into(),
        ));
    }

    let mut doy = day;
    for &dim in &days_in_month[..(month as usize - 1)] {
        doy += dim;
    }

    Ok(doy)
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── Julian Date ──────────────────────────────────────────────────

    #[test]
    fn j2000_epoch() {
        // J2000.0 = 2000 Jan 1, 12:00 TT → JD 2451545.0
        let jd = calendar_to_jd(2000, 1, 1.5).unwrap();
        assert!(
            (jd - J2000_JD).abs() < 1e-6,
            "J2000: {} vs {}",
            jd,
            J2000_JD
        );
    }

    #[test]
    fn known_date_sputnik() {
        // Sputnik launch: 1957 Oct 4 → JD ~2436116.0
        let jd = calendar_to_jd(1957, 10, 4.0).unwrap();
        assert!((jd - 2_436_115.5).abs() < 0.5, "Sputnik JD: {}", jd);
    }

    #[test]
    fn calendar_roundtrip() {
        let jd = calendar_to_jd(2024, 7, 15.75).unwrap();
        let (y, m, d) = jd_to_calendar(jd);
        assert_eq!(y, 2024);
        assert_eq!(m, 7);
        assert!((d - 15.75).abs() < 1e-10, "day: {d}");
    }

    #[test]
    fn calendar_roundtrip_leap() {
        let jd = calendar_to_jd(2024, 2, 29.0).unwrap();
        let (y, m, d) = jd_to_calendar(jd);
        assert_eq!(y, 2024);
        assert_eq!(m, 2);
        assert!((d - 29.0).abs() < 1e-10);
    }

    #[test]
    fn calendar_invalid_month() {
        assert!(calendar_to_jd(2024, 0, 1.0).is_err());
        assert!(calendar_to_jd(2024, 13, 1.0).is_err());
    }

    // ── MJD ──────────────────────────────────────────────────────────

    #[test]
    fn mjd_roundtrip() {
        let jd = 2_460_000.5;
        let mjd = jd_to_mjd(jd);
        let jd2 = mjd_to_jd(mjd);
        assert!((jd - jd2).abs() < 1e-15);
    }

    // ── Unix ─────────────────────────────────────────────────────────

    #[test]
    fn unix_epoch() {
        let jd = unix_to_jd(0.0);
        assert!((jd - UNIX_EPOCH_JD).abs() < 1e-10, "unix epoch JD: {}", jd);
    }

    #[test]
    fn unix_roundtrip() {
        let ts = 1_700_000_000.0; // ~2023-11-14
        let jd = unix_to_jd(ts);
        let ts2 = jd_to_unix(jd);
        // Large values lose some f64 precision
        assert!((ts - ts2).abs() < 0.01, "unix roundtrip: {ts} vs {ts2}");
    }

    // ── GMST ─────────────────────────────────────────────────────────

    #[test]
    fn gmst_range() {
        // GMST should always be in [0, 2π)
        for day_offset in 0..365 {
            let jd = J2000_JD + day_offset as f64;
            let g = gmst(jd);
            assert!(
                (0.0..std::f64::consts::TAU).contains(&g),
                "GMST out of range at JD {jd}: {g}"
            );
        }
    }

    #[test]
    fn gmst_monotonic_over_day() {
        // GMST should increase over the course of a day (modulo 2π wrap)
        let jd_base = J2000_JD + 100.0; // arbitrary day
        let g1 = gmst(jd_base);
        let g2 = gmst(jd_base + 0.25); // 6 hours later
        // 6 hours ≈ π/2 radians of sidereal rotation
        let diff = (g2 - g1).rem_euclid(std::f64::consts::TAU);
        assert!(
            (diff - std::f64::consts::FRAC_PI_2).abs() < 0.05,
            "6h should ≈ π/2 rad: diff={diff}"
        );
    }

    #[test]
    fn gmst_j2000() {
        // At J2000.0 (2000 Jan 1.5 UT1), GMST = 18h 41m 50.55s = 280.4606°
        // IAU 1982: at 0h UT1 on Jan 1 2000, GMST = 24110.54841s = 6h 41m 50.5s
        // Plus 12h of sidereal rotation for the half-day → 280.46°
        let g = gmst(J2000_JD);
        let g_deg = g.to_degrees();
        assert!((g_deg - 280.46).abs() < 0.1, "GMST at J2000: {g_deg}°");
    }

    // ── Day of year ──────────────────────────────────────────────────

    #[test]
    fn doy_jan1() {
        assert_eq!(day_of_year(2024, 1, 1).unwrap(), 1);
    }

    #[test]
    fn doy_leap_dec31() {
        assert_eq!(day_of_year(2024, 12, 31).unwrap(), 366);
    }

    #[test]
    fn doy_non_leap_dec31() {
        assert_eq!(day_of_year(2023, 12, 31).unwrap(), 365);
    }

    #[test]
    fn doy_march1_leap() {
        // 2024 is leap: Jan(31) + Feb(29) + 1 = 61
        assert_eq!(day_of_year(2024, 3, 1).unwrap(), 61);
    }

    #[test]
    fn doy_invalid() {
        assert!(day_of_year(2024, 0, 1).is_err());
        assert!(day_of_year(2024, 1, 0).is_err());
    }

    #[test]
    fn doy_feb30_rejected() {
        assert!(day_of_year(2024, 2, 30).is_err());
        // But Feb 29 is valid in leap year
        assert!(day_of_year(2024, 2, 29).is_ok());
        // And rejected in non-leap year
        assert!(day_of_year(2023, 2, 29).is_err());
    }

    #[test]
    fn doy_apr31_rejected() {
        assert!(day_of_year(2024, 4, 31).is_err());
        assert!(day_of_year(2024, 4, 30).is_ok());
    }

    // ── Julian centuries ─────────────────────────────────────────────

    #[test]
    fn centuries_at_j2000() {
        assert!((julian_centuries_since_j2000(J2000_JD)).abs() < 1e-15);
    }

    #[test]
    fn centuries_one_century() {
        let jd = J2000_JD + DAYS_PER_JULIAN_CENTURY;
        assert!((julian_centuries_since_j2000(jd) - 1.0).abs() < 1e-12);
    }
}
