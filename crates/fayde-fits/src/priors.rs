//! Priors from an image header.
//!
//! The rules (design doc §1): scale from instrument metadata, then a unit-checked linear WCS,
//! rejecting headers that carry both CD and PC matrices; pointing from the header WCS or pointing
//! keywords, when the header does not declare a frame other than ICRS/J2000; time from JD/MJD
//! keywords before ISO dates, when it falls in a plausible range. A header WCS is only ever a
//! *prior*: it is never treated as a solution.

use fayde_core::Portable;
use fayde_core::{Epoch, SkyPos};
use fayde_wcs::{Header, Value};

/// What a header says about a frame. Each field records where it came from.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct HeaderPriors {
    /// Pointing and its origin keyword(s).
    pub pointing: Option<(SkyPos, &'static str)>,
    /// Pixel scale (arcsec/px) and its origin.
    pub scale_arcsec: Option<(f64, &'static str)>,
    /// Mid-exposure epoch and its origin.
    pub epoch: Option<(Epoch, &'static str)>,
    /// `NAXIS1`, `NAXIS2` (or `IMAGEW`, `IMAGEH`).
    pub size: Option<(u32, u32)>,
    /// Problems found (ambiguous or invalid keywords that were ignored).
    pub warnings: Vec<String>,
}

fn num(h: &Header, k: &str) -> Option<f64> {
    match h.get(k)? {
        Value::Str(s) => s.trim().parse().ok(),
        v => v.as_f64(),
    }
    .filter(|v: &f64| v.is_finite())
}

fn text(h: &Header, k: &str) -> Option<String> {
    h.get(k)
        .and_then(|v| v.as_str())
        .map(|s| s.trim().to_string())
}

/// Parse `"hh mm ss.s"`, `"hh:mm:ss.s"` or `"-dd mm ss"` into a decimal value of the leading unit.
pub fn parse_sexagesimal(s: &str) -> Option<f64> {
    let s = s.trim();
    let neg = s.starts_with('-');
    let parts: Vec<f64> = s
        .trim_start_matches(['+', '-'])
        .split([' ', ':'])
        .filter(|p| !p.is_empty())
        .map(str::parse)
        .collect::<Result<_, _>>()
        .ok()?;
    if parts.is_empty() || parts.len() > 3 {
        return None;
    }
    // Minutes and seconds are unsigned and below 60; the sign belongs to the whole value.
    if parts[1..].iter().any(|v| !(0.0..60.0).contains(v)) || parts[0] < 0.0 {
        return None;
    }
    let v = parts[0] + parts.get(1).unwrap_or(&0.0) / 60.0 + parts.get(2).unwrap_or(&0.0) / 3600.0;
    Some(if neg { -v } else { v })
}

/// Modified Julian Date from an ISO-8601 date or date-time (`YYYY-MM-DD[Thh:mm:ss[.sss]]`), UTC.
pub fn mjd_from_iso(s: &str) -> Option<f64> {
    let s = s.trim();
    let (date, time) = match s.split_once('T') {
        Some((d, t)) => (d, Some(t)),
        None => (s, None),
    };
    let d: Vec<i64> = date
        .split('-')
        .map(str::parse)
        .collect::<Result<_, _>>()
        .ok()?;
    let [y, m, day] = d[..] else { return None };
    // Bounded before the day arithmetic below, which overflows for an absurd year.
    if !(1..=9999).contains(&y) {
        return None;
    }
    let leap = (y % 4 == 0 && y % 100 != 0) || y % 400 == 0;
    let days_in_month = [
        31,
        if leap { 29 } else { 28 },
        31,
        30,
        31,
        30,
        31,
        31,
        30,
        31,
        30,
        31,
    ];
    if !(1..=12).contains(&m) || day < 1 || day > days_in_month[(m - 1) as usize] {
        return None;
    }
    // Days from 1970-01-01 (H. Hinnant, "chrono-compatible low-level date algorithms").
    let y2 = if m <= 2 { y - 1 } else { y };
    let era = y2.div_euclid(400);
    let yoe = y2 - era * 400;
    let mp = (m + 9) % 12;
    let doy = (153 * mp + 2) / 5 + day - 1;
    let doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    let days = era * 146_097 + doe - 719_468;
    let mut frac = 0.0;
    if let Some(t) = time {
        let p: Vec<f64> = t
            .trim_end_matches('Z')
            .split(':')
            .map(str::parse)
            .collect::<Result<_, _>>()
            .ok()?;
        let (hh, mm, ss) = (
            *p.first()?,
            *p.get(1).unwrap_or(&0.0),
            *p.get(2).unwrap_or(&0.0),
        );
        // 60 ≤ ss < 61 allows a leap second.
        if p.len() > 3
            || !(0.0..24.0).contains(&hh)
            || !(0.0..60.0).contains(&mm)
            || !(0.0..61.0).contains(&ss)
        {
            return None;
        }
        frac = (hh + mm / 60.0 + ss / 3600.0) / 24.0;
    }
    Some(days as f64 + 40_587.0 + frac)
}

/// Why the header's coordinates are not ICRS/J2000, when it says so: a `RADESYS` other than ICRS
/// or FK5, or an `EQUINOX` other than 2000. The older `EPOCH` keyword is not consulted: it is
/// also written for the observation date, and a prior should not be lost to that ambiguity.
fn foreign_frame(h: &Header) -> Option<String> {
    let sys = text(h, "RADESYS")
        .or_else(|| text(h, "RADECSYS"))
        .map(|s| s.to_ascii_uppercase());
    if let Some(s) = sys.as_deref().filter(|s| !matches!(*s, "ICRS" | "FK5")) {
        return Some(format!("RADESYS = '{s}'"));
    }
    let equinox = match h.get("EQUINOX")? {
        Value::Str(s) => s.trim().trim_start_matches(['J', 'j']).parse::<f64>().ok(),
        v => v.as_f64(),
    };
    match equinox {
        Some(e) if (e - 2000.0).abs() < 1e-6 => None,
        Some(e) => Some(format!("EQUINOX = {e}")),
        None => Some("EQUINOX is not a number".into()),
    }
}

/// Read priors from an image header.
pub fn header_priors(h: &Header) -> HeaderPriors {
    let mut out = HeaderPriors::default();

    // Size.
    let dim = |a: &str, b: &str| Some((num(h, a)? as u32, num(h, b)? as u32));
    out.size = dim("NAXIS1", "NAXIS2")
        .or_else(|| dim("IMAGEW", "IMAGEH"))
        .filter(|(w, h)| *w > 0 && *h > 0);

    // Scale: instrument metadata first, then the linear WCS, then explicit scale keywords.
    let has_cd = ["CD1_1", "CD1_2", "CD2_1", "CD2_2"]
        .iter()
        .any(|k| h.contains(k));
    let has_pc = ["PC1_1", "PC1_2", "PC2_1", "PC2_2"]
        .iter()
        .any(|k| h.contains(k));
    let units_ok = ["CUNIT1", "CUNIT2"]
        .iter()
        .all(|k| text(h, k).is_none_or(|u| fayde_wcs::header::is_degree_unit(&u)));
    if has_cd && has_pc {
        out.warnings
            .push("both CDi_j and PCi_j present: linear WCS ignored for the scale".into());
    }
    if !units_ok {
        out.warnings
            .push("CUNIT is not deg: linear WCS ignored for the scale".into());
    }
    let from_optics = match (num(h, "XPIXSZ"), num(h, "FOCALLEN")) {
        (Some(px_um), Some(f_mm)) if px_um > 0.0 && f_mm > 0.0 => Some(206.264_806 * px_um / f_mm),
        _ => None,
    };
    let linear = if has_cd && has_pc || !units_ok {
        None
    } else if has_cd {
        let g = |k| num(h, k).unwrap_or(0.0);
        let det = g("CD1_1") * g("CD2_2") - g("CD1_2") * g("CD2_1");
        (det != 0.0).then(|| det.abs().sqrt() * 3600.0)
    } else {
        match (num(h, "CDELT1"), num(h, "CDELT2")) {
            (Some(c1), Some(c2)) if c1 != 0.0 && c2 != 0.0 => {
                let g = |k, d| num(h, k).unwrap_or(d);
                let det = g("PC1_1", 1.0) * g("PC2_2", 1.0) - g("PC1_2", 0.0) * g("PC2_1", 0.0);
                (det != 0.0).then(|| (c1 * c2 * det).abs().sqrt() * 3600.0)
            }
            _ => None,
        }
    };
    // Placeholder values (e.g. CDELT = 1 on frames without a real WCS) imply absurd fields.
    let plausible = |s: f64| match out.size {
        Some((w, hgt)) => 0.5 * f64::from(w).phypot(f64::from(hgt)) * s / 3600.0 < 60.0,
        None => s < 600.0,
    };
    let linear = match linear {
        Some(s) if !plausible(s) => {
            out.warnings.push(format!(
                "header linear WCS scale {s:.1} arcsec/px is implausible: ignored"
            ));
            None
        }
        other => other,
    };
    let explicit = num(h, "PIXSCALE")
        .or_else(|| num(h, "SECPIX"))
        .filter(|v| *v > 0.0);
    out.scale_arcsec = from_optics
        .map(|s| (s, "XPIXSZ/FOCALLEN"))
        .or(linear.map(|s| (s, "header linear WCS")))
        .or(explicit.map(|s| (s, "PIXSCALE/SECPIX")));

    // Pointing: header WCS reference point, else pointing keywords.
    let tan = text(h, "CTYPE1").is_some_and(|c| c.starts_with("RA---"))
        && text(h, "CTYPE2").is_some_and(|c| c.starts_with("DEC--"));
    // CRVAL is in CUNIT (default degrees): radians are converted; any other unit leaves the
    // WCS reference point unused, with a warning, rather than misread as degrees.
    let crval_factor = {
        let units: Vec<Option<String>> = ["CUNIT1", "CUNIT2"].iter().map(|k| text(h, k)).collect();
        let rad = |u: &Option<String>| {
            u.as_deref()
                .is_some_and(|u| u.trim().eq_ignore_ascii_case("rad"))
        };
        if units
            .iter()
            .all(|u| u.as_deref().is_none_or(fayde_wcs::header::is_degree_unit))
        {
            Some(1.0)
        } else if units.iter().all(rad) {
            Some(180.0 / std::f64::consts::PI)
        } else {
            if tan {
                out.warnings
                    .push("CUNIT is neither deg nor rad: CRVAL not used for the pointing".into());
            }
            None
        }
    };
    let wcs_ptg = (tan && crval_factor.is_some())
        .then(|| {
            let f = crval_factor?;
            // A longitude beyond one turn either way is not a coordinate; it is not wrapped.
            let ra = Some(num(h, "CRVAL1")? * f).filter(|ra| ra.abs() <= 360.0)?;
            SkyPos::from_degrees(ra, num(h, "CRVAL2")? * f)
        })
        .flatten();
    let kw = |ra_k: &str, dec_k: &str| -> Option<SkyPos> {
        let ra: f64 = match h.get(ra_k)? {
            Value::Str(s) if s.contains([' ', ':']) => parse_sexagesimal(s)? * 15.0,
            Value::Str(s) => s.trim().parse().ok()?,
            v => v.as_f64()?,
        };
        // Out of range (25 hours, 400°) is an error in the header, not a value to wrap.
        if !(0.0..360.0).contains(&ra) {
            return None;
        }
        let dec = match h.get(dec_k)? {
            Value::Str(s) if s.contains([' ', ':']) => parse_sexagesimal(s)?,
            Value::Str(s) => s.trim().parse().ok()?,
            v => v.as_f64()?,
        };
        SkyPos::from_degrees(ra, dec)
    };
    let found = wcs_ptg
        .map(|p| (p, "CRVAL1/CRVAL2"))
        .or_else(|| kw("RA", "DEC").map(|p| (p, "RA/DEC")))
        .or_else(|| kw("OBJCTRA", "OBJCTDEC").map(|p| (p, "OBJCTRA/OBJCTDEC")));
    for (a, b) in [("RA", "DEC"), ("OBJCTRA", "OBJCTDEC")] {
        if h.contains(a) && h.contains(b) && kw(a, b).is_none() {
            out.warnings
                .push(format!("{a}/{b} is not a valid position: ignored"));
        }
    }
    // The pointing is used as ICRS. Coordinates the header declares in another frame or equinox
    // (FK4/B1950, the equinox of date) are up to about 0.7° away from that, so they are not used.
    out.pointing = match (found, foreign_frame(h)) {
        (Some((_, origin)), Some(why)) => {
            out.warnings.push(format!(
                "{origin} not used for the pointing: the header gives {why}, not ICRS/J2000"
            ));
            None
        }
        (found, _) => found,
    };

    // Time: numeric clocks before ISO strings; DATE-OBS is the exposure start, so add half the
    // exposure time for mid-exposure.
    let half_exp_days = num(h, "EXPTIME")
        .or_else(|| num(h, "EXPOSURE"))
        .filter(|v| *v > 0.0)
        .map_or(0.0, |e| e / 2.0 / 86_400.0);
    let date_obs = || {
        let d = text(h, "DATE-OBS")?;
        let iso = match (d.contains('T'), text(h, "TIME-OBS")) {
            (false, Some(t)) => format!("{d}T{t}"),
            _ => d,
        };
        Epoch::from_mjd(mjd_from_iso(&iso)? + half_exp_days)
    };
    let clocks: [(&str, Option<Epoch>, &'static str); 3] = [
        (
            "JD-OBS",
            num(h, "JD-OBS").and_then(|jd| Epoch::from_jd(jd + half_exp_days)),
            "JD-OBS + exposure/2",
        ),
        (
            "MJD-OBS",
            num(h, "MJD-OBS").and_then(|m| Epoch::from_mjd(m + half_exp_days)),
            "MJD-OBS + exposure/2",
        ),
        ("DATE-OBS", date_obs(), "DATE-OBS + exposure/2"),
    ];
    if clocks[2].1.is_none() && text(h, "DATE-OBS").is_some() {
        out.warnings
            .push("DATE-OBS is not a valid date: ignored".into());
    }
    for (keyword, epoch, origin) in clocks {
        let Some(epoch) = epoch else { continue };
        if epoch.is_plausible() {
            out.epoch = Some((epoch, origin));
            break;
        }
        // A clock error (a Julian Date in MJD-OBS, an unset camera clock): say so and try the
        // next keyword rather than move every star by centuries of proper motion.
        out.warnings.push(format!(
            "{keyword} gives year {:.1}, outside {:.0} to {:.0}: ignored",
            epoch.year(),
            Epoch::PLAUSIBLE_YEARS.start(),
            Epoch::PLAUSIBLE_YEARS.end()
        ));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn crval_follows_cunit() {
        let tan = ["CTYPE1  = 'RA---TAN'", "CTYPE2  = 'DEC--TAN'"];
        let rad = Header::parse_cards(tan.into_iter().chain([
            "CUNIT1  = 'rad     '",
            "CUNIT2  = 'rad     '",
            "CRVAL1  =                  1.0",
            "CRVAL2  =                  0.5",
        ]));
        let (p, _) = header_priors(&rad).pointing.unwrap();
        let (ra, dec) = p.to_degrees();
        assert!((ra - 57.295_779_5).abs() < 1e-6 && (dec - 28.647_889_8).abs() < 1e-6);
        let arcmin = Header::parse_cards(tan.into_iter().chain([
            "CUNIT1  = 'arcmin  '",
            "CUNIT2  = 'arcmin  '",
            "CRVAL1  =                  1.0",
            "CRVAL2  =                  0.5",
        ]));
        // Not misread as degrees; no RA/DEC keywords either, so no pointing at all.
        assert!(header_priors(&arcmin).pointing.is_none());
    }

    #[test]
    fn sexagesimal_and_dates() {
        assert!((parse_sexagesimal("05 22 57.30").unwrap() - 5.382_583_333).abs() < 1e-9);
        assert!((parse_sexagesimal("-00:11:24.3").unwrap() + 0.190_083_333).abs() < 1e-9);
        assert!(parse_sexagesimal("a b").is_none());
        assert_eq!(mjd_from_iso("1858-11-17"), Some(0.0));
        assert_eq!(mjd_from_iso("2000-01-01T12:00:00"), Some(51_544.5));
        assert!((mjd_from_iso("2020-11-21T16:46:40.074").unwrap() - 59_174.699_075).abs() < 1e-6);
        assert!(mjd_from_iso("2020-13-01").is_none());
        assert!(mjd_from_iso("2021-02-29").is_none());
        assert!(mjd_from_iso("2020-02-29").is_some());
        assert!(mjd_from_iso("2020-04-31").is_none());
        assert!(mjd_from_iso("2020-01-01T25:00:00").is_none());
        assert!(mjd_from_iso("2020-01-01T12:61:00").is_none());
        assert!(mjd_from_iso("2016-12-31T23:59:60.5").is_some());
        assert!((parse_sexagesimal("- 00 11 24.3").unwrap() + 0.190_083_333).abs() < 1e-9);
        assert!(parse_sexagesimal("-00 -11 -24.3").is_none());
        assert!(parse_sexagesimal("05 61 00").is_none());
    }
}
