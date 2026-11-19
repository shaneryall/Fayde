use fayde_core::Portable;
use fayde_core::{Epoch, SkyPos};
use fayde_wcs::{Header, Value};

#[derive(Debug, Clone, PartialEq, Default)]
pub struct HeaderPriors {
    pub pointing: Option<(SkyPos, &'static str)>,
    pub scale_arcsec: Option<(f64, &'static str)>,
    pub epoch: Option<(Epoch, &'static str)>,
    pub size: Option<(u32, u32)>,
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

pub fn parse_sexagesimal(s: &str) -> Option<f64> {
    let s = s.trim();
    let neg = s.starts_with('-');
    let unsigned = s.strip_prefix(['+', '-']).unwrap_or(s);
    if unsigned.starts_with(['+', '-']) {
        return None;
    }
    let parts: Vec<f64> = unsigned
        .split([' ', ':'])
        .filter(|p| !p.is_empty())
        .map(str::parse)
        .collect::<Result<_, _>>()
        .ok()?;
    if parts.is_empty() || parts.len() > 3 {
        return None;
    }
    let last = parts.len() - 1;
    if parts.iter().enumerate().skip(1).any(|(i, v)| {
        if i == last {
            !(0.0..=60.0).contains(v)
        } else {
            !(0.0..60.0).contains(v)
        }
    }) || !parts[0].is_finite()
        || parts[0] < 0.0
    {
        return None;
    }
    let v = parts[0] + parts.get(1).unwrap_or(&0.0) / 60.0 + parts.get(2).unwrap_or(&0.0) / 3600.0;
    Some(if neg { -v } else { v })
}

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

pub fn header_priors(h: &Header) -> HeaderPriors {
    header_priors_for_size(h, (None, None))
}

pub fn header_priors_for_size(h: &Header, size: (Option<u32>, Option<u32>)) -> HeaderPriors {
    let mut out = HeaderPriors::default();

    let dim = |key| num(h, key).and_then(crate::pixel_dimension);
    out.size = dim("NAXIS1")
        .zip(dim("NAXIS2"))
        .or_else(|| dim("IMAGEW").zip(dim("IMAGEH")));

    let size = size
        .0
        .or(out.size.map(|s| s.0))
        .zip(size.1.or(out.size.map(|s| s.1)));

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
    let plausible = |s: f64| match size {
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

    let tan = text(h, "CTYPE1").is_some_and(|c| c.starts_with("RA---"))
        && text(h, "CTYPE2").is_some_and(|c| c.starts_with("DEC--"));
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
            let ra = Some(num(h, "CRVAL1")? * f).filter(|ra| ra.abs() <= 360.0)?;
            SkyPos::from_degrees(ra, num(h, "CRVAL2")? * f)
        })
        .flatten();
    let kw = |ra_k: &str, dec_k: &str| -> Option<SkyPos> {
        let ra: f64 = match h.get(ra_k)? {
            Value::Str(s) if s.trim().contains([' ', ':']) => {
                let ra = parse_sexagesimal(s)? * 15.0;
                if ra == 360.0 {
                    0.0
                } else {
                    ra
                }
            }
            Value::Str(s) => s.trim().parse().ok()?,
            v => v.as_f64()?,
        };
        if !(0.0..360.0).contains(&ra) {
            return None;
        }
        let dec = match h.get(dec_k)? {
            Value::Str(s) if s.trim().contains([' ', ':']) => parse_sexagesimal(s)?,
            Value::Str(s) => s.trim().parse().ok()?,
            v => v.as_f64()?,
        };
        SkyPos::from_degrees(ra, dec)
    };
    let centre = size.and_then(|(w, hgt)| {
        linear?;
        wcs_ptg?;
        let wcs = fayde_wcs::Wcs::from_header(h).ok()?;
        let pixel =
            fayde_core::PixelPos::new((f64::from(w) - 1.0) / 2.0, (f64::from(hgt) - 1.0) / 2.0);
        if pixel == wcs.crpix() {
            return wcs_ptg.map(|p| (p, "CRVAL1/CRVAL2"));
        }
        Some((wcs.pixel_to_sky(pixel), "WCS image centre"))
    });
    let found = centre
        .or_else(|| wcs_ptg.map(|p| (p, "CRVAL1/CRVAL2")))
        .or_else(|| kw("RA", "DEC").map(|p| (p, "RA/DEC")))
        .or_else(|| kw("OBJCTRA", "OBJCTDEC").map(|p| (p, "OBJCTRA/OBJCTDEC")));
    for (a, b) in [("RA", "DEC"), ("OBJCTRA", "OBJCTDEC")] {
        if h.contains(a) && h.contains(b) && kw(a, b).is_none() {
            out.warnings
                .push(format!("{a}/{b} is not a valid position: ignored"));
        }
    }
    out.pointing = match (found, foreign_frame(h)) {
        (Some((_, origin)), Some(why)) => {
            out.warnings.push(format!(
                "{origin} not used for the pointing: the header gives {why}, not ICRS/J2000"
            ));
            None
        }
        (found, _) => found,
    };

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
    fn invalid_axis_sizes_fall_back_to_image_dimensions() {
        for bad in ["NAXIS1  = 0", "NAXIS1  = 99.5", "NAXIS1  = 1e12"] {
            let h = Header::parse_cards([bad, "NAXIS2  = 20", "IMAGEW  = 100", "IMAGEH  = 50"]);
            assert_eq!(header_priors(&h).size, Some((100, 50)));
        }
        for input in ["--12:00", "+-2:30", "NaN", "inf:00"] {
            assert!(parse_sexagesimal(input).is_none(), "{input}");
        }
    }

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

    #[test]
    fn sexagesimal_seconds_rounded_up_to_sixty() {
        assert!((parse_sexagesimal("20:29:60.0").unwrap() - 20.5).abs() < 1e-12);
        assert!((parse_sexagesimal("-09 43 60.0").unwrap() + (9.0 + 44.0 / 60.0)).abs() < 1e-12);
        assert!(parse_sexagesimal("20:29:60.5").is_none());
        assert!(parse_sexagesimal("20:60:00").is_none());
        assert!((parse_sexagesimal("20 60.0").unwrap() - 21.0).abs() < 1e-12);
    }

    #[test]
    fn header_pointing_with_seconds_rounded_up_to_sixty() {
        let h = Header::parse_cards(["RA      = '23:59:60.0'", "DEC     = '-06:56:60.0'"]);
        let (p, _) = header_priors(&h).pointing.expect("pointing kept");
        let (ra, dec) = p.to_degrees();
        assert!(ra.abs() < 1e-9 || (ra - 360.0).abs() < 1e-9, "{ra}");
        assert!((dec + (6.0 + 57.0 / 60.0)).abs() < 1e-9, "{dec}");
    }

    #[test]
    fn padded_decimal_coordinates_are_degrees_not_sexagesimal() {
        let h = Header::parse_cards(["RA      = ' 24.0'", "DEC     = '  2.5 '"]);
        let (p, _) = header_priors(&h).pointing.expect("pointing kept");
        let (ra, dec) = p.to_degrees();
        assert!((ra - 24.0).abs() < 1e-9, "{ra}");
        assert!((dec - 2.5).abs() < 1e-9, "{dec}");
        let h = Header::parse_cards(["RA      = ' 12.5'", "DEC     = '-1.0'"]);
        let (p, _) = header_priors(&h).pointing.expect("pointing kept");
        assert!((p.to_degrees().0 - 12.5).abs() < 1e-9);
    }
}
