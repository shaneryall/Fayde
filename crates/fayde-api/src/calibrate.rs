//! An instrument profile from solved frames (`fayde calibrate`).
//!
//! A profile holds what frames of one instrument share: the pixel scale, the parity, the rule
//! that gives the roll from a header keyword, and how far the header's pointing is from the
//! truth. All four can be read off a set of solutions, so a profile needs no other solver: solve
//! some frames with the pointing and a nominal scale alone, and calibrate from those.
//!
//! Every value is a median, and every spread comes from the median absolute deviation (Hampel
//! 1974; the factor 1.4826 makes it the standard deviation of a Gaussian), so a few wrong or
//! poor solutions among the frames do not move the result.

use fayde_core::{to_arcsec, PixelPos, Portable};
use fayde_fits::header_priors;
use fayde_wcs::{header::Value as Card, Header, Parity, Wcs};
use serde_json::{json, Value};

/// Fewest solved frames a calibration is made from.
pub const CALIBRATE_MIN_FRAMES: usize = 5;
/// The scale's 1σ is no tighter than this fraction of the scale. The solver's hard limit is
/// 3σ on the scale it fits to a frame, and a fit to a dozen stars scatters by 0.1 to 0.2%
/// whatever the instrument does: a σ from well-populated frames alone rejects the sparse ones.
const SCALE_SIGMA_FLOOR: f64 = 0.001;
/// The scale's 1σ as a multiple of the robust spread of the fitted scales.
const SCALE_SIGMA_FACTOR: f64 = 1.5;
/// The roll's 1σ as a multiple of the robust spread about the rule, and its floor, degrees.
const ROLL_SIGMA_FACTOR: f64 = 2.0;
const ROLL_SIGMA_FLOOR_DEG: f64 = 0.05;
/// The roll's hard limit as a multiple of its 1σ, and its floor, degrees.
const ROLL_LIMIT_FACTOR: f64 = 10.0;
const ROLL_LIMIT_FLOOR_DEG: f64 = 1.0;
/// A keyword whose rule leaves more scatter than this, degrees, does not give the roll.
const ROLL_MAX_SPREAD_DEG: f64 = 2.0;
/// Median of a Rayleigh distribution in units of its σ, `√(2 ln 2)`: the radial offset of a
/// pointing whose two axes each have a Gaussian error σ.
const RAYLEIGH_MEDIAN: f64 = 1.177_410_022_515_474_7;

/// One solved frame.
#[derive(Debug, Clone)]
pub struct SolvedFrame {
    /// The frame's header.
    pub header: Header,
    /// Its solution.
    pub wcs: Wcs,
    /// Image width, pixels.
    pub width: u32,
    /// Image height, pixels.
    pub height: u32,
}

/// What [`calibrate_profile`] is told about the instrument's headers.
#[derive(Debug, Clone, Default)]
pub struct CalibrateOptions {
    /// Header keyword holding the instrument's own roll angle, degrees. Without it the profile
    /// has no roll rule.
    pub roll_keyword: Option<String>,
    /// Header keywords of the binning factors along x and y (the x factor when y is absent).
    /// With them the profile's scale is the unbinned one.
    pub binning: Option<(String, Option<String>)>,
    /// The profile's name.
    pub name: Option<String>,
    /// Measure the header pointing without the header's WCS, and have the profile tell the
    /// solver to read it the same way (`header_wcs: false`).
    pub no_header_wcs: bool,
}

/// Median of `v` (the upper of the two middle values of an even count). `v` is not empty.
fn median(v: &[f64]) -> f64 {
    let mut v = v.to_vec();
    v.sort_by(f64::total_cmp);
    v[v.len() / 2]
}

/// Median and robust standard deviation (1.4826 × the median absolute deviation).
fn robust(v: &[f64]) -> (f64, f64) {
    let m = median(v);
    let dev: Vec<f64> = v.iter().map(|x| (x - m).abs()).collect();
    (m, 1.4826 * median(&dev))
}

/// An angle in degrees wrapped into [−180, 180).
fn wrap(a: f64) -> f64 {
    (a + 180.0).rem_euclid(360.0) - 180.0
}

fn number(h: &Header, key: &str) -> Option<f64> {
    match h.get(key)? {
        Card::Str(s) => s.trim().parse().ok(),
        v => v.as_f64(),
    }
    .filter(|v: &f64| v.is_finite())
}

/// `v` to `decimals` decimal places: the number that prints with that many and no more.
fn rounded(v: f64, decimals: usize) -> f64 {
    format!("{v:.decimals$}").parse().unwrap_or(v)
}

/// An instrument profile from solved frames, as the JSON `--instrument` reads.
///
/// - `scale_arcsec`: the median fitted scale (unbinned when binning keywords are given);
///   `scale_sigma_arcsec`: 1.5 × its robust spread, at least 0.1% of the scale.
/// - `parity`: the solutions' parity, when they all agree.
/// - `roll`: with a roll keyword, the rule `roll = offset + sign × keyword` that fits best
///   (both signs are tried; the offset is the median), its 1σ twice the robust spread about
///   the rule (at least 0.05°) and its hard limit ten times that (at least 1°).
/// - `pointing_sigma_deg`: the per-axis error of the header's pointing, from the median
///   distance between it and the solved frame centre.
///
/// The object also carries a `calibration` member with the numbers behind it (frames used,
/// spreads, worst deviations), which the profile reader ignores.
pub fn calibrate_profile(frames: &[SolvedFrame], o: &CalibrateOptions) -> Result<Value, String> {
    if frames.len() < CALIBRATE_MIN_FRAMES {
        return Err(format!(
            "calibration needs at least {CALIBRATE_MIN_FRAMES} solved frames, got {}",
            frames.len()
        ));
    }
    let frame = |i: usize| format!("solved frame {}", i + 1);

    // Scale, unbinned.
    let mut scales = Vec::with_capacity(frames.len());
    for (i, f) in frames.iter().enumerate() {
        let bin = match &o.binning {
            None => 1.0,
            Some((x_key, y_key)) => {
                let bx = number(&f.header, x_key);
                let by = y_key.as_deref().map_or(bx, |k| number(&f.header, k));
                match (bx, by) {
                    (Some(x), Some(y)) if x >= 1.0 && y >= 1.0 => (x * y).sqrt(),
                    _ => {
                        return Err(format!(
                            "{}: no valid binning in {x_key}{}",
                            frame(i),
                            y_key.as_deref().map_or(String::new(), |k| format!("/{k}"))
                        ))
                    }
                }
            }
        };
        let scale = f.wcs.pixel_scale_arcsec() / bin;
        if !(scale.is_finite() && scale > 0.0) {
            return Err(format!("{}: its WCS has no scale", frame(i)));
        }
        scales.push(scale);
    }
    let (scale, scale_spread) = robust(&scales);
    let scale_sigma = (SCALE_SIGMA_FACTOR * scale_spread).max(SCALE_SIGMA_FLOOR * scale);

    // Parity: one camera has one.
    let positive = frames
        .iter()
        .filter(|f| f.wcs.parity() == Parity::Positive)
        .count();
    let parity = match positive {
        0 => Some("negative"),
        n if n == frames.len() => Some("positive"),
        _ => None,
    };

    // Roll: the position angle of the image +y axis, east of north, against the keyword.
    let mut roll = None;
    let mut roll_stats = Value::Null;
    if let Some(key) = &o.roll_keyword {
        let mut pairs = Vec::with_capacity(frames.len());
        for (i, f) in frames.iter().enumerate() {
            let v = number(&f.header, key)
                .ok_or_else(|| format!("{}: no number in {key}", frame(i)))?;
            let cd = f.wcs.cd();
            pairs.push((cd[0][1].patan2(cd[1][1]).to_degrees(), v));
        }
        // offset = roll − sign × keyword, for the sign whose offsets agree best. The offsets
        // are centred on their mean direction first, so that values near ±180° do not split.
        let fit = |sign: f64| {
            let offsets: Vec<f64> = pairs.iter().map(|(r, v)| wrap(r - sign * v)).collect();
            let (sn, cs) = offsets.iter().fold((0.0, 0.0), |(s, c), a| {
                let (sin, cos) = a.to_radians().psin_cos();
                (s + sin, c + cos)
            });
            let centre = sn.patan2(cs).to_degrees();
            let about: Vec<f64> = offsets.iter().map(|a| wrap(a - centre)).collect();
            let (m, spread) = robust(&about);
            let offset = wrap(centre + m);
            let worst = offsets
                .iter()
                .map(|a| wrap(a - offset).abs())
                .fold(0.0, f64::max);
            (sign, offset, spread, worst)
        };
        let (plus, minus) = (fit(1.0), fit(-1.0));
        let (sign, offset, spread, worst) = if minus.2 < plus.2 { minus } else { plus };
        if spread > ROLL_MAX_SPREAD_DEG {
            return Err(format!(
                "{key} does not give the roll: the best rule leaves a scatter of {spread:.1}°"
            ));
        }
        let sigma = (ROLL_SIGMA_FACTOR * spread).max(ROLL_SIGMA_FLOOR_DEG);
        roll = Some(json!({
            "keyword": key,
            "offset_deg": rounded(offset, 3),
            "sign": sign,
            "sigma_deg": rounded(sigma, 3),
            "max_dev_deg": rounded((ROLL_LIMIT_FACTOR * sigma).max(ROLL_LIMIT_FLOOR_DEG), 3),
        }));
        roll_stats = json!({"robust_sigma_deg": spread, "worst_deg": worst});
    }

    // Pointing: the header's against the solved frame centre.
    let mut offsets = Vec::new();
    for f in frames {
        let header = if o.no_header_wcs {
            header_priors(&f.header.without_wcs())
        } else {
            header_priors(&f.header)
        };
        if let Some((pointing, _)) = header.pointing {
            let centre = f.wcs.pixel_to_sky(PixelPos::new(
                (f64::from(f.width) - 1.0) / 2.0,
                (f64::from(f.height) - 1.0) / 2.0,
            ));
            offsets.push(to_arcsec(pointing.separation(centre)) / 3600.0);
        }
    }
    let mut pointing_stats = Value::Null;
    let pointing_sigma = (offsets.len() >= CALIBRATE_MIN_FRAMES).then(|| {
        let mut sorted = offsets.clone();
        sorted.sort_by(f64::total_cmp);
        let at = |q: f64| sorted[((sorted.len() as f64 - 1.0) * q).round() as usize];
        pointing_stats = json!({
            "frames": sorted.len(),
            "offset_median_deg": at(0.5),
            "offset_p90_deg": at(0.9),
            "offset_max_deg": at(1.0),
        });
        // Not below 0.36″: a pointing that good is beyond what the offsets can show.
        rounded((median(&offsets) / RAYLEIGH_MEDIAN).max(1e-4), 5)
    });

    let mut profile = serde_json::Map::new();
    let mut put = |k: &str, v: Value| {
        profile.insert(k.into(), v);
    };
    put(
        "name",
        json!(o
            .name
            .clone()
            .unwrap_or_else(|| "calibrated profile".into())),
    );
    put(
        "source",
        json!(format!(
            "Calibrated by fayde {} from {} solved frames: medians, with spreads from the median \
             absolute deviation.",
            env!("CARGO_PKG_VERSION"),
            frames.len()
        )),
    );
    put("scale_arcsec", json!(rounded(scale, 5)));
    put("scale_sigma_arcsec", json!(rounded(scale_sigma, 5)));
    if let Some(p) = parity {
        put("parity", json!(p));
    }
    if let Some(r) = roll {
        put("roll", r);
    }
    if let Some((x, y)) = &o.binning {
        put("binning", json!({"x_keyword": x, "y_keyword": y}));
    }
    if let Some(s) = pointing_sigma {
        put("pointing_sigma_deg", json!(s));
    }
    if o.no_header_wcs {
        put("header_wcs", json!(false));
    }
    put(
        "calibration",
        json!({
            "frames": frames.len(),
            "scale_robust_sigma_arcsec": scale_spread,
            "scale_min_arcsec": scales.iter().copied().fold(f64::INFINITY, f64::min),
            "scale_max_arcsec": scales.iter().copied().fold(f64::NEG_INFINITY, f64::max),
            "parity_positive_frames": positive,
            "roll": roll_stats,
            "pointing": pointing_stats,
        }),
    );
    Ok(Value::Object(profile))
}

#[cfg(test)]
mod tests {
    use super::*;
    use fayde_core::SkyPos;
    use fayde_fits::InstrumentProfile;
    use fayde_wcs::TanWcs;

    /// A frame of a camera with scale 2.5″/px (unbinned), positive parity and
    /// `roll = 33 − ROLLKEY`, whose header pointing is `off` degrees south of the truth.
    fn frame(k: usize, bin: f64, scale_err: f64, roll_err: f64, off: f64) -> SolvedFrame {
        let key = -170.0 + 37.0 * k as f64; // crosses ±180° in the offsets
        let roll: f64 = (33.0 - key + roll_err).to_radians();
        let s = 2.5 * bin * (1.0 + scale_err) / 3600.0;
        let (sn, cs) = roll.psin_cos();
        let (ra, dec) = (20.0 + 30.0 * k as f64, -40.0 + 9.0 * k as f64);
        let wcs = TanWcs::new(
            SkyPos::from_degrees(ra, dec).unwrap(),
            PixelPos::new(255.5, 255.5),
            [[s * cs, s * sn], [-s * sn, s * cs]],
        )
        .unwrap();
        assert_eq!(wcs.parity(), Parity::Positive);
        let mut header = Header::default();
        header.set("ROLLKEY", Card::Float(key));
        header.set("XBINNING", Card::Int(bin as i64));
        header.set("YBINNING", Card::Int(bin as i64));
        header.set("RA", Card::Float(ra));
        header.set("DEC", Card::Float(dec - off));
        // A WCS in the header that is far from the truth: only read when allowed.
        header.set("CTYPE1", Card::Str("RA---TAN".into()));
        header.set("CTYPE2", Card::Str("DEC--TAN".into()));
        header.set("CRVAL1", Card::Float(ra));
        header.set("CRVAL2", Card::Float(dec - 1.0));
        SolvedFrame {
            header,
            wcs: wcs.into(),
            width: 512,
            height: 512,
        }
    }

    fn options() -> CalibrateOptions {
        CalibrateOptions {
            roll_keyword: Some("ROLLKEY".into()),
            binning: Some(("XBINNING".into(), Some("YBINNING".into()))),
            name: Some("test camera".into()),
            no_header_wcs: true,
        }
    }

    #[test]
    fn recovers_scale_parity_roll_rule_and_pointing_error() {
        // Nine frames, one of them binned, one with a wild scale and roll (a poor solution).
        let frames: Vec<SolvedFrame> = (0..9)
            .map(|k| {
                let bin = if k == 3 { 2.0 } else { 1.0 };
                let (se, re) = if k == 5 {
                    (0.03, 4.0)
                } else {
                    (1e-4 * (k as f64 - 4.0), 0.01 * (k as f64 - 4.0))
                };
                frame(k, bin, se, re, 0.01 + 0.002 * k as f64)
            })
            .collect();
        let json = calibrate_profile(&frames, &options()).unwrap();
        // The profile reader accepts it.
        let p = InstrumentProfile::from_json(&json.to_string()).unwrap();
        assert_eq!(p.name.as_deref(), Some("test camera"));
        assert!((p.scale_arcsec.unwrap() - 2.5).abs() < 5e-4, "{json}");
        // The floor: 0.1% of the scale, though the frames agree far better.
        assert!(
            (p.scale_sigma_arcsec.unwrap() - 0.0025).abs() < 1e-9,
            "{json}"
        );
        assert_eq!(p.parity.as_deref(), Some("positive"));
        let roll = p.roll.unwrap();
        assert_eq!((roll.keyword.as_str(), roll.sign), ("ROLLKEY", -1.0));
        assert!((roll.offset_deg - 33.0).abs() < 0.03, "{json}");
        // Spread about the rule 0.03° × 1.4826, doubled; the limit is at its 1° floor.
        assert!(
            (roll.sigma_deg - 0.089).abs() < 0.005 && roll.max_dev_deg == 1.0,
            "{json}"
        );
        // Offsets of 0.010° to 0.026°: median 0.018°, so σ = 0.018 / 1.1774.
        assert!(
            (p.pointing_sigma_deg.unwrap() - 0.018 / 1.1774).abs() < 1e-4,
            "{json}"
        );
        assert_eq!(p.header_wcs, Some(false));
        assert_eq!(json["binning"]["x_keyword"], "XBINNING");
        assert_eq!(json["calibration"]["frames"], 9);
        assert!((json["calibration"]["roll"]["worst_deg"].as_f64().unwrap() - 4.0).abs() < 0.05);

        // With the header WCS allowed the pointing is the header's CRVAL, 1° off.
        let with_wcs = CalibrateOptions {
            no_header_wcs: false,
            ..options()
        };
        let json = calibrate_profile(&frames, &with_wcs).unwrap();
        assert!((json["pointing_sigma_deg"].as_f64().unwrap() - 1.0 / 1.1774).abs() < 1e-3);
        assert!(json.get("header_wcs").is_none());
    }

    #[test]
    fn says_what_it_cannot_calibrate() {
        let frames: Vec<SolvedFrame> = (0..6).map(|k| frame(k, 1.0, 0.0, 0.0, 0.01)).collect();
        // Too few frames.
        let err = calibrate_profile(&frames[..4], &options()).unwrap_err();
        assert!(err.contains("at least 5"), "{err}");
        // A keyword the roll does not follow.
        let mut shuffled = frames.clone();
        for (k, f) in shuffled.iter_mut().enumerate() {
            f.header.set("ROLLKEY", Card::Float(71.0 * (k * k) as f64));
        }
        let err = calibrate_profile(&shuffled, &options()).unwrap_err();
        assert!(err.contains("does not give the roll"), "{err}");
        // A frame without the keyword, or without its binning.
        let mut missing = frames.clone();
        missing[2].header = Header::default();
        let err = calibrate_profile(&missing, &options()).unwrap_err();
        assert!(err.contains("solved frame 3"), "{err}");
        // Mixed parity gives no parity; no roll keyword, no roll rule; no pointing in the
        // headers, no pointing error.
        let mut mixed = frames.clone();
        let cd = mixed[1].wcs.cd();
        mixed[1].wcs = TanWcs::new(
            mixed[1].wcs.crval(),
            mixed[1].wcs.crpix(),
            [[-cd[0][0], cd[0][1]], [-cd[1][0], cd[1][1]]],
        )
        .unwrap()
        .into();
        for f in &mut mixed {
            let mut h = Header::default();
            h.set("XBINNING", Card::Int(1));
            f.header = h;
        }
        let plain = CalibrateOptions {
            roll_keyword: None,
            binning: Some(("XBINNING".into(), None)),
            ..CalibrateOptions::default()
        };
        let json = calibrate_profile(&mixed, &plain).unwrap();
        assert!(json.get("parity").is_none() && json.get("roll").is_none());
        assert!(json.get("pointing_sigma_deg").is_none());
        assert_eq!(json["calibration"]["parity_positive_frames"], 5);
        InstrumentProfile::from_json(&json.to_string()).unwrap();
    }
}
