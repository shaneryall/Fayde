use fayde_core::{math::round_decimal, to_arcsec, PixelPos, Portable};
use fayde_fits::header_priors_for_size;
use fayde_wcs::{header::Value as Card, Header, Parity, Wcs};
use serde_json::{json, Value};

pub const CALIBRATE_MIN_FRAMES: usize = 5;
const SCALE_SIGMA_FLOOR: f64 = 0.001;
const SCALE_SIGMA_FACTOR: f64 = 1.5;
const ROLL_SIGMA_FACTOR: f64 = 2.0;
const ROLL_SIGMA_FLOOR_DEG: f64 = 0.05;
const ROLL_LIMIT_FACTOR: f64 = 10.0;
const ROLL_LIMIT_FLOOR_DEG: f64 = 1.0;
const ROLL_MAX_SPREAD_DEG: f64 = 2.0;
const RAYLEIGH_MEDIAN: f64 = 1.177_410_022_515_474_7;

#[derive(Debug, Clone)]
pub struct SolvedFrame {
    pub header: Header,
    pub wcs: Wcs,
    pub width: u32,
    pub height: u32,
}

#[derive(Debug, Clone, Default)]
pub struct CalibrateOptions {
    pub roll_keyword: Option<String>,
    pub binning: Option<(String, Option<String>)>,
    pub name: Option<String>,
    pub no_header_wcs: bool,
}

fn median(v: &[f64]) -> f64 {
    let mut v = v.to_vec();
    v.sort_by(f64::total_cmp);
    v[v.len() / 2]
}

fn robust(v: &[f64]) -> (f64, f64) {
    let m = median(v);
    let dev: Vec<f64> = v.iter().map(|x| (x - m).abs()).collect();
    (m, 1.4826 * median(&dev))
}

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

pub fn calibrate_profile(frames: &[SolvedFrame], o: &CalibrateOptions) -> Result<Value, String> {
    if frames.len() < CALIBRATE_MIN_FRAMES {
        return Err(format!(
            "calibration needs at least {CALIBRATE_MIN_FRAMES} solved frames, got {}",
            frames.len()
        ));
    }
    let frame = |i: usize| format!("solved frame {}", i + 1);

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

    let positive = frames
        .iter()
        .filter(|f| f.wcs.parity() == Parity::Positive)
        .count();
    let parity = match positive {
        0 => Some("negative"),
        n if n == frames.len() => Some("positive"),
        _ => None,
    };

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
            "offset_deg": round_decimal(offset, 3),
            "sign": sign,
            "sigma_deg": round_decimal(sigma, 3),
            "max_dev_deg": round_decimal((ROLL_LIMIT_FACTOR * sigma).max(ROLL_LIMIT_FLOOR_DEG), 3),
        }));
        roll_stats = json!({"robust_sigma_deg": spread, "worst_deg": worst});
    }

    let mut offsets = Vec::new();
    for f in frames {
        let size = (Some(f.width), Some(f.height));
        let header = if o.no_header_wcs {
            header_priors_for_size(&f.header.without_wcs(), size)
        } else {
            header_priors_for_size(&f.header, size)
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
        round_decimal((median(&offsets) / RAYLEIGH_MEDIAN).max(1e-4), 5)
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
    put("scale_arcsec", json!(round_decimal(scale, 5)));
    put("scale_sigma_arcsec", json!(round_decimal(scale_sigma, 5)));
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

#[derive(Debug, Clone, Copy)]
pub struct TransferStar {
    pub x: f64,
    pub y: f64,
    pub flux: f64,
    pub sigma: f64,
    pub residual: [f64; 2],
}

#[derive(Debug, Clone, Default)]
pub struct TransferFrame {
    pub stars: Vec<TransferStar>,
    pub width: usize,
    pub height: usize,
    pub epoch: Option<f64>,
}

pub fn transfer_frame(
    solver: &crate::Solver,
    image: &fayde_fits::Image,
    wcs: &Wcs,
    detect: &crate::DetectOptions,
) -> Result<TransferFrame, String> {
    let raw = crate::DetectOptions {
        charge_transfer: None,
        ..*detect
    };
    let sources = crate::detect_image(image, &raw, crate::DetectionReport::None)?.sources;
    let opts = crate::RequestOptions {
        header: Some(image.header.clone()),
        no_priors: true,
        epoch: detect.epoch,
        ..crate::RequestOptions::default()
    };
    let size = (
        u32::try_from(image.width).ok(),
        u32::try_from(image.height).ok(),
    );
    let built = crate::build_request(sources, size, &opts)?;
    let req = &built.request;
    let v = fayde_solve::verify_wcs(req, solver.catalog(), wcs).map_err(|e| e.to_string())?;
    let stars = v
        .matches
        .iter()
        .filter(|m| m.inlier_prob > 0.5)
        .filter_map(|m| {
            let s = &req.sources[m.source];
            Some(TransferStar {
                x: s.x,
                y: s.y,
                flux: s.flux?,
                sigma: s.sigma_px.unwrap_or(req.policy.default_sigma_px),
                residual: m.residual_px,
            })
        })
        .collect();
    Ok(TransferFrame {
        stars,
        width: image.width,
        height: image.height,
        epoch: req.policy.epoch.map(|e| e.year()),
    })
}

pub const TRANSFER_MIN_STARS: usize = 200;
const TRANSFER_MIN_SIGNIFICANCE: f64 = 5.0;
const TRANSFER_SIGMA_FLOOR: f64 = 0.03;
const TRANSFER_CLIP: f64 = 5.0;
const TRANSFER_MIN_SPAN_YEARS: f64 = 1.0;

fn solve_small(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    let n = b.len();
    for c in 0..n {
        let p = (c..n).max_by(|&i, &j| a[i][c].abs().total_cmp(&a[j][c].abs()))?;
        if a[p][c].abs() < 1e-12 * (1.0 + a[c][c].abs()) {
            return None;
        }
        a.swap(c, p);
        b.swap(c, p);
        let pivot = a[c].clone();
        for r in 0..n {
            if r != c {
                let f = a[r][c] / pivot[c];
                for (ark, pk) in a[r][c..].iter_mut().zip(&pivot[c..]) {
                    *ark -= f * pk;
                }
                b[r] -= f * b[c];
            }
        }
    }
    Some((0..n).map(|i| b[i] / a[i][i]).collect())
}

struct TransferFit {
    axis: &'static str,
    edge: &'static str,
    shift: (f64, f64),
    per_year: Option<(f64, f64)>,
    improvement: f64,
    stars: usize,
    frames: usize,
}

type Keep = Vec<Vec<bool>>;

fn fit_transfer(
    frames: &[TransferFrame],
    flux_ref: f64,
    axis: fayde_fits::TransferAxis,
    edge: fayde_fits::ReadoutEdge,
    t0: Option<f64>,
    fixed: Option<&Keep>,
) -> Option<(TransferFit, Keep)> {
    use fayde_fits::{ChargeTransfer, ReadoutEdge, TransferAxis};
    let unit = ChargeTransfer {
        axis,
        readout_edge: edge,
        shift_px: 1.0,
        flux_ref,
        per_year: 0.0,
        epoch: None,
    };
    let along = |s: &TransferStar| match axis {
        TransferAxis::X => s.residual[0],
        TransferAxis::Y => s.residual[1],
    };
    let weight =
        |s: &TransferStar| 1.0 / (s.sigma * s.sigma + TRANSFER_SIGMA_FLOOR * TRANSFER_SIGMA_FLOOR);
    let np = if t0.is_some() { 2 } else { 1 };
    let dt = |f: &TransferFrame| match (t0, f.epoch) {
        (Some(t0), Some(t)) => Some(t - t0),
        (Some(_), None) => None,
        (None, _) => Some(0.0),
    };
    let mut keep: Keep = match fixed {
        Some(k) => k.clone(),
        None => frames.iter().map(|f| vec![true; f.stars.len()]).collect(),
    };
    let rounds = if fixed.is_some() { 1 } else { 3 };
    let mut result = None;
    for round in 0..rounds {
        let mut gg = vec![vec![0.0; np]; np];
        let mut gr = vec![0.0; np];
        let mut rr = 0.0;
        let (mut used, mut frames_used) = (0usize, 0usize);
        let mut per_frame = Vec::with_capacity(frames.len());
        for (f, k) in frames.iter().zip(&keep) {
            let Some(dt) = dt(f) else {
                per_frame.push(None);
                continue;
            };
            let rows: Vec<([f64; 3], Vec<f64>, f64, f64)> = f
                .stars
                .iter()
                .zip(k)
                .filter(|(_, &k)| k)
                .map(|(s, _)| {
                    let g = unit.displacement(s.x, s.y, f.width, f.height, s.flux, 1.0);
                    let gv = if np == 2 { vec![g, g * dt] } else { vec![g] };
                    ([1.0, s.x / 1000.0, s.y / 1000.0], gv, along(s), weight(s))
                })
                .collect();
            if rows.len() < 6 {
                per_frame.push(None);
                continue;
            }
            let mut m = vec![vec![0.0; 3]; 3];
            let mut b = vec![vec![0.0; np]; 3];
            let mut c = vec![0.0; 3];
            let mut fgg = vec![vec![0.0; np]; np];
            let mut fgr = vec![0.0; np];
            let mut frr = 0.0;
            for (n, g, r, w) in &rows {
                for i in 0..3 {
                    for j in 0..3 {
                        m[i][j] += w * n[i] * n[j];
                    }
                    for (bij, gj) in b[i].iter_mut().zip(g) {
                        *bij += w * n[i] * gj;
                    }
                    c[i] += w * n[i] * r;
                }
                for i in 0..np {
                    for j in 0..np {
                        fgg[i][j] += w * g[i] * g[j];
                    }
                    fgr[i] += w * g[i] * r;
                }
                frr += w * r * r;
            }
            let Some(mc) = solve_small(m.clone(), c.clone()) else {
                per_frame.push(None);
                continue;
            };
            let Some(mb) = (0..np)
                .map(|j| solve_small(m.clone(), (0..3).map(|i| b[i][j]).collect()))
                .collect::<Option<Vec<Vec<f64>>>>()
            else {
                per_frame.push(None);
                continue;
            };
            for i in 0..np {
                for (j, mbj) in mb.iter().enumerate() {
                    gg[i][j] += fgg[i][j] - (0..3).map(|k| b[k][i] * mbj[k]).sum::<f64>();
                }
                gr[i] += fgr[i] - (0..3).map(|k| b[k][i] * mc[k]).sum::<f64>();
            }
            rr += frr - (0..3).map(|k| c[k] * mc[k]).sum::<f64>();
            used += rows.len();
            frames_used += 1;
            per_frame.push(Some((mc, mb, dt)));
        }
        let beta = solve_small(gg.clone(), gr.clone())?;
        let chi2 = rr - (0..np).map(|i| beta[i] * gr[i]).sum::<f64>();
        let dof = used.saturating_sub(3 * frames_used + np).max(1) as f64;
        let s2 = (chi2 / dof).max(1e-12);
        let sig = |i: usize| {
            let mut e = vec![0.0; np];
            e[i] = 1.0;
            solve_small(gg.clone(), e).map(|col| (col[i] * s2).sqrt())
        };
        result = Some(TransferFit {
            axis: match axis {
                TransferAxis::X => "x",
                TransferAxis::Y => "y",
            },
            edge: match edge {
                ReadoutEdge::Low => "low",
                ReadoutEdge::High => "high",
            },
            shift: (beta[0], sig(0)?),
            per_year: (np == 2).then(|| (beta[1], sig(1).unwrap_or(f64::NAN))),
            improvement: (0..np).map(|i| beta[i] * gr[i]).sum::<f64>(),
            stars: used,
            frames: frames_used,
        });
        if round + 1 == rounds {
            break;
        }
        let mut z: Vec<(usize, usize, f64)> = Vec::new();
        for (fi, (f, pf)) in frames.iter().zip(&per_frame).enumerate() {
            let Some((mc, mb, dt)) = pf else { continue };
            let own: Vec<f64> = (0..3)
                .map(|i| mc[i] - (0..np).map(|j| mb[j][i] * beta[j]).sum::<f64>())
                .collect();
            for (si, s) in f.stars.iter().enumerate().filter(|(si, _)| keep[fi][*si]) {
                let g = unit.displacement(s.x, s.y, f.width, f.height, s.flux, 1.0);
                let model = own[0]
                    + own[1] * s.x / 1000.0
                    + own[2] * s.y / 1000.0
                    + beta[0] * g
                    + if np == 2 { beta[1] * g * dt } else { 0.0 };
                z.push((fi, si, (along(s) - model) * weight(s).sqrt()));
            }
        }
        let values: Vec<f64> = z.iter().map(|v| v.2).collect();
        if values.is_empty() {
            break;
        }
        let (_, spread) = robust(&values);
        for (fi, si, v) in z {
            keep[fi][si] = v.abs() <= TRANSFER_CLIP * spread.max(1e-9);
        }
    }
    result.map(|r| (r, keep))
}

pub fn calibrate_charge_transfer(
    frames: &[TransferFrame],
) -> Result<(Option<Value>, Value), String> {
    use fayde_fits::{ReadoutEdge, TransferAxis};
    let mut fluxes: Vec<f64> = frames
        .iter()
        .flat_map(|f| f.stars.iter().map(|s| s.flux))
        .filter(|f| f.is_finite() && *f > 0.0)
        .collect();
    if fluxes.len() < TRANSFER_MIN_STARS {
        return Err(format!(
            "a charge-transfer calibration needs at least {TRANSFER_MIN_STARS} matched stars, got {}",
            fluxes.len()
        ));
    }
    fluxes.sort_by(f64::total_cmp);
    let flux_ref = rounded_sig(
        fluxes[((fluxes.len() - 1) as f64 * 0.99).round() as usize],
        3,
    );
    let epochs: Vec<f64> = frames
        .iter()
        .filter(|f| f.stars.len() >= 6)
        .filter_map(|f| f.epoch)
        .collect();
    let span = epochs.iter().copied().fold(f64::NEG_INFINITY, f64::max)
        - epochs.iter().copied().fold(f64::INFINITY, f64::min);
    let t0 = (epochs.len() >= CALIBRATE_MIN_FRAMES && span >= TRANSFER_MIN_SPAN_YEARS)
        .then(|| round_decimal(median(&epochs), 2));
    let geometries = [
        (TransferAxis::Y, ReadoutEdge::Low),
        (TransferAxis::Y, ReadoutEdge::High),
        (TransferAxis::X, ReadoutEdge::Low),
        (TransferAxis::X, ReadoutEdge::High),
    ];
    let fit = |axis, edge, fixed: Option<&Keep>| {
        fit_transfer(frames, flux_ref, axis, edge, t0, fixed)
            .map(|r| (r, t0))
            .or_else(|| fit_transfer(frames, flux_ref, axis, edge, None, fixed).map(|r| (r, None)))
    };
    let qualifies =
        |f: &TransferFit| f.shift.0 > 0.0 && f.shift.0 >= TRANSFER_MIN_SIGNIFICANCE * f.shift.1;
    let better = |a: &&((TransferFit, Keep), Option<f64>),
                  b: &&((TransferFit, Keep), Option<f64>)| {
        (qualifies(&a.0 .0), a.0 .0.improvement)
            .partial_cmp(&(qualifies(&b.0 .0), b.0 .0.improvement))
            .unwrap_or(std::cmp::Ordering::Equal)
    };
    let first: Vec<_> = geometries
        .iter()
        .filter_map(|&(a, e)| fit(a, e, None))
        .collect();
    let lead = first
        .iter()
        .max_by(better)
        .ok_or("the charge-transfer fit failed: too few stars per frame")?;
    let common = lead.0 .1.clone();
    let fits: Vec<_> = geometries
        .iter()
        .filter_map(|&(a, e)| fit(a, e, Some(&common)))
        .collect();
    let best = fits
        .iter()
        .max_by(better)
        .ok_or("the charge-transfer fit failed: too few stars per frame")?;
    let stats = json!({
        "stars": best.0 .0.stars,
        "frames": best.0 .0.frames,
        "flux_ref": flux_ref,
        "fits": fits.iter().map(|((f, _), _)| json!({
            "axis": f.axis,
            "readout_edge": f.edge,
            "shift_px": f.shift.0,
            "shift_sigma_px": f.shift.1,
            "per_year": f.per_year.map(|p| p.0),
            "per_year_sigma": f.per_year.map(|p| p.1),
            "chi2_improvement": f.improvement,
        })).collect::<Vec<_>>(),
    });
    let ((best, _), best_t0) = best;
    if !qualifies(best) {
        return Ok((None, stats));
    }
    let mut term = json!({
        "axis": best.axis,
        "readout_edge": best.edge,
        "shift_px": round_decimal(best.shift.0, 4),
        "flux_ref": flux_ref,
    });
    if let (Some((rate, _)), Some(t0)) = (best.per_year, best_t0) {
        term["per_year"] = json!(round_decimal(rate, 4));
        term["epoch"] = json!(t0);
    }
    Ok((Some(term), stats))
}

fn rounded_sig(v: f64, digits: i32) -> f64 {
    if v == 0.0 || !v.is_finite() {
        return v;
    }
    let scale = 10f64.ppowi(digits - 1 - v.abs().plog10().floor() as i32);
    (v * scale).round() / scale
}

#[cfg(test)]
mod tests {
    use super::*;
    use fayde_core::SkyPos;
    use fayde_fits::InstrumentProfile;
    use fayde_wcs::TanWcs;

    fn frame(k: usize, bin: f64, scale_err: f64, roll_err: f64, off: f64) -> SolvedFrame {
        let key = -170.0 + 37.0 * k as f64;
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
        let p = InstrumentProfile::from_json(&json.to_string()).unwrap();
        assert_eq!(p.name.as_deref(), Some("test camera"));
        assert!((p.scale_arcsec.unwrap() - 2.5).abs() < 5e-4, "{json}");
        assert!(
            (p.scale_sigma_arcsec.unwrap() - 0.0025).abs() < 1e-9,
            "{json}"
        );
        assert_eq!(p.parity.as_deref(), Some("positive"));
        let roll = p.roll.unwrap();
        assert_eq!((roll.keyword.as_str(), roll.sign), ("ROLLKEY", -1.0));
        assert!((roll.offset_deg - 33.0).abs() < 0.03, "{json}");
        assert!(
            (roll.sigma_deg - 0.089).abs() < 0.005 && roll.max_dev_deg == 1.0,
            "{json}"
        );
        assert!(
            (p.pointing_sigma_deg.unwrap() - 0.018 / 1.1774).abs() < 1e-4,
            "{json}"
        );
        assert_eq!(p.header_wcs, Some(false));
        assert_eq!(json["binning"]["x_keyword"], "XBINNING");
        assert_eq!(json["calibration"]["frames"], 9);
        assert!((json["calibration"]["roll"]["worst_deg"].as_f64().unwrap() - 4.0).abs() < 0.05);

        let with_wcs = CalibrateOptions {
            no_header_wcs: false,
            ..options()
        };
        let json = calibrate_profile(&frames, &with_wcs).unwrap();
        assert!((json["pointing_sigma_deg"].as_f64().unwrap() - 1.0 / 1.1774).abs() < 1e-3);
        assert!(json.get("header_wcs").is_none());
    }

    fn smeared(
        axis: fayde_fits::TransferAxis,
        edge: fayde_fits::ReadoutEdge,
        shift: f64,
        per_year: f64,
    ) -> Vec<TransferFrame> {
        use fayde_detect::synth::Rng;
        use fayde_fits::{ChargeTransfer, TransferAxis};
        let model = ChargeTransfer {
            axis,
            readout_edge: edge,
            shift_px: shift,
            flux_ref: 100_000.0,
            per_year,
            epoch: Some(2021.0),
        };
        let mut rng = Rng(7);
        (0..10)
            .map(|k| {
                let t = 2020.0 + 0.4 * k as f64;
                let s = model.shift_at(Some(t));
                let (a, b, c) = (0.3 * rng.normal(), 0.2 * rng.normal(), 0.2 * rng.normal());
                let stars = (0..300)
                    .map(|_| {
                        let (x, y) = (1024.0 * rng.uniform(), 768.0 * rng.uniform());
                        let flux = 300.0 * (400.0f64).ppowf(rng.uniform());
                        let sigma = (0.4 * (500.0 / flux).sqrt()).max(0.05);
                        let d = model.displacement(x, y, 1024, 768, flux, s);
                        let wcs = a + b * x / 1000.0 + c * y / 1000.0;
                        let (on, off) = (d + wcs + sigma * rng.normal(), sigma * rng.normal());
                        let residual = match axis {
                            TransferAxis::X => [on, off],
                            TransferAxis::Y => [off, on],
                        };
                        TransferStar {
                            x,
                            y,
                            flux,
                            sigma,
                            residual,
                        }
                    })
                    .collect();
                TransferFrame {
                    stars,
                    width: 1024,
                    height: 768,
                    epoch: Some(t),
                }
            })
            .collect()
    }

    #[test]
    fn recovers_the_charge_transfer_smear_and_its_geometry() {
        use fayde_fits::{ReadoutEdge, TransferAxis};
        let frames = smeared(TransferAxis::Y, ReadoutEdge::Low, 0.2, 0.02);
        let (term, stats) = calibrate_charge_transfer(&frames).unwrap();
        let term = term.expect("a smear this strong is found");
        assert_eq!(
            (term["axis"].as_str(), term["readout_edge"].as_str()),
            (Some("y"), Some("low")),
            "{stats}"
        );
        assert!(
            (term["shift_px"].as_f64().unwrap()
                - 0.2
                - 0.02 * (term["epoch"].as_f64().unwrap() - 2021.0))
                .abs()
                < 0.01,
            "{term}"
        );
        assert!(
            (term["per_year"].as_f64().unwrap() - 0.02).abs() < 0.006,
            "{term}"
        );
        let p =
            InstrumentProfile::from_json(&json!({ "charge_transfer": term }).to_string()).unwrap();
        assert!(p.charge_transfer.is_some());

        let (term, _) =
            calibrate_charge_transfer(&smeared(TransferAxis::X, ReadoutEdge::High, 0.2, 0.0))
                .unwrap();
        let term = term.unwrap();
        assert_eq!(
            (term["axis"].as_str(), term["readout_edge"].as_str()),
            (Some("x"), Some("high")),
            "{term}"
        );

        let (term, stats) =
            calibrate_charge_transfer(&smeared(TransferAxis::Y, ReadoutEdge::Low, 0.0, 0.0))
                .unwrap();
        assert!(term.is_none(), "{stats}");
        assert_eq!(stats["fits"].as_array().unwrap().len(), 4);

        let base = calibrate_charge_transfer(&frames).unwrap().0.unwrap();
        let mut extra = frames.clone();
        let mut line = frames[0].clone();
        for s in &mut line.stars {
            s.y = 100.0;
            s.residual[1] += 5.0 * s.x / 1000.0;
        }
        extra.push(line);
        extra.extend((0..100).map(|_| TransferFrame {
            stars: Vec::new(),
            ..frames[0].clone()
        }));
        let mut undated = smeared(TransferAxis::Y, ReadoutEdge::Low, 0.6, 0.0)[0].clone();
        undated.epoch = None;
        extra.push(undated);
        let (term, stats) = calibrate_charge_transfer(&extra).unwrap();
        let term = term.unwrap();
        let at = |t: &Value, year: f64| {
            t["shift_px"].as_f64().unwrap()
                + t["per_year"].as_f64().unwrap() * (year - t["epoch"].as_f64().unwrap())
        };
        assert!(
            (at(&term, 2022.0) - at(&base, 2022.0)).abs() < 1e-3,
            "{term} against {base}; {stats}"
        );
        assert_eq!(term["per_year"], base["per_year"], "{term} against {base}");
        assert_eq!(stats["frames"], 10, "{stats}");

        let mut wrong = frames.clone();
        let towards = fayde_fits::ChargeTransfer {
            axis: TransferAxis::X,
            readout_edge: ReadoutEdge::Low,
            shift_px: 0.5,
            flux_ref: 100_000.0,
            per_year: 0.0,
            epoch: None,
        };
        for f in &mut wrong {
            for s in &mut f.stars {
                s.residual[0] -= towards.displacement(s.x, s.y, f.width, f.height, s.flux, 0.5);
            }
        }
        let (term, stats) = calibrate_charge_transfer(&wrong).unwrap();
        let best = stats["fits"]
            .as_array()
            .unwrap()
            .iter()
            .max_by(|a, b| {
                a["chi2_improvement"]
                    .as_f64()
                    .unwrap()
                    .total_cmp(&b["chi2_improvement"].as_f64().unwrap())
            })
            .unwrap();
        assert!(best["shift_px"].as_f64().unwrap() < 0.0, "{stats}");
        let term = term.unwrap_or_else(|| panic!("{stats}"));
        assert!(term["shift_px"].as_f64().unwrap() > 0.0, "{term}");

        let few = vec![TransferFrame {
            stars: frames[0].stars[..50].to_vec(),
            ..frames[0].clone()
        }];
        let err = calibrate_charge_transfer(&few).unwrap_err();
        assert!(err.contains("at least 200"), "{err}");
    }

    #[test]
    fn says_what_it_cannot_calibrate() {
        let frames: Vec<SolvedFrame> = (0..6).map(|k| frame(k, 1.0, 0.0, 0.0, 0.01)).collect();
        let err = calibrate_profile(&frames[..4], &options()).unwrap_err();
        assert!(err.contains("at least 5"), "{err}");
        let mut shuffled = frames.clone();
        for (k, f) in shuffled.iter_mut().enumerate() {
            f.header.set("ROLLKEY", Card::Float(71.0 * (k * k) as f64));
        }
        let err = calibrate_profile(&shuffled, &options()).unwrap_err();
        assert!(err.contains("does not give the roll"), "{err}");
        let mut missing = frames.clone();
        missing[2].header = Header::default();
        let err = calibrate_profile(&missing, &options()).unwrap_err();
        assert!(err.contains("solved frame 3"), "{err}");
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
