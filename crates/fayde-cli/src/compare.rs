use fayde_core::Portable;
use std::time::Instant;

use fayde_api::{build_request, RequestOptions};
use fayde_catalog::ucac4::Ucac4Catalog;
use fayde_core::{to_arcsec, PixelConvention, PixelPos};
use fayde_fits::InstrumentProfile;
use fayde_index::Index;
use fayde_solve::{solve_auto, Outcome, Priors, SolveRequest, Source};
use fayde_wcs::{Header, Wcs};
use serde::Deserialize;
use serde_json::{Map, Value};

use crate::{index_paths, num, options, seconds};

#[derive(Deserialize)]
struct FrameInput<'a> {
    width: f64,
    height: f64,
    sources_fits: Vec<[f64; 2]>,
    header: Option<Map<String, Value>>,
    #[serde(borrow)]
    solve_field_wcs_cards: Option<Vec<&'a str>>,
}

struct Run {
    wcs: Option<Wcs>,
    agree_arcsec: Option<f64>,
    seconds: f64,
}

fn run_one(
    req: &SolveRequest,
    cat: &Ucac4Catalog,
    indexes: &[Index],
    reference: Option<&Wcs>,
    w: f64,
    h: f64,
) -> Result<(Run, String), String> {
    let t0 = Instant::now();
    let out = solve_auto(req, cat, indexes).map_err(|e| e.to_string())?;
    let seconds = t0.elapsed().as_secs_f64();
    Ok(match out {
        Outcome::Solved(sol) => {
            let agree = reference.map(|r| max_disagreement(&sol.wcs, r, w, h));
            let rms_of = |wcs: &Wcs| {
                let mut sum = 0.0;
                let mut count = 0;
                for m in sol.matches.iter().filter(|m| m.inlier_prob > 0.5) {
                    let p = wcs.sky_to_pixel(m.star_pos)?;
                    let src = &req.sources[m.source];
                    sum += (src.x - p.x).ppowi(2) + (src.y - p.y).ppowi(2);
                    count += 1;
                }
                (count > 0).then(|| (sum / f64::from(count)).sqrt())
            };
            let text = format!(
                "{:?} order {} | {} matches rms {:.2}px escale {:.2} σ(centre) {:.2}″ NFA {:.0e} | {}",
                sol.method,
                sol.wcs.sip_order(),
                sol.matches.iter().filter(|m| m.inlier_prob > 0.5).count(),
                sol.rms_px,
                sol.error_scale,
                sol.covariance.radial_sigma_arcsec(sol.wcs.crpix()),
                sol.nfa,
                match (agree, reference) {
                    (Some(a), Some(reference)) => match (rms_of(&sol.wcs), rms_of(reference)) {
                        (Some(f), Some(r)) => format!("Δ vs ref {a:.2}″; same matches: fayde rms {f:.2}px, ref rms {r:.2}px"),
                        _ => format!("Δ vs ref {a:.2}″; matched-pair RMS unavailable"),
                    },
                    _ => "no reference".into(),
                },
            );
            (
                Run {
                    wcs: Some(sol.wcs.clone()),
                    agree_arcsec: agree,
                    seconds,
                },
                text,
            )
        }
        Outcome::Abstained {
            reason,
            hypotheses_tested,
            ..
        } => (
            Run {
                wcs: None,
                agree_arcsec: None,
                seconds,
            },
            format!("abstained: {reason} ({hypotheses_tested} hyp)"),
        ),
    })
}

fn max_disagreement(wcs: &Wcs, reference: &Wcs, w: f64, h: f64) -> f64 {
    let mut worst = 0.0f64;
    for i in 0..5 {
        for j in 0..5 {
            let p = PixelPos::new(-0.5 + w * i as f64 / 4.0, -0.5 + h * j as f64 / 4.0);
            worst = worst.max(to_arcsec(
                wcs.pixel_to_sky(p).separation(reference.pixel_to_sky(p)),
            ));
        }
    }
    worst
}

fn median(sorted: &[f64]) -> Option<f64> {
    let last = sorted.len().checked_sub(1)?;
    Some(sorted[last / 2].midpoint(sorted[sorted.len() / 2]))
}

pub(crate) fn run(args: &[&str]) -> Result<(), String> {
    let o = options(
        args,
        &[
            "ucac4",
            "index",
            "frames",
            "limit",
            "deadline",
            "instrument",
        ],
        &["no-header-wcs", "hinted-only"],
    )?;
    let hinted_only = o.contains_key("hinted-only");
    let no_header_wcs = o.contains_key("no-header-wcs");
    let instrument = o
        .get("instrument")
        .map(|p| InstrumentProfile::load(p).map_err(|e| format!("{p}: {e}")))
        .transpose()?;
    let deadline = Some(seconds(&o, "deadline")?.unwrap_or(std::time::Duration::from_secs(60)));
    let cat = Ucac4Catalog::open(o.get("ucac4").ok_or("--ucac4 is required")?)
        .map_err(|e| e.to_string())?;
    let paths = index_paths(&o);
    if paths.is_empty() && !hinted_only {
        return Err("--index is required (the blind run needs one; or give --hinted-only)".into());
    }
    let indexes: Vec<Index> = paths
        .iter()
        .map(|p| Index::read(p).map_err(|e| e.to_string()))
        .collect::<Result<Vec<Index>, String>>()?;
    let text = std::fs::read_to_string(o.get("frames").ok_or("--frames is required")?)
        .map_err(|e| e.to_string())?;
    let data: Value = serde_json::from_str(&text).map_err(|e| e.to_string())?;
    let frames = data["frames"].as_array().ok_or("no frames")?;
    let limit: usize = num(&o, "limit")?.unwrap_or(usize::MAX);

    let mut stats = serde_json::Map::new();
    let (mut n, mut ref_solved) = (0usize, 0usize);
    let (mut hp_solved, mut hp_agree, mut hp_disagree, mut hp_only) = (0, 0, 0, 0);
    let (mut bl_solved, mut bl_agree, mut bl_disagree, mut bl_only) = (0, 0, 0, 0);
    let (mut ref_only_hp, mut ref_only_bl) = (0, 0);
    let (mut both, mut both_agree) = (0usize, 0usize);
    let mut own_agreements = Vec::new();
    let mut agreements = Vec::new();
    let mut times = (Vec::new(), Vec::new());
    for fr in frames.iter().take(limit) {
        let name = fr["name"].as_str().unwrap_or("unnamed frame");
        let frame = FrameInput::deserialize(fr).map_err(|e| format!("{name}: {e}"))?;
        let dimension = |v: f64, key: &str| {
            (v.fract() == 0.0 && (1.0..=f64::from(u32::MAX)).contains(&v))
                .then_some(v as u32)
                .ok_or_else(|| format!("{name}: {key} must be a positive integer image dimension"))
        };
        let (width, height) = (
            dimension(frame.width, "width")?,
            dimension(frame.height, "height")?,
        );
        let (w, h) = (f64::from(width), f64::from(height));
        let sources = frame
            .sources_fits
            .into_iter()
            .map(|[x, y]| {
                let pos = PixelConvention::Fits.to_internal(x, y);
                Source::at(pos.x, pos.y)
            })
            .collect();
        let reference = frame
            .solve_field_wcs_cards
            .filter(|cards| !cards.is_empty())
            .map(|cards| Wcs::from_header(&Header::parse_cards(cards)))
            .transpose()
            .map_err(|e| format!("{name}: reference WCS: {e}"))?;
        let header = frame
            .header
            .as_ref()
            .map(fayde_api::header_from_json)
            .transpose()?;
        let built = build_request(
            sources,
            (Some(width), Some(height)),
            &RequestOptions {
                header,
                instrument: instrument.clone(),
                no_header_wcs,
                deadline,
                ..RequestOptions::default()
            },
        )?;
        let hinted = built.request;
        n += 1;
        ref_solved += usize::from(reference.is_some());
        let (a, ta) = run_one(&hinted, &cat, &indexes, reference.as_ref(), w, h)?;
        let (b, tb) = if hinted_only {
            (None, "skipped (--hinted-only)".to_string())
        } else {
            let base = SolveRequest {
                priors: Priors::default(),
                ..hinted.clone()
            };
            let (run, text) = run_one(&base, &cat, &indexes, reference.as_ref(), w, h)?;
            (Some(run), text)
        };
        times.0.push(a.seconds);
        if let Some(b) = &b {
            times.1.push(b.seconds);
        }
        let own = match (&a.wcs, b.as_ref().and_then(|b| b.wcs.as_ref())) {
            (Some(x), Some(y)) => {
                let d = max_disagreement(x, y, w, h);
                both += 1;
                both_agree += usize::from(d < 5.0);
                own_agreements.push(d);
                format!("with priors and blind {d:.2}″ apart")
            }
            _ => "with priors and blind: not both solved".to_string(),
        };
        for (r, solved, agree, disagree, only, ref_only) in [
            (
                Some(&a),
                &mut hp_solved,
                &mut hp_agree,
                &mut hp_disagree,
                &mut hp_only,
                &mut ref_only_hp,
            ),
            (
                b.as_ref(),
                &mut bl_solved,
                &mut bl_agree,
                &mut bl_disagree,
                &mut bl_only,
                &mut ref_only_bl,
            ),
        ] {
            let Some(r) = r else { continue };
            *solved += usize::from(r.wcs.is_some());
            match (r.wcs.is_some(), reference.is_some(), r.agree_arcsec) {
                (true, true, Some(d)) if d < 5.0 => *agree += 1,
                (true, true, _) => *disagree += 1,
                (true, false, _) => *only += 1,
                (false, true, _) => *ref_only += 1,
                _ => {}
            }
        }
        if let Some(d) = a.agree_arcsec {
            agreements.push(d);
        }
        let blind_text = match b {
            Some(b) => format!("{:6.2}s  {tb}", b.seconds),
            None => tb,
        };
        eprintln!(
            "{:<42} {:>3} src  reference {:<6}\n    with priors   {:6.2}s  {ta}\n    blind         {blind_text}\n    {own}",
            name,
            hinted.sources.len(),
            if reference.is_some() { "given" } else { "none" },
            a.seconds,
        );
    }
    agreements.sort_by(f64::total_cmp);
    own_agreements.sort_by(f64::total_cmp);
    times.0.sort_by(f64::total_cmp);
    times.1.sort_by(f64::total_cmp);
    for (k, v) in [
        ("frames", n),
        ("priors_and_blind_both_solved", both),
        ("priors_and_blind_agree_lt_5arcsec", both_agree),
        ("priors_and_blind_disagree", both - both_agree),
        ("solve_field_solved", ref_solved),
        ("header_priors_solved", hp_solved),
        ("header_priors_agree_lt_5arcsec", hp_agree),
        ("header_priors_disagree", hp_disagree),
        ("header_priors_solved_where_solve_field_failed", hp_only),
        ("header_priors_failed_where_solve_field_solved", ref_only_hp),
        ("blind_solved", bl_solved),
        ("blind_agree_lt_5arcsec", bl_agree),
        ("blind_disagree", bl_disagree),
        ("blind_solved_where_solve_field_failed", bl_only),
        ("blind_failed_where_solve_field_solved", ref_only_bl),
    ] {
        stats.insert(k.into(), v.into());
    }
    stats.insert(
        "priors_and_blind_arcsec_median".into(),
        median(&own_agreements).into(),
    );
    stats.insert(
        "priors_and_blind_arcsec_max".into(),
        own_agreements.last().copied().into(),
    );
    stats.insert("agreement_arcsec_median".into(), median(&agreements).into());
    stats.insert(
        "agreement_arcsec_max".into(),
        agreements.last().copied().into(),
    );
    stats.insert(
        "header_priors_seconds_median".into(),
        median(&times.0).into(),
    );
    stats.insert("blind_seconds_median".into(), median(&times.1).into());
    crate::out(&serde_json::to_string_pretty(&Value::Object(stats)).expect("json"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reported_median_uses_both_middle_values() {
        assert_eq!(median(&[]), None);
        assert_eq!(median(&[1.0, 7.0]), Some(4.0));
        assert_eq!(median(&[1.0, 5.0, 7.0]), Some(5.0));
    }
}
