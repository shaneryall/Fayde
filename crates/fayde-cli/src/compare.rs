//! `fayde compare`: solve exported producer source lists frame by frame, with priors and blind,
//! and compare the two solutions; where the file carries another solver's WCS, compare with
//! that too.
//!
//! Input: JSON with a `frames` array. Each frame has `name`, `width`, `height`, `sources_fits`
//! (`[x, y]` pairs, FITS 1-based, brightest first), `header` (keyword → value, every valued card
//! of the frame header) and, optionally, `solve_field_wcs_cards` (another solver's WCS as card
//! images). Each frame is solved twice:
//! - **with priors:** as `fayde solve --header` does (`fayde_api::build_request`): pointing,
//!   scale and epoch from the frame header, or from `--instrument` where the profile has them;
//!   pointing 1σ 0.5° and hard limit 2°, scale 1σ 2%, unless the profile says otherwise;
//! - **blind:** no priors, pattern index only (skipped with `--hinted-only`).
//!
//! Fayde is checked against itself: where both runs solve, their WCS are compared on a 5×5
//! grid over the frame. The other solver's WCS, where given, is a second opinion, not an
//! input: with `--no-header-wcs` nothing a plate solver wrote into the header is read either
//! (the pointing then comes from `RA`/`DEC` or `OBJCTRA`/`OBJCTDEC`).

use fayde_core::Portable;
use std::time::Instant;

use fayde_api::{build_request, RequestOptions};
use fayde_catalog::ucac4::Ucac4Catalog;
use fayde_core::{to_arcsec, PixelConvention, PixelPos};
use fayde_fits::InstrumentProfile;
use fayde_index::Index;
use fayde_solve::{solve_auto, Outcome, Priors, SolveRequest, Source};
use fayde_wcs::{Header, Wcs};
use serde_json::Value;

use crate::{index_paths, num, options, seconds};

struct Run {
    /// The solution's WCS.
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
            // Arbiter: both WCS against the same confident matches (source pixel vs catalogue star).
            let rms_of = |wcs: &Wcs| {
                let (mut s, mut n) = (0.0f64, 0.0f64);
                for m in sol.matches.iter().filter(|m| m.inlier_prob > 0.5) {
                    if let Some(p) = wcs.sky_to_pixel(m.star_pos) {
                        let src = &req.sources[m.source];
                        s += (src.x - p.x).ppowi(2) + (src.y - p.y).ppowi(2);
                        n += 1.0;
                    }
                }
                (s / n.max(1.0)).sqrt()
            };
            let ref_rms = reference.map(rms_of);
            let text = format!(
                "{:?} order {} | {} matches rms {:.2}px escale {:.2} σ(centre) {:.2}″ NFA {:.0e} | {}",
                sol.method,
                sol.wcs.sip_order(),
                sol.matches.iter().filter(|m| m.inlier_prob > 0.5).count(),
                sol.rms_px,
                sol.error_scale,
                sol.covariance.radial_sigma_arcsec(sol.wcs.crpix()),
                sol.nfa,
                match (agree, ref_rms) {
                    (Some(a), Some(r)) => format!("Δ vs ref {a:.2}″; same matches: fayde rms {:.2}px, ref rms {r:.2}px", rms_of(&sol.wcs)),
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

/// Largest sky separation between two WCS over a 5×5 grid, arcsec.
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

fn f(v: &Value, k: &str) -> Option<f64> {
    v.get(k).and_then(Value::as_f64)
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
    // Fayde against itself: frames both runs solve, and how far apart the two solutions are.
    let (mut both, mut both_agree) = (0usize, 0usize);
    let mut own_agreements = Vec::new();
    let mut agreements = Vec::new();
    let mut times = (Vec::new(), Vec::new());
    for fr in frames.iter().take(limit) {
        let (w, h) = (
            f(fr, "width").unwrap_or(0.0),
            f(fr, "height").unwrap_or(0.0),
        );
        let hdr = &fr["header"];
        let sources: Vec<Source> = fr["sources_fits"]
            .as_array()
            .ok_or("sources")?
            .iter()
            .map(|p| {
                let pos = PixelConvention::Fits.to_internal(
                    p[0].as_f64().unwrap_or(f64::NAN),
                    p[1].as_f64().unwrap_or(f64::NAN),
                );
                Source::at(pos.x, pos.y)
            })
            .collect();
        let reference = fr["solve_field_wcs_cards"].as_array().and_then(|cards| {
            let cards: Vec<&str> = cards.iter().filter_map(Value::as_str).collect();
            Wcs::from_header(&Header::parse_cards(cards)).ok()
        });
        // The request `fayde solve --header [--instrument] [--no-header-wcs]` builds.
        let header = hdr
            .as_object()
            .map(fayde_api::header_from_json)
            .transpose()?;
        let built = build_request(
            sources,
            (Some(w as u32), Some(h as u32)),
            &RequestOptions {
                header,
                instrument: instrument.clone(),
                no_header_wcs,
                deadline,
                ..RequestOptions::default()
            },
        )?;
        let hinted = built.request;
        // Blind: the same frame, sources and epoch, and no prior.
        let base = SolveRequest {
            priors: Priors::default(),
            ..hinted.clone()
        };
        n += 1;
        ref_solved += usize::from(reference.is_some());
        let (a, ta) = run_one(&hinted, &cat, &indexes, reference.as_ref(), w, h)?;
        let (b, tb) = if hinted_only {
            (
                Run {
                    wcs: None,
                    agree_arcsec: None,
                    seconds: 0.0,
                },
                "skipped (--hinted-only)".to_string(),
            )
        } else {
            run_one(&base, &cat, &indexes, reference.as_ref(), w, h)?
        };
        times.0.push(a.seconds);
        if !hinted_only {
            times.1.push(b.seconds);
        }
        let own = match (&a.wcs, &b.wcs) {
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
                &a,
                &mut hp_solved,
                &mut hp_agree,
                &mut hp_disagree,
                &mut hp_only,
                &mut ref_only_hp,
            ),
            (
                &b,
                &mut bl_solved,
                &mut bl_agree,
                &mut bl_disagree,
                &mut bl_only,
                &mut ref_only_bl,
            ),
        ] {
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
        eprintln!(
            "{:<42} {:>3} src  reference {:<6}\n    with priors   {:6.2}s  {ta}\n    blind         {:6.2}s  {tb}\n    {own}",
            fr["name"].as_str().unwrap_or("?"),
            base.sources.len(),
            if reference.is_some() { "given" } else { "none" },
            a.seconds,
            b.seconds,
        );
    }
    agreements.sort_by(f64::total_cmp);
    own_agreements.sort_by(f64::total_cmp);
    let med = |v: &mut Vec<f64>| {
        v.sort_by(f64::total_cmp);
        v.get(v.len() / 2).copied()
    };
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
        own_agreements.get(own_agreements.len() / 2).copied().into(),
    );
    stats.insert(
        "priors_and_blind_arcsec_max".into(),
        own_agreements.last().copied().into(),
    );
    stats.insert(
        "agreement_arcsec_median".into(),
        agreements.get(agreements.len() / 2).copied().into(),
    );
    stats.insert(
        "agreement_arcsec_max".into(),
        agreements.last().copied().into(),
    );
    stats.insert(
        "header_priors_seconds_median".into(),
        med(&mut times.0).into(),
    );
    stats.insert("blind_seconds_median".into(), med(&mut times.1).into());
    crate::out(&serde_json::to_string_pretty(&Value::Object(stats)).expect("json"))
}
