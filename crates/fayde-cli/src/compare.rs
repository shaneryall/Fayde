//! `fayde compare`: solve exported producer source lists and compare with an existing solver's
//! WCS (e.g. solve-field), frame by frame.
//!
//! Input: JSON with a `frames` array. Each frame has `name`, `width`, `height`, `sources_fits`
//! (`[x, y]` pairs, FITS 1-based, brightest first), `header` (keyword → value, every valued card
//! of the frame header) and, where the other solver solved it, `solve_field_wcs_cards` (its WCS
//! as card images). Each frame is solved twice:
//! - **header priors:** pointing, scale and epoch read from the frame header as `fayde solve
//!   --header` reads them (`fayde_fits::header_priors`); pointing 1σ 0.5°, hard limit 2°, scale
//!   1σ 2%, each used when present — the information a pipeline has;
//! - **blind:** no priors, pattern index only.
//!
//! Where the other solver produced a WCS, the two are compared on a 5×5 grid over the frame.

use fayde_core::Portable;
use std::time::Instant;

use fayde_catalog::ucac4::Ucac4Catalog;
use fayde_core::{to_arcsec, PixelConvention, PixelPos};
use fayde_fits::header_priors;
use fayde_index::Index;
use fayde_solve::{
    solve_auto, Frame, Outcome, PointingPrior, Policy, Priors, ScalePrior, Solution, SolveRequest,
    Source,
};
use fayde_wcs::{Header, Wcs};
use serde_json::Value;

use crate::{index_paths, num, options, seconds};

struct Run {
    solved: bool,
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
            let agree = reference.map(|r| max_disagreement(&sol, r, w, h));
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
                    solved: true,
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
                solved: false,
                agree_arcsec: None,
                seconds,
            },
            format!("abstained: {reason} ({hypotheses_tested} hyp)"),
        ),
    })
}

/// Largest sky separation between the two WCS over a 5×5 grid, arcsec.
fn max_disagreement(sol: &Solution, reference: &Wcs, w: f64, h: f64) -> f64 {
    let mut worst = 0.0f64;
    for i in 0..5 {
        for j in 0..5 {
            let p = PixelPos::new(-0.5 + w * i as f64 / 4.0, -0.5 + h * j as f64 / 4.0);
            worst = worst.max(to_arcsec(
                sol.wcs
                    .pixel_to_sky(p)
                    .separation(reference.pixel_to_sky(p)),
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
        &["ucac4", "index", "frames", "limit", "deadline"],
        &[],
    )?;
    let deadline = Some(seconds(&o, "deadline")?.unwrap_or(std::time::Duration::from_secs(60)));
    let cat = Ucac4Catalog::open(o.get("ucac4").ok_or("--ucac4 is required")?)
        .map_err(|e| e.to_string())?;
    let paths = index_paths(&o);
    if paths.is_empty() {
        return Err("--index is required".into());
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
        // The same header reading as `fayde solve --header`: pointing (CRVAL with CTYPE, or
        // RA/DEC), scale (CD, PC+CDELT, CDELT, optics; placeholders ignored), mid-exposure
        // epoch (JD/MJD/DATE-OBS + exposure/2), reference-frame checks.
        let hp = hdr
            .as_object()
            .map(fayde_api::header_from_json)
            .transpose()?
            .map(|h| header_priors(&h))
            .unwrap_or_default();
        let epoch = hp.epoch.map(|(e, _)| e);
        let scale = hp.scale_arcsec.map(|(s, _)| s);
        let pointing = hp.pointing.map(|(p, _)| p);
        let base = SolveRequest {
            frame: Frame {
                width: w as u32,
                height: h as u32,
                trail: None,
            },
            sources,
            priors: Priors::default(),
            policy: Policy {
                epoch,
                deadline,
                ..Policy::default()
            },
        };
        let mut hinted = base.clone();
        // Each prior stands on its own: a placeholder scale must not discard a usable pointing.
        if let Some(c) = pointing {
            hinted.priors.pointing = Some(PointingPrior {
                centre: c,
                sigma: 0.5f64.to_radians(),
                max_radius: 2f64.to_radians(),
            });
        }
        if let Some(s) = scale {
            hinted.priors.scale = Some(ScalePrior {
                arcsec_per_px: s,
                sigma: 0.02 * s,
            });
        }
        n += 1;
        ref_solved += usize::from(reference.is_some());
        let (a, ta) = run_one(&hinted, &cat, &indexes, reference.as_ref(), w, h)?;
        let (b, tb) = run_one(&base, &cat, &indexes, reference.as_ref(), w, h)?;
        times.0.push(a.seconds);
        times.1.push(b.seconds);
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
            *solved += usize::from(r.solved);
            match (r.solved, reference.is_some(), r.agree_arcsec) {
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
            "{:<42} {:>3} src  solve-field {:<6}\n    header priors {:6.2}s  {ta}\n    blind         {:6.2}s  {tb}",
            fr["name"].as_str().unwrap_or("?"),
            base.sources.len(),
            if reference.is_some() { "solved" } else { "FAILED" },
            a.seconds,
            b.seconds,
        );
    }
    agreements.sort_by(f64::total_cmp);
    let med = |v: &mut Vec<f64>| {
        v.sort_by(f64::total_cmp);
        v.get(v.len() / 2).copied()
    };
    for (k, v) in [
        ("frames", n),
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
