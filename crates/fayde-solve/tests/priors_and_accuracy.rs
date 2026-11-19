use fayde_core::Portable;
mod common;

use common::{field, Field, Rng};
use fayde_catalog::memory::MemoryCatalog;
use fayde_core::{to_arcsec, PixelPos};
use fayde_solve::{
    refine_wcs, solve, solve_auto, verify_wcs, Outcome, RollPrior, ScalePrior, Solution,
    SolveMethod, Trail,
};
use fayde_wcs::{Parity, Poly, Sip, TanWcs, Wcs};

fn solved(out: Outcome) -> Box<Solution> {
    match out {
        Outcome::Solved(s) => s,
        o => panic!("{o:?}"),
    }
}

fn grid_rms(wcs: &Wcs, truth: &Wcs) -> f64 {
    let mut s = 0.0;
    for i in 0..9 {
        for j in 0..9 {
            let p = PixelPos::new(1023.0 * i as f64 / 8.0, 1023.0 * j as f64 / 8.0);
            s += to_arcsec(wcs.pixel_to_sky(p).separation(truth.pixel_to_sky(p))).ppowi(2);
        }
    }
    (s / 81.0).sqrt()
}

fn roll_of(w: &Wcs) -> f64 {
    let cd = w.cd();
    cd[0][1].patan2(cd[1][1])
}

#[test]
fn accuracy_is_reported_where_it_is_worst_and_leave_one_out_tracks_noise() {
    let (cat, f) = field((150.0, 2.0), 37.0, Parity::Negative, 31);
    let sol = solved(solve(&f.request, &cat).unwrap());
    let centre = sol.covariance.radial_sigma_arcsec(sol.wcs.crpix());
    let loo = sol.loo_rms_px.expect("enough matches");
    eprintln!(
        "centre σ {centre:.3}″, worst σ {:.3}″, rms {:.3} px, leave-one-out rms {loo:.3} px",
        sol.worst_sigma_arcsec, sol.rms_px
    );
    assert!(
        sol.worst_sigma_arcsec > 1.5 * centre,
        "corners are less certain than the centre"
    );
    assert!(
        loo >= sol.rms_px && loo < 1.3 * sol.rms_px,
        "{loo} vs {}",
        sol.rms_px
    );
    assert!((loo / (0.3 * 2f64.sqrt()) - 1.0).abs() < 0.25);
}

fn truncate_edges(f: &mut Field, trail: Trail, mark: bool) -> usize {
    let (ca, sa) = (trail.angle.pcos(), trail.angle.psin());
    let mut n = 0;
    for (k, s) in f.request.sources.iter_mut().enumerate() {
        let edge = s.x.min(s.y).min(1023.0 - s.x).min(1023.0 - s.y);
        if f.source_star[k].is_some() && edge < 0.5 * trail.length_px {
            let cut = 0.5 * trail.length_px - edge;
            let sign = if (s.x - 511.5) * ca + (s.y - 511.5) * sa > 0.0 {
                -1.0
            } else {
                1.0
            };
            s.x += sign * 0.5 * cut * ca;
            s.y += sign * 0.5 * cut * sa;
            s.truncated = mark;
            n += 1;
        }
    }
    n
}

#[test]
fn truncated_sources_are_kept_as_evidence_but_out_of_the_fit() {
    let trail = Trail {
        length_px: 90.0,
        angle: 0.6,
    };
    let (cat, mut plain) = field((150.0, 2.0), 37.0, Parity::Negative, 32);
    let mut flagged = Field {
        request: plain.request.clone(),
        truth: plain.truth.clone(),
        source_star: plain.source_star.clone(),
    };
    let n = truncate_edges(&mut plain, trail, false);
    truncate_edges(&mut flagged, trail, true);
    flagged.request.frame.trail = Some(trail);
    plain.request.frame.trail = Some(trail);
    let a = solved(solve(&plain.request, &cat).unwrap());
    let b = solved(solve(&flagged.request, &cat).unwrap());
    let (ea, eb) = (
        grid_rms(&a.wcs, &plain.truth),
        grid_rms(&b.wcs, &flagged.truth),
    );
    let kept = b
        .matches
        .iter()
        .filter(|m| flagged.request.sources[m.source].truncated && m.inlier_prob > 0.5)
        .count();
    eprintln!(
        "{n} truncated sources: grid rms unflagged {ea:.3}″, flagged {eb:.3}″; \
         {kept} truncated sources matched as evidence"
    );
    assert!(n >= 5, "the test needs edge sources");
    assert!(eb < ea, "flagging must help");
    assert!(eb < 0.5, "flagged solution is noise-limited: {eb}″");
    assert!(kept >= n / 2, "truncated sources still count as matches");
}

#[test]
fn verify_accepts_the_truth_and_rejects_a_rotated_wcs() {
    let (cat, f) = field((150.0, 2.0), 37.0, Parity::Negative, 33);
    let good = verify_wcs(&f.request, &cat, &f.truth).unwrap();
    let t = f.truth.tan;
    let (s, c) = 2f64.to_radians().psin_cos();
    let cd = t.cd();
    let rot = [
        [c * cd[0][0] - s * cd[1][0], c * cd[0][1] - s * cd[1][1]],
        [s * cd[0][0] + c * cd[1][0], s * cd[0][1] + c * cd[1][1]],
    ];
    let bad_wcs: Wcs = TanWcs::new(t.crval(), t.crpix(), rot).unwrap().into();
    let bad = verify_wcs(&f.request, &cat, &bad_wcs).unwrap();
    eprintln!(
        "verify: truth log-odds {:.1} ({} matches, rms {:.2} px); rotated 2°: {:.1}",
        good.log_odds,
        good.matches.len(),
        good.rms_px,
        bad.log_odds
    );
    assert!(good.accepted && good.matches.len() > 50);
    assert!(!bad.accepted);
}

#[test]
fn verify_ignores_search_priors() {
    let (cat, mut f) = field((150.0, 2.0), 37.0, Parity::Negative, 33);
    let expected = verify_wcs(&f.request, &cat, &f.truth).unwrap();
    f.request.priors.scale = Some(ScalePrior {
        arcsec_per_px: 90_000.0,
        sigma: 1.0,
    });
    f.request.priors.roll = Some(RollPrior {
        angle: 0.0,
        sigma: 1.0,
        max_dev: -1.0,
    });
    assert_eq!(verify_wcs(&f.request, &cat, &f.truth).unwrap(), expected);
    assert!(solve(&f.request, &cat).is_err());
}

#[test]
fn refit_repairs_a_perturbed_overfitted_wcs() {
    let (cat, f) = field((150.0, 2.0), 37.0, Parity::Negative, 34);
    let t = f.truth.tan;
    let (s, c) = 0.3f64.to_radians().psin_cos();
    let cd = t.cd();
    let rot = [
        [c * cd[0][0] - s * cd[1][0], c * cd[0][1] - s * cd[1][1]],
        [s * cd[0][0] + c * cd[1][0], s * cd[0][1] + c * cd[1][1]],
    ];
    let shifted = t.pixel_to_sky(PixelPos::new(t.crpix().x + 4.0, t.crpix().y - 2.0));
    let empty = Poly {
        order: 0,
        terms: vec![],
    };
    let stored = Wcs {
        tan: TanWcs::new(shifted, t.crpix(), rot).unwrap(),
        sip: Some(Sip {
            a: Poly {
                order: 2,
                terms: vec![(2, 0, 8e-6), (1, 1, -5e-6)],
            },
            b: Poly {
                order: 2,
                terms: vec![(0, 2, 6e-6)],
            },
            ap: empty.clone(),
            bp: empty,
        }),
    };
    let before = grid_rms(&stored, &f.truth);
    let sol = solved(refine_wcs(&f.request, &cat, &stored).unwrap());
    let after = grid_rms(&sol.wcs, &f.truth);
    eprintln!(
        "refit: grid rms {before:.2}″ → {after:.3}″, SIP order {}",
        sol.wcs.sip_order()
    );
    assert_eq!(sol.method, SolveMethod::Refit);
    assert_eq!(sol.hypotheses_tested, 1);
    assert!(before > 10.0 && after < 0.5);
    assert_eq!(sol.wcs.sip_order(), 1, "no distortion in the truth");
}

fn sparse(seed: u64, n: usize) -> (MemoryCatalog, Field) {
    let (cat, mut f) = field((150.0, 2.0), 37.0, Parity::Negative, seed);
    let mut real: Vec<usize> = (0..f.request.sources.len())
        .filter(|&k| f.source_star[k].is_some())
        .collect();
    real.sort_by(|&a, &b| {
        f.request.sources[b]
            .flux
            .unwrap()
            .total_cmp(&f.request.sources[a].flux.unwrap())
    });
    let mut keep: Vec<usize> = real[..n].to_vec();
    keep.extend(
        (0..f.request.sources.len())
            .filter(|&k| f.source_star[k].is_none())
            .take(3),
    );
    let mut rng = Rng(seed);
    keep.sort_by_key(|_| (rng.uniform() * 1e9) as u64);
    f.request.sources = keep.iter().map(|&k| f.request.sources[k]).collect();
    f.source_star = keep.iter().map(|&k| f.source_star[k]).collect();
    f.request.priors.scale = Some(ScalePrior {
        arcsec_per_px: 3.0,
        sigma: 0.003,
    });
    (cat, f)
}

#[test]
fn roll_and_parity_priors_on_sparse_frames() {
    for n in [4usize, 5, 6] {
        let (mut solved_without, mut solved_with, mut hyp_without, mut hyp_with) =
            (0, 0, 0u64, 0u64);
        for seed in 0..10u64 {
            let (cat, mut f) = sparse(100 + seed, n);
            let tested = |o: &Outcome| match o {
                Outcome::Solved(s) => s.hypotheses_tested,
                Outcome::Abstained {
                    hypotheses_tested, ..
                } => *hypotheses_tested,
            };
            let check = |o: &Outcome, f: &Field| match o {
                Outcome::Solved(s) => {
                    assert!(
                        grid_rms(&s.wcs, &f.truth) < 5.0,
                        "wrong solution with {n} stars"
                    );
                    1
                }
                Outcome::Abstained { .. } => 0,
            };
            let without = solve(&f.request, &cat).unwrap();
            f.request.priors.roll = Some(RollPrior {
                angle: roll_of(&f.truth),
                sigma: 0.1f64.to_radians(),
                max_dev: 1f64.to_radians(),
            });
            f.request.priors.parity = Some(Parity::Negative);
            let with = solve(&f.request, &cat).unwrap();
            solved_without += check(&without, &f);
            solved_with += check(&with, &f);
            hyp_without += tested(&without);
            hyp_with += tested(&with);
        }
        eprintln!(
            "{n} stars + 3 false, 10 fields: without roll/parity {solved_without}/10 solved \
             ({hyp_without} hypotheses), with {solved_with}/10 ({hyp_with} hypotheses)"
        );
        assert!(solved_with >= solved_without);
        assert!(hyp_with <= hyp_without);
        assert!(
            solved_with >= 8,
            "{n} stars: only {solved_with}/10 solved with roll/parity"
        );
    }
}

#[test]
fn a_wrong_hard_limit_costs_a_retry_not_the_solution() {
    let (cat, mut f) = field((150.0, 2.0), 37.0, Parity::Negative, 35);
    f.request.priors.roll = Some(RollPrior {
        angle: roll_of(&f.truth) + 10f64.to_radians(),
        sigma: 0.1f64.to_radians(),
        max_dev: 1f64.to_radians(),
    });
    f.request.priors.parity = Some(Parity::Negative);
    let strict = solve_auto(&f.request, &cat, &[]).unwrap();
    let Outcome::Abstained { reason, .. } = &strict else {
        panic!("{strict:?}");
    };
    assert!(reason.contains("outside the prior hard limits"), "{reason}");

    f.request.policy.relax_priors_on_failure = true;
    let sol = solved(solve_auto(&f.request, &cat, &[]).unwrap());
    let err = grid_rms(&sol.wcs, &f.truth);
    assert!(sol.elapsed.as_secs_f64() >= sol.search.iter().map(|s| s.seconds).sum::<f64>());
    eprintln!(
        "relaxed retry: grid rms {err:.3}″, violated {:?}, {} hypotheses over {} stages",
        sol.violated_priors,
        sol.hypotheses_tested,
        sol.search.len()
    );
    assert!(err < 1.0);
    assert!(sol.prior_violated);
    assert_eq!(sol.violated_priors, vec!["roll"]);
    assert!(sol
        .search
        .iter()
        .any(|s| s.stage.starts_with("relaxed priors: ")));
    assert!(
        (sol.nfa - sol.hypotheses_tested as f64 * (-sol.log_odds).pexp()).abs()
            <= 1e-9 * sol.nfa.max(1e-300)
    );

    f.request.priors.roll = Some(RollPrior {
        angle: roll_of(&f.truth),
        sigma: 0.1f64.to_radians(),
        max_dev: 1f64.to_radians(),
    });
    let sol = solved(solve_auto(&f.request, &cat, &[]).unwrap());
    assert!(!sol.prior_violated && sol.violated_priors.is_empty());
    assert!(!sol.search.iter().any(|s| s.stage.starts_with("relaxed")));
}

#[test]
fn instrument_priors_hold_sparse_fits_near_the_instrument() {
    let (mut loose_err, mut tight_err, mut loose_ds, mut tight_ds) =
        (vec![], vec![], vec![], vec![]);
    for seed in 0..20u64 {
        let (cat, mut f) =
            common::field_with_noise((150.0, 2.0), 37.0, Parity::Negative, 200 + seed, 1.5);
        let mut real: Vec<usize> = (0..f.request.sources.len())
            .filter(|&k| f.source_star[k].is_some())
            .collect();
        real.sort_by(|&a, &b| {
            f.request.sources[b]
                .flux
                .unwrap()
                .total_cmp(&f.request.sources[a].flux.unwrap())
        });
        let mut keep: Vec<usize> = real[..6].to_vec();
        keep.extend(
            (0..f.request.sources.len())
                .filter(|&k| f.source_star[k].is_none())
                .take(3),
        );
        f.request.sources = keep.iter().map(|&k| f.request.sources[k]).collect();
        f.source_star = keep.iter().map(|&k| f.source_star[k]).collect();
        f.request.priors.parity = Some(Parity::Negative);
        let truth_scale = f.truth.pixel_scale_arcsec();
        let mut run = |scale_sigma: f64, roll_sigma_deg: f64| {
            f.request.priors.scale = Some(ScalePrior {
                arcsec_per_px: truth_scale,
                sigma: scale_sigma * truth_scale,
            });
            f.request.priors.roll = Some(RollPrior {
                angle: roll_of(&f.truth),
                sigma: roll_sigma_deg.to_radians(),
                max_dev: (10.0 * roll_sigma_deg).to_radians(),
            });
            match solve(&f.request, &cat).unwrap() {
                Outcome::Solved(s) => Some((
                    grid_rms(&s.wcs, &f.truth),
                    (s.wcs.pixel_scale_arcsec() / truth_scale - 1.0).abs(),
                )),
                Outcome::Abstained { .. } => None,
            }
        };
        if let (Some(l), Some(t)) = (run(0.02, 3.0), run(0.001, 0.1)) {
            loose_err.push(l.0);
            tight_err.push(t.0);
            loose_ds.push(l.1);
            tight_ds.push(t.1);
        }
    }
    let mean = |v: &[f64]| v.iter().sum::<f64>() / v.len().max(1) as f64;
    let max = |v: &[f64]| v.iter().copied().fold(0.0, f64::max);
    eprintln!(
        "{} fields solved by both: grid rms loose {:.2}″ (max {:.2}″), tight {:.2}″ (max {:.2}″); \
         |scale error| loose {:.3}% (max {:.3}%), tight {:.3}% (max {:.3}%)",
        loose_err.len(),
        mean(&loose_err),
        max(&loose_err),
        mean(&tight_err),
        max(&tight_err),
        100.0 * mean(&loose_ds),
        100.0 * max(&loose_ds),
        100.0 * mean(&tight_ds),
        100.0 * max(&tight_ds),
    );
    assert!(loose_err.len() >= 15);
    assert!(mean(&tight_err) < 0.8 * mean(&loose_err));
    assert!(max(&tight_ds) < 0.004, "tight scale prior held to 0.4%");
}
