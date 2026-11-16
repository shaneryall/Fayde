//! Regression tests for the code-review findings fixed together: input validation, hard limits
//! on returned WCSs, references for off-centre CRPIX, sparse catalogues and truncated sources.

use fayde_core::Portable;
mod common;

use common::field;
use fayde_catalog::memory::MemoryCatalog;
use fayde_core::{to_arcsec, PixelPos};
use fayde_solve::{refine_wcs, solve, solve_auto, Outcome, RollPrior, SolveError, Trail};
use fayde_wcs::{Parity, TanWcs, Wcs};

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

#[test]
fn non_finite_priors_are_input_errors() {
    let (cat, mut f) = field((150.0, 2.0), 37.0, Parity::Negative, 41);
    f.request.priors.roll = Some(RollPrior {
        angle: f64::NAN,
        sigma: 0.1,
        max_dev: 0.1,
    });
    assert!(matches!(
        solve_auto(&f.request, &cat, &[]),
        Err(SolveError::InvalidInput(_))
    ));
    assert!(matches!(
        refine_wcs(&f.request, &cat, &f.truth),
        Err(SolveError::InvalidInput(_))
    ));
}

#[test]
fn refit_respects_hard_limits() {
    let (cat, mut f) = field((150.0, 2.0), 37.0, Parity::Negative, 42);
    f.request.priors.parity = Some(Parity::Positive);
    let out = refine_wcs(&f.request, &cat, &f.truth).unwrap();
    assert!(matches!(out, Outcome::Abstained { .. }), "{out:?}");
    f.request.priors.parity = Some(Parity::Negative);
    assert!(matches!(
        refine_wcs(&f.request, &cat, &f.truth).unwrap(),
        Outcome::Solved(_)
    ));
}

#[test]
fn references_follow_the_frame_not_crval() {
    // The same linear WCS re-expressed with its reference pixel far outside the frame.
    let (cat, f) = field((150.0, 2.0), 37.0, Parity::Negative, 43);
    let t = f.truth.tan;
    let crpix = PixelPos::new(-500.0, -500.0);
    let off: Wcs = TanWcs::new(t.pixel_to_sky(crpix), crpix, t.cd())
        .unwrap()
        .into();
    let Outcome::Solved(sol) = refine_wcs(&f.request, &cat, &off).unwrap() else {
        panic!("refit around an off-centre CRPIX found no references");
    };
    let err = grid_rms(&sol.wcs, &f.truth);
    assert!(err < 1.0, "{err}″");
}

#[test]
fn a_catalogue_sparser_than_the_seed_count_does_not_panic() {
    let (cat, f) = field((150.0, 2.0), 37.0, Parity::Negative, 44);
    let centre = f.request.priors.pointing.unwrap().centre;
    let mut stars = cat.stars().to_vec();
    stars.sort_by(|a, b| {
        a.pos
            .separation(centre)
            .total_cmp(&b.pos.separation(centre))
    });
    stars.truncate(6);
    let sparse = MemoryCatalog::new(stars);
    assert!(f.request.policy.seed_sources > 6);
    // Any outcome is fine; a panic is not.
    let _ = solve(&f.request, &sparse).unwrap();
}

#[test]
fn truncated_matches_do_not_buy_distortion_terms() {
    // All but six real sources truncated: they are evidence, not fit points, so the fit has
    // six points and must stay linear.
    let (cat, mut f) = field((150.0, 2.0), 37.0, Parity::Negative, 45);
    let mut kept = 0;
    for (k, s) in f.request.sources.iter_mut().enumerate() {
        if f.source_star[k].is_some() {
            if kept < 6 {
                kept += 1;
            } else {
                s.truncated = true;
            }
        }
    }
    f.request.frame.trail = Some(Trail {
        length_px: 4.0,
        angle: 0.3,
    });
    if let Outcome::Solved(sol) = solve_auto(&f.request, &cat, &[]).unwrap() {
        assert_eq!(sol.wcs.sip_order(), 1);
    }
}

#[test]
fn nothing_to_search_is_an_abstention_not_an_error() {
    // No pointing and no index (e.g. the scale prior was dropped for unknown binning): a valid
    // request with nothing to search answers with a reason.
    let (cat, mut f) = field((150.0, 2.0), 37.0, Parity::Negative, 46);
    f.request.priors.pointing = None;
    match solve_auto(&f.request, &cat, &[]).unwrap() {
        Outcome::Abstained { reason, .. } => assert!(reason.contains("nothing to search")),
        o => panic!("{o:?}"),
    }
}

#[test]
fn a_huge_deadline_is_no_deadline() {
    let mut req = fayde_solve::SolveRequest {
        frame: fayde_solve::Frame {
            width: 100,
            height: 100,
            trail: None,
        },
        sources: vec![fayde_solve::Source::at(10.0, 10.0)],
        priors: fayde_solve::Priors {
            pointing: Some(fayde_solve::PointingPrior {
                centre: fayde_core::SkyPos::from_degrees(10.0, 20.0).unwrap(),
                sigma: 0.01,
                max_radius: 0.02,
            }),
            scale: Some(fayde_solve::ScalePrior {
                arcsec_per_px: 2.0,
                sigma: 0.002,
            }),
            ..fayde_solve::Priors::default()
        },
        policy: fayde_solve::Policy::default(),
    };
    req.policy.deadline = Some(std::time::Duration::MAX);
    let cat = MemoryCatalog::new(vec![]);
    // Must not panic adding the deadline to the start time.
    let _ = fayde_solve::solve_auto(&req, &cat, &[]);
}

/// The abstention says which budget stopped the search and in which
/// stage; each stage records how it ended; a refit honours the deadline.
#[test]
fn budget_stops_name_the_stage_and_refit_honours_the_deadline() {
    let (cat, f) = field((150.0, 2.0), 37.0, Parity::Negative, 44);
    let mut req = f.request.clone();
    // Noise sources only: nothing can verify, so the search runs until a budget stops it.
    let mut seed = 7u64;
    for s in &mut req.sources {
        seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        s.x = (seed >> 11) as f64 / (1u64 << 53) as f64 * f64::from(req.frame.width);
        seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        s.y = (seed >> 11) as f64 / (1u64 << 53) as f64 * f64::from(req.frame.height);
    }
    req.policy.max_hypotheses = 50;
    let Outcome::Abstained { reason, search, .. } = solve_auto(&req, &cat, &[]).unwrap() else {
        panic!("noise must not solve");
    };
    assert!(
        reason.starts_with("hypothesis limit (50) reached during pair seeds"),
        "{reason}"
    );
    assert_eq!(search.last().unwrap().outcome, "hypothesis limit");

    req.policy.max_hypotheses = 5_000_000;
    req.policy.deadline = Some(std::time::Duration::ZERO);
    let Outcome::Abstained { reason, search, .. } = solve_auto(&req, &cat, &[]).unwrap() else {
        panic!("a zero deadline must not solve");
    };
    assert!(
        reason.starts_with("deadline reached during pair seeds"),
        "{reason}"
    );
    assert_eq!(search.last().unwrap().outcome, "deadline");

    let truth = f.truth.clone();
    let mut good = f.request.clone();
    good.policy.deadline = Some(std::time::Duration::ZERO);
    let Outcome::Abstained { reason, .. } = refine_wcs(&good, &cat, &truth).unwrap() else {
        panic!("a refit past its deadline must abstain");
    };
    assert!(
        reason.starts_with("deadline reached during refit"),
        "{reason}"
    );
    good.policy.deadline = None;
    assert!(matches!(
        refine_wcs(&good, &cat, &truth).unwrap(),
        Outcome::Solved(_)
    ));
}
