mod common;

use common::{field, field_full, field_with_noise, Field};
use fayde_catalog::memory::MemoryCatalog;
use fayde_core::{to_arcsec, PixelPos};
use fayde_solve::{solve, Outcome, Solution, SolveRequest};
use fayde_wcs::Parity;

fn solved(req: &SolveRequest, cat: &MemoryCatalog) -> Box<Solution> {
    let out = solve(req, cat).unwrap();
    match out {
        Outcome::Solved(s) => s,
        _ => panic!("expected a solution: {out:?}"),
    }
}

#[test]
fn tan_with_underreported_sigma_keeps_order_one() {
    let (cat, mut f) = field((150.0, 2.0), 37.0, Parity::Negative, 1);
    for s in &mut f.request.sources {
        s.sigma_px = Some(0.1);
    }
    let sol = solved(&f.request, &cat);
    eprintln!(
        "underreported sigma: order {}, scale {:.4}",
        sol.wcs.sip_order(),
        sol.error_scale
    );
    assert_eq!(sol.wcs.sip_order(), 1);
    assert!(sol.error_scale > 2.0);
}

#[test]
fn barrel_with_default_sigma_selects_sip() {
    let (cat, mut f) = field_full((150.0, 2.0), 37.0, Parity::Negative, 8, 0.2, 2e-8);
    for s in &mut f.request.sources {
        s.sigma_px = None;
    }
    let sol = solved(&f.request, &cat);
    eprintln!(
        "default sigma: order {}, scale {:.4}, RMS {:.4}",
        sol.wcs.sip_order(),
        sol.error_scale,
        sol.rms_px
    );
    assert!(sol.wcs.sip_order() >= 2);
    assert!(sol.error_scale < 0.4);
}

#[test]
fn pair_prior_rejections_do_not_count_as_hypotheses() {
    let (cat, mut f) = field((150.0, 2.0), 37.0, Parity::Negative, 1);
    f.request.priors.roll = Some(fayde_solve::RollPrior {
        angle: 1.23456789,
        sigma: 0.1,
        max_dev: 0.0,
    });
    let out = solve(&f.request, &cat).unwrap();
    let Outcome::Abstained {
        hypotheses_tested, ..
    } = out
    else {
        panic!("{out:?}")
    };
    eprintln!("pair prior rejects counted: {hypotheses_tested}");
    assert_eq!(hypotheses_tested, 0);
}

#[test]
fn pair_search_starts_at_scale_prior_mean() {
    let (cat, f) = field((150.0, 2.0), 37.0, Parity::Negative, 1);
    let sol = solved(&f.request, &cat);
    eprintln!("scale-prior ordering: {} hypotheses", sol.hypotheses_tested);
    const LOW_EDGE_HYPOTHESES: u64 = 78;
    assert!(2 * sol.hypotheses_tested < LOW_EDGE_HYPOTHESES);
}

#[test]
fn low_probability_sources_preserve_log_odds() {
    let (cat, mut f) = field((150.0, 2.0), 37.0, Parity::Negative, 1);
    let baseline = solved(&f.request, &cat);
    let mut rng = common::Rng(123);
    for _ in 0..120 {
        f.request.sources.push(fayde_solve::Source {
            x: rng.uniform() * 1023.0,
            y: rng.uniform() * 1023.0,
            flux: Some(0.01),
            sigma_px: Some(0.3),
            p_star: Some(0.05),
            truncated: false,
        });
    }
    let augmented = solved(&f.request, &cat);
    let loss = baseline.log_odds - augmented.log_odds;
    eprintln!(
        "low-p sources: log odds {:.4} -> {:.4}, loss {loss:.4}",
        baseline.log_odds, augmented.log_odds
    );
    assert!(loss < 15.0, "log-odds loss {loss}");
}

#[test]
fn a_hint_scale_outside_its_limits_gives_no_bent_solution() {
    for seed in [11u64, 12, 13, 14, 15, 16] {
        let (cat, mut f) = field_with_noise((150.0, 2.0), 37.0, Parity::Negative, seed, 1.0);
        let s = 3.0 * 1.03;
        f.request.priors.scale = Some(fayde_solve::ScalePrior {
            arcsec_per_px: s,
            sigma: s * 0.02 / 3.0,
        });
        if let Outcome::Solved(sol) = solve(&f.request, &cat).unwrap() {
            let err = worst_error(&sol, &f);
            assert!(
                err < 3.0,
                "seed {seed}: {err}″ off, order {}",
                sol.wcs.sip_order()
            );
        }
    }
}

fn worst_error(sol: &Solution, f: &Field) -> f64 {
    [
        (0.0, 0.0),
        (1023.0, 0.0),
        (0.0, 1023.0),
        (1023.0, 1023.0),
        (511.5, 511.5),
    ]
    .iter()
    .map(|&(x, y)| {
        let p = PixelPos::new(x, y);
        to_arcsec(sol.wcs.pixel_to_sky(p).separation(f.truth.pixel_to_sky(p)))
    })
    .fold(0.0, f64::max)
}
