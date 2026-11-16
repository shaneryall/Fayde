//! End-to-end hinted solves on synthetic skies (no external data needed).
//!
//! Each field: stars drawn in a cap with a realistic magnitude distribution, projected through a
//! known TAN WCS; 20% of real stars dropped, 30% extra false sources added, 0.3 px centroid
//! noise, pointing prior offset from the truth.

use fayde_core::Portable;
mod common;

use common::{field, field_full, field_with_noise, Field, Rng};
use fayde_catalog::memory::MemoryCatalog;
use fayde_core::{arcsec, to_arcsec, PixelPos};
use fayde_solve::{solve, Outcome, RollPrior};
use fayde_wcs::Parity;

fn check_solves(name: &str, cat: &MemoryCatalog, f: &Field) {
    check_solves_within(name, cat, f, 1.0);
}

/// Errors at the frame corners and centre must lie within `max_err` arcsec and within 5× the
/// solver's own predicted radial 1σ (from the fit covariance).
fn check_solves_within(name: &str, cat: &MemoryCatalog, f: &Field, max_err: f64) {
    let out = solve(&f.request, cat).unwrap();
    let Outcome::Solved(sol) = out else {
        panic!("{name}: {out:?}");
    };
    let mut worst = 0.0f64;
    let mut worst_z = 0.0f64;
    for (x, y) in [
        (0.0, 0.0),
        (1023.0, 0.0),
        (0.0, 1023.0),
        (1023.0, 1023.0),
        (511.5, 511.5),
    ] {
        let p = PixelPos::new(x, y);
        let err = to_arcsec(sol.wcs.pixel_to_sky(p).separation(f.truth.pixel_to_sky(p)));
        worst = worst.max(err);
        worst_z = worst_z.max(err / sol.covariance.radial_sigma_arcsec(p).max(1e-9));
    }
    let confident: Vec<_> = sol.matches.iter().filter(|m| m.inlier_prob > 0.5).collect();
    let correct = confident
        .iter()
        .filter(|m| f.source_star[m.source] == Some(m.star_id))
        .count();
    eprintln!(
        "{name}: worst sky error {worst:.3}″ ({worst_z:.1}σ), SIP order {}, {}/{} confident matches correct, log-odds {:.1}, \
         NFA {:.1e}, {} hypotheses, rms {:.2} px, error scale {:.2}, {:?}",
        sol.wcs.sip_order(),
        correct,
        confident.len(),
        sol.log_odds,
        sol.nfa,
        sol.hypotheses_tested,
        sol.rms_px,
        sol.error_scale,
        sol.elapsed
    );
    assert!(worst < max_err, "{name}: {worst}″");
    if max_err > 1e-3 {
        assert!(
            worst_z < 5.0,
            "{name}: error {worst_z:.1}σ of the reported covariance"
        );
    }
    assert!(
        confident.len() >= 20,
        "{name}: only {} matches",
        confident.len()
    );
    assert!(
        correct as f64 >= 0.98 * confident.len() as f64,
        "{name}: wrong matches"
    );
    assert!(sol.nfa < 1e-9, "{name}: NFA {}", sol.nfa);
    assert_eq!(sol.wcs.parity(), f.truth.parity());
}

#[test]
fn solves_equatorial_field() {
    let (cat, f) = field((150.0, 2.0), 37.0, Parity::Negative, 1);
    check_solves("equator", &cat, &f);
}

#[test]
fn solves_near_pole() {
    let (cat, f) = field((12.0, 85.0), -120.0, Parity::Negative, 2);
    check_solves("dec 85°", &cat, &f);
}

#[test]
fn solves_across_ra_zero() {
    let (cat, f) = field((0.1, -30.0), 5.0, Parity::Negative, 3);
    check_solves("RA 0.1°", &cat, &f);
}

#[test]
fn solves_mirrored_image() {
    let (cat, f) = field((210.0, 45.0), 200.0, Parity::Positive, 4);
    check_solves("flipped parity", &cat, &f);
}

#[test]
fn roll_prior_restricts_and_still_solves() {
    let (cat, mut f) = field((150.0, 2.0), 37.0, Parity::Negative, 5);
    f.request.priors.roll = Some(RollPrior {
        angle: 37f64.to_radians(),
        sigma: 0.5f64.to_radians(),
        max_dev: 2f64.to_radians(),
    });
    f.request.priors.parity = Some(Parity::Negative);
    check_solves("roll + parity priors", &cat, &f);

    // A wrong roll prior must not produce a (wrong) solution.
    f.request.priors.roll = Some(RollPrior {
        angle: 137f64.to_radians(),
        sigma: 0.5f64.to_radians(),
        max_dev: 2f64.to_radians(),
    });
    assert!(matches!(
        solve(&f.request, &cat).unwrap(),
        Outcome::Abstained { .. }
    ));
}

#[test]
fn random_sources_are_never_solved() {
    let (cat, mut f) = field((150.0, 2.0), 37.0, Parity::Negative, 6);
    let mut rng = Rng(99);
    for s in &mut f.request.sources {
        s.x = rng.uniform() * 1023.0;
        s.y = rng.uniform() * 1023.0;
    }
    match solve(&f.request, &cat).unwrap() {
        Outcome::Abstained {
            hypotheses_tested,
            best_log_odds,
            ..
        } => {
            eprintln!("random: abstained after {hypotheses_tested} hypotheses, best log-odds {best_log_odds:?}");
            assert!(hypotheses_tested > 1000);
        }
        Outcome::Solved(s) => panic!("random sources solved: log-odds {}", s.log_odds),
    }
}

#[test]
fn rejects_blind_request() {
    let (cat, mut f) = field((150.0, 2.0), 0.0, Parity::Negative, 7);
    f.request.priors.pointing = None;
    assert!(matches!(
        solve(&f.request, &cat),
        Err(fayde_solve::SolveError::Unsupported(_))
    ));
    let _ = arcsec(1.0);
}

#[test]
fn noiseless_field_is_recovered_exactly() {
    // Without centroid noise the fit must reproduce the true WCS to numerical precision: any
    // bias in seeding, matching or fitting would show here.
    let (cat, f) = field_with_noise((210.0, 45.0), 200.0, Parity::Positive, 4, 0.0);
    check_solves_within("noiseless", &cat, &f, 1e-4);
}

/// RMS sky error over a 9×9 grid covering the frame, arcsec.
fn grid_rms(wcs: &fayde_wcs::Wcs, truth: &fayde_wcs::Wcs) -> f64 {
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
fn distorted_field_gets_sip() {
    // ~1.4 px barrel distortion at the corners. TAN alone leaves a systematic error that shows in
    // the error scale; SIP (chosen by BIC) removes it. A per-frame cubic is weakly constrained at
    // the extreme corners, which the reported covariance says honestly (design doc §5.5: fit
    // distortion per instrument across frames, then only small per-frame corrections).
    let (cat, mut f) = field_full((150.0, 2.0), 37.0, Parity::Negative, 8, 0.3, 2e-8);
    check_solves_within("barrel, SIP", &cat, &f, 5.0);
    let Outcome::Solved(sol) = solve(&f.request, &cat).unwrap() else {
        panic!()
    };
    assert!(
        sol.wcs.sip_order() >= 3,
        "BIC kept order {}",
        sol.wcs.sip_order()
    );

    f.request.policy.max_sip_order = 1;
    let Outcome::Solved(tan) = solve(&f.request, &cat).unwrap() else {
        panic!()
    };
    let (tan_rms, sip_rms) = (grid_rms(&tan.wcs, &f.truth), grid_rms(&sol.wcs, &f.truth));
    eprintln!(
        "grid RMS: TAN {tan_rms:.2}″ (error scale {:.2}), SIP {sip_rms:.2}″ (error scale {:.2})",
        tan.error_scale, sol.error_scale
    );
    assert!(tan_rms > 2.0 * sip_rms);
    assert!(
        tan.error_scale > 1.5,
        "TAN misfit must show in the error scale"
    );
    assert!(sol.error_scale < 1.2);
}

#[test]
fn one_missing_flux_keeps_flux_ranking() {
    // Previously a single source without flux switched the whole list to input order (which is
    // shuffled here); the solve must still start from the brightest sources.
    let (cat, mut f) = field((150.0, 2.0), 37.0, Parity::Negative, 9);
    f.request.sources[0].flux = None;
    let Outcome::Solved(sol) = solve(&f.request, &cat).unwrap() else {
        panic!()
    };
    assert!(
        sol.hypotheses_tested < 5_000,
        "{} hypotheses",
        sol.hypotheses_tested
    );
}

#[test]
fn duplicate_fragments_are_demoted() {
    // Every one of the 15 brightest sources gets a fragment 1 px away with 80% of its flux.
    let (cat, mut f) = field((150.0, 2.0), 37.0, Parity::Negative, 10);
    let mut by_flux: Vec<usize> = (0..f.request.sources.len()).collect();
    by_flux.sort_by(|&a, &b| {
        f.request.sources[b]
            .flux
            .unwrap()
            .total_cmp(&f.request.sources[a].flux.unwrap())
    });
    for &i in by_flux.iter().take(15) {
        let mut d = f.request.sources[i];
        d.x += 1.0;
        d.flux = d.flux.map(|v| v * 0.8);
        f.request.sources.push(d);
        f.source_star.push(None);
    }
    let Outcome::Solved(sol) = solve(&f.request, &cat).unwrap() else {
        panic!()
    };
    let n0 = f.request.sources.len() - 15;
    let dup_matches = sol
        .matches
        .iter()
        .filter(|m| m.source >= n0 && m.inlier_prob > 0.5)
        .count();
    eprintln!(
        "duplicates: {} hypotheses, {dup_matches} confident matches on fragments",
        sol.hypotheses_tested
    );
    assert_eq!(dup_matches, 0);
    check_solves_within("with fragments", &cat, &f, 1.0);
}
