mod common;

use common::{field, Field};
use fayde_catalog::{memory::MemoryCatalog, Catalog};
use fayde_core::{to_arcsec, PixelPos, SkyPos};
use fayde_index::{build, BuildParams, Index};
use fayde_solve::{solve_auto, Outcome, PointingPrior, Priors, ScalePrior, Solution, SolveMethod};
use fayde_wcs::Parity;

const CENTRE: (f64, f64) = (150.0, 2.0);

fn setup() -> (MemoryCatalog, Field, Vec<Index>) {
    let (cat, f) = field(CENTRE, 37.0, Parity::Negative, 1);
    let index = build(
        cat.stars(),
        cat.id(),
        &BuildParams {
            fov_deg: 1.3,
            min_edge_deg: 0.05,
            stars_per_fov: 40.0,
            region: Some([CENTRE.0, CENTRE.1, 1.6]),
            ..BuildParams::default()
        },
    )
    .unwrap();
    (cat, f, vec![index])
}

fn stages(sol: &Solution) -> Vec<String> {
    sol.search.iter().map(|s| s.stage.clone()).collect()
}

fn accurate(sol: &Solution, f: &Field) {
    let p = PixelPos::new(0.0, 0.0);
    let err = to_arcsec(sol.wcs.pixel_to_sky(p).separation(f.truth.pixel_to_sky(p)));
    assert!(err < 5.0 * sol.covariance.radial_sigma_arcsec(p), "{err}″");
}

fn solved(out: Outcome) -> Box<Solution> {
    match out {
        Outcome::Solved(s) => s,
        o => panic!("{o:?}"),
    }
}

#[test]
fn tight_priors_start_with_pair_seeds() {
    let (cat, f, idx) = setup();
    let sol = solved(solve_auto(&f.request, &cat, &idx).unwrap());
    eprintln!("tight: {:?}", stages(&sol));
    assert_eq!(sol.method, SolveMethod::PairSeeds);
    assert!(stages(&sol)[0].starts_with("pair seeds"));
    assert_eq!(sol.search.len(), 1);
    accurate(&sol, &f);
}

#[test]
fn pointing_with_loose_scale_searches_the_inner_ring_first() {
    let (cat, mut f, idx) = setup();
    f.request.priors.scale = Some(ScalePrior {
        arcsec_per_px: 3.2,
        sigma: 0.5,
    });
    let mut pt = f.request.priors.pointing.unwrap();
    pt.max_radius = 2f64.to_radians();
    f.request.priors.pointing = Some(pt);
    let sol = solved(solve_auto(&f.request, &cat, &idx).unwrap());
    let s = stages(&sol);
    eprintln!("loose scale: {s:?}");
    assert_eq!(sol.method, SolveMethod::PatternHash);
    assert!(s[0].ends_with("r<0.40°"), "{s:?}");
    assert_eq!(s.len(), 1, "solved in the inner ring");
    assert!(!sol.prior_violated);
    accurate(&sol, &f);
}

#[test]
fn pointing_without_scale_uses_rings_not_all_sky() {
    let (cat, mut f, idx) = setup();
    f.request.priors.scale = None;
    let sol = solved(solve_auto(&f.request, &cat, &idx).unwrap());
    let s = stages(&sol);
    eprintln!("no scale: {s:?}");
    assert!(s
        .iter()
        .all(|x| !x.contains("all sky") && !x.contains("outside")));
    accurate(&sol, &f);
}

#[test]
fn no_metadata_is_blind() {
    let (cat, mut f, idx) = setup();
    f.request.priors = Priors::default();
    let sol = solved(solve_auto(&f.request, &cat, &idx).unwrap());
    let s = stages(&sol);
    eprintln!("blind: {s:?}");
    assert!(s[0].ends_with("all sky"));
    accurate(&sol, &f);
}

#[test]
fn wrong_pointing_abstains_unless_outside_search_is_allowed() {
    let (cat, mut f, idx) = setup();
    f.request.priors.pointing = Some(PointingPrior {
        centre: SkyPos::from_degrees(CENTRE.0 + 1.2, CENTRE.1).unwrap(),
        sigma: 0.2f64.to_radians(),
        max_radius: 0.5f64.to_radians(),
    });
    match solve_auto(&f.request, &cat, &idx).unwrap() {
        Outcome::Abstained { search, .. } => {
            let s: Vec<_> = search.iter().map(|x| x.stage.as_str()).collect();
            eprintln!("wrong pointing: abstained after {s:?}");
            assert!(s.iter().all(|x| !x.contains("outside")));
        }
        Outcome::Solved(s) => panic!(
            "solved outside its prior without permission: {:?}",
            stages(&s)
        ),
    }

    f.request.policy.search_outside_pointing = true;
    let sol = solved(solve_auto(&f.request, &cat, &idx).unwrap());
    let s = stages(&sol);
    eprintln!("wrong pointing, outside allowed: {s:?}");
    assert!(sol.prior_violated);
    assert!(s.last().unwrap().contains("outside pointing limit"));
    accurate(&sol, &f);
}

#[test]
fn indexes_that_cannot_cover_the_pointing_are_skipped() {
    let (cat, mut f, idx) = setup();
    f.request.priors.scale = Some(ScalePrior {
        arcsec_per_px: 3.2,
        sigma: 0.5,
    });
    f.request.priors.pointing = Some(PointingPrior {
        centre: SkyPos::from_degrees(CENTRE.0 + 20.0, CENTRE.1).unwrap(),
        sigma: 0.2f64.to_radians(),
        max_radius: 0.5f64.to_radians(),
    });
    match solve_auto(&f.request, &cat, &idx).unwrap() {
        Outcome::Abstained {
            search,
            hypotheses_tested,
            reason,
            ..
        } => {
            assert!(search.is_empty(), "{search:?}");
            assert_eq!(reason, "no index covers the pointing prior");
            assert_eq!(hypotheses_tested, 0);
        }
        Outcome::Solved(_) => panic!("solved far outside the index region"),
    }
}

#[test]
fn a_small_pointing_sigma_searches_close_to_the_pointing_first() {
    let (cat, mut f, _) = setup();
    let mut pt = f.request.priors.pointing.unwrap();
    pt.max_radius = 2f64.to_radians();
    pt.sigma = 0.02f64.to_radians();
    f.request.priors.pointing = Some(pt);
    let sol = solved(solve_auto(&f.request, &cat, &[]).unwrap());
    let s = stages(&sol);
    eprintln!("small sigma: {s:?}");
    assert_eq!(s, ["pair seeds r<0.06°"]);
    assert_eq!(sol.method, SolveMethod::PairSeeds);
    assert!(!sol.prior_violated);
    accurate(&sol, &f);

    pt.sigma = 0.4f64.to_radians();
    f.request.priors.pointing = Some(pt);
    let wide = solved(solve_auto(&f.request, &cat, &[]).unwrap());
    assert_eq!(stages(&wide), ["pair seeds r<2.00°"]);
    accurate(&wide, &f);
    let p = PixelPos::new(0.0, 0.0);
    let apart = to_arcsec(sol.wcs.pixel_to_sky(p).separation(wide.wcs.pixel_to_sky(p)));
    assert!(apart < 0.5, "{apart}″ between the two solutions");
    pt.sigma = 0.0;
    f.request.priors.pointing = Some(pt);
    let flat = solved(solve_auto(&f.request, &cat, &[]).unwrap());
    assert_eq!(stages(&flat), ["pair seeds r<2.00°"]);
}

#[test]
fn a_pointing_in_the_tail_is_found_by_the_search_to_the_limit() {
    let (cat, mut f, _) = setup();
    f.request.priors.pointing = Some(PointingPrior {
        centre: SkyPos::from_degrees(CENTRE.0 + 1.6, CENTRE.1).unwrap(),
        sigma: 0.05f64.to_radians(),
        max_radius: 2f64.to_radians(),
    });
    let sol = solved(solve_auto(&f.request, &cat, &[]).unwrap());
    let s = stages(&sol);
    eprintln!("pointing in the tail: {s:?}");
    assert_eq!(s, ["pair seeds r<0.15°", "pair seeds r<2.00°"]);
    assert_eq!(sol.search[0].outcome, "exhausted");
    assert_eq!(sol.search[1].outcome, "solved");
    assert!(!sol.prior_violated);
    accurate(&sol, &f);
}
