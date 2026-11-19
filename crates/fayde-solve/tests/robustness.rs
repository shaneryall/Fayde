use fayde_core::Portable;
mod common;

use common::{field, field_full, Field, Rng};
use fayde_core::{to_arcsec, PixelPos};
use fayde_solve::{solve, Outcome, Solution, Source};
use fayde_wcs::Parity;

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

fn sparse(f: &mut Field, n_real: usize) {
    let src = &f.request.sources;
    let mut real: Vec<usize> = (0..src.len())
        .filter(|&i| f.source_star[i].is_some())
        .collect();
    real.sort_by(|&a, &b| src[b].flux.unwrap().total_cmp(&src[a].flux.unwrap()));
    real.truncate(n_real);
    let fake = (0..src.len())
        .filter(|&i| f.source_star[i].is_none())
        .take((n_real * 3).div_ceil(10));
    let keep: Vec<usize> = real.into_iter().chain(fake).collect();
    let sources: Vec<Source> = keep.iter().map(|&i| f.request.sources[i]).collect();
    f.source_star = keep.iter().map(|&i| f.source_star[i]).collect();
    f.request.sources = sources;
}

#[test]
fn sparse_frames_solve_without_an_interpolating_polynomial() {
    for n_real in [6usize, 8, 10, 14] {
        for seed in 1100..1106u64 {
            let (cat, mut f) = field(
                (20.0 + 3.0 * seed as f64 % 300.0, -40.0 + seed as f64 % 80.0),
                17.0 * seed as f64,
                if seed % 2 == 0 {
                    Parity::Negative
                } else {
                    Parity::Positive
                },
                seed,
            );
            sparse(&mut f, n_real);
            let out = solve(&f.request, &cat).unwrap();
            let Outcome::Solved(sol) = out else {
                panic!("{n_real} stars, seed {seed}: {out:?}")
            };
            assert_eq!(sol.wcs.sip_order(), 1, "{n_real} stars, seed {seed}");
            let err = worst_error(&sol, &f);
            assert!(err < 15.0, "{n_real} stars, seed {seed}: {err}″");
        }
    }
}

#[test]
fn strong_distortion_with_tight_centroids_solves() {
    for seed in [300u64, 301, 302, 303] {
        let (cat, f) = field_full((150.0, 2.0), 37.0, Parity::Negative, seed, 0.1, 1e-7);
        let out = solve(&f.request, &cat).unwrap();
        let Outcome::Solved(sol) = out else {
            panic!("seed {seed}: {out:?}")
        };
        assert!(
            sol.wcs.sip_order() >= 3,
            "seed {seed}: order {}",
            sol.wcs.sip_order()
        );
        let err = worst_error(&sol, &f);
        assert!(err < 2.0, "seed {seed}: {err}″");
    }
}

#[test]
fn faint_false_sources_do_not_block_a_solve() {
    let (cat, mut f) = field((150.0, 2.0), 37.0, Parity::Negative, 1);
    let clean = match solve(&f.request, &cat).unwrap() {
        Outcome::Solved(s) => s,
        out => panic!("{out:?}"),
    };
    let faintest = f
        .request
        .sources
        .iter()
        .filter_map(|s| s.flux)
        .fold(f64::INFINITY, f64::min);
    let mut rng = Rng(77);
    for _ in 0..400 {
        f.request.sources.push(Source {
            x: rng.uniform() * 1023.0,
            y: rng.uniform() * 1023.0,
            flux: Some(faintest * (0.2 + 0.7 * rng.uniform())),
            sigma_px: Some(0.3),
            p_star: None,
            truncated: false,
        });
        f.source_star.push(None);
    }
    let out = solve(&f.request, &cat).unwrap();
    let Outcome::Solved(sol) = out else {
        panic!("{out:?}")
    };
    assert!(worst_error(&sol, &f) < 1.0);
    assert!(
        sol.log_odds > 0.8 * clean.log_odds,
        "{} vs {}",
        sol.log_odds,
        clean.log_odds
    );
}

#[test]
fn short_random_lists_are_never_solved() {
    for seed in 900..940u64 {
        let (cat, mut f) = field((150.0, 2.0), 37.0, Parity::Negative, seed);
        let mut rng = Rng(seed * 3 + 11);
        f.request.sources.truncate(10);
        for s in &mut f.request.sources {
            s.x = rng.uniform() * 1023.0;
            s.y = rng.uniform() * 1023.0;
        }
        match solve(&f.request, &cat).unwrap() {
            Outcome::Abstained { best_log_odds, .. } => {
                assert!(
                    best_log_odds.is_none_or(|b| b < 1e9f64.pln()),
                    "seed {seed}: random list scored {best_log_odds:?}"
                );
            }
            Outcome::Solved(s) => {
                panic!("seed {seed}: random list solved, log-odds {}", s.log_odds)
            }
        }
    }
}

#[test]
fn distortion_tune_respects_last_hypothesis() {
    let (cat, mut f) = field_full((150., 2.), 37., Parity::Negative, 301, 0.1, 1e-7);
    for cap in [10, 11] {
        f.request.policy.max_hypotheses = cap;
        let tested = match solve(&f.request, &cat).unwrap() {
            Outcome::Solved(s) => s.hypotheses_tested,
            Outcome::Abstained {
                hypotheses_tested, ..
            } => hypotheses_tested,
        };
        assert!(tested <= cap, "{tested} hypotheses with cap {cap}");
    }
}

#[test]
fn high_order_refits_emit_reusable_wcs() {
    use fayde_catalog::{memory::MemoryCatalog, Star};
    use fayde_solve::{refine_wcs, Frame, Policy, SolveRequest};
    use fayde_wcs::{Header, Wcs};

    for (order, coeff) in [(8, 2e-22), (9, 4e-25)] {
        let cards = [
            "CTYPE1  = 'RA---TAN-SIP'",
            "CTYPE2  = 'DEC--TAN-SIP'",
            "CRVAL1  = 42.",
            "CRVAL2  = 15.",
            "CRPIX1  = 512.5",
            "CRPIX2  = 512.5",
            "CDELT1  = -0.0004",
            "CDELT2  = 0.0004",
        ];
        let terms = [
            format!("A_ORDER = {order}"),
            format!("B_ORDER = {order}"),
            format!("A_{order}_0   = {coeff:e}"),
            format!("B_0_{order}   = {:e}", -coeff),
        ];
        let truth = Wcs::from_header(&Header::parse_cards(
            cards.into_iter().chain(terms.iter().map(String::as_str)),
        ))
        .unwrap();
        let mut req = SolveRequest {
            frame: Frame::new(1024, 1024),
            sources: vec![],
            policy: Policy {
                max_sip_order: order,
                ..Default::default()
            },
            priors: Default::default(),
        };
        let mut stars = Vec::new();
        for i in 0..28 {
            for j in 0..28 {
                let (x, y) = (7. + 37. * f64::from(i), 7. + 37. * f64::from(j));
                let mut src = Source::at(x, y);
                src.sigma_px = Some(0.3);
                req.sources.push(src);
                stars.push(Star {
                    id: stars.len() as u64 + 1,
                    pos: truth.pixel_to_sky(PixelPos::new(x, y)),
                    mag: Some(10.),
                    pos_err_mas: [1.; 2],
                    pm_mas_yr: [None; 2],
                });
            }
        }
        let Outcome::Solved(sol) = refine_wcs(&req, &MemoryCatalog::new(stars), &truth).unwrap()
        else {
            panic!("order {order} refit abstained");
        };
        assert_eq!(sol.wcs.sip_order(), order);
        let restored = Wcs::from_header(&Header::parse_cards(
            sol.wcs.to_cards().iter().map(String::as_str),
        ))
        .unwrap();
        for x in [0., 255., 511.5, 1023.] {
            for y in [0., 511.5, 1023.] {
                let back = restored
                    .sky_to_pixel(truth.pixel_to_sky(PixelPos::new(x, y)))
                    .unwrap();
                assert!(
                    (back.x - x).abs().max((back.y - y).abs()) < 0.02,
                    "{back:?}"
                );
            }
        }
    }
}
