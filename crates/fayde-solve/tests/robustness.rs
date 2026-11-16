//! Regression fields for sparse frames, lens distortion and contaminated source lists.

use fayde_core::Portable;
mod common;

use common::{field, field_full, Field, Rng};
use fayde_core::{to_arcsec, PixelPos};
use fayde_solve::{solve, Outcome, Solution, Source};
use fayde_wcs::Parity;

/// Worst sky error over the frame corners and centre, arcsec.
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

/// Keep the `n_real` brightest real sources and 30% as many false ones.
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
    // With exactly 6 or 10 matches, order 2 or 3 interpolates the stars. The profiled BIC used
    // to choose it every time and return a WCS wrong by arcminutes at the corners.
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
    // 5% barrel distortion (27 px per axis at the corners) with 0.1 px centroids: no linear
    // model puts the outer stars inside a 0.5 px gate, so acceptance needs the tune-up.
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
    // 400 false sources fainter than every star, about five per real star: summed over the
    // whole list their penalties outweigh the matches.
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
    // The best prefix ends near the last star, so the junk costs little.
    assert!(
        sol.log_odds > 0.8 * clean.log_odds,
        "{} vs {}",
        sol.log_odds,
        clean.log_odds
    );
}

#[test]
fn short_random_lists_are_never_solved() {
    // Few sources leave little to outvote a chance alignment; the seed stars must not count.
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
