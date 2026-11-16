//! Pattern-hash solving: loose or missing priors, using a regional pattern index built from the
//! same synthetic sky.

use fayde_core::Portable;
mod common;

use common::{field, Field, Rng};
use fayde_catalog::{memory::MemoryCatalog, Catalog};
use fayde_core::{to_arcsec, PixelPos};
use fayde_index::{build, BuildParams, Index};
use fayde_solve::{solve_with_index, solve_with_indexes, Outcome, Priors, ScalePrior, SolveMethod};
use fayde_wcs::Parity;

/// Index over the synthetic sky around `centre`, for frames up to 1.3° across.
fn index_for(cat: &MemoryCatalog, centre: (f64, f64)) -> Index {
    let params = BuildParams {
        fov_deg: 1.3,
        min_edge_deg: 0.05,
        stars_per_fov: 40.0,
        region: Some([centre.0, centre.1, 1.6]),
        ..BuildParams::default()
    };
    build(cat.stars(), cat.id(), &params).unwrap()
}

fn check(name: &str, cat: &MemoryCatalog, index: &Index, f: &Field) {
    let out = solve_with_index(&f.request, cat, index).unwrap();
    let Outcome::Solved(sol) = out else {
        panic!("{name}: {out:?}")
    };
    let mut worst = 0.0f64;
    for (x, y) in [
        (0.0, 0.0),
        (1023.0, 0.0),
        (0.0, 1023.0),
        (1023.0, 1023.0),
        (511.5, 511.5),
    ] {
        let p = PixelPos::new(x, y);
        worst = worst.max(to_arcsec(
            sol.wcs.pixel_to_sky(p).separation(f.truth.pixel_to_sky(p)),
        ));
    }
    let conf: Vec<_> = sol.matches.iter().filter(|m| m.inlier_prob > 0.5).collect();
    let correct = conf
        .iter()
        .filter(|m| f.source_star[m.source] == Some(m.star_id))
        .count();
    eprintln!(
        "{name}: {:?}, worst sky error {worst:.3}″, {correct}/{} correct, NFA {:.1e}, {} hypotheses, {:?}",
        sol.method,
        conf.len(),
        sol.nfa,
        sol.hypotheses_tested,
        sol.elapsed
    );
    assert_eq!(sol.method, SolveMethod::PatternHash, "{name}");
    // 5σ noise bound (see ucac4_solve.rs): σ = 0.3 px, 3″/px.
    let bound = 5.0 * 0.3 * (14.0 / conf.len() as f64).sqrt() * 3.0;
    assert!(worst < bound, "{name}: {worst}″ > {bound}″");
    assert!(correct == conf.len() && conf.len() >= 20, "{name}");
    assert!(sol.nfa < 1e-9);
}

#[test]
fn loose_scale_prior_uses_pattern_hash() {
    let (cat, mut f) = field((150.0, 2.0), 37.0, Parity::Negative, 1);
    let index = index_for(&cat, (150.0, 2.0));
    f.request.priors.scale = Some(ScalePrior {
        arcsec_per_px: 3.3,
        sigma: 0.5,
    }); // ±17%, off by 10%
    check("loose scale", &cat, &index, &f);
}

#[test]
fn solves_with_no_priors_at_all() {
    let (cat, mut f) = field((210.0, 45.0), 200.0, Parity::Positive, 4);
    let index = index_for(&cat, (210.0, 45.0));
    f.request.priors = Priors::default();
    check("no priors", &cat, &index, &f);
}

#[test]
fn tight_priors_still_use_pair_seeds() {
    let (cat, f) = field((150.0, 2.0), 37.0, Parity::Negative, 1);
    let index = index_for(&cat, (150.0, 2.0));
    let Outcome::Solved(sol) = solve_with_index(&f.request, &cat, &index).unwrap() else {
        panic!()
    };
    assert_eq!(sol.method, SolveMethod::PairSeeds);
}

#[test]
fn random_sources_without_priors_abstain() {
    let (cat, mut f) = field((150.0, 2.0), 37.0, Parity::Negative, 6);
    let index = index_for(&cat, (150.0, 2.0));
    f.request.priors = Priors::default();
    let mut rng = Rng(5);
    for s in &mut f.request.sources {
        s.x = rng.uniform() * 1023.0;
        s.y = rng.uniform() * 1023.0;
    }
    match solve_with_index(&f.request, &cat, &index).unwrap() {
        Outcome::Abstained {
            hypotheses_tested,
            best_log_odds,
            ..
        } => {
            eprintln!(
                "random/no priors: abstained after {hypotheses_tested}, best {best_log_odds:?}"
            )
        }
        Outcome::Solved(s) => panic!("random sources solved: log-odds {}", s.log_odds),
    }
}

/// Indexes in √2 field-size steps over the field's synthetic sky.
fn bands(cat: &MemoryCatalog, centre: (f64, f64)) -> Vec<Index> {
    [0.65, 0.92, 1.3]
        .iter()
        .map(|&fov| {
            let params = BuildParams {
                fov_deg: fov,
                min_edge_deg: fov / 25.0,
                stars_per_fov: 40.0,
                region: Some([centre.0, centre.1, 1.6]),
                ..BuildParams::default()
            };
            build(cat.stars(), cat.id(), &params).unwrap()
        })
        .collect()
}

#[test]
fn blind_multi_band_solves_different_scales() {
    for (scale, seed) in [(1.5, 21u64), (3.0, 22)] {
        let (cat, mut f) = common::field_scaled((60.0, -20.0), 75.0, Parity::Negative, seed, scale);
        let idx = bands(&cat, (60.0, -20.0));
        let fov_of = |h: u64| {
            idx.iter()
                .find(|i| i.content_hash() == h)
                .map(|i| i.meta().params.fov_deg)
        };
        let real = f.source_star.iter().filter(|s| s.is_some()).count();

        // Blind: no priors at all.
        f.request.priors = Priors::default();
        let out = solve_with_indexes(&f.request, &cat, &idx).unwrap();
        let Outcome::Solved(sol) = out else {
            panic!("{scale}″/px blind ({real} real sources): {out:?}")
        };
        let p = PixelPos::new(0.0, 0.0);
        let err = to_arcsec(sol.wcs.pixel_to_sky(p).separation(f.truth.pixel_to_sky(p)));
        let z = err / sol.covariance.radial_sigma_arcsec(p);
        eprintln!(
            "{scale}″/px blind: band {:?}°, corner error {err:.2}″ ({z:.1}σ), {} hypotheses, {:?}",
            sol.index_hash.and_then(fov_of),
            sol.hypotheses_tested,
            sol.elapsed
        );
        assert!(z < 5.0);
        assert!((sol.wcs.pixel_scale_arcsec() / scale - 1.0).abs() < 0.01);

        // Loose scale prior: the best-matched band is searched first and solves.
        f.request.priors.scale = Some(ScalePrior {
            arcsec_per_px: scale * 1.1,
            sigma: scale * 0.2,
        });
        let Outcome::Solved(sol) = solve_with_indexes(&f.request, &cat, &idx).unwrap() else {
            panic!()
        };
        let diag = 1024.0 * std::f64::consts::SQRT_2 * scale * 1.1 / 3600.0;
        let nearest = idx
            .iter()
            .min_by(|a, b| {
                let d = |i: &Index| (i.meta().params.fov_deg / diag).pln().abs();
                d(a).total_cmp(&d(b))
            })
            .unwrap();
        assert_eq!(sol.index_hash, Some(nearest.content_hash()), "{scale}″/px");
    }
}

#[test]
fn pattern_prior_rejections_do_not_count_as_hypotheses() {
    let (cat, mut f) = field((150.0, 2.0), 37.0, Parity::Negative, 1);
    let index = index_for(&cat, (150.0, 2.0));
    f.request.priors.scale = None;
    f.request.priors.roll = Some(fayde_solve::RollPrior {
        angle: 1.23456789,
        sigma: 0.1,
        max_dev: 0.0,
    });
    let out = solve_with_index(&f.request, &cat, &index).unwrap();
    let Outcome::Abstained {
        hypotheses_tested, ..
    } = out
    else {
        panic!("{out:?}")
    };
    eprintln!("pattern prior rejects counted: {hypotheses_tested}");
    assert_eq!(hypotheses_tested, 0);
}

#[test]
fn near_rectangular_pattern_survives_canonical_order_swap() {
    use fayde_catalog::Star;
    use fayde_core::SkyPos;
    use fayde_index::{
        image_shape,
        pattern::{angle, shape_of},
    };
    use fayde_solve::{Frame, Policy, SolveRequest, Source};
    use fayde_wcs::TanWcs;

    for parity in [Parity::Negative, Parity::Positive] {
        let sign = if parity == Parity::Negative {
            -1.0
        } else {
            1.0
        };
        let truth = TanWcs::new(
            SkyPos::from_degrees(150.0, 2.0).unwrap(),
            PixelPos::new(511.5, 511.5),
            [[sign * 3.0 / 3600.0, 0.0], [0.0, 3.0 / 3600.0]],
        )
        .unwrap();
        let pixels = [
            [200.0, 250.0],
            [800.0, 250.04],
            [800.0, 750.0],
            [200.0, 750.2],
            [350.0, 400.0],
            [650.0, 330.0],
            [475.0, 650.0],
            [850.0, 550.0],
            [120.0, 480.0],
        ];
        let stars: Vec<Star> = pixels
            .iter()
            .enumerate()
            .map(|(i, p)| Star {
                id: i as u64 + 1,
                pos: truth.pixel_to_sky(PixelPos::new(p[0], p[1])),
                mag: Some(8.0 + i as f32),
                pos_err_mas: [0.0; 2],
                pm_mas_yr: [Some(0.0); 2],
            })
            .collect();
        let mut sources: Vec<Source> = pixels
            .iter()
            .enumerate()
            .map(|(i, p)| Source {
                x: p[0],
                y: p[1],
                flux: Some(100.0 - i as f64),
                sigma_px: Some(0.3),
                p_star: None,
                truncated: false,
            })
            .collect();
        sources[0].x -= 0.4;
        sources[2].y += 0.3;
        let image = std::array::from_fn(|i| [sources[i].x, sources[i].y]);
        let sky = std::array::from_fn(|i| stars[i].pos.to_unit_vector().0);
        let (_, image_order) = image_shape(&image);
        let (_, sky_order) = shape_of(&sky, angle);
        eprintln!("{parity:?}: image order {image_order:?}, catalogue order {sky_order:?}");
        assert_ne!(
            image_order, sky_order,
            "fixture must exercise the centroid-order swap"
        );
        let cat = MemoryCatalog::new(stars);
        let index = build(
            cat.stars(),
            cat.id(),
            &BuildParams {
                fov_deg: 1.3,
                min_edge_deg: 0.01,
                stars_per_fov: 100.0,
                pattern_stars_per_fov: 100.0,
                per_cell: 16,
                region: Some([150.0, 2.0, 1.0]),
                patterns_per_field: 200,
                ..BuildParams::default()
            },
        )
        .unwrap();
        let req = SolveRequest {
            frame: Frame {
                width: 1024,
                height: 1024,
                trail: None,
            },
            sources,
            priors: Priors::default(),
            policy: Policy {
                pattern_sources: 4,
                max_sip_order: 1,
                ..Policy::default()
            },
        };
        let out = solve_with_index(&req, &cat, &index).unwrap();
        let Outcome::Solved(sol) = out else {
            panic!("canonical swap: {out:?}")
        };
        eprintln!(
            "{parity:?}: {} hypotheses, RMS {:.4}",
            sol.hypotheses_tested, sol.rms_px
        );
        assert_eq!(sol.method, SolveMethod::PatternHash);
        assert_eq!(sol.wcs.parity(), parity);
        assert!(sol.rms_px < 0.3);
        assert_eq!(sol.matches.len(), 9);
        assert!(sol.matches.iter().all(|m| m.star_id == m.source as u64 + 1));

        // With only the seed stars, every surviving assignment fails the extra-hit
        // pre-filter, but each must still consume the hypothesis budget.
        let mut seed_only = req.clone();
        seed_only.sources.truncate(4);
        seed_only.policy.max_hypotheses = 2;
        let out = solve_with_index(&seed_only, &cat, &index).unwrap();
        let Outcome::Abstained {
            reason,
            hypotheses_tested,
            ..
        } = out
        else {
            panic!("seed-only field must fail the pre-filter: {out:?}")
        };
        assert!(
            reason.starts_with("hypothesis limit (2) reached during pattern hash"),
            "{reason}"
        );
        assert_eq!(hypotheses_tested, 2);
    }
}
