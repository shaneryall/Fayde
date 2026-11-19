use fayde_core::Portable;
use std::path::PathBuf;

use fayde_catalog::{ucac4::Ucac4Catalog, Catalog, StarQuery};
use fayde_core::{to_arcsec, Epoch, PixelPos, SkyPos};
use fayde_index::{build, BuildParams};
use fayde_solve::{
    solve, solve_with_index, Frame, Outcome, PointingPrior, Policy, Priors, ScalePrior,
    SolveMethod, SolveRequest, Source,
};
use fayde_wcs::TanWcs;

struct RealField {
    cat: Ucac4Catalog,
    req: SolveRequest,
    truth: TanWcs,
    truth_ids: Vec<Option<u64>>,
}

fn real_field() -> Option<RealField> {
    let dir = std::env::var_os("FAYDE_UCAC4_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data/ucac4"));
    if !(458..=463).all(|z| dir.join(format!("u4b/z{z}")).is_file()) {
        panic!("UCAC4 zones z458–z463 not at {}", dir.display());
    }
    let cat = Ucac4Catalog::open(&dir).unwrap();
    let epoch = Epoch::julian_year(2025.0);
    let centre = SkyPos::from_degrees(150.05, 1.98).unwrap();
    let s = 3.0 / 3600.0;
    let (sr, cr) = 23f64.to_radians().psin_cos();
    let truth = TanWcs::new(
        centre,
        PixelPos::new(199.5, 199.5),
        [[-s * cr, s * sr], [s * sr, s * cr]],
    )
    .unwrap();

    let q = StarQuery {
        epoch,
        mag_max: Some(15.5),
        ..StarQuery::default()
    };
    let stars = cat.cone(centre, 0.3f64.to_radians(), &q).unwrap();
    let mut seed = 11u64;
    let mut rnd = || {
        seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((seed >> 11) as f64 + 0.5) / (1u64 << 53) as f64
    };
    let gauss = |r: &mut dyn FnMut() -> f64| {
        (-2.0 * r().pln()).sqrt() * (2.0 * std::f64::consts::PI * r()).pcos()
    };
    let mut sources = Vec::new();
    let mut truth_ids = Vec::new();
    for st in &stars {
        let Some(p) = truth.sky_to_pixel(st.pos) else {
            continue;
        };
        if !(p.x >= 0.0 && p.y >= 0.0 && p.x < 399.0 && p.y < 399.0) || rnd() < 0.2 {
            continue;
        }
        let flux = 10f64.ppowf(-0.4 * (f64::from(st.mag.unwrap()) - 20.0));
        sources.push(Source {
            x: p.x + 0.3 * gauss(&mut rnd),
            y: p.y + 0.3 * gauss(&mut rnd),
            flux: Some(flux),
            sigma_px: Some(0.3),
            p_star: None,
            truncated: false,
        });
        truth_ids.push(Some(st.id));
    }
    let n_false = sources.len() * 3 / 10;
    for k in 0..n_false {
        let f = sources[k % sources.len()].flux.unwrap() * (0.5 + rnd());
        sources.push(Source {
            x: rnd() * 399.0,
            y: rnd() * 399.0,
            flux: Some(f),
            sigma_px: Some(0.3),
            p_star: None,
            truncated: false,
        });
        truth_ids.push(None);
    }
    let req = SolveRequest {
        frame: Frame {
            width: 400,
            height: 400,
            trail: None,
        },
        sources,
        priors: Priors {
            pointing: Some(PointingPrior {
                centre: SkyPos::from_degrees(150.0, 2.05).unwrap(),
                sigma: 0.1f64.to_radians(),
                max_radius: 0.2f64.to_radians(),
            }),
            scale: Some(ScalePrior {
                arcsec_per_px: 3.02,
                sigma: 0.03,
            }),
            roll: None,
            parity: None,
        },
        policy: Policy {
            epoch,
            ..Policy::default()
        },
    };
    Some(RealField {
        cat,
        req,
        truth,
        truth_ids,
    })
}

fn check(name: &str, f: &RealField, out: Outcome) -> SolveMethod {
    let Outcome::Solved(sol) = out else {
        panic!("{name}: did not solve: {out:?}")
    };
    let mut worst = 0.0f64;
    for (x, y) in [
        (0.0, 0.0),
        (399.0, 0.0),
        (0.0, 399.0),
        (399.0, 399.0),
        (199.5, 199.5),
    ] {
        let p = PixelPos::new(x, y);
        worst = worst.max(to_arcsec(
            sol.wcs.pixel_to_sky(p).separation(f.truth.pixel_to_sky(p)),
        ));
    }
    let conf: Vec<_> = sol.matches.iter().filter(|m| m.inlier_prob > 0.5).collect();
    let correct = conf
        .iter()
        .filter(|m| f.truth_ids[m.source] == Some(m.star_id))
        .count();
    eprintln!(
        "{name}: {:?}, worst sky error {worst:.3}″, {correct}/{} correct, NFA {:.1e}, {} hypotheses, {:?}",
        sol.method,
        conf.len(),
        sol.nfa,
        sol.hypotheses_tested,
        sol.elapsed
    );
    let bound = 5.0 * 0.3 * (14.0 / conf.len() as f64).sqrt() * 3.0;
    assert!(worst < bound, "{name}: {worst}″ > {bound}″");
    assert!(conf.len() >= 10 && correct == conf.len(), "{name}");
    assert_eq!(sol.catalog.name, "UCAC4");
    sol.method
}

#[test]
#[ignore = "needs UCAC4 data: set FAYDE_UCAC4_DIR, run with --include-ignored"]
fn solves_real_ucac4_field() {
    let Some(f) = real_field() else { return };
    let out = solve(&f.req, &f.cat).unwrap();
    assert_eq!(check("UCAC4 pair seeds", &f, out), SolveMethod::PairSeeds);
}

#[test]
#[ignore = "needs UCAC4 data: set FAYDE_UCAC4_DIR, run with --include-ignored"]
fn solves_real_ucac4_field_without_priors() {
    let Some(mut f) = real_field() else { return };
    let params = BuildParams {
        fov_deg: 0.5,
        min_edge_deg: 0.02,
        stars_per_fov: 30.0,
        mag_limit: Some(15.5),
        epoch_year: 2025.0,
        region: Some([150.0, 2.0, 0.38]),
        ..BuildParams::default()
    };
    let q = StarQuery {
        epoch: Epoch::julian_year(2025.0),
        ..StarQuery::default()
    };
    let centre = SkyPos::from_degrees(150.0, 2.0).unwrap();
    let stars = f.cat.cone(centre, 0.38f64.to_radians(), &q).unwrap();
    let index = build(&stars, f.cat.id(), &params).unwrap();
    f.req.priors = Priors::default();
    let out = solve_with_index(&f.req, &f.cat, &index).unwrap();
    assert_eq!(check("UCAC4 no priors", &f, out), SolveMethod::PatternHash);
}
