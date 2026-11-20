use fayde_api::{
    build_request, detect_sources, detected_trail, retry_as_points, DetectOptions, RequestOptions,
};
use fayde_catalog::{memory::MemoryCatalog, Star};
use fayde_core::{to_arcsec, PixelPos, Portable, SkyPos, UnitVector};
use fayde_detect::synth::{Psf, Rng, Scene, Star as Point, Streak, Trail};
use fayde_fits::Image;
use fayde_solve::{solve, Outcome};
use fayde_wcs::Header;
use fayde_wcs::TanWcs;

const SIZE: usize = 1024;
const SCALE: f64 = 3.0;

fn sky(centre: SkyPos, radius_deg: f64, rng: &mut Rng) -> Vec<Star> {
    let r = radius_deg.to_radians();
    let area =
        2.0 * std::f64::consts::PI * (1.0 - r.pcos()) * (180.0 / std::f64::consts::PI).ppowi(2);
    let c = centre.to_unit_vector().0;
    let e = {
        let (s, co) = centre.ra.psin_cos();
        [-s, co, 0.0]
    };
    let n = [
        c[1] * e[2] - c[2] * e[1],
        c[2] * e[0] - c[0] * e[2],
        c[0] * e[1] - c[1] * e[0],
    ];
    (0..(400.0 * area) as usize)
        .map(|i| {
            let cos_t = 1.0 - rng.uniform() * (1.0 - r.pcos());
            let sin_t = (1.0 - cos_t * cos_t).sqrt();
            let phi = 2.0 * std::f64::consts::PI * rng.uniform();
            let v =
                [0, 1, 2].map(|k| c[k] * cos_t + sin_t * (phi.pcos() * e[k] + phi.psin() * n[k]));
            let (lo, hi) = (10f64.ppowf(0.3 * 6.0), 10f64.ppowf(0.3 * 16.0));
            let mag = ((lo + rng.uniform() * (hi - lo)).plog10() / 0.3) as f32;
            Star {
                id: i as u64 + 1,
                pos: UnitVector(v).to_sky(),
                mag: Some(mag),
                pos_err_mas: [20.0, 20.0],
                pm_mas_yr: [Some(0.0), Some(0.0)],
            }
        })
        .collect()
}

struct Frame {
    catalog: MemoryCatalog,
    truth: TanWcs,
    scene: Scene,
    pixels: Vec<f32>,
}

fn frame(ra: f64, dec: f64, roll_deg: f64, psf: Psf, seed: u64) -> Frame {
    let mut rng = Rng(seed);
    let centre = SkyPos::from_degrees(ra, dec).unwrap();
    let stars = sky(centre, 1.0, &mut rng);
    let s = SCALE / 3600.0;
    let (sr, cr) = roll_deg.to_radians().psin_cos();
    let mid = (SIZE as f64 - 1.0) / 2.0;
    let truth = TanWcs::new(
        centre,
        PixelPos::new(mid, mid),
        [[-s * cr, s * sr], [s * sr, s * cr]],
    )
    .unwrap();
    let mut scene = Scene::new(SIZE, SIZE);
    scene.psf = psf;
    scene.gradient = (0.03, -0.02);
    scene.saturation = Some(65_535.0);
    scene.quantise = true;
    let faintest = 4.0 * scene.sky_noise() * 2.0 * std::f64::consts::PI.sqrt() * scene.psf_sigma();
    for st in &stars {
        let Some(p) = truth.sky_to_pixel(st.pos) else {
            continue;
        };
        if p.x < -20.0 || p.y < -20.0 || p.x > SIZE as f64 + 20.0 || p.y > SIZE as f64 + 20.0 {
            continue;
        }
        let flux = faintest * 10f64.ppowf(-0.4 * (f64::from(st.mag.unwrap()) - 16.0));
        scene.stars.push(Point {
            x: p.x,
            y: p.y,
            flux,
        });
    }
    let noise = scene.sky_noise();
    scene.streaks.push(Streak {
        x0: -30.0,
        y0: 310.0,
        x1: 1060.0,
        y1: 705.0,
        flux_per_px: 25.0 * noise,
    });
    scene.streaks.push(Streak {
        x0: 640.0,
        y0: 150.0,
        x1: 720.0,
        y1: 118.0,
        flux_per_px: 60.0 * noise,
    });
    for _ in 0..40 {
        scene.hot.push((
            (rng.uniform() * SIZE as f64) as usize,
            (rng.uniform() * SIZE as f64) as usize,
            noise * (20.0 + 500.0 * rng.uniform()),
        ));
    }
    let pixels = scene.render(seed ^ 0x5eed);
    Frame {
        catalog: MemoryCatalog::new(stars),
        truth,
        scene,
        pixels,
    }
}

fn solve_frame(f: &Frame, ra: f64, dec: f64, name: &str) -> f64 {
    let image = Image {
        width: SIZE,
        height: SIZE,
        data: f.pixels.clone(),
        header: Header::default(),
        hdu: 0,
    };
    let (sources, report) = detect_sources(&image, &DetectOptions::default()).unwrap();
    let streaks = report["streaks"].as_array().unwrap();
    eprintln!(
        "{name}: {} stars of {} detections, FWHM {:.3} px (true {:.3}), noise {:.2} (true {:.2}), \
         {} streaks, bin {}, {:.1} ms",
        sources.len(),
        report["sources"].as_array().unwrap().len(),
        report["fwhm_px"].as_f64().unwrap(),
        f.scene.psf_sigma() * 2.354_82,
        report["noise"].as_f64().unwrap(),
        f.scene.sky_noise(),
        streaks.len(),
        report["bin"],
        report["seconds"].as_f64().unwrap() * 1e3,
    );
    assert_eq!(report["count"], sources.len());
    assert_eq!(streaks.len(), 2, "{name}: {streaks:?}");
    let on_a_streak = |x: f64, y: f64| {
        f.scene.streaks.iter().any(|t| {
            let (dx, dy) = (t.x1 - t.x0, t.y1 - t.y0);
            let len = dx.phypot(dy);
            let along = ((x - t.x0) * dx + (y - t.y0) * dy) / len;
            let across = (-(x - t.x0) * dy + (y - t.y0) * dx) / len;
            across.abs() < 2.0 && along > 0.0 && along < len
        })
    };
    let near_a_star = |x: f64, y: f64| {
        f.scene
            .stars
            .iter()
            .any(|s| (s.x - x).phypot(s.y - y) < 3.0)
    };
    let fragments = sources
        .iter()
        .filter(|s| on_a_streak(s.x, s.y) && !near_a_star(s.x, s.y))
        .count();
    assert_eq!(
        fragments, 0,
        "{name}: streak fragments among the solver's sources"
    );
    let hot =
        sources
            .iter()
            .filter(|s| {
                f.scene.hot.iter().any(|&(hx, hy, _)| {
                    (hx as f64 - s.x).abs() < 1.0 && (hy as f64 - s.y).abs() < 1.0
                }) && !near_a_star(s.x, s.y)
            })
            .count();
    assert_eq!(hot, 0, "{name}: hot pixels among the solver's sources");

    let opts = RequestOptions {
        ra: Some(ra + 0.2),
        dec: Some(dec - 0.2),
        radius: Some(1.0),
        scale: Some(SCALE * 1.004),
        scale_sigma: Some(0.01),
        ..RequestOptions::default()
    };
    let built = build_request(sources, (Some(SIZE as u32), Some(SIZE as u32)), &opts).unwrap();
    let out = solve(&built.request, &f.catalog).unwrap();
    let Outcome::Solved(sol) = out else {
        panic!("{name}: {out:?}");
    };
    let mut worst = 0.0f64;
    let edge = SIZE as f64 - 1.0;
    for (x, y) in [
        (0.0, 0.0),
        (edge, 0.0),
        (0.0, edge),
        (edge, edge),
        (edge / 2.0, edge / 2.0),
    ] {
        let p = PixelPos::new(x, y);
        worst = worst.max(to_arcsec(
            sol.wcs.pixel_to_sky(p).separation(f.truth.pixel_to_sky(p)),
        ));
    }
    eprintln!(
        "{name}: solved, worst sky error {worst:.3}″, {} matches, rms {:.3} px, error scale {:.2}, \
         log-odds {:.0}, {} hypotheses",
        sol.matches.len(),
        sol.rms_px,
        sol.error_scale,
        sol.log_odds,
        sol.hypotheses_tested
    );
    worst
}

#[test]
fn a_frame_with_streaks_solves_from_its_pixels() {
    let f = frame(81.0, 12.0, 25.0, Psf::Gaussian(1.4), 11);
    let worst = solve_frame(&f, 81.0, 12.0, "gaussian 1.4 px");
    assert!(worst < 0.3, "{worst}″");
}

#[test]
fn seeing_limited_and_oversampled_stars_solve_too() {
    let f = frame(
        200.0,
        -35.0,
        160.0,
        Psf::Moffat {
            fwhm: 7.0,
            beta: 3.0,
        },
        12,
    );
    let worst = solve_frame(&f, 200.0, -35.0, "moffat fwhm 7 px");
    assert!(worst < 0.6, "{worst}″");
}

#[test]
fn undersampled_stars_solve_too() {
    let f = frame(310.0, 60.0, 270.0, Psf::Gaussian(0.6), 13);
    let worst = solve_frame(&f, 310.0, 60.0, "gaussian 0.6 px");
    assert!(worst < 0.3, "{worst}″");
}

#[test]
fn a_trailed_frame_solves_from_its_pixels() {
    let mut f = frame(140.0, -20.0, 70.0, Psf::Gaussian(1.3), 14);
    f.scene.trail = Some(Trail {
        length: 30.0,
        angle: 0.6,
    });
    f.scene.movers.push(Point {
        x: 512.3,
        y: 511.8,
        flux: 50_000.0,
    });
    let pixels = f.scene.render(14 ^ 0x5eed);
    let image = Image {
        width: SIZE,
        height: SIZE,
        data: pixels,
        header: Header::default(),
        hdu: 0,
    };
    let (sources, report) = detect_sources(&image, &DetectOptions::default()).unwrap();
    assert_eq!(report["mode"], "trails", "{}", report["census"]);
    let (len, angle) = detected_trail(&report).expect("a trail");
    assert!((len / 30.0 - 1.0).abs() < 0.03, "{len}");
    assert!((angle - 0.6f64.to_degrees()).abs() < 0.3, "{angle}");
    assert!(report["classes"]["point"].as_u64().unwrap() >= 1);
    assert!(sources
        .iter()
        .all(|s| (s.x - 512.3).phypot(s.y - 511.8) > 3.0));
    let opts = RequestOptions {
        ra: Some(140.2),
        dec: Some(-20.2),
        radius: Some(1.0),
        scale: Some(SCALE * 1.004),
        scale_sigma: Some(0.01),
        trail_length: Some(len),
        trail_angle: Some(angle),
        ..RequestOptions::default()
    };
    let built = build_request(sources, (Some(SIZE as u32), Some(SIZE as u32)), &opts).unwrap();
    let out = solve(&built.request, &f.catalog).unwrap();
    let Outcome::Solved(sol) = out else {
        panic!("{out:?}");
    };
    let edge = SIZE as f64 - 1.0;
    let worst = [
        (0.0, 0.0),
        (edge, 0.0),
        (0.0, edge),
        (edge, edge),
        (edge / 2.0, edge / 2.0),
    ]
    .iter()
    .map(|&(x, y)| {
        let p = PixelPos::new(x, y);
        to_arcsec(sol.wcs.pixel_to_sky(p).separation(f.truth.pixel_to_sky(p)))
    })
    .fold(0.0f64, f64::max);
    eprintln!(
        "trailed: {} matches, rms {:.3} px, worst sky error {worst:.3}″",
        sol.matches.len(),
        sol.rms_px
    );
    assert!(worst < 1.0, "{worst}″");
}
#[test]
fn a_curved_frame_solves_from_exposure_averaged_positions() {
    let mut f = frame(140.0, -20.0, 70.0, Psf::Gaussian(0.9), 14);
    let (sin, cos) = 0.40_f64.psin_cos();
    for star in std::mem::take(&mut f.scene.stars) {
        for i in 0..40 {
            let t = (i as f64 + 0.5) / 40.0;
            let u = 40.0 * t.min(0.5) - 15.0;
            let v = 40.0 * (t - 0.5).max(0.0) - 5.0;
            f.scene.stars.push(Point {
                x: star.x + cos * u - sin * v,
                y: star.y + sin * u + cos * v,
                flux: 3.0 * star.flux / 40.0,
            });
        }
    }
    let image = Image {
        width: SIZE,
        height: SIZE,
        data: f.scene.render(42),
        header: Header::default(),
        hdu: 0,
    };
    let (sources, report) = detect_sources(&image, &DetectOptions::default()).unwrap();
    assert_eq!(report["mode"], "kernel", "{}", report["census"]);
    assert_eq!(report["kernel"]["regions"], 1);
    assert!(report["kernel"]["bend_px"].as_f64().unwrap() > 3.0);
    assert!(detected_trail(&report).is_none());
    assert!(retry_as_points(&report));

    let opts = RequestOptions {
        ra: Some(140.2),
        dec: Some(-20.2),
        radius: Some(1.0),
        scale: Some(SCALE),
        scale_sigma: Some(0.01),
        ..RequestOptions::default()
    };
    let built = build_request(sources, (Some(SIZE as u32), Some(SIZE as u32)), &opts).unwrap();
    let out = solve(&built.request, &f.catalog).unwrap();
    let Outcome::Solved(sol) = out else {
        panic!("{out:?}");
    };
    let edge = SIZE as f64 - 1.0;
    for (x, y) in [
        (0.0, 0.0),
        (edge, 0.0),
        (0.0, edge),
        (edge, edge),
        (edge / 2.0, edge / 2.0),
    ] {
        let p = PixelPos::new(x, y);
        let error = to_arcsec(sol.wcs.pixel_to_sky(p).separation(f.truth.pixel_to_sky(p)));
        assert!(error < 0.3 * SCALE, "{error} arcsec at ({x}, {y})");
    }
}
