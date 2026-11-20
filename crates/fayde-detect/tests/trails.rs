use fayde_core::Portable;
use fayde_detect::synth::{Cosmic, Psf, Rng, Scene, Star, Streak, Trail};
use fayde_detect::{detect, Class, DetectOptions, Detected, Detection, FrameMode};

fn trailed(n: usize, len: f64, angle: f64, snr: (f64, f64), seed: u64) -> Scene {
    let mut sc = Scene::new(768, 768);
    sc.psf = Psf::Gaussian(1.3);
    sc.gain = 2.0;
    sc.trail = Some(Trail { length: len, angle });
    let mut rng = Rng(seed);
    let unit = sc.sky_noise() * 2.0 * std::f64::consts::PI.sqrt() * sc.psf_sigma();
    let mut tries = 0;
    while sc.stars.len() < n && tries < 100_000 {
        tries += 1;
        let (x, y) = (4.0 + rng.uniform() * 760.0, 4.0 + rng.uniform() * 760.0);
        if sc
            .stars
            .iter()
            .any(|s| (s.x - x).phypot(s.y - y) < len + 8.0)
        {
            continue;
        }
        let s = snr.0 * (snr.1 / snr.0).ppowf(rng.uniform());
        sc.stars.push(Star {
            x,
            y,
            flux: s * unit,
        });
    }
    sc
}

fn nearest(det: &Detection, x: f64, y: f64, r: f64) -> Option<&Detected> {
    det.sources
        .iter()
        .filter(|d| (d.x - x).phypot(d.y - y) < r)
        .min_by(|a, b| {
            (a.x - x)
                .phypot(a.y - y)
                .total_cmp(&(b.x - x).phypot(b.y - y))
        })
}

fn robust(v: &[f64]) -> f64 {
    let mut a: Vec<f64> = v.iter().map(|x| x.abs()).collect();
    a.sort_by(f64::total_cmp);
    1.4826 * a[a.len() / 2]
}

#[test]
fn narrow_trails_keep_hot_pixels_masked() {
    let mut sc = trailed(80, 24.0, 0.35, (100.0, 1000.0), 62);
    sc.psf = Psf::Gaussian(0.3);
    let mut rng = Rng(63);
    while sc.hot.len() < 40 {
        let x = (8.0 + rng.uniform() * 752.0) as usize;
        let y = (8.0 + rng.uniform() * 752.0) as usize;
        if sc
            .stars
            .iter()
            .all(|s| (s.x - x as f64).phypot(s.y - y as f64) > 16.0)
        {
            sc.hot.push((x, y, 4000.0));
        }
    }
    let opt = DetectOptions {
        psf_sigma: Some(0.4),
        ..DetectOptions::default()
    };
    let det = detect(&sc.render(64), sc.width, sc.height, &opt).unwrap();
    assert_eq!(det.census.mode, FrameMode::Trails, "{:?}", det.census);
    for &(x, y, _) in &sc.hot {
        let source = nearest(&det, x as f64, y as f64, 1.0);
        assert_eq!(source.map(|s| s.class), Some(Class::Sharp), "{x}, {y}");
    }
}

#[test]
fn bias_subtracted_trails_with_a_clipped_sky() {
    let scene = trailed(100, 40.0, -0.13, (12.0, 300.0), 29);
    let mut image = scene.render(36);
    let pedestal = (scene.sky + 1.4 * scene.sky_noise()) as f32;
    for value in &mut image {
        *value = (*value - pedestal).round().max(0.0);
    }
    let zeros = image.iter().filter(|&&v| v == 0.0).count();
    assert!(zeros > 9 * image.len() / 10);
    let opt = DetectOptions {
        gain: Some(scene.gain),
        ..DetectOptions::default()
    };
    let det = detect(&image, scene.width, scene.height, &opt).unwrap();
    assert_eq!(det.census.mode, FrameMode::Trails, "{:?}", det.census);
    let stars: Vec<_> = det
        .sources
        .iter()
        .filter(|d| d.class == Class::Star)
        .collect();
    let matches = stars
        .iter()
        .filter(|d| {
            scene
                .stars
                .iter()
                .any(|s| (d.x - s.x).phypot(d.y - s.y) < 3.0)
        })
        .count();
    let recovered = scene
        .stars
        .iter()
        .filter(|s| nearest(&det, s.x, s.y, 3.0).is_some())
        .count();
    assert!(
        matches * 10 >= 9 * stars.len(),
        "{matches}/{} matches",
        stars.len()
    );
    assert!(recovered >= 70, "{recovered}/100 recovered");
}

#[test]
fn trailed_frame_is_measured_as_trails() {
    for (len, deg, seed) in [(24.0, 20.0, 1u64), (36.0, -55.0, 2), (60.0, 80.0, 3)] {
        let angle = f64::to_radians(deg);
        let mut sc = trailed(150, len, angle, (8.0, 400.0), seed);
        let unit = sc.sky_noise() * 2.0 * std::f64::consts::PI.sqrt() * sc.psf_sigma();
        sc.movers.push(Star {
            x: 384.3,
            y: 383.6,
            flux: 200.0 * unit,
        });
        let mut rng = Rng(seed + 100);
        for _ in 0..40 {
            let (x, y) = (
                (rng.uniform() * 767.0) as usize,
                (rng.uniform() * 767.0) as usize,
            );
            sc.hot.push((x, y, 400.0));
        }
        for _ in 0..10 {
            sc.cosmics.push(Cosmic {
                x: rng.uniform() * 760.0,
                y: rng.uniform() * 760.0,
                length: rng.uniform() * 8.0,
                angle: rng.uniform() * 3.0,
                adu_per_px: 600.0,
            });
        }
        sc.streaks.push(Streak {
            x0: 20.0,
            y0: 700.0,
            x1: 740.0,
            y1: 60.0,
            flux_per_px: 6.0 * sc.sky_noise(),
        });
        let img = sc.render(seed + 7);
        let opt = DetectOptions {
            gain: Some(sc.gain),
            sigma_floor: 0.0,
            ..DetectOptions::default()
        };
        let det = detect(&img, 768, 768, &opt).unwrap();
        assert_eq!(
            det.census.mode,
            FrameMode::Trails,
            "{len} px: {:?}",
            det.census
        );
        let t = det.trail.unwrap();
        println!(
            "L {len} at {deg}°: measured {:.2} ± {:.2} px at {:.2}° ± {:.3}°, width {:.2}, fitted {}; census {:?}",
            t.length,
            t.sigma_length,
            t.angle.to_degrees(),
            t.sigma_angle.to_degrees(),
            t.width,
            t.fitted,
            det.census
        );
        assert!((t.length / len - 1.0).abs() < 0.05, "{len}: {}", t.length);
        let da = (t.angle - angle).psin().abs();
        assert!(da < 0.02, "{deg}: {}", t.angle.to_degrees());

        let (sn, cs) = angle.psin_cos();
        let s_psf = sc.psf_sigma();
        let to_trail = 2.0 * std::f64::consts::PI.sqrt() * s_psf
            / (4.0 * std::f64::consts::PI * s_psf * s_psf
                + 2.0 * std::f64::consts::PI.sqrt() * s_psf * len)
                .sqrt();
        let (mut ea, mut ec, mut za, mut zc) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
        let (mut found, mut strong) = (0, 0);
        for s in &sc.stars {
            let edge = 0.5 * len + 6.0;
            if s.x < edge || s.y < edge || s.x > 767.0 - edge || s.y > 767.0 - edge {
                continue;
            }
            let snr = s.flux / unit * to_trail;
            let candidate = det
                .sources
                .iter()
                .filter(|d| d.class == Class::Star)
                .filter(|d| {
                    let (dx, dy) = (d.x - s.x, d.y - s.y);
                    let (a, c) = (dx * cs + dy * sn, -dx * sn + dy * cs);
                    a.abs() <= (5.0 * d.sigma_along).clamp(2.0, 0.25 * len)
                        && c.abs() <= (5.0 * d.sigma_across).max(1.0)
                })
                .min_by(|a, b| {
                    (a.x - s.x)
                        .phypot(a.y - s.y)
                        .total_cmp(&(b.x - s.x).phypot(b.y - s.y))
                });
            if snr >= 10.0 {
                strong += 1;
            }
            let Some(d) = candidate else {
                if snr >= 10.0 {
                    let near = nearest(&det, s.x, s.y, 0.5 * len + 5.0);
                    println!(
                        "    missed {:.1},{:.1} trail snr {snr:.0}: nearest {:?}",
                        s.x,
                        s.y,
                        near.map(|d| (d.x, d.y, d.class))
                    );
                }
                continue;
            };
            if snr >= 10.0 {
                found += 1;
            }
            let (dx, dy) = (d.x - s.x, d.y - s.y);
            let (a, c) = (dx * cs + dy * sn, -dx * sn + dy * cs);
            ea.push(a);
            ec.push(c);
            za.push(a / d.sigma_along);
            zc.push(c / d.sigma_across);
        }
        let med = |v: &[f64]| {
            let mut a = v.to_vec();
            a.sort_by(f64::total_cmp);
            a[a.len() / 2]
        };
        let (bias_a, bias_c) = (med(&ea), med(&ec));
        println!(
            "  noise measured {:.2}, true {:.2}; rms along {:.3} across {:.3}",
            det.stats.noise,
            sc.sky_noise(),
            robust(&ea),
            robust(&ec)
        );
        println!(
            "  stars: {found} of {strong} above trail SNR 10 found; bias {bias_a:.3} / {bias_c:.3} \
             px; error over claim {:.2} along, {:.2} across",
            robust(&za),
            robust(&zc)
        );
        assert!(found as f64 >= 0.95 * strong as f64, "{found} of {strong}");
        let se = |v: &[f64]| 3.0 * 1.2533 * robust(v) / (v.len() as f64).sqrt();
        assert!(bias_a.abs() < se(&ea).max(0.02), "{bias_a}");
        assert!(bias_c.abs() < se(&ec).max(0.02), "{bias_c}");
        assert!((0.65..1.35).contains(&robust(&za)), "{}", robust(&za));
        assert!((0.65..1.35).contains(&robust(&zc)), "{}", robust(&zc));

        let m = nearest(&det, 384.3, 383.6, 3.0).expect("tracked object found");
        assert_eq!(m.class, Class::Point);
        let false_stars = det
            .stars()
            .filter(|d| {
                !sc.stars
                    .iter()
                    .any(|s| (s.x - d.x).phypot(s.y - d.y) < 0.5 * len + 3.0)
            })
            .count();
        println!(
            "  false stars: {false_stars}; streaks kept: {}",
            det.streaks.len()
        );
        for t in &det.streaks {
            println!(
                "    streak {:.0},{:.0} to {:.0},{:.0}: {:.0} px at {:.1}° cr {}",
                t.x0,
                t.y0,
                t.x1,
                t.y1,
                t.length,
                t.angle.to_degrees(),
                t.cosmic_ray
            );
        }
        for d in det.stars().filter(|d| {
            !sc.stars
                .iter()
                .any(|s| (s.x - d.x).phypot(s.y - d.y) < 0.5 * len + 3.0)
        }) {
            let (x0, y0, x1, y1) = (20.0, 700.0, 740.0, 60.0);
            let l = (x1 - x0).phypot(y1 - y0);
            let dist = ((d.x - x0) * (y1 - y0) - (d.y - y0) * (x1 - x0)).abs() / l;
            println!(
                "    false at {:.1},{:.1}: snr {:.1} p {:.2} flux {:.0}; {:.1} px from the satellite",
                d.x, d.y, d.snr, d.p_star, d.flux, dist
            );
        }
        assert!(false_stars <= 3, "{false_stars}");
    }
}

#[test]
fn point_frames_are_untouched() {
    for seed in 1..=3u64 {
        let mut sc = trailed(300, 0.0, 0.0, (8.0, 400.0), seed);
        sc.trail = None;
        sc.streaks.push(Streak {
            x0: 10.0,
            y0: 30.0,
            x1: 700.0,
            y1: 500.0,
            flux_per_px: 30.0,
        });
        let img = sc.render(seed);
        let on = detect(&img, 768, 768, &DetectOptions::default()).unwrap();
        let off = detect(
            &img,
            768,
            768,
            &DetectOptions {
                trails: false,
                ..DetectOptions::default()
            },
        )
        .unwrap();
        assert_eq!(on.census.mode, FrameMode::Points, "{:?}", on.census);
        assert_eq!(on.sources, off.sources);
        assert_eq!(on.streaks, off.streaks);
    }
}

#[test]
#[ignore]
fn error_claims_cover_the_truth() {
    for (len, deg) in [(24.0, 20.0), (36.0, -55.0), (60.0, 80.0)] {
        let angle = f64::to_radians(deg);
        let (sn, cs) = angle.psin_cos();
        let (mut ea, mut ec, mut za, mut zc) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
        for seed in 0..12u64 {
            let sc = trailed(120, len, angle, (8.0, 400.0), 1000 + seed);
            let img = sc.render(2000 + seed);
            let opt = DetectOptions {
                gain: Some(sc.gain),
                sigma_floor: 0.0,
                ..DetectOptions::default()
            };
            let det = detect(&img, 768, 768, &opt).unwrap();
            if det.census.mode != FrameMode::Trails {
                println!("  seed {seed}: {:?}", det.census);
                continue;
            }
            for s in &sc.stars {
                let edge = 0.5 * len + 6.0;
                if s.x < edge || s.y < edge || s.x > 767.0 - edge || s.y > 767.0 - edge {
                    continue;
                }
                let Some(d) = det
                    .stars()
                    .filter(|d| {
                        let (dx, dy) = (d.x - s.x, d.y - s.y);
                        let reach = (5.0 * d.sigma_along).clamp(2.0, 0.25 * len);
                        (dx * cs + dy * sn).abs() <= reach
                            && (-dx * sn + dy * cs).abs() <= (5.0 * d.sigma_across).max(1.0)
                    })
                    .min_by(|a, b| {
                        (a.x - s.x)
                            .phypot(a.y - s.y)
                            .total_cmp(&(b.x - s.x).phypot(b.y - s.y))
                    })
                else {
                    continue;
                };
                let (dx, dy) = (d.x - s.x, d.y - s.y);
                let (a, c) = (dx * cs + dy * sn, -dx * sn + dy * cs);
                ea.push(a);
                ec.push(c);
                za.push(a / d.sigma_along);
                zc.push(c / d.sigma_across);
            }
        }
        let within =
            |z: &[f64], k: f64| z.iter().filter(|v| v.abs() <= k).count() as f64 / z.len() as f64;
        let mean = |v: &[f64]| {
            let mut a = v.to_vec();
            a.sort_by(f64::total_cmp);
            a[a.len() / 2]
        };
        for (z, axis) in [(&za, "along"), (&zc, "across")] {
            let (one, two) = (within(z, 1.0), within(z, 2.0));
            assert!(
                (0.58..=0.78).contains(&one),
                "L {len} {axis}: {one} within 1σ"
            );
            assert!(two >= 0.86, "L {len} {axis}: {two} within 2σ");
        }
        println!(
            "L {len}: {} stars; median error {:.3} / {:.3} px; within 1σ {:.2} / {:.2}, 2σ {:.2} / {:.2}; spread {:.2} / {:.2}",
            ea.len(),
            mean(&ea),
            mean(&ec),
            within(&za, 1.0),
            within(&zc, 1.0),
            within(&za, 2.0),
            within(&zc, 2.0),
            robust(&za),
            robust(&zc)
        );
    }
}

#[test]
fn short_trails_stay_with_the_point_detector() {
    let sc = trailed(150, 12.0, 0.4, (8.0, 400.0), 9);
    let img = sc.render(9);
    let det = detect(&img, 768, 768, &DetectOptions::default()).unwrap();
    assert_eq!(det.census.mode, FrameMode::Points, "{:?}", det.census);
    assert!(
        det.census.reason.contains("shorter"),
        "{}",
        det.census.reason
    );
    assert!(det.trail.is_none());
}

#[test]
fn curved_trails_are_refused() {
    const LENGTH: f64 = 40.0;
    const SAMPLES: u32 = 80;
    let mut sc = trailed(150, LENGTH, 0.0, (30.0, 1000.0), 18);
    sc.psf = Psf::Gaussian(0.85);
    sc.trail = None;
    let turn = std::f64::consts::FRAC_PI_2;
    let radius = LENGTH / turn;
    sc.stars = sc
        .stars
        .iter()
        .flat_map(|s| {
            (0..SAMPLES).map(move |i| {
                let angle = turn * ((f64::from(i) + 0.5) / f64::from(SAMPLES) - 0.5);
                let (sn, cs) = angle.psin_cos();
                Star {
                    x: s.x + radius * sn,
                    y: s.y + radius * (1.0 - cs),
                    flux: s.flux / f64::from(SAMPLES),
                }
            })
        })
        .collect();
    let img = sc.render(18);
    let det = detect(&img, sc.width, sc.height, &DetectOptions::default()).unwrap();
    assert!(
        matches!(det.census.mode, FrameMode::Refused | FrameMode::Kernel),
        "{:?}",
        det.census
    );
    assert!(det.trail.is_none());
}

#[test]
fn aligned_segments_do_not_make_point_stars_trails() {
    let mut sc = trailed(150, 0.0, 0.0, (8.0, 400.0), 21);
    sc.trail = None;
    for k in 0..8 {
        let (x, y) = (60.0 + 80.0 * k as f64, 100.0 + 70.0 * k as f64);
        sc.streaks.push(Streak {
            x0: x,
            y0: y,
            x1: x + 40.0 * 0.9,
            y1: y + 40.0 * 0.44,
            flux_per_px: 40.0,
        });
    }
    let img = sc.render(21);
    let on = detect(&img, 768, 768, &DetectOptions::default()).unwrap();
    let off = detect(
        &img,
        768,
        768,
        &DetectOptions {
            trails: false,
            ..DetectOptions::default()
        },
    )
    .unwrap();
    assert_ne!(on.census.mode, FrameMode::Trails, "{:?}", on.census);
    assert_eq!(on.sources, off.sources);
}

#[test]
fn undefined_pixels_are_left_out() {
    let sc = trailed(150, 30.0, 0.5, (8.0, 400.0), 22);
    let mut img = sc.render(22);
    for y in 300..340 {
        for x in 300..420 {
            img[y * 768 + x] = f32::NAN;
        }
    }
    let det = detect(&img, 768, 768, &DetectOptions::default()).unwrap();
    assert_eq!(det.census.mode, FrameMode::Trails, "{:?}", det.census);
    for d in &det.sources {
        assert!(d.x.is_finite() && d.y.is_finite() && d.sigma_px.is_finite() && d.flux.is_finite());
    }
    let t = det.trail.unwrap();
    assert!(t.length.is_finite() && t.angle.is_finite() && t.width.is_finite());
}

#[test]
fn an_explicit_threshold_holds_for_trails() {
    let sc = trailed(150, 30.0, 0.5, (8.0, 400.0), 23);
    let img = sc.render(23);
    let n = |t: Option<f64>| {
        detect(
            &img,
            768,
            768,
            &DetectOptions {
                threshold: t,
                ..DetectOptions::default()
            },
        )
        .unwrap()
        .stars()
        .count()
    };
    let (default, high) = (n(None), n(Some(60.0)));
    assert!(high < default / 2, "{high} of {default}");
}

#[test]
fn long_trails_at_any_angle() {
    for (deg, seed) in [(42.0f64, 31u64), (47.0, 32), (-83.0, 33)] {
        let angle = deg.to_radians();
        let sc = trailed(40, 160.0, angle, (300.0, 3000.0), seed);
        let img = sc.render(seed);
        let det = detect(&img, 768, 768, &DetectOptions::default()).unwrap();
        assert_eq!(
            det.census.mode,
            FrameMode::Trails,
            "{deg}: {:?}",
            det.census
        );
        let t = det.trail.unwrap();
        assert!(
            (t.length / 160.0 - 1.0).abs() < 0.05,
            "{deg}: {} (fitted {})",
            t.length,
            t.fitted
        );
        assert!(
            (t.angle - angle).psin().abs() < 0.01,
            "{deg}: {}",
            t.angle.to_degrees()
        );
    }
}

#[test]
fn two_trail_populations_are_refused() {
    let mut sc = trailed(140, 30.0, 0.3, (20.0, 400.0), 41);
    let unit = sc.sky_noise() * 2.0 * std::f64::consts::PI.sqrt() * sc.psf_sigma();
    let mut rng = Rng(42);
    for _ in 0..140 {
        let (x, y) = (20.0 + rng.uniform() * 728.0, 20.0 + rng.uniform() * 728.0);
        let flux = 20.0 * (20.0f64).ppowf(rng.uniform()) * unit;
        let (sn, cs) = 1.6f64.psin_cos();
        sc.streaks.push(Streak {
            x0: x - 15.0 * cs,
            y0: y - 15.0 * sn,
            x1: x + 15.0 * cs,
            y1: y + 15.0 * sn,
            flux_per_px: flux / 30.0,
        });
    }
    let img = sc.render(41);
    let on = detect(&img, 768, 768, &DetectOptions::default()).unwrap();
    assert_ne!(on.census.mode, FrameMode::Trails, "{:?}", on.census);
    let chosen = detect(
        &img,
        768,
        768,
        &DetectOptions {
            trail_angle: Some(0.3),
            ..DetectOptions::default()
        },
    )
    .unwrap();
    if chosen.census.mode == FrameMode::Trails {
        let a = chosen.trail.unwrap().angle;
        assert!((a - 0.3).abs() < 0.02, "{}", a.to_degrees());
    }
}

#[test]
fn trails_longer_than_the_frame_are_refused() {
    let mut sc = Scene::new(512, 512);
    sc.psf = Psf::Gaussian(1.3);
    sc.trail = Some(Trail {
        length: 1200.0,
        angle: 0.2,
    });
    let unit = sc.sky_noise() * 2.0 * std::f64::consts::PI.sqrt() * sc.psf_sigma();
    let mut rng = Rng(51);
    for _ in 0..30 {
        sc.stars.push(Star {
            x: rng.uniform() * 512.0,
            y: rng.uniform() * 512.0,
            flux: 3000.0 * unit,
        });
    }
    let img = sc.render(51);
    let det = detect(&img, 512, 512, &DetectOptions::default()).unwrap();
    assert_ne!(det.census.mode, FrameMode::Trails, "{:?}", det.census);
    assert!(det.trail.is_none());
}

#[test]
fn trails_are_pinned_to_the_bit() {
    let mut sc = trailed(80, 28.0, 0.7, (10.0, 400.0), 61);
    sc.quantise = true;
    let img = sc.render(61);
    let det = detect(&img, 768, 768, &DetectOptions::default()).unwrap();
    assert_eq!(det.census.mode, FrameMode::Trails);
    let mut h = 0xcbf2_9ce4_8422_2325u64;
    let mut put = |v: f64| {
        for b in v.to_bits().to_le_bytes() {
            h = (h ^ u64::from(b)).wrapping_mul(0x0100_0000_01b3);
        }
    };
    let t = det.trail.unwrap();
    for v in [t.length, t.angle, t.width, t.sigma_length, t.sigma_angle] {
        put(v);
    }
    for s in &det.sources {
        for v in [
            s.x,
            s.y,
            s.sigma_along,
            s.sigma_across,
            s.flux,
            s.snr,
            s.p_star,
        ] {
            put(v);
        }
    }
    assert_eq!(
        h,
        TRAILS_PINNED,
        "fingerprint {h:#018x} ({} sources)",
        det.sources.len()
    );
}

const TRAILS_PINNED: u64 = 0x9554_0a8f_8220_e592;

fn assert_recovered_trails(sc: &Scene, img: &[f32]) {
    let det = detect(img, sc.width, sc.height, &DetectOptions::default()).unwrap();
    assert_eq!(det.census.mode, FrameMode::Trails, "{:?}", det.census);
    let truth = sc.trail.unwrap();
    let trail = det.trail.unwrap();
    assert!((trail.angle - truth.angle).psin().abs() < 0.03, "{trail:?}");
    assert!((trail.length / truth.length - 1.0).abs() < 0.1, "{trail:?}");
    let recovered = sc
        .stars
        .iter()
        .filter(|s| nearest(&det, s.x, s.y, 3.0).is_some_and(|d| d.class == Class::Star))
        .count();
    assert!(
        recovered * 5 >= sc.stars.len() * 4,
        "{recovered} of {} stars recovered",
        sc.stars.len()
    );
}

#[test]
fn cosmic_ray_tracks_do_not_outvote_stars() {
    let mut sc = trailed(100, 40.0, -0.48, (40.0, 400.0), 72);
    sc.psf = Psf::Gaussian(0.85);
    let mut rng = Rng(92);
    for _ in 0..480 {
        sc.cosmics.push(Cosmic {
            x: rng.uniform() * 767.0,
            y: rng.uniform() * 767.0,
            angle: rng.uniform() * std::f64::consts::PI,
            length: 3.0 + rng.uniform() * 9.0,
            adu_per_px: 400.0,
        });
    }
    assert_recovered_trails(&sc, &sc.render(10));
}

/// A radiation storm: hits, short tracks and hot pixels outnumber the stars twenty to one, and
/// none of the stars is round. Most hits spread their charge into the neighbouring pixels (σ 0.4
/// to 0.6 px), so they share a size, as stars do: on real frames that mode, not the stars, was
/// measured as the frame's width. The width must still be the stars' across their trail.
fn storm(seed: u64) -> (Scene, Vec<f32>) {
    let mut sc = trailed(100, 35.0, -0.48, (40.0, 400.0), seed);
    sc.psf = Psf::Gaussian(1.1);
    let mut rng = Rng(seed + 1);
    for _ in 0..300 {
        sc.cosmics.push(Cosmic {
            x: rng.uniform() * 767.0,
            y: rng.uniform() * 767.0,
            angle: rng.uniform() * std::f64::consts::PI,
            length: 3.0 + rng.uniform() * 12.0,
            adu_per_px: 300.0 + rng.uniform() * 2700.0,
        });
    }
    while sc.hot.len() < 400 {
        let (x, y) = (
            (rng.uniform() * 767.0) as usize,
            (rng.uniform() * 767.0) as usize,
        );
        sc.hot.push((x, y, 200.0 + rng.uniform() * 4000.0));
    }
    let mut img = sc.render(seed + 2);
    let (w, h) = (sc.width, sc.height);
    for _ in 0..2500 {
        let (cx, cy) = (2.0 + rng.uniform() * 763.0, 2.0 + rng.uniform() * 763.0);
        let s = 0.4 + rng.uniform() * 0.2;
        let total = 1500.0 + rng.uniform() * 15000.0;
        let norm = total / (std::f64::consts::TAU * s * s);
        for y in (cy as usize).saturating_sub(3)..((cy as usize) + 4).min(h) {
            for x in (cx as usize).saturating_sub(3)..((cx as usize) + 4).min(w) {
                let r2 = (x as f64 - cx).ppowi(2) + (y as f64 - cy).ppowi(2);
                img[y * w + x] += (norm * (-0.5 * r2 / (s * s)).pexp()) as f32;
            }
        }
    }
    (sc, img)
}

#[test]
fn radiation_storm_does_not_set_the_star_width() {
    for seed in [140, 150] {
        let (sc, img) = storm(seed);
        let det = detect(&img, sc.width, sc.height, &DetectOptions::default()).unwrap();
        let sigma = det.stats.psf_sigma;
        assert!(
            (sigma / sc.psf_sigma() - 1.0).abs() < 0.25,
            "seed {seed}: measured σ {sigma:.2} px for stars of σ {:.2} px",
            sc.psf_sigma()
        );
        assert_recovered_trails(&sc, &img);
    }
}

#[test]
fn smooth_background_does_not_set_the_census_direction() {
    let sc = trailed(100, 40.0, -0.48, (40.0, 400.0), 72);
    let mut img = sc.render(10);
    let amp = 2.0 * sc.sky_noise();
    let tau = std::f64::consts::TAU;
    for (i, v) in img.iter_mut().enumerate() {
        let (x, y) = ((i % sc.width) as f64, (i / sc.width) as f64);
        let field =
            (tau * (x + 0.3 * y) / 130.0).pcos() + 0.7 * (tau * (y - 0.4 * x) / 170.0).psin();
        *v += (amp * field) as f32;
    }
    assert_recovered_trails(&sc, &img);
}

#[test]
fn a_bright_tracked_object_is_measured_once_where_it_is() {
    // The trail-matched filter spreads a bright point into a ridge of seeds about two trail
    // lengths long. Two of them used to survive, and each fit took the other as a neighbour
    // carrying the whole flux: the object was lost, or reported pixels away with a fraction
    // of its flux. A bright object was also flagged as a hit beside its point.
    for (dim_sky, snr, seed) in [
        (false, 200.0, 71u64),
        (false, 600.0, 72),
        (false, 2000.0, 73),
        (true, 200.0, 74),
        (true, 700.0, 75),
        (true, 2000.0, 76),
    ] {
        let mut sc = trailed(120, 48.0, f64::to_radians(32.0), (8.0, 400.0), seed);
        if dim_sky {
            sc.psf = Psf::Gaussian(1.6);
            sc.sky = 200.0;
            sc.gain = 1.0;
            sc.read_noise = 6.0;
        }
        let unit = sc.sky_noise() * 2.0 * std::f64::consts::PI.sqrt() * sc.psf_sigma();
        let (x, y) = (384.3, 383.6);
        sc.stars.retain(|s| (s.x - x).phypot(s.y - y) > 60.0);
        let flux = snr * unit;
        sc.movers.push(Star { x, y, flux });
        let det = detect(
            &sc.render(seed),
            sc.width,
            sc.height,
            &DetectOptions::default(),
        )
        .unwrap();
        assert_eq!(det.census.mode, FrameMode::Trails, "{:?}", det.census);
        let near: Vec<&Detected> = det
            .sources
            .iter()
            .filter(|d| (d.x - x).phypot(d.y - y) < 30.0)
            .collect();
        assert_eq!(near.len(), 1, "SNR {snr}, seed {seed}: {near:?}");
        let m = near[0];
        assert_eq!(m.class, Class::Point, "SNR {snr}, seed {seed}");
        assert!((m.x - x).phypot(m.y - y) < 0.3, "SNR {snr}: {m:?}");
        assert!(
            (m.flux / flux - 1.0).abs() < 0.05,
            "SNR {snr}: {} vs {flux}",
            m.flux
        );
    }
}

/// A crowded track-rate field: 48 px trails that may overlap, faint to bright, radiation hits
/// a few pixels long, on a dim sky.
fn crowded(n: usize, seed: u64) -> Scene {
    let mut sc = Scene::new(1024, 1024);
    sc.psf = Psf::Gaussian(1.6);
    sc.sky = 200.0;
    sc.gain = 1.0;
    sc.read_noise = 6.0;
    sc.trail = Some(Trail {
        length: 48.0,
        angle: f64::to_radians(32.0),
    });
    let mut rng = Rng(seed);
    let unit = sc.sky_noise() * 2.0 * std::f64::consts::PI.sqrt() * sc.psf_sigma();
    for _ in 0..n {
        let snr = 10.0 * 300f64.ppowf(rng.uniform().ppowf(1.6));
        sc.stars.push(Star {
            x: -20.0 + rng.uniform() * 1064.0,
            y: -20.0 + rng.uniform() * 1064.0,
            flux: snr * unit,
        });
    }
    for _ in 0..120 {
        sc.cosmics.push(Cosmic {
            x: 2.0 + rng.uniform() * 1020.0,
            y: 2.0 + rng.uniform() * 1020.0,
            length: rng.uniform() * 5.0,
            angle: rng.uniform() * std::f64::consts::PI,
            adu_per_px: 800.0 + 5200.0 * rng.uniform(),
        });
    }
    sc
}

#[test]
fn hits_in_a_crowded_trailed_frame_are_flagged_once() {
    // The round peaks of a trailed frame are hits, and their thinness was taken for the stars':
    // no hit was sharper than that, and hits came back as tracked objects instead.
    for (n, seed) in [(260usize, 81u64), (350, 82)] {
        let sc = crowded(n, seed);
        let det = detect(
            &sc.render(seed),
            sc.width,
            sc.height,
            &DetectOptions::default(),
        )
        .unwrap();
        assert_eq!(det.census.mode, FrameMode::Trails, "{:?}", det.census);
        let near = |c: &Cosmic, class: Class, r: f64| {
            det.sources
                .iter()
                .filter(|s| s.class == class && (s.x - c.x).phypot(s.y - c.y) < r)
                .count()
        };
        let flagged = sc
            .cosmics
            .iter()
            .filter(|c| near(c, Class::Sharp, 6.0) > 0)
            .count();
        assert!(flagged >= 114, "{n} trails: {flagged} of 120 hits flagged");
        let as_points: usize = sc.cosmics.iter().map(|c| near(c, Class::Point, 6.0)).sum();
        assert_eq!(as_points, 0, "{n} trails: hits reported as tracked objects");
        let sharp = det
            .sources
            .iter()
            .filter(|s| s.class == Class::Sharp)
            .count();
        assert!(
            sharp <= 126,
            "{n} trails: {sharp} sharp sources for 120 hits"
        );
    }
}

#[test]
fn chains_of_trails_are_not_streaks_but_satellites_are() {
    // The streak tracer joined roughly collinear trails into streaks nearly along them; the
    // trails on such a chain were then lost with it. Satellites crossing nearly along the
    // trails must still be found.
    let mut sc = crowded(350, 83);
    let noise = sc.sky_noise();
    let lines = [
        (-20.0, 300.0, 1040.0, 728.0, 12.0),
        (100.0, 1000.0, 1000.0, 489.0, 8.0),
        (300.0, 80.0, 700.0, 350.0, 20.0),
    ];
    for &(x0, y0, x1, y1, snr) in &lines {
        sc.streaks.push(Streak {
            x0,
            y0,
            x1,
            y1,
            flux_per_px: snr * noise * (2.0 * std::f64::consts::PI).sqrt() * sc.psf_sigma(),
        });
    }
    let det = detect(
        &sc.render(83),
        sc.width,
        sc.height,
        &DetectOptions::default(),
    )
    .unwrap();
    assert_eq!(det.census.mode, FrameMode::Trails, "{:?}", det.census);
    let off_line = |s: &fayde_detect::Streak, l: &(f64, f64, f64, f64, f64)| {
        let (dx, dy) = (l.2 - l.0, l.3 - l.1);
        let len = dx.phypot(dy);
        [(s.x0, s.y0), (s.x1, s.y1)]
            .iter()
            .map(|&(x, y)| ((x - l.0) * dy - (y - l.1) * dx).abs() / len)
            .fold(0.0, f64::max)
    };
    for l in &lines {
        assert!(
            det.streaks.iter().any(|s| off_line(s, l) < 3.0),
            "satellite {l:?} not found: {:?}",
            det.streaks
        );
    }
    let spurious: Vec<_> = det
        .streaks
        .iter()
        .filter(|s| lines.iter().all(|l| off_line(s, l) >= 3.0))
        .collect();
    assert!(
        spurious.len() <= 2,
        "{} spurious: {spurious:?}",
        spurious.len()
    );
}
