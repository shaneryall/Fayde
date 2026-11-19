use fayde_core::Portable;
use fayde_detect::synth::{Cosmic, Psf, Rng, Scene, Star, Stray, Streak};
use fayde_detect::{detect, Class, DetectOptions, Detected, Detection};

fn flux_for(scene: &Scene, snr: f64) -> f64 {
    snr * scene.sky_noise() * 2.0 * std::f64::consts::PI.sqrt() * scene.psf_sigma()
}

fn scatter(scene: &mut Scene, n: usize, snr: (f64, f64), rng: &mut Rng) {
    let (w, h) = (scene.width as f64, scene.height as f64);
    while scene.stars.len() < n {
        let (x, y) = (
            10.0 + rng.uniform() * (w - 20.0),
            10.0 + rng.uniform() * (h - 20.0),
        );
        if scene.stars.iter().any(|s| (s.x - x).phypot(s.y - y) < 12.0) {
            continue;
        }
        let s = snr.0 * (snr.1 / snr.0).ppowf(rng.uniform());
        let flux = flux_for(scene, s);
        scene.stars.push(Star { x, y, flux });
    }
}

fn nearest(det: &Detection, x: f64, y: f64, r: f64) -> Option<&Detected> {
    det.stars()
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
fn positions_are_as_good_as_claimed_and_as_the_bound_allows() {
    for (psf, count, seed) in [(0.9, 500, 1u64), (1.5, 400, 2), (2.6, 150, 3)] {
        let mut scene = Scene::new(768, 768);
        scene.psf = Psf::Gaussian(psf);
        scene.gain = 4.0;
        scatter(&mut scene, count, (10.0, 2000.0), &mut Rng(seed));
        let img = scene.render(seed + 50);
        let opt = DetectOptions {
            gain: Some(scene.gain),
            sigma_floor: 0.0,
            ..DetectOptions::default()
        };
        let det = detect(&img, 768, 768, &opt).unwrap();
        let width = (psf * psf + 1.0 / 12.0).sqrt();
        assert!(
            (det.stats.psf_sigma / width - 1.0).abs() < 0.05,
            "psf {psf}: measured {}",
            det.stats.psf_sigma
        );
        let noise = det.stats.noise / scene.sky_noise();
        assert!(
            (0.97..1.12).contains(&noise),
            "psf {psf}: noise ratio {noise}"
        );
        let (mut z, mut faint) = (Vec::new(), Vec::new());
        let mut missed = 0;
        for st in &scene.stars {
            let Some(d) = nearest(&det, st.x, st.y, 3.0) else {
                missed += 1;
                continue;
            };
            let snr = st.flux / flux_for(&scene, 1.0);
            for e in [d.x - st.x, d.y - st.y] {
                z.push(e / d.sigma_px);
                if snr < 40.0 {
                    faint.push(e * snr / (std::f64::consts::SQRT_2 * width));
                }
            }
        }
        assert!(missed <= 2, "psf {psf}: {missed} of {count} stars missed");
        let spread = robust(&z);
        assert!(
            (0.8..1.25).contains(&spread),
            "psf {psf}: error / claimed = {spread}"
        );
        let bound = robust(&faint);
        assert!(
            (0.8..1.3).contains(&bound),
            "psf {psf}: error / bound = {bound}"
        );
    }
}

#[test]
fn empty_frames_give_about_the_false_stars_asked_for() {
    let mut total = 0;
    for seed in 0..6 {
        let mut scene = Scene::new(640, 640);
        scene.gradient = (0.05, -0.03);
        let det = detect(&scene.render(seed), 640, 640, &DetectOptions::default()).unwrap();
        assert!(det.streaks.is_empty());
        for s in det.stars() {
            assert!(s.p_star < 0.9, "a noise peak with p_star {}", s.p_star);
            total += 1;
        }
    }
    assert!(total <= 14, "{total} false stars in 6 frames");
    let scene = Scene::new(640, 640);
    let img = scene.render(77);
    let many = DetectOptions {
        false_alarms: 100.0,
        ..DetectOptions::default()
    };
    let n = detect(&img, 640, 640, &many).unwrap().sources.len();
    assert!((30..300).contains(&n), "{n}");
    let none = DetectOptions {
        threshold: Some(7.0),
        ..DetectOptions::default()
    };
    assert!(detect(&img, 640, 640, &none).unwrap().sources.is_empty());
}

fn has_streak(det: &Detection, truth: (f64, f64, f64, f64), tol: f64) -> bool {
    det.streaks.iter().any(|t| {
        let d = |a: (f64, f64), b: (f64, f64)| (a.0 - b.0).phypot(a.1 - b.1);
        let (p, q) = ((t.x0, t.y0), (t.x1, t.y1));
        let (a, b) = ((truth.0, truth.1), (truth.2, truth.3));
        (d(p, a) < tol && d(q, b) < tol) || (d(p, b) < tol && d(q, a) < tol)
    })
}

#[test]
fn streaks_are_traced_end_to_end_and_kept_out_of_the_stars() {
    let mut scene = Scene::new(900, 700);
    let noise = scene.sky_noise();
    scatter(&mut scene, 150, (8.0, 3000.0), &mut Rng(4));
    let streaks = [
        (120.0, 80.0, 760.0, 560.0, 15.0),
        (860.0, -30.0, 520.0, 730.0, 45.0),
        (90.0, 600.0, 118.0, 612.0, 35.0),
        (300.0, 660.0, 650.0, 640.0, 6.0),
    ];
    for &(x0, y0, x1, y1, f) in &streaks {
        scene.streaks.push(Streak {
            x0,
            y0,
            x1,
            y1,
            flux_per_px: f * noise,
        });
    }
    scene.stars.push(Star {
        x: 120.0 - 11.0,
        y: 80.0 - 8.25,
        flux: flux_for(&scene, 1500.0),
    });
    let det = detect(&scene.render(31), 900, 700, &DetectOptions::default()).unwrap();
    assert_eq!(det.streaks.len(), 4, "{:#?}", det.streaks);
    assert!(
        has_streak(&det, (120.0, 80.0, 760.0, 560.0), 2.0),
        "{:#?}",
        det.streaks
    );
    assert!(has_streak(&det, (90.0, 600.0, 118.0, 612.0), 2.0));
    assert!(
        has_streak(&det, (300.0, 660.0, 650.0, 640.0), 10.0),
        "{:#?}",
        det.streaks
    );
    let long = det
        .streaks
        .iter()
        .find(|t| t.touches_edge)
        .expect("a streak at the edge");
    assert!(long.length > 750.0, "{long:?}");
    let slope = (long.y1 - long.y0) / (long.x1 - long.x0);
    assert!((slope - 760.0 / -340.0).abs() < 0.01, "{slope}");
    for t in &det.streaks {
        assert!((t.width / 1.53 - 1.0).abs() < 0.25, "width {}", t.width);
    }
    assert!(nearest(&det, 109.0, 71.75, 1.0).is_some());
    let stray = det
        .stars()
        .filter(|d| {
            let on_streak = streaks.iter().any(|&(x0, y0, x1, y1, _)| {
                let (dx, dy) = (x1 - x0, y1 - y0);
                let len = dx.phypot(dy);
                let along = ((d.x - x0) * dx + (d.y - y0) * dy) / len;
                let across = (-(d.x - x0) * dy + (d.y - y0) * dx) / len;
                across.abs() < 2.5 && along > -2.0 && along < len + 2.0
            });
            on_streak
                && !scene
                    .stars
                    .iter()
                    .any(|s| (s.x - d.x).phypot(s.y - d.y) < 3.0)
        })
        .count();
    assert_eq!(stray, 0);
    assert!(
        det.sources
            .iter()
            .filter(|d| d.class == Class::Streak)
            .count()
            > 200
    );
    let off = DetectOptions {
        streaks: false,
        ..DetectOptions::default()
    };
    let plain = detect(&scene.render(31), 900, 700, &off).unwrap();
    assert!(plain.streaks.is_empty());
    assert!(
        plain
            .sources
            .iter()
            .filter(|d| d.class == Class::Extended)
            .count()
            > 150
    );
}

#[test]
fn a_row_of_stars_is_not_a_streak() {
    let mut scene = Scene::new(400, 300);
    for i in 0..12 {
        let f = flux_for(&scene, 300.0);
        scene.stars.push(Star {
            x: 40.0 + 9.0 * i as f64,
            y: 100.0 + 4.5 * i as f64,
            flux: f,
        });
    }
    let det = detect(&scene.render(3), 400, 300, &DetectOptions::default()).unwrap();
    assert!(det.streaks.is_empty(), "{:#?}", det.streaks);
    assert_eq!(det.stars().count(), 12);
}

#[test]
fn hot_pixels_and_cosmic_rays_are_not_stars() {
    for psf in [0.7, 1.5, 3.0] {
        let mut scene = Scene::new(512, 512);
        scene.psf = Psf::Gaussian(psf);
        scatter(&mut scene, 80, (15.0, 800.0), &mut Rng(9));
        let noise = scene.sky_noise();
        let mut rng = Rng(10);
        while scene.hot.len() < 120 {
            let (x, y) = (
                (6.0 + rng.uniform() * 500.0) as usize,
                (6.0 + rng.uniform() * 500.0) as usize,
            );
            if scene
                .stars
                .iter()
                .all(|s| (s.x - x as f64).phypot(s.y - y as f64) > 12.0)
            {
                scene
                    .hot
                    .push((x, y, noise * (15.0 + 400.0 * rng.uniform())));
            }
        }
        let hot: Vec<(usize, usize, f64)> = scene.hot.clone();
        for i in 0..5 {
            scene.hot.push((40 + i, 300 + i / 2, 200.0 * noise));
        }
        let det = detect(&scene.render(11), 512, 512, &DetectOptions::default()).unwrap();
        let passed = hot
            .iter()
            .filter(|&&(x, y, _)| nearest(&det, x as f64, y as f64, 1.5).is_some())
            .count();
        assert!(
            passed <= 6,
            "psf {psf}: {passed} of 120 hot pixels passed as stars"
        );
        assert!(
            nearest(&det, 42.0, 301.0, 3.0).is_none(),
            "psf {psf}: cosmic ray passed"
        );
        let found = scene
            .stars
            .iter()
            .filter(|s| nearest(&det, s.x, s.y, 2.0).is_some())
            .count();
        assert!(found >= 78, "psf {psf}: {found} of 80 stars");
        let width = (psf * psf + 1.0 / 12.0).sqrt();
        assert!(
            (det.stats.psf_sigma / width - 1.0).abs() < 0.1,
            "{}",
            det.stats.psf_sigma
        );
    }
}

#[test]
fn undersampled_stars_survive_sharpness_checks() {
    let mut scene = Scene::new(512, 512);
    scene.psf = Psf::Gaussian(0.5 / 2.354_820_045);
    scatter(&mut scene, 100, (15.0, 3000.0), &mut Rng(9));
    for (x, y) in [(5, 80), (5, 180), (5, 280)] {
        scene.hot.push((x, y, 500.0 * scene.sky_noise()));
    }
    let det = detect(&scene.render(11), 512, 512, &DetectOptions::default()).unwrap();
    let found = scene
        .stars
        .iter()
        .filter(|s| nearest(&det, s.x, s.y, 1.0).is_some())
        .count();
    assert!(found >= 95, "{found} of 100 undersampled stars found");
    assert!(
        det.stars().count() <= found + 10,
        "{} detections without a matching star",
        det.stars().count() - found
    );
}

#[test]
fn stray_light_makes_no_false_stars() {
    let mut scene = Scene::new(768, 768);
    scene.stray = Some(Stray {
        x: -40.0,
        y: -60.0,
        peak: 20_000.0,
        half_px: 120.0,
    });
    scatter(&mut scene, 200, (10.0, 1000.0), &mut Rng(14));
    let det = detect(&scene.render(15), 768, 768, &DetectOptions::default()).unwrap();
    assert!(det.streaks.is_empty(), "{:#?}", det.streaks);
    let false_stars = det
        .stars()
        .filter(|d| {
            !scene
                .stars
                .iter()
                .any(|s| (s.x - d.x).phypot(s.y - d.y) < 3.0)
        })
        .count();
    assert!(
        false_stars <= 3,
        "{false_stars} false stars under stray light"
    );
    let model = scene.model();
    let (mut bright, mut found) = (0, 0);
    for st in &scene.stars {
        let i = (st.y as usize).min(759) * 768 + (st.x as usize + 8).min(767);
        let noise = (model[i] / scene.gain + scene.read_noise * scene.read_noise).sqrt();
        if st.flux > 12.0 * noise * 2.0 * std::f64::consts::PI.sqrt() * scene.psf_sigma() {
            bright += 1;
            found += usize::from(nearest(&det, st.x, st.y, 2.0).is_some());
        }
    }
    assert!(bright > 100 && found + 2 >= bright, "{found} of {bright}");
    assert!(det.stats.background > 1000.0);
}

#[test]
fn cosmic_rays_by_the_hundred_are_not_stars() {
    let mut scene = Scene::new(512, 512);
    scatter(&mut scene, 60, (30.0, 800.0), &mut Rng(16));
    let noise = scene.sky_noise();
    let mut rng = Rng(17);
    while scene.cosmics.len() < 900 {
        let (x, y) = (rng.uniform() * 512.0, rng.uniform() * 512.0);
        let length = match scene.cosmics.len() % 3 {
            0 => 0.0,
            1 => 1.0 + 3.0 * rng.uniform(),
            _ => 5.0 + 55.0 * rng.uniform() * rng.uniform(),
        };
        let angle = std::f64::consts::PI * rng.uniform();
        scene.cosmics.push(Cosmic {
            x,
            y,
            length,
            angle,
            adu_per_px: noise * (8.0 + 200.0 * rng.uniform()),
        });
    }
    let det = detect(&scene.render(18), 512, 512, &DetectOptions::default()).unwrap();
    assert!(
        (det.stats.psf_sigma / 1.528 - 1.0).abs() < 0.1,
        "{}",
        det.stats.psf_sigma
    );
    assert_eq!(det.stats.saturation, None);
    let on_hit = |d: &Detected| {
        scene.cosmics.iter().any(|c| {
            let (sn, cs) = c.angle.psin_cos();
            let (dx, dy) = (d.x - c.x, d.y - c.y);
            let (along, across) = (dx * cs + dy * sn, -dx * sn + dy * cs);
            across.abs() < 1.5 && along > -1.5 && along < c.length + 1.5
        })
    };
    let near_star = |d: &Detected| {
        scene
            .stars
            .iter()
            .any(|s| (s.x - d.x).phypot(s.y - d.y) < 3.0)
    };
    let passed = det.stars().filter(|d| on_hit(d) && !near_star(d)).count();
    assert!(passed <= 27, "{passed} of 900 hits passed as stars");
    let (mut clear, mut found) = (0, 0);
    for st in &scene.stars {
        let hit_nearby = scene.cosmics.iter().any(|c| {
            let (sn, cs) = c.angle.psin_cos();
            let (dx, dy) = (st.x - c.x, st.y - c.y);
            let (along, across) = (dx * cs + dy * sn, -dx * sn + dy * cs);
            across.abs() < 6.0 && along > -6.0 && along < c.length + 6.0
        });
        if !hit_nearby {
            clear += 1;
            found += usize::from(nearest(&det, st.x, st.y, 1.5).is_some());
        }
    }
    assert!(clear >= 20 && found == clear, "{found} of {clear}");
    assert!(
        det.streaks.iter().all(|t| t.cosmic_ray),
        "{:#?}",
        det.streaks
    );
}

#[test]
fn stars_smeared_by_a_drift_are_still_stars() {
    let mut scene = Scene::new(640, 640);
    scene.psf = Psf::Gaussian(0.8);
    let mut centres = Scene::new(640, 640);
    centres.psf = Psf::Gaussian(0.8);
    scatter(&mut centres, 150, (20.0, 3000.0), &mut Rng(23));
    let (sn, cs) = 0.5f64.psin_cos();
    for st in &centres.stars {
        for k in 0..9 {
            let t = 0.6 * (k as f64 - 4.0);
            scene.stars.push(Star {
                x: st.x + t * cs,
                y: st.y + t * sn,
                flux: st.flux / 9.0,
            });
        }
    }
    let noise = scene.sky_noise();
    let mut rng = Rng(24);
    while scene.hot.len() < 40 {
        let (x, y) = (
            (6.0 + rng.uniform() * 628.0) as usize,
            (6.0 + rng.uniform() * 628.0) as usize,
        );
        if centres
            .stars
            .iter()
            .all(|s| (s.x - x as f64).phypot(s.y - y as f64) > 12.0)
        {
            scene
                .hot
                .push((x, y, noise * (30.0 + 400.0 * rng.uniform())));
        }
    }
    let det = detect(&scene.render(25), 640, 640, &DetectOptions::default()).unwrap();
    assert!(det.streaks.is_empty(), "{:#?}", det.streaks);
    let mut lost = Vec::new();
    for st in &centres.stars {
        match nearest(&det, st.x, st.y, 2.0) {
            Some(d) => {
                let err = (d.x - st.x).phypot(d.y - st.y);
                assert!(err < (4.0 * d.sigma_px).max(0.35), "{d:?} for {st:?}");
            }
            None => lost.push(st.flux / flux_for(&centres, 1.0)),
        }
    }
    assert!(lost.is_empty(), "smeared stars lost, at SNR {lost:?}");
    let mut by_flux: Vec<&Star> = centres.stars.iter().collect();
    by_flux.sort_by(|a, b| b.flux.total_cmp(&a.flux));
    for st in &by_flux[..20] {
        let rank = det
            .stars()
            .position(|d| (d.x - st.x).phypot(d.y - st.y) < 2.0)
            .unwrap();
        assert!(rank < 30, "a top-20 star at rank {rank}");
    }
    let passed = scene
        .hot
        .iter()
        .filter(|&&(x, y, _)| nearest(&det, x as f64, y as f64, 1.5).is_some())
        .count();
    assert!(passed <= 2, "{passed} of 40 hot pixels passed as stars");
}

#[test]
fn saturated_stars_keep_their_centre_with_or_without_a_bleed_trail() {
    let (w, h) = (600usize, 600usize);
    let mut scene = Scene::new(w, h);
    scene.psf = Psf::Gaussian(1.8);
    scene.saturation = Some(60_000.0);
    scene.quantise = true;
    scatter(&mut scene, 60, (20.0, 500.0), &mut Rng(21));
    let bright: Vec<(f64, f64, f64)> = (0..8)
        .map(|i| {
            (
                70.3 + 65.0 * i as f64,
                90.6 + 55.0 * i as f64,
                3.0 * 2f64.ppowi(i),
            )
        })
        .collect();
    for &(x, y, over) in &bright {
        let flux = over * 59_000.0 * 2.0 * std::f64::consts::PI * 1.8 * 1.8;
        scene.stars.push(Star { x, y, flux });
    }
    let mut img = scene.render(8);
    for &(x, y, _) in &bright[5..] {
        let (cx, cy) = (x.round() as usize, y.round() as usize);
        for yy in cy.saturating_sub(150)..(cy + 60).min(h) {
            for xx in cx - 1..=cx + 1 {
                img[yy * w + xx] = 60_000.0;
            }
        }
    }
    for saturation in [None, Some(60_000.0f32)] {
        let opt = DetectOptions {
            saturation,
            ..DetectOptions::default()
        };
        let det = detect(&img, w, h, &opt).unwrap();
        assert!((det.stats.saturation.unwrap() / 60_000.0 - 1.0).abs() < 0.03);
        for &(x, y, over) in &bright {
            let d = nearest(&det, x, y, 6.0).unwrap_or_else(|| panic!("star ×{over} lost"));
            assert_eq!(d.class, Class::Saturated, "star ×{over}");
            let err = (d.x - x).phypot(d.y - y);
            assert!(err < 0.6, "star ×{over}: {err} px off");
            assert!(
                err < 5.0 * d.sigma_px,
                "star ×{over}: {err} px, claimed {}",
                d.sigma_px
            );
            let near = det
                .sources
                .iter()
                .filter(|s| (s.x - x).phypot(s.y - y) < 10.0)
                .count();
            assert_eq!(near, 1, "star ×{over}");
        }
        let ranks: Vec<usize> = bright
            .iter()
            .map(|&(x, y, _)| {
                det.sources
                    .iter()
                    .position(|s| (s.x - x).phypot(s.y - y) < 6.0)
                    .unwrap()
            })
            .collect();
        assert!(ranks.iter().all(|r| *r < 8), "{ranks:?}");
        assert!(ranks.windows(2).all(|p| p[0] > p[1]), "{ranks:?}");
        let on_trails = det
            .stars()
            .filter(|s| {
                bright[5..].iter().any(|&(x, y, _)| {
                    (s.x - x).abs() < 4.0
                        && (s.y - y).abs() > 12.0
                        && s.y > y - 155.0
                        && s.y < y + 65.0
                })
            })
            .count();
        assert_eq!(on_trails, 0);
    }
}

#[test]
fn wide_stars_are_binned_and_a_given_width_is_used() {
    let mut scene = Scene::new(900, 900);
    scene.psf = Psf::Moffat {
        fwhm: 9.0,
        beta: 3.0,
    };
    scatter(&mut scene, 40, (30.0, 600.0), &mut Rng(5));
    let img = scene.render(6);
    let det = detect(&img, 900, 900, &DetectOptions::default()).unwrap();
    assert!(
        det.stats.bin >= 2 && det.stats.passes >= 2,
        "{:?}",
        det.stats
    );
    assert!(det.stats.psf_measured);
    assert!(
        (det.stats.psf_sigma * 2.3548 / 9.0 - 1.0).abs() < 0.25,
        "{}",
        det.stats.psf_sigma
    );
    let found = scene
        .stars
        .iter()
        .filter(|s| nearest(&det, s.x, s.y, 3.0).is_some())
        .count();
    assert!(found >= 39, "{found} of 40");
    let told = DetectOptions {
        psf_sigma: Some(9.0 / 2.3548),
        ..DetectOptions::default()
    };
    let det2 = detect(&img, 900, 900, &told).unwrap();
    assert!(
        det2.stats.passes == 1 && det2.stats.bin >= 2,
        "{:?}",
        det2.stats
    );
    assert!(!det2.stats.psf_measured);
    assert!(scene
        .stars
        .iter()
        .all(|s| nearest(&det2, s.x, s.y, 3.0).is_some()));
}

#[test]
fn undefined_pixels_and_borders_are_handled() {
    let mut scene = Scene::new(400, 400);
    scatter(&mut scene, 40, (30.0, 500.0), &mut Rng(2));
    let f = flux_for(&scene, 300.0);
    scene.stars.push(Star {
        x: 150.0,
        y: 150.0,
        flux: f,
    });
    scene.stars.push(Star {
        x: 3.0,
        y: 200.0,
        flux: f,
    });
    let mut img = scene.render(4);
    for y in 130..170 {
        for x in 130..170 {
            img[y * 400 + x] = f32::NAN;
        }
    }
    img[20] = f32::INFINITY;
    let det = detect(&img, 400, 400, &DetectOptions::default()).unwrap();
    assert_eq!(det.stats.undefined, 1601);
    assert!(det
        .sources
        .iter()
        .all(|s| !(128.0..172.0).contains(&s.x) || !(128.0..172.0).contains(&s.y)));
    let edge = nearest(&det, 3.0, 200.0, 2.0).expect("the star at the edge");
    assert!(edge.truncated);
    let inner = det.stars().filter(|s| !s.truncated).count();
    assert!(inner >= 38, "{inner}");
    assert!(det
        .sources
        .iter()
        .all(|s| s.x.is_finite() && s.y.is_finite() && s.sigma_px > 0.0 && s.flux.is_finite()));

    assert!(detect(&img, 400, 399, &DetectOptions::default()).is_err());
    assert!(detect(&[0.0; 100], 10, 10, &DetectOptions::default()).is_err());
    assert!(detect(&vec![f32::NAN; 1024], 32, 32, &DetectOptions::default()).is_err());
    let bad = DetectOptions {
        false_alarms: 0.0,
        ..DetectOptions::default()
    };
    assert!(detect(&img, 400, 400, &bad).is_err());
    let flat = detect(&vec![100.0; 64 * 64], 64, 64, &DetectOptions::default()).unwrap();
    assert!(flat.sources.is_empty());
}

#[test]
fn hostile_settings_and_frames_end_quickly() {
    let started = std::time::Instant::now();
    for (w, h) in [(16usize, 16usize), (20, 40), (33, 17)] {
        let scene = Scene::new(w, h);
        let wide = DetectOptions {
            psf_sigma: Some(21.0),
            ..DetectOptions::default()
        };
        let det = detect(&scene.render(1), w, h, &wide).unwrap();
        assert_eq!(det.stats.bin, 1, "{w} × {h}");
    }
    let mut scene = Scene::new(512, 512);
    scatter(&mut scene, 60, (20.0, 2000.0), &mut Rng(31));
    let img = scene.render(32);
    let plain = detect(&img, 512, 512, &DetectOptions::default()).unwrap();
    let low = DetectOptions {
        saturation: Some(scene.sky as f32),
        ..DetectOptions::default()
    };
    let det = detect(&img, 512, 512, &low).unwrap();
    assert_eq!(det.stats.saturation, None);
    assert_eq!(det.stars().count(), plain.stars().count());
    let most = DetectOptions {
        saturation: Some((scene.sky + 25.0 * scene.sky_noise()) as f32),
        ..DetectOptions::default()
    };
    let det = detect(&img, 512, 512, &most).unwrap();
    assert!(det.stars().filter(|d| d.class == Class::Saturated).count() >= 20);
    let mut blocked = img.clone();
    for y in 100..400 {
        for x in 150..420 {
            blocked[y * 512 + x] = 60_000.0;
        }
    }
    for i in 0..512 {
        blocked[i * 512 + i] = 60_000.0;
    }
    let sat = DetectOptions {
        saturation: Some(60_000.0),
        ..DetectOptions::default()
    };
    let det = detect(&blocked, 512, 512, &sat).unwrap();
    assert!(det.stars().filter(|d| d.class == Class::Saturated).count() <= 2);
    let flat = vec![60_000f32; 700 * 700];
    assert!(detect(&flat, 700, 700, &sat).unwrap().sources.is_empty());
    let too_low = DetectOptions {
        threshold: Some(0.2),
        ..DetectOptions::default()
    };
    assert!(detect(&Scene::new(512, 512).render(5), 512, 512, &too_low).is_err());
    assert!(started.elapsed().as_secs() < 60, "{:?}", started.elapsed());
}

#[test]
fn the_list_is_brightest_first_and_can_be_cut() {
    let mut scene = Scene::new(500, 500);
    scatter(&mut scene, 120, (10.0, 2000.0), &mut Rng(8));
    let img = scene.render(1);
    let all = detect(&img, 500, 500, &DetectOptions::default()).unwrap();
    assert!(all.sources.windows(2).all(|p| p[0].flux >= p[1].flux));
    let cut = DetectOptions {
        max_sources: 25,
        ..DetectOptions::default()
    };
    let few = detect(&img, 500, 500, &cut).unwrap();
    assert_eq!(few.stars().count(), 25);
    assert_eq!(few.sources[..], all.sources[..few.sources.len()]);
    for st in scene
        .stars
        .iter()
        .filter(|s| s.flux > flux_for(&scene, 100.0))
    {
        let d = nearest(&all, st.x, st.y, 2.0).unwrap();
        assert!(
            (d.flux / st.flux - 1.0).abs() < 0.1,
            "{} vs {}",
            d.flux,
            st.flux
        );
    }
}

fn fingerprint(det: &Detection) -> u64 {
    let mut h = 0xcbf2_9ce4_8422_2325u64;
    let mut put = |v: f64| {
        for b in v.to_bits().to_le_bytes() {
            h = (h ^ u64::from(b)).wrapping_mul(0x0100_0000_01b3);
        }
    };
    for s in &det.sources {
        for v in [
            s.x,
            s.y,
            s.sigma_px,
            s.flux,
            s.snr,
            s.p_star,
            s.sigma_minor,
            s.theta,
        ] {
            put(v);
        }
    }
    for t in &det.streaks {
        for v in [
            t.x0, t.y0, t.x1, t.y1, t.length, t.angle, t.width, t.flux, t.snr,
        ] {
            put(v);
        }
    }
    put(det.stats.noise);
    put(det.stats.psf_sigma);
    h
}

#[test]
fn results_are_pinned_to_the_bit() {
    let mut scene = Scene::new(512, 384);
    scene.gradient = (0.04, 0.02);
    scene.saturation = Some(40_000.0);
    scene.quantise = true;
    scatter(&mut scene, 90, (6.0, 4000.0), &mut Rng(42));
    let noise = scene.sky_noise();
    scene.streaks.push(Streak {
        x0: 30.0,
        y0: 300.0,
        x1: 470.0,
        y1: 60.0,
        flux_per_px: 20.0 * noise,
    });
    scene.hot.push((100, 100, 3000.0));
    let img = scene.render(43);
    let a = detect(&img, 512, 384, &DetectOptions::default()).unwrap();
    assert_eq!(
        a,
        detect(&img, 512, 384, &DetectOptions::default()).unwrap()
    );
    assert_eq!(a.streaks.len(), 1);
    let got = fingerprint(&a);
    assert_eq!(
        got,
        PINNED,
        "fingerprint {got:#018x} ({} sources)",
        a.sources.len()
    );
}

const PINNED: u64 = 0xb2d2_5006_567b_e21e;

#[test]
fn bias_subtracted_points_with_a_clipped_sky() {
    let mut scene = Scene::new(512, 512);
    scatter(&mut scene, 100, (12.0, 300.0), &mut Rng(29));
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
    let det = detect(&image, 512, 512, &opt).unwrap();
    let stars: Vec<_> = det.stars().collect();
    let matches = stars
        .iter()
        .filter(|d| {
            scene
                .stars
                .iter()
                .any(|s| (d.x - s.x).phypot(d.y - s.y) < 2.0)
        })
        .count();
    let recovered = scene
        .stars
        .iter()
        .filter(|s| nearest(&det, s.x, s.y, 2.0).is_some())
        .count();
    assert!(
        matches * 10 >= 9 * stars.len(),
        "{matches}/{} matches",
        stars.len()
    );
    assert!(recovered >= 70, "{recovered}/100 recovered");
}
