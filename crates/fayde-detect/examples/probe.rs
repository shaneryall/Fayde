use std::time::Instant;

use fayde_core::Portable;
use fayde_detect::synth::{Cosmic, Psf, Rng, Scene, Star, Stray, Streak};
use fayde_detect::{detect, Class, DetectOptions};

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let arg = |i: usize, d: f64| args.get(i).and_then(|v| v.parse().ok()).unwrap_or(d);
    let size = arg(0, 2048.0) as usize;
    let psf = arg(1, 1.5);
    let n_stars = arg(2, 400.0) as usize;
    let seed = arg(3, 1.0) as u64;
    let mode = args.get(4).map(String::as_str).unwrap_or("full");

    let mut rng = Rng(seed);
    let mut sc = Scene::new(size, size);
    sc.psf = if mode == "moffat" {
        Psf::Moffat {
            fwhm: psf * 2.3548,
            beta: 3.0,
        }
    } else {
        Psf::Gaussian(psf)
    };
    sc.gradient = (0.02, -0.01);
    let noise = sc.sky_noise();
    for _ in 0..n_stars {
        let snr = 3.0 * 1000f64.ppowf(rng.uniform());
        let flux = snr * noise * 2.0 * std::f64::consts::PI.sqrt() * psf;
        sc.stars.push(Star {
            x: 8.0 + rng.uniform() * (size as f64 - 16.0),
            y: 8.0 + rng.uniform() * (size as f64 - 16.0),
            flux,
        });
    }
    if mode != "stars" && mode != "noise" {
        let s = size as f64;
        sc.streaks.push(Streak {
            x0: 0.1 * s,
            y0: 0.2 * s,
            x1: 0.7 * s,
            y1: 0.55 * s,
            flux_per_px: 12.0 * noise,
        });
        sc.streaks.push(Streak {
            x0: 0.8 * s,
            y0: -20.0,
            x1: 0.55 * s,
            y1: s + 20.0,
            flux_per_px: 40.0 * noise,
        });
        sc.streaks.push(Streak {
            x0: 0.3 * s,
            y0: 0.8 * s,
            x1: 0.3 * s + 25.0,
            y1: 0.8 * s + 9.0,
            flux_per_px: 30.0 * noise,
        });
        for _ in 0..60 {
            sc.hot.push((
                (rng.uniform() * size as f64) as usize,
                (rng.uniform() * size as f64) as usize,
                noise * (10.0 + 300.0 * rng.uniform()),
            ));
        }
    }
    if mode == "noise" {
        sc.stars.clear();
    }
    if let Ok(v) = std::env::var("STRAY") {
        let p: Vec<f64> = v.split(',').filter_map(|x| x.parse().ok()).collect();
        sc.stray = Some(Stray {
            x: -40.0,
            y: -60.0,
            peak: p[0],
            half_px: p[1],
        });
    }
    if let Ok(v) = std::env::var("COSMICS") {
        let n: usize = v.parse().unwrap_or(0);
        for i in 0..n {
            let length = match i % 3 {
                0 => 0.0,
                1 => 1.0 + 3.0 * rng.uniform(),
                _ => 5.0 + 60.0 * rng.uniform() * rng.uniform(),
            };
            sc.cosmics.push(Cosmic {
                x: rng.uniform() * size as f64,
                y: rng.uniform() * size as f64,
                length,
                angle: std::f64::consts::PI * rng.uniform(),
                adu_per_px: noise * (8.0 + 200.0 * rng.uniform()),
            });
        }
    }
    if let Ok(v) = std::env::var("SAT") {
        sc.saturation = v.parse().ok();
        sc.quantise = true;
    }
    let img = sc.render(seed + 100);

    let mut opt = DetectOptions::default();
    if let Ok(t) = std::env::var("THRESH") {
        opt.threshold = t.parse().ok();
    }
    if let Ok(t) = std::env::var("PSF") {
        opt.psf_sigma = t.parse().ok();
    }
    if std::env::var("GAIN").is_ok() {
        opt.gain = Some(sc.gain);
        opt.sigma_floor = 0.0;
    }
    let t0 = Instant::now();
    let det = detect(&img, size, size, &opt).unwrap();
    let dt = t0.elapsed().as_secs_f64();
    let mut best = dt;
    for _ in 0..4 {
        let t = Instant::now();
        let _ = detect(&img, size, size, &opt).unwrap();
        best = best.min(t.elapsed().as_secs_f64());
    }
    println!(
        "{size}×{size}  detect {:.1} ms ({:.2} ms/Mpx)  stats {:?}",
        best * 1e3,
        best * 1e3 / (size * size) as f64 * 1e6,
        det.stats
    );
    let count = |c: Class| det.sources.iter().filter(|s| s.class == c).count();
    println!(
        "sources {}: star {} saturated {} extended {} sharp {} streak {}",
        det.sources.len(),
        count(Class::Star),
        count(Class::Saturated),
        count(Class::Extended),
        count(Class::Sharp),
        count(Class::Streak)
    );

    let stellar: Vec<_> = det.stars().collect();
    let mut used = vec![false; stellar.len()];
    let mut bins = [(0usize, 0usize); 6];
    let bin_of = |snr: f64| match snr {
        s if s < 5.0 => 0,
        s if s < 8.0 => 1,
        s if s < 15.0 => 2,
        s if s < 50.0 => 3,
        s if s < 500.0 => 4,
        _ => 5,
    };
    let mut chi = Vec::new();
    let mut size_z = vec![Vec::new(); 6];
    let mut sat_err = Vec::new();
    let mut sat_chi = Vec::new();
    let mut err_by_bin = vec![Vec::new(); 6];
    let mut sig_by_bin = vec![Vec::new(); 6];
    for st in &sc.stars {
        let snr = st.flux / (noise * 2.0 * std::f64::consts::PI.sqrt() * psf);
        let b = bin_of(snr);
        bins[b].0 += 1;
        let mut best: Option<(usize, f64)> = None;
        for (i, d) in stellar.iter().enumerate() {
            let r = (d.x - st.x).phypot(d.y - st.y);
            if r < 2.0 + 1.5 * psf && !used[i] && best.is_none_or(|(_, br)| r < br) {
                best = Some((i, r));
            }
        }
        if let Some((i, _)) = best {
            used[i] = true;
            bins[b].1 += 1;
            let d = stellar[i];
            if d.class == Class::Saturated {
                sat_err.push(d.x - st.x);
                sat_err.push(d.y - st.y);
                sat_chi.push((d.x - st.x) / d.sigma_px);
                sat_chi.push((d.y - st.y) / d.sigma_px);
                continue;
            }
            let k2 = det.stats.kernel_sigma * det.stats.kernel_sigma;
            let truth = psf * psf + 1.0 / 12.0;
            size_z[b].push((d.sigma_minor * d.sigma_minor - truth) / (truth + k2) * d.snr);
            chi.push((d.x - st.x) / d.sigma_px);
            chi.push((d.y - st.y) / d.sigma_px);
            err_by_bin[b].push(d.x - st.x);
            err_by_bin[b].push(d.y - st.y);
            sig_by_bin[b].push(d.sigma_px);
        }
    }
    let rms = |v: &[f64]| (v.iter().map(|x| x * x).sum::<f64>() / v.len().max(1) as f64).sqrt();
    let robust = |v: &[f64]| {
        let mut a: Vec<f64> = v.iter().map(|x| x.abs()).collect();
        a.sort_by(f64::total_cmp);
        if a.is_empty() {
            0.0
        } else {
            1.4826 * a[a.len() / 2]
        }
    };
    let labels = ["<5", "5-8", "8-15", "15-50", "50-500", ">500"];
    for (i, (n, f)) in bins.iter().enumerate() {
        let mean_sig = sig_by_bin[i].iter().sum::<f64>() / sig_by_bin[i].len().max(1) as f64;
        println!(
            "  SNR {:>7}: {f:4}/{n:4} found   rms err {:.4} px (robust {:.4})   mean claimed {:.4} px",
            labels[i],
            rms(&err_by_bin[i]),
            robust(&err_by_bin[i]),
            mean_sig
        );
    }
    println!(
        "  error / claimed sigma: rms {:.3}, robust {:.3}",
        rms(&chi),
        robust(&chi)
    );
    let med = |v: &[f64]| {
        let mut a = v.to_vec();
        a.sort_by(f64::total_cmp);
        if a.is_empty() {
            0.0
        } else {
            a[a.len() / 2]
        }
    };
    println!(
        "  minor-axis variance error × SNR / total, by SNR bin: spread {:?}, median {:?}",
        size_z
            .iter()
            .map(
                |v| (robust(&v.iter().map(|x| x - med(v)).collect::<Vec<_>>()) * 100.0).round()
                    / 100.0
            )
            .collect::<Vec<_>>(),
        size_z
            .iter()
            .map(|v| (med(v) * 100.0).round() / 100.0)
            .collect::<Vec<_>>()
    );
    if !sat_err.is_empty() {
        println!(
            "  saturated: {} stars, rms err {:.3} px (robust {:.3}), error / claimed: rms {:.2}, robust {:.2}",
            sat_err.len() / 2,
            rms(&sat_err),
            robust(&sat_err),
            rms(&sat_chi),
            robust(&sat_chi)
        );
    }
    let false_stellar: Vec<_> = stellar
        .iter()
        .zip(&used)
        .filter(|(_, u)| !**u)
        .map(|(d, _)| *d)
        .collect();
    println!("  unmatched stellar detections: {}", false_stellar.len());
    if !sc.cosmics.is_empty() {
        let hit = |d: &fayde_detect::Detected| {
            sc.cosmics.iter().any(|c| {
                let (sn, cs) = c.angle.psin_cos();
                let (dx, dy) = (d.x - c.x, d.y - c.y);
                let along = dx * cs + dy * sn;
                let across = -dx * sn + dy * cs;
                across.abs() < 1.5 && along > -1.5 && along < c.length + 1.5
            })
        };
        println!(
            "  on a cosmic-ray hit: {} of the unmatched; streaks narrower than half the PSF: {} of {}",
            false_stellar.iter().filter(|d| hit(d)).count(),
            det.streaks.iter().filter(|t| t.width < 0.5 * det.stats.psf_sigma).count(),
            det.streaks.len()
        );
    }
    for d in false_stellar.iter().take(12) {
        println!(
            "    ({:7.1}, {:7.1}) snr {:6.1} p_star {:.2} {:?} maj {:.2} min {:.2}",
            d.x, d.y, d.snr, d.p_star, d.class, d.sigma_major, d.sigma_minor
        );
    }
    for t in &det.streaks {
        println!(
            "  streak ({:.1}, {:.1}) → ({:.1}, {:.1}) len {:.1} angle {:.2}° width {:.2} snr {:.0} cover {:.2} edge {}",
            t.x0, t.y0, t.x1, t.y1, t.length, t.angle.to_degrees(), t.width, t.snr, t.coverage, t.touches_edge
        );
    }
    for t in &sc.streaks {
        println!(
            "  true streak ({:.1}, {:.1}) → ({:.1}, {:.1})",
            t.x0, t.y0, t.x1, t.y1
        );
        let (dx, dy) = (t.x1 - t.x0, t.y1 - t.y0);
        let len = dx.phypot(dy);
        for st in &sc.stars {
            let along = ((st.x - t.x0) * dx + (st.y - t.y0) * dy) / len;
            let across = (-(st.x - t.x0) * dy + (st.y - t.y0) * dx) / len;
            if across.abs() < 4.0
                && (along > -60.0 && along < 0.0 || along > len && along < len + 60.0)
            {
                let snr = st.flux / (noise * 2.0 * std::f64::consts::PI.sqrt() * psf);
                println!(
                    "    star near its line: along {along:.1} across {across:.1} snr {snr:.0}"
                );
            }
        }
    }
}
