use fayde_core::Portable;
use fayde_detect::synth::{Cosmic, Psf, Rng, Scene, Star};
use fayde_detect::{detect, DetectOptions, Detection, FrameMode};

const SIDE: usize = 640;
const STEPS: usize = 48;

fn scene(path: impl Fn(f64, f64, f64) -> (f64, f64)) -> (Scene, Vec<Star>) {
    let mut scene = Scene::new(SIDE, SIDE);
    scene.sky = 40.0;
    scene.gain = 2.0;
    scene.psf = Psf::Gaussian(0.9);
    scene.gradient = (0.013, 0.025);
    let mut rng = Rng(430);
    let mut stars = Vec::new();
    for row in 0..9 {
        for col in 0..9 {
            let x = 55.0 + 65.0 * col as f64 + 10.0 * rng.uniform();
            let y = 52.0 + 65.0 * row as f64 + 10.0 * rng.uniform();
            let flux = 600.0 * 30.0_f64.ppowf(rng.uniform());
            let mut sum = (0.0, 0.0);
            for i in 0..STEPS {
                let t = (i as f64 + 0.5) / STEPS as f64;
                let (dx, dy) = path(x, y, t);
                sum.0 += dx;
                sum.1 += dy;
                scene.stars.push(Star {
                    x: x + dx,
                    y: y + dy,
                    flux: flux / STEPS as f64,
                });
            }
            stars.push(Star {
                x: x + sum.0 / STEPS as f64,
                y: y + sum.1 / STEPS as f64,
                flux,
            });
        }
    }
    (scene, stars)
}

fn check(detection: &Detection, truth: &[Star]) {
    let sources: Vec<_> = detection.stars().collect();
    let distance = |x: f64, y: f64, s: &Star| (x - s.x).phypot(y - s.y);
    let mut errors: Vec<_> = truth
        .iter()
        .filter_map(|s| {
            sources
                .iter()
                .map(|d| distance(d.x, d.y, s))
                .filter(|&e| e <= 3.0)
                .min_by(f64::total_cmp)
        })
        .collect();
    errors.sort_by(f64::total_cmp);
    let purity = sources
        .iter()
        .filter(|d| truth.iter().any(|s| distance(d.x, d.y, s) <= 3.0))
        .count();
    assert!(
        errors.len() >= truth.len() * 9 / 10,
        "found {}",
        errors.len()
    );
    assert!(
        purity * 10 >= sources.len() * 9,
        "{purity}/{}",
        sources.len()
    );
    let p90 = errors[errors.len() * 9 / 10];
    assert!(p90 < 0.8, "{p90} px");
    for d in sources {
        assert!(d.flux > 0.0 && d.sigma_px.is_finite() && d.sigma_px >= 0.05);
        assert!((0.0..=1.0).contains(&d.p_star));
    }
}

#[test]
fn shared_right_angle_path_with_defects() {
    let (mut scene, truth) = scene(|_, _, t| {
        let (u, v) = (40.0 * t.min(0.5), 40.0 * (t - 0.5).max(0.0));
        let (sin, cos) = 0.48_f64.psin_cos();
        (cos * u - sin * v, sin * u + cos * v)
    });
    let mut rng = Rng(159);
    for _ in 0..35 {
        scene.hot.push((
            (SIDE as f64 * rng.uniform()) as usize,
            (SIDE as f64 * rng.uniform()) as usize,
            800.0,
        ));
        scene.cosmics.push(Cosmic {
            x: SIDE as f64 * rng.uniform(),
            y: SIDE as f64 * rng.uniform(),
            length: 1.0 + 5.0 * rng.uniform(),
            angle: 3.0 * rng.uniform(),
            adu_per_px: 600.0,
        });
    }
    let image = scene.render(99);
    let options = DetectOptions {
        gain: Some(scene.gain),
        ..DetectOptions::default()
    };
    let detection = detect(&image, SIDE, SIDE, &options).unwrap();
    assert_eq!(
        detection.census.mode,
        FrameMode::Kernel,
        "{:?}",
        detection.census
    );
    let k = detection.kernel.unwrap();
    assert_eq!(k.regions, 1);
    assert!((5..=64).contains(&k.stacked));
    assert!((34.0..46.0).contains(&k.path_length), "{k:?}");
    assert!(k.bend > 3.0, "{k:?}");
    assert!(detection.trail.is_none());
    check(&detection, &truth);

    let points = DetectOptions {
        trails: false,
        ..options
    };
    let detection = detect(&image, SIDE, SIDE, &points).unwrap();
    assert_eq!(detection.census.mode, FrameMode::Points);
    assert!(detection.kernel.is_none());
}

#[test]
fn rotation_changes_the_template_across_the_field() {
    let (scene, truth) = scene(|x, y, t| {
        let (x, y) = (x - 328.0, y - 301.0);
        let (sin, cos) = (0.067 * (t - 0.5)).psin_cos();
        (x * (cos - 1.0) - y * sin, y * (cos - 1.0) + x * sin)
    });
    let image = scene.render(190);
    let detection = detect(&image, SIDE, SIDE, &DetectOptions::default()).unwrap();
    assert_eq!(
        detection.census.mode,
        FrameMode::Kernel,
        "{:?}",
        detection.census
    );
    let k = detection.kernel.unwrap();
    assert!(k.regions > 1 && k.stacked >= 5, "{k:?}");
    check(&detection, &truth);
}
