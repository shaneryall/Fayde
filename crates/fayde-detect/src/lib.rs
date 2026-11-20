//! Star and streak detection in astronomical images. [`detect`] accepts a row-major pixel array.

#![forbid(unsafe_code)]

mod filter;
mod kernel;
mod mesh;
pub mod pattern;
mod peaks;
mod saturate;
mod stats;
mod streak;
pub mod synth;
mod trail;

use std::borrow::Cow;

use fayde_core::Portable;
use filter::Kernel;
use mesh::{ClippedBackground, Mesh};
use peaks::Peak;

pub use kernel::FrameKernel;
pub use streak::Streak;
pub use trail::{Census, FrameMode, FrameTrail};

pub const MIN_SIDE: usize = 16;

#[derive(Debug, Clone, PartialEq)]
pub struct DetectOptions {
    /// Peak threshold in units of local filtered noise. Defaults to the threshold for [`Self::false_alarms`].
    pub threshold: Option<f64>,
    pub false_alarms: f64,
    pub psf_sigma: Option<f64>,
    pub saturation: Option<f32>,
    /// Maximum number of stellar detections to retain. Zero removes the cap.
    pub max_sources: usize,
    pub mesh: usize,
    pub streaks: bool,
    pub sigma_floor: f64,
    pub gain: Option<f64>,
    pub trails: bool,
    pub trail_angle: Option<f64>,
}

impl Default for DetectOptions {
    fn default() -> Self {
        Self {
            threshold: None,
            false_alarms: 1.0,
            psf_sigma: None,
            saturation: None,
            max_sources: 0,
            mesh: 64,
            streaks: true,
            sigma_floor: 0.05,
            gain: None,
            trails: true,
            trail_angle: None,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Class {
    Star,
    Saturated,
    Extended,
    Sharp,
    Streak,
    Point,
}

impl Class {
    pub fn name(self) -> &'static str {
        match self {
            Self::Star => "star",
            Self::Saturated => "saturated",
            Self::Extended => "extended",
            Self::Sharp => "sharp",
            Self::Streak => "streak",
            Self::Point => "point",
        }
    }

    pub fn is_stellar(self) -> bool {
        matches!(self, Self::Star | Self::Saturated | Self::Extended)
    }
}

/// Source position and shape in the original image, with zero-based pixel centres.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Detected {
    pub x: f64,
    pub y: f64,
    /// One-sigma centroid uncertainty in pixels, using the least precise axis for a trailed source.
    pub sigma_px: f64,
    pub sigma_along: f64,
    pub sigma_across: f64,
    pub flux: f64,
    pub snr: f64,
    pub p_star: f64,
    pub sigma_major: f64,
    pub sigma_minor: f64,
    pub theta: f64,
    pub class: Class,
    pub streak: Option<usize>,
    pub truncated: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub struct FrameStats {
    pub background: f64,
    pub noise: f64,
    pub psf_sigma: f64,
    pub psf_measured: bool,
    pub bin: usize,
    pub kernel_sigma: f64,
    pub threshold: f64,
    pub passes: usize,
    pub candidates: usize,
    pub saturation: Option<f64>,
    pub undefined: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Detection {
    pub sources: Vec<Detected>,
    pub streaks: Vec<Streak>,
    pub stats: FrameStats,
    pub census: Census,
    pub trail: Option<FrameTrail>,
    pub kernel: Option<FrameKernel>,
}

impl Detection {
    /// Iterate detections eligible for plate solving, including saturated and extended sources.
    pub fn stars(&self) -> impl Iterator<Item = &Detected> {
        self.sources.iter().filter(|s| s.class.is_stellar())
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DetectError {
    BadShape(String),
    BadOption(String),
    Empty,
}

impl std::fmt::Display for DetectError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::BadShape(s) | Self::BadOption(s) => write!(f, "{s}"),
            Self::Empty => write!(f, "the image has no defined pixel"),
        }
    }
}

impl std::error::Error for DetectError {}

pub const MIN_THRESHOLD: f64 = 1.5;

const MAX_PEAKS: usize = 200_000;

const FIRST_SIGMA: f64 = 1.2;
const MIN_EFFICIENCY: f64 = 0.93;
const KERNEL_RANGE: (f64, f64) = (0.6, 2.0);
const BIN_ABOVE: f64 = 2.0;

struct Pass<'a> {
    bin: usize,
    kernel: Kernel,
    w: usize,
    h: usize,
    raw: Cow<'a, [f32]>,
    s: Vec<f32>,
    mesh: Mesh,
    clipped: Option<&'a ClippedBackground>,
    peaks: Vec<Peak>,
    threshold: f64,
}

/// Detect stars and streaks in a row-major image. Non-finite pixels are masked.
pub fn detect(
    data: &[f32],
    width: usize,
    height: usize,
    opt: &DetectOptions,
) -> Result<Detection, DetectError> {
    if width.checked_mul(height) != Some(data.len()) {
        return Err(DetectError::BadShape(format!(
            "{} pixels given for a {width} × {height} image",
            data.len()
        )));
    }
    if width < MIN_SIDE || height < MIN_SIDE {
        return Err(DetectError::BadShape(format!(
            "image is {width} × {height}: each side must be at least {MIN_SIDE} pixels"
        )));
    }
    let bad = |name: &str| Err(DetectError::BadOption(format!("{name} is out of range")));
    if opt
        .threshold
        .is_some_and(|t| !(t.is_finite() && t >= MIN_THRESHOLD))
    {
        return Err(DetectError::BadOption(format!(
            "threshold must be at least {MIN_THRESHOLD} (σ of the filtered image)"
        )));
    }
    if !(opt.false_alarms.is_finite() && opt.false_alarms > 0.0) {
        return bad("false_alarms");
    }
    if opt.psf_sigma.is_some_and(|s| !(s.is_finite() && s > 0.0)) {
        return bad("psf_sigma");
    }
    if !(opt.sigma_floor.is_finite() && opt.sigma_floor >= 0.0) {
        return bad("sigma_floor");
    }
    if opt.mesh < 8 {
        return bad("mesh");
    }
    if opt.gain.is_some_and(|g| !(g.is_finite() && g > 0.0)) {
        return bad("gain");
    }
    if opt.trail_angle.is_some_and(|a| !a.is_finite()) {
        return bad("trail_angle");
    }

    let clipped = ClippedBackground::estimate(data, width, height, opt.mesh);
    let valid = clipped
        .as_ref()
        .and_then(|sky| sky.valid_pixels(data, width))
        .or_else(|| {
            data.iter()
                .any(|v| !v.is_finite())
                .then(|| data.iter().map(|v| v.is_finite()).collect())
        });
    let undefined = valid
        .as_ref()
        .map_or(0, |mask| mask.iter().filter(|&&ok| !ok).count());
    let image = if let Some(mask) = &valid {
        if undefined == data.len() {
            return Err(DetectError::Empty);
        }
        let step = (data.len() / 4096).max(1);
        let mut sample: Vec<f64> = data
            .iter()
            .zip(mask)
            .step_by(step)
            .filter_map(|(&v, &ok)| ok.then_some(f64::from(v)))
            .collect();
        let fill = match stats::median(&mut sample) {
            Some(m) => m as f32,
            None => data
                .iter()
                .zip(mask)
                .find_map(|(&v, &ok)| ok.then_some(v))
                .expect("a defined pixel"),
        };
        Cow::Owned(
            data.iter()
                .zip(mask)
                .map(|(&v, &ok)| if ok { v } else { fill })
                .collect::<Vec<_>>(),
        )
    } else {
        Cow::Borrowed(data)
    };

    let max_bin = (width.min(height) / MIN_SIDE).max(1);
    let (mut bin, mut sigma_k) = match opt.psf_sigma {
        Some(s) => plan_for(s, 1, max_bin),
        None => (1, FIRST_SIGMA),
    };
    let run = |bin, sigma| {
        run_pass(
            &image,
            valid.as_deref(),
            (width, height),
            bin,
            sigma,
            opt,
            clipped.as_ref(),
        )
    };
    let mut pass = run(bin, sigma_k);
    let mut passes = 1;
    let mut first: Option<Pass> = None;
    let mut psf = opt.psf_sigma.map(|s| scaled_psf(s, pass.bin));
    let mut measured = false;
    if opt.psf_sigma.is_none() {
        while passes < 3 {
            let Some(v) = measure_psf(&pass) else { break };
            psf = Some(v);
            measured = true;
            let (b, sk) = plan_for(v.sqrt() * pass.bin as f64, pass.bin, max_bin);
            let ks = pass.kernel.sigma;
            let eff = 2.0 * sk * ks / (sk * sk + ks * ks);
            if b == pass.bin && eff >= MIN_EFFICIENCY {
                break;
            }
            if width / b < MIN_SIDE || height / b < MIN_SIDE {
                break;
            }
            (bin, sigma_k) = (b, sk);
            let next = run(bin, sigma_k);
            let old = std::mem::replace(&mut pass, next);
            if opt.trails && first.is_none() && old.bin == 1 && old.kernel.sigma == FIRST_SIGMA {
                first = Some(old);
            }
            passes += 1;
        }
        if passes > 1 {
            if let Some(v) = measure_psf(&pass) {
                psf = Some(v);
            }
        }
    }
    let psf_var = psf.unwrap_or(pass.kernel.sigma * pass.kernel.sigma);
    let undersampled = pass.bin == 1
        && psf_var <= UNDERSAMPLED_VAR
        && (opt.psf_sigma.is_some() || (measured && psf_var > PSF_VAR_FLOOR));
    let mut det = finish(&mut pass, psf_var, measured, passes, undefined, opt);
    if opt.trails {
        if pass.bin == 1 {
            trails(&pass, &mut det, psf_var, opt, valid.as_deref(), true);
        } else {
            let full = match first.take() {
                Some(f) => f,
                None => run(1, FIRST_SIGMA),
            };
            let b = pass.bin as f64;
            let var = (psf_var * b * b - (b * b - 1.0) / 12.0).max(0.05);
            trails(&full, &mut det, var, opt, valid.as_deref(), false);
        }
    }
    if matches!(det.census.mode, FrameMode::Points | FrameMode::Refused) {
        if undersampled && det.sources.iter().any(|s| s.class == Class::Sharp) {
            for s in det.sources.iter_mut().filter(|s| s.class == Class::Sharp) {
                s.class = Class::Star;
            }
            score_points(&mut det.sources, &pass);
        }
        det.sources
            .truncate(source_limit(&det.sources, opt.max_sources));
    }
    Ok(det)
}

fn source_limit(sources: &[Detected], max_sources: usize) -> usize {
    if max_sources == 0 {
        return sources.len();
    }
    sources
        .iter()
        .enumerate()
        .filter(|(_, s)| s.class.is_stellar())
        .nth(max_sources - 1)
        .map_or(sources.len(), |(i, _)| i + 1)
}

fn trails(
    pass: &Pass,
    det: &mut Detection,
    psf_var: f64,
    opt: &DetectOptions,
    valid: Option<&[bool]>,
    point_classes: bool,
) {
    debug_assert_eq!(pass.bin, 1, "the trail pass runs at full resolution");
    let (w, h) = (pass.w, pass.h);
    let point_sources = &det.sources[..source_limit(&det.sources, opt.max_sources)];
    let mut mask: Vec<bool> = match valid {
        Some(v) if v.len() == w * h => v.iter().map(|ok| !ok).collect(),
        _ => vec![false; w * h],
    };
    let mut set = |x: f64, y: f64, r: f64| {
        let (x0, x1) = (
            (x - r).floor().max(0.0) as usize,
            (x + r).ceil().max(0.0) as usize,
        );
        let (y0, y1) = (
            (y - r).floor().max(0.0) as usize,
            (y + r).ceil().max(0.0) as usize,
        );
        for yy in y0..=y1.min(h - 1) {
            for xx in x0..=x1.min(w - 1) {
                mask[yy * w + xx] = true;
            }
        }
    };
    if point_classes {
        for s in point_sources.iter().filter(|s| s.class == Class::Sharp) {
            set(s.x, s.y, 1.5);
        }
        for t in det.streaks.iter().filter(|t| t.cosmic_ray) {
            let n = t.length.ceil().max(1.0) as usize;
            for k in 0..=n {
                let f = k as f64 / n as f64;
                set(t.x0 + f * (t.x1 - t.x0), t.y0 + f * (t.y1 - t.y0), 1.5);
            }
        }
    }
    let k2 = pass.kernel.sigma * pass.kernel.sigma;
    let mut across: Vec<f64> = pass
        .peaks
        .iter()
        .filter(|p| p.fitted && p.snr() >= 10.0 && p.elongation() > PSF_ELONGATION)
        .map(|p| (p.minor2 - k2).max(0.05))
        .collect();
    let minor = stats::median(&mut across)
        .unwrap_or(psf_var)
        .min(psf_var.max(0.05));
    let mut raw_hits: Vec<Detected> = Vec::new();
    for p in &pass.peaks {
        if sharp_in_raw(p, pass, minor) {
            set(p.x, p.y, 1.5);
            raw_hits.push(hit_source(p, k2, opt.sigma_floor));
        }
    }
    let input = trail::Input {
        raw: &pass.raw,
        orig: &pass.raw,
        s: &pass.s,
        w,
        h,
        mesh_tile: opt.mesh,
        mesh: &pass.mesh,
        clipped: pass.clipped,
        kernel_sigma: pass.kernel.sigma,
        kernel_gain: pass.kernel.noise_gain(),
        mask,
        peaks: pass
            .peaks
            .iter()
            .filter(|p| p.fitted)
            .map(|p| (p.x, p.y, p.minor2, p.snr()))
            .collect(),
        psf_var,
        saturation: det.stats.saturation.map(|v| v as f32),
        gain: opt.gain,
        false_alarms: opt.false_alarms,
        threshold: opt.threshold,
        sigma_floor: opt.sigma_floor,
        prior: opt.trail_angle,
        streaks: det
            .streaks
            .iter()
            .map(|t| {
                (
                    [t.x0, t.y0, t.x1, t.y1],
                    t.angle,
                    t.length,
                    t.cosmic_ray,
                    t.width,
                    t.flux,
                )
            })
            .collect(),
    };
    let out = trail::run(&input, valid);
    det.census = out.census;
    if let Some(out) = out.kernel {
        det.stats.bin = 1;
        det.stats.psf_sigma = out.psf_sigma;
        det.stats.psf_measured = true;
        det.stats.passes += 2;
        det.stats.background = out.background;
        det.stats.noise = out.noise;
        det.stats.kernel_sigma = out.equivalent_sigma;
        det.stats.threshold = out.threshold;
        det.stats.candidates = out.candidates;
        det.sources
            .truncate(source_limit(&det.sources, opt.max_sources));
        retain_streaks(det, &out.keep_streak);
        add_hits(&mut det.sources, raw_hits);
        det.sources.extend(out.sources);
        sort_sources(&mut det.sources, opt.max_sources);
        det.kernel = Some(out.info);
        return;
    }
    let Some(t) = out.trail else { return };
    det.stats.psf_sigma = t.width;
    det.stats.psf_measured = true;
    det.stats.kernel_sigma = out.kernel_sigma;
    det.stats.threshold = out.threshold;
    det.stats.passes += usize::from(det.stats.bin > 1);
    det.stats.bin = pass.bin;
    let trail = t;
    let keep: Vec<bool> = det
        .streaks
        .iter()
        .map(|s| {
            !trail::star_like(s.angle, s.length, s.cosmic_ray, &trail)
                && !out.chains.contains(&[s.x0, s.y0, s.x1, s.y1])
        })
        .collect();
    det.sources
        .truncate(source_limit(&det.sources, opt.max_sources));
    retain_streaks(det, &keep);
    add_hits(&mut det.sources, raw_hits);
    let sources = &mut det.sources;
    let major = (trail.length * trail.length / 12.0 + trail.width * trail.width).sqrt();
    // Hits the first pass flagged are masked for the trail fit, but the rest of a hit's track
    // can still come back as a point: one hit, reported once.
    let hits: Vec<(f64, f64)> = sources
        .iter()
        .filter(|s| s.class == Class::Sharp)
        .map(|s| (s.x, s.y))
        .collect();
    for s in &out.sources {
        let point = s.kind != trail::Kind::Trail;
        if point
            && hits
                .iter()
                .any(|h| (h.0 - s.x).phypot(h.1 - s.y) < HIT_TRACK_REACH)
        {
            continue;
        }
        let (along, across) = (s.sigma_along, s.sigma_across);
        sources.push(Detected {
            x: s.x,
            y: s.y,
            sigma_px: along.max(across),
            sigma_along: along,
            sigma_across: across,
            flux: s.flux,
            snr: s.peak,
            p_star: s.p_star,
            sigma_major: if point { trail.width } else { major },
            sigma_minor: trail.width,
            theta: trail.angle,
            class: match s.kind {
                trail::Kind::Trail => Class::Star,
                trail::Kind::Point => Class::Point,
                trail::Kind::Hit => Class::Sharp,
            },
            streak: None,
            truncated: s.truncated,
        });
    }
    sort_sources(sources, opt.max_sources);
    det.trail = Some(trail);
}

/// A hit that [`sharp_in_raw`] found at full resolution while the trail fit was prepared.
fn hit_source(p: &Peak, k2: f64, sigma_floor: f64) -> Detected {
    let (minor2, major2) = if p.fitted {
        (p.minor2, p.major2)
    } else {
        (k2, k2)
    };
    let snr_amp = p.amp / p.noise;
    let sigma = (minor2 + major2) / (std::f64::consts::SQRT_2 * k2.sqrt() * snr_amp);
    let sigma_px = sigma.max(sigma_floor);
    Detected {
        x: p.x,
        y: p.y,
        sigma_px,
        sigma_along: sigma_px,
        sigma_across: sigma_px,
        flux: 2.0 * std::f64::consts::PI * p.amp * (minor2 * major2.min(100.0 * minor2)).sqrt(),
        snr: p.snr(),
        p_star: 0.0,
        sigma_major: (major2 - k2).max(0.0).sqrt(),
        sigma_minor: (minor2 - k2).max(0.0).sqrt(),
        theta: p.theta,
        class: Class::Sharp,
        streak: None,
        truncated: false,
    }
}

/// Reports the hits found at full resolution that the first pass did not flag: on a binned
/// pass a short track blurs into a star, and next to a trail a hit is taken into it. The trail
/// fit masks them either way; this keeps them in the report.
fn add_hits(sources: &mut Vec<Detected>, hits: Vec<Detected>) {
    let mut known: Vec<(f64, f64)> = sources
        .iter()
        .filter(|s| s.class == Class::Sharp)
        .map(|s| (s.x, s.y))
        .collect();
    for h in hits {
        if known
            .iter()
            .all(|k| (k.0 - h.x).phypot(k.1 - h.y) >= HIT_TRACK_REACH)
        {
            known.push((h.x, h.y));
            sources.push(h);
        }
    }
}

fn retain_streaks(det: &mut Detection, keep: &[bool]) {
    let mut map = vec![None; det.streaks.len()];
    let mut streaks = Vec::new();
    for (i, s) in det.streaks.iter().enumerate() {
        if keep[i] {
            map[i] = Some(streaks.len());
            let mut s = s.clone();
            s.id = streaks.len();
            streaks.push(s);
        }
    }
    det.sources = det
        .sources
        .iter()
        .filter(|s| match s.class {
            Class::Sharp => true,
            Class::Streak => s.streak.is_some_and(|i| keep[i]),
            _ => false,
        })
        .map(|s| Detected {
            streak: s.streak.and_then(|i| map[i]),
            ..*s
        })
        .collect();
    det.streaks = streaks;
}

fn sort_sources(sources: &mut Vec<Detected>, max_sources: usize) {
    sources.sort_by(|a, b| {
        b.flux
            .total_cmp(&a.flux)
            .then(a.y.total_cmp(&b.y))
            .then(a.x.total_cmp(&b.x))
    });
    sources.truncate(source_limit(sources, max_sources));
}

fn plan_for(psf: f64, at_least: usize, at_most: usize) -> (usize, f64) {
    let mut bin = at_least.max(1);
    let at = |b: usize| ((psf * psf + (b * b - 1) as f64 / 12.0).sqrt()) / b as f64;
    while at(bin) > BIN_ABOVE && bin < at_most.min(64) {
        bin += 1;
    }
    (bin, at(bin).clamp(KERNEL_RANGE.0, KERNEL_RANGE.1))
}

fn scaled_psf(psf: f64, bin: usize) -> f64 {
    (psf * psf + (bin * bin - 1) as f64 / 12.0) / (bin * bin) as f64
}

fn run_pass<'a>(
    image: &'a [f32],
    valid: Option<&[bool]>,
    (width, height): (usize, usize),
    bin: usize,
    sigma_k: f64,
    opt: &DetectOptions,
    clipped: Option<&'a ClippedBackground>,
) -> Pass<'a> {
    let (raw, w, h): (Cow<[f32]>, usize, usize) = if bin == 1 {
        (Cow::Borrowed(image), width, height)
    } else {
        let (b, bw, bh) = filter::bin(image, width, height, bin);
        (Cow::Owned(b), bw, bh)
    };
    let valid: Option<Cow<[bool]>> = valid.map(|v| {
        if bin == 1 {
            Cow::Borrowed(v)
        } else {
            Cow::Owned(
                (0..w * h)
                    .map(|i| {
                        let (bx, by) = (i % w, i / w);
                        (0..bin)
                            .all(|dy| (0..bin).all(|dx| v[(by * bin + dy) * width + bx * bin + dx]))
                    })
                    .collect(),
            )
        }
    });
    let kernel = Kernel::gaussian(sigma_k);
    let (mut tmp, mut s) = (Vec::new(), Vec::new());
    filter::smooth(&raw, w, h, &kernel, &mut tmp, &mut s);
    drop(tmp);
    let tile = (opt.mesh / bin).max(24);
    let mesh = if let Some(sky) = clipped {
        let gain = kernel.noise_gain();
        let cubes: f64 = kernel.taps.iter().map(|&v| f64::from(v).ppowi(3)).sum();
        sky.filtered(
            w,
            h,
            gain / bin as f64,
            cubes * cubes / gain.ppowi(3) / bin as f64,
        )
    } else {
        Mesh::build(&s, w, h, tile, kernel.radius + 1, valid.as_deref())
    };
    let threshold = opt
        .threshold
        .unwrap_or_else(|| stats::threshold_for(opt.false_alarms, (w * h) as f64, sigma_k));
    let mut found = peaks::find(&s, w, h, &mesh, threshold as f32, valid.as_deref());
    if found.len() > MAX_PEAKS {
        let height = |&(x, y): &(usize, usize)| {
            let (bg, noise) = mesh.at(x as f64, y as f64);
            (s[y * w + x] - bg) / noise
        };
        found.select_nth_unstable_by(MAX_PEAKS, |a, b| {
            height(b)
                .total_cmp(&height(a))
                .then(a.1.cmp(&b.1))
                .then(a.0.cmp(&b.0))
        });
        found.truncate(MAX_PEAKS);
        found.sort_unstable_by_key(|&(x, y)| (y, x));
    }
    let peaks = found
        .into_iter()
        .map(|(x, y)| {
            let (bg, noise) = mesh.at(x as f64, y as f64);
            peaks::measure(&s, w, x, y, f64::from(bg), f64::from(noise))
        })
        .collect();
    Pass {
        bin,
        kernel,
        w,
        h,
        raw,
        s,
        mesh,
        clipped,
        peaks,
        threshold,
    }
}

fn measure_psf(pass: &Pass) -> Option<f64> {
    psf_of(measure_round_psf(pass), measure_trail_width(pass))
}

/// [`measure_psf`]'s choice between the round peaks' width and the trails'.
fn psf_of(round: Option<f64>, trail: Option<f64>) -> Option<f64> {
    match (round, trail) {
        (Some(r), Some(t)) if r <= HIT_MODE_VAR && t > TRAIL_WIDTH_DISAGREEMENT * r => Some(t),
        _ => round,
    }
}

/// The stars' variance across their trail, kernel removed, px² of the pass's image: the median
/// narrow-axis variance of the bright, elongated peaks that share the frame's dominant
/// direction. `None` unless that direction stands out: at least [`MIN_TRAIL_PEAKS`] peaks within
/// 7.5° of it, three times what peaks of random directions would put there, and spread over at
/// least [`MIN_TRAIL_LINES`] parallel lines so that one streak, or a few, cannot pass for it.
/// Cosmic-ray tracks point anywhere; stars smeared by one motion share it.
fn measure_trail_width(pass: &Pass) -> Option<f64> {
    trail_width(&pass.peaks, pass.kernel.sigma * pass.kernel.sigma)
}

/// [`measure_trail_width`] over a pass's peaks; `k2` is the kernel's variance.
fn trail_width(peaks: &[Peak], k2: f64) -> Option<f64> {
    const BINS: usize = 36;
    let pick = |min_snr: f64| -> Vec<&Peak> {
        peaks
            .iter()
            .filter(|p| {
                p.fitted
                    && p.snr() >= min_snr
                    && p.minor2.is_finite()
                    && p.elongation() > PSF_ELONGATION
            })
            .collect()
    };
    let mut v = pick(20.0);
    if v.len() < MIN_TRAIL_PEAKS {
        v = pick(8.0);
    }
    if v.len() < MIN_TRAIL_PEAKS {
        return None;
    }
    // Directions are axial: θ in (−π/2, π/2], bins of 5° that wrap around.
    let bin_of = |t: f64| {
        let u = (t + std::f64::consts::FRAC_PI_2) / std::f64::consts::PI;
        ((u * BINS as f64) as usize).min(BINS - 1)
    };
    let mut count = [0usize; BINS];
    for p in &v {
        count[bin_of(p.theta)] += 1;
    }
    let (best, in_window) = (0..BINS)
        .map(|i| {
            let n = count[(i + BINS - 1) % BINS] + count[i] + count[(i + 1) % BINS];
            (i, n)
        })
        .max_by_key(|&(i, n)| (n, std::cmp::Reverse(i)))?;
    let random = 3.0 * v.len() as f64 / BINS as f64;
    if in_window < MIN_TRAIL_PEAKS || (in_window as f64) < 3.0 * random {
        return None;
    }
    let centre =
        -std::f64::consts::FRAC_PI_2 + (best as f64 + 0.5) * std::f64::consts::PI / BINS as f64;
    let half = 7.5f64.to_radians().psin();
    let members: Vec<&Peak> = v
        .into_iter()
        .filter(|p| (p.theta - centre).psin().abs() <= half)
        .collect();
    // Offsets across the direction, in bins a few star widths wide: distinct lines.
    let (s, c) = centre.psin_cos();
    let mut lines: Vec<i64> = members
        .iter()
        .map(|p| ((c * p.y - s * p.x) / TRAIL_LINE_BIN).floor() as i64)
        .collect();
    lines.sort_unstable();
    lines.dedup();
    if lines.len() < MIN_TRAIL_LINES {
        return None;
    }
    let mut w: Vec<f64> = members.iter().map(|p| p.minor2 - k2).collect();
    stats::median(&mut w).map(|m| m.max(PSF_VAR_FLOOR))
}

/// Least number of elongated peaks sharing the dominant direction for [`measure_trail_width`].
const MIN_TRAIL_PEAKS: usize = 12;
/// Least number of distinct parallel lines those peaks lie on.
const MIN_TRAIL_LINES: usize = 8;
/// Width of a line across the dominant direction, pixels of the pass's image.
const TRAIL_LINE_BIN: f64 = 8.0;
/// Variance ratio of the trails' width to the round peaks' beyond which the trails' is the
/// stars' ([`measure_psf`]).
const TRAIL_WIDTH_DISAGREEMENT: f64 = 2.0;
/// Round-peak variance, px² of the pass's image, at or below which the round mode may be hits
/// rather than stars: a FWHM of 1.5 px. Hits whose charge spreads into the neighbouring pixels
/// measure 0.9 to 1.4 px.
const HIT_MODE_VAR: f64 = (1.5 / FWHM_PER_SIGMA) * (1.5 / FWHM_PER_SIGMA);

/// The widest distinct size mode of the bright, round peaks ([`measure_psf`]).
fn measure_round_psf(pass: &Pass) -> Option<f64> {
    let k2 = pass.kernel.sigma * pass.kernel.sigma;
    let pick = |min_snr: f64| -> Vec<f64> {
        pass.peaks
            .iter()
            .filter(|p| p.fitted && p.snr() >= min_snr && p.elongation() <= PSF_ELONGATION)
            .map(|p| p.size2() - k2)
            .collect()
    };
    let mut v = pick(20.0);
    if v.len() < 5 {
        v = pick(8.0);
    }
    if v.len() < 3 {
        return None;
    }
    v.sort_by(f64::total_cmp);
    let middle = |s: &[f64]| s[s.len() / 2];
    if v.len() >= 12 {
        let x: Vec<f64> = v.iter().map(|s| (s.max(0.0) + SIZE_OFFSET).pln()).collect();
        let (lo, hi) = (x[0], x[x.len() - 1]);
        let bins = (((hi - lo) / SIZE_BIN) as usize + 1).min(400);
        let mut count = vec![0f64; bins + 2];
        for xi in &x {
            count[(((xi - lo) / SIZE_BIN) as usize).min(bins - 1) + 1] += 1.0;
        }
        let smooth: Vec<f64> = (1..=bins)
            .map(|i| count[i - 1] + 2.0 * count[i] + count[i + 1])
            .collect();
        let tallest = smooth.iter().copied().fold(0.0, f64::max);
        for i in (0..bins).rev() {
            let s = smooth[i];
            let is_max = (i + 1 >= bins || s > smooth[i + 1]) && (i == 0 || s >= smooth[i - 1]);
            if !is_max || s < MODE_MIN || s < 0.05 * tallest {
                continue;
            }
            let dip = smooth[..i].iter().rev().any(|d| *d <= 0.5 * s);
            let anything_higher_below = smooth[..i].iter().any(|d| *d > s);
            let nearly_tallest = s >= MODE_NEARLY_TALLEST * tallest;
            if dip || !anything_higher_below || nearly_tallest {
                let centre = lo + (i as f64 + 0.5) * SIZE_BIN;
                let members: Vec<f64> = v
                    .iter()
                    .zip(&x)
                    .filter(|(_, xi)| (*xi - centre).abs() <= 1.75 * SIZE_BIN)
                    .map(|(s, _)| *s)
                    .collect();
                if !members.is_empty() {
                    return Some(middle(&members).max(PSF_VAR_FLOOR));
                }
            }
        }
    }
    let resolved: Vec<f64> = v.iter().copied().filter(|x| *x > RESOLVED).collect();
    let picked = if resolved.len() >= 5 {
        &resolved[..]
    } else {
        &v[..]
    };
    Some(middle(picked).max(PSF_VAR_FLOOR))
}

const PSF_ELONGATION: f64 = 2.0;
/// Distance, pixels, within which a point found by the trail fit belongs to a flagged hit.
const HIT_TRACK_REACH: f64 = 4.0;

#[derive(Debug, Clone, Copy)]
struct StarShape {
    minor: f64,
    stretch: f64,
}

fn measure_shape(pass: &Pass, psf_var: f64) -> StarShape {
    let k2 = pass.kernel.sigma * pass.kernel.sigma;
    let centre = (psf_var + SIZE_OFFSET).pln();
    let pick = |min_snr: f64| -> (Vec<f64>, Vec<f64>) {
        pass.peaks
            .iter()
            .filter(|p| p.fitted && p.snr() >= min_snr && p.elongation() <= PSF_ELONGATION)
            .filter(|p| {
                let x = ((p.size2() - k2).max(0.0) + SIZE_OFFSET).pln();
                (x - centre).abs() <= 1.75 * SIZE_BIN
            })
            .map(|p| (p.minor2 - k2, p.elongation()))
            .unzip()
    };
    let (mut minor, mut elong) = pick(20.0);
    if minor.len() < 5 {
        (minor, elong) = pick(8.0);
    }
    match (stats::median(&mut minor), stats::median(&mut elong)) {
        (Some(m), Some(e)) if minor.len() >= 3 => StarShape {
            minor: m.clamp(0.02, psf_var),
            stretch: (e / 1.1).max(1.0),
        },
        _ => StarShape {
            minor: psf_var,
            stretch: 1.0,
        },
    }
}

/// The star shape that [`classify`] judges sharpness against, on a frame whose stars trail.
///
/// On such a frame no star is round, and the round peaks [`measure_shape`] reads are whatever
/// else the frame holds. Hit tracks a few pixels long are round enough and as large as stars
/// while a pixel thin: their narrow-axis variance was taken for the stars' (0.06 px² for
/// trails 1.8 px² across), and no hit was narrower than that. On a binned pass, fragments of
/// trails make the round mode as wide as a trail is long (3.6 px² for trails 0.5 px² across),
/// and a tracked object, round and as wide as a star, was flagged as a hit. Across a trail the
/// profile is the star's, so when the trails ([`measure_trail_width`]) and the round peaks
/// disagree by more than [`TRAIL_WIDTH_DISAGREEMENT`] either way, sharpness is judged against
/// the trails. The stars' width ([`measure_psf`]) is left as it is. Frames without a shared
/// trail direction keep their shape.
fn trail_shape(shape: StarShape, trail: Option<f64>) -> StarShape {
    match trail {
        Some(t)
            if t > TRAIL_WIDTH_DISAGREEMENT * shape.minor
                || shape.minor > TRAIL_WIDTH_DISAGREEMENT * t =>
        {
            StarShape {
                minor: t,
                stretch: 1.0,
            }
        }
        _ => shape,
    }
}

const SIZE_BIN: f64 = 0.25;
const MODE_MIN: f64 = 12.0;
const MODE_NEARLY_TALLEST: f64 = 0.75;

const SIZE_OFFSET: f64 = 0.05;

const RESOLVED: f64 = 0.12;

const PSF_VAR_FLOOR: f64 = 0.02;

const UNDERSAMPLED_VAR: f64 = 1.0 / (FWHM_PER_SIGMA * FWHM_PER_SIGMA);

const SATURATED_AREA: usize = 6;
const SATURATED_SIGMA: (f64, f64) = (0.1, 0.03);

pub(crate) const FWHM_PER_SIGMA: f64 = 2.354_820_045_030_949;

fn stand_out(pass: &mut Pass, psf_var: f64, saturation: Option<f32>) {
    let k2 = pass.kernel.sigma * pass.kernel.sigma;
    let rho = 3.0 * (psf_var + k2).sqrt();
    let mut kept = Vec::with_capacity(pass.peaks.len());
    for p in &pass.peaks {
        if saturation.is_some_and(|level| flat_top(pass, p).0 >= level) {
            kept.push(*p);
            continue;
        }
        let ring = peaks::ring_median(&pass.s, pass.w, pass.h, p.ix as f64, p.iy as f64, rho);
        match ring {
            None => kept.push(*p),
            Some(r) if (r - p.bg).abs() <= 2.0 * peaks::RING_NOISE * p.noise => kept.push(*p),
            Some(r) => {
                let height = f64::from(pass.s[p.iy * pass.w + p.ix]) - r;
                if height > (pass.threshold - 1.0) * p.noise {
                    kept.push(peaks::measure(&pass.s, pass.w, p.ix, p.iy, r, p.noise));
                }
            }
        }
    }
    pass.peaks = kept;
}

fn finish(
    pass: &mut Pass,
    psf_var: f64,
    psf_measured: bool,
    passes: usize,
    undefined: usize,
    opt: &DetectOptions,
) -> Detection {
    let k2 = pass.kernel.sigma * pass.kernel.sigma;
    let b = pass.bin as f64;
    let (typ_bg, typ_noise) = pass.mesh.typical();
    let gain = pass.kernel.noise_gain();
    let raw_noise = f64::from(typ_noise) / gain * b;

    let floor = f64::from(typ_bg) + 20.0 * f64::from(typ_noise) / gain;
    let saturation = find_saturation(pass, opt).filter(|level| f64::from(*level) > floor);
    let candidates = pass.peaks.len();
    let shape = trail_shape(measure_shape(pass, psf_var), measure_trail_width(pass));
    stand_out(pass, psf_var, saturation);
    let mut classes: Vec<Class> = pass
        .peaks
        .iter()
        .map(|p| classify(p, pass, psf_var, shape, saturation))
        .collect();
    let n = pass.peaks.len();
    let mut plateaus: Vec<Option<saturate::Plateau>> = vec![None; n];
    let mut absorbed = vec![false; n];
    if let Some(level) = saturation {
        let mut order: Vec<usize> = (0..n).filter(|&i| classes[i] == Class::Saturated).collect();
        order.sort_by(|&a, &b| {
            pass.peaks[b]
                .height
                .total_cmp(&pass.peaks[a].height)
                .then(a.cmp(&b))
        });
        if !order.is_empty() {
            let mut regions = saturate::Regions::new(pass.w, pass.h);
            let mut owner: Vec<Option<usize>> = Vec::new();
            for i in order {
                let p = &pass.peaks[i];
                let Some(k) = regions.at(&pass.raw, p.ix, p.iy, level) else {
                    continue;
                };
                owner.resize(regions.info.len(), None);
                let pl = regions.info[k];
                if pl.linear || pl.truncated {
                    classes[i] = Class::Extended;
                } else if owner[k].is_some() {
                    absorbed[i] = true;
                } else {
                    owner[k] = Some(i);
                    plateaus[i] = Some(pl);
                }
            }
            let pad = (1.0 + (2.0 * (psf_var + k2).sqrt()).ceil()) as usize;
            let star = |k: usize| owner.get(k).copied().flatten().is_some();
            for (j, q) in pass.peaks.iter().enumerate() {
                if plateaus[j].is_none()
                    && !absorbed[j]
                    && regions.near(q.x, q.y, pad, star).is_some_and(star)
                {
                    absorbed[j] = true;
                }
            }
        }
    }

    let mut streaks: Vec<Streak> = Vec::new();
    let mut on_streak: Vec<Option<usize>> = vec![None; n];
    if opt.streaks {
        let traced = streak::trace_all(
            &pass.s,
            pass.w,
            pass.h,
            &pass.mesh,
            &pass.peaks,
            &classes,
            &absorbed,
            k2,
            psf_var,
            (shape.minor, shape.stretch),
        );
        for (i, member) in traced.members.iter().enumerate() {
            if let Some(si) = member {
                classes[i] = Class::Streak;
                on_streak[i] = Some(*si);
            }
        }
        streaks = traced.streaks;
    }

    let to_input = |v: f64| b * v + 0.5 * (b - 1.0);
    let margin = 1.0 + (2.0 * (psf_var + k2).sqrt()).ceil();
    let mut sources: Vec<Detected> = pass
        .peaks
        .iter()
        .enumerate()
        .filter(|(i, _)| !absorbed[*i])
        .map(|(i, p)| {
            let class = classes[i];
            let total2 = if p.fitted { p.size2() } else { psf_var + k2 };
            let snr_amp = p.amp / p.noise;
            let mut sigma = total2 / (std::f64::consts::SQRT_2 * pass.kernel.sigma * snr_amp);
            let minor2 = if p.fitted { p.minor2 } else { psf_var + k2 };
            let major2 = if p.fitted { p.major2 } else { psf_var + k2 };
            let mut flux = 2.0
                * std::f64::consts::PI
                * p.amp
                * (minor2 * major2.min(100.0 * minor2)).sqrt()
                * b
                * b;
            if let Some(gain) = opt.gain {
                let sp2 = (total2 - k2).max(1.0 / 12.0);
                let electrons = (flux * gain).max(1.0);
                let t8 = total2 * total2 * total2 * total2;
                let d = k2 + 2.0 * sp2;
                sigma = (sigma * sigma + sp2 * t8 / (electrons * k2 * k2 * d * d)).sqrt();
            }
            let (mut px, mut py) = (p.x, p.y);
            if let (Some(pl), Some(level)) = (plateaus[i], saturation) {
                if pl.area >= SATURATED_AREA {
                    (px, py) = (pl.x, pl.y);
                    sigma = SATURATED_SIGMA.0 + SATURATED_SIGMA.1 * pl.radius;
                } else {
                    sigma = sigma.max(SATURATED_SIGMA.0);
                }
                let top = (f64::from(level) - p.bg).max(p.amp);
                flux =
                    flux.max(top * (pl.area as f64 + 2.0 * std::f64::consts::PI * total2) * b * b);
            }
            let sigma_px = (sigma * b).max(opt.sigma_floor);
            let near_edge = px < margin
                || py < margin
                || px > pass.w as f64 - 1.0 - margin
                || py > pass.h as f64 - 1.0 - margin;
            Detected {
                x: to_input(px),
                y: to_input(py),
                sigma_px,
                sigma_along: sigma_px,
                sigma_across: sigma_px,
                flux,
                snr: p.snr(),
                p_star: 0.0,
                sigma_major: b * (major2 - k2).max(0.0).sqrt(),
                sigma_minor: b * (minor2 - k2).max(0.0).sqrt(),
                theta: p.theta,
                class,
                streak: on_streak[i],
                truncated: near_edge,
            }
        })
        .collect();
    score_points(&mut sources, pass);
    sources.sort_by(|a, b| {
        b.flux
            .total_cmp(&a.flux)
            .then(a.y.total_cmp(&b.y))
            .then(a.x.total_cmp(&b.x))
    });
    for s in &mut streaks {
        s.rescale(b);
    }
    streaks.sort_by(|a, b| b.length.total_cmp(&a.length));
    if !streaks.is_empty() {
        let mut map = vec![0usize; streaks.len()];
        for (new, s) in streaks.iter().enumerate() {
            map[s.id] = new;
        }
        for s in &mut sources {
            s.streak = s.streak.map(|i| map[i]);
        }
        for (i, s) in streaks.iter_mut().enumerate() {
            s.id = i;
        }
    }

    Detection {
        sources,
        streaks,
        stats: FrameStats {
            background: f64::from(typ_bg),
            noise: raw_noise,
            psf_sigma: (psf_var * b * b - (b * b - 1.0) / 12.0).max(0.0).sqrt(),
            psf_measured,
            bin: pass.bin,
            kernel_sigma: pass.kernel.sigma,
            threshold: pass.threshold,
            passes,
            candidates,
            saturation: saturation.map(f64::from),
            undefined,
        },
        census: Census {
            mode: FrameMode::Points,
            reason: "not examined".into(),
            votes: 0,
            mode_votes: 0,
            mode_share: 0.0,
            second_share: 0.0,
            elongated_share: 0.0,
            shapes: None,
            autocov: None,
            regions: (0, 0.0, 0.0),
        },
        trail: None,
        kernel: None,
    }
}

fn score_points(sources: &mut [Detected], pass: &Pass) {
    let tests = (pass.w * pass.h) as f64;
    let mut heights: Vec<f64> = sources
        .iter()
        .filter(|s| s.class.is_stellar())
        .map(|s| s.snr)
        .collect();
    heights.sort_by(f64::total_cmp);
    for s in sources {
        let p_real = 1.0
            - stats::local_false_discovery_rate(
                s.snr,
                pass.threshold,
                pass.kernel.sigma,
                tests,
                &heights,
            )
            .min(1.0);
        let shape = match s.class {
            Class::Star => 1.0,
            Class::Saturated => 0.9,
            Class::Extended => 0.5,
            Class::Sharp => 0.02,
            Class::Streak => 0.01,
            Class::Point => 0.0,
        };
        s.p_star = (p_real * shape).clamp(0.0, 1.0);
    }
}

fn find_saturation(pass: &Pass, opt: &DetectOptions) -> Option<f32> {
    if let Some(s) = opt.saturation {
        return Some(s);
    }
    let tops: Vec<(f32, usize)> = pass
        .peaks
        .iter()
        .filter(|p| p.snr() > 50.0)
        .map(|p| flat_top(pass, p))
        .collect();
    let highest = tops.iter().map(|t| t.0).fold(f32::NEG_INFINITY, f32::max);
    let (bg, _) = pass.mesh.typical();
    let mut flat: Vec<f64> = tops
        .iter()
        .filter(|(top, count)| *count >= 4 && *top - bg >= 0.9 * (highest - bg))
        .map(|(top, _)| f64::from(*top))
        .collect();
    stats::median(&mut flat).map(|m| (0.98 * m) as f32)
}

fn flat_top(pass: &Pass, p: &Peak) -> (f32, usize) {
    let (w, h) = (pass.w as isize, pass.h as isize);
    let at = |dx: isize, dy: isize| -> Option<f32> {
        let (x, y) = (p.ix as isize + dx, p.iy as isize + dy);
        (x >= 0 && y >= 0 && x < w && y < h).then(|| pass.raw[(y * w + x) as usize])
    };
    let mut top = f32::NEG_INFINITY;
    for dy in -1..=1 {
        for dx in -1..=1 {
            if let Some(v) = at(dx, dy) {
                top = top.max(v);
            }
        }
    }
    let tol = 0.005 * (top - p.bg as f32).abs();
    let mut count = 0;
    for dy in -2..=2 {
        for dx in -2..=2 {
            if at(dx, dy).is_some_and(|v| v >= top - tol) {
                count += 1;
            }
        }
    }
    (top, count)
}

fn classify(
    p: &Peak,
    pass: &Pass,
    psf_var: f64,
    shape: StarShape,
    saturation: Option<f32>,
) -> Class {
    let k2 = pass.kernel.sigma * pass.kernel.sigma;
    if !p.fitted {
        return if p.snr() > 2.0 * pass.threshold {
            Class::Sharp
        } else {
            Class::Star
        };
    }
    let rel = 2.0 / p.snr();
    let minor = p.minor2 - k2;
    let tol = 2.5 * rel * p.minor2;
    if (shape.minor >= RESOLVED && minor < 0.4 * shape.minor - tol)
        || sharp_in_raw(p, pass, shape.minor)
    {
        return Class::Sharp;
    }
    if let Some(level) = saturation {
        if flat_top(pass, p).0 >= level {
            return Class::Saturated;
        }
    }
    let elong = p.elongation();
    if elong > 1.5 * shape.stretch + 6.0 * rel {
        return Class::Extended;
    }
    let size = p.size2() - k2;
    if size > 2.0 * psf_var + 2.5 * rel * p.size2() {
        return Class::Extended;
    }
    Class::Star
}

fn sharp_in_raw(p: &Peak, pass: &Pass, psf_minor: f64) -> bool {
    let (w, h) = (pass.w, pass.h);
    if p.ix < 2 || p.iy < 2 || p.ix + 2 >= w || p.iy + 2 >= h {
        return false;
    }
    let (mut bx, mut by, mut top) = (p.ix, p.iy, f32::NEG_INFINITY);
    for y in p.iy - 1..=p.iy + 1 {
        for x in p.ix - 1..=p.ix + 1 {
            let v = pass.raw[y * w + x];
            if v > top {
                (bx, by, top) = (x, y, v);
            }
        }
    }
    let bg = p.bg;
    let c = f64::from(top) - bg;
    let mut around = [(bx - 1, by), (bx + 1, by), (bx, by - 1), (bx, by + 1)]
        .map(|(x, y)| f64::from(pass.raw[y * w + x]) - bg);
    around.sort_unstable_by(f64::total_cmp);
    let expected = (-1.0 / psf_minor.max(0.1)).pexp();
    let raw_noise = p.noise / pass.kernel.noise_gain();
    c > 10.0 * raw_noise && around[1] < 0.5 * expected * c - 2.0 * raw_noise
}

#[cfg(test)]
mod tests {
    use super::*;
    use synth::Rng;

    /// A bright peak at (x, y) with variances along and across `theta`, kernel included.
    fn peak(x: f64, y: f64, major2: f64, minor2: f64, theta: f64) -> Peak {
        Peak {
            ix: x as usize,
            iy: y as usize,
            x,
            y,
            height: 500.0,
            amp: 500.0,
            noise: 10.0,
            bg: 0.0,
            major2,
            minor2,
            theta,
            fitted: true,
        }
    }

    const K2: f64 = 1.44;

    /// Trails across the frame in one direction, minor variance 1.2 px² beyond the kernel.
    fn trails(n: usize, rng: &mut Rng) -> Vec<Peak> {
        (0..n)
            .map(|_| {
                let (x, y) = (rng.uniform() * 760.0, rng.uniform() * 760.0);
                let theta = 0.6 + (rng.uniform() - 0.5) * 0.04;
                peak(
                    x,
                    y,
                    f64::INFINITY,
                    K2 + 1.2 + 0.1 * (rng.uniform() - 0.5),
                    theta,
                )
            })
            .collect()
    }

    #[test]
    fn trail_width_comes_from_the_shared_direction() {
        let mut rng = Rng(5);
        let mut peaks = trails(60, &mut rng);
        // Ten times as many elongated hit tracks in random directions, narrower than stars.
        for _ in 0..600 {
            let (x, y) = (rng.uniform() * 760.0, rng.uniform() * 760.0);
            let theta = (rng.uniform() - 0.5) * std::f64::consts::PI;
            peaks.push(peak(x, y, K2 + 3.0, K2 + 0.1, theta));
        }
        let w = trail_width(&peaks, K2).expect("a shared direction");
        assert!((w - 1.2).abs() < 0.1, "{w}");
    }

    #[test]
    fn random_directions_or_one_streak_give_no_trail_width() {
        let mut rng = Rng(6);
        let random: Vec<Peak> = (0..300)
            .map(|_| {
                let (x, y) = (rng.uniform() * 760.0, rng.uniform() * 760.0);
                let theta = (rng.uniform() - 0.5) * std::f64::consts::PI;
                peak(x, y, f64::INFINITY, K2 + 1.0, theta)
            })
            .collect();
        assert_eq!(trail_width(&random, K2), None);
        // Forty peaks along one satellite streak: one line, not stars.
        let streak: Vec<Peak> = (0..40)
            .map(|i| {
                let t = i as f64 * 15.0;
                peak(
                    50.0 + t * 0.8,
                    80.0 + t * 0.6,
                    f64::INFINITY,
                    K2 + 2.0,
                    0.6435,
                )
            })
            .collect();
        assert_eq!(trail_width(&streak, K2), None);
    }

    #[test]
    fn the_trail_width_replaces_only_a_mode_of_hits() {
        // Hits taken for stars: a narrow round mode, trails more than twice as wide.
        assert_eq!(psf_of(Some(0.27), Some(1.27)), Some(1.27));
        assert_eq!(psf_of(Some(0.02), Some(0.9)), Some(0.9));
        // A narrow round mode the trails agree with (undersampled stars): kept.
        assert_eq!(psf_of(Some(0.1), Some(0.15)), Some(0.1));
        // Round stars of an ordinary width: kept whatever the trails say, bit for bit.
        assert_eq!(psf_of(Some(1.1), Some(1.3)), Some(1.1));
        assert_eq!(psf_of(Some(0.5), Some(1.3)), Some(0.5));
        assert_eq!(psf_of(Some(6.2), Some(1.75)), Some(6.2));
        // No round mode, or no trails: unchanged.
        assert_eq!(psf_of(None, Some(1.2)), None);
        assert_eq!(psf_of(Some(0.9), None), Some(0.9));
        assert_eq!(psf_of(None, None), None);
    }
}
