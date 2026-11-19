//! Star and streak detection in astronomical images. [`detect`] accepts a row-major pixel array.

#![forbid(unsafe_code)]

mod filter;
mod mesh;
mod peaks;
mod saturate;
mod stats;
mod streak;
pub mod synth;

use std::borrow::Cow;

use fayde_core::Portable;
use filter::Kernel;
use mesh::{ClippedBackground, Mesh};
use peaks::Peak;

pub use streak::Streak;

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
}

impl Class {
    pub fn name(self) -> &'static str {
        match self {
            Self::Star => "star",
            Self::Saturated => "saturated",
            Self::Extended => "extended",
            Self::Sharp => "sharp",
            Self::Streak => "streak",
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
    pub sigma_px: f64,
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
    peaks: Vec<Peak>,
    threshold: f64,
}

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
    let mut pass = run_pass(
        &image,
        valid.as_deref(),
        (width, height),
        bin,
        sigma_k,
        opt,
        clipped.as_ref(),
    )?;
    let mut passes = 1;
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
            pass = run_pass(
                &image,
                valid.as_deref(),
                (width, height),
                bin,
                sigma_k,
                opt,
                clipped.as_ref(),
            )?;
            passes += 1;
        }
        if passes > 1 {
            if let Some(v) = measure_psf(&pass) {
                psf = Some(v);
            }
        }
    }
    let psf_var = psf.unwrap_or(pass.kernel.sigma * pass.kernel.sigma);
    Ok(finish(pass, psf_var, measured, passes, undefined, opt))
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
    clipped: Option<&ClippedBackground>,
) -> Result<Pass<'a>, DetectError> {
    let (raw, w, h): (Cow<[f32]>, usize, usize) = if bin == 1 {
        (Cow::Borrowed(image), width, height)
    } else {
        let (b, bw, bh) = filter::bin(image, width, height, bin);
        (Cow::Owned(b), bw, bh)
    };
    // A binned pixel is usable when all the pixels it averages are.
    let valid: Option<Vec<bool>> = valid.map(|v| {
        if bin == 1 {
            v.to_vec()
        } else {
            (0..w * h)
                .map(|i| {
                    let (bx, by) = (i % w, i / w);
                    (0..bin).all(|dy| (0..bin).all(|dx| v[(by * bin + dy) * width + bx * bin + dx]))
                })
                .collect()
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
    Ok(Pass {
        bin,
        kernel,
        w,
        h,
        raw,
        s,
        mesh,
        peaks,
        threshold,
    })
}

fn measure_psf(pass: &Pass) -> Option<f64> {
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

const SIZE_BIN: f64 = 0.25;
const MODE_MIN: f64 = 12.0;

const SIZE_OFFSET: f64 = 0.05;

const RESOLVED: f64 = 0.12;

pub const MIN_THRESHOLD: f64 = 1.5;
const MODE_NEARLY_TALLEST: f64 = 0.75;
const PSF_VAR_FLOOR: f64 = 0.02;
const FWHM_PER_SIGMA: f64 = 2.354_820_045_030_949;
const UNDERSAMPLED_VAR: f64 = 1.0 / (FWHM_PER_SIGMA * FWHM_PER_SIGMA);
const SATURATED_AREA: usize = 6;
const SATURATED_SIGMA: (f64, f64) = (0.1, 0.03);

/// Half-width of the height window of the local false discovery rate, in σ.
const FDR_WINDOW: f64 = 0.25;

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
    mut pass: Pass,
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
    let saturation = find_saturation(&pass, opt).filter(|level| f64::from(*level) > floor);
    let candidates = pass.peaks.len();
    let shape = measure_shape(&pass, psf_var);
    let undersampled = pass.bin == 1
        && psf_var <= UNDERSAMPLED_VAR
        && (opt.psf_sigma.is_some() || (psf_measured && psf_var > PSF_VAR_FLOOR));
    stand_out(&mut pass, psf_var, saturation);
    let mut classes: Vec<Class> = pass
        .peaks
        .iter()
        .map(|p| classify(p, &pass, psf_var, shape, saturation, undersampled))
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

    // Probability that a peak of this height is not a noise peak: one minus the local false
    // discovery rate (Efron 2008), the expected number of noise peaks with a height within a
    // quarter σ of this one over the number of peaks seen there (the expected number from
    // `stats::noise_peaks_above`). Peaks already known to be artefacts are not counted as seen.
    let tests = (pass.w * pass.h) as f64;
    let mut heights: Vec<f64> = (0..n)
        .filter(|&i| classes[i].is_stellar() && !absorbed[i])
        .map(|i| pass.peaks[i].snr())
        .collect();
    heights.sort_by(f64::total_cmp);
    let p_real: Vec<f64> = pass
        .peaks
        .iter()
        .map(|p| {
            let s = p.snr();
            // Nothing was kept below the threshold less one (`stand_out`).
            let (lo, hi) = ((s - FDR_WINDOW).max(pass.threshold - 1.0), s + FDR_WINDOW);
            let seen = heights.partition_point(|h| *h <= hi) - heights.partition_point(|h| *h < lo);
            let sk = pass.kernel.sigma;
            let noise_peaks =
                tests * (stats::noise_peaks_above(lo, sk) - stats::noise_peaks_above(hi, sk));
            1.0 - (noise_peaks / seen.max(1) as f64).min(1.0)
        })
        .collect();

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
            let shape = match class {
                Class::Star => 1.0,
                Class::Saturated => 0.9,
                Class::Extended => 0.5,
                Class::Sharp => 0.02,
                Class::Streak => 0.01,
            };
            let near_edge = px < margin
                || py < margin
                || px > pass.w as f64 - 1.0 - margin
                || py > pass.h as f64 - 1.0 - margin;
            Detected {
                x: to_input(px),
                y: to_input(py),
                sigma_px,
                flux,
                snr: p.snr(),
                p_star: (p_real[i] * shape).clamp(0.0, 1.0),
                sigma_major: b * (major2 - k2).max(0.0).sqrt(),
                sigma_minor: b * (minor2 - k2).max(0.0).sqrt(),
                theta: p.theta,
                class,
                streak: on_streak[i],
                truncated: near_edge,
            }
        })
        .collect();
    sources.sort_by(|a, b| {
        b.flux
            .total_cmp(&a.flux)
            .then(a.y.total_cmp(&b.y))
            .then(a.x.total_cmp(&b.x))
    });
    if opt.max_sources > 0 {
        // The limit counts what a solver would use; artefacts ranked above the last kept star
        // stay in the list.
        let mut kept = 0;
        let mut cut = sources.len();
        for (i, s) in sources.iter().enumerate() {
            if s.class.is_stellar() {
                kept += 1;
                if kept == opt.max_sources {
                    cut = i + 1;
                    break;
                }
            }
        }
        sources.truncate(cut);
    }
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
    undersampled: bool,
) -> Class {
    let k2 = pass.kernel.sigma * pass.kernel.sigma;
    if !p.fitted {
        return if !undersampled && p.snr() > 2.0 * pass.threshold {
            Class::Sharp
        } else {
            Class::Star
        };
    }
    let rel = 2.0 / p.snr();
    let minor = p.minor2 - k2;
    let tol = 2.5 * rel * p.minor2;
    if !undersampled
        && ((shape.minor >= RESOLVED && minor < 0.4 * shape.minor - tol)
            || sharp_in_raw(p, pass, shape.minor))
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
