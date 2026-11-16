//! Star and streak detection for Fayde.
//!
//! [`detect`] takes an image of point-like stars (a frame tracked at the sidereal rate), with
//! possibly a satellite streak or two, and returns the stars with what a plate solver needs
//! from each: a position, its 1σ error, a flux for ranking and the probability that it is a
//! real star. Streaks are traced and reported separately, and the detections that lie on them
//! are kept out of the star list.
//!
//! The method, in order (sources in `docs/references.md`):
//!
//! 1. **Matched filter.** The image is smoothed with a Gaussian the width of the stars
//!    (`filter.rs`). No width is asked for: the first pass measures it and a second pass is run
//!    only if the first kernel was a poor match. Images whose stars span many pixels are binned
//!    first, which is where most of the speed on large frames comes from.
//! 2. **Background and noise** of the filtered image on a coarse mesh (`mesh.rs`).
//! 3. **Peaks**: local maxima above a threshold set from the number of false detections the
//!    caller accepts per frame (`peaks.rs`). Each must also stand that far above the median of
//!    a ring around it, so that what a coarse background map misses (stray light, halos) does
//!    not pass for sources.
//! 4. **One closed-form fit per peak**: a quadratic in the logarithm of the 3 × 3 pixels around
//!    it gives position, size, elongation, orientation and amplitude. The same numbers classify
//!    the peak (star, hot pixel or cosmic ray, extended object, ridge) and give the position
//!    error.
//! 5. **Streaks**: every ridge-like peak is traced along its direction (`streak.rs`).
//! 6. **Probabilities**: each source's `p_star` combines the chance that a noise peak of its
//!    height exists (false discovery rate) with its shape class.

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
use mesh::Mesh;
use peaks::Peak;

pub use streak::Streak;

/// Smallest image side accepted, pixels.
pub const MIN_SIDE: usize = 16;

/// What [`detect`] is asked to do. `Default` suits a frame about which nothing is known.
#[derive(Debug, Clone, PartialEq)]
pub struct DetectOptions {
    /// Peak threshold in units of the local noise of the filtered image. `None`: set from
    /// `false_alarms`.
    pub threshold: Option<f64>,
    /// Number of pure-noise peaks accepted per frame when `threshold` is `None` (default 1).
    pub false_alarms: f64,
    /// Gaussian σ of the stars, pixels, if known (FWHM / 2.355). `None`: measured.
    pub psf_sigma: Option<f64>,
    /// Level at which pixels saturate. `None`: found from flat-topped stars.
    pub saturation: Option<f32>,
    /// Keep at most this many sources, brightest first (0: all).
    pub max_sources: usize,
    /// Background mesh size, pixels of the input image. On a binned image the mesh covers
    /// the same area, down to 24 binned pixels.
    pub mesh: usize,
    /// Trace streaks (default true). Without it ridge-like peaks are only marked elongated.
    pub streaks: bool,
    /// Smallest position error reported, pixels: the floor set by what the noise model leaves
    /// out (pixel response, PSF asymmetry, and the star's own photon noise when `gain` is not
    /// given).
    pub sigma_floor: f64,
    /// Gain, electrons per input unit. When given, the photon noise of each star is added to
    /// its position error.
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

/// What a detection is taken to be.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Class {
    /// A point source of the frame's PSF width.
    Star,
    /// A star with saturated pixels: real, with a less certain position and flux.
    Saturated,
    /// Wider than the PSF or mildly elongated: a blend, a galaxy, a defocused star.
    Extended,
    /// Narrower than the PSF: a hot pixel or a cosmic-ray hit.
    Sharp,
    /// On a traced streak.
    Streak,
}

impl Class {
    /// Lower-case name, as written in reports.
    pub fn name(self) -> &'static str {
        match self {
            Self::Star => "star",
            Self::Saturated => "saturated",
            Self::Extended => "extended",
            Self::Sharp => "sharp",
            Self::Streak => "streak",
        }
    }

    /// Whether detections of this class go to a plate solver.
    pub fn is_stellar(self) -> bool {
        matches!(self, Self::Star | Self::Saturated | Self::Extended)
    }
}

/// One detection. Positions are zero-based pixels of the input image.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Detected {
    /// Column.
    pub x: f64,
    /// Row.
    pub y: f64,
    /// 1σ position error per axis, pixels.
    pub sigma_px: f64,
    /// Flux above the background, input units × pixels.
    pub flux: f64,
    /// Peak of the filtered image over its local noise.
    pub snr: f64,
    /// Probability that this is a real star, 0–1.
    pub p_star: f64,
    /// Gaussian σ along the major axis, pixels, kernel removed (infinite on a ridge).
    pub sigma_major: f64,
    /// Gaussian σ along the minor axis, pixels, kernel removed.
    pub sigma_minor: f64,
    /// Direction of the major axis, radians from +x towards +y.
    pub theta: f64,
    /// Class.
    pub class: Class,
    /// Index into [`Detection::streaks`] for [`Class::Streak`].
    pub streak: Option<usize>,
    /// The fit window reached the frame edge or an undefined pixel: the position is biased.
    pub truncated: bool,
}

/// Frame-wide measurements.
#[derive(Debug, Clone, PartialEq)]
pub struct FrameStats {
    /// Typical background, input units.
    pub background: f64,
    /// Typical noise per pixel of the input image, input units (1σ).
    pub noise: f64,
    /// Gaussian σ of the stars, pixels of the input image.
    pub psf_sigma: f64,
    /// Whether `psf_sigma` was measured from this frame's stars (else given or assumed).
    pub psf_measured: bool,
    /// Binning used for detection (1: none).
    pub bin: usize,
    /// σ of the matched filter, pixels of the binned image.
    pub kernel_sigma: f64,
    /// Peak threshold used, in units of the filtered image's noise.
    pub threshold: f64,
    /// Filter passes run (1 to 3).
    pub passes: usize,
    /// Peaks above the threshold, before classification.
    pub candidates: usize,
    /// Saturation level used, if any.
    pub saturation: Option<f64>,
    /// Undefined (non-finite) pixels in the input.
    pub undefined: usize,
}

/// The result of [`detect`].
#[derive(Debug, Clone, PartialEq)]
pub struct Detection {
    /// Detections, brightest first. Every class is listed; [`Detection::stars`] gives those a
    /// plate solver should use.
    pub sources: Vec<Detected>,
    /// Traced streaks, longest first.
    pub streaks: Vec<Streak>,
    /// Frame-wide measurements.
    pub stats: FrameStats,
}

impl Detection {
    /// The detections to hand to a plate solver: stars, saturated stars and extended sources.
    pub fn stars(&self) -> impl Iterator<Item = &Detected> {
        self.sources.iter().filter(|s| s.class.is_stellar())
    }
}

/// Why [`detect`] refused an input.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DetectError {
    /// The pixel count does not match `width × height`, or a side is below [`MIN_SIDE`].
    BadShape(String),
    /// An option is out of range.
    BadOption(String),
    /// No pixel is defined.
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

/// Most peaks measured in one pass: far more than any frame of stars holds.
const MAX_PEAKS: usize = 200_000;

/// σ of the first filter pass, pixels: within 7% of the best signal-to-noise ratio for stars
/// of σ 0.8 to 1.75 px (FWHM 1.9 to 4.1 px), the range most frames fall in.
const FIRST_SIGMA: f64 = 1.2;
/// A kernel is kept when its signal-to-noise ratio is at least this fraction of the matched
/// one's, `2·σk·σp / (σk² + σp²)`.
const MIN_EFFICIENCY: f64 = 0.93;
/// Kernel σ range, pixels of the binned image.
const KERNEL_RANGE: (f64, f64) = (0.6, 2.0);
/// Stars wider than this (σ, binned pixels) are binned further.
const BIN_ABOVE: f64 = 2.0;

/// One filter pass at a given binning and kernel.
struct Pass<'a> {
    bin: usize,
    kernel: Kernel,
    w: usize,
    h: usize,
    /// The (binned) image, undefined pixels filled.
    raw: Cow<'a, [f32]>,
    /// The filtered image.
    s: Vec<f32>,
    mesh: Mesh,
    peaks: Vec<Peak>,
    threshold: f64,
}

/// Detect stars and streaks in `data` (`width × height`, row-major, any linear unit).
/// Non-finite pixels are treated as undefined.
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
    if opt.threshold.is_some_and(|t| !(t.is_finite() && t > 0.0)) {
        return bad("threshold");
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

    // Undefined pixels are filled with a typical value so the filter has something to read,
    // and remembered so that nothing is detected on them.
    let undefined = data
        .iter()
        .fold(0usize, |a, v| a + usize::from(!v.is_finite()));
    let (image, valid): (Cow<[f32]>, Option<Vec<bool>>) = if undefined == 0 {
        (Cow::Borrowed(data), None)
    } else {
        if undefined == data.len() {
            return Err(DetectError::Empty);
        }
        let step = (data.len() / 4096).max(1);
        let mut sample: Vec<f64> = data
            .iter()
            .step_by(step)
            .filter(|v| v.is_finite())
            .map(|v| f64::from(*v))
            .collect();
        let fill = match stats::median(&mut sample) {
            Some(m) => m as f32,
            None => *data.iter().find(|v| v.is_finite()).expect("a finite pixel"),
        };
        let valid: Vec<bool> = data.iter().map(|v| v.is_finite()).collect();
        let filled: Vec<f32> = data
            .iter()
            .map(|v| if v.is_finite() { *v } else { fill })
            .collect();
        (Cow::Owned(filled), Some(valid))
    };

    // The filter scale: given, or measured by a first pass and corrected if it was off.
    // Binning must leave an image worth detecting on, whatever width is given or measured.
    let max_bin = (width.min(height) / MIN_SIDE).max(1);
    let (mut bin, mut sigma_k) = match opt.psf_sigma {
        Some(s) => plan_for(s, 1, max_bin),
        None => (1, FIRST_SIGMA),
    };
    let mut pass = run_pass(&image, valid.as_deref(), width, height, bin, sigma_k, opt)?;
    let mut passes = 1;
    let mut psf = opt.psf_sigma.map(|s| scaled_psf(s, pass.bin));
    let mut measured = false;
    if opt.psf_sigma.is_none() {
        // At most two corrections: the first estimate of a very wide star is rough.
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
            // Binning must leave an image worth detecting on.
            if width / b < MIN_SIDE || height / b < MIN_SIDE {
                break;
            }
            (bin, sigma_k) = (b, sk);
            pass = run_pass(&image, valid.as_deref(), width, height, bin, sigma_k, opt)?;
            passes += 1;
        }
        if passes > 1 {
            if let Some(v) = measure_psf(&pass) {
                psf = Some(v);
            }
        }
    }
    // With no star to measure, the PSF is taken to be the kernel.
    let psf_var = psf.unwrap_or(pass.kernel.sigma * pass.kernel.sigma);
    Ok(finish(pass, psf_var, measured, passes, undefined, opt))
}

/// Binning and kernel σ for stars of σ `psf` pixels of the input image, keeping at least the
/// binning `at_least` already in use and at most `at_most`.
fn plan_for(psf: f64, at_least: usize, at_most: usize) -> (usize, f64) {
    let mut bin = at_least.max(1);
    // Pixel integration: binning by b adds the variance of a box of b input pixels.
    let at = |b: usize| ((psf * psf + (b * b - 1) as f64 / 12.0).sqrt()) / b as f64;
    while at(bin) > BIN_ABOVE && bin < at_most.min(64) {
        bin += 1;
    }
    (bin, at(bin).clamp(KERNEL_RANGE.0, KERNEL_RANGE.1))
}

/// The variance of a star of σ `psf` input pixels, in pixels² of an image binned by `bin`.
fn scaled_psf(psf: f64, bin: usize) -> f64 {
    (psf * psf + (bin * bin - 1) as f64 / 12.0) / (bin * bin) as f64
}

fn run_pass<'a>(
    image: &'a [f32],
    valid: Option<&[bool]>,
    width: usize,
    height: usize,
    bin: usize,
    sigma_k: f64,
    opt: &DetectOptions,
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
    let mesh = Mesh::build(&s, w, h, tile, kernel.radius + 1, valid.as_deref());
    // As many noise peaks per frame as the caller accepts, unless the threshold is given.
    let threshold = opt
        .threshold
        .unwrap_or_else(|| stats::threshold_for(opt.false_alarms, (w * h) as f64, sigma_k));
    let mut found = peaks::find(&s, w, h, &mesh, threshold as f32, valid.as_deref());
    if found.len() > MAX_PEAKS {
        // A threshold far down in the noise, or a frame that is all structure: the highest
        // peaks are kept, so that the work after this point is bounded.
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

/// The variance of the stars' profile with the kernel removed, in pixels² of the pass's image,
/// from the bright, round peaks. `None` when there are too few to tell.
///
/// The peaks are not all stars: hot pixels, cosmic-ray hits and charge-transfer tails are
/// compact too, and on a warm or irradiated detector (a space telescope crossing the radiation
/// belts) they outnumber the stars ten to one or more. But nothing that came through the
/// optics is sharper than the optics allow, and the stars all have one width where the
/// artefacts have many. So the stars are the *widest distinct mode* of the size distribution:
/// in a histogram of ln(size), the widest local maximum that holds at least 5% of the tallest
/// one's count and is separated from the narrower sizes by a dip to half its height or less.
/// The width is the median of the sizes within that mode. With one mode this is close to the
/// median of all; blends and galaxies, wider than the stars but with no common size, form a
/// tail and not a mode. The rule is Fayde's own.
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
        // Histogram of ln(size + offset), bins of SIZE_BIN, smoothed 1-2-1. The offset keeps
        // unresolved peaks (size about zero) finite.
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
        // From the widest sizes down: the first local maximum that is high enough and has a
        // dip below it (or nothing below it at all).
        for i in (0..bins).rev() {
            let s = smooth[i];
            let is_max = (i + 1 >= bins || s > smooth[i + 1]) && (i == 0 || s >= smooth[i - 1]);
            if !is_max || s < MODE_MIN || s < 0.05 * tallest {
                continue;
            }
            let dip = smooth[..i].iter().rev().any(|d| *d <= 0.5 * s);
            let anything_higher_below = smooth[..i].iter().any(|d| *d > s);
            if dip || !anything_higher_below {
                let centre = lo + (i as f64 + 0.5) * SIZE_BIN;
                let members: Vec<f64> = v
                    .iter()
                    .zip(&x)
                    .filter(|(_, xi)| (*xi - centre).abs() <= 1.75 * SIZE_BIN)
                    .map(|(s, _)| *s)
                    .collect();
                if !members.is_empty() {
                    return Some(middle(&members).max(0.02));
                }
            }
        }
    }
    // Too few peaks for a histogram. Hot pixels have no width of their own: when enough
    // peaks are clearly resolved, the stars are those.
    let resolved: Vec<f64> = v.iter().copied().filter(|x| *x > RESOLVED).collect();
    let picked = if resolved.len() >= 5 {
        &resolved[..]
    } else {
        &v[..]
    };
    Some(middle(picked).max(0.02))
}

/// Peaks more elongated than this (axis ratio, kernel included) are not used to measure the
/// stars: ridges and blends. Stars smeared by a drift during the exposure stay below it.
const PSF_ELONGATION: f64 = 2.0;

/// The shape the frame's stars share, beyond their mean width.
#[derive(Debug, Clone, Copy)]
struct StarShape {
    /// Variance along the stars' narrow axis, kernel removed, px².
    minor: f64,
    /// How much more elongated the stars are than round ones measure (1 for round stars): the
    /// median axis ratio over 1.1, the ratio noise alone gives a round star.
    stretch: f64,
}

/// The narrow-axis width and the elongation common to the frame's stars: medians over the
/// bright peaks whose mean size is the stars' (`psf_var`, within the width of a mode of
/// [`measure_psf`]).
///
/// A drift during the exposure smears every star alike, two to one or more. Each source's
/// narrow axis is then well under the stars' *mean* width, and a test for "narrower than a
/// star" against the mean throws out the brightest stars as cosmic rays (their sizes are the
/// best measured, so no tolerance saves them). Narrowness and elongation are judged against
/// what the frame's stars show.
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
        // Too few stars to tell: round.
        _ => StarShape {
            minor: psf_var,
            stretch: 1.0,
        },
    }
}

/// Bin width of the histogram of ln(size) in [`measure_psf`]: 25% in variance, the scatter of
/// a star's measured size at a signal-to-noise ratio of 20.
const SIZE_BIN: f64 = 0.25;
/// Least smoothed count (four times a bin's own) for a mode: about four peaks.
const MODE_MIN: f64 = 12.0;

/// Added to a size, px², before its logarithm is taken in [`measure_psf`].
const SIZE_OFFSET: f64 = 0.05;

/// Intrinsic variance, px², above which a peak is certainly wider than a single pixel: a star
/// integrated over pixels has at least 1/12 px² plus its optical width.
const RESOLVED: f64 = 0.12;

/// A saturated region of at least this many pixels has no maximum left to fit.
const SATURATED_AREA: usize = 6;
/// Position error of a saturated star, pixels: a constant plus a fraction of the radius of its
/// saturated core (measured on synthetic frames; see `tests/synthetic.rs`).
const SATURATED_SIGMA: (f64, f64) = (0.1, 0.03);

/// Half-width of the height window of the local false discovery rate, in σ.
const FDR_WINDOW: f64 = 0.25;

/// Keep the peaks that stand out against their own surroundings, and measure on its local
/// background each peak whose surroundings differ from the mesh.
///
/// The background mesh follows the frame on the scale of its tiles. Stray light, the halo of a
/// bright star or a nebula change faster than that, and what the mesh misses is left in the
/// filtered image as broad bumps that pass a threshold on height above the mesh. A star is
/// narrow: at three times its width from its centre nothing of it is left, so its height above
/// the median of a ring of that radius is its height. A broad bump stands hardly above its own
/// ring.
///
/// Where the ring agrees with the mesh (within twice the ring's own noise) the mesh, which is
/// the less noisy of the two, stands. Elsewhere the ring is the background, and the peak must
/// rise above it by the threshold less one σ: the allowance keeps a faint star whose ring
/// happens to read high, and a bump two or more star-widths across still falls short.
/// Saturated stars are exempt: their flat tops are wider than the ring.
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
            // No ring this close to the frame edge, or one that says what the mesh says.
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

/// Classify the last pass's peaks, trace streaks and build the result.
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
    // Noise of one input pixel, from the filtered noise and the kernel's noise gain (white
    // noise assumed): the filter multiplies the noise by Σk², binning by 1/b.
    let gain = pass.kernel.noise_gain();
    let raw_noise = f64::from(typ_noise) / gain * b;

    // A saturation level down in the background is not one (a header value in other units, a
    // mistaken option): every pixel would count as saturated.
    let floor = f64::from(typ_bg) + 20.0 * f64::from(typ_noise) / gain;
    let saturation = find_saturation(&pass, opt).filter(|level| f64::from(*level) > floor);
    let candidates = pass.peaks.len();
    let shape = measure_shape(&pass, psf_var);
    stand_out(&mut pass, psf_var, saturation);
    let mut classes: Vec<Class> = pass
        .peaks
        .iter()
        .map(|p| classify(p, &pass, psf_var, shape, saturation))
        .collect();

    // Saturated stars, brightest first: the centre of each comes from its saturated region,
    // and the other peaks on that region (its rim, a bleed trail) are part of it.
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
            // The peak each region's star is reported at.
            let mut owner: Vec<Option<usize>> = Vec::new();
            for i in order {
                let p = &pass.peaks[i];
                let Some(k) = regions.at(&pass.raw, p.ix, p.iy, level) else {
                    continue;
                };
                owner.resize(regions.info.len(), None);
                let pl = regions.info[k];
                if pl.linear || pl.truncated {
                    // A saturated streak, left to the streak tracer; or part of a frame
                    // saturated over a wide area, which is no star.
                    classes[i] = Class::Extended;
                } else if owner[k].is_some() {
                    absorbed[i] = true;
                } else {
                    owner[k] = Some(i);
                    plateaus[i] = Some(pl);
                }
            }
            // Every other peak within reach of a star's saturated pixels is part of that star.
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

    // Streaks: trace from every ridge-like peak, brightest first.
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
            // Position error of the maximum of a filtered image: noise gradient over signal
            // curvature, σ = σt² / (√2 · σk · SNR) with σt² the profile's variance after
            // filtering. Equal to the Cramér–Rao bound √2·σpsf / SNR when the kernel matches.
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
                // The star's own photon noise moves the maximum of the filtered image by
                // σp²·σt⁸ / (N·σk⁴·(σk² + 2σp²)²) in variance, N electrons: (16/9)·σp²/N for
                // a matched kernel (derivation in docs/design/detection.md).
                let sp2 = (total2 - k2).max(1.0 / 12.0);
                let electrons = (flux * gain).max(1.0);
                let t8 = total2 * total2 * total2 * total2;
                let d = k2 + 2.0 * sp2;
                sigma = (sigma * sigma + sp2 * t8 / (electrons * k2 * k2 * d * d)).sqrt();
            }
            let (mut px, mut py) = (p.x, p.y);
            if let (Some(pl), Some(level)) = (plateaus[i], saturation) {
                // A top of a few pixels still has a maximum to fit; a wider one does not, and
                // the centre of its saturated region is used. The flux is a lower bound that
                // grows with the region: enough to rank saturated stars among themselves and
                // above the rest.
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
    // Streak indices changed with the sort: map them through.
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

/// The saturation level: the caller's, else the level of flat-topped bright peaks, if any.
fn find_saturation(pass: &Pass, opt: &DetectOptions) -> Option<f32> {
    if let Some(s) = opt.saturation {
        return Some(s);
    }
    // A saturated star has several pixels at (almost) the same value at its top; an
    // unsaturated one, however bright, has one. Saturation is also the highest level a pixel
    // can hold: flat tops well below the brightest peak of the frame (a cosmic-ray track of
    // even charge, a column stuck at a value) are something else.
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

/// The highest unfiltered value within one pixel of the peak, and how many of the 5 × 5 pixels
/// around the peak lie within 0.5% of it (relative to the background).
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

/// Class of a peak from its fitted shape, before streak tracing.
fn classify(
    p: &Peak,
    pass: &Pass,
    psf_var: f64,
    shape: StarShape,
    saturation: Option<f32>,
) -> Class {
    let k2 = pass.kernel.sigma * pass.kernel.sigma;
    if !p.fitted {
        // Too sharp for the fit (a neighbour at the background): a single hot pixel under a
        // narrow kernel; or a noise peak barely above the threshold, which the false discovery
        // rate handles.
        return if p.snr() > 2.0 * pass.threshold {
            Class::Sharp
        } else {
            Class::Star
        };
    }
    // 1σ relative error of a fitted variance, from the noise in a second difference of the
    // filtered image over the peak height (measured on synthetic frames: about 2 / SNR).
    let rel = 2.0 / p.snr();
    let minor = p.minor2 - k2;
    let tol = 2.5 * rel * p.minor2;
    // Narrower across its narrowest axis than the stars are across theirs, by more than the
    // noise allows. Tested before saturation: a hot pixel at the saturation level is still a
    // hot pixel.
    if shape.minor >= RESOLVED && minor < 0.4 * shape.minor - tol {
        return Class::Sharp;
    }
    if sharp_in_raw(p, pass, shape.minor) {
        return Class::Sharp;
    }
    if let Some(level) = saturation {
        if flat_top(pass, p).0 >= level {
            return Class::Saturated;
        }
    }
    // More elongated than the frame's stars are.
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

/// Hot-pixel and cosmic-ray test on the unfiltered image (the idea of DAOFIND's sharpness
/// statistic, Stetson 1987, applied to single neighbours). Light through the optics spreads to
/// all four neighbours of its brightest pixel: for a Gaussian star the faintest of them, the
/// one opposite the direction the star is displaced in, still holds `exp(−1/σ²)` of the peak
/// at the worst displacement of half a pixel (σ along the stars' narrow axis). Charge deposited in the detector does not: a
/// single pixel has no bright neighbour, a track one or two. The test is on the second
/// faintest neighbour, so that noise on one of them decides nothing.
fn sharp_in_raw(p: &Peak, pass: &Pass, psf_minor: f64) -> bool {
    let (w, h) = (pass.w, pass.h);
    if p.ix < 2 || p.iy < 2 || p.ix + 2 >= w || p.iy + 2 >= h {
        return false;
    }
    // The brightest raw pixel within one pixel of the filtered maximum.
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
    // The faintest neighbours of a smeared star lie across the smear: the stars' narrow axis.
    let expected = (-1.0 / psf_minor.max(0.1)).pexp();
    let raw_noise = p.noise / pass.kernel.noise_gain();
    // Half the least a star would show, less twice the noise of a pixel.
    c > 10.0 * raw_noise && around[1] < 0.5 * expected * c - 2.0 * raw_noise
}
