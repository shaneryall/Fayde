use fayde_core::Portable;

use crate::filter::{self, Kernel};
use crate::mesh::{ClippedBackground, Mesh};
use crate::peaks;
use crate::stats;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FrameMode {
    Points,
    Trails,
    Kernel,
    Refused,
}

impl FrameMode {
    pub const POINT_RETRY_MODES: [Self; 2] = [Self::Trails, Self::Kernel];

    pub fn name(self) -> &'static str {
        match self {
            Self::Points => "points",
            Self::Trails => "trails",
            Self::Kernel => "kernel",
            Self::Refused => "refused",
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Census {
    pub mode: FrameMode,
    pub reason: String,
    pub votes: usize,
    pub mode_votes: usize,
    pub mode_share: f64,
    pub second_share: f64,
    pub elongated_share: f64,
    pub shapes: Option<(f64, f64)>,
    pub autocov: Option<(f64, f64)>,
    pub regions: (usize, f64, f64),
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FrameTrail {
    pub length: f64,
    pub angle: f64,
    pub width: f64,
    pub sigma_length: f64,
    pub sigma_angle: f64,
    pub fitted: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Kind {
    Trail,
    Point,
    /// Point-like, but sharper than a star: a cosmic-ray hit.
    Hit,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct TrailSource {
    pub x: f64,
    pub y: f64,
    pub sigma_along: f64,
    pub sigma_across: f64,
    pub flux: f64,
    pub snr: f64,
    pub peak: f64,
    pub p_star: f64,
    pub kind: Kind,
    pub truncated: bool,
}

pub(crate) struct Input<'a> {
    pub raw: &'a [f32],
    pub orig: &'a [f32],
    pub s: &'a [f32],
    pub w: usize,
    pub h: usize,
    pub mesh_tile: usize,
    pub mesh: &'a Mesh,
    pub clipped: Option<&'a ClippedBackground>,
    pub kernel_sigma: f64,
    pub kernel_gain: f64,
    pub mask: Vec<bool>,
    pub peaks: Vec<(f64, f64, f64, f64)>,
    pub psf_var: f64,
    pub saturation: Option<f32>,
    pub gain: Option<f64>,
    pub false_alarms: f64,
    pub threshold: Option<f64>,
    pub sigma_floor: f64,
    pub prior: Option<f64>,
    pub streaks: Vec<([f64; 4], f64, f64, bool, f64, f64)>,
}

pub(crate) struct Output {
    pub census: Census,
    pub trail: Option<FrameTrail>,
    pub kernel: Option<crate::kernel::Output>,
    pub sources: Vec<TrailSource>,
    pub threshold: f64,
    pub kernel_sigma: f64,
    /// Endpoints of streak candidates found to be chains of overlapping star trails.
    pub chains: Vec<[f64; 4]>,
}

const COMPONENT_SIGMA: f32 = 4.0;
const VOTE_AXIS_RATIO: f64 = 1.5;
const OPTICAL_AXIS_RATIO: f64 = 4.0;
const HIT_VARIANCE_MARGIN: f64 = 1.0 / 3.0;
const HIT_LENGTH_WIDTHS: f64 = 12.0;
const VOTE_CAP: f64 = 30.0;
const BRIGHT: f64 = 6.0;
const MIN_MODE_VOTES: usize = 5;
const MIN_MODE_SHARE: f64 = 0.5;
const MIN_ELONGATED_SHARE: f64 = 0.05;
const MAX_SECOND_SHARE: f64 = 0.5;
const MIN_LENGTH_FWHM: f64 = 7.0;
const LENGTH_AGREEMENT: f64 = 0.2;
const CENSUS_AGREEMENT: f64 = 3.0;
const CENSUS_FLOOR: f64 = 0.08;
const ANGLE_AGREEMENT: f64 = 2.0;
const ANGLE_FLOOR: f64 = 0.05;
const REFINE_TRAILS: usize = 15;
const REFINE_CANDIDATES: usize = 300;
const MAX_SEEDS: usize = 5000;
const MAX_ROUND_SHARE: f64 = 4.0;
const AUTOCOV_RATIO: f64 = 2.0;
const DIRECTION_BINS: usize = 36;
const POINT_CHECK_SHARE: f64 = 0.4;
const POINT_CHECK_CHI2: f64 = 2.0;
const ON_STREAK_CHI2: f64 = 1.5;
/// First-pass fits closer than this many trail widths are one source found from two seeds.
const SAME_SOURCE_WIDTHS: f64 = 1.5;
/// A point source whose brightest pixel holds more than this many times the share of its flux
/// that a star's profile puts in one pixel is a hit, not a tracked object.
const HIT_PEAK_SHARE: f64 = 2.5;
/// Streak candidates within this angle of the trails, radians, may be chains of trails
/// ([`measure_trails`]): the chains seen in crowded synthetic fields lay within 10.4°.
const CHAIN_ANGLE: f64 = 12.0 * std::f64::consts::PI / 180.0;
/// Robust spread of a candidate's light along it, over its median, beyond which it is a chain
/// of trails: satellite streaks crossing crowded trailed fields measured 0.02 to 0.13, chains
/// 0.6 to 1.9.
const CHAIN_UNEVENNESS: f64 = 0.5;
const MAX_COMPONENTS: usize = 50_000;
const MAX_COMPONENT_AREA: usize = 40_000;

#[derive(Debug, Clone, Copy)]
pub(crate) struct Component {
    pub(crate) x: f64,
    pub(crate) y: f64,
    pub(crate) major: f64,
    pub(crate) minor: f64,
    pub(crate) theta: f64,
    pub(crate) peak: f64,
    pub(crate) edge: bool,
    pub(crate) masked: bool,
    pub(crate) saturated: bool,
    pub(crate) bounds: [usize; 4],
    pub(crate) first_pixel: usize,
    pub(crate) end_pixel: usize,
}

impl Component {
    pub(crate) fn length(&self) -> f64 {
        (12.0 * (self.major - self.minor)).max(0.0).sqrt()
    }

    pub(crate) fn elongated(&self) -> bool {
        self.major >= VOTE_AXIS_RATIO * VOTE_AXIS_RATIO * self.minor && self.length() >= 3.0
    }

    fn narrow_track(&self, kernel_var: f64) -> bool {
        self.minor <= kernel_var + HIT_VARIANCE_MARGIN
            && self.length() <= HIT_LENGTH_WIDTHS * self.minor.sqrt()
    }

    fn weight(&self) -> f64 {
        self.peak.min(VOTE_CAP)
    }
}

fn axial(a: f64) -> f64 {
    let pi = std::f64::consts::PI;
    let mut a = a % pi;
    if a <= -0.5 * pi {
        a += pi;
    } else if a > 0.5 * pi {
        a -= pi;
    }
    a
}

pub(crate) fn axial_diff(a: f64, b: f64) -> f64 {
    axial(a - b).abs()
}

pub(crate) fn star_like(angle: f64, length: f64, cosmic: bool, t: &FrameTrail) -> bool {
    let tol = (2.0 * t.width / t.length).max(ANGLE_FLOOR);
    axial_diff(angle, t.angle) <= 3.0 * tol && length <= 3.0 * t.length && !cosmic
}

fn axial_mean(items: &[(f64, f64)]) -> Option<f64> {
    let (mut c, mut s) = (0.0, 0.0);
    for &(a, w) in items {
        let (sn, cs) = (2.0 * a).psin_cos();
        c += w * cs;
        s += w * sn;
    }
    (c * c + s * s > 0.0).then(|| axial(0.5 * s.patan2(c)))
}

fn weighted_median(items: &mut [(f64, f64)]) -> Option<f64> {
    if items.is_empty() {
        return None;
    }
    items.sort_by(|a, b| a.0.total_cmp(&b.0));
    let total: f64 = items.iter().map(|v| v.1).sum();
    let mut acc = 0.0;
    for &(v, w) in items.iter() {
        acc += w;
        if acc >= 0.5 * total {
            return Some(v);
        }
    }
    items.last().map(|v| v.0)
}

fn components(inp: &Input, collect_pixels: bool) -> (Vec<Component>, Vec<usize>) {
    let (w, h) = (inp.w, inp.h);
    let mut seen = vec![false; w * h];
    let mut out = Vec::new();
    let mut pixels = Vec::new();
    let mut stack: Vec<usize> = Vec::new();
    let (nx, ny) = inp.mesh.tiles();
    let mut floor = vec![0.0f32; w * h];
    for ty in 0..ny {
        let (y0, y1) = inp.mesh.span(ty, false, h);
        for tx in 0..nx {
            let (x0, x1) = inp.mesh.span(tx, true, w);
            let f = inp.mesh.floor(tx, ty, COMPONENT_SIGMA, w, h);
            for y in y0..y1 {
                floor[y * w + x0..y * w + x1].fill(f);
            }
        }
    }
    let above = |i: usize| -> bool {
        if inp.s[i] <= floor[i] {
            return false;
        }
        let (bg, noise) = inp.mesh.at_px(i % w, i / w);
        inp.s[i] - bg > COMPONENT_SIGMA * noise
    };
    const BLOCK: usize = 64;
    let starts = (0..w * h).step_by(BLOCK).flat_map(|b0| {
        let b1 = (b0 + BLOCK).min(w * h);
        let any = inp.s[b0..b1]
            .iter()
            .zip(&floor[b0..b1])
            .fold(0u32, |a, (v, f)| a + u32::from(v > f));
        if any == 0 {
            b1..b1
        } else {
            b0..b1
        }
    });
    for start in starts {
        if seen[start] || out.len() >= MAX_COMPONENTS {
            continue;
        }
        if !above(start) {
            continue;
        }
        let (mut sw, mut sx, mut sy, mut sxx, mut syy, mut sxy) = (0.0f64, 0.0, 0.0, 0.0, 0.0, 0.0);
        let (mut peak, mut area, mut edge, mut masked, mut sat) =
            (0.0f64, 0usize, false, false, false);
        let mut best = (0.0f64, start);
        let first = pixels.len();
        let mut bounds = [w, h, 0, 0];
        seen[start] = true;
        stack.clear();
        stack.push(start);
        while let Some(i) = stack.pop() {
            let (px, py) = (i % w, i / w);
            let (bg, noise) = inp.mesh.at_px(px, py);
            let v = f64::from(inp.s[i] - bg);
            let z = v / f64::from(noise);
            if z > best.0 {
                best = (z, i);
            }
            peak = peak.max(z);
            area += 1;
            if area <= MAX_COMPONENT_AREA {
                if collect_pixels {
                    pixels.push(i);
                }
                bounds[0] = bounds[0].min(px);
                bounds[1] = bounds[1].min(py);
                bounds[2] = bounds[2].max(px);
                bounds[3] = bounds[3].max(py);
                let (fx, fy) = (px as f64, py as f64);
                sw += v;
                sx += v * fx;
                sy += v * fy;
                sxx += v * fx * fx;
                syy += v * fy * fy;
                sxy += v * fx * fy;
            }
            if px <= 1 || py <= 1 || px + 2 >= w || py + 2 >= h {
                edge = true;
            }
            if inp.saturation.is_some_and(|l| inp.orig[i] >= l) {
                sat = true;
            }
            for dy in -1i64..=1 {
                for dx in -1i64..=1 {
                    let (nx, ny) = (px as i64 + dx, py as i64 + dy);
                    if nx < 0 || ny < 0 || nx >= w as i64 || ny >= h as i64 {
                        continue;
                    }
                    let j = ny as usize * w + nx as usize;
                    if seen[j] {
                        continue;
                    }
                    if above(j) {
                        seen[j] = true;
                        stack.push(j);
                    }
                }
            }
        }
        masked |= inp.mask[best.1];
        if area > MAX_COMPONENT_AREA || sw <= 0.0 {
            pixels.truncate(first);
            continue;
        }
        let (mx, my) = (sx / sw, sy / sw);
        let cxx = sxx / sw - mx * mx + 1.0 / 12.0;
        let cyy = syy / sw - my * my + 1.0 / 12.0;
        let cxy = sxy / sw - mx * my;
        let tr = 0.5 * (cxx + cyy);
        let d = (0.25 * (cxx - cyy) * (cxx - cyy) + cxy * cxy).sqrt();
        out.push(Component {
            x: mx,
            y: my,
            major: tr + d,
            minor: (tr - d).max(1e-6),
            theta: axial(0.5 * (2.0 * cxy).patan2(cxx - cyy)),
            peak,
            edge,
            masked,
            saturated: sat,
            first_pixel: first,
            end_pixel: pixels.len(),
            bounds,
        });
    }
    (out, pixels)
}

fn modes(votes: &[&Component]) -> Vec<(f64, f64, usize)> {
    let nb = DIRECTION_BINS;
    let pi = std::f64::consts::PI;
    let bin_of = |a: f64| (((a + 0.5 * pi) / pi * nb as f64) as usize).min(nb - 1);
    let mut hist = vec![0.0f64; nb];
    let mut count = vec![0usize; nb];
    for c in votes {
        let b = bin_of(c.theta);
        hist[b] += c.weight();
        count[b] += 1;
    }
    let sm: Vec<f64> = (0..nb)
        .map(|i| hist[(i + nb - 1) % nb] + 2.0 * hist[i] + hist[(i + 1) % nb])
        .collect();
    let mut found: Vec<(f64, f64, usize)> = Vec::new();
    for i in 0..nb {
        let (l, r) = (sm[(i + nb - 1) % nb], sm[(i + 1) % nb]);
        if sm[i] > 0.0 && sm[i] > l && sm[i] >= r {
            let span = [(i + nb - 1) % nb, i, (i + 1) % nb];
            let wsum: f64 = span.iter().map(|&k| hist[k]).sum();
            let n: usize = span.iter().map(|&k| count[k]).sum();
            let centre = -0.5 * pi + (i as f64 + 0.5) * pi / nb as f64;
            found.push((centre, wsum, n));
        }
    }
    found.sort_by(|a, b| b.1.total_cmp(&a.1).then(a.0.total_cmp(&b.0)));
    found
}

fn autocovariance(
    inp: &Input,
    l0: f64,
    theta0: f64,
    sigma_t: f64,
    selected: Option<(&[&Component], &[usize])>,
) -> Option<(f64, f64)> {
    let (w, h) = (inp.w, inp.h);
    let mut off = inp.mask.clone();
    let shapes_trail = FrameTrail {
        length: l0,
        angle: theta0,
        width: sigma_t,
        sigma_length: 0.0,
        sigma_angle: 0.0,
        fitted: 0,
    };
    for &(e, angle, length, cosmic, sw, _) in &inp.streaks {
        if star_like(angle, length, cosmic, &shapes_trail) {
            continue;
        }
        let reach = 2.0 * sw.max(1.0) + 2.0;
        let (dx, dy) = (e[2] - e[0], e[3] - e[1]);
        let n = dx.phypot(dy).ceil().max(1.0) as usize;
        for k in 0..=n {
            let f = k as f64 / n as f64;
            let (cx, cy) = (e[0] + f * dx, e[1] + f * dy);
            let (x0, x1) = ((cx - reach).floor().max(0.0), (cx + reach).ceil());
            let (y0, y1) = ((cy - reach).floor().max(0.0), (cy + reach).ceil());
            for yy in y0 as usize..=(y1 as usize).min(h - 1) {
                for xx in x0 as usize..=(x1 as usize).min(w - 1) {
                    off[yy * w + xx] = true;
                }
            }
        }
    }
    if let Some(level) = inp.saturation {
        let r = 3usize;
        for y in 0..h {
            for x in 0..w {
                if inp.orig[y * w + x] >= 0.8 * level {
                    for yy in y.saturating_sub(r)..=(y + r).min(h - 1) {
                        for xx in x.saturating_sub(r)..=(x + r).min(w - 1) {
                            off[yy * w + xx] = true;
                        }
                    }
                }
            }
        }
    }
    let z = |x: f64, y: f64| -> f64 {
        let (ix, iy) = (x.round(), y.round());
        if ix < 0.0 || iy < 0.0 || ix >= w as f64 || iy >= h as f64 {
            return 0.0;
        }
        if off[iy as usize * w + ix as usize] {
            return 0.0;
        }
        match peaks::sample(inp.s, w, h, x, y) {
            Some(v) => {
                let (bg, n) = inp.mesh.at(x, y);
                ((v - f64::from(bg)) / f64::from(n)).max(0.0)
            }
            None => 0.0,
        }
    };
    let mut anchors: Vec<(f64, f64, f64)> = Vec::new();
    if let Some((votes, pixels)) = selected {
        for c in votes.iter().filter(|c| !c.saturated && !c.edge) {
            for &i in &pixels[c.first_pixel..c.end_pixel] {
                if !off[i] {
                    let (x, y) = (i % w, i / w);
                    let (bg, n) = inp.mesh.at_px(x, y);
                    let v = (inp.s[i] - bg) / n;
                    anchors.push((x as f64, y as f64, f64::from(v)));
                }
            }
        }
    } else {
        let (mut bg, mut noise, mut zrow) = (vec![0.0f32; w], vec![0.0f32; w], vec![0.0f32; w]);
        for y in 0..h {
            inp.mesh.row(y, &mut bg, &mut noise);
            let srow = &inp.s[y * w..(y + 1) * w];
            for (((z, s), b), n) in zrow.iter_mut().zip(srow).zip(&bg).zip(&noise) {
                *z = (s - b) / n;
            }
            for (x, &v) in zrow.iter().enumerate() {
                if !off[y * w + x] && v > 3.0 {
                    anchors.push((x as f64, y as f64, f64::from(v)));
                }
            }
        }
    }
    if anchors.len() < 50 {
        return None;
    }
    let stride = (anchors.len() / 5_000).max(1);
    let anchors: Vec<(f64, f64, f64)> = anchors.into_iter().step_by(stride).collect();
    let corr = |a: f64, d: f64| -> f64 {
        let (sn, cs) = a.psin_cos();
        anchors
            .iter()
            .map(|&(x, y, v)| (z(x + d * cs, y + d * sn) / v).min(AUTOCOV_RATIO))
            .sum::<f64>()
    };
    let pi = std::f64::consts::PI;
    let n_dir = 36;
    let half_trail = (0.5 * l0).max(2.0 * sigma_t);
    let mut lag = (4.0 * sigma_t).min(half_trail);
    let parabola = |l: f64, c: f64, r: f64| -> f64 {
        let den = l - 2.0 * c + r;
        if den < 0.0 {
            (0.5 * (l - r) / den).clamp(-0.5, 0.5)
        } else {
            0.0
        }
    };
    let vals: Vec<f64> = (0..n_dir)
        .map(|k| corr(-0.5 * pi + k as f64 * pi / n_dir as f64, lag))
        .collect();
    let (kb, _) = vals
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.total_cmp(b.1).then(b.0.cmp(&a.0)))?;
    let off = parabola(
        vals[(kb + n_dir - 1) % n_dir],
        vals[kb],
        vals[(kb + 1) % n_dir],
    );
    let mut theta = axial(-0.5 * pi + (kb as f64 + off) * pi / n_dir as f64);
    while lag < half_trail {
        lag = (2.0 * lag).min(half_trail);
        let step = 0.25 * sigma_t / lag;
        for _ in 0..2 {
            let vals: Vec<f64> = (-8..=8)
                .map(|k| corr(theta + k as f64 * step, lag))
                .collect();
            let (kb, _) = vals
                .iter()
                .enumerate()
                .max_by(|a, b| a.1.total_cmp(b.1).then(b.0.cmp(&a.0)))?;
            let edge = kb == 0 || kb == vals.len() - 1;
            let off = if edge {
                0.0
            } else {
                parabola(vals[kb - 1], vals[kb], vals[kb + 1])
            };
            theta = axial(theta + (kb as f64 - 8.0 + off) * step);
            if !edge {
                break;
            }
        }
    }
    let step = (l0 / 16.0).max(0.5);
    let n_lag = ((2.0 * l0) / step).ceil() as usize + 1;
    let first = (2.0 * sigma_t / step).ceil() as usize;
    if n_lag <= first + 4 {
        return None;
    }
    let lags: Vec<f64> = (first..n_lag).map(|i| i as f64 * step).collect();
    let prof: Vec<f64> = lags.iter().map(|&d| corr(theta, d)).collect();
    let blur = std::f64::consts::SQRT_2 * sigma_t;
    let taps: Vec<(f64, f64)> = {
        let n = 24;
        (0..=2 * n)
            .map(|k| {
                let t = (k as f64 - n as f64) / n as f64 * 3.0 * blur;
                (t, (-0.5 * t * t / (blur * blur)).pexp())
            })
            .collect()
    };
    let model = |d: f64, len: f64| -> f64 {
        let mut acc = 0.0;
        let mut wsum = 0.0;
        for &(t, g) in &taps {
            acc += g * (1.0 - (d - t).abs() / len).max(0.0);
            wsum += g;
        }
        acc / wsum
    };
    let fit = |len: f64| -> (f64, f64) {
        let (mut s1, mut sm, mut smm, mut sp, mut smp, mut spp) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        for (&d, &p) in lags.iter().zip(&prof) {
            let m = model(d, len);
            s1 += 1.0;
            sm += m;
            smm += m * m;
            sp += p;
            smp += m * p;
            spp += p * p;
        }
        let det = s1 * smm - sm * sm;
        if det.abs() < 1e-12 {
            return (f64::INFINITY, 0.0);
        }
        let a = (s1 * smp - sm * sp) / det;
        let c = (sp - a * sm) / s1;
        let ssr = spp - 2.0 * a * smp - 2.0 * c * sp + a * a * smm + 2.0 * a * c * sm + c * c * s1;
        (ssr, a)
    };
    let mut best: Option<(f64, f64)> = None;
    let mut len = 0.5 * l0;
    while len <= 2.0 * l0 {
        let (ssr, a) = fit(len);
        if a > 0.0 && best.is_none_or(|b| ssr < b.1) {
            best = Some((len, ssr));
        }
        len += 0.25;
    }
    best.map(|(l2, _)| (theta, l2))
}

fn census(
    inp: &Input,
    comps: &[Component],
    pixels: Option<&[usize]>,
    reject_narrow: bool,
) -> (Census, Option<(f64, f64, f64)>) {
    let k2 = inp.kernel_sigma * inp.kernel_sigma;
    let usable: Vec<&Component> = comps.iter().filter(|c| !c.masked).collect();
    let bright_weight: f64 = usable
        .iter()
        .filter(|c| c.peak >= BRIGHT)
        .map(|c| c.weight())
        .sum();
    let optical = |c: &Component| {
        !c.narrow_track(k2)
            && c.major >= OPTICAL_AXIS_RATIO * OPTICAL_AXIS_RATIO * c.minor
            && c.length()
                >= MIN_LENGTH_FWHM * crate::FWHM_PER_SIGMA * (c.minor - k2).max(0.0).sqrt()
    };
    let votes: Vec<&Component> = usable
        .iter()
        .copied()
        .filter(|c| c.peak >= BRIGHT && c.elongated() && (!reject_narrow || optical(c)))
        .collect();
    let elongated_weight: f64 = votes.iter().map(|c| c.weight()).sum();
    let mut cen = Census {
        mode: FrameMode::Points,
        reason: String::new(),
        votes: votes.len(),
        mode_votes: 0,
        mode_share: 0.0,
        second_share: 0.0,
        elongated_share: if bright_weight > 0.0 {
            elongated_weight / bright_weight
        } else {
            0.0
        },
        shapes: None,
        autocov: None,
        regions: (0, 0.0, 0.0),
    };
    let found = modes(&votes);
    let Some(&first) = found.first() else {
        cen.reason = "no elongated sources".into();
        return (cen, None);
    };
    let mut chosen = first;
    let mut prior_matched = false;
    if let Some(pa) = inp.prior {
        if let Some(m) = found
            .iter()
            .find(|m| axial_diff(m.0, pa) <= 3.0 * std::f64::consts::PI / DIRECTION_BINS as f64)
        {
            chosen = *m;
            prior_matched = true;
        }
    }
    cen.mode_votes = chosen.2;
    cen.mode_share = chosen.1 / elongated_weight.max(f64::MIN_POSITIVE);
    let second = found
        .iter()
        .filter(|m| axial_diff(m.0, chosen.0) > 3.0 * std::f64::consts::PI / DIRECTION_BINS as f64)
        .map(|m| (m.1, m.2))
        .fold((0.0f64, 0usize), |a, b| if b.0 > a.0 { b } else { a });
    cen.second_share = second.0 / chosen.1.max(f64::MIN_POSITIVE);
    if chosen.2 < MIN_MODE_VOTES
        || (cen.mode_share < MIN_MODE_SHARE && !prior_matched)
        || cen.elongated_share < MIN_ELONGATED_SHARE
    {
        cen.reason = "no common trail direction".into();
        return (cen, None);
    }
    let near = |c: &Component| {
        axial_diff(c.theta, chosen.0) <= 2.5 * std::f64::consts::PI / DIRECTION_BINS as f64
    };
    let in_mode: Vec<&Component> = votes.iter().copied().filter(|c| near(c)).collect();
    let theta1 = axial_mean(
        &in_mode
            .iter()
            .map(|c| (c.theta, c.weight()))
            .collect::<Vec<_>>(),
    )
    .unwrap_or(chosen.0);
    let mut lens: Vec<(f64, f64)> = in_mode
        .iter()
        .filter(|c| !c.edge && !c.saturated && c.peak >= 3.0 * BRIGHT)
        .map(|c| (c.length(), c.weight()))
        .collect();
    if lens.len() < 3 {
        lens = in_mode
            .iter()
            .filter(|c| !c.edge)
            .map(|c| (c.length(), c.weight()))
            .collect();
    }
    let whole = lens.len();
    let width_var = if pixels.is_some() {
        let mut widths: Vec<_> = in_mode
            .iter()
            .filter(|c| c.peak >= 3.0 * BRIGHT && !c.edge && !c.saturated)
            .map(|c| (c.minor - k2).max(0.05))
            .collect();
        stats::median(&mut widths)
    } else {
        let mut across: Vec<_> = inp
            .peaks
            .iter()
            .filter(|p| p.3 >= 10.0)
            .filter(|p| {
                in_mode
                    .iter()
                    .any(|c| (p.0 - c.x).phypot(p.1 - c.y) <= 0.5 * c.length() + 2.0)
            })
            .map(|p| (p.2 - k2).max(0.05))
            .collect();
        stats::median(&mut across)
    }
    .unwrap_or(inp.psf_var);
    let width = width_var.sqrt();
    let sigma_t = (width_var + k2).sqrt();
    let star_var = width_var + k2;
    let round_stars = usable
        .iter()
        .filter(|c| c.peak >= BRIGHT && !c.elongated() && !c.saturated)
        .filter(|c| c.minor >= 0.6 * star_var && c.major <= 2.0 * star_var)
        .count();
    let Some(l1) = weighted_median(&mut lens) else {
        cen.mode = FrameMode::Refused;
        cen.reason = "trails longer than the frame".into();
        return (cen, None);
    };
    cen.shapes = Some((theta1, l1));
    let fwhm = crate::FWHM_PER_SIGMA * width;
    if l1 < MIN_LENGTH_FWHM * fwhm {
        cen.reason = format!("trails shorter than {MIN_LENGTH_FWHM} FWHM");
        return (cen, None);
    }
    cen.mode = FrameMode::Refused;
    if round_stars as f64 > (MAX_ROUND_SHARE * in_mode.len() as f64).max(3.0) {
        cen.reason = format!(
            "point stars as well as trails ({round_stars} round sources, {} trails)",
            in_mode.len()
        );
        return (cen, None);
    }
    if second.1 >= MIN_MODE_VOTES && cen.second_share >= MAX_SECOND_SHARE && !prior_matched {
        cen.reason = "two trail populations".into();
        return (cen, None);
    }
    if whole < 3 && in_mode.iter().filter(|c| c.edge).count() * 4 >= in_mode.len() * 3 {
        cen.reason = "trails longer than the frame".into();
        return (cen, None);
    }
    let selected = pixels.map(|pixels| (&votes[..], pixels));
    let Some((theta2, l2)) = autocovariance(inp, l1, theta1, sigma_t, selected) else {
        cen.reason = "too little light for the second estimate".into();
        return (cen, None);
    };
    cen.autocov = Some((theta2, l2));
    let tol = (CENSUS_AGREEMENT * sigma_t / l1).max(CENSUS_FLOOR);
    if axial_diff(theta1, theta2) > tol || (l2 / l1 - 1.0).abs() > LENGTH_AGREEMENT {
        cen.reason = "the two estimates disagree".into();
        return (cen, None);
    }
    let (mut measured, mut worst_a, mut worst_l) = (0usize, 0.0f64, 0.0f64);
    for ry in 0..3 {
        for rx in 0..3 {
            let inside = |c: &&&Component| {
                let (fx, fy) = (c.x / inp.w as f64 * 3.0, c.y / inp.h as f64 * 3.0);
                fx as usize == rx && fy as usize == ry
            };
            let here: Vec<&&Component> = in_mode.iter().filter(inside).collect();
            if here.len() < 4 {
                continue;
            }
            measured += 1;
            let a = axial_mean(
                &here
                    .iter()
                    .map(|c| (c.theta, c.weight()))
                    .collect::<Vec<_>>(),
            )
            .unwrap_or(theta1);
            let mut ls: Vec<(f64, f64)> = here
                .iter()
                .filter(|c| !c.edge && !c.saturated && c.peak >= 3.0 * BRIGHT)
                .map(|c| (c.length(), c.weight()))
                .collect();
            worst_a = worst_a.max(axial_diff(a, theta1));
            if ls.len() >= 5 {
                if let Some(lr) = weighted_median(&mut ls) {
                    worst_l = worst_l.max((lr / l1 - 1.0).abs());
                }
            }
        }
    }
    cen.regions = (measured, worst_a, worst_l);
    if worst_a > 2.0 * tol || worst_l > 0.5 {
        cen.reason = "the trails change across the frame".into();
        return (cen, None);
    }
    cen.mode = FrameMode::Trails;
    cen.reason = "one straight trail".into();
    (cen, Some((theta1, l1, width)))
}

struct Tables {
    gauss: Vec<f64>,
    erf: Vec<f64>,
}

const TABLE_STEP: f64 = 1.0 / 512.0;
const GAUSS_REACH: f64 = 9.0;
const ERF_REACH: f64 = 6.0;

fn tables() -> &'static Tables {
    static T: std::sync::OnceLock<Tables> = std::sync::OnceLock::new();
    T.get_or_init(|| {
        let n = |r: f64| (r / TABLE_STEP).ceil() as usize + 2;
        Tables {
            gauss: (0..n(GAUSS_REACH))
                .map(|i| {
                    let t = i as f64 * TABLE_STEP;
                    (-0.5 * t * t).pexp()
                })
                .collect(),
            erf: (0..n(ERF_REACH))
                .map(|i| 1.0 - (i as f64 * TABLE_STEP).perfc())
                .collect(),
        }
    })
}

#[inline]
fn lookup(table: &[f64], t: f64, beyond: f64) -> f64 {
    let x = t.abs() / TABLE_STEP;
    let i = x as usize;
    if i + 1 >= table.len() {
        return beyond;
    }
    let f = x - i as f64;
    table[i] + f * (table[i + 1] - table[i])
}

#[inline]
fn gauss(t: f64) -> f64 {
    lookup(&tables().gauss, t, 0.0)
}

#[inline]
fn erf_t(t: f64) -> f64 {
    let v = lookup(&tables().erf, t, 1.0);
    if t < 0.0 {
        -v
    } else {
        v
    }
}

#[inline(always)]
fn profile(u: f64, v: f64, len: f64, s: f64) -> [f64; 5] {
    let sqrt_2pi = (2.0 * std::f64::consts::PI).sqrt();
    let g = gauss(v / s) / (sqrt_2pi * s);
    let dg_dv = -v / (s * s) * g;
    let dg_ds = g * (v * v / (s * s * s) - 1.0 / s);
    if len < 1e-6 {
        let a = gauss(u / s) / (sqrt_2pi * s);
        let da_du = -u / (s * s) * a;
        let da_ds = a * (u * u / (s * s * s) - 1.0 / s);
        return [a * g, da_du * g, a * dg_dv, 0.0, da_ds * g + a * dg_ds];
    }
    let k = 1.0 / (std::f64::consts::SQRT_2 * s);
    let (p, m) = (u + 0.5 * len, u - 0.5 * len);
    let r2 = std::f64::consts::SQRT_2;
    let (ea, eb) = (gauss(r2 * p * k), gauss(r2 * m * k));
    let a = (erf_t(p * k) - erf_t(m * k)) / (2.0 * len);
    let two_over_sqrt_pi = 2.0 / std::f64::consts::PI.sqrt();
    let da_du = two_over_sqrt_pi * k * (ea - eb) / (2.0 * len);
    let da_dl = two_over_sqrt_pi * k * 0.5 * (ea + eb) / (2.0 * len) - a / len;
    let da_ds = two_over_sqrt_pi * (-k / s) * (ea * p - eb * m) / (2.0 * len);
    [
        a * g,
        da_du * g,
        a * dg_dv,
        da_dl * g,
        da_ds * g + a * dg_ds,
    ]
}

#[derive(Debug, Clone, Copy)]
struct Fit {
    x: f64,
    y: f64,
    flux: f64,
    sigma_flux: f64,
    cov: [f64; 3],
    chi2: f64,
    dof: usize,
    len: f64,
    theta: f64,
    width: f64,
    converged: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Free {
    Held,
    Geometry,
}

struct Pixels<'a> {
    raw: &'a [f32],
    orig: &'a [f32],
    extra: Option<&'a [f32]>,
    w: usize,
    h: usize,
    mesh: &'a Mesh,
    kernel_gain: f64,
    mask: &'a [bool],
    saturation: Option<f32>,
    gain: Option<f64>,
    clipped: bool,
}

impl Pixels<'_> {
    fn var(&self, x: usize, y: usize, signal: f64) -> f64 {
        let (_, n) = self.mesh.at(x as f64, y as f64);
        let pix = f64::from(n) / self.kernel_gain;
        let removed = self.extra.map_or(0.0, |e| f64::from(e[y * self.w + x]));
        pix * pix
            + self
                .gain
                .map_or(0.0, |g| (signal.max(0.0) + removed.max(0.0)) / g)
    }
}

fn cholesky_solve(a: &[f64], b: &[f64], n: usize) -> Option<(Vec<f64>, Vec<f64>)> {
    let mut l = vec![0.0f64; n * n];
    for i in 0..n {
        for j in 0..=i {
            let mut sum = a[i * n + j];
            for k in 0..j {
                sum -= l[i * n + k] * l[j * n + k];
            }
            if i == j {
                if sum <= 0.0 {
                    return None;
                }
                l[i * n + i] = sum.sqrt();
            } else {
                l[i * n + j] = sum / l[j * n + j];
            }
        }
    }
    let solve = |rhs: &[f64]| -> Vec<f64> {
        let mut y = vec![0.0f64; n];
        for i in 0..n {
            let mut s = rhs[i];
            for k in 0..i {
                s -= l[i * n + k] * y[k];
            }
            y[i] = s / l[i * n + i];
        }
        let mut x = vec![0.0f64; n];
        for i in (0..n).rev() {
            let mut s = y[i];
            for k in i + 1..n {
                s -= l[k * n + i] * x[k];
            }
            x[i] = s / l[i * n + i];
        }
        x
    };
    let x = solve(b);
    let mut inv = vec![0.0f64; n * n];
    for c in 0..n {
        let mut e = vec![0.0f64; n];
        e[c] = 1.0;
        let col = solve(&e);
        for r in 0..n {
            inv[r * n + c] = col[r];
        }
    }
    Some((x, inv))
}

#[derive(Debug, Clone, Copy)]
struct StampPx {
    xf: f64,
    yf: f64,
    bx: f64,
    by: f64,
    d: f64,
    sub: f64,
    pix2: f64,
    removed: f64,
}

#[inline(always)]
fn normal_equations<const NP: usize>(
    pix: &[StampPx],
    p: &[f64; 9],
    (sn, cs): (f64, f64),
    s: f64,
    gain: Option<f64>,
    ata_out: &mut [f64],
    atb_out: &mut [f64],
) -> f64 {
    let mut ata = [[0.0f64; NP]; NP];
    let mut atb = [0.0f64; NP];
    let mut chi2 = 0.0;
    let mut jrow = [0.0f64; NP];
    for q in pix {
        let (dx, dy) = (q.xf - p[0], q.yf - p[1]);
        let (u, v) = (dx * cs + dy * sn, -dx * sn + dy * cs);
        let [t, tu, tv, tl, ts] = profile(u, v, p[6], s);
        let m = p[3] + p[4] * q.bx + p[5] * q.by + p[2] * t;
        let signal = m - p[3] + q.sub;
        let wt = 1.0 / (q.pix2 + gain.map_or(0.0, |g| (signal.max(0.0) + q.removed) / g));
        let r = q.d - m;
        chi2 += wt * r * r;
        jrow[0] = p[2] * (-cs * tu + sn * tv);
        jrow[1] = p[2] * (-sn * tu - cs * tv);
        jrow[2] = t;
        jrow[3] = 1.0;
        jrow[4] = q.bx;
        jrow[5] = q.by;
        if NP == 9 {
            jrow[6 % NP] = p[2] * tl;
            jrow[7 % NP] = p[2] * (tu * v - tv * u);
            jrow[8 % NP] = p[2] * ts * p[8] / s;
        }
        for a in 0..NP {
            let wa = wt * jrow[a];
            atb[a] += wa * r;
            for (sum, jb) in ata[a].iter_mut().zip(&jrow).take(a + 1) {
                *sum += wa * jb;
            }
        }
    }
    for a in 0..NP {
        atb_out[a] = atb[a];
        ata_out[a * NP..a * NP + a + 1].copy_from_slice(&ata[a][..=a]);
    }
    chi2
}

#[derive(Debug, Clone, Copy)]
struct Spec {
    x0: f64,
    y0: f64,
    len: f64,
    theta: f64,
    width: f64,
    free: Free,
    stamp_len: f64,
    stamp_width: f64,
    max_iter: usize,
}

type Neighbour = (f64, f64, f64, f64);

fn fit(px: &Pixels, spec: &Spec, others: &[Neighbour]) -> Option<Fit> {
    let Spec {
        x0,
        y0,
        len: len0,
        theta: theta0,
        width: width0,
        free,
        stamp_len,
        stamp_width,
        max_iter,
    } = *spec;
    let (w, h) = (px.w, px.h);
    let (sn0, cs0) = theta0.psin_cos();
    let half_u = 0.5 * stamp_len + 3.0 * stamp_width + 3.0;
    let half_v = 3.0 * stamp_width + 1.5;
    let reach = half_u.max(half_v).ceil() as i64 + 1;
    let mut pix: Vec<StampPx> = Vec::new();
    let (cx, cy) = (x0.round() as i64, y0.round() as i64);
    for yy in (cy - reach).max(0)..=(cy + reach).min(h as i64 - 1) {
        let dy = yy as f64 - y0;
        let (mut lo, mut hi) = (f64::NEG_INFINITY, f64::INFINITY);
        let mut clip = |c: f64, centre: f64, half: f64| {
            if c.abs() > 1e-9 {
                let (a, b) = ((centre - half) / c, (centre + half) / c);
                lo = lo.max(a.min(b));
                hi = hi.min(a.max(b));
            }
        };
        clip(cs0, -dy * sn0, half_u);
        clip(-sn0, -dy * cs0, half_v);
        let xa = (cx - reach).max(0).max((x0 + lo).floor() as i64 - 1);
        let xb = (cx + reach)
            .min(w as i64 - 1)
            .min((x0 + hi).ceil() as i64 + 1);
        for xx in xa..=xb {
            let (dx, dy) = (xx as f64 - x0, yy as f64 - y0);
            let (u, v) = (dx * cs0 + dy * sn0, -dx * sn0 + dy * cs0);
            if u.abs() > half_u || v.abs() > half_v {
                continue;
            }
            let i = yy as usize * w + xx as usize;
            if px.mask[i] || px.saturation.is_some_and(|l| px.orig[i] >= l) {
                continue;
            }
            let mut d = f64::from(px.raw[i]);
            let mut sub = 0.0;
            let s_n = (width0 * width0 + 1.0 / 12.0).sqrt();
            for &(ox, oy, of, ol) in others {
                let (ex, ey) = (xx as f64 - ox, yy as f64 - oy);
                let (u, v) = (ex * cs0 + ey * sn0, -ex * sn0 + ey * cs0);
                if u.abs() < 0.5 * ol + 5.0 * s_n && v.abs() < 5.0 * s_n {
                    sub += of * profile(u, v, ol, s_n)[0];
                }
            }
            d -= sub;
            if !d.is_finite() {
                continue;
            }
            let (_, noise) = px.mesh.at_px(xx as usize, yy as usize);
            let pn = f64::from(noise) / px.kernel_gain;
            let removed = px.extra.map_or(0.0, |e| f64::from(e[i]));
            pix.push(StampPx {
                xf: xx as f64,
                yf: yy as f64,
                bx: xx as f64 - x0,
                by: yy as f64 - y0,
                d,
                sub,
                pix2: pn * pn,
                removed: removed.max(0.0),
            });
        }
    }
    let np = if free == Free::Geometry { 9 } else { 6 };
    if pix.len() < 3 * np {
        return None;
    }
    let (bg0, noise0) = px.mesh.at(x0, y0);
    let background_weight = if px.clipped {
        let sigma = 0.25 * f64::from(noise0) / px.kernel_gain;
        1.0 / (sigma * sigma)
    } else {
        0.0
    };
    let mut p = [x0, y0, 0.0, f64::from(bg0), 0.0, 0.0, len0, theta0, width0];
    p[2] = pix.iter().map(|q| q.d - p[3]).sum::<f64>().max(1.0);
    let s_of = |wd: f64| (wd * wd + 1.0 / 12.0).sqrt();
    let mut converged = false;
    let mut settled = false;
    let mut last: Option<(Vec<f64>, f64)> = None;
    let mut damp = 1.0f64;
    let mut prev_step = (0.0f64, 0.0f64);
    for it in 0..=max_iter {
        let (sn, cs) = p[7].psin_cos();
        let s = s_of(p[8]);
        let mut ata = vec![0.0f64; np * np];
        let mut atb = vec![0.0f64; np];
        let mut chi2 = if np == 9 {
            normal_equations::<9>(&pix, &p, (sn, cs), s, px.gain, &mut ata, &mut atb)
        } else {
            normal_equations::<6>(&pix, &p, (sn, cs), s, px.gain, &mut ata, &mut atb)
        };
        if px.clipped {
            for (k, centre, weight) in [
                (3, f64::from(bg0), background_weight),
                (4, 0.0, background_weight * half_u * half_u),
                (5, 0.0, background_weight * half_v * half_v),
            ] {
                let r = centre - p[k];
                chi2 += weight * r * r;
                ata[k * np + k] += weight;
                atb[k] += weight * r;
            }
        }
        for a in 0..np {
            for b in 0..a {
                ata[b * np + a] = ata[a * np + b];
            }
        }
        let (step, inv) = cholesky_solve(&ata, &atb, np)?;
        let pos_err = (inv[0] + inv[np + 1]).max(0.0).sqrt();
        let (angle_err, len_err, width_err) = if np == 9 {
            (
                inv[7 * np + 7].max(0.0).sqrt(),
                inv[6 * np + 6].max(0.0).sqrt(),
                inv[8 * np + 8].max(0.0).sqrt(),
            )
        } else {
            (0.0, 0.0, 0.0)
        };
        last = Some((inv, chi2));
        if converged || it == max_iter {
            break;
        }
        if step[0] * prev_step.0 + step[1] * prev_step.1 < 0.0 {
            damp *= 0.5;
        }
        let mut scale: f64 = damp;
        let pos = step[0].phypot(step[1]);
        let max_pos = (0.1 * p[6]).max(1.0);
        if pos > max_pos {
            scale = scale.min(max_pos / pos);
        }
        prev_step = (scale * step[0], scale * step[1]);
        if np == 9 {
            if step[7].abs() > 0.05 {
                scale = scale.min(0.05 / step[7].abs());
            }
            if step[6].abs() > 0.2 * p[6].max(1.0) {
                scale = scale.min(0.2 * p[6].max(1.0) / step[6].abs());
            }
            if step[8].abs() > 0.3 * p[8] {
                scale = scale.min(0.3 * p[8] / step[8].abs());
            }
        }
        for k in 0..np {
            p[k] += scale * step[k];
        }
        if np == 9 {
            p[6] = p[6].max(0.5);
            p[7] = axial(p[7]);
            p[8] = p[8].clamp(0.3, 20.0);
        }
        if (p[0] - x0).abs() > half_u || (p[1] - y0).abs() > half_u {
            return None;
        }
        let geometry_done = np == 6
            || (scale * step[7].abs() < (0.1 * angle_err).clamp(1e-5, 1e-3)
                && scale * step[6].abs() < (0.1 * len_err).clamp(1e-3, 0.2)
                && scale * step[8].abs() < (0.1 * width_err).clamp(1e-4, 0.01));
        if scale * pos < (0.1 * pos_err).clamp(1e-3, 0.05) && geometry_done {
            converged = true;
        }
        settled = converged
            || (scale * pos < 0.5 * pos_err
                && (np == 6
                    || (scale * step[7].abs() < 0.5 * angle_err
                        && scale * step[6].abs() < 0.5 * len_err
                        && scale * step[8].abs() < 0.5 * width_err)));
    }
    let (inv, chi2) = last?;
    if !(p.iter().all(|v| v.is_finite()) && chi2.is_finite() && inv.iter().all(|v| v.is_finite())) {
        return None;
    }
    Some(Fit {
        x: p[0],
        y: p[1],
        flux: p[2],
        sigma_flux: inv[2 * np + 2].max(0.0).sqrt(),
        cov: [inv[0], inv[1], inv[np + 1]],
        chi2,
        dof: pix.len() - np,
        len: p[6],
        theta: p[7],
        width: p[8],
        converged: settled,
    })
}

fn refine(comps: &[Component], px: &Pixels, start: (f64, f64, f64)) -> Option<FrameTrail> {
    let (theta1, l1, width1) = start;
    let mut cand: Vec<&Component> = comps
        .iter()
        .filter(|c| !c.masked && !c.edge && !c.saturated && c.elongated())
        .filter(|c| axial_diff(c.theta, theta1) <= 0.15 && (c.length() / l1 - 1.0).abs() < 0.4)
        .collect();
    cand.sort_by(|a, b| b.peak.total_cmp(&a.peak).then(a.y.total_cmp(&b.y)));
    cand.truncate(REFINE_CANDIDATES);
    let isolated = |c: &Component| {
        !comps
            .iter()
            .any(|o| o.peak >= BRIGHT && !std::ptr::eq(o, c) && (o.x - c.x).phypot(o.y - c.y) < l1)
    };
    let mut ls: Vec<f64> = Vec::new();
    let mut ts: Vec<f64> = Vec::new();
    let mut ws: Vec<f64> = Vec::new();
    let mut attempted = 0;
    let mut shortened = 0;
    for c in cand.iter().filter(|c| isolated(c)).take(REFINE_TRAILS) {
        attempted += 1;
        let spec = Spec {
            x0: c.x,
            y0: c.y,
            len: l1,
            theta: theta1,
            width: width1,
            free: Free::Geometry,
            stamp_len: 1.3 * l1,
            stamp_width: width1,
            max_iter: 30,
        };
        let Some(f) = fit(px, &spec, &[]) else {
            continue;
        };
        if f.flux > 0.0 && f.len < (1.0 - LENGTH_AGREEMENT) * l1.min(c.length()) {
            shortened += 1;
        }
        if !f.converged || f.chi2 > 4.0 * f.dof as f64 || f.flux <= 0.0 {
            continue;
        }
        ls.push(f.len);
        ts.push(theta1 + axial(f.theta - theta1));
        ws.push(f.width);
    }
    if shortened >= 3 && 2 * shortened > attempted {
        return None;
    }
    let n = ls.len();
    if n == 0 {
        return Some(FrameTrail {
            length: l1,
            angle: theta1,
            width: width1,
            sigma_length: 0.1 * l1,
            sigma_angle: width1 / l1,
            fitted: 0,
        });
    }
    let spread = |v: &mut Vec<f64>, m: f64, floor: f64| -> f64 {
        if n < 3 {
            return floor;
        }
        let mut d: Vec<f64> = v.iter().map(|x| (x - m).abs()).collect();
        1.4826 * stats::median(&mut d).unwrap_or(0.0) * 1.2533 / (n as f64).sqrt()
    };
    let lm = stats::median(&mut ls.clone()).unwrap_or(l1);
    let tm = stats::median(&mut ts.clone()).unwrap_or(theta1);
    let wm = stats::median(&mut ws.clone()).unwrap_or(width1);
    Some(FrameTrail {
        length: lm,
        angle: axial(tm),
        width: wm,
        sigma_length: spread(&mut ls, lm, 0.02 * lm),
        sigma_angle: spread(&mut ts, tm, wm / lm),
        fitted: n,
    })
}

pub(crate) fn aimed_box(img: &[f32], w: usize, h: usize, theta: f64, len: f64) -> Vec<f32> {
    let (sn, cs) = theta.psin_cos();
    if cs.abs() >= sn.abs() {
        line_sums(img, w, h, sn / cs, 0.5 * len * cs.abs(), 1.0 / cs.abs())
    } else {
        let t = transpose(img, w, h);
        let s = line_sums(&t, h, w, cs / sn, 0.5 * len * sn.abs(), 1.0 / sn.abs());
        transpose(&s, h, w)
    }
}

fn transpose(img: &[f32], w: usize, h: usize) -> Vec<f32> {
    const B: usize = 32;
    let mut t = vec![0.0f32; w * h];
    for y0 in (0..h).step_by(B) {
        for x0 in (0..w).step_by(B) {
            for y in y0..(y0 + B).min(h) {
                for x in x0..(x0 + B).min(w) {
                    t[x * h + y] = img[y * w + x];
                }
            }
        }
    }
    t
}

fn line_sums(img: &[f32], w: usize, h: usize, slope: f64, half: f64, ds: f64) -> Vec<f32> {
    let span = slope * (w as f64 - 1.0);
    let c0 = (-span.max(0.0)).floor() - 1.0;
    let c_max = (h as f64 - 1.0) - span.min(0.0);
    let nc = (c_max - c0).ceil() as usize + 2;
    let (n, fr, centre) = if half >= 0.5 {
        let n = (half - 0.5).floor();
        (n as usize, half - 0.5 - n, 1.0)
    } else {
        (0, 0.0, 2.0 * half)
    };
    let reach = n + 1;
    let active = |c: f64| -> (usize, usize) {
        if slope == 0.0 {
            return if c > -1.0 && c < h as f64 {
                (0, w)
            } else {
                (0, 0)
            };
        }
        let (p, q) = ((-1.0 - c) / slope, (h as f64 - c) / slope);
        let (lo, hi) = (p.min(q), p.max(q));
        let a = lo.floor().clamp(0.0, w as f64) as usize;
        let b = (hi.ceil() + 1.0).clamp(0.0, w as f64) as usize;
        (a, b.max(a))
    };
    let sample = |x: usize, yf: f64| -> f64 {
        let y0 = yf.floor();
        let f = yf - y0;
        let y0 = y0 as i64;
        if y0 >= 0 && y0 + 1 < h as i64 {
            let i = y0 as usize * w + x;
            return (1.0 - f) * f64::from(img[i]) + f * f64::from(img[i + w]);
        }
        let get = |y: i64| -> f64 {
            if y < 0 || y >= h as i64 {
                0.0
            } else {
                f64::from(img[y as usize * w + x])
            }
        };
        (1.0 - f) * get(y0) + f * get(y0 + 1)
    };
    let mut vals = vec![0.0f64; w];
    let mut pre = vec![0.0f64; w + 1];
    let mut row = |ci: usize, sums: &mut [f64], old: (usize, usize)| -> (usize, usize) {
        sums[old.0..old.1].fill(0.0);
        let c = c0 + ci as f64;
        let (a, b) = active(c);
        if a >= b {
            return (0, 0);
        }
        pre[a] = 0.0;
        let mut acc = 0.0f64;
        for x in a..b {
            let v = sample(x, c + slope * x as f64);
            vals[x] = v;
            acc += v;
            pre[x + 1] = acc;
        }
        let val = |x: usize| if x >= a && x < b { vals[x] } else { 0.0 };
        let cum = |x: usize| pre[x.clamp(a, b)];
        let (ea, eb) = (a.saturating_sub(reach), (b + reach).min(w));
        let (lo, hi) = if centre < 1.0 {
            (eb, eb)
        } else {
            let lo = (a + n + 1).clamp(ea, eb);
            (lo, b.saturating_sub(n + 1).clamp(lo, eb))
        };
        for x in lo..hi {
            let mut v = pre[x + n + 1] - pre[x - n];
            if fr > 0.0 {
                v += fr * vals[x - n - 1];
                v += fr * vals[x + n + 1];
            }
            sums[x] = v * ds;
        }
        for x in (ea..lo).chain(hi..eb) {
            let out = &mut sums[x];
            let mut v = if centre < 1.0 {
                centre * val(x)
            } else {
                cum((x + n + 1).min(w)) - cum(x.saturating_sub(n))
            };
            if fr > 0.0 {
                if x > n {
                    v += fr * val(x - n - 1);
                }
                if x + n + 1 < w {
                    v += fr * val(x + n + 1);
                }
            }
            *out = v * ds;
        }
        (ea, eb)
    };
    let mut out = vec![0.0f32; w * h];
    let mut prev = vec![0.0f64; w];
    let mut cur = vec![0.0f64; w];
    let mut prev_range = row(0, &mut prev, (0, 0));
    let mut cur_range = (0, 0);
    for k in 0..nc.saturating_sub(1) {
        cur_range = row(k + 1, &mut cur, cur_range);
        let (a0, b0) = active(c0 + k as f64);
        let (a1, b1) = active(c0 + (k + 1) as f64);
        let (xa, xb) = if a0 >= b0 {
            (a1, b1)
        } else if a1 >= b1 {
            (a0, b0)
        } else {
            (a0.min(a1), b0.max(b1))
        };
        for x in xa..xb {
            let base = (k as f64 + c0 + slope * x as f64).ceil() as i64;
            if base >= 0 && base < h as i64 {
                let g = base as f64 - slope * x as f64 - c0 - k as f64;
                if (1e-6..=1.0 - 1e-6).contains(&g) {
                    out[base as usize * w + x] = ((1.0 - g) * prev[x] + g * cur[x]) as f32;
                    continue;
                }
            }
            for y in base - 1..=base + 1 {
                if y < 0 || y >= h as i64 {
                    continue;
                }
                let c = y as f64 - slope * x as f64 - c0;
                if c.floor() as i64 != k as i64 {
                    continue;
                }
                let g = c - k as f64;
                out[y as usize * w + x] = ((1.0 - g) * prev[x] + g * cur[x]) as f32;
            }
        }
        std::mem::swap(&mut prev, &mut cur);
        std::mem::swap(&mut prev_range, &mut cur_range);
    }
    out
}

pub(crate) fn equivalent_sigma(len: f64, sk: f64) -> f64 {
    let d = 0.25;
    let sqrt2 = std::f64::consts::SQRT_2;
    let erf = |x: f64| 1.0 - x.perfc();
    let along = |u: f64| {
        if len < 1e-6 {
            (-0.5 * u * u / (sk * sk)).pexp()
        } else {
            erf((u + 0.5 * len) / (sqrt2 * sk)) - erf((u - 0.5 * len) / (sqrt2 * sk))
        }
    };
    let across = |v: f64| (-0.5 * v * v / (sk * sk)).pexp();
    let ratio = |f: &dyn Fn(f64) -> f64, reach: f64| -> f64 {
        let n = (reach / d).ceil() as i64;
        let (mut s2, mut g2) = (0.0, 0.0);
        for i in -n..=n {
            let t = i as f64 * d;
            let v = f(t);
            let g = (f(t + 0.5 * d) - f(t - 0.5 * d)) / d;
            s2 += v * v;
            g2 += g * g;
        }
        g2 / s2
    };
    let lu = ratio(&along, 0.5 * len + 6.0 * sk);
    let lv = ratio(&across, 6.0 * sk);
    (1.0 / (2.0 * (lu * lv).sqrt())).sqrt()
}

fn crowded_by_hits(comps: &[Component], kernel_sigma: f64) -> bool {
    let kernel_var = kernel_sigma * kernel_sigma;
    let mut narrow = Vec::new();
    let mut total = 0;
    for c in comps
        .iter()
        .filter(|c| c.peak >= 3.0 * BRIGHT && c.elongated() && !c.saturated)
    {
        total += 1;
        if c.narrow_track(kernel_var) {
            narrow.push(c);
        }
    }
    narrow.len() >= MIN_MODE_VOTES
        && 5 * narrow.len() >= total
        && modes(&narrow)
            .first()
            .is_some_and(|mode| 2 * mode.2 <= narrow.len())
}

fn mask_hits(raw: &mut [f32], w: usize, h: usize, mesh: &Mesh, mask: &mut [bool]) {
    let mut hits = vec![false; w * h];
    for y in 1..h - 1 {
        for x in 1..w - 1 {
            let i = y * w + x;
            let v = raw[i];
            let n = mesh.at_px(x, y).1;
            if v <= 8.0 * n {
                continue;
            }
            let mut neighbors = [raw[i - 1], raw[i + 1], raw[i - w], raw[i + w]];
            neighbors.sort_unstable_by(f32::total_cmp);
            hits[i] = neighbors[1] < 0.25 * v - 1.5 * n;
        }
    }
    for y in 1..h - 1 {
        for x in 1..w - 1 {
            let i = y * w + x;
            let n = mesh.at_px(x, y).1;
            if hits[i]
                || raw[i] > 4.0 * n && (hits[i - 1] || hits[i + 1] || hits[i - w] || hits[i + w])
            {
                mask[i] = true;
                raw[i] = 0.0;
            }
        }
    }
}

fn clean_census(
    inp: &Input,
    valid: Option<&[bool]>,
    cosmic: bool,
    original_regions: usize,
) -> Option<Output> {
    let (w, h) = (inp.w, inp.h);
    let background = Mesh::build(inp.raw, w, h, inp.mesh_tile.min(16), 3, valid);
    let mut raw: Vec<f32> = inp
        .raw
        .iter()
        .enumerate()
        .map(|(i, &v)| {
            if valid.is_none_or(|v| v[i]) {
                v - background.at_px(i % w, i / w).0
            } else {
                0.0
            }
        })
        .collect();
    let mut mask = inp.mask.clone();
    if cosmic {
        for (i, m) in mask.iter_mut().enumerate() {
            *m = valid.is_some_and(|v| !v[i]);
        }
        mask_hits(&mut raw, w, h, &background, &mut mask);
    }
    let frequencies: Vec<_> = crate::pattern::find(&[(&raw, w, h)])
        .into_iter()
        .filter(|p| {
            p.share >= 0.02
                && p.frequency[0].phypot(p.frequency[1])
                    >= 2.0 / crate::pattern::PATTERN_TILE as f64
        })
        .map(|p| p.frequency)
        .collect();
    crate::pattern::remove(&mut raw, w, h, &frequencies);
    let kernel = Kernel::gaussian(if cosmic {
        inp.kernel_sigma.max(1.0)
    } else {
        inp.kernel_sigma
    });
    let mut filtered = Vec::new();
    filter::smooth(&raw, w, h, &kernel, &mut Vec::new(), &mut filtered);
    let mesh = Mesh::build(
        &filtered,
        w,
        h,
        inp.mesh_tile.min(24),
        kernel.radius + 1,
        valid,
    );
    let streaks = inp
        .streaks
        .iter()
        .copied()
        .filter(|&(e, _, len, hit, _, _)| {
            if !cosmic {
                return true;
            }
            if hit {
                return false;
            }
            let n = len.ceil().max(1.0) as usize;
            let supported = (0..=n)
                .filter(|&i| {
                    let t = i as f64 / n as f64;
                    let x = (e[0] + t * (e[2] - e[0]))
                        .round()
                        .clamp(0.0, (w - 1) as f64) as usize;
                    let y = (e[1] + t * (e[3] - e[1]))
                        .round()
                        .clamp(0.0, (h - 1) as f64) as usize;
                    let (bg, noise) = mesh.at_px(x, y);
                    filtered[y * w + x] > bg + 3.0 * noise
                })
                .count();
            4 * supported >= 3 * (n + 1)
        })
        .collect();
    let prepared = Input {
        raw: &raw,
        s: &filtered,
        mesh: &mesh,
        kernel_sigma: kernel.sigma,
        kernel_gain: kernel.noise_gain(),
        mask,
        peaks: Vec::new(),
        streaks,
        ..*inp
    };
    let (comps, pixels) = components(&prepared, true);
    let (mut census, start) = census(&prepared, &comps, Some(&pixels), cosmic);
    if let Some(out) = kernel_output(&prepared, &comps, &mut census) {
        return Some(out);
    }
    let start = start?;
    if census.regions.0 < original_regions {
        return None;
    }
    let trail = refine(&comps, &geometry_pixels(&prepared), start)?;
    if trail.fitted < 3 {
        return None;
    }
    Some(measure_trails(&prepared, census, trail))
}

pub(crate) fn run(inp: &Input, valid: Option<&[bool]>) -> Output {
    let (comps, _) = components(inp, false);
    let cosmic = crowded_by_hits(&comps, inp.kernel_sigma);
    let (mut census, start) = census(inp, &comps, None, false);
    let retry = cosmic
        || (start.is_none()
            && comps
                .iter()
                .filter(|c| {
                    !c.masked
                        && !c.saturated
                        && c.peak >= BRIGHT
                        && c.major >= OPTICAL_AXIS_RATIO * OPTICAL_AXIS_RATIO * c.minor
                })
                .take(MIN_MODE_VOTES)
                .count()
                == MIN_MODE_VOTES);
    if retry {
        let regions = if start.is_some() { census.regions.0 } else { 0 };
        if let Some(out) = clean_census(inp, valid, cosmic, regions) {
            return out;
        }
    }
    if let Some(out) = kernel_output(inp, &comps, &mut census) {
        return out;
    }
    if let Some(start) = start {
        if let Some(trail) = refine(&comps, &geometry_pixels(inp), start) {
            return measure_trails(inp, census, trail);
        }
        census.mode = FrameMode::Refused;
        census.reason = "the trails do not fit a straight segment".into();
        if let Some(out) = kernel_output(inp, &comps, &mut census) {
            return out;
        }
    }
    Output {
        chains: Vec::new(),
        census,
        trail: None,
        kernel: None,
        sources: Vec::new(),
        threshold: 0.0,
        kernel_sigma: 0.0,
    }
}

fn kernel_output(inp: &Input, comps: &[Component], census: &mut Census) -> Option<Output> {
    let kernel = crate::kernel::run(inp, comps, census)?;
    census.mode = FrameMode::Kernel;
    census.reason = kernel.reason.into();
    Some(Output {
        chains: Vec::new(),
        census: census.clone(),
        trail: None,
        kernel: Some(kernel),
        sources: Vec::new(),
        threshold: 0.0,
        kernel_sigma: 0.0,
    })
}

fn geometry_pixels<'a>(inp: &'a Input<'a>) -> Pixels<'a> {
    Pixels {
        raw: inp.raw,
        orig: inp.orig,
        extra: None,
        w: inp.w,
        h: inp.h,
        mesh: inp.mesh,
        kernel_gain: inp.kernel_gain,
        mask: &inp.mask,
        saturation: inp.saturation,
        gain: inp.gain,
        clipped: inp.clipped.is_some(),
    }
}

fn measure_trails(inp: &Input, census: Census, trail: FrameTrail) -> Output {
    let others: Vec<StreakInput> = inp
        .streaks
        .iter()
        .filter(|s| !star_like(s.1, s.2, s.3, &trail))
        .copied()
        .collect();
    // The streak tracer can join a few roughly collinear trails, across the faint stars and
    // noise between them, into one streak nearly along the trails that star_like does not
    // catch. Subtracted as a streak, it took the trails' light with it, and the trails on it
    // were lost. A satellite's streak is even along its length; such a chain alternates bright
    // trails and faint gaps. A candidate nearly along the trails and uneven is a chain: it is
    // neither subtracted nor reported.
    let (chains, real): (Vec<StreakInput>, Vec<StreakInput>) = others.into_iter().partition(|s| {
        axial_diff(s.1, trail.angle) <= CHAIN_ANGLE
            && unevenness(inp, s.0, s.2, &trail).is_some_and(|u| u > CHAIN_UNEVENNESS)
    });
    let mut out = measure_with(inp, census, trail, &real);
    out.chains = chains.into_iter().map(|s| s.0).collect();
    out
}

/// A streak candidate as the trail stage gets it: endpoints, angle, length, flagged as a hit
/// track, width, flux.
type StreakInput = ([f64; 4], f64, f64, bool, f64, f64);

/// How uneven a streak is along its length: the robust spread of the light in bins half a
/// trail long, across ±2 trail widths, over their median. `None` for fewer than three bins
/// or no light.
fn unevenness(inp: &Input, e: [f64; 4], length: f64, t: &FrameTrail) -> Option<f64> {
    let bin = 0.5 * t.length;
    let bins = (length / bin).floor() as usize;
    if bins < 3 || length <= 0.0 {
        return None;
    }
    let (ux, uy) = ((e[2] - e[0]) / length, (e[3] - e[1]) / length);
    let half = (2.0 * t.width).ceil() as i64;
    let mut light: Vec<f64> = (0..bins)
        .map(|b| {
            let mut sum = 0.0;
            for k in 0..bin.round() as usize {
                let u = b as f64 * bin + k as f64;
                for v in -half..=half {
                    let x = (e[0] + u * ux - v as f64 * uy).round();
                    let y = (e[1] + u * uy + v as f64 * ux).round();
                    if x < 0.0 || y < 0.0 || x >= inp.w as f64 || y >= inp.h as f64 {
                        continue;
                    }
                    let i = y as usize * inp.w + x as usize;
                    if inp.mask[i] {
                        continue;
                    }
                    sum += f64::from(inp.raw[i]) - f64::from(inp.mesh.at(x, y).0);
                }
            }
            sum
        })
        .collect();
    let median = stats::median(&mut light.clone())?;
    if median.is_nan() || median <= 0.0 {
        return None;
    }
    for l in &mut light {
        *l = (*l - median).abs();
    }
    Some(1.4826 * stats::median(&mut light)? / median)
}

fn measure_with(
    inp: &Input,
    mut census: Census,
    trail: FrameTrail,
    others: &[StreakInput],
) -> Output {
    let (w, h) = (inp.w, inp.h);
    let clipped = inp.clipped.is_some();
    let px = geometry_pixels(inp);
    let (len, theta, width) = (trail.length, trail.angle, trail.width);
    if census.regions.0 < 2 {
        let refuse = |census: &mut Census, why: &str| {
            census.mode = FrameMode::Refused;
            census.reason = why.into();
        };
        if trail.fitted < 3 {
            refuse(
                &mut census,
                "too few trails to check that they are straight",
            );
        } else {
            let spread = trail.sigma_angle * (trail.fitted as f64).sqrt() / 1.2533;
            let tol = (ANGLE_AGREEMENT * width / len).max(ANGLE_FLOOR);
            if spread > tol {
                refuse(&mut census, "the trails change across the frame");
            }
        }
        if census.mode == FrameMode::Refused {
            return Output {
                chains: Vec::new(),
                census,
                trail: None,
                kernel: None,
                sources: Vec::new(),
                threshold: 0.0,
                kernel_sigma: 0.0,
            };
        }
    }
    let mut residual: Vec<f32> = inp.raw.to_vec();
    let mut removed = vec![0.0f32; w * h];
    for &(e, _, l, _, sw, flux) in others {
        if l <= 0.0 || flux <= 0.0 {
            continue;
        }
        let (cx, cy) = (0.5 * (e[0] + e[2]), 0.5 * (e[1] + e[3]));
        let (dx, dy) = (e[2] - e[0], e[3] - e[1]);
        let (ux, uy) = (dx / l, dy / l);
        let s_w = (sw * sw + 1.0 / 12.0).sqrt();
        let pad = 0.5 * l + 5.0 * s_w;
        let (x0, x1) = (
            (cx - pad).floor().max(0.0) as usize,
            (cx + pad).ceil().min(w as f64 - 1.0),
        );
        let (y0, y1) = (
            (cy - pad).floor().max(0.0) as usize,
            (cy + pad).ceil().min(h as f64 - 1.0),
        );
        if x1 < 0.0 || y1 < 0.0 {
            continue;
        }
        for yy in y0..=y1 as usize {
            for xx in x0..=x1 as usize {
                let (px_, py_) = (xx as f64 - cx, yy as f64 - cy);
                let (u, v) = (px_ * ux + py_ * uy, -px_ * uy + py_ * ux);
                if v.abs() > 5.0 * s_w || u.abs() > pad {
                    continue;
                }
                let light = (flux * profile(u, v, l, s_w)[0]) as f32;
                residual[yy * w + xx] -= light;
                removed[yy * w + xx] += light;
            }
        }
    }
    let lag = (len + 6.0 * width).ceil() as usize + 1;
    let raw_mesh = inp.clipped.map_or_else(
        || {
            Mesh::build(
                &residual,
                w,
                h,
                inp.mesh_tile.max(24).max(2 * lag),
                lag,
                Some(&inp.mask.iter().map(|m| !m).collect::<Vec<bool>>()),
            )
        },
        |sky| sky.filtered(w, h, 1.0, 0.0),
    );
    let px = Pixels {
        raw: &residual,
        extra: Some(&removed),
        mesh: &raw_mesh,
        kernel_gain: 1.0,
        ..px
    };

    let mut flat = vec![0.0f32; w * h];
    let (mut bg, mut noise) = (vec![0.0f32; w], vec![0.0f32; w]);
    for y in 0..h {
        inp.mesh.row(y, &mut bg, &mut noise);
        let r = y * w..(y + 1) * w;
        for (((f, v), b), m) in flat[r.clone()]
            .iter_mut()
            .zip(&residual[r.clone()])
            .zip(&bg)
            .zip(&inp.mask[r])
        {
            if !*m {
                *f = v - b;
            }
        }
    }
    let boxed = aimed_box(&flat, w, h, theta, len);
    drop(flat);
    let sk = width.clamp(0.6, 4.0);
    let kernel = Kernel::gaussian(sk);
    let (mut tmp, mut f) = (Vec::new(), Vec::new());
    filter::smooth(&boxed, w, h, &kernel, &mut tmp, &mut f);
    drop(boxed);
    drop(tmp);
    let lag = (len + 3.0 * sk).ceil() as usize + 1;
    let tile = inp.mesh_tile.max(24).max(2 * lag);
    let mut fmesh = Mesh::build(&f, w, h, tile, lag, None);
    if let Some(sky) = inp.clipped {
        let across = (sk * sk + 1.0 / 3.0).sqrt();
        let gain = (len / (2.0 * std::f64::consts::PI.sqrt() * across)).sqrt();
        let skew_gain = 2.0 * gain / (3.0_f64.sqrt() * len);
        fmesh.raise_noise(&sky.filtered(w, h, gain, skew_gain));
    }
    let sigma_eff = equivalent_sigma(len, sk);
    let threshold = inp.threshold.unwrap_or_else(|| {
        let z = stats::threshold_for(inp.false_alarms, (w * h) as f64, sigma_eff);
        if clipped {
            z.max(5.0)
        } else {
            z
        }
    });
    let mut seeds: Vec<(usize, usize, f64)> = peaks::find(&f, w, h, &fmesh, threshold as f32, None)
        .into_iter()
        .map(|(x, y)| {
            let (bg, n) = fmesh.at(x as f64, y as f64);
            (x, y, f64::from((f[y * w + x] - bg) / n))
        })
        .collect();
    seeds.sort_by(|a, b| b.2.total_cmp(&a.2).then(a.1.cmp(&b.1)).then(a.0.cmp(&b.0)));
    seeds.truncate(MAX_SEEDS);
    let (sn, cs) = theta.psin_cos();
    let mut kept: Vec<(usize, usize, f64)> = Vec::new();
    for s in seeds {
        let dup = kept.iter().any(|k| {
            let (dx, dy) = (s.0 as f64 - k.0 as f64, s.1 as f64 - k.1 as f64);
            (dx * cs + dy * sn).abs() < 0.5 * len && (-dx * sn + dy * cs).abs() < 1.5 * sk
        });
        if !dup {
            kept.push(s);
        }
    }

    let tests = (w * h) as f64;
    let mut heights: Vec<f64> = kept.iter().map(|s| s.2).collect();
    heights.sort_by(f64::total_cmp);
    let lfdr =
        |z| stats::local_false_discovery_rate(z, threshold, sigma_eff, tests, &heights).min(1.0);

    let start = |sx: f64, sy: f64| -> (f64, f64, f64) {
        let reach = 0.5 * len + 3.0 * width;
        let n = reach.ceil() as i64;
        let (mut m0, mut m1) = (0.0f64, 0.0f64);
        let mut cols: Vec<(f64, f64)> = Vec::with_capacity(2 * n as usize + 1);
        for k in -n..=n {
            let u = k as f64;
            let mut col = 0.0;
            for j in -2..=2 {
                let v = 0.5 * j as f64 * width;
                let (x, y) = (sx + u * cs - v * sn, sy + u * sn + v * cs);
                if let Some(val) = peaks::sample(&residual, w, h, x, y) {
                    let (bg, _) = inp.mesh.at(x, y);
                    col += val - f64::from(bg);
                }
            }
            let c = col.max(0.0);
            m0 += c;
            m1 += c * u;
            cols.push((u, c));
        }
        if m0 > 0.0 {
            let u = (m1 / m0).clamp(-reach, reach) * if clipped { 0.5 } else { 1.0 };
            let peak = cols
                .iter()
                .max_by(|a, b| a.1.total_cmp(&b.1).then(b.0.total_cmp(&a.0)))
                .map_or(0.0, |c| c.0);
            let near: f64 = cols
                .iter()
                .filter(|c| (c.0 - peak).abs() <= 2.0 * width + 1.0)
                .map(|c| c.1)
                .sum();
            (sx + u * cs, sy + u * sn, near / m0)
        } else {
            (sx, sy, 1.0)
        }
    };
    let on_streak = |x: f64, y: f64| {
        others.iter().map(|s| s.0).any(|e| {
            let (dx, dy) = (e[2] - e[0], e[3] - e[1]);
            let l = dx.phypot(dy).max(1e-9);
            let (ux, uy) = (dx / l, dy / l);
            let (px_, py_) = (x - e[0], y - e[1]);
            let along = px_ * ux + py_ * uy;
            let across = (-px_ * uy + py_ * ux).abs();
            across <= 2.0 * width + 1.5 && along >= -0.5 * len && along <= l + 0.5 * len
        })
    };
    let refit = |x0: f64, y0: f64, len_m: f64, others: &[Neighbour]| -> Option<Fit> {
        let (mut x, mut y) = (x0, y0);
        let mut out = None;
        for _ in 0..2 {
            let spec = Spec {
                x0: x,
                y0: y,
                len: len_m,
                theta,
                width,
                free: Free::Held,
                stamp_len: len,
                stamp_width: width,
                max_iter: if clipped { 24 } else { 8 },
            };
            let f = fit(&px, &spec, others)?;
            let moved = (f.x - x).phypot(f.y - y);
            (x, y) = (f.x, f.y);
            out = Some(f);
            if moved <= 1.0 {
                break;
            }
        }
        out
    };
    #[derive(Clone, Copy)]
    struct Measured {
        z: f64,
        conc: f64,
        trail: Fit,
        point: Fit,
    }
    let measure = |x0: f64, y0: f64, z: f64, conc: f64, others: &[Neighbour]| -> Option<Measured> {
        let trail = refit(x0, y0, len, others)?;
        if !(trail.flux > 0.0 && trail.sigma_flux > 0.0) {
            return None;
        }
        let poor = trail.chi2 > POINT_CHECK_CHI2 * trail.dof as f64;
        let point = if conc >= POINT_CHECK_SHARE || poor {
            refit(x0, y0, 0.0, others)?
        } else {
            Fit {
                chi2: f64::INFINITY,
                ..trail
            }
        };
        Some(Measured {
            z,
            conc,
            trail,
            point,
        })
    };
    let as_neighbour = |m: &Measured| -> Neighbour {
        if m.point.chi2 < m.trail.chi2 {
            (m.point.x, m.point.y, m.point.flux.max(0.0), 0.0)
        } else {
            (m.trail.x, m.trail.y, m.trail.flux, len)
        }
    };
    let mut first: Vec<Measured> = Vec::new();
    for &(sx, sy, z) in &kept {
        let (x0, y0, conc) = start(sx as f64, sy as f64);
        if let Some(m) = measure(x0, y0, z, conc, &[]) {
            let best = if m.point.chi2 < m.trail.chi2 {
                &m.point
            } else {
                &m.trail
            };
            if on_streak(best.x, best.y) && best.chi2 > ON_STREAK_CHI2 * best.dof as f64 {
                continue;
            }
            first.push(m);
        }
    }
    // The trail-matched filter spreads one source into a ridge of seeds about two trail
    // lengths long, so the duplicate test above can keep two seeds of one source. Both fits
    // converge on the same position. Keep one, or the neighbour refit below would give each
    // copy the other as a neighbour carrying its whole flux, and both fits would collapse.
    // Two seeds of a blend also converge together, but their single-source fit is poor, and
    // the neighbour refit is what separates the blend: those are kept.
    let best = |m: &Measured| {
        if m.point.chi2 < m.trail.chi2 {
            m.point
        } else {
            m.trail
        }
    };
    let same = SAME_SOURCE_WIDTHS * width;
    let mut unique: Vec<Measured> = Vec::with_capacity(first.len());
    for m in first {
        let b = best(&m);
        let explained = |f: &Fit| f.chi2 <= POINT_CHECK_CHI2 * f.dof as f64;
        let twin = explained(&b)
            .then(|| {
                unique.iter_mut().find(|u| {
                    let o = best(u);
                    explained(&o) && (o.x - b.x).phypot(o.y - b.y) < same
                })
            })
            .flatten();
        match twin {
            Some(u) if m.z > u.z => *u = m,
            Some(_) => {}
            None => unique.push(m),
        }
    }
    let first = unique;
    let all: Vec<Neighbour> = first.iter().map(as_neighbour).collect();
    let fwhm = crate::FWHM_PER_SIGMA * width;
    let reliable: Vec<bool> = first
        .iter()
        .map(|m| {
            let f = if m.point.chi2 < m.trail.chi2 {
                m.point
            } else {
                m.trail
            };
            !clipped
                || (f.converged
                    && f.flux > 5.0 * f.sigma_flux
                    && (f.cov[0] + f.cov[2]).sqrt() < 2.0 * fwhm)
        })
        .collect();
    let reach_u = len + 6.0 * width + 4.0;
    let reach_v = 6.0 * width + 4.0;
    let mut measured: Vec<Measured> = Vec::new();
    for (i, m) in first.iter().enumerate() {
        let (x, y) = (all[i].0, all[i].1);
        let others: Vec<Neighbour> = all
            .iter()
            .enumerate()
            .filter(|(j, o)| {
                let (dx, dy) = (o.0 - x, o.1 - y);
                *j != i
                    && reliable[*j]
                    && (dx * cs + dy * sn).abs() < reach_u
                    && (-dx * sn + dy * cs).abs() < reach_v
            })
            .map(|(_, o)| *o)
            .collect();
        if others.is_empty() {
            measured.push(*m);
        } else {
            measured.push(measure(x, y, m.z, m.conc, &others).unwrap_or(*m));
        }
    }

    let chi2_on = |stamp: (f64, f64), at: (f64, f64), len_m: f64| -> Option<f64> {
        let s_n = (width * width + 1.0 / 12.0).sqrt();
        let half_u = 0.5 * len + 3.0 * width + 3.0;
        let half_v = 3.0 * width + 1.5;
        let reach = half_u.max(half_v).ceil() as i64 + 1;
        let (cx, cy) = (stamp.0.round() as i64, stamp.1.round() as i64);
        let mut ata = [0.0f64; 16];
        let mut atb = [0.0f64; 4];
        let mut rows: Vec<([f64; 4], f64, f64)> = Vec::new();
        for yy in (cy - reach).max(0)..=(cy + reach).min(h as i64 - 1) {
            for xx in (cx - reach).max(0)..=(cx + reach).min(w as i64 - 1) {
                let (dx, dy) = (xx as f64 - stamp.0, yy as f64 - stamp.1);
                if (dx * cs + dy * sn).abs() > half_u || (-dx * sn + dy * cs).abs() > half_v {
                    continue;
                }
                let i = yy as usize * w + xx as usize;
                if inp.mask[i] || inp.saturation.is_some_and(|l| inp.orig[i] >= l) {
                    continue;
                }
                let (ex, ey) = (xx as f64 - at.0, yy as f64 - at.1);
                let (u, v) = (ex * cs + ey * sn, -ex * sn + ey * cs);
                let t = profile(u, v, len_m, s_n)[0];
                let d = f64::from(residual[i]);
                let (bg, _) = px.mesh.at(xx as f64, yy as f64);
                let wt = 1.0 / px.var(xx as usize, yy as usize, d - f64::from(bg));
                let j = [t, 1.0, dx, dy];
                for a in 0..4 {
                    atb[a] += wt * j[a] * d;
                    for b in 0..4 {
                        ata[a * 4 + b] += wt * j[a] * j[b];
                    }
                }
                rows.push((j, d, wt));
            }
        }
        let (sol, _) = cholesky_solve(&ata, &atb, 4)?;
        let chi2 = rows
            .iter()
            .map(|(j, d, wt)| {
                let m: f64 = (0..4).map(|a| sol[a] * j[a]).sum();
                wt * (d - m) * (d - m)
            })
            .sum::<f64>();
        chi2.is_finite().then_some(chi2)
    };

    // A tracked object has the stars' profile across their trails; a hit the first pass could
    // not see (a short track blurred by binning) puts most of its charge in a pixel or two.
    // The share of a centred Gaussian's flux in its brightest pixel bounds a star's.
    let star_share = {
        let e = erf_t(0.5 / (std::f64::consts::SQRT_2 * width.max(0.3)));
        e * e
    };
    let too_sharp = |x: f64, y: f64, flux: f64| -> bool {
        if flux.is_nan() || flux <= 0.0 {
            return false;
        }
        let (cx, cy) = (x.round() as i64, y.round() as i64);
        let mut best: Option<(f64, usize, usize)> = None;
        for yy in (cy - 1).max(0)..=(cy + 1).min(h as i64 - 1) {
            for xx in (cx - 1).max(0)..=(cx + 1).min(w as i64 - 1) {
                let (xu, yu) = (xx as usize, yy as usize);
                let i = yu * w + xu;
                if inp.mask[i] {
                    continue;
                }
                let (bg, _) = inp.mesh.at(xx as f64, yy as f64);
                let v = f64::from(inp.raw[i]) - f64::from(bg);
                if best.is_none_or(|b| v > b.0) {
                    best = Some((v, xu, yu));
                }
            }
        }
        best.is_some_and(|(v, xu, yu)| {
            v - 3.0 * px.var(xu, yu, v).sqrt() > HIT_PEAK_SHARE * star_share * flux
        })
    };
    let mut sources = Vec::new();
    for m in &measured {
        let (t, pt, z) = (m.trail, m.point, m.z);
        if !(t.converged || pt.chi2.is_finite() && pt.converged) {
            continue;
        }
        let d = if pt.chi2.is_finite() {
            match (
                chi2_on((t.x, t.y), (t.x, t.y), len),
                chi2_on((t.x, t.y), (pt.x, pt.y), 0.0),
            ) {
                (Some(ct), Some(cp)) => 0.5 * (cp - ct),
                _ => 0.5 * (pt.chi2 - t.chi2),
            }
        } else {
            f64::INFINITY
        };
        let p_trail = 1.0 / (1.0 + (-d.clamp(-50.0, 50.0)).pexp());
        let kind = if p_trail >= 0.5 {
            Kind::Trail
        } else {
            Kind::Point
        };
        let f = if kind == Kind::Trail { t } else { pt };
        if !f.converged {
            continue;
        }
        let kind = if kind == Kind::Point && too_sharp(f.x, f.y, f.flux) {
            Kind::Hit
        } else {
            kind
        };
        // A trail may run off the frame with its centre outside; a point cannot be measured
        // from outside the image, and one fitted there is the visible end of something else.
        let inside = f.x >= -0.5 && f.y >= -0.5 && f.x <= w as f64 - 0.5 && f.y <= h as f64 - 0.5;
        if kind != Kind::Trail && !inside {
            continue;
        }
        let [cxx, cxy, cyy] = f.cov;
        let along = (cs * cs * cxx + 2.0 * sn * cs * cxy + sn * sn * cyy)
            .max(0.0)
            .sqrt();
        let across = (sn * sn * cxx - 2.0 * sn * cs * cxy + cs * cs * cyy)
            .max(0.0)
            .sqrt();
        if clipped && !(f.x >= 0.0 && f.x <= (w - 1) as f64 && f.y >= 0.0 && f.y <= (h - 1) as f64)
        {
            continue;
        }
        let along_limit = if clipped && kind == Kind::Trail {
            fwhm
        } else if kind == Kind::Trail {
            0.5 * len
        } else {
            2.0 * width
        };
        let fits = along <= along_limit && across <= 2.0 * width;
        if !((fits || kind == Kind::Hit) && f.flux > 0.0) {
            continue;
        }
        let half = if kind == Kind::Trail { 0.5 * len } else { 0.0 };
        let (ex, ey) = (half * cs.abs() + 2.0 * width, half * sn.abs() + 2.0 * width);
        let truncated =
            f.x < ex || f.y < ey || f.x > w as f64 - 1.0 - ex || f.y > h as f64 - 1.0 - ey;
        sources.push(TrailSource {
            x: f.x,
            y: f.y,
            sigma_along: along.max(inp.sigma_floor),
            sigma_across: across.max(inp.sigma_floor),
            flux: f.flux,
            snr: f.flux / f.sigma_flux,
            peak: z,
            p_star: ((1.0 - lfdr(z)) * p_trail).clamp(0.0, 1.0),
            kind,
            truncated,
        });
    }
    Output {
        chains: Vec::new(),
        census,
        trail: Some(trail),
        kernel: None,
        sources,
        threshold,
        kernel_sigma: sk,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn profile_derivatives_match_differences() {
        for &(u, v, l, s) in &[
            (1.3, -0.7, 12.0, 1.4),
            (-6.2, 0.4, 12.0, 1.1),
            (0.3, 1.1, 0.0, 1.3),
        ] {
            let p = profile(u, v, l, s);
            let e = 1e-3;
            let du = (profile(u + e, v, l, s)[0] - profile(u - e, v, l, s)[0]) / (2.0 * e);
            let dv = (profile(u, v + e, l, s)[0] - profile(u, v - e, l, s)[0]) / (2.0 * e);
            let ds = (profile(u, v, l, s + e)[0] - profile(u, v, l, s - e)[0]) / (2.0 * e);
            assert!((p[1] - du).abs() < 1e-4, "du {} {}", p[1], du);
            assert!((p[2] - dv).abs() < 1e-4, "dv {} {}", p[2], dv);
            assert!((p[4] - ds).abs() < 1e-4, "ds {} {}", p[4], ds);
            if l > 0.0 {
                let dl = (profile(u, v, l + e, s)[0] - profile(u, v, l - e, s)[0]) / (2.0 * e);
                assert!((p[3] - dl).abs() < 1e-4, "dl {} {}", p[3], dl);
            }
        }
    }

    #[test]
    fn profile_has_unit_flux() {
        for &l in &[0.0, 5.0, 30.0] {
            let mut sum = 0.0;
            for i in -400..=400 {
                for j in -60..=60 {
                    sum += profile(i as f64 * 0.1, j as f64 * 0.1, l, 1.3)[0] * 0.01;
                }
            }
            assert!((sum - 1.0).abs() < 1e-3, "{l}: {sum}");
        }
    }

    #[test]
    fn equivalent_sigma_of_a_point_is_its_sigma() {
        for &sk in &[0.8, 1.3, 2.0] {
            let e = equivalent_sigma(0.0, sk);
            assert!((e / sk - 1.0).abs() < 0.01, "{sk}: {e}");
        }
        assert!(equivalent_sigma(30.0, 1.3) > 2.0 * 1.3);
    }

    #[test]
    fn aimed_box_on_a_wide_frame() {
        let (w, h) = (20_000usize, 24usize);
        let mut img = vec![0.0f32; w * h];
        img[12 * w + 10_000] = 1.0;
        let out = aimed_box(&img, w, h, std::f64::consts::FRAC_PI_4, 10.0);
        let total: f64 = out.iter().map(|v| f64::from(*v)).sum();
        assert!((total / 10.0 - 1.0).abs() < 0.01, "{total}");
    }

    #[test]
    fn aimed_box_sums_along_the_line() {
        let (w, h) = (64usize, 64usize);
        let mut img = vec![0.0f32; w * h];
        img[32 * w + 30] = 1.0;
        for &deg in &[0.0f64, 17.0, 45.0, 63.0, 90.0, -30.0] {
            let th = deg.to_radians();
            let len = 15.5;
            let out = aimed_box(&img, w, h, th, len);
            let total: f64 = out.iter().map(|v| f64::from(*v)).sum();
            assert!((total / len - 1.0).abs() < 0.01, "{deg}: {total}");
            let (sn, cs) = th.psin_cos();
            let (mut mu, mut sw) = (0.0, 0.0);
            for (i, v) in out.iter().enumerate() {
                let (dx, dy) = ((i % w) as f64 - 30.0, (i / w) as f64 - 32.0);
                mu += f64::from(*v) * (dx * cs + dy * sn);
                sw += f64::from(*v);
            }
            assert!((mu / sw).abs() < 0.05, "{deg}: {}", mu / sw);
        }
    }
}
