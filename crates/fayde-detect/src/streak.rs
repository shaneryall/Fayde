//! Streaks: straight linear features traced along their ridge.
//!
//! A peak on a streak is a ridge: the filtered image curves down across the streak and hardly
//! at all along it, which the 3 × 3 fit reports as a large elongation and a direction (the
//! eigenvector of the weaker curvature; Steger 1998, IEEE TPAMI 20, 113, uses the same Hessian
//! for line points). From each such peak the ridge is followed in both directions, one pixel at
//! a time, along a straight line that is refitted by least squares as points are added, and
//! stops where the ridge stays below the trace threshold. Holding the trace to a line keeps it
//! on course through a bright star beside the streak or another streak crossing it. The final
//! line is a total-least-squares fit, and each end is put where the ridge falls to half its
//! level. A trace is kept as a streak only if most of the peaks on it are themselves ridges,
//! which a row of stars is not. The tracing and these tests are Fayde's own.
//!
//! Only streaks bright enough to raise a peak above the detection threshold are found. Fainter
//! ones do not disturb the star list, which is what the plate solver needs; finding them takes
//! an integration along the line (Nir, Zackay & Ofek 2018, AJ 156, 229) that is not done here.

use fayde_core::Portable;

use crate::mesh::Mesh;
use crate::peaks::{sample, Peak};
use crate::Class;

/// A traced streak. Positions are zero-based pixels of the input image.
#[derive(Debug, Clone, PartialEq)]
pub struct Streak {
    /// Position in [`crate::Detection::streaks`].
    pub id: usize,
    /// One end.
    pub x0: f64,
    /// One end.
    pub y0: f64,
    /// The other end.
    pub x1: f64,
    /// The other end.
    pub y1: f64,
    /// Length, pixels.
    pub length: f64,
    /// Direction, radians from +x towards +y, in `(−π/2, π/2]`.
    pub angle: f64,
    /// Gaussian σ of the cross-section, pixels, filter removed.
    pub width: f64,
    /// Flux above the background, input units × pixels.
    pub flux: f64,
    /// Signal-to-noise ratio of the whole streak.
    pub snr: f64,
    /// Fraction of the length where the ridge is above the trace threshold (below 1 for a
    /// tumbling object or one crossing gaps).
    pub coverage: f64,
    /// An end is at the frame edge: the streak continues outside.
    pub touches_edge: bool,
    /// The cross-section is narrower than the stars': the track of a cosmic ray through the
    /// detector, not something imaged through the optics.
    pub cosmic_ray: bool,
    /// Half-width of the band whose detections belong to the streak, pixels.
    pub(crate) band: f64,
    /// Median ridge height above the background, filtered-image units.
    pub(crate) level: f64,
}

impl Streak {
    /// Convert from pixels of an image binned by `b` to pixels of the input image.
    pub(crate) fn rescale(&mut self, b: f64) {
        let p = |v: f64| b * v + 0.5 * (b - 1.0);
        (self.x0, self.y0, self.x1, self.y1) = (p(self.x0), p(self.y0), p(self.x1), p(self.y1));
        self.length *= b;
        self.width *= b;
        self.band *= b;
        self.flux *= b * b;
    }

    /// Distance of `(x, y)` from the streak's axis and position along it from the first end.
    fn frame(&self, x: f64, y: f64) -> (f64, f64) {
        let (sn, cs) = self.direction();
        let (dx, dy) = (x - self.x0, y - self.y0);
        ((-dx * sn + dy * cs).abs(), dx * cs + dy * sn)
    }

    fn direction(&self) -> (f64, f64) {
        let l = self.length.max(1e-9);
        ((self.y1 - self.y0) / l, (self.x1 - self.x0) / l)
    }

    /// Whether `other` lies along this streak's line: both its ends within the band, and no
    /// further from this streak along the line than a gap a tumbling object could leave.
    fn same_line(&self, other: &Streak) -> bool {
        let (a0, t0) = self.frame(other.x0, other.y0);
        let (a1, t1) = self.frame(other.x1, other.y1);
        let band = self.band.max(other.band);
        let (lo, hi) = (t0.min(t1), t0.max(t1));
        let gap = (lo - self.length).max(-hi).max(0.0);
        a0 <= band && a1 <= band && gap <= 0.25 * (self.length + other.length)
    }

    /// Extend this streak to cover `other`, which lies on the same line.
    fn absorb(&mut self, other: &Streak) {
        let (sn, cs) = self.direction();
        let t = [
            0.0,
            self.length,
            self.frame(other.x0, other.y0).1,
            self.frame(other.x1, other.y1).1,
        ];
        let lo = t.iter().copied().fold(f64::INFINITY, f64::min);
        let hi = t.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let added = ((hi - lo) - self.length).max(0.0);
        // The part of `other` outside this streak brings its share of flux and of the
        // traced length.
        let share = (added / other.length.max(1e-9)).min(1.0);
        let traced = self.coverage * self.length + share * other.coverage * other.length;
        self.flux += share * other.flux;
        self.snr = (self.snr * self.snr + share * other.snr * other.snr).sqrt();
        let (x0, y0) = (self.x0, self.y0);
        (self.x0, self.y0) = (x0 + lo * cs, y0 + lo * sn);
        (self.x1, self.y1) = (x0 + hi * cs, y0 + hi * sn);
        self.length = hi - lo;
        self.coverage = (traced / self.length.max(1e-9)).min(1.0);
        self.touches_edge |= other.touches_edge;
    }

    fn contains(&self, x: f64, y: f64) -> bool {
        let (across, along) = self.frame(x, y);
        across <= self.band && along >= -2.0 * self.band && along <= self.length + 2.0 * self.band
    }
}

/// Streaks and, for each peak, the streak it lies on.
pub(crate) struct Traced {
    pub streaks: Vec<Streak>,
    pub members: Vec<Option<usize>>,
}

/// Ridge points accepted while tracing must exceed this many σ of the filtered noise. Lower
/// than the peak threshold: a point is tested only where a ridge is already being followed.
const TRACE_SIGMA: f64 = 3.0;
/// A peak this many times more elongated (major over minor axis) than the frame's stars
/// starts a trace.
const SEED_ELONGATION: f64 = 2.0;
/// Fewest pixels of length for a streak, and the same in units of the cross-section σ.
const MIN_LENGTH: f64 = 10.0;
const MIN_LENGTH_SIGMAS: f64 = 8.0;
/// Most streaks traced in a frame.
const MAX_STREAKS: usize = 1000;
/// Least fraction of the length above the trace threshold.
const MIN_COVERAGE: f64 = 0.5;
/// A traced point this close to the predicted line, pixels, refines the line.
const CLEAN: f64 = 1.25;
/// A peak that is not a ridge counts against a trace only above this signal-to-noise ratio,
/// where the error of its fitted axis ratio (about 2 / SNR per axis) leaves no doubt.
const ROUND_SNR: f64 = 12.0;
/// Least fraction of the peaks on a trace that must be ridges.
const MIN_RIDGE_FRACTION: f64 = 0.5;
/// Least signal-to-noise ratio of the whole streak.
const MIN_SNR: f64 = 10.0;

#[derive(Clone, Copy)]
struct Point {
    x: f64,
    y: f64,
    value: f64,
    noise: f64,
    /// On the predicted line: used for the line fit.
    clean: bool,
}

/// Follow the ridge from `start` in direction `dir`. Returns the accepted points in order and
/// whether the trace ran off the frame while still on the ridge.
fn walk(
    s: &[f32],
    w: usize,
    h: usize,
    mesh: &Mesh,
    start: (f64, f64),
    dir: (f64, f64),
) -> (Vec<Point>, bool) {
    let mut pts: Vec<Point> = Vec::new();
    let (ux, uy) = dir;
    let (nx, ny) = (-uy, ux);
    // The line, in the frame of the seed: lateral offset d = a + b·t at distance t along
    // `dir`. Sums for its least-squares fit over the clean points.
    let (mut a, mut b) = (0.0f64, 0.0f64);
    let (mut n, mut st, mut stt, mut sd, mut std) = (0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64);
    let mut t = 0.0f64;
    let mut gap = 0usize;
    loop {
        t += 1.0;
        let d0 = a + b * t;
        let (x, y) = (start.0 + t * ux + d0 * nx, start.1 + t * uy + d0 * ny);
        if x < 2.0 || y < 2.0 || x > w as f64 - 3.0 || y > h as f64 - 3.0 {
            return (pts, gap <= 2);
        }
        // The brightest of five samples across the track, refined with a parabola.
        let mut v = [f64::NEG_INFINITY; 5];
        for (i, val) in v.iter_mut().enumerate() {
            let o = i as f64 - 2.0;
            if let Some(q) = sample(s, w, h, x + o * nx, y + o * ny) {
                *val = q;
            }
        }
        let mut best = 2;
        for i in 0..5 {
            if v[i] > v[best] {
                best = i;
            }
        }
        let mut off = best as f64 - 2.0;
        let mut top = v[best];
        if (1..4).contains(&best) && v[best - 1].is_finite() && v[best + 1].is_finite() {
            let curve = v[best - 1] - 2.0 * v[best] + v[best + 1];
            if curve < 0.0 {
                let d = (0.5 * (v[best - 1] - v[best + 1]) / curve).clamp(-0.5, 0.5);
                off += d;
                top -= 0.25 * (v[best - 1] - v[best + 1]) * d;
            }
        }
        let (bg, noise) = mesh.at(x, y);
        let (bg, noise) = (f64::from(bg), f64::from(noise));
        if top - bg > TRACE_SIGMA * noise {
            // Near the seed its 3 × 3 direction is all there is, and it can be several degrees
            // off: every point found counts until the line has a lever arm.
            let clean = off.abs() <= CLEAN || t <= 24.0;
            pts.push(Point {
                x: x + off * nx,
                y: y + off * ny,
                value: top - bg,
                noise,
                clean,
            });
            gap = 0;
            if clean {
                let d = d0 + off;
                n += 1.0;
                st += t;
                stt += t * t;
                sd += d;
                std += t * d;
                let det = n * stt - st * st;
                if n >= 4.0 && det > 0.0 {
                    b = (n * std - st * sd) / det;
                    a = (sd - b * st) / n;
                } else {
                    a = sd / n;
                }
            }
        } else {
            gap += 1;
            // Gaps (a tumbling object, a dip of a faint streak under the noise) are bridged in
            // proportion to what has been traced, so a short feature cannot leap to a
            // neighbour.
            let allowed = (8 + pts.len() / 10).min(30);
            if gap > allowed {
                return (pts, false);
            }
        }
    }
}

/// Drop accepted points at the far end of a walk that stand alone: a noise sample above the
/// trace threshold inside the gap allowance would otherwise extend the streak by up to that
/// allowance. The end is where at least 5 of the last 8 pixels were accepted.
fn trim(pts: &mut Vec<Point>) {
    while let Some(last) = pts.last().copied() {
        let near = pts
            .iter()
            .rev()
            .take_while(|p| (p.x - last.x).phypot(p.y - last.y) <= 8.0)
            .count();
        if near >= 5 || pts.len() < 5 {
            break;
        }
        pts.pop();
    }
}

/// Remove from each end of an ordered trace what is not the streak:
///
/// - a group of points beyond a gap that is shorter than the gap (or than `reach` pixels): a
///   star or a noise clump the gap allowance reached, where a continuing streak would run on
///   for at least as long as it was interrupted;
/// - a star touching the end: points more than twice as bright as the ridge within half of
///   `reach` of the end (a star's own radius; a brighter streak crossing this one further in
///   is left alone), with everything outside them, up to the dip before the streak.
fn peel(pts: &mut Vec<Point>, reach: f64) {
    let dist = |a: &Point, b: &Point| (a.x - b.x).phypot(a.y - b.y);
    for _ in 0..2 {
        // Groups beyond a gap, from this end inwards.
        loop {
            let n = pts.len();
            let Some(i) = (1..n / 2).find(|&i| dist(&pts[i - 1], &pts[i]) > 2.5) else {
                break;
            };
            let extent = dist(&pts[0], &pts[i - 1]);
            let gap = dist(&pts[i - 1], &pts[i]);
            if extent >= gap.max(reach) {
                break;
            }
            pts.drain(..i);
        }
        if pts.len() >= 16 {
            let mut values: Vec<f64> = pts.iter().map(|p| p.value).collect();
            let level = crate::stats::median(&mut values).unwrap_or(0.0);
            // The innermost bright point near the end, then the faintest point within `reach`
            // further in: the cut.
            let bright = (0..pts.len() / 2)
                .take_while(|&i| dist(&pts[0], &pts[i]) <= 0.5 * reach)
                .filter(|&i| pts[i].value > 2.0 * level)
                .last();
            if let Some(core) = bright {
                let mut cut = core + 1;
                for i in core + 1..pts.len() / 2 {
                    if dist(&pts[core], &pts[i]) > reach {
                        break;
                    }
                    if pts[i].value < pts[cut].value {
                        cut = i;
                    }
                }
                pts.drain(..cut);
            }
        }
        pts.reverse();
    }
}

/// Fit a streak to traced points. `None` if they do not make one.
#[allow(clippy::too_many_arguments)]
fn fit(
    s: &[f32],
    w: usize,
    h: usize,
    mesh: &Mesh,
    pts: &[Point],
    edge: (bool, bool),
    k2: f64,
    psf_var: f64,
) -> Option<Streak> {
    if pts.len() < 6 {
        return None;
    }
    // Total least squares over the points on the line: through their centroid, along the
    // principal axis.
    let on_line = || pts.iter().filter(|p| p.clean);
    let nc = on_line().count();
    if nc < 4 {
        return None;
    }
    let (cx, cy) = (
        on_line().map(|p| p.x).sum::<f64>() / nc as f64,
        on_line().map(|p| p.y).sum::<f64>() / nc as f64,
    );
    let (mut sxx, mut sxy, mut syy) = (0.0, 0.0, 0.0);
    for p in on_line() {
        let (dx, dy) = (p.x - cx, p.y - cy);
        sxx += dx * dx;
        sxy += dx * dy;
        syy += dy * dy;
    }
    let mut theta = 0.5 * (2.0 * sxy).patan2(sxx - syy);
    // Point from the first traced point to the last, so that `edge` keeps its meaning.
    let (first, last) = (pts[0], pts[pts.len() - 1]);
    {
        let (sn, cs) = theta.psin_cos();
        if (last.x - first.x) * cs + (last.y - first.y) * sn < 0.0 {
            theta += std::f64::consts::PI;
        }
    }
    let (sn, cs) = theta.psin_cos();
    let along = |p: &Point| (p.x - cx) * cs + (p.y - cy) * sn;
    let (lo, hi) = (along(&first), along(&last));

    // Cross-section variance: three samples across the ridge, in the logarithm (exact for a
    // Gaussian profile), at up to 32 points; the median.
    let mut widths: Vec<f64> = Vec::new();
    let step = (pts.len() / 32).max(1);
    for p in pts.iter().step_by(step) {
        let (bg, _) = mesh.at(p.x, p.y);
        let v = |o: f64| sample(s, w, h, p.x - o * sn, p.y + o * cs).map(|q| q - f64::from(bg));
        if let (Some(a), Some(c), Some(b)) = (v(-1.0), v(0.0), v(1.0)) {
            if a > 0.0 && b > 0.0 && c > 0.0 {
                let d = a.pln() + b.pln() - 2.0 * c.pln();
                if d < 0.0 {
                    widths.push(-1.0 / d);
                }
            }
        }
    }
    let total2 = crate::stats::median(&mut widths).unwrap_or(psf_var + k2);
    let sigma_c = total2.sqrt();
    let mut values: Vec<f64> = pts.iter().map(|p| p.value).collect();
    let level = crate::stats::median(&mut values).unwrap_or(0.0);

    // An end of a uniform streak blurred by the PSF and the filter is where the ridge falls to
    // half its level; the trace, which stops at the trace threshold, runs past that on a
    // bright streak and short of it on a faint one. From a point inside, step outwards half a
    // pixel at a time to the crossing.
    let ridge = |t: f64| -> Option<f64> {
        let (x, y) = (cx + t * cs, cy + t * sn);
        let (bg, _) = mesh.at(x, y);
        sample(s, w, h, x, y).map(|q| q - f64::from(bg))
    };
    let reach = 3.0 * sigma_c + 2.0;
    let end = |t_end: f64, out: f64| -> f64 {
        let mut t = t_end - out * reach;
        let mut prev = match ridge(t) {
            Some(v) if v >= 0.5 * level => v,
            _ => return t_end,
        };
        while (t - t_end) * out < reach {
            t += out * 0.5;
            match ridge(t) {
                Some(v) if v < 0.5 * level => {
                    return t - out * 0.5 * (0.5 * level - v) / (prev - v);
                }
                Some(v) => prev = v,
                None => return t_end,
            }
        }
        t_end
    };
    let lo = if edge.0 { lo } else { end(lo, -1.0) };
    let hi = if edge.1 { hi } else { end(hi, 1.0) };
    let length = hi - lo;

    let min_len = MIN_LENGTH.max(MIN_LENGTH_SIGMAS * sigma_c);
    let nf = pts.len() as f64;
    let coverage = (nf / (length + 1.0)).min(1.0);
    // Along the track the filtered noise is correlated over 2√π·σk pixels, so the points hold
    // that many times fewer independent ones.
    let corr = (2.0 * std::f64::consts::PI.sqrt() * k2.sqrt()).max(1.0);
    let snr = pts.iter().map(|p| p.value / p.noise).sum::<f64>() / (nf * corr).sqrt();
    if !(length >= min_len && coverage >= MIN_COVERAGE && snr >= MIN_SNR) {
        return None;
    }
    let flux =
        pts.iter().map(|p| p.value).sum::<f64>() * (2.0 * std::f64::consts::PI * total2).sqrt();
    let mut angle = theta;
    while angle > std::f64::consts::FRAC_PI_2 {
        angle -= std::f64::consts::PI;
    }
    while angle <= -std::f64::consts::FRAC_PI_2 {
        angle += std::f64::consts::PI;
    }
    Some(Streak {
        id: 0,
        x0: cx + lo * cs,
        y0: cy + lo * sn,
        x1: cx + hi * cs,
        y1: cy + hi * sn,
        length,
        angle,
        width: (total2 - k2).max(0.0).sqrt(),
        flux,
        snr,
        coverage,
        touches_edge: edge.0 || edge.1,
        // Set by `trace_all`, which knows the stars' narrow-axis width.
        cosmic_ray: false,
        // Out to where the cross-section falls to the noise: peaks there are the streak's
        // flanks, not stars.
        band: {
            let mut noise: Vec<f64> = pts.iter().map(|p| p.noise).collect();
            let noise = crate::stats::median(&mut noise).unwrap_or(1.0);
            let ratio = (level / noise.max(1e-12)).max(7.4);
            (sigma_c * (2.0 * ratio.pln()).sqrt()).max(1.5) + 1.0
        },
        level,
    })
}

/// Walk both ways from `seed` along `dir`, clean the ends and fit a streak.
#[allow(clippy::too_many_arguments)]
fn trace_from(
    s: &[f32],
    w: usize,
    h: usize,
    mesh: &Mesh,
    seed: Point,
    dir: (f64, f64),
    k2: f64,
    psf_var: f64,
) -> Option<Streak> {
    let (mut fwd, edge_f) = walk(s, w, h, mesh, (seed.x, seed.y), dir);
    let (mut back, edge_b) = walk(s, w, h, mesh, (seed.x, seed.y), (-dir.0, -dir.1));
    trim(&mut fwd);
    trim(&mut back);
    let mut pts: Vec<Point> = back.into_iter().rev().collect();
    pts.push(seed);
    pts.extend(fwd);
    let before = (pts[0], pts[pts.len() - 1]);
    peel(&mut pts, 6.0 * (psf_var + k2).sqrt() + 2.0);
    // An end that was cut back no longer reaches the frame edge.
    let same = |a: &Point, b: &Point| a.x == b.x && a.y == b.y;
    let edge = (
        edge_b && same(&before.0, &pts[0]),
        edge_f && same(&before.1, &pts[pts.len() - 1]),
    );
    fit(s, w, h, mesh, &pts, edge, k2, psf_var)
}

/// Trace a streak from every ridge-like peak that is not already on one.
#[allow(clippy::too_many_arguments)]
pub(crate) fn trace_all(
    s: &[f32],
    w: usize,
    h: usize,
    mesh: &Mesh,
    peaks: &[Peak],
    classes: &[Class],
    skip: &[bool],
    k2: f64,
    psf_var: f64,
    (psf_minor, stretch): (f64, f64),
) -> Traced {
    // A ridge is that much more elongated than the frame's stars, which a drift during the
    // exposure may have stretched (`stretch`, 1 for round stars).
    let ridge = |i: usize| peaks[i].fitted && peaks[i].elongation() >= SEED_ELONGATION * stretch;
    let mut seeds: Vec<usize> = (0..peaks.len())
        .filter(|&i| !skip[i] && ridge(i) && matches!(classes[i], Class::Extended | Class::Star))
        .collect();
    seeds.sort_by(|&a, &b| peaks[b].snr().total_cmp(&peaks[a].snr()).then(a.cmp(&b)));

    // Whether peak `i` belongs to streak `t`: on it, and not a round peak far brighter than the
    // ridge (a star the streak crosses).
    let member = |t: &Streak, i: usize| {
        let p = &peaks[i];
        let star_on_top = !ridge(i)
            && matches!(classes[i], Class::Star | Class::Saturated)
            && p.height > 4.0 * t.level;
        t.contains(p.x, p.y) && !star_on_top
    };

    let mut streaks: Vec<Streak> = Vec::new();
    for &i in &seeds {
        // Far more than any frame holds; the cost of the loops over streaks stays bounded.
        if streaks.len() >= MAX_STREAKS {
            break;
        }
        let p = &peaks[i];
        if streaks.iter().any(|t| t.contains(p.x, p.y)) {
            continue;
        }
        let (sn, cs) = p.theta.psin_cos();
        let seed = Point {
            x: p.x,
            y: p.y,
            value: p.height,
            noise: p.noise,
            clean: true,
        };
        let Some(mut t) = trace_from(s, w, h, mesh, seed, (cs, sn), k2, psf_var) else {
            continue;
        };
        // Once more from the middle of the fitted line: its direction, from the whole trace,
        // is far better than the seed's 3 × 3 fit, whose error can end a walk early on a faint
        // streak. The longer of the two traces is kept.
        let (dy, dx) = t.direction();
        let mid = (0.5 * (t.x0 + t.x1), 0.5 * (t.y0 + t.y1));
        let (bg, noise) = mesh.at(mid.0, mid.1);
        if let Some(value) = sample(s, w, h, mid.0, mid.1) {
            let start = Point {
                x: mid.0,
                y: mid.1,
                value: value - f64::from(bg),
                noise: f64::from(noise),
                clean: true,
            };
            if let Some(again) = trace_from(s, w, h, mesh, start, (dx, dy), k2, psf_var) {
                if again.length > t.length {
                    t = again;
                }
            }
        }
        // A streak's peaks are ridges; a row of stars, or two blended ones, has round peaks
        // with at most a ridge-like saddle between them. Only peaks bright enough for their
        // shape to be told count as round: on a faint streak the noise rounds many.
        let (mut ridges, mut others) = (0usize, 0usize);
        for j in (0..peaks.len()).filter(|&j| !skip[j] && member(&t, j)) {
            if ridge(j) {
                ridges += 1;
            } else if classes[j] != Class::Sharp && peaks[j].snr() >= ROUND_SNR {
                others += 1;
            }
        }
        if (ridges as f64) < MIN_RIDGE_FRACTION * (ridges + others) as f64 {
            continue;
        }
        // The same line found again (from a peak beside the first trace's band, or past a gap
        // it did not bridge): extend the first instead of reporting two.
        if let Some(first) = streaks.iter_mut().find(|f| f.same_line(&t)) {
            first.absorb(&t);
            continue;
        }
        t.id = streaks.len();
        streaks.push(t);
    }
    // Narrower than the stars across their narrow axis: the criterion of the "sharp" class.
    for t in &mut streaks {
        t.cosmic_ray = t.width * t.width < 0.4 * psf_minor;
    }

    let members = (0..peaks.len())
        .map(|i| {
            (!skip[i])
                .then(|| streaks.iter().position(|t| member(t, i)))
                .flatten()
        })
        .collect();
    Traced { streaks, members }
}
