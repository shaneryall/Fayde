use fayde_core::Portable;

use crate::mesh::Mesh;
use crate::peaks::{sample, Peak};
use crate::Class;

/// Visible segment of a streak, with endpoints in zero-based image pixels.
#[derive(Debug, Clone, PartialEq)]
pub struct Streak {
    pub id: usize,
    pub x0: f64,
    pub y0: f64,
    pub x1: f64,
    pub y1: f64,
    pub length: f64,
    /// Axis orientation in radians from +x toward +y, modulo π.
    pub angle: f64,
    /// Cross-track Gaussian sigma in pixels, corrected for the detection filter.
    pub width: f64,
    pub flux: f64,
    pub snr: f64,
    /// A crossing star was masked. Flux uses the median ridge height for masked samples. SNR uses
    /// only unmasked samples.
    pub photometry_masked: bool,
    pub coverage: f64,
    pub touches_edge: bool,
    pub cosmic_ray: bool,
    pub(crate) band: f64,
    pub(crate) level: f64,
}

impl Streak {
    pub(crate) fn rescale(&mut self, b: f64) {
        let p = |v: f64| b * v + 0.5 * (b - 1.0);
        (self.x0, self.y0, self.x1, self.y1) = (p(self.x0), p(self.y0), p(self.x1), p(self.y1));
        self.length *= b;
        self.width *= b;
        self.band *= b;
        self.flux *= b * b;
    }

    fn frame(&self, x: f64, y: f64) -> (f64, f64) {
        let (sn, cs) = self.direction();
        let (dx, dy) = (x - self.x0, y - self.y0);
        ((-dx * sn + dy * cs).abs(), dx * cs + dy * sn)
    }

    fn direction(&self) -> (f64, f64) {
        let l = self.length.max(1e-9);
        ((self.y1 - self.y0) / l, (self.x1 - self.x0) / l)
    }

    fn same_line(&self, other: &Streak) -> bool {
        let (a0, t0) = self.frame(other.x0, other.y0);
        let (a1, t1) = self.frame(other.x1, other.y1);
        let band = self.band.max(other.band);
        let (lo, hi) = (t0.min(t1), t0.max(t1));
        let gap = (lo - self.length).max(-hi).max(0.0);
        a0 <= band && a1 <= band && gap <= 0.25 * (self.length + other.length)
    }

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
        let share = (added / other.length.max(1e-9)).min(1.0);
        let traced = self.coverage * self.length + share * other.coverage * other.length;
        self.photometry_masked |= share > 0.0 && other.photometry_masked;
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

pub(crate) struct Traced {
    pub streaks: Vec<Streak>,
    pub members: Vec<Option<usize>>,
}

const TRACE_SIGMA: f64 = 3.0;
const SEED_ELONGATION: f64 = 2.0;
const CROSSING_CONTRAST: f64 = 4.0;
const MIN_LENGTH: f64 = 10.0;
const MIN_LENGTH_SIGMAS: f64 = 8.0;
const MAX_STREAKS: usize = 1000;
const MIN_COVERAGE: f64 = 0.5;
const CLEAN: f64 = 1.25;
const ROUND_SNR: f64 = 12.0;
const MIN_RIDGE_FRACTION: f64 = 0.5;
const MIN_SNR: f64 = 10.0;

#[derive(Clone, Copy)]
struct Point {
    x: f64,
    y: f64,
    value: f64,
    noise: f64,
    clean: bool,
}

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
            let allowed = (8 + pts.len() / 10).min(30);
            if gap > allowed {
                return (pts, false);
            }
        }
    }
}

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

fn peel(pts: &mut Vec<Point>, reach: f64) {
    let dist = |a: &Point, b: &Point| (a.x - b.x).phypot(a.y - b.y);
    for _ in 0..2 {
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
    stars: &[&Peak],
) -> Option<Streak> {
    if pts.len() < 6 {
        return None;
    }
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
    let corr = (2.0 * std::f64::consts::PI.sqrt() * k2.sqrt()).max(1.0);
    let mut masked = vec![false; pts.len()];
    for star in stars
        .iter()
        .filter(|p| p.height > CROSSING_CONTRAST * level)
    {
        let variance = star.major2.min(100.0 * star.minor2).max(psf_var + k2);
        let radius2 = 2.0 * variance * star.snr().max(1.0).pln();
        let across = -(star.x - cx) * sn + (star.y - cy) * cs;
        let t = (star.x - cx) * cs + (star.y - cy) * sn;
        if across * across > radius2 || t < lo - radius2.sqrt() || t > hi + radius2.sqrt() {
            continue;
        }
        for (p, mask) in pts.iter().zip(&mut masked) {
            let (dx, dy) = (p.x - star.x, p.y - star.y);
            *mask |= dx * dx + dy * dy <= radius2;
        }
    }
    let clean = || {
        pts.iter()
            .zip(&masked)
            .filter_map(|(p, &m)| (!m).then_some(p))
    };
    let n_clean = clean().count();
    if n_clean == 0 {
        return None;
    }
    let snr = clean().map(|p| p.value / p.noise).sum::<f64>() / (n_clean as f64 * corr).sqrt();
    if !(length >= min_len && coverage >= MIN_COVERAGE && snr >= MIN_SNR) {
        return None;
    }
    let flux = (clean().map(|p| p.value).sum::<f64>() + (pts.len() - n_clean) as f64 * level)
        * (2.0 * std::f64::consts::PI * total2).sqrt();
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
        photometry_masked: n_clean != pts.len(),
        coverage,
        touches_edge: edge.0 || edge.1,
        cosmic_ray: false,
        band: {
            let mut noise: Vec<f64> = pts.iter().map(|p| p.noise).collect();
            let noise = crate::stats::median(&mut noise).unwrap_or(1.0);
            let ratio = (level / noise.max(1e-12)).max(7.4);
            (sigma_c * (2.0 * ratio.pln()).sqrt()).max(1.5) + 1.0
        },
        level,
    })
}

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
    stars: &[&Peak],
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
    let same = |a: &Point, b: &Point| a.x == b.x && a.y == b.y;
    let edge = (
        edge_b && same(&before.0, &pts[0]),
        edge_f && same(&before.1, &pts[pts.len() - 1]),
    );
    fit(s, w, h, mesh, &pts, edge, k2, psf_var, stars)
}

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
    let ridge = |i: usize| peaks[i].fitted && peaks[i].elongation() >= SEED_ELONGATION * stretch;
    let stellar = |i: usize| !ridge(i) && matches!(classes[i], Class::Star | Class::Saturated);
    let stars: Vec<_> = peaks
        .iter()
        .enumerate()
        .filter_map(|(i, p)| stellar(i).then_some(p))
        .collect();
    let mut seeds: Vec<usize> = (0..peaks.len())
        .filter(|&i| !skip[i] && ridge(i) && matches!(classes[i], Class::Extended | Class::Star))
        .collect();
    seeds.sort_by(|&a, &b| peaks[b].snr().total_cmp(&peaks[a].snr()).then(a.cmp(&b)));

    let member = |t: &Streak, i: usize| {
        let p = &peaks[i];
        let star_on_top = stellar(i) && p.height > CROSSING_CONTRAST * t.level;
        t.contains(p.x, p.y) && !star_on_top
    };

    let mut streaks: Vec<Streak> = Vec::new();
    for &i in &seeds {
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
        let Some(mut t) = trace_from(s, w, h, mesh, seed, (cs, sn), k2, psf_var, &stars) else {
            continue;
        };
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
            if let Some(again) = trace_from(s, w, h, mesh, start, (dx, dy), k2, psf_var, &stars) {
                if again.length > t.length {
                    t = again;
                }
            }
        }
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
        if let Some(first) = streaks.iter_mut().find(|f| f.same_line(&t)) {
            first.absorb(&t);
            continue;
        }
        t.id = streaks.len();
        streaks.push(t);
    }
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
