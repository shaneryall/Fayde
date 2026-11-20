mod rotation;

use fayde_core::Portable;
use rotation::Rotation;

use crate::mesh::Mesh;
use crate::trail::{axial_diff, Component, Input};
use crate::{peaks, stats, Census, Class, Detected, FrameMode};

const MIN_STAMPS: usize = 5;
const MAX_STAMPS: usize = 64;
const MAX_RADIUS: usize = 64;
const MARGIN: usize = 5;
const MIN_BEND: f64 = 0.5;
const OVERRIDE_STAMPS: usize = 15;
const MIN_STACKED: usize = 10;
const CORE_FRACTION: f32 = 0.02;
const OVERRIDE_BEND: f64 = 1.0;
const MAX_SEEDS: usize = 5000;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FrameKernel {
    pub width: usize,
    pub height: usize,
    pub stacked: usize,
    pub path_length: f64,
    pub bend: f64,
    pub regions: usize,
}

pub(crate) struct Output {
    pub reason: &'static str,
    pub info: FrameKernel,
    pub sources: Vec<Detected>,
    pub threshold: f64,
    pub equivalent_sigma: f64,
    pub background: f64,
    pub noise: f64,
    pub candidates: usize,
    pub psf_sigma: f64,
    pub keep_streak: Vec<bool>,
}

struct Stamp {
    pixels: Vec<f32>,
    origin: (f64, f64),
    noise: f64,
}

#[derive(Clone, Copy)]
struct Tap {
    x: isize,
    y: isize,
    value: f32,
}

struct Template {
    pixels: Vec<f32>,
    side: usize,
    taps: Vec<Tap>,
    power: f64,
    peak: f64,
    support: Vec<(isize, isize)>,
    gradient: [f64; 3],
    theta: f64,
    major: f64,
    minor: f64,
    psf_sigma: f64,
    info: FrameKernel,
}

struct Image<'a> {
    pixels: Vec<f32>,
    input: &'a Input<'a>,
    mesh: Mesh,
    noise: f64,
    background: f64,
}

impl<'a> Image<'a> {
    fn new(input: &'a Input<'a>) -> Self {
        let (w, h) = (input.w, input.h);
        let valid: Vec<bool> = input.mask.iter().map(|&m| !m).collect();
        let mesh = Mesh::build(input.raw, w, h, input.mesh_tile, 3, Some(&valid));
        let (background, noise) = mesh.typical();
        let mut pixels = vec![0.0; w * h];
        let (mut bg, mut sig) = (vec![0.0; w], vec![0.0; w]);
        for (y, row) in pixels.chunks_exact_mut(w).enumerate() {
            mesh.row(y, &mut bg, &mut sig);
            let start = y * w;
            for (x, ((out, &b), &ok)) in row
                .iter_mut()
                .zip(&bg)
                .zip(&valid[start..start + w])
                .enumerate()
            {
                if ok {
                    *out = input.raw[start + x] - b;
                }
            }
        }
        Self {
            pixels,
            input,
            mesh,
            noise: f64::from(noise),
            background: f64::from(background),
        }
    }
}

fn training_stars<'a>(
    comps: &'a [Component],
    theta: f64,
    length: f64,
    input: &Input,
) -> Vec<&'a Component> {
    let mut stars: Vec<_> = comps
        .iter()
        .filter(|c| {
            !c.edge
                && !c.masked
                && !c.saturated
                && c.peak >= 15.0
                && c.elongated()
                && axial_diff(c.theta, theta) <= 0.15
                && (c.length() / length - 1.0).abs() <= 0.22
        })
        .collect();
    stars.sort_by(|a, b| b.peak.total_cmp(&a.peak));
    stars
        .into_iter()
        .take(4 * MAX_STAMPS)
        .filter(|c| isolated(c, comps, input) && stamp_radius(c) <= MAX_RADIUS)
        .take(MAX_STAMPS)
        .collect()
}

fn isolated(c: &Component, comps: &[Component], input: &Input) -> bool {
    let [x0, y0, x1, y1] = c.bounds;
    if x0 < MARGIN || y0 < MARGIN || x1 + MARGIN >= input.w || y1 + MARGIN >= input.h {
        return false;
    }
    !comps.iter().any(|other| {
        if other.peak < 8.0 || std::ptr::eq(c, other) {
            return false;
        }
        let [a, b, d, e] = other.bounds;
        a <= x1 + MARGIN && d + MARGIN >= x0 && b <= y1 + MARGIN && e + MARGIN >= y0
    })
}

fn stamp_radius(c: &Component) -> usize {
    let [x0, y0, x1, y1] = c.bounds.map(|n| n as f64);
    (c.x - x0).max(x1 - c.x).max(c.y - y0).max(y1 - c.y).ceil() as usize + MARGIN
}

fn square(radius: usize) -> impl Iterator<Item = (isize, isize)> {
    let r = radius as isize;
    (-r..=r).flat_map(move |y| (-r..=r).map(move |x| (x, y)))
}

fn cubic(p: [f64; 4], t: f64) -> (f64, f64) {
    let a = 0.5 * (p[2] - p[0]);
    let b = p[0] - 2.5 * p[1] + 2.0 * p[2] - 0.5 * p[3];
    let c = 0.5 * (p[3] - p[0]) + 1.5 * (p[1] - p[2]);
    (
        p[1] + t * (a + t * (b + t * c)),
        a + t * (2.0 * b + t * 3.0 * c),
    )
}

fn cubic_value(p: [f64; 4], t: f64) -> f64 {
    let a = 0.5 * (p[2] - p[0]);
    let b = p[0] - 2.5 * p[1] + 2.0 * p[2] - 0.5 * p[3];
    let c = 0.5 * (p[3] - p[0]) + 1.5 * (p[1] - p[2]);
    p[1] + t * (a + t * (b + t * c))
}

fn pixel(grid: &[f32], side: usize, x: isize, y: isize) -> f64 {
    if (0..side as isize).contains(&x) && (0..side as isize).contains(&y) {
        f64::from(grid[y as usize * side + x as usize])
    } else {
        0.0
    }
}

fn sample(grid: &[f32], side: usize, x: f64, y: f64) -> (f64, f64, f64) {
    let (x0, y0) = (x.floor() as isize, y.floor() as isize);
    let (tx, ty) = (x - x0 as f64, y - y0 as f64);
    let mut values = [0.0; 4];
    let mut slopes = [0.0; 4];
    for (j, (v, d)) in values.iter_mut().zip(&mut slopes).enumerate() {
        let row =
            std::array::from_fn(|i| pixel(grid, side, x0 + i as isize - 1, y0 + j as isize - 1));
        (*v, *d) = cubic(row, tx);
    }
    let (v, dy) = cubic(values, ty);
    (v, cubic(slopes, ty).0, dy)
}

fn sample_value(grid: &[f32], side: usize, x: f64, y: f64) -> f64 {
    let (x0, y0) = (x.floor() as isize, y.floor() as isize);
    let (tx, ty) = (x - x0 as f64, y - y0 as f64);
    let values = std::array::from_fn(|j| {
        let row =
            std::array::from_fn(|i| pixel(grid, side, x0 + i as isize - 1, y0 + j as isize - 1));
        cubic_value(row, tx)
    });
    cubic_value(values, ty)
}

fn cutout(image: &Image, c: &Component, radius: usize) -> Option<Stamp> {
    let (cx, cy) = (c.x.round() as isize, c.y.round() as isize);
    let [x0, y0, x1, y1] = c.bounds;
    let (w, h) = (image.input.w as isize, image.input.h as isize);
    let side = 2 * radius + 1;
    let mut pixels = vec![0.0; side * side];
    let (mut sum, mut mx, mut my) = (0.0, 0.0, 0.0);
    for ((x, y), out) in square(radius).zip(&mut pixels) {
        let (ix, iy) = (cx + x, cy + y);
        if !(0..w).contains(&ix) || !(0..h).contains(&iy) {
            continue;
        }
        let i = (iy * w + ix) as usize;
        if ix >= x0 as isize - 1
            && ix <= x1 as isize + 1
            && iy >= y0 as isize - 1
            && iy <= y1 as isize + 1
            && (image.input.mask[i]
                || image
                    .input
                    .saturation
                    .is_some_and(|s| image.input.raw[i] >= s))
        {
            return None;
        }
        *out = image.pixels[i];
        let (b, n) = image.input.mesh.at_px(ix as usize, iy as usize);
        if image.input.s[i] > b + 3.0 * n {
            let v = f64::from(*out);
            sum += v;
            mx += v * x as f64;
            my += v * y as f64;
        }
    }
    if !sum.is_finite() || sum <= 0.0 || (mx / sum).abs() > 3.0 || (my / sum).abs() > 3.0 {
        return None;
    }
    pixels
        .iter_mut()
        .for_each(|v| *v = (f64::from(*v) / sum) as f32);
    let (_, noise) = image.mesh.at(c.x, c.y);
    Some(Stamp {
        pixels,
        origin: (radius as f64 + mx / sum, radius as f64 + my / sum),
        noise: f64::from(noise) / sum,
    })
}

fn stack(stamps: &[Stamp], side: usize) -> Vec<f32> {
    let radius = side / 2;
    let mut values = Vec::with_capacity(stamps.len());
    let mut scratch = Vec::with_capacity(stamps.len());
    square(radius)
        .map(|(x, y)| {
            values.clear();
            values.extend(stamps.iter().map(|s| {
                sample_value(
                    &s.pixels,
                    side,
                    s.origin.0 + x as f64,
                    s.origin.1 + y as f64,
                )
            }));
            scratch.clear();
            scratch.extend(&values);
            let median = stats::median(&mut scratch).expect("training stars");
            let (mut light, mut weight) = (0.0, 0.0);
            for (&v, s) in values.iter().zip(stamps) {
                let scale = s.noise + 0.05 * median.abs();
                if (v - median).abs() <= 4.0 * scale {
                    let w = 1.0 / (scale * scale);
                    light += w * v;
                    weight += w;
                }
            }
            if weight > 0.0 && light / weight > 2.0 / weight.sqrt() {
                (light / weight) as f32
            } else {
                0.0
            }
        })
        .collect()
}

fn quadratic(f: impl Fn(isize, isize) -> f64) -> (f64, f64) {
    let z = f(0, 0);
    let (a, b, c, d) = (f(-1, 0), f(1, 0), f(0, -1), f(0, 1));
    let (xx, yy) = (a + b - 2.0 * z, c + d - 2.0 * z);
    let xy = 0.25 * (f(1, 1) + f(-1, -1) - f(-1, 1) - f(1, -1));
    let det = xx * yy - xy * xy;
    if xx >= 0.0 || yy >= 0.0 || det <= 0.0 {
        return (0.0, 0.0);
    }
    let (gx, gy) = (0.5 * (b - a), 0.5 * (d - c));
    (
        ((xy * gy - yy * gx) / det).clamp(-0.75, 0.75),
        ((xy * gx - xx * gy) / det).clamp(-0.75, 0.75),
    )
}

fn train(image: &Image, stars: &[&Component]) -> Option<Template> {
    let radius = stars.iter().map(|&c| stamp_radius(c)).max()?;
    let side = 2 * radius + 1;
    let mut stamps: Vec<_> = stars
        .iter()
        .filter_map(|c| cutout(image, c, radius))
        .collect();
    if stamps.len() < MIN_STAMPS {
        return None;
    }
    let mut mean = stack(&stamps, side);
    for _ in 0..2 {
        let taps: Vec<_> = square(radius)
            .zip(&mean)
            .filter(|(_, v)| **v > 0.0)
            .collect();
        for stamp in &mut stamps {
            let (dx, dy) = quadratic(|dx, dy| {
                taps.iter()
                    .map(|&((x, y), v)| {
                        f64::from(*v)
                            * sample_value(
                                &stamp.pixels,
                                side,
                                stamp.origin.0 + (x + dx) as f64,
                                stamp.origin.1 + (y + dy) as f64,
                            )
                    })
                    .sum()
            });
            stamp.origin.0 += dx;
            stamp.origin.1 += dy;
        }
        mean = stack(&stamps, side);
    }
    let peak = mean.iter().fold(0.0f32, |m, &v| m.max(v));
    let (mut sum, mut mx, mut my) = (0.0, 0.0, 0.0);
    for ((x, y), &v) in square(radius).zip(&mean) {
        if v < CORE_FRACTION * peak {
            continue;
        }
        let v = f64::from(v);
        sum += v;
        mx += v * x as f64;
        my += v * y as f64;
    }
    if !sum.is_finite() || sum <= 0.0 {
        return None;
    }
    let pixels = square(radius)
        .map(|(x, y)| {
            sample_value(
                &mean,
                side,
                (x + radius as isize) as f64 + mx / sum,
                (y + radius as isize) as f64 + my / sum,
            )
            .max(0.0) as f32
        })
        .collect();
    Template::new(pixels, side, stamps.len(), None)
}

impl Template {
    fn new(
        mut pixels: Vec<f32>,
        side: usize,
        stacked: usize,
        geometry: Option<PathGeometry>,
    ) -> Option<Self> {
        let sum: f64 = pixels.iter().map(|&v| f64::from(v)).sum();
        if !sum.is_finite() || sum <= 0.0 {
            return None;
        }
        let radius = side / 2;
        let mut taps: Vec<_> = square(radius)
            .zip(&pixels)
            .filter(|(_, v)| **v > 0.0)
            .map(|((x, y), &value)| Tap { x, y, value })
            .collect();
        taps.sort_by(|a, b| b.value.total_cmp(&a.value));
        let mut light = 0.0;
        taps.retain(|t| {
            let keep = light < 0.99 * sum;
            light += f64::from(t.value);
            keep
        });
        let peak = f64::from(taps.first()?.value) / sum;
        let kept: f64 = taps.iter().map(|t| f64::from(t.value)).sum();
        taps.sort_by_key(|t| (t.y, t.x));
        let [mut x0, mut y0, mut x1, mut y1] = [0; 4];
        let [mut xx, mut xy, mut yy] = [0.0; 3];
        let mut power = 0.0;
        for t in &mut taps {
            let v = f64::from(t.value) / kept;
            t.value = v as f32;
            let (x, y) = (t.x as f64, t.y as f64);
            xx += x * x * v;
            xy += x * y * v;
            yy += y * y * v;
            power += v * v;
            x0 = x0.min(t.x);
            y0 = y0.min(t.y);
            x1 = x1.max(t.x);
            y1 = y1.max(t.y);
        }
        for v in &mut pixels {
            *v = (f64::from(*v) / sum) as f32;
        }
        let theta = 0.5 * (2.0 * xy).patan2(xx - yy);
        let scale = 0.5 * (xx - yy).phypot(2.0 * xy);
        let major = (0.5 * (xx + yy) + scale).sqrt();
        let minor = (0.5 * (xx + yy) - scale).max(0.0).sqrt();
        let geometry = geometry.unwrap_or_else(|| path(&taps, theta, radius));
        let mut gradient = [0.0; 3];
        for (x, y) in square(radius) {
            let (x, y) = (x + radius as isize, y + radius as isize);
            let gx = 0.5 * (pixel(&pixels, side, x + 1, y) - pixel(&pixels, side, x - 1, y));
            let gy = 0.5 * (pixel(&pixels, side, x, y + 1) - pixel(&pixels, side, x, y - 1));
            gradient[0] += gx * gx;
            gradient[1] += gx * gy;
            gradient[2] += gy * gy;
        }
        let r = radius as isize;
        let mut support = vec![false; side * side];
        for t in &taps {
            for sy in (t.y - 1).max(-r)..=(t.y + 1).min(r) {
                for sx in (t.x - 1).max(-r)..=(t.x + 1).min(r) {
                    support[(sy + r) as usize * side + (sx + r) as usize] = true;
                }
            }
        }
        let support = square(radius)
            .zip(support)
            .filter_map(|(offset, used)| used.then_some(offset))
            .collect();
        Some(Self {
            pixels,
            side,
            taps,
            power,
            peak,
            support,
            gradient,
            theta,
            major,
            minor,
            psf_sigma: geometry.width,
            info: FrameKernel {
                width: (x1 - x0 + 1) as usize,
                height: (y1 - y0 + 1) as usize,
                stacked,
                path_length: geometry.length,
                bend: geometry.bend,
                regions: 1,
            },
        })
    }

    fn equivalent_sigma(&self) -> f64 {
        let [xx, xy, yy] = self.gradient;
        (0.5 * self.power / (xx * yy - xy * xy).sqrt()).sqrt()
    }

    fn sample(&self, x: f64, y: f64) -> (f64, f64, f64) {
        let r = (self.side / 2) as f64;
        sample(&self.pixels, self.side, x + r, y + r)
    }

    fn sample_value(&self, x: f64, y: f64) -> f64 {
        let r = (self.side / 2) as f64;
        sample_value(&self.pixels, self.side, x + r, y + r)
    }
}

struct PathGeometry {
    length: f64,
    bend: f64,
    width: f64,
}

fn path(taps: &[Tap], theta: f64, radius: usize) -> PathGeometry {
    let (sin, cos) = theta.psin_cos();
    let mut slices = vec![[0.0; 6]; 2 * radius + 3];
    for t in taps {
        let (x, y) = (t.x as f64, t.y as f64);
        let (u, v) = (cos * x + sin * y, cos * y - sin * x);
        let i = (0.5 * u).round() as isize + radius as isize + 1;
        let a = &mut slices[i as usize];
        let light = f64::from(t.value);
        a[0] += light;
        a[1] += light * u;
        a[2] += light * v;
        a[3] += light * u * u;
        a[4] += light * u * v;
        a[5] += light * v * v;
    }
    let max_light = slices.iter().map(|s| s[0]).fold(0.0, f64::max);
    let points: Vec<_> = slices
        .into_iter()
        .filter(|s| s[0] > 0.25 * max_light)
        .map(|[w, u, v, uu, uv, vv]| (w, u / w, v / w, uu / w, uv / w, vv / w))
        .collect();
    let mut width2 = Vec::new();
    let (mut sum, mut bend, mut length) = (0.0, 0.0, 2.0);
    for (i, &(w, u, v, uu, uv, vv)) in points.iter().enumerate() {
        sum += w;
        bend += w * v * v;
        if i > 0 {
            length += (u - points[i - 1].1).phypot(v - points[i - 1].2);
        }
        if i > 0 && i + 1 < points.len() {
            let (lo, hi) = (points[i - 1], points[i + 1]);
            let slope = (hi.2 - lo.2) / (hi.1 - lo.1);
            let var = (vv - v * v - 2.0 * slope * (uv - u * v) + slope * slope * (uu - u * u))
                / (1.0 + slope * slope);
            width2.push(var);
        }
    }
    PathGeometry {
        length,
        bend: (bend / sum).sqrt(),
        width: (stats::median(&mut width2).unwrap_or(1.0) - 1.0 / 12.0)
            .max(0.05)
            .sqrt(),
    }
}

#[derive(Clone, Copy)]
struct Region {
    x0: usize,
    y0: usize,
    x1: usize,
    y1: usize,
}

fn correlate(image: &Image, template: &Template, region: Region, filtered: &mut [f32]) {
    let (w, h) = (image.input.w, image.input.h);
    for &Tap { x, y, value } in &template.taps {
        let x0 = region.x0.max((-x).max(0) as usize);
        let x1 = region.x1.min(w.saturating_sub(x.max(0) as usize));
        let y0 = region.y0.max((-y).max(0) as usize);
        let y1 = region.y1.min(h.saturating_sub(y.max(0) as usize));
        if x0 >= x1 {
            continue;
        }
        let n = x1 - x0;
        for row in y0..y1 {
            let a = (row as isize + y) as usize * w + (x0 as isize + x) as usize;
            let b = row * w + x0;
            for (out, &v) in filtered[b..b + n].iter_mut().zip(&image.pixels[a..a + n]) {
                *out += v * value;
            }
        }
    }
}

#[derive(Clone, Copy)]
struct Fit {
    x: f64,
    y: f64,
    flux: f64,
    snr: f64,
    sigma_px: f64,
    sigma_along: f64,
    sigma_across: f64,
    truncated: bool,
    saturated: bool,
    shape: [f64; 3],
}

struct FitPixel {
    i: usize,
    x: f64,
    y: f64,
}

fn fit(template: &Template, image: &Image, start: (f64, f64), bg: f64, noise: f64) -> Option<Fit> {
    let input = image.input;
    let (ix, iy) = (start.0.round() as isize, start.1.round() as isize);
    let (mut dx, mut dy) = (start.0 - ix as f64, start.1 - iy as f64);
    let mut samples = Vec::with_capacity(template.support.len());
    let (mut truncated, mut saturated) = (false, false);
    for &(x, y) in &template.support {
        let (px, py) = (ix + x, iy + y);
        if !(0..input.w as isize).contains(&px) || !(0..input.h as isize).contains(&py) {
            truncated = true;
            continue;
        }
        let i = py as usize * input.w + px as usize;
        if input.mask[i] {
            truncated = true;
            continue;
        }
        if input.saturation.is_some_and(|s| input.raw[i] >= s) {
            saturated = true;
            continue;
        }
        samples.push(FitPixel {
            i,
            x: x as f64,
            y: y as f64,
        });
    }
    let mut model = Vec::with_capacity(samples.len());
    for iter in 0..=5 {
        let (mut dot, mut norm) = (0.0, 0.0);
        model.clear();
        for s in &samples {
            let (k, gx, gy) = template.sample(s.x - dx, s.y - dy);
            let k = k.max(0.0);
            let v = f64::from(image.pixels[s.i]) - bg;
            norm += k * k;
            dot += k * v;
            model.push((v, k, -gx, -gy));
        }
        if norm < 0.30 * template.power || dot <= 0.0 {
            return None;
        }
        let flux = dot / norm;
        let snr = dot / norm.sqrt() / noise;
        let [mut xx, mut xy, mut yy, mut bx, mut by] = [0.0; 5];
        let mut chi2 = 0.0;
        let mut n = 0;
        for &(v, k, gx, gy) in &model {
            let expected = flux * k;
            let residual = v - expected;
            let photon = input.gain.map_or(0.0, |gain| expected / gain);
            let variance = noise * noise + photon + (0.055 * expected) * (0.055 * expected);
            let weight = 1.0 / variance;
            let (gx, gy) = (flux * gx, flux * gy);
            xx += weight * gx * gx;
            xy += weight * gx * gy;
            yy += weight * gy * gy;
            bx += weight * gx * residual;
            by += weight * gy * residual;
            if k > 0.15 * template.peak {
                chi2 += weight * residual * residual;
                n += 1;
            }
        }
        let det = xx * yy - xy * xy;
        if !det.is_finite() || det <= 0.0 {
            return None;
        }
        let sx = ((yy * bx - xy * by) / det).clamp(-0.75, 0.75);
        let sy = ((xx * by - xy * bx) / det).clamp(-0.75, 0.75);
        if iter == 5 || sx.abs().max(sy.abs()) < 0.025 {
            if chi2 > 3.5 * n.max(1) as f64 {
                return None;
            }
            let mean = 0.5 * (xx + yy);
            let diff = 0.5 * (xx - yy).phypot(2.0 * xy);
            let floor = input.sigma_floor * input.sigma_floor;
            let (sn, cs) = template.theta.psin_cos();
            let cross = 2.0 * sn * cs * xy;
            let along = (cs * cs * yy - cross + sn * sn * xx) / det;
            let across = (sn * sn * yy + cross + cs * cs * xx) / det;
            return Some(Fit {
                x: ix as f64 + dx,
                y: iy as f64 + dy,
                flux,
                snr,
                sigma_px: (floor + (mean + diff) / det).sqrt(),
                sigma_along: (floor + along).sqrt(),
                sigma_across: (floor + across).sqrt(),
                truncated,
                saturated,
                shape: [template.major, template.minor, template.theta],
            });
        }
        dx = (dx + sx).clamp(-2.5, 2.5);
        dy = (dy + sy).clamp(-2.5, 2.5);
    }
    None
}

fn subtract(image: &mut Image, template: &Template, fit: Fit) {
    let (w, h) = (image.input.w as isize, image.input.h as isize);
    let (cx, cy) = (fit.x.round() as isize, fit.y.round() as isize);
    for (x, y) in square(template.side / 2) {
        let (ix, iy) = (cx + x, cy + y);
        if (0..w).contains(&ix) && (0..h).contains(&iy) {
            let i = (iy * w + ix) as usize;
            if !image.input.mask[i] {
                let v = template
                    .sample_value(ix as f64 - fit.x, iy as f64 - fit.y)
                    .max(0.0);
                let mut model = (v * fit.flux) as f32;
                if image
                    .input
                    .saturation
                    .is_some_and(|s| image.input.raw[i] >= s)
                {
                    model = model.min(image.pixels[i].max(0.0));
                }
                image.pixels[i] -= model;
            }
        }
    }
}

fn path_streaks(image: &Image, template: &Template, fit: Fit, keep: &mut [bool]) {
    for (streak, keep) in image.input.streaks.iter().zip(keep) {
        let ([x0, y0, x1, y1], _, len, cosmic, _, _) = *streak;
        let radius = template.side as f64 * 0.5;
        if *keep
            && !cosmic
            && len <= template.side as f64
            && (0.5 * (x0 + x1) - fit.x).abs() < radius
            && (0.5 * (y0 + y1) - fit.y).abs() < radius
        {
            *keep = !(0..=4).all(|i| {
                let t = 0.25 * i as f64;
                let x = x0 + t * (x1 - x0) - fit.x;
                let y = y0 + t * (y1 - y0) - fit.y;
                template.sample_value(x, y) > 0.06 * template.peak
            });
        }
    }
}

fn measure(
    mut image: Image,
    regions: &[(Region, Template)],
    rotation: Option<&Rotation>,
    reason: &'static str,
) -> Option<Output> {
    let input = image.input;
    let (w, h) = (input.w, input.h);
    let median = |value: fn(&Template) -> f64| {
        stats::median(&mut regions.iter().map(|(_, t)| value(t)).collect::<Vec<_>>())
            .expect("at least one region")
    };
    let sigma = median(Template::equivalent_sigma);
    let tests = (w * h) as f64;
    let threshold = input
        .threshold
        .unwrap_or_else(|| stats::threshold_for(input.false_alarms, tests, sigma));
    let mut info = regions.first()?.1.info;
    if regions.iter().any(|(_, t)| t.info.stacked < MIN_STACKED) {
        return None;
    }
    info.width = regions.iter().map(|(_, t)| t.info.width).max()?;
    info.height = regions.iter().map(|(_, t)| t.info.height).max()?;
    info.path_length = median(|t| t.info.path_length);
    info.bend = median(|t| t.info.bend);
    info.regions = regions.len();
    let lag = info.width.max(info.height) + 2;
    let mut detections: Vec<Fit> = Vec::new();
    let mut candidates = 0;
    let mut keep_streak = vec![true; input.streaks.len()];
    for _ in 0..2 {
        let mut filtered = vec![0.0; w * h];
        for (region, template) in regions {
            correlate(&image, template, *region, &mut filtered);
        }
        let mesh = Mesh::build(&filtered, w, h, input.mesh_tile.max(2 * lag), lag, None);
        let mut seeds: Vec<_> = peaks::find(&filtered, w, h, &mesh, threshold as f32, None)
            .into_iter()
            .map(|(x, y)| {
                let (bg, _) = mesh.at_px(x, y);
                (x, y, filtered[y * w + x] - bg)
            })
            .collect();
        seeds.sort_by(|a, b| b.2.total_cmp(&a.2).then(a.1.cmp(&b.1)).then(a.0.cmp(&b.0)));
        candidates += seeds.len();
        for (x, y, _) in seeds.into_iter().take(MAX_SEEDS) {
            if detections
                .iter()
                .any(|f| (f.x - x as f64).abs() < 2.0 && (f.y - y as f64).abs() < 2.0)
            {
                continue;
            }
            let (bx, by) = quadratic(|dx, dy| {
                f64::from(filtered[(y as isize + dy) as usize * w + (x as isize + dx) as usize])
            });
            let start = (x as f64 + bx, y as f64 + by);
            let local;
            let template = if let Some(rotation) = rotation {
                let Some(t) = rotation.template(start.0, start.1) else {
                    continue;
                };
                local = t;
                &local
            } else {
                &regions[0].1
            };
            let (bg, filtered_noise) = mesh.at_px(x, y);
            let noise = (f64::from(filtered_noise) / template.power.sqrt())
                .max(f64::from(image.mesh.at_px(x, y).1));
            let Some(f) = fit(template, &image, start, f64::from(bg), noise) else {
                continue;
            };
            if f.snr >= threshold
                && !detections
                    .iter()
                    .any(|d| (d.x - f.x).phypot(d.y - f.y) < 2.0)
            {
                subtract(&mut image, template, f);
                path_streaks(&image, template, f, &mut keep_streak);
                detections.push(f);
            }
        }
    }
    if detections.len() < MIN_STAMPS {
        return None;
    }
    let mut heights: Vec<_> = detections.iter().map(|f| f.snr).collect();
    heights.sort_by(f64::total_cmp);
    let sources = detections
        .into_iter()
        .map(|f| {
            let lfdr = stats::local_false_discovery_rate(f.snr, threshold, sigma, tests, &heights);
            Detected {
                x: f.x,
                y: f.y,
                sigma_px: f.sigma_px,
                sigma_along: f.sigma_along,
                sigma_across: f.sigma_across,
                flux: f.flux,
                snr: f.snr,
                p_star: 1.0 - lfdr.clamp(0.0, 1.0),
                sigma_major: f.shape[0],
                sigma_minor: f.shape[1],
                theta: f.shape[2],
                class: if f.saturated {
                    Class::Saturated
                } else {
                    Class::Star
                },
                streak: None,
                truncated: f.truncated,
            }
        })
        .collect();
    Some(Output {
        reason,
        info,
        sources,
        threshold,
        equivalent_sigma: sigma,
        background: image.background,
        noise: image.noise,
        candidates,
        psf_sigma: median(|t| t.psf_sigma),
        keep_streak,
    })
}

pub(crate) fn run(inp: &Input, components: &[Component], census: &Census) -> Option<Output> {
    let mut built = None;
    if census.mode != FrameMode::Trails
        && census.elongated_share >= 0.5
        && census.mode_share < 0.55
        && census.votes >= 12
    {
        let image = Image::new(inp);
        match Rotation::learn(components, &image)
            .and_then(|r| r.regions(inp.w, inp.h).map(|g| (r, g)))
        {
            Some((rotation, regions)) => {
                if let Some(output) = measure(
                    image,
                    &regions,
                    Some(&rotation),
                    "star paths vary with field rotation",
                ) {
                    return Some(output);
                }
            }
            None => built = Some(image),
        }
    }
    if census.mode == FrameMode::Points || census.mode_share < 0.5 {
        return None;
    }
    let (theta, length) = census.shapes?;
    let stars = training_stars(components, theta, length, inp);
    let needed = if census.mode == FrameMode::Trails {
        OVERRIDE_STAMPS
    } else {
        MIN_STAMPS
    };
    if stars.len() < needed.max(MIN_STAMPS) {
        return None;
    }
    let image = built.unwrap_or_else(|| Image::new(inp));
    let template = train(&image, &stars)?;
    if template.info.bend < MIN_BEND {
        return None;
    }
    if census.mode == FrameMode::Trails
        && (template.info.stacked < OVERRIDE_STAMPS || template.info.bend < OVERRIDE_BEND)
    {
        return None;
    }
    let region = Region {
        x0: 0,
        x1: inp.w,
        y0: 0,
        y1: inp.h,
    };
    measure(
        image,
        &[(region, template)],
        None,
        "shared curved star path",
    )
}
