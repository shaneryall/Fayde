use std::ops::Range;

use crate::stats::{normal_density, normal_quantile_upper, normal_tail};

pub(crate) struct ClippedBackground {
    floor: f32,
    latent: Mesh,
    padding: Vec<(Range<usize>, Range<usize>)>,
}

fn censored_normal(buf: &mut [f32], floor: f32) -> Option<(f32, f32)> {
    let n = buf.len();
    buf.sort_unstable_by(f32::total_cmp);
    let c = buf.partition_point(|&v| v <= floor);
    if c < 8 || n - c < 16 {
        return None;
    }
    let upper = buf[c + (n - c) / 2];
    let k = buf.partition_point(|&v| v <= upper);
    if k == n {
        return None;
    }
    let mid = |a, b| (f64::from(a) + f64::from(b)) / 2.0;
    let x0 = mid(floor, buf[c]);
    let x1 = mid(upper, buf[k]);
    let quantile = |rank| {
        if rank <= n / 2 {
            -normal_quantile_upper(rank as f64 / n as f64)
        } else {
            normal_quantile_upper((n - rank) as f64 / n as f64)
        }
    };
    let (z0, z1) = (quantile(c), quantile(k));
    let sigma = (x1 - x0) / (z1 - z0);
    let mu = x0 - sigma * z0;
    let (mu, sigma) = (mu as f32, sigma as f32);
    (mu.is_finite() && sigma.is_finite() && sigma > 0.0).then_some((mu, sigma))
}

fn floor_is_interspersed(
    raw: &[f32],
    w: usize,
    xs: Range<usize>,
    ys: Range<usize>,
    floor: f32,
    fraction: f64,
) -> bool {
    let (mut pairs, mut transitions) = (0, 0);
    let mut pair = |a: f32, b: f32| {
        if a.is_finite() && b.is_finite() {
            pairs += 1;
            transitions += usize::from((a == floor) != (b == floor));
        }
    };
    for y in ys.clone() {
        let row = &raw[y * w + xs.start..y * w + xs.end];
        for p in row.windows(2) {
            pair(p[0], p[1]);
        }
        if y > ys.start {
            let prev = &raw[(y - 1) * w + xs.start..(y - 1) * w + xs.end];
            for (&a, &b) in row.iter().zip(prev) {
                pair(a, b);
            }
        }
    }
    transitions as f64 > 0.5 * fraction * (1.0 - fraction) * pairs as f64
}

impl ClippedBackground {
    pub fn estimate(raw: &[f32], w: usize, h: usize, tile: usize) -> Option<Self> {
        let (mut floor, mut count, mut n) = (f32::INFINITY, 0usize, 0usize);
        for &v in raw.iter().filter(|v| v.is_finite()) {
            n += 1;
            if v < floor {
                floor = v;
                count = 1;
            } else if v == floor {
                count += 1;
            }
        }
        if count < n.div_ceil(4) || n - count < 16 {
            return None;
        }
        let (nx, ny, tw, th) = Mesh::grid(w, h, tile.max(64));
        let mut bg = vec![f32::NAN; nx * ny];
        let mut noise = bg.clone();
        let mut buf = Vec::new();
        let mut padding = Vec::new();
        let (mut sky_tiles, mut fits) = (0, 0);
        for ty in 0..ny {
            let (y0, y1) = Mesh::bounds(ty, th, h);
            for tx in 0..nx {
                let (x0, x1) = Mesh::bounds(tx, tw, w);
                buf.clear();
                for y in y0..y1 {
                    buf.extend(raw[y * w + x0..y * w + x1].iter().filter(|v| v.is_finite()));
                }
                let n = buf.len();
                if n == 0 {
                    continue;
                }
                let count = buf.iter().filter(|&&v| v == floor).count();
                if count >= n - n / 100
                    || (count >= 8
                        && !floor_is_interspersed(
                            raw,
                            w,
                            x0..x1,
                            y0..y1,
                            floor,
                            count as f64 / n as f64,
                        ))
                {
                    padding.push((x0..x1, y0..y1));
                    continue;
                }
                if n < 24 {
                    continue;
                }
                sky_tiles += 1;
                if count < n.div_ceil(4) {
                    continue;
                }
                if let Some((mu, sigma)) = censored_normal(&mut buf, floor) {
                    bg[ty * nx + tx] = mu;
                    noise[ty * nx + tx] = sigma;
                    fits += 1;
                }
            }
        }
        if fits <= sky_tiles / 2 {
            return None;
        }
        let latent = Mesh::from_tiles(w, h, nx, ny, bg, noise);
        Some(Self {
            floor,
            latent,
            padding,
        })
    }

    pub fn valid_pixels(&self, raw: &[f32], w: usize) -> Option<Vec<bool>> {
        if self.padding.is_empty() {
            return None;
        }
        let mut valid: Vec<bool> = raw.iter().map(|v| v.is_finite()).collect();
        for (xs, ys) in &self.padding {
            for y in ys.clone() {
                let row = y * w + xs.start..y * w + xs.end;
                for (ok, &v) in valid[row.clone()].iter_mut().zip(&raw[row]) {
                    *ok &= v != self.floor;
                }
            }
        }
        Some(valid)
    }

    pub fn filtered(&self, w: usize, h: usize, noise_gain: f64, skew_gain: f64) -> Mesh {
        let mut bg = Vec::with_capacity(self.latent.bg.len());
        let mut noise = Vec::with_capacity(self.latent.noise.len());
        for (&mu, &sigma) in self.latent.bg.iter().zip(&self.latent.noise) {
            let (mu, sigma, floor) = (f64::from(mu), f64::from(sigma), f64::from(self.floor));
            let a = (floor - mu) / sigma;
            let (q, phi) = (normal_tail(a), normal_density(a));
            let m1 = phi - a * q;
            let m2 = q - a * m1;
            let m3 = 2.0 * m1 - a * m2;
            let sd = (m2 - m1 * m1).sqrt();
            let skew = ((m3 - 3.0 * m2 * m1 + 2.0 * m1 * m1 * m1) / (sd * sd * sd)).max(0.0);
            bg.push((floor + sigma * m1) as f32);
            noise.push((sigma * sd * noise_gain * (1.0 + 0.8 * skew_gain * skew)) as f32);
        }
        Mesh::from_tiles(w, h, self.latent.nx, self.latent.ny, bg, noise)
    }
}

#[derive(Debug, Clone)]
pub(crate) struct Mesh {
    nx: usize,
    ny: usize,
    tw: f64,
    th: f64,
    bg: Vec<f32>,
    noise: Vec<f32>,
    cols: Vec<(u32, u32, f32)>,
    rows: Vec<(u32, u32, f32)>,
}

fn clipped_median(buf: &mut [f32]) -> Option<f32> {
    let n = buf.len();
    if n < 8 {
        return None;
    }
    let cmp = |a: &f32, b: &f32| a.total_cmp(b);
    let m = n / 2;
    let med = *buf.select_nth_unstable_by(m, cmp).1;
    let q = *buf[..m]
        .select_nth_unstable_by(((0.1587 * n as f64) as usize).min(m - 1), cmp)
        .1;
    let hi = med + 3.0 * (med - q);
    let kept = buf.iter().filter(|v| **v <= hi).count();
    let m2 = (kept / 2).min(m);
    if m2 == m {
        return Some(med);
    }
    Some(*buf[..m].select_nth_unstable_by(m2, cmp).1)
}

fn pair_noise(dif: &mut [f32]) -> Option<f32> {
    let n = dif.len();
    if n < 8 {
        return None;
    }
    let cmp = |a: &f32, b: &f32| a.total_cmp(b);
    let med = *dif.select_nth_unstable_by(n / 2, cmp).1;
    for v in dif.iter_mut() {
        *v = (*v - med).abs();
    }
    let mut top = n / 2;
    let mut mad = *dif.select_nth_unstable_by(top, cmp).1;
    let mut factor = 1.4826f32;
    for _ in 0..3 {
        let hi = 3.0 * factor * mad;
        let kept = dif.iter().filter(|v| **v <= hi).count();
        let m = kept / 2;
        if m >= top || m == 0 {
            break;
        }
        mad = *dif[..top].select_nth_unstable_by(m, cmp).1;
        top = m;
        factor = 1.4872;
    }
    Some(factor * mad / std::f32::consts::SQRT_2)
}

fn despike(v: &[f32], nx: usize, ny: usize, excess: impl Fn(f32) -> f32) -> Vec<f32> {
    let mut out = v.to_vec();
    if nx < 3 || ny < 3 {
        return out;
    }
    for y in 1..ny - 1 {
        for x in 1..nx - 1 {
            let at = |dx: isize, dy: isize| {
                v[(y as isize + dy) as usize * nx + (x as isize + dx) as usize]
            };
            let predicted = (0.5 * (at(-1, 0) + at(1, 0)))
                .max(0.5 * (at(0, -1) + at(0, 1)))
                .max(0.5 * (at(-1, -1) + at(1, 1)))
                .max(0.5 * (at(-1, 1) + at(1, -1)));
            if v[y * nx + x] > predicted + excess(predicted) {
                out[y * nx + x] = predicted;
            }
        }
    }
    out
}

impl Mesh {
    fn grid(w: usize, h: usize, tile: usize) -> (usize, usize, f64, f64) {
        let nx = ((w + tile / 2) / tile).max(1);
        let ny = ((h + tile / 2) / tile).max(1);
        (nx, ny, w as f64 / nx as f64, h as f64 / ny as f64)
    }

    fn bounds(t: usize, pitch: f64, limit: usize) -> (usize, usize) {
        (
            (t as f64 * pitch) as usize,
            (((t + 1) as f64 * pitch) as usize).min(limit),
        )
    }

    pub fn build(
        s: &[f32],
        w: usize,
        h: usize,
        tile: usize,
        lag: usize,
        valid: Option<&[bool]>,
    ) -> Self {
        let tile = tile.max(8);
        let (nx, ny, tw, th) = Self::grid(w, h, tile);
        let stride = (tile / 20).max(1);
        let mut bg = vec![f32::NAN; nx * ny];
        let mut noise = vec![f32::NAN; nx * ny];
        let mut buf: Vec<f32> = Vec::with_capacity((tile / stride + 2).pow(2));
        let mut dif: Vec<f32> = Vec::with_capacity(buf.capacity());
        let lag = lag.clamp(1, (w.min(h) - 1).max(1));
        for ty in 0..ny {
            let (y0, y1) = ((ty as f64 * th) as usize, (((ty + 1) as f64) * th) as usize);
            for tx in 0..nx {
                let (x0, x1) = ((tx as f64 * tw) as usize, (((tx + 1) as f64) * tw) as usize);
                buf.clear();
                dif.clear();
                let mut y = y0 + stride / 2;
                while y < y1.min(h) {
                    let mut x = x0 + stride / 2;
                    while x < x1.min(w) {
                        let i = y * w + x;
                        if valid.is_none_or(|v| v[i]) {
                            buf.push(s[i]);
                            let px = if x + lag < w { x + lag } else { x - lag.min(x) };
                            let py = if y + lag < h { y + lag } else { y - lag.min(y) };
                            let j = py * w + px;
                            if j != i && valid.is_none_or(|v| v[j]) {
                                dif.push(s[i] - s[j]);
                            }
                        }
                        x += stride;
                    }
                    y += stride;
                }
                if let Some(m) = clipped_median(&mut buf) {
                    bg[ty * nx + tx] = m;
                }
                if let Some(sg) = pair_noise(&mut dif) {
                    noise[ty * nx + tx] = sg;
                }
            }
        }
        Self::from_tiles(w, h, nx, ny, bg, noise)
    }

    fn from_tiles(
        w: usize,
        h: usize,
        nx: usize,
        ny: usize,
        mut bg: Vec<f32>,
        mut noise: Vec<f32>,
    ) -> Self {
        let (tw, th) = (w as f64 / nx as f64, h as f64 / ny as f64);
        let typical = |v: &[f32], positive: bool| -> f32 {
            let mut ok: Vec<f32> = v
                .iter()
                .copied()
                .filter(|x| x.is_finite() && (!positive || *x > 0.0))
                .collect();
            if ok.is_empty() {
                return if positive { 1.0 } else { 0.0 };
            }
            let m = ok.len() / 2;
            *ok.select_nth_unstable_by(m, f32::total_cmp).1
        };
        let (tb, tn) = (typical(&bg, false), typical(&noise, true));
        for (b, n) in bg.iter_mut().zip(noise.iter_mut()) {
            if !b.is_finite() {
                *b = tb;
            }
            if !(n.is_finite() && *n > 0.0) {
                *n = tn;
            }
        }
        let bg = despike(&bg, nx, ny, |_| 0.5 * tn);
        let noise = despike(&noise, nx, ny, |p| 0.3 * p);
        let axes = |n: usize, pitch: f64, tiles: usize| -> Vec<(u32, u32, f32)> {
            (0..n)
                .map(|v| {
                    let (a, b, t) = Self::axis(v as f64, pitch, tiles);
                    (a as u32, b as u32, t)
                })
                .collect()
        };
        Self {
            bg,
            noise,
            nx,
            ny,
            tw,
            th,
            cols: axes(w, tw, nx),
            rows: axes(h, th, ny),
        }
    }

    #[inline]
    pub fn at_px(&self, x: usize, y: usize) -> (f32, f32) {
        let (x0, x1, tx) = self.cols[x];
        let (y0, y1, ty) = self.rows[y];
        let (x0, x1, y0, y1) = (x0 as usize, x1 as usize, y0 as usize, y1 as usize);
        let lerp = |v: &[f32], tx: f32, ty: f32| {
            let a = v[y0 * self.nx + x0] * (1.0 - tx) + v[y0 * self.nx + x1] * tx;
            let b = v[y1 * self.nx + x0] * (1.0 - tx) + v[y1 * self.nx + x1] * tx;
            a * (1.0 - ty) + b * ty
        };
        (
            lerp(&self.bg, tx, ty),
            lerp(&self.noise, tx.clamp(0.0, 1.0), ty.clamp(0.0, 1.0)),
        )
    }

    pub fn row(&self, y: usize, bg: &mut [f32], noise: &mut [f32]) {
        let (y0, y1, ty) = self.rows[y];
        let (y0, y1) = (y0 as usize, y1 as usize);
        let tyc = ty.clamp(0.0, 1.0);
        let w = bg.len();
        let mut x = 0;
        while x < w {
            let (x0, x1, _) = self.cols[x];
            let mut end = x + 1;
            while end < w && self.cols[end].0 == x0 {
                end += 1;
            }
            let (x0, x1) = (x0 as usize, x1 as usize);
            let corners = |v: &[f32]| {
                (
                    v[y0 * self.nx + x0],
                    v[y0 * self.nx + x1],
                    v[y1 * self.nx + x0],
                    v[y1 * self.nx + x1],
                )
            };
            let (b00, b01, b10, b11) = corners(&self.bg);
            for (o, c) in bg[x..end].iter_mut().zip(&self.cols[x..end]) {
                let tx = c.2;
                let a = b00 * (1.0 - tx) + b01 * tx;
                let b = b10 * (1.0 - tx) + b11 * tx;
                *o = a * (1.0 - ty) + b * ty;
            }
            let (n00, n01, n10, n11) = corners(&self.noise);
            for (o, c) in noise[x..end].iter_mut().zip(&self.cols[x..end]) {
                let tx = c.2.clamp(0.0, 1.0);
                let a = n00 * (1.0 - tx) + n01 * tx;
                let b = n10 * (1.0 - tx) + n11 * tx;
                *o = a * (1.0 - tyc) + b * tyc;
            }
            x = end;
        }
    }

    #[inline]
    fn axis(v: f64, pitch: f64, n: usize) -> (usize, usize, f32) {
        if n == 1 {
            return (0, 0, 0.0);
        }
        let f = (v + 0.5) / pitch - 0.5;
        let i = (f.max(0.0) as usize).min(n - 2);
        (i, i + 1, (f - i as f64) as f32)
    }

    #[inline]
    pub fn at(&self, x: f64, y: f64) -> (f32, f32) {
        let (x0, x1, tx) = Self::axis(x, self.tw, self.nx);
        let (y0, y1, ty) = Self::axis(y, self.th, self.ny);
        let lerp = |v: &[f32], tx: f32, ty: f32| {
            let a = v[y0 * self.nx + x0] * (1.0 - tx) + v[y0 * self.nx + x1] * tx;
            let b = v[y1 * self.nx + x0] * (1.0 - tx) + v[y1 * self.nx + x1] * tx;
            a * (1.0 - ty) + b * ty
        };
        (
            lerp(&self.bg, tx, ty),
            lerp(&self.noise, tx.clamp(0.0, 1.0), ty.clamp(0.0, 1.0)),
        )
    }

    pub fn tiles(&self) -> (usize, usize) {
        (self.nx, self.ny)
    }

    pub fn span(&self, t: usize, horizontal: bool, limit: usize) -> (usize, usize) {
        let p = if horizontal { self.tw } else { self.th };
        Self::bounds(t, p, limit)
    }

    pub fn floor(&self, tx: usize, ty: usize, k: f32, w: usize, h: usize) -> f32 {
        let (x0, x1) = self.span(tx, true, w);
        let (y0, y1) = self.span(ty, false, h);
        let xs = [
            x0 as f64,
            (tx as f64 + 0.5) * self.tw - 0.5,
            x1 as f64 - 1.0,
        ];
        let ys = [
            y0 as f64,
            (ty as f64 + 0.5) * self.th - 0.5,
            y1 as f64 - 1.0,
        ];
        let mut lo = f32::INFINITY;
        let mut noise = f32::INFINITY;
        for y in ys {
            for x in xs {
                let (b, n) = self.at(x, y);
                lo = lo.min(b + k * n);
                noise = noise.min(n);
            }
        }
        lo - 1e-3 * noise
    }

    pub fn raise_noise(&mut self, model: &Self) {
        for ty in 0..self.ny {
            for tx in 0..self.nx {
                let x = (tx as f64 + 0.5) * self.tw - 0.5;
                let y = (ty as f64 + 0.5) * self.th - 0.5;
                let i = ty * self.nx + tx;
                self.noise[i] = self.noise[i].max(model.at(x, y).1);
            }
        }
    }

    pub fn typical(&self) -> (f32, f32) {
        let med = |v: &[f32]| {
            let mut c = v.to_vec();
            let m = c.len() / 2;
            *c.select_nth_unstable_by(m, f32::total_cmp).1
        };
        (med(&self.bg), med(&self.noise))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fayde_core::Portable;

    pub(crate) fn normals(n: usize, seed: u64) -> Vec<f32> {
        use fayde_core::Portable;
        let mut s = seed;
        let mut uni = move || {
            s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((s >> 11) as f64 + 0.5) / (1u64 << 53) as f64
        };
        (0..n)
            .map(|_| {
                let (u, v) = (uni(), uni());
                ((-2.0 * u.pln()).sqrt() * (2.0 * std::f64::consts::PI * v).pcos()) as f32
            })
            .collect()
    }

    #[test]
    fn robust_statistics_ignore_a_bright_tail() {
        let mut v = normals(4000, 3);
        for x in v.iter_mut().take(400) {
            *x += 50.0;
        }
        let med = clipped_median(&mut v).unwrap();
        assert!(med.abs() < 0.1, "{med}");
        let a = normals(4000, 5);
        let b = normals(4000, 6);
        let mut d: Vec<f32> = a.iter().zip(&b).map(|(x, y)| x - y).collect();
        for x in d.iter_mut().take(400) {
            *x += 50.0;
        }
        let sig = pair_noise(&mut d).unwrap();
        assert!((sig - 1.0).abs() < 0.1, "{sig}");
    }

    #[test]
    fn censored_sky_recovers_the_level_below_the_floor() {
        let (w, h) = (256, 256);
        let normal = normals(w * h, 41);
        for (floor, below) in [(0.0f32, 2.0), (-20.0, 7.0)] {
            let raw: Vec<_> = normal
                .iter()
                .map(|&z| (floor - below + 5.0 * z).round())
                .collect();
            assert!(ClippedBackground::estimate(&raw, w, h, 64).is_none());
            let mut clipped: Vec<_> = raw.iter().map(|&v| v.max(floor)).collect();
            let sky: Vec<_> = clipped.clone();
            for value in clipped.iter_mut().step_by(200) {
                *value += 50.0;
            }
            for value in clipped.iter_mut().skip(1).step_by(100) {
                *value = f32::NAN;
            }
            let model = ClippedBackground::estimate(&clipped, w, h, 64).unwrap();
            let (mu, sigma) = model.latent.typical();
            assert!((mu - (floor - below)).abs() < 0.8, "{mu}");
            assert!((sigma - 5.0).abs() < 0.7, "{sigma}");

            let (bg, noise) = model.filtered(w, h, 1.0, 0.0).typical();
            let mean = sky.iter().sum::<f32>() / sky.len() as f32;
            let variance = sky.iter().map(|v| (v - mean).ppowi(2)).sum::<f32>() / sky.len() as f32;
            assert!((bg - mean).abs() < 0.2, "{bg} vs {mean}");
            assert!(
                (noise - variance.sqrt()).abs() < 0.3,
                "{noise} vs {variance}"
            );
        }
    }

    #[test]
    fn contiguous_padding_is_distinct_from_clipped_sky() {
        let (w, h) = (256, 192);
        let normal = normals(w * h, 43);
        let mut padded: Vec<_> = normal.iter().map(|&z| 500.0 + 5.0 * z).collect();
        for row in padded.chunks_mut(w) {
            row[..89].fill(0.0);
        }
        assert!(ClippedBackground::estimate(&padded, w, h, 64).is_none());

        let mut clipped: Vec<_> = normal.iter().map(|&z| (5.0 * z - 2.0).max(0.0)).collect();
        let model = ClippedBackground::estimate(&clipped, w, h, 64).unwrap();
        assert!(model.valid_pixels(&clipped, w).is_none());

        for row in clipped.chunks_mut(w) {
            row[..89].fill(0.0);
        }
        let model = ClippedBackground::estimate(&clipped, w, h, 64).unwrap();
        let valid = model.valid_pixels(&clipped, w).unwrap();
        for row in valid.chunks(w) {
            assert!(row[..64].iter().all(|&ok| !ok));
            assert!(row[128..].iter().all(|&ok| ok));
        }
    }

    #[test]
    fn a_gradient_is_followed_and_noise_recovered() {
        let (w, h) = (512, 384);
        let n = normals(w * h, 11);
        let img: Vec<f32> = (0..w * h)
            .map(|i| 100.0 + 0.2 * (i % w) as f32 + 3.0 * n[i])
            .collect();
        let m = Mesh::build(&img, w, h, 64, 5, None);
        for (x, y) in [(40.0, 40.0), (256.0, 190.0), (470.0, 340.0), (511.0, 2.0)] {
            let (b, s) = m.at(x, y);
            assert!((b - (100.0 + 0.2 * x as f32)).abs() < 1.5, "{x}: {b}");
            assert!((2.75..3.25).contains(&s), "{s}");
        }
        let (nx, ny) = m.tiles();
        for ty in 0..ny {
            for tx in 0..nx {
                let (x0, x1) = m.span(tx, true, w);
                let (y0, y1) = m.span(ty, false, h);
                let f = m.floor(tx, ty, 5.0, w, h);
                for y in y0..y1 {
                    for x in x0..x1 {
                        let (b, s) = m.at(x as f64, y as f64);
                        assert!(f <= b + 5.0 * s, "tile ({tx}, {ty}) pixel ({x}, {y})");
                    }
                }
            }
        }
    }

    #[test]
    fn a_bright_tile_is_replaced_and_a_curved_background_is_kept() {
        let (nx, ny) = (8, 6);
        let glow =
            |x: usize, y: usize| 1000.0 * 0.5f32.ppowf(((x * x + y * y) as f32).sqrt() / 2.0);
        let mut v: Vec<f32> = (0..nx * ny).map(|i| glow(i % nx, i / nx)).collect();
        v[3 * nx + 4] += 40.0;
        let out = despike(&v, nx, ny, |_| 2.0);
        for (i, got) in out.iter().enumerate() {
            let want = glow(i % nx, i / nx);
            let tol = if i == 3 * nx + 4 { 25.0 } else { 1e-3 };
            assert!((got - want).abs() <= tol, "tile {i}: {got} vs {want}");
        }
        assert!(out[3 * nx + 4] < v[3 * nx + 4] - 15.0);
    }

    #[test]
    fn constant_and_tiny_images_do_not_break() {
        let m = Mesh::build(&[5.0; 100], 10, 10, 64, 5, None);
        let (b, s) = m.at(3.0, 3.0);
        assert_eq!(b, 5.0);
        assert!(s > 0.0);
        let m = Mesh::build(&[1.0, 2.0, 3.0, 4.0], 2, 2, 64, 5, None);
        assert!(m.at(0.0, 0.0).1 > 0.0);
    }
}
