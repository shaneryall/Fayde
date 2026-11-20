use fayde_core::Portable;

pub const PATTERN_TILE: usize = 32;
const PATTERN_CLIP: f32 = 3.0;
const PATTERN_MIN_PER_TERM: usize = 8;
const PEAK_OVER_MEDIAN: f64 = 30.0;
const PEAK_MIN_SHARE: f64 = 0.005;
const MAX_FOUND: usize = 4;
const FIND_MAX_SIDE: usize = 512;

#[derive(Debug, Clone, Default, PartialEq)]
pub struct Removed {
    pub rms: f64,
    pub amplitude: Vec<f64>,
    pub tiles: usize,
    pub skipped: usize,
}

fn median_sigma(v: &mut [f32]) -> Option<(f32, f32)> {
    if v.is_empty() {
        return None;
    }
    let mid = v.len() / 2;
    let med = *v.select_nth_unstable_by(mid, f32::total_cmp).1;
    for x in v.iter_mut() {
        *x = (*x - med).abs();
    }
    let mad = *v.select_nth_unstable_by(mid, f32::total_cmp).1;
    Some((med, 1.4826 * mad))
}

fn cholesky_solve(a: &mut [f64], b: &mut [f64], n: usize) -> Option<()> {
    for j in 0..n {
        let mut d = a[j * n + j];
        for k in 0..j {
            d -= a[j * n + k] * a[j * n + k];
        }
        if d.is_nan() || d <= 1e-12 * a[j * n + j].abs().max(1e-300) {
            return None;
        }
        let d = d.sqrt();
        a[j * n + j] = d;
        for i in j + 1..n {
            let mut s = a[i * n + j];
            for k in 0..j {
                s -= a[i * n + k] * a[j * n + k];
            }
            a[i * n + j] = s / d;
        }
    }
    for i in 0..n {
        let mut s = b[i];
        for k in 0..i {
            s -= a[i * n + k] * b[k];
        }
        b[i] = s / a[i * n + i];
    }
    for i in (0..n).rev() {
        let mut s = b[i];
        for k in i + 1..n {
            s -= a[k * n + i] * b[k];
        }
        b[i] = s / a[i * n + i];
    }
    Some(())
}

fn on_grid(f: f64) -> bool {
    (2.0 * f).fract() == 0.0
}

fn phase_table(f: f64, n: usize) -> (Vec<f64>, Vec<f64>) {
    if on_grid(f) {
        let alt = (2.0 * f) as i64 % 2 != 0;
        return (0..n)
            .map(|i| (if alt && i % 2 == 1 { -1.0 } else { 1.0 }, 0.0))
            .unzip();
    }
    let tau = 2.0 * std::f64::consts::PI;
    (0..n)
        .map(|i| (tau * f * i as f64).psin_cos())
        .map(|(s, c)| (c, s))
        .unzip()
}

struct Waves {
    x: Vec<(Vec<f64>, Vec<f64>)>,
    y: Vec<(Vec<f64>, Vec<f64>)>,
    cols: Vec<usize>,
}

impl Waves {
    fn new(freqs: &[[f64; 2]], w: usize, h: usize) -> Self {
        let cols = freqs
            .iter()
            .enumerate()
            .flat_map(|(j, f)| {
                let has_sine = !(on_grid(f[0]) && on_grid(f[1]));
                std::iter::once(2 * j).chain(has_sine.then_some(2 * j + 1))
            })
            .collect();
        Self {
            x: freqs.iter().map(|f| phase_table(f[0], w)).collect(),
            y: freqs.iter().map(|f| phase_table(f[1], h)).collect(),
            cols,
        }
    }

    fn len(&self) -> usize {
        self.x.len()
    }

    fn at(&self, j: usize, x: usize, y: usize) -> (f64, f64) {
        let (cx, sx) = (self.x[j].0[x], self.x[j].1[x]);
        let (cy, sy) = (self.y[j].0[y], self.y[j].1[y]);
        (cx * cy - sx * sy, sx * cy + cx * sy)
    }
}

fn tile_bounds(tile: usize, size: usize) -> (usize, usize) {
    (tile * PATTERN_TILE, ((tile + 1) * PATTERN_TILE).min(size))
}

fn source_mask(data: &[f32], w: usize, h: usize) -> Vec<bool> {
    let (nx, ny) = (w.div_ceil(PATTERN_TILE), h.div_ceil(PATTERN_TILE));
    let mut source = vec![false; w * h];
    let mut buf: Vec<f32> = Vec::with_capacity(PATTERN_TILE * PATTERN_TILE);
    for ty_ in 0..ny {
        for tx_ in 0..nx {
            let (x0, x1) = tile_bounds(tx_, w);
            let (y0, y1) = tile_bounds(ty_, h);
            buf.clear();
            for y in y0..y1 {
                buf.extend(
                    data[y * w + x0..y * w + x1]
                        .iter()
                        .filter(|v| v.is_finite()),
                );
            }
            let Some((med, sig)) = median_sigma(&mut buf) else {
                continue;
            };
            for y in y0..y1 {
                for x in x0..x1 {
                    let v = data[y * w + x];
                    if !v.is_finite() || (v - med).abs() > PATTERN_CLIP * sig {
                        source[y * w + x] = true;
                    }
                }
            }
        }
    }
    let mut masked = source.clone();
    for y in 0..h {
        for x in 0..w {
            if source[y * w + x] {
                for yy in y.saturating_sub(2)..=(y + 2).min(h - 1) {
                    for xx in x.saturating_sub(2)..=(x + 2).min(w - 1) {
                        masked[yy * w + xx] = true;
                    }
                }
            }
        }
    }

    masked
}

fn fit_tile_amplitudes(
    data: &[f32],
    masked: &[bool],
    w: usize,
    h: usize,
    waves: &Waves,
) -> Vec<Option<Vec<f64>>> {
    let (nx, ny) = (w.div_ceil(PATTERN_TILE), h.div_ceil(PATTERN_TILE));
    let k = waves.len();
    let np = 1 + waves.cols.len();
    let mut amp: Vec<Option<Vec<f64>>> = Vec::with_capacity(nx * ny);
    let mut basis = vec![0.0f64; np];
    let mut wave_terms = vec![0.0f64; 2 * k];
    for ty_ in 0..ny {
        for tx_ in 0..nx {
            let (x0, x1) = tile_bounds(tx_, w);
            let (y0, y1) = tile_bounds(ty_, h);
            let mut a = vec![0.0f64; np * np];
            let mut b = vec![0.0f64; np];
            let mut n = 0usize;
            for y in y0..y1 {
                for x in x0..x1 {
                    if masked[y * w + x] {
                        continue;
                    }
                    basis[0] = 1.0;
                    for j in 0..k {
                        let (c, s) = waves.at(j, x, y);
                        wave_terms[2 * j] = c;
                        wave_terms[2 * j + 1] = s;
                    }
                    for (bc, &col) in basis[1..].iter_mut().zip(&waves.cols) {
                        *bc = wave_terms[col];
                    }
                    let v = f64::from(data[y * w + x]);
                    for r in 0..np {
                        for c in 0..=r {
                            a[r * np + c] += basis[r] * basis[c];
                        }
                        b[r] += basis[r] * v;
                    }
                    n += 1;
                }
            }
            for r in 0..np {
                for c in r + 1..np {
                    a[r * np + c] = a[c * np + r];
                }
            }
            let fitted =
                n >= PATTERN_MIN_PER_TERM * np && cholesky_solve(&mut a, &mut b, np).is_some();
            amp.push(fitted.then(|| {
                let mut t = vec![0.0f64; 2 * k];
                for (&col, v) in waves.cols.iter().zip(&b[1..]) {
                    t[col] = *v;
                }
                t
            }));
        }
    }
    amp
}

fn subtract_interpolated(
    data: &mut [f32],
    w: usize,
    h: usize,
    amp: &[Vec<f64>],
    waves: &Waves,
) -> f64 {
    let (nx, ny) = (w.div_ceil(PATTERN_TILE), h.div_ceil(PATTERN_TILE));
    let k = waves.len();
    let centre = |t: usize, n: usize| {
        (t * PATTERN_TILE) as f64 + 0.5 * ((PATTERN_TILE.min(n - t * PATTERN_TILE)) as f64 - 1.0)
    };
    let locate = |p: usize, n_t: usize, size: usize| -> (usize, usize, f64) {
        let pf = p as f64;
        let first = centre(0, size);
        let last = centre(n_t - 1, size);
        if n_t == 1 || pf <= first {
            return (0, 0, 0.0);
        }
        if pf >= last {
            return (n_t - 1, n_t - 1, 0.0);
        }
        let mut t = p / PATTERN_TILE;
        if pf < centre(t, size) {
            t -= 1;
        }
        let (c0, c1) = (centre(t, size), centre(t + 1, size));
        (t, t + 1, (pf - c0) / (c1 - c0))
    };
    let (mut ss, mut count) = (0.0f64, 0usize);
    let mut terms = vec![0.0f64; 2 * k];
    for y in 0..h {
        let (ya, yb, fy) = locate(y, ny, h);
        for x in 0..w {
            let i = y * w + x;
            if !data[i].is_finite() {
                continue;
            }
            let (xa, xb, fx) = locate(x, nx, w);
            let (a00, a01, a10, a11) = (
                &amp[ya * nx + xa],
                &amp[ya * nx + xb],
                &amp[yb * nx + xa],
                &amp[yb * nx + xb],
            );
            for (t, term) in terms.iter_mut().enumerate() {
                *term = (1.0 - fy) * ((1.0 - fx) * a00[t] + fx * a01[t])
                    + fy * ((1.0 - fx) * a10[t] + fx * a11[t]);
            }
            let mut p = 0.0;
            for j in 0..k {
                let (c, s) = waves.at(j, x, y);
                p += terms[2 * j] * c + terms[2 * j + 1] * s;
            }
            data[i] -= p as f32;
            ss += p * p;
            count += 1;
        }
    }
    (ss / count.max(1) as f64).sqrt()
}

pub fn remove(data: &mut [f32], w: usize, h: usize, freqs: &[[f64; 2]]) -> Removed {
    let k = freqs.len();
    if k == 0 || w == 0 || h == 0 {
        return Removed::default();
    }
    let waves = Waves::new(freqs, w, h);
    let amp = {
        let masked = source_mask(data, w, h);
        fit_tile_amplitudes(data, &masked, w, h, &waves)
    };
    let total_tiles = amp.len();
    let tiles = amp.iter().filter(|a| a.is_some()).count();
    if tiles == 0 {
        return Removed {
            skipped: total_tiles,
            ..Removed::default()
        };
    }
    let mean: Vec<f64> = (0..2 * k)
        .map(|i| amp.iter().flatten().map(|a| a[i]).sum::<f64>() / tiles as f64)
        .collect();
    let amp: Vec<Vec<f64>> = amp
        .into_iter()
        .map(|a| a.unwrap_or_else(|| mean.clone()))
        .collect();

    let rms = subtract_interpolated(data, w, h, &amp, &waves);
    Removed {
        rms,
        amplitude: (0..k)
            .map(|j| {
                amp.iter()
                    .map(|a| a[2 * j].phypot(a[2 * j + 1]))
                    .sum::<f64>()
                    / amp.len() as f64
            })
            .collect(),
        tiles,
        skipped: total_tiles - tiles,
    }
}

fn fft(re: &mut [f64], im: &mut [f64]) {
    let n = re.len();
    let mut j = 0usize;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j |= bit;
        if i < j {
            re.swap(i, j);
            im.swap(i, j);
        }
    }
    let mut len = 2;
    while len <= n {
        let (s, c) = (-2.0 * std::f64::consts::PI / len as f64).psin_cos();
        for start in (0..n).step_by(len) {
            let (mut wr, mut wi) = (1.0f64, 0.0f64);
            for k in 0..len / 2 {
                let (a, b) = (start + k, start + k + len / 2);
                let (tr, ti) = (re[b] * wr - im[b] * wi, re[b] * wi + im[b] * wr);
                re[b] = re[a] - tr;
                im[b] = im[a] - ti;
                re[a] += tr;
                im[a] += ti;
                (wr, wi) = (wr * c - wi * s, wr * s + wi * c);
            }
        }
        len <<= 1;
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Found {
    pub frequency: [f64; 2],
    pub share: f64,
}

pub fn find(frames: &[(&[f32], usize, usize)]) -> Vec<Found> {
    let mut f = Finder::default();
    for &(data, w, h) in frames {
        f.add(data, w, h);
    }
    f.finish().0
}

#[derive(Debug, Clone, Default)]
pub struct Finder {
    n: usize,
    mean: Vec<f64>,
    used: usize,
}

impl Finder {
    pub fn add(&mut self, data: &[f32], w: usize, h: usize) {
        if self.n == 0 {
            let s = w.min(h);
            if s < 64 {
                return;
            }
            self.n = 1usize << (usize::BITS - 1 - s.min(FIND_MAX_SIDE).leading_zeros());
            self.mean = vec![0.0; self.n * self.n];
        }
        let n = self.n;
        if w.min(h) < n {
            return;
        }
        let (x0, y0) = ((w - n) / 2, (h - n) / 2);
        let mut crop: Vec<f32> = (0..n)
            .flat_map(|y| {
                data[(y0 + y) * w + x0..(y0 + y) * w + x0 + n]
                    .iter()
                    .copied()
            })
            .collect();
        let mut tmp: Vec<f32> = crop.iter().copied().filter(|v| v.is_finite()).collect();
        let Some((med, sig)) = median_sigma(&mut tmp) else {
            return;
        };
        for v in crop.iter_mut() {
            *v = if v.is_finite() && (*v - med).abs() <= PATTERN_CLIP * sig {
                *v - med
            } else {
                0.0
            };
        }
        let mut re: Vec<f64> = crop.iter().map(|&v| f64::from(v)).collect();
        let mut im = vec![0.0f64; n * n];
        let (mut rr, mut ri) = (vec![0.0; n], vec![0.0; n]);
        for y in 0..n {
            rr.copy_from_slice(&re[y * n..(y + 1) * n]);
            ri.copy_from_slice(&im[y * n..(y + 1) * n]);
            fft(&mut rr, &mut ri);
            re[y * n..(y + 1) * n].copy_from_slice(&rr);
            im[y * n..(y + 1) * n].copy_from_slice(&ri);
        }
        for x in 0..n {
            for y in 0..n {
                rr[y] = re[y * n + x];
                ri[y] = im[y * n + x];
            }
            fft(&mut rr, &mut ri);
            for y in 0..n {
                re[y * n + x] = rr[y];
                im[y * n + x] = ri[y];
            }
        }
        let power: Vec<f64> = re.iter().zip(&im).map(|(a, b)| a * a + b * b).collect();
        let total: f64 = power.iter().skip(1).sum();
        if total <= 0.0 {
            return;
        }
        for (m, p) in self.mean.iter_mut().zip(&power) {
            *m += p / total;
        }
        self.used += 1;
    }

    pub fn finish(self) -> (Vec<Found>, usize) {
        let (n, used, mut mean) = (self.n, self.used, self.mean);
        if used == 0 {
            return (Vec::new(), 0);
        }
        for m in mean.iter_mut() {
            *m /= used as f64;
        }
        mean[0] = 0.0;
        let signed = |i: usize| {
            if i <= n / 2 {
                i as i64
            } else {
                i as i64 - n as i64
            }
        };
        let at = |kx: i64, ky: i64| {
            mean[(ky.rem_euclid(n as i64) as usize) * n + kx.rem_euclid(n as i64) as usize]
        };
        let mut sorted: Vec<f64> = mean.iter().copied().filter(|&v| v > 0.0).collect();
        let mid = sorted.len() / 2;
        let median = if sorted.is_empty() {
            0.0
        } else {
            *sorted.select_nth_unstable_by(mid, f64::total_cmp).1
        };
        let min_k = (0.02 * n as f64).ceil() as i64;
        let mut found: Vec<Found> = Vec::new();
        for iy in 0..n {
            for ix in 0..n {
                let (kx, ky) = (signed(ix), signed(iy));
                if !(kx > 0 || (kx == 0 && ky > 0)) || kx * kx + ky * ky < min_k * min_k {
                    continue;
                }
                let p = at(kx, ky);
                if p < PEAK_OVER_MEDIAN * median {
                    continue;
                }
                let local_max = (-1..=1)
                    .all(|dy| (-1..=1).all(|dx| (dx == 0 && dy == 0) || at(kx + dx, ky + dy) < p));
                if !local_max {
                    continue;
                }
                let share = 2.0
                    * (-1..=1)
                        .map(|dy| (-1..=1).map(|dx| at(kx + dx, ky + dy)).sum::<f64>())
                        .sum::<f64>();
                if share < PEAK_MIN_SHARE {
                    continue;
                }
                let refine = |l: f64, c: f64, r: f64| {
                    let (l, c, r) = (
                        l.max(1e-300).pln(),
                        c.max(1e-300).pln(),
                        r.max(1e-300).pln(),
                    );
                    let den = l - 2.0 * c + r;
                    if den < 0.0 {
                        (0.5 * (l - r) / den).clamp(-0.5, 0.5)
                    } else {
                        0.0
                    }
                };
                let ox = refine(at(kx - 1, ky), p, at(kx + 1, ky));
                let oy = refine(at(kx, ky - 1), p, at(kx, ky + 1));
                found.push(Found {
                    frequency: [(kx as f64 + ox) / n as f64, (ky as f64 + oy) / n as f64],
                    share,
                });
            }
        }
        found.sort_by(|a, b| b.share.total_cmp(&a.share));
        found.truncate(MAX_FOUND);
        (found, used)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::synth::Rng;

    fn frame(w: usize, h: usize, f: [f64; 2], amp: f64, seed: u64) -> (Vec<f32>, Vec<f32>) {
        let mut rng = Rng(seed);
        let tau = 2.0 * std::f64::consts::PI;
        let mut clean = vec![0.0f32; w * h];
        let mut data = vec![0.0f32; w * h];
        let stars: Vec<(f64, f64, f64)> = (0..30)
            .map(|_| {
                (
                    rng.uniform() * w as f64,
                    rng.uniform() * h as f64,
                    3000.0 * rng.uniform(),
                )
            })
            .collect();
        for y in 0..h {
            for x in 0..w {
                let mut v = 100.0 + 10.0 * rng.normal();
                for &(sx, sy, fl) in &stars {
                    let d2 = (x as f64 - sx).ppowi(2) + (y as f64 - sy).ppowi(2);
                    if d2 < 64.0 {
                        v += fl / (2.0 * std::f64::consts::PI * 2.0) * (-d2 / 4.0).pexp();
                    }
                }
                clean[y * w + x] = v as f32;
                let drift = 0.3 * (x as f64 / w as f64);
                let p = amp * (tau * (f[0] * x as f64 + f[1] * y as f64) + drift).pcos();
                data[y * w + x] = (v + p) as f32;
            }
        }
        (data, clean)
    }

    #[test]
    fn a_wave_is_found_and_removed() {
        let (w, h) = (300, 260);
        let f = [0.1706, 0.1608];
        let frames: Vec<(Vec<f32>, Vec<f32>)> =
            (0..3).map(|s| frame(w, h, f, 8.0, 40 + s)).collect();
        let views: Vec<(&[f32], usize, usize)> =
            frames.iter().map(|(d, _)| (d.as_slice(), w, h)).collect();
        let found = find(&views);
        assert!(!found.is_empty(), "nothing found");
        let g = found[0].frequency;
        assert!(
            (g[0] - f[0]).abs() < 0.002 && (g[1] - f[1]).abs() < 0.002,
            "{found:?}"
        );
        assert!(found[0].share > 0.2, "{found:?}");

        let (mut data, clean) = frames[0].clone();
        let r = remove(&mut data, w, h, &[g]);
        assert!((r.amplitude[0] - 8.0).abs() < 0.8, "{r:?}");
        let resid: Vec<f64> = data
            .iter()
            .zip(&clean)
            .map(|(a, b)| f64::from(a - b))
            .collect();
        let rms = (resid.iter().map(|v| v * v).sum::<f64>() / resid.len() as f64).sqrt();
        assert!(rms < 1.5, "residual pattern rms {rms}");

        let (mut data, clean) = frames[1].clone();
        for (i, v) in data.iter_mut().enumerate() {
            *v += if (i / w) % 2 == 1 { 3.0 } else { -3.0 };
        }
        let r = remove(&mut data, w, h, &[g, [0.0, 0.5]]);
        assert!((r.amplitude[1] - 3.0).abs() < 0.5, "{r:?}");
        let rms = (data
            .iter()
            .zip(&clean)
            .map(|(a, b)| f64::from(a - b).ppowi(2))
            .sum::<f64>()
            / data.len() as f64)
            .sqrt();
        assert!(rms < 1.5, "residual rms {rms}");

        let (plain, _) = frame(w, h, f, 0.0, 7);
        assert!(find(&[(plain.as_slice(), w, h)]).is_empty());
        let mut copy = plain.clone();
        let r = remove(&mut copy, w, h, &[f]);
        assert!(r.rms < 1.0, "{r:?}");
    }

    #[test]
    fn fft_matches_a_direct_transform() {
        let n = 16;
        let x: Vec<f64> = (0..n).map(|i| ((i * 7 + 3) % 11) as f64 - 5.0).collect();
        let (mut re, mut im) = (x.clone(), vec![0.0; n]);
        fft(&mut re, &mut im);
        for k in 0..n {
            let (mut sr, mut si) = (0.0, 0.0);
            for (t, v) in x.iter().enumerate() {
                let (s, c) = (-2.0 * std::f64::consts::PI * (k * t) as f64 / n as f64).psin_cos();
                sr += v * c;
                si += v * s;
            }
            assert!((re[k] - sr).abs() < 1e-9 && (im[k] - si).abs() < 1e-9);
        }
    }
}
