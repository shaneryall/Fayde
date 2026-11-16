//! Background and noise of the filtered image on a coarse mesh.
//!
//! The image is cut into tiles. Each tile gives a background level: the median of a sample of
//! its pixels, after one clip of the bright tail that stars add. It also gives a noise level,
//! taken from the *differences* of pixel pairs a fixed lag apart: a difference cancels the
//! background and any gradient across the tile, which a spread of values would count as noise.
//! The lag is longer than the filter kernel, so the two pixels of a pair carry independent
//! noise and the noise is `1.4826 × MAD / √2` of the differences. Tiles sitting on a bright
//! star or a galaxy are replaced by what their neighbours predict (`despike`), and values
//! between tile centres are interpolated bilinearly.
//!
//! The mesh is SExtractor's (Bertin & Arnouts 1996, A&AS 117, 393, §3);
//! noise from pixel-pair differences is astrometry.net's (`dsigma`, Lang et al. 2010); the
//! median absolute deviation with its factor 1.4826 is Hampel 1974 (JASA 69, 383).
//!
//! The statistics are taken on the *filtered* image, so the noise is that of the quantity the
//! peak test thresholds, whatever the noise correlation of the input (binned, resampled or
//! debayered frames included).

/// Background and noise on a tile grid.
#[derive(Debug, Clone)]
pub(crate) struct Mesh {
    nx: usize,
    ny: usize,
    /// Tile pitch, pixels (not necessarily an integer).
    tw: f64,
    th: f64,
    bg: Vec<f32>,
    noise: Vec<f32>,
}

/// Median of `buf` after one clip of its bright tail at median + 3 × (median − 15.87th
/// percentile). `None` for fewer than 8 values.
fn clipped_median(buf: &mut [f32]) -> Option<f32> {
    let n = buf.len();
    if n < 8 {
        return None;
    }
    let cmp = |a: &f32, b: &f32| a.total_cmp(b);
    let m = n / 2;
    let med = *buf.select_nth_unstable_by(m, cmp).1;
    // Everything below index m is ≤ the median, in no particular order.
    let q = *buf[..m]
        .select_nth_unstable_by(((0.1587 * n as f64) as usize).min(m - 1), cmp)
        .1;
    let hi = med + 3.0 * (med - q);
    // The k-th smallest value of the clipped set is the k-th smallest overall, and lies at or
    // below the first median.
    let kept = buf.iter().filter(|v| **v <= hi).count();
    let m2 = (kept / 2).min(m);
    if m2 == m {
        return Some(med);
    }
    Some(*buf[..m].select_nth_unstable_by(m2, cmp).1)
}

/// 1σ of the noise in each member of the pairs whose differences are `dif`: `1.4826 × MAD / √2`
/// after up to three clips at 3σ (pairs with a star in one member; each clip tightens the next
/// in a crowded tile). `None` for fewer than 8 values.
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
    // `top` is one past the index of the current median: the values below it are the smallest
    // `top`, in no particular order, and each later median lies among them.
    let mut top = n / 2;
    let mut mad = *dif.select_nth_unstable_by(top, cmp).1;
    let mut factor = 1.4826f32;
    for _ in 0..3 {
        // 3σ of a difference is 3 × 1.4826 × MAD.
        let hi = 3.0 * factor * mad;
        let kept = dif.iter().filter(|v| **v <= hi).count();
        let m = kept / 2;
        if m >= top || m == 0 {
            break;
        }
        mad = *dif[..top].select_nth_unstable_by(m, cmp).1;
        top = m;
        // The clip takes 0.27% of a Gaussian's values, which lowers the median of the rest
        // by 0.31%: the factor is raised to match.
        factor = 1.4872;
    }
    Some(factor * mad / std::f32::consts::SQRT_2)
}

/// Replace tiles that stand above their surroundings (a bright star, a galaxy, a crowded
/// patch: all add light) by what the surroundings predict.
///
/// The prediction for a tile is the largest, over the four lines through it (row, column and
/// both diagonals), of the mean of its two neighbours on that line. Each is exact where the
/// background is a plane. The glow of stray light curves upwards along the direction to its
/// source and downwards across it, so the lines disagree; the one nearest that direction
/// predicts at least the tile's value, and taking the largest keeps a gradient or a glow from
/// being mistaken for an outlier. Tiles on the border, which lack a neighbour on one side, are
/// kept as measured. A tile is replaced when it exceeds the prediction by more than
/// `excess(prediction)`.
///
/// SExtractor filters its mesh with a plain 3 × 3 median, which flattens the corner of a frame
/// lit from outside; this one-sided filter is Fayde's own.
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
    /// Measure `s` (`w × h`) on tiles of about `tile` pixels. `lag` is the pixel offset, in
    /// both axes, between the two members of a noise pair: more than the reach of the filter
    /// that made `s`. `valid` marks usable pixels.
    pub fn build(
        s: &[f32],
        w: usize,
        h: usize,
        tile: usize,
        lag: usize,
        valid: Option<&[bool]>,
    ) -> Self {
        let tile = tile.max(8);
        let nx = ((w + tile / 2) / tile).max(1);
        let ny = ((h + tile / 2) / tile).max(1);
        let (tw, th) = (w as f64 / nx as f64, h as f64 / ny as f64);
        // About 20 samples a side per tile: enough for a 5% noise estimate, a small fraction
        // of the pixels.
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
                            // The partner lies towards the frame's interior.
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
        // Tiles without a measurement (masked, or constant: zero noise) take the frame's
        // typical values, so nothing is detected in them on the strength of a zero threshold.
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
        // A tile's median has an error of about a tenth of the noise; half the noise above the
        // prediction is light that does not belong to the background. The noise itself is
        // replaced where it is a third above its surroundings.
        let bg = despike(&bg, nx, ny, |_| 0.5 * tn);
        let noise = despike(&noise, nx, ny, |p| 0.3 * p);
        Self {
            bg,
            noise,
            nx,
            ny,
            tw,
            th,
        }
    }

    /// The two tiles on either side of pixel coordinate `v` along an axis of `n` tiles of pitch
    /// `pitch`, and the weight of the second. Outside the outermost tile centres the weight
    /// leaves `[0, 1]`: a straight-line extrapolation.
    #[inline]
    fn axis(v: f64, pitch: f64, n: usize) -> (usize, usize, f32) {
        if n == 1 {
            return (0, 0, 0.0);
        }
        let f = (v + 0.5) / pitch - 0.5;
        let i = (f.max(0.0) as usize).min(n - 2);
        (i, i + 1, (f - i as f64) as f32)
    }

    /// `(background, noise)` at pixel `(x, y)`: bilinear between tile centres. Beyond the
    /// outermost centres the background is extrapolated in a straight line, so a gradient
    /// running off the frame is followed to its edge; the noise is held at the edge value.
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

    /// Tile grid size.
    pub fn tiles(&self) -> (usize, usize) {
        (self.nx, self.ny)
    }

    /// Pixel range `[start, end)` of tile column `tx` (or row, with the height pitch).
    pub fn span(&self, t: usize, horizontal: bool, limit: usize) -> (usize, usize) {
        let p = if horizontal { self.tw } else { self.th };
        (
            (t as f64 * p) as usize,
            (((t + 1) as f64 * p) as usize).min(limit),
        )
    }

    /// A lower bound of `background + k × noise` over the pixels of tile `(tx, ty)`. The
    /// interpolated value is bilinear on each side of the tile's centre lines, so its least
    /// value over the tile is at a corner, the middle of an edge or the centre.
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
        // Rounding in the interpolation must not put the bound above the value.
        lo - 1e-3 * noise
    }

    /// Median background and noise over the tiles.
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

    /// Deterministic unit normals.
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
            *x += 50.0; // 10% of the samples are "stars"
        }
        let med = clipped_median(&mut v).unwrap();
        assert!(med.abs() < 0.1, "{med}");
        // Differences of pairs: 10% of them have a star in one member.
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
    fn a_gradient_is_followed_and_noise_recovered() {
        let (w, h) = (512, 384);
        let n = normals(w * h, 11);
        let img: Vec<f32> = (0..w * h)
            .map(|i| 100.0 + 0.2 * (i % w) as f32 + 3.0 * n[i])
            .collect();
        let m = Mesh::build(&img, w, h, 64, 5, None);
        // The last point is at the frame edge, beyond the outermost tile centre.
        for (x, y) in [(40.0, 40.0), (256.0, 190.0), (470.0, 340.0), (511.0, 2.0)] {
            let (b, s) = m.at(x, y);
            assert!((b - (100.0 + 0.2 * x as f32)).abs() < 1.5, "{x}: {b}");
            // The gradient across a tile (13 units, four times the noise) does not count.
            assert!((2.75..3.25).contains(&s), "{s}");
        }
        // The floor is a lower bound of the interpolated threshold.
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
        // A glow that halves every two tiles from the corner, and one tile raised by a "star".
        let (nx, ny) = (8, 6);
        let glow =
            |x: usize, y: usize| 1000.0 * 0.5f32.ppowf(((x * x + y * y) as f32).sqrt() / 2.0);
        let mut v: Vec<f32> = (0..nx * ny).map(|i| glow(i % nx, i / nx)).collect();
        v[3 * nx + 4] += 40.0;
        let out = despike(&v, nx, ny, |_| 2.0);
        for (i, got) in out.iter().enumerate() {
            let want = glow(i % nx, i / nx);
            // The raised tile comes down to the largest line mean, which lies above the glow
            // where it curves; every other tile is untouched.
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
