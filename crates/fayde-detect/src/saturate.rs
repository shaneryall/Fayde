//! Saturated stars: the centre of a clipped top.
//!
//! A saturated star has a flat top with no maximum to fit, and often a bleed trail along one
//! detector axis. Its centre is taken from the saturated region itself: the mean position of
//! the region's *core*, the pixels whose runs of saturated pixels, along the row and along the
//! column, are each at least half the longest such run in the region. For a disc that is a
//! symmetric part of the disc; a bleed trail, a few pixels wide, has short runs across it and
//! contributes only where it crosses the disc. The method is Fayde's own.
//!
//! The work is bounded by the frame, whatever the image and the saturation level: every
//! saturated pixel is flooded at most once ([`Regions`] remembers which region owns it), and a
//! region's measurements cost its pixel count, not its bounding box.

/// A connected region of saturated pixels.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Plateau {
    /// Centre of the region's core.
    pub x: f64,
    /// Centre of the region's core.
    pub y: f64,
    /// Saturated pixels in the region.
    pub area: usize,
    /// Radius of a disc with the area of the core, pixels.
    pub radius: f64,
    /// The region is a line of even width (a saturated streak), not a star: over three times
    /// as long as wide, with no part wider than the rest.
    pub linear: bool,
    /// The flood stopped at [`MAX_AREA`] pixels: part of a larger region (a frame saturated
    /// over a wide area), not a star.
    pub truncated: bool,
}

/// A region is followed for at most this many pixels.
pub(crate) const MAX_AREA: usize = 200_000;

/// The saturated regions found so far in a frame, and the region each pixel belongs to.
pub(crate) struct Regions {
    w: usize,
    h: usize,
    /// 0: not in a region; `k + 1`: in `info[k]`.
    label: Vec<u32>,
    /// One entry per region, in the order found.
    pub info: Vec<Plateau>,
}

impl Regions {
    /// No region yet, for a frame of `w × h` pixels.
    pub fn new(w: usize, h: usize) -> Self {
        Self {
            w,
            h,
            label: vec![0; w * h],
            info: Vec::new(),
        }
    }

    /// The saturated region (pixels ≥ `level`, 8-connected) touching the 5 × 5 pixels around
    /// `(ix, iy)`: its index in `info`. A region not seen before is flooded now.
    pub fn at(&mut self, raw: &[f32], ix: usize, iy: usize, level: f32) -> Option<usize> {
        let (w, h) = (self.w, self.h);
        let mut seed = None;
        for y in iy.saturating_sub(2)..(iy + 3).min(h) {
            for x in ix.saturating_sub(2)..(ix + 3).min(w) {
                let l = self.label[y * w + x];
                if l != 0 {
                    return Some(l as usize - 1);
                }
                if seed.is_none() && raw[y * w + x] >= level {
                    seed = Some((x, y));
                }
            }
        }
        let seed = seed?;
        let id = u32::try_from(self.info.len() + 1).ok()?;
        // Flood fill. A pixel is labelled when it is pushed, so it is never pushed twice, in
        // this region or a later one.
        let mut pixels: Vec<(u32, u32)> = Vec::new();
        let mut stack = vec![seed];
        self.label[seed.1 * w + seed.0] = id;
        let mut truncated = false;
        while let Some((x, y)) = stack.pop() {
            pixels.push((y as u32, x as u32));
            if pixels.len() + stack.len() >= MAX_AREA {
                truncated = true;
                continue;
            }
            for yy in y.saturating_sub(1)..(y + 2).min(h) {
                for xx in x.saturating_sub(1)..(x + 2).min(w) {
                    let i = yy * w + xx;
                    if self.label[i] == 0 && raw[i] >= level {
                        self.label[i] = id;
                        stack.push((xx, yy));
                    }
                }
            }
        }
        self.info.push(measure(&mut pixels, truncated));
        Some(id as usize - 1)
    }

    /// A region with a pixel within `pad` pixels of `(x, y)` (in either axis), preferring one
    /// that `keep` accepts.
    pub fn near(&self, x: f64, y: f64, pad: usize, keep: impl Fn(usize) -> bool) -> Option<usize> {
        let (cx, cy) = (x.round().max(0.0) as usize, y.round().max(0.0) as usize);
        let mut other = None;
        for yy in cy.saturating_sub(pad)..(cy + pad + 1).min(self.h) {
            for xx in cx.saturating_sub(pad)..(cx + pad + 1).min(self.w) {
                let l = self.label[yy * self.w + xx];
                if l != 0 {
                    let k = l as usize - 1;
                    if keep(k) {
                        return Some(k);
                    }
                    other = Some(k);
                }
            }
        }
        other
    }
}

/// Centre, size and shape of a region from its pixels `(y, x)`. Reorders `pixels`.
fn measure(pixels: &mut [(u32, u32)], truncated: bool) -> Plateau {
    let n = pixels.len();
    // Runs along rows: with the pixels in (y, x) order a run is a stretch of consecutive x.
    pixels.sort_unstable();
    let mut hrun = vec![0u32; n];
    let mut start = 0;
    for i in 1..=n {
        let breaks = i == n || pixels[i].0 != pixels[i - 1].0 || pixels[i].1 != pixels[i - 1].1 + 1;
        if breaks {
            for v in &mut hrun[start..i] {
                *v = (i - start) as u32;
            }
            start = i;
        }
    }
    // Runs along columns: the same in (x, y) order, through an index.
    let mut by_col: Vec<u32> = (0..n as u32).collect();
    by_col.sort_unstable_by_key(|&i| (pixels[i as usize].1, pixels[i as usize].0));
    let mut vrun = vec![0u32; n];
    start = 0;
    for i in 1..=n {
        let breaks = i == n || {
            let (a, b) = (pixels[by_col[i] as usize], pixels[by_col[i - 1] as usize]);
            a.1 != b.1 || a.0 != b.0 + 1
        };
        if breaks {
            for &j in &by_col[start..i] {
                vrun[j as usize] = (i - start) as u32;
            }
            start = i;
        }
    }
    let longest_h = hrun.iter().copied().max().unwrap_or(0);
    // The longest column run among pixels of long rows: a bleed trail's own length (rows of
    // which are short) does not set the scale.
    let longest_v = (0..n)
        .filter(|&i| 2 * hrun[i] >= longest_h)
        .map(|i| vrun[i])
        .max()
        .unwrap_or(0);
    // A bleed trail is long in one axis: the disc's diameter is the shorter of the two.
    let scale = longest_h.min(longest_v);
    let (mut sx, mut sy, mut core) = (0.0f64, 0.0f64, 0usize);
    for i in 0..n {
        if 2 * hrun[i] >= scale && 2 * vrun[i] >= scale {
            sx += f64::from(pixels[i].1);
            sy += f64::from(pixels[i].0);
            core += 1;
        }
    }
    // Axis ratio from the second moments of the whole region. A star with a bleed trail is
    // long too, but most of its pixels are not core.
    let area = n as f64;
    let (mx, my) = (
        pixels.iter().map(|p| f64::from(p.1)).sum::<f64>() / area,
        pixels.iter().map(|p| f64::from(p.0)).sum::<f64>() / area,
    );
    let (mut sxx, mut sxy, mut syy) = (0.0f64, 0.0f64, 0.0f64);
    for p in pixels.iter() {
        let (dx, dy) = (f64::from(p.1) - mx, f64::from(p.0) - my);
        sxx += dx * dx;
        sxy += dx * dy;
        syy += dy * dy;
    }
    let mean = 0.5 * (sxx + syy);
    let root = (0.25 * (sxx - syy) * (sxx - syy) + sxy * sxy).sqrt();
    // Variances of a pixel region include 1/12 px² per axis from the pixels themselves.
    let (major, minor) = (mean + root, (mean - root).max(0.0) + area / 12.0);
    // The core is never empty (the pixel of the longest runs is in it); the mean of the whole
    // region stands in if it were.
    let (cx, cy) = if core > 0 {
        (sx / core as f64, sy / core as f64)
    } else {
        (mx, my)
    };
    let core = core.max(1);
    Plateau {
        x: cx,
        y: cy,
        area: n,
        radius: (core as f64 / std::f64::consts::PI).sqrt(),
        linear: major > 9.0 * minor && core as f64 > 0.6 * area,
        truncated,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn disc(w: usize, h: usize, cx: f64, cy: f64, r: f64) -> Vec<f32> {
        (0..w * h)
            .map(|i| {
                let (dx, dy) = ((i % w) as f64 - cx, (i / w) as f64 - cy);
                if dx * dx + dy * dy <= r * r {
                    100.0
                } else {
                    1.0
                }
            })
            .collect()
    }

    fn plateau(img: &[f32], w: usize, h: usize, ix: usize, iy: usize) -> Option<Plateau> {
        let mut r = Regions::new(w, h);
        r.at(img, ix, iy, 50.0).map(|k| r.info[k])
    }

    #[test]
    fn the_centre_of_a_disc_is_found() {
        let img = disc(60, 50, 30.3, 24.6, 6.0);
        let p = plateau(&img, 60, 50, 28, 23).unwrap();
        assert!(
            (p.x - 30.3).abs() < 0.4 && (p.y - 24.6).abs() < 0.4,
            "{p:?}"
        );
        assert!(p.area > 100 && p.area < 125);
        assert!(!p.linear && !p.truncated);
        assert!(plateau(&img, 60, 50, 5, 5).is_none());
    }

    #[test]
    fn a_saturated_streak_is_not_a_star() {
        let (w, h) = (80, 80);
        // A diagonal band 4 pixels wide.
        let img: Vec<f32> = (0..w * h)
            .map(|i| {
                let (x, y) = ((i % w) as f64, (i / w) as f64);
                let across = (y - 0.5 * x - 20.0) / 1.118;
                if across.abs() <= 2.0 && (10.0..70.0).contains(&x) {
                    100.0
                } else {
                    1.0
                }
            })
            .collect();
        let p = plateau(&img, w, h, 40, 40).unwrap();
        assert!(p.linear, "{p:?}");
    }

    #[test]
    fn a_bleed_trail_does_not_move_the_centre() {
        let (w, h) = (60, 120);
        let mut img = disc(w, h, 30.0, 80.0, 6.0);
        // A trail 3 pixels wide, running 70 pixels up from the star and 10 down.
        for y in 10..=90 {
            for x in 29..=31 {
                img[y * w + x] = 100.0;
            }
        }
        let p = plateau(&img, w, h, 30, 80).unwrap();
        assert!(
            (p.x - 30.0).abs() < 0.5 && (p.y - 80.0).abs() < 0.8,
            "{p:?}"
        );
        assert!(!p.linear);
    }

    #[test]
    fn every_pixel_is_flooded_once_and_a_huge_region_is_cut_short() {
        // A frame saturated all over: larger than one region may be.
        let (w, h) = (600, 500);
        let img = vec![100f32; w * h];
        let mut r = Regions::new(w, h);
        let first = r.at(&img, 10, 10, 50.0).unwrap();
        assert!(r.info[first].truncated && r.info[first].area <= MAX_AREA);
        // The same place again is the same region, with no new flood.
        assert_eq!(r.at(&img, 11, 9, 50.0), Some(first));
        assert_eq!(r.info.len(), 1);
        // Asking everywhere labels each pixel at most once: the regions' areas add up to no
        // more than the frame.
        for y in (0..h).step_by(7) {
            for x in (0..w).step_by(7) {
                r.at(&img, x, y, 50.0);
            }
        }
        let total: usize = r.info.iter().map(|p| p.area).sum();
        assert!(total <= w * h, "{total}");
        assert!(r.near(300.0, 250.0, 3, |_| true).is_some());
    }
}
