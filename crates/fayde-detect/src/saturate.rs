#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Plateau {
    pub x: f64,
    pub y: f64,
    pub area: usize,
    pub radius: f64,
    pub linear: bool,
    pub truncated: bool,
}

pub(crate) const MAX_AREA: usize = 200_000;

pub(crate) struct Regions {
    w: usize,
    h: usize,
    label: Vec<u32>,
    pub info: Vec<Plateau>,
}

impl Regions {
    pub fn new(w: usize, h: usize) -> Self {
        Self {
            w,
            h,
            label: vec![0; w * h],
            info: Vec::new(),
        }
    }

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

fn measure(pixels: &mut [(u32, u32)], truncated: bool) -> Plateau {
    let n = pixels.len();
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
    let longest_v = (0..n)
        .filter(|&i| 2 * hrun[i] >= longest_h)
        .map(|i| vrun[i])
        .max()
        .unwrap_or(0);
    let scale = longest_h.min(longest_v);
    let (mut sx, mut sy, mut core) = (0.0f64, 0.0f64, 0usize);
    for i in 0..n {
        if 2 * hrun[i] >= scale && 2 * vrun[i] >= scale {
            sx += f64::from(pixels[i].1);
            sy += f64::from(pixels[i].0);
            core += 1;
        }
    }
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
    let (major, minor) = (mean + root, (mean - root).max(0.0) + area / 12.0);
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
        let (w, h) = (600, 500);
        let img = vec![100f32; w * h];
        let mut r = Regions::new(w, h);
        let first = r.at(&img, 10, 10, 50.0).unwrap();
        assert!(r.info[first].truncated && r.info[first].area <= MAX_AREA);
        assert_eq!(r.at(&img, 11, 9, 50.0), Some(first));
        assert_eq!(r.info.len(), 1);
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
