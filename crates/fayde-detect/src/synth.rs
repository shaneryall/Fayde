use fayde_core::Portable;

#[derive(Debug, Clone)]
pub struct Rng(pub u64);

impl Rng {
    pub fn uniform(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64 + 0.5) / (1u64 << 53) as f64
    }

    pub fn normal(&mut self) -> f64 {
        let (u, v) = (self.uniform(), self.uniform());
        (-2.0 * u.pln()).sqrt() * (2.0 * std::f64::consts::PI * v).pcos()
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Star {
    pub x: f64,
    pub y: f64,
    pub flux: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Streak {
    pub x0: f64,
    pub y0: f64,
    pub x1: f64,
    pub y1: f64,
    pub flux_per_px: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Trail {
    pub length: f64,
    pub angle: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Stray {
    pub x: f64,
    pub y: f64,
    pub peak: f64,
    pub half_px: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Cosmic {
    pub x: f64,
    pub y: f64,
    pub length: f64,
    pub angle: f64,
    pub adu_per_px: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Psf {
    Gaussian(f64),
    Moffat { fwhm: f64, beta: f64 },
}

#[derive(Debug, Clone, PartialEq)]
pub struct Scene {
    pub width: usize,
    pub height: usize,
    pub sky: f64,
    pub gradient: (f64, f64),
    pub read_noise: f64,
    pub gain: f64,
    pub saturation: Option<f64>,
    pub psf: Psf,
    pub stars: Vec<Star>,
    pub trail: Option<Trail>,
    pub movers: Vec<Star>,
    pub streaks: Vec<Streak>,
    pub hot: Vec<(usize, usize, f64)>,
    pub stray: Option<Stray>,
    pub cosmics: Vec<Cosmic>,
    pub quantise: bool,
}

impl Scene {
    pub fn new(width: usize, height: usize) -> Self {
        Self {
            width,
            height,
            sky: 1000.0,
            gradient: (0.0, 0.0),
            read_noise: 5.0,
            gain: 1.0,
            saturation: None,
            psf: Psf::Gaussian(1.5),
            stars: Vec::new(),
            trail: None,
            movers: Vec::new(),
            streaks: Vec::new(),
            hot: Vec::new(),
            stray: None,
            cosmics: Vec::new(),
            quantise: false,
        }
    }

    pub fn psf_sigma(&self) -> f64 {
        match self.psf {
            Psf::Gaussian(s) => s,
            Psf::Moffat { fwhm, .. } => fwhm / 2.354_820_045,
        }
    }

    pub fn sky_noise(&self) -> f64 {
        (self.sky / self.gain + self.read_noise * self.read_noise).sqrt()
    }

    pub fn model(&self) -> Vec<f64> {
        let (w, h) = (self.width, self.height);
        let (cx, cy) = (0.5 * (w as f64 - 1.0), 0.5 * (h as f64 - 1.0));
        let mut img: Vec<f64> = (0..w * h)
            .map(|i| {
                let (x, y) = ((i % w) as f64, (i / w) as f64);
                let glow = self.stray.map_or(0.0, |s| {
                    let d = (x - s.x).phypot(y - s.y);
                    s.peak * (-std::f64::consts::LN_2 * d / s.half_px).pexp()
                });
                self.sky + self.gradient.0 * (x - cx) + self.gradient.1 * (y - cy) + glow
            })
            .collect();
        let floor = 0.02 * self.sky_noise();
        let sigma = self.psf_sigma();
        let trailed = self.trail.filter(|t| t.length > 0.0);
        for s in &self.stars {
            if let Some(t) = trailed {
                let (sn, cs) = t.angle.psin_cos();
                let (hx, hy) = (0.5 * t.length * cs, 0.5 * t.length * sn);
                let seg = Streak {
                    x0: s.x - hx,
                    y0: s.y - hy,
                    x1: s.x + hx,
                    y1: s.y + hy,
                    flux_per_px: s.flux / t.length,
                };
                add_streak(&mut img, w, h, &seg, sigma);
                continue;
            }
            match self.psf {
                Psf::Gaussian(sigma) => add_gaussian(&mut img, w, h, s, sigma, floor),
                Psf::Moffat { fwhm, beta } => add_moffat(&mut img, w, h, s, fwhm, beta, floor),
            }
        }
        for s in &self.movers {
            match self.psf {
                Psf::Gaussian(sigma) => add_gaussian(&mut img, w, h, s, sigma, floor),
                Psf::Moffat { fwhm, beta } => add_moffat(&mut img, w, h, s, fwhm, beta, floor),
            }
        }
        for t in &self.streaks {
            add_streak(&mut img, w, h, t, sigma);
        }
        img
    }

    pub fn render(&self, seed: u64) -> Vec<f32> {
        let mut rng = Rng(seed);
        let mut img = self.model();
        let rn2 = self.read_noise * self.read_noise;
        for v in &mut img {
            let var = v.max(0.0) / self.gain + rn2;
            *v += var.sqrt() * rng.normal();
        }
        for &(x, y, adu) in &self.hot {
            if x < self.width && y < self.height {
                img[y * self.width + x] += adu;
            }
        }
        for c in &self.cosmics {
            let (sn, cs) = c.angle.psin_cos();
            let steps = (c.length * 10.0).ceil().max(1.0) as usize;
            let each = c.adu_per_px * c.length.max(1.0) / steps as f64;
            for k in 0..steps {
                let t = c.length * (k as f64 + 0.5) / steps as f64;
                let (x, y) = ((c.x + t * cs).round(), (c.y + t * sn).round());
                if x >= 0.0 && y >= 0.0 && (x as usize) < self.width && (y as usize) < self.height {
                    img[y as usize * self.width + x as usize] += each;
                }
            }
        }
        img.iter()
            .map(|&v| {
                let v = self.saturation.map_or(v, |s| v.min(s));
                if self.quantise {
                    v.round() as f32
                } else {
                    v as f32
                }
            })
            .collect()
    }
}

fn erf(x: f64) -> f64 {
    1.0 - x.perfc()
}

fn pixel_fraction(i: isize, c: f64, sigma: f64) -> f64 {
    let k = 1.0 / (std::f64::consts::SQRT_2 * sigma);
    0.5 * (erf((i as f64 + 0.5 - c) * k) - erf((i as f64 - 0.5 - c) * k))
}

fn add_gaussian(img: &mut [f64], w: usize, h: usize, s: &Star, sigma: f64, floor: f64) {
    let peak = s.flux / (2.0 * std::f64::consts::PI * sigma * sigma);
    let r = if peak > floor {
        sigma * (2.0 * (peak / floor).pln()).sqrt()
    } else {
        0.0
    };
    let r = r.max(3.0 * sigma).ceil() as isize + 1;
    let (ix, iy) = (s.x.round() as isize, s.y.round() as isize);
    let fx: Vec<f64> = (ix - r..=ix + r)
        .map(|i| pixel_fraction(i, s.x, sigma))
        .collect();
    for y in (iy - r).max(0)..=(iy + r).min(h as isize - 1) {
        let fy = s.flux * pixel_fraction(y, s.y, sigma);
        for x in (ix - r).max(0)..=(ix + r).min(w as isize - 1) {
            img[y as usize * w + x as usize] += fy * fx[(x - ix + r) as usize];
        }
    }
}

fn add_moffat(img: &mut [f64], w: usize, h: usize, s: &Star, fwhm: f64, beta: f64, floor: f64) {
    let alpha = fwhm / (2.0 * (2f64.ppowf(1.0 / beta) - 1.0).sqrt());
    let peak = s.flux * (beta - 1.0) / (std::f64::consts::PI * alpha * alpha);
    let r = if peak > floor {
        alpha * ((peak / floor).ppowf(1.0 / beta) - 1.0).max(0.0).sqrt()
    } else {
        0.0
    };
    let r = r.clamp(3.0 * fwhm, 60.0 * fwhm).ceil() as isize + 1;
    let (ix, iy) = (s.x.round() as isize, s.y.round() as isize);
    const SUB: usize = 5;
    for y in (iy - r).max(0)..=(iy + r).min(h as isize - 1) {
        for x in (ix - r).max(0)..=(ix + r).min(w as isize - 1) {
            let mut acc = 0.0;
            for sy in 0..SUB {
                for sx in 0..SUB {
                    let dx = x as f64 + (sx as f64 + 0.5) / SUB as f64 - 0.5 - s.x;
                    let dy = y as f64 + (sy as f64 + 0.5) / SUB as f64 - 0.5 - s.y;
                    acc += (1.0 + (dx * dx + dy * dy) / (alpha * alpha)).ppowf(-beta);
                }
            }
            img[y as usize * w + x as usize] += peak * acc / (SUB * SUB) as f64;
        }
    }
}

fn add_streak(img: &mut [f64], w: usize, h: usize, t: &Streak, sigma: f64) {
    let (dx, dy) = (t.x1 - t.x0, t.y1 - t.y0);
    let len = dx.phypot(dy);
    if len <= 0.0 {
        return;
    }
    let (ux, uy) = (dx / len, dy / len);
    let s2 = sigma * sigma + 1.0 / 12.0;
    let s = s2.sqrt();
    let pad = 5.0 * s + 1.0;
    let x_lo = (t.x0.min(t.x1) - pad).floor().max(0.0) as usize;
    let x_hi = (t.x0.max(t.x1) + pad).ceil().min(w as f64 - 1.0);
    let y_lo = (t.y0.min(t.y1) - pad).floor().max(0.0) as usize;
    let y_hi = (t.y0.max(t.y1) + pad).ceil().min(h as f64 - 1.0);
    if x_hi < 0.0 || y_hi < 0.0 {
        return;
    }
    let k = 1.0 / (std::f64::consts::SQRT_2 * s);
    let norm = t.flux_per_px / ((2.0 * std::f64::consts::PI).sqrt() * s);
    for y in y_lo..=y_hi as usize {
        for x in x_lo..=x_hi as usize {
            let (px, py) = (x as f64 - t.x0, y as f64 - t.y0);
            let along = px * ux + py * uy;
            let across = -px * uy + py * ux;
            if across.abs() > pad {
                continue;
            }
            let ends = 0.5 * (erf(along * k) - erf((along - len) * k));
            if ends <= 0.0 {
                continue;
            }
            img[y * w + x] += norm * (-0.5 * across * across / s2).pexp() * ends;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn flux_is_conserved() {
        let mut sc = Scene::new(64, 64);
        sc.sky = 0.0;
        sc.stars.push(Star {
            x: 30.3,
            y: 31.7,
            flux: 1e4,
        });
        let sum: f64 = sc.model().iter().sum();
        assert!((sum / 1e4 - 1.0).abs() < 1e-3, "{sum}");
        sc.psf = Psf::Moffat {
            fwhm: 3.0,
            beta: 3.0,
        };
        let sum: f64 = sc.model().iter().sum();
        assert!((sum / 1e4 - 1.0).abs() < 0.02, "{sum}");
        sc.stars.clear();
        sc.streaks.push(Streak {
            x0: 10.0,
            y0: 12.0,
            x1: 50.0,
            y1: 42.0,
            flux_per_px: 100.0,
        });
        let sum: f64 = sc.model().iter().sum();
        assert!((sum / 5000.0 - 1.0).abs() < 1e-3, "{sum}");
        sc.streaks.clear();
        sc.psf = Psf::Gaussian(1.3);
        sc.trail = Some(Trail {
            length: 20.0,
            angle: 0.7,
        });
        sc.stars.push(Star {
            x: 31.2,
            y: 30.6,
            flux: 1e4,
        });
        let model = sc.model();
        let sum: f64 = model.iter().sum();
        assert!((sum / 1e4 - 1.0).abs() < 1e-3, "{sum}");
        let (mut mx, mut my) = (0.0, 0.0);
        for (i, v) in model.iter().enumerate() {
            mx += v * (i % 64) as f64;
            my += v * (i / 64) as f64;
        }
        assert!((mx / sum - 31.2).abs() < 0.01 && (my / sum - 30.6).abs() < 0.01);
    }

    #[test]
    fn noise_matches_the_model() {
        let sc = Scene::new(128, 128);
        let img = sc.render(7);
        let n = img.len() as f64;
        let mean = img.iter().map(|v| f64::from(*v)).sum::<f64>() / n;
        let var = img
            .iter()
            .map(|v| (f64::from(*v) - mean) * (f64::from(*v) - mean))
            .sum::<f64>()
            / n;
        assert!((mean - 1000.0).abs() < 1.0, "{mean}");
        assert!((var.sqrt() / sc.sky_noise() - 1.0).abs() < 0.03);
        assert_eq!(img, sc.render(7));
    }
}
