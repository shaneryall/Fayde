use fayde_core::Portable;

use super::{
    cubic, cutout, isolated, sample_value, square, stamp_radius, Image, PathGeometry, Region,
    Template, MAX_RADIUS, MAX_STAMPS, MIN_STAMPS,
};
use crate::stats;
use crate::trail::{axial_diff, Component};

pub(super) struct Rotation {
    centre: (f64, f64),
    rate: f64,
    width: f64,
    profile: [f64; 17],
    stacked: usize,
}

fn rotation_centre(stars: &[&Component]) -> Option<(f64, f64)> {
    let [mut xx, mut xy, mut yy, mut bx, mut by] = [0.0; 5];
    for c in stars {
        let (sin, cos) = c.theta.psin_cos();
        let projection = cos * c.x + sin * c.y;
        xx += cos * cos;
        xy += cos * sin;
        yy += sin * sin;
        bx += cos * projection;
        by += sin * projection;
    }
    let det = xx * yy - xy * xy;
    (det > 0.07 * (xx + yy) * (xx + yy))
        .then(|| ((yy * bx - xy * by) / det, (xx * by - xy * bx) / det))
}

impl Rotation {
    pub(super) fn learn(comps: &[Component], image: &Image) -> Option<Self> {
        let input = image.input;
        let mut stars: Vec<_> = comps
            .iter()
            .filter(|c| {
                c.peak >= 16.0
                    && !c.masked
                    && !c.saturated
                    && !c.edge
                    && c.elongated()
                    && c.length() >= 7.0
            })
            .collect();
        if stars.len() < 12 {
            return None;
        }
        let n = stars.len();
        let mut centre = rotation_centre(&stars)?;
        for _ in 0..2 {
            stars.retain(|c| {
                let (x, y) = (c.x - centre.0, c.y - centre.1);
                axial_diff(c.theta, x.patan2(-y)) <= 0.18
            });
            if stars.len() * 4 < 3 * n {
                return None;
            }
            centre = rotation_centre(&stars)?;
        }
        if !(0.0..input.w as f64).contains(&centre.0) || !(0.0..input.h as f64).contains(&centre.1)
        {
            return None;
        }
        let distance = |c: &Component| (c.x - centre.0).phypot(c.y - centre.1);
        let mut rates: Vec<_> = stars.iter().map(|c| c.length() / distance(c)).collect();
        let rate = stats::median(&mut rates)?;
        if rate <= 0.0 || rate > 0.15 {
            return None;
        }
        stars.retain(|c| (c.length() / distance(c) / rate - 1.0).abs() < 0.22);
        if stars.len() * 3 < n * 2 {
            return None;
        }
        stars.sort_by(|a, b| b.peak.total_cmp(&a.peak));
        let mut profiles = Vec::new();
        for c in stars.into_iter().take(4 * MAX_STAMPS) {
            let (nx, ny, len) = rotation_axes(centre, rate, c.x, c.y);
            let radius = stamp_radius(c);
            if radius > MAX_RADIUS || !isolated(c, comps, input) {
                continue;
            }
            let Some(stamp) = cutout(image, c, radius) else {
                continue;
            };
            let mut profile = [0.0; 17];
            for (j, val) in profile.iter_mut().enumerate() {
                let v = (j as f64 - 8.0) * 0.5;
                for i in -3..=3 {
                    let u = len * i as f64 / 10.0;
                    *val += sample_value(
                        &stamp.pixels,
                        2 * radius + 1,
                        stamp.origin.0 - ny * u + nx * v,
                        stamp.origin.1 + nx * u + ny * v,
                    );
                }
                *val = val.max(0.0);
            }
            let sum: f64 = profile.iter().sum();
            if sum > 0.0 {
                profile.iter_mut().for_each(|v| *v /= sum);
                profiles.push(profile);
            }
            if profiles.len() == MAX_STAMPS {
                break;
            }
        }
        if profiles.len() < MIN_STAMPS {
            return None;
        }
        let mut profile = [0.0; 17];
        for (j, p) in profile.iter_mut().enumerate() {
            *p = stats::median(
                &mut profiles
                    .iter()
                    .map(|s| 0.5 * (s[j] + s[16 - j]))
                    .collect::<Vec<_>>(),
            )?;
        }
        let power: f64 = profile.iter().sum();
        if power <= 0.0 {
            return None;
        }
        let var: f64 = profile
            .iter()
            .enumerate()
            .map(|(j, v)| {
                let x = (j as f64 - 8.0) * 0.5;
                x * x * v
            })
            .sum::<f64>()
            / power;
        let stacked = profiles.len();
        Some(Self {
            centre,
            rate,
            width: var.sqrt().clamp(0.55, 2.0),
            profile,
            stacked,
        })
    }
}

fn rotation_axes(centre: (f64, f64), rate: f64, x: f64, y: f64) -> (f64, f64, f64) {
    let (dx, dy) = (x - centre.0, y - centre.1);
    let radius = dx.phypot(dy);
    if radius <= f64::EPSILON {
        (1.0, 0.0, 0.0)
    } else {
        (dx / radius, dy / radius, rate * radius)
    }
}

impl Rotation {
    pub(super) fn template(&self, x: f64, y: f64) -> Option<Template> {
        let (nx, ny, len) = rotation_axes(self.centre, self.rate, x, y);
        let radius = (0.5 * len + 5.0 * self.width).ceil() as usize;
        if radius > MAX_RADIUS {
            return None;
        }
        let k = std::f64::consts::SQRT_2 * self.width;
        let profile = |v: f64| {
            let v = 2.0 * v.abs() + 8.0;
            let i = v.floor() as usize;
            let p = |j| self.profile.get(j).copied().unwrap_or(0.0);
            cubic([p(i - 1), p(i), p(i + 1), p(i + 2)], v - i as f64)
                .0
                .max(0.0)
        };
        let pixels = square(radius)
            .map(|(x, y)| {
                let (x, y) = (x as f64, y as f64);
                let (u, v) = (-ny * x + nx * y, nx * x + ny * y);
                let sag = if len > 1e-6 {
                    self.rate * (len / 24.0 - 0.5 * u * u / len)
                } else {
                    0.0
                };
                let across = profile(v - sag);
                if across == 0.0 || u.abs() > 0.5 * len + 4.0 * self.width {
                    return 0.0;
                }
                let along = if len > 1e-6 {
                    (((u - 0.5 * len) / k).perfc() - ((u + 0.5 * len) / k).perfc()) / len
                } else {
                    (-0.5 * u * u / (self.width * self.width)).pexp()
                };
                (across * along) as f32
            })
            .collect();
        let geometry = PathGeometry {
            length: len,
            bend: self.rate * len / 720.0_f64.sqrt(),
            width: (self.width * self.width - 1.0 / 12.0).sqrt(),
        };
        Template::new(pixels, 2 * radius + 1, self.stacked, Some(geometry))
    }

    pub(super) fn regions(&self, w: usize, h: usize) -> Option<Vec<(Region, Template)>> {
        let step = (2.5 / self.rate).round().clamp(32.0, 128.0) as usize;
        let mut regions = Vec::new();
        let (nx, ny) = (w.div_ceil(step), h.div_ceil(step));
        for j in 0..ny {
            for i in 0..nx {
                let region = Region {
                    x0: w * i / nx,
                    x1: w * (i + 1) / nx,
                    y0: h * j / ny,
                    y1: h * (j + 1) / ny,
                };
                let x = 0.5 * (region.x0 + region.x1 - 1) as f64;
                let y = 0.5 * (region.y0 + region.y1 - 1) as f64;
                regions.push((region, self.template(x, y)?));
            }
        }
        Some(regions)
    }
}
