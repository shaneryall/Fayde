use fayde_core::Portable;

use crate::mesh::Mesh;

#[derive(Debug, Clone, Copy)]
pub(crate) struct Peak {
    pub ix: usize,
    pub iy: usize,
    pub x: f64,
    pub y: f64,
    pub height: f64,
    pub amp: f64,
    pub noise: f64,
    pub bg: f64,
    pub major2: f64,
    pub minor2: f64,
    pub theta: f64,
    pub fitted: bool,
}

impl Peak {
    pub fn snr(&self) -> f64 {
        self.height / self.noise
    }

    pub fn size2(&self) -> f64 {
        0.5 * (self.major2.min(100.0 * self.minor2) + self.minor2)
    }

    pub fn elongation(&self) -> f64 {
        (self.major2 / self.minor2).sqrt()
    }
}

pub(crate) fn find(
    s: &[f32],
    w: usize,
    h: usize,
    mesh: &Mesh,
    k: f32,
    valid: Option<&[bool]>,
) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    if w < 3 || h < 3 {
        return out;
    }
    let (nx, ny) = mesh.tiles();
    for ty in 0..ny {
        let (y0, y1) = mesh.span(ty, false, h);
        let floors: Vec<(usize, usize, f32)> = (0..nx)
            .map(|tx| {
                let (x0, x1) = mesh.span(tx, true, w);
                (x0.max(1), x1.min(w - 1), mesh.floor(tx, ty, k, w, h))
            })
            .collect();
        for y in y0.max(1)..y1.min(h - 1) {
            let row = &s[y * w..(y + 1) * w];
            for &(x0, x1, floor) in &floors {
                if x0 >= x1 {
                    continue;
                }
                let seg = &row[x0..x1];
                let above = seg.iter().fold(0u32, |a, &v| a + u32::from(v > floor));
                if above == 0 {
                    continue;
                }
                for (dx, &v) in seg.iter().enumerate() {
                    if v <= floor {
                        continue;
                    }
                    let x = x0 + dx;
                    let (up, dn) = (&s[(y - 1) * w..y * w], &s[(y + 1) * w..(y + 2) * w]);
                    let is_max = v > up[x - 1]
                        && v > up[x]
                        && v > up[x + 1]
                        && v > row[x - 1]
                        && v >= row[x + 1]
                        && v >= dn[x - 1]
                        && v >= dn[x]
                        && v >= dn[x + 1];
                    if !is_max {
                        continue;
                    }
                    let (bg, noise) = mesh.at(x as f64, y as f64);
                    if v - bg <= k * noise {
                        continue;
                    }
                    if let Some(ok) = valid {
                        let clean =
                            (y - 1..=y + 1).all(|yy| (x - 1..=x + 1).all(|xx| ok[yy * w + xx]));
                        if !clean {
                            continue;
                        }
                    }
                    out.push((x, y));
                }
            }
        }
    }
    out
}

pub(crate) fn ring_median(s: &[f32], w: usize, h: usize, x: f64, y: f64, rho: f64) -> Option<f64> {
    const DIRS: [(f64, f64); 8] = [
        (0.923_879_532_511_286_7, 0.382_683_432_365_089_8),
        (0.382_683_432_365_089_8, 0.923_879_532_511_286_7),
        (-0.382_683_432_365_089_8, 0.923_879_532_511_286_7),
        (-0.923_879_532_511_286_7, 0.382_683_432_365_089_8),
        (-0.923_879_532_511_286_7, -0.382_683_432_365_089_8),
        (-0.382_683_432_365_089_8, -0.923_879_532_511_286_7),
        (0.382_683_432_365_089_8, -0.923_879_532_511_286_7),
        (0.923_879_532_511_286_7, -0.382_683_432_365_089_8),
    ];
    let mut v = [0f64; 8];
    let mut n = 0;
    for (cs, sn) in DIRS {
        if let Some(q) = sample(s, w, h, x + rho * cs, y + rho * sn) {
            v[n] = q;
            n += 1;
        }
    }
    if n < 6 {
        return None;
    }
    let v = &mut v[..n];
    v.sort_unstable_by(f64::total_cmp);
    Some(v[(n - 1) / 2])
}

pub(crate) const RING_NOISE: f64 = 0.44;

pub(crate) fn measure(s: &[f32], w: usize, ix: usize, iy: usize, bg: f64, noise: f64) -> Peak {
    let v = |dx: isize, dy: isize| -> f64 {
        let i = (iy as isize + dy) as usize * w + (ix as isize + dx) as usize;
        f64::from(s[i]) - bg
    };
    let height = v(0, 0);
    let mut peak = Peak {
        ix,
        iy,
        x: ix as f64,
        y: iy as f64,
        height,
        amp: height,
        noise,
        bg,
        major2: f64::INFINITY,
        minor2: f64::INFINITY,
        theta: 0.0,
        fitted: false,
    };
    let mut z = [[0f64; 3]; 3];
    for (j, row) in z.iter_mut().enumerate() {
        for (i, c) in row.iter_mut().enumerate() {
            let val = v(i as isize - 1, j as isize - 1);
            if val <= 0.02 * height {
                return peak;
            }
            *c = val.pln();
        }
    }
    let col = |i: usize| z[0][i] + z[1][i] + z[2][i];
    let row = |j: usize| z[j][0] + z[j][1] + z[j][2];
    let total = col(0) + col(1) + col(2);
    let gx = (col(2) - col(0)) / 6.0;
    let gy = (row(2) - row(0)) / 6.0;
    let hxx = (col(0) + col(2) - 2.0 * col(1)) / 3.0;
    let hyy = (row(0) + row(2) - 2.0 * row(1)) / 3.0;
    let hxy = (z[2][2] - z[2][0] - z[0][2] + z[0][0]) / 4.0;
    let p0 = (5.0 * total - 3.0 * (col(0) + col(2) + row(0) + row(2))) / 9.0;

    let mean = 0.5 * (hxx + hyy);
    let diff = 0.5 * (hxx - hyy);
    let root = diff.phypot(hxy);
    let (lo, hi) = (mean - root, mean + root);
    if lo >= 0.0 {
        return peak;
    }
    peak.fitted = true;
    peak.minor2 = -1.0 / lo;
    peak.major2 = if hi < 0.0 { -1.0 / hi } else { f64::INFINITY };
    peak.theta = 0.5 * (2.0 * hxy).patan2(hxx - hyy);
    let det = hxx * hyy - hxy * hxy;
    if hi < 0.0 && det > 0.0 {
        let dx = -(hyy * gx - hxy * gy) / det;
        let dy = -(hxx * gy - hxy * gx) / det;
        if dx.abs() <= 1.0 && dy.abs() <= 1.0 {
            peak.x += dx;
            peak.y += dy;
            let top = p0 + 0.5 * (gx * dx + gy * dy);
            peak.amp = top.pexp().max(height);
        }
    } else {
        let (sn, cs) = peak.theta.psin_cos();
        let g_across = -gx * sn + gy * cs;
        let d = (-g_across / lo).clamp(-1.0, 1.0);
        peak.x += -sn * d;
        peak.y += cs * d;
    }
    peak
}

pub(crate) fn sample(s: &[f32], w: usize, h: usize, x: f64, y: f64) -> Option<f64> {
    if !(x >= 0.0 && y >= 0.0 && x <= (w - 1) as f64 && y <= (h - 1) as f64) {
        return None;
    }
    let (x0, y0) = ((x as usize).min(w - 2), (y as usize).min(h - 2));
    let (fx, fy) = (x - x0 as f64, y - y0 as f64);
    let p = |xx: usize, yy: usize| f64::from(s[yy * w + xx]);
    let a = p(x0, y0) * (1.0 - fx) + p(x0 + 1, y0) * fx;
    let b = p(x0, y0 + 1) * (1.0 - fx) + p(x0 + 1, y0 + 1) * fx;
    Some(a * (1.0 - fy) + b * fy)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn gaussian(w: usize, h: usize, x0: f64, y0: f64, a2: f64, b2: f64, th: f64) -> Vec<f32> {
        let (sn, cs) = th.psin_cos();
        (0..w * h)
            .map(|i| {
                let (dx, dy) = ((i % w) as f64 - x0, (i / w) as f64 - y0);
                let (u, v) = (dx * cs + dy * sn, -dx * sn + dy * cs);
                (10.0 + 500.0 * (-0.5 * (u * u / a2 + v * v / b2)).pexp()) as f32
            })
            .collect()
    }

    fn flat_mesh(w: usize, h: usize) -> Mesh {
        let mut corner: Vec<f32> = vec![10.0; w * h];
        for (i, v) in corner.iter_mut().enumerate() {
            *v += ((i * 7919) % 1000) as f32 / 1000.0 - 0.5;
        }
        Mesh::build(&corner, w, h, w.max(h), 3, None)
    }

    #[test]
    fn the_fit_is_exact_for_an_elliptical_gaussian() {
        let (w, h) = (41, 37);
        for &(x0, y0, a2, b2, th) in &[
            (20.3, 18.6, 2.5, 2.5, 0.0),
            (19.5, 17.5, 6.0, 2.0, 0.6),
            (20.0, 18.0, 4.0, 1.5, -1.2),
            (20.45, 18.49, 9.0, 3.0, 1.5),
        ] {
            let img = gaussian(w, h, x0, y0, a2, b2, th);
            let mesh = flat_mesh(w, h);
            let found = find(&img, w, h, &mesh, 5.0, None);
            assert_eq!(found.len(), 1, "{found:?}");
            let (b, n) = mesh.at(found[0].0 as f64, found[0].1 as f64);
            let p = measure(&img, w, found[0].0, found[0].1, f64::from(b), f64::from(n));
            assert!(p.fitted);
            let (bg, _) = mesh.at(20.0, 18.0);
            let tol = 0.02 + f64::from((bg - 10.0).abs()) * 0.01;
            assert!((p.x - x0).abs() < tol && (p.y - y0).abs() < tol, "{p:?}");
            assert!((p.major2 / a2 - 1.0).abs() < 0.03, "{} vs {a2}", p.major2);
            assert!((p.minor2 / b2 - 1.0).abs() < 0.03, "{} vs {b2}", p.minor2);
            assert!((p.amp / 500.0 - 1.0).abs() < 0.01, "{}", p.amp);
            if a2 != b2 {
                let d = (p.theta - th).psin().abs();
                assert!(d < 0.02, "theta {} vs {th}", p.theta);
            }
        }
    }

    #[test]
    fn a_plateau_yields_one_maximum_and_borders_none() {
        let (w, h) = (20, 20);
        let mut img = vec![10f32; w * h];
        for y in 8..11 {
            for x in 8..11 {
                img[y * w + x] = 100.0;
            }
        }
        img[0] = 500.0;
        let mesh = flat_mesh(w, h);
        assert_eq!(find(&img, w, h, &mesh, 5.0, None), vec![(8, 8)]);
        let mut valid = vec![true; w * h];
        valid[7 * w + 7] = false;
        assert!(find(&img, w, h, &mesh, 5.0, Some(&valid)).is_empty());
    }

    #[test]
    fn a_ridge_has_an_infinite_major_axis_along_its_direction() {
        let (w, h) = (40, 40);
        let th = 30f64.to_radians();
        let (sn, cs) = th.psin_cos();
        let mut img: Vec<f32> = (0..w * h)
            .map(|i| {
                let (dx, dy) = ((i % w) as f64 - 20.0, (i / w) as f64 - 20.0);
                let c = -dx * sn + dy * cs;
                (10.0 + 200.0 * (-0.5 * c * c / 2.0).pexp()) as f32
            })
            .collect();
        img[20 * w + 20] += 0.5;
        let mesh = flat_mesh(w, h);
        let (b, n) = mesh.at(20.0, 20.0);
        let p = measure(&img, w, 20, 20, f64::from(b), f64::from(n));
        assert!(p.fitted);
        assert!(p.elongation() > 5.0, "{p:?}");
        assert!((p.minor2 / 2.0 - 1.0).abs() < 0.1, "{}", p.minor2);
        assert!((p.theta - th).psin().abs() < 0.05, "{}", p.theta);
    }

    #[test]
    fn the_ring_median_ignores_a_gradient_a_neighbour_and_a_streak() {
        let (w, h) = (60, 60);
        let img: Vec<f32> = (0..w * h)
            .map(|i| {
                let (x, y) = ((i % w) as f64, (i / w) as f64);
                let (dx, dy) = (x - 35.0, y - 31.0);
                let star = 900.0 * (-0.25 * (dx * dx + dy * dy)).pexp();
                let across = y - 30.0 - 0.3 * (x - 30.0);
                let line = 300.0 * (-0.5 * across * across).pexp();
                (50.0 + 0.5 * x - 0.3 * y + star + line) as f32
            })
            .collect();
        let m = ring_median(&img, w, h, 30.0, 30.0, 5.0).unwrap();
        assert!((m - 56.0).abs() < 2.5, "{m}");
        assert!(ring_median(&img, w, h, 1.0, 1.0, 5.0).is_none());
    }

    #[test]
    fn bilinear_sampling() {
        let img = [0f32, 1.0, 2.0, 3.0];
        assert_eq!(sample(&img, 2, 2, 0.5, 0.5), Some(1.5));
        assert_eq!(sample(&img, 2, 2, 1.0, 1.0), Some(3.0));
        assert_eq!(sample(&img, 2, 2, 1.1, 0.0), None);
    }
}
