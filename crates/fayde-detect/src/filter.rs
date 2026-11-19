use fayde_core::Portable;

#[derive(Debug, Clone)]
pub(crate) struct Kernel {
    pub taps: Vec<f32>,
    pub radius: usize,
    pub sigma: f64,
}

impl Kernel {
    pub fn gaussian(sigma: f64) -> Self {
        let radius = (3.0 * sigma).ceil().max(1.0) as usize;
        let raw: Vec<f64> = (0..=2 * radius)
            .map(|i| {
                let d = i as f64 - radius as f64;
                (-(d * d) / (2.0 * sigma * sigma)).pexp()
            })
            .collect();
        let sum: f64 = raw.iter().sum();
        Self {
            taps: raw.iter().map(|v| (v / sum) as f32).collect(),
            radius,
            sigma,
        }
    }

    pub fn noise_gain2(&self) -> f64 {
        let g: f64 = self
            .taps
            .iter()
            .map(|t| f64::from(*t) * f64::from(*t))
            .sum();
        g * g
    }

    pub fn noise_gain(&self) -> f64 {
        self.noise_gain2().sqrt()
    }
}

#[inline]
fn reflect(i: isize, n: usize) -> usize {
    if n <= 1 {
        return 0;
    }
    let n = n as isize;
    let mut i = i;
    loop {
        if i < 0 {
            i = -i;
        } else if i >= n {
            i = 2 * (n - 1) - i;
        } else {
            return i as usize;
        }
    }
}

pub(crate) fn smooth(
    src: &[f32],
    w: usize,
    h: usize,
    k: &Kernel,
    tmp: &mut Vec<f32>,
    dst: &mut Vec<f32>,
) {
    let r = k.radius;
    tmp.clear();
    tmp.resize(w * h, 0.0);
    dst.clear();
    dst.resize(w * h, 0.0);

    for y in 0..h {
        let row = &src[y * w..(y + 1) * w];
        let out = &mut tmp[y * w..(y + 1) * w];
        if w > 2 * r {
            let inner = &mut out[r..w - r];
            for (j, &kj) in k.taps.iter().enumerate() {
                for (o, &v) in inner.iter_mut().zip(&row[j..j + w - 2 * r]) {
                    *o += kj * v;
                }
            }
        }
        let edge = r.min(w);
        for x in (0..edge).chain(w.saturating_sub(r).max(edge)..w) {
            let mut acc = 0f32;
            for (j, &kj) in k.taps.iter().enumerate() {
                acc += kj * row[reflect(x as isize + j as isize - r as isize, w)];
            }
            out[x] = acc;
        }
    }
    for y in 0..h {
        let out = &mut dst[y * w..(y + 1) * w];
        for (j, &kj) in k.taps.iter().enumerate() {
            let yy = reflect(y as isize + j as isize - r as isize, h);
            for (o, &v) in out.iter_mut().zip(&tmp[yy * w..(yy + 1) * w]) {
                *o += kj * v;
            }
        }
    }
}

pub(crate) fn bin(src: &[f32], w: usize, h: usize, b: usize) -> (Vec<f32>, usize, usize) {
    let (bw, bh) = (w / b, h / b);
    let mut out = vec![0f32; bw * bh];
    let scale = 1.0 / (b * b) as f32;
    for by in 0..bh {
        let o = &mut out[by * bw..(by + 1) * bw];
        for dy in 0..b {
            let row = &src[(by * b + dy) * w..(by * b + dy) * w + bw * b];
            for (v, c) in o.iter_mut().zip(row.chunks_exact(b)) {
                let mut s = 0f32;
                for &p in c {
                    s += p;
                }
                *v += s;
            }
        }
        for v in o.iter_mut() {
            *v *= scale;
        }
    }
    (out, bw, bh)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kernel_is_normalised_and_symmetric() {
        for sigma in [0.6, 1.0, 1.2, 2.4] {
            let k = Kernel::gaussian(sigma);
            assert_eq!(k.taps.len(), 2 * k.radius + 1);
            let sum: f32 = k.taps.iter().sum();
            assert!((sum - 1.0).abs() < 1e-6, "{sigma}: {sum}");
            for i in 0..k.radius {
                assert_eq!(k.taps[i], k.taps[2 * k.radius - i]);
            }
        }
    }

    #[test]
    fn a_constant_image_is_unchanged_and_a_point_becomes_the_kernel() {
        let (w, h) = (23, 17);
        let k = Kernel::gaussian(1.2);
        let (mut tmp, mut out) = (Vec::new(), Vec::new());
        smooth(&vec![7.0; w * h], w, h, &k, &mut tmp, &mut out);
        assert!(out.iter().all(|v| (v - 7.0).abs() < 1e-5));

        let mut img = vec![0f32; w * h];
        img[8 * w + 11] = 1.0;
        smooth(&img, w, h, &k, &mut tmp, &mut out);
        let r = k.radius;
        for dy in 0..=2 * r {
            for dx in 0..=2 * r {
                let want = k.taps[dx] * k.taps[dy];
                let got = out[(8 + dy - r) * w + 11 + dx - r];
                assert!((got - want).abs() < 1e-7);
            }
        }
        let sum: f32 = out.iter().sum();
        assert!((sum - 1.0).abs() < 1e-5);
    }

    #[test]
    fn images_narrower_than_the_kernel_are_handled() {
        let k = Kernel::gaussian(2.0);
        let (mut tmp, mut out) = (Vec::new(), Vec::new());
        smooth(
            &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            3,
            2,
            &k,
            &mut tmp,
            &mut out,
        );
        assert!(out.iter().all(|v| v.is_finite() && *v > 1.0 && *v < 6.0));
    }

    #[test]
    fn a_single_row_or_column_is_filtered_without_hanging() {
        let k = Kernel::gaussian(1.2);
        let (mut tmp, mut out) = (Vec::new(), Vec::new());
        smooth(&[1.0, 2.0, 3.0], 1, 3, &k, &mut tmp, &mut out);
        assert!(out.iter().all(|v| v.is_finite()));
        smooth(&[1.0, 2.0, 3.0], 3, 1, &k, &mut tmp, &mut out);
        assert!(out.iter().all(|v| v.is_finite()));
        smooth(&[4.0], 1, 1, &k, &mut tmp, &mut out);
        assert_eq!(out, [4.0]);
        assert_eq!(reflect(5, 1), 0);
    }

    #[test]
    fn binning_averages_blocks() {
        let src: Vec<f32> = (0..30).map(|v| v as f32).collect();
        let (out, w, h) = bin(&src, 6, 5, 2);
        assert_eq!((w, h), (3, 2));
        assert_eq!(out[0], 3.5);
        assert_eq!(out[5], (16.0 + 17.0 + 22.0 + 23.0) / 4.0);
    }
}
