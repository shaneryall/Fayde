//! Normal tail probabilities and the false-alarm threshold.

use fayde_core::Portable;

/// Upper tail of the standard normal distribution, `P(Z > z)`.
pub(crate) fn normal_tail(z: f64) -> f64 {
    0.5 * (z / std::f64::consts::SQRT_2).perfc()
}

/// `z` with `P(Z > z) = p`, by bisection on [`normal_tail`]: `p` in `(0, 0.5]`.
#[cfg(test)]
pub(crate) fn normal_quantile_upper(p: f64) -> f64 {
    let (mut lo, mut hi) = (0.0f64, 40.0f64);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if normal_tail(mid) > p {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// Expected number of local maxima above `u` (in units of the noise) per pixel of white noise
/// smoothed by a Gaussian of standard deviation `sigma_k` pixels.
///
/// Two bounds, of which the smaller is used. For a continuous field the density of maxima
/// above `u` tends to `u·φ(u) / (4π σk²)` (Bond & Efstathiou 1987, MNRAS 226, 655, eq. A1.12;
/// Cheng & Schwartzman 2017, Ann. Stat. 45, 529, with κ = 1 for a Gaussian kernel). On pixels a
/// maximum above `u` is a pixel above `u`, so there are no more of them than `P(Z > u)` per
/// pixel, which is the smaller figure for narrow kernels. Counted on empty synthetic frames,
/// the smaller bound is within 45% of the number of peaks found for σk from 0.7 to 2 px.
pub(crate) fn noise_peaks_above(u: f64, sigma_k: f64) -> f64 {
    let phi = (-0.5 * u * u).pexp() / (2.0 * std::f64::consts::PI).sqrt();
    let peaks = u.max(0.0) * phi / (4.0 * std::f64::consts::PI * sigma_k * sigma_k);
    normal_tail(u).min(peaks.max(0.0))
}

/// The threshold above which `count` noise peaks are expected among `pixels` pixels.
pub(crate) fn threshold_for(count: f64, pixels: f64, sigma_k: f64) -> f64 {
    // The density falls with u above 1: bisect there.
    let (mut lo, mut hi) = (1.0f64, 40.0f64);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if pixels * noise_peaks_above(mid, sigma_k) > count {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// Median of `v` (the upper middle for an even count); `None` when empty. Reorders `v`.
pub(crate) fn median(v: &mut [f64]) -> Option<f64> {
    if v.is_empty() {
        return None;
    }
    let m = v.len() / 2;
    Some(*v.select_nth_unstable_by(m, f64::total_cmp).1)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tail_and_quantile_agree_with_tables() {
        assert!((normal_tail(0.0) - 0.5).abs() < 1e-15);
        assert!((normal_tail(1.0) - 0.158_655_253_931_457).abs() < 1e-12);
        assert!((normal_tail(5.0) - 2.866_515_718_791_94e-7).abs() < 1e-18);
        for p in [0.4, 1e-3, 1e-6, 2.5e-7, 1e-12] {
            let z = normal_quantile_upper(p);
            assert!((normal_tail(z) / p - 1.0).abs() < 1e-9, "{p}");
        }
        assert!((normal_quantile_upper(1e-6) - 4.753_424).abs() < 1e-5);
    }

    #[test]
    fn noise_peak_density_and_threshold() {
        // A narrow kernel: the pixel bound. A wide one: the peak density, which is lower.
        assert_eq!(noise_peaks_above(5.0, 0.7), normal_tail(5.0));
        let wide = noise_peaks_above(5.0, 2.0);
        assert!((wide / 1.479e-7 - 1.0).abs() < 1e-3, "{wide}");
        for (count, s) in [(1.0, 0.7), (1.0, 2.0), (10.0, 1.2)] {
            let u = threshold_for(count, 4e6, s);
            assert!((4e6 * noise_peaks_above(u, s) / count - 1.0).abs() < 1e-6);
        }
        assert!((threshold_for(1.0, 4e6, 1.2) - normal_quantile_upper(2.5e-7)).abs() < 1e-6);
        assert!(threshold_for(1.0, 4e6, 2.0) < threshold_for(1.0, 4e6, 1.2));
    }
}
