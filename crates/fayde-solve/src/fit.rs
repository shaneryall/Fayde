//! Weighted WCS fits with `CRPIX` fixed (design doc §5.5).
//!
//! The model is a polynomial from pixel offsets `(u, v) = pixel − CRPIX` to tangent-plane
//! coordinates at `CRVAL`: `ξ = Σ c_pq u^p v^q`, `η = Σ d_pq u^p v^q`, `0 ≤ p+q ≤ order`.
//! Order 1 is TAN (affine); order ≥ 2 is decomposed into a CD matrix (linear terms) and SIP
//! polynomials `A = CD⁻¹·(higher terms)` (Shupe et al. 2005), the same decomposition
//! astrometry.net uses. Each iteration re-centres the tangent point on the constant term.
//!
//! Weights come from the caller (IRLS: inlier probability / σ², Holland & Welsch 1977; Hogg, Bovy
//! & Lang 2010). The inverse normal matrix, scaled by the fitted error scale², is the parameter
//! covariance; [`WcsCovariance`] turns it into predicted position errors. Order is chosen by BIC
//! (Schwarz 1978), not by an F-test.
//!
//! Instrument priors on the scale and roll can enter the fit as two pseudo-observations on the
//! linear part ([`LinearPrior`]): a Gaussian prior on a linear-model parameter is equivalent to
//! one more weighted observation of it (Gelman et al. 2013, *Bayesian Data Analysis* 3rd ed.,
//! ch. 14, prior information as additional data points). They hold sparse fits near what the
//! instrument is known to do and barely move fits the data determine.

use fayde_core::Portable;
use fayde_core::{PixelPos, SkyPos};
use fayde_wcs::{
    gnomonic, gnomonic_inverse,
    sip::{invert_dense, monomials, solve_dense},
    Parity, Poly, Sip, TanWcs, Wcs,
};

/// One weighted correspondence.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Pair {
    pub px: [f64; 2],
    pub sky: SkyPos,
    pub w: f64,
}

/// Parameter covariance of a fitted WCS.
#[derive(Debug, Clone, PartialEq)]
pub struct WcsCovariance {
    /// Monomials `(p, q)` of the fitted polynomial, in parameter order.
    pub terms: Vec<(usize, usize)>,
    /// Inverse normal matrix (row-major, `terms.len()²`), shared by ξ and η.
    pub inv_normal: Vec<f64>,
    /// Error scale² applied to `inv_normal` (may be < 1 for default centroid errors).
    pub scale2: f64,
    /// Pixel scale, radians per pixel. Weights are 1/σ² in pixels while the fitted values are
    /// tangent-plane radians, so the parameter covariance is `inv_normal × scale2 × rad_per_px²`.
    pub rad_per_px: f64,
    /// Reference pixel.
    pub crpix: PixelPos,
    /// Offset normalization `S` used in the fit: `inv_normal` refers to monomials of
    /// `(u/S, v/S)`.
    pub norm: f64,
}

impl WcsCovariance {
    /// 1σ error of the predicted sky position of pixel `p`, arcsec, per tangent-plane axis
    /// `[ξ, η]`. At `CRPIX` this is the uncertainty of `CRVAL`.
    pub fn position_sigma_arcsec(&self, p: PixelPos) -> f64 {
        (self.basis_quadratic(p) * self.scale2).sqrt() * self.rad_per_px.to_degrees() * 3600.0
    }

    /// fᵀ (AᵀWA)⁻¹ f for the monomial vector f at pixel `p`, without the error scale. Times a
    /// pair's fit weight this is that pair's leverage (the hat-matrix diagonal).
    pub(crate) fn basis_quadratic(&self, p: PixelPos) -> f64 {
        let (u, v) = (
            (p.x - self.crpix.x) / self.norm,
            (p.y - self.crpix.y) / self.norm,
        );
        let f: Vec<f64> = self
            .terms
            .iter()
            .map(|&(a, b)| u.ppowi(a as i32) * v.ppowi(b as i32))
            .collect();
        let k = f.len();
        let mut q = 0.0;
        for i in 0..k {
            for j in 0..k {
                q += f[i] * self.inv_normal[i * k + j] * f[j];
            }
        }
        q.max(0.0)
    }

    /// Worst radial 1σ over a 5×5 grid covering a `width × height` frame, arcsec.
    pub fn worst_radial_sigma_arcsec(&self, width: f64, height: f64) -> f64 {
        let mut worst = 0.0f64;
        for i in 0..5 {
            for j in 0..5 {
                let p = PixelPos::new(
                    -0.5 + width * i as f64 / 4.0,
                    -0.5 + height * j as f64 / 4.0,
                );
                worst = worst.max(self.radial_sigma_arcsec(p));
            }
        }
        worst
    }

    /// Radial 1σ error (√2 × per-axis) of pixel `p`, arcsec.
    pub fn radial_sigma_arcsec(&self, p: PixelPos) -> f64 {
        std::f64::consts::SQRT_2 * self.position_sigma_arcsec(p)
    }
}

/// Gaussian priors on the linear part of a fit, from the instrument: scale and roll, as two
/// pseudo-observations linearised at the prior CD matrix
/// `CD₀ = s₀ · [[f cos θ₀, sin θ₀], [−f sin θ₀, cos θ₀]]` (f = −1 for negative parity).
/// The scale is observed along CD₀, the roll along ∂CD₀/∂θ; the two remaining directions of the
/// CD matrix (skew and pixel anisotropy) stay free.
#[derive(Debug, Clone, Copy)]
pub(crate) struct LinearPrior {
    /// Scale and 1σ, radians per pixel.
    pub scale: Option<(f64, f64)>,
    /// Roll `atan2(CD1_2, CD2_2)` and 1σ, radians.
    pub roll: Option<(f64, f64)>,
    /// Parity, when known; otherwise that of the data-only fit.
    pub parity: Option<Parity>,
    /// The factor between the caller's pair weights and true inverse variances
    /// (`w = w_true × weight_scale2`), so the priors carry their stated weight.
    pub weight_scale2: f64,
}

/// Re-solve the linear system of [`fit_poly`] for ξ and η jointly with the prior rows. `cx`, `cy`
/// are the data-only solution (fit units: coefficients of normalised offsets `u/norm`).
#[allow(clippy::too_many_arguments)]
fn solve_with_prior(
    normal: &[f64],
    vx: &[f64],
    vy: &[f64],
    k: usize,
    cx: &[f64],
    cy: &[f64],
    norm: f64,
    pr: &LinearPrior,
) -> Option<(Vec<f64>, Vec<f64>)> {
    let (c1, c2, d1, d2) = (cx[1], cx[2], cy[1], cy[2]);
    let det = c1 * d2 - c2 * d1;
    let flip = match pr.parity {
        Some(Parity::Negative) => -1.0,
        Some(Parity::Positive) => 1.0,
        None if det < 0.0 => -1.0,
        None => 1.0,
    };
    // Reference scale (fit units) and roll of the linearisation point.
    let s_ref = pr.scale.map_or(det.abs().sqrt(), |(s, _)| s * norm);
    let theta = pr.roll.map_or(c2.patan2(d2), |(r, _)| r);
    if !(s_ref > 0.0 && s_ref.is_finite() && theta.is_finite()) {
        return None;
    }
    let (sn, cs) = theta.psin_cos();
    // Weights: the data weights are true inverse variances (in pixels) × weight_scale2, and a
    // fit-unit coefficient error is rad_per_px × norm per pixel.
    let rad_per_px = s_ref / norm;
    let w0 = rad_per_px * rad_per_px * pr.weight_scale2;
    // Rows over (c1, c2, d1, d2): (row, target, weight).
    let mut rows: Vec<([f64; 4], f64, f64)> = Vec::new();
    if let Some((s0, sigma)) = pr.scale.filter(|&(_, sg)| sg > 0.0) {
        // ŝ = ⟨m, CD₀⟩ / |CD₀|², times s_ref: exact at CD₀, first order around it.
        rows.push((
            [flip * cs / 2.0, sn / 2.0, -flip * sn / 2.0, cs / 2.0],
            s0 * norm,
            w0 / (sigma * norm).ppowi(2),
        ));
    }
    if let Some((_, sigma)) = pr.roll.filter(|&(_, sg)| sg > 0.0) {
        // δθ = ⟨m, ∂CD₀/∂θ⟩ / |∂CD₀/∂θ|²; ⟨CD₀, ∂CD₀/∂θ⟩ = 0, so δθ = 0 at CD₀.
        let t = 1.0 / (2.0 * s_ref);
        rows.push((
            [-flip * sn * t, cs * t, -flip * cs * t, -sn * t],
            0.0,
            w0 / (sigma * sigma),
        ));
    }
    if rows.is_empty() {
        return None;
    }
    let n = 2 * k;
    let mut m = vec![0.0; n * n];
    let mut v = vec![0.0; n];
    for i in 0..k {
        for j in 0..k {
            m[i * n + j] = normal[i * k + j];
            m[(k + i) * n + (k + j)] = normal[i * k + j];
        }
        v[i] = vx[i];
        v[k + i] = vy[i];
    }
    let idx = [1, 2, k + 1, k + 2];
    for (r, target, w) in rows {
        for a in 0..4 {
            v[idx[a]] += w * r[a] * target;
            for b in 0..4 {
                m[idx[a] * n + idx[b]] += w * r[a] * r[b];
            }
        }
    }
    let x = solve_dense(m, v)?;
    Some((x[..k].to_vec(), x[k..].to_vec()))
}

/// Result of one polynomial fit.
pub(crate) struct PolyFit {
    pub wcs: Wcs,
    pub cov: WcsCovariance,
}

/// Fit a polynomial WCS of `order` (1 = TAN). Needs at least as many weighted pairs as parameters
/// per axis (exactly as many gives the interpolating solution). With `prior`, the linear part is
/// also held to the instrument's scale and roll. The covariance stays that of the data alone
/// (conservative: the priors only add information).
pub(crate) fn fit_poly(
    crval0: SkyPos,
    crpix: PixelPos,
    pairs: &[Pair],
    order: usize,
    prior: Option<&LinearPrior>,
) -> Option<PolyFit> {
    let terms = monomials(0, order);
    let k = terms.len();
    if pairs.iter().filter(|p| p.w > 0.0).count() < k {
        return None;
    }
    // Conditioning: fit in normalized offsets (u/S, v/S) with S the largest |offset|, so every
    // monomial is O(1) and the normal matrix stays well conditioned at any order and detector
    // size (raw u³ on a 2048 px chip spans ~1e18 against the constant term). Coefficients are
    // rescaled by S^-(p+q) afterwards.
    let norm = pairs
        .iter()
        .filter(|p| p.w > 0.0)
        .map(|p| (p.px[0] - crpix.x).abs().max((p.px[1] - crpix.y).abs()))
        .fold(1.0f64, f64::max);
    let mut crval = crval0;
    let mut cx = vec![0.0; k];
    let mut cy = vec![0.0; k];
    let mut normal = vec![0.0; k * k];
    for _ in 0..5 {
        normal = vec![0.0; k * k];
        let mut vx = vec![0.0; k];
        let mut vy = vec![0.0; k];
        for p in pairs.iter().filter(|p| p.w > 0.0) {
            let (xi, eta) = gnomonic(crval, p.sky)?;
            let (u, v) = ((p.px[0] - crpix.x) / norm, (p.px[1] - crpix.y) / norm);
            let f: Vec<f64> = terms
                .iter()
                .map(|&(a, b)| u.ppowi(a as i32) * v.ppowi(b as i32))
                .collect();
            for i in 0..k {
                for j in 0..k {
                    normal[i * k + j] += p.w * f[i] * f[j];
                }
                vx[i] += p.w * f[i] * xi;
                vy[i] += p.w * f[i] * eta;
            }
        }
        cx = solve_dense(normal.clone(), vx.clone())?;
        cy = solve_dense(normal.clone(), vy.clone())?;
        if let Some(pr) = prior {
            if let Some((px, py)) = solve_with_prior(&normal, &vx, &vy, k, &cx, &cy, norm, pr) {
                (cx, cy) = (px, py);
            }
        }
        // Term 0 is the constant: the tangent-point offset.
        crval = gnomonic_inverse(crval, cx[0], cy[0]);
        if cx[0].phypot(cy[0]) < 1e-13 {
            cx[0] = 0.0;
            cy[0] = 0.0;
            break;
        }
    }
    // Back to pixel units.
    for (i, &(a, b)) in terms.iter().enumerate() {
        let s = norm.ppowi((a + b) as i32);
        cx[i] /= s;
        cy[i] /= s;
    }
    // Terms 1 and 2 are (1,0) and (0,1): the linear part.
    let cd_rad = [[cx[1], cx[2]], [cy[1], cy[2]]];
    let tan = TanWcs::new(crval, crpix, cd_rad.map(|r| r.map(f64::to_degrees))).ok()?;
    let sip = (order >= 2).then(|| {
        let det = cd_rad[0][0] * cd_rad[1][1] - cd_rad[0][1] * cd_rad[1][0];
        let inv = [
            [cd_rad[1][1] / det, -cd_rad[0][1] / det],
            [-cd_rad[1][0] / det, cd_rad[0][0] / det],
        ];
        let mut a = Vec::new();
        let mut b = Vec::new();
        for (i, &(p, q)) in terms.iter().enumerate().skip(3) {
            a.push((p, q, inv[0][0] * cx[i] + inv[0][1] * cy[i]));
            b.push((p, q, inv[1][0] * cx[i] + inv[1][1] * cy[i]));
        }
        Sip {
            a: Poly { order, terms: a },
            b: Poly { order, terms: b },
            ap: Poly {
                order: 0,
                terms: vec![],
            },
            bp: Poly {
                order: 0,
                terms: vec![],
            },
        }
    });
    let inv_normal = invert_dense(&normal, k)?;
    Some(PolyFit {
        wcs: Wcs { tan, sip },
        cov: WcsCovariance {
            terms,
            inv_normal,
            scale2: 1.0,
            rad_per_px: (cd_rad[0][0] * cd_rad[1][1] - cd_rad[0][1] * cd_rad[1][0])
                .abs()
                .sqrt(),
            crpix,
            norm,
        },
    })
}

/// Affine TAN fit (order 1).
#[cfg(test)]
pub(crate) fn fit_tan(crval0: SkyPos, crpix: PixelPos, pairs: &[Pair]) -> Option<TanWcs> {
    fit_poly(crval0, crpix, pairs, 1, None).map(|f| f.wcs.tan)
}

/// Weighted χ² of `wcs` over `pairs`, in pixels (residual / σ through the weights).
pub(crate) fn chi2(wcs: &Wcs, pairs: &[Pair]) -> Option<f64> {
    let mut s = 0.0;
    for p in pairs.iter().filter(|p| p.w > 0.0) {
        let q = wcs.sky_to_pixel(p.sky)?;
        s += p.w * ((p.px[0] - q.x).ppowi(2) + (p.px[1] - q.y).ppowi(2));
    }
    Some(s)
}

/// Observations (2 per match) per fitted parameter a distortion order needs to be a BIC
/// candidate in a fit that is returned: n ≥ 4k, i.e. 24 matches for order 2 and 40 for order 3.
/// See [`fit_best_order`].
pub(crate) const FINAL_OBS_PER_PARAM: f64 = 4.0;
/// The same limit while a tune-up is still collecting matches (n ≥ 2k: 12 and 20 matches). Its
/// result is refitted under [`FINAL_OBS_PER_PARAM`] before it can be returned.
pub(crate) const GROWTH_OBS_PER_PARAM: f64 = 2.0;

/// Fit orders `1..=max_order` and keep the lowest BIC = n·ln(χ²/n) + k·ln(n), with k = 2 × terms and
/// n = 2 × (sum of inlier probabilities). Returns the fit and its χ².
///
/// An order above 1 is a candidate only with `obs_per_param` observations per parameter
/// (n ≥ `obs_per_param` × k). The profiled BIC has no floor: as k approaches n the fit
/// interpolates, χ² → 0 and n·ln(χ²/n) → −∞, so without a rule a frame with exactly 6 or 10
/// matches always picks the interpolating polynomial, which is wrong by many pixels away from
/// its stars. Short of interpolation it is still too liberal for a few dozen matches with
/// outliers among them: on 1,200 undistorted frames of 12–50 stars, n ≥ 2k chose distortion 16
/// times (corner errors up to 107″), n ≥ 3k 5 times (84″), n ≥ 4k once (3″). A BIC margin
/// instead of a count did not help (still 3 false choices at a margin of 50). Fayde's own rule.
pub(crate) fn fit_best_order(
    crval0: SkyPos,
    crpix: PixelPos,
    pairs: &[Pair],
    inlier_sum: f64,
    max_order: usize,
    obs_per_param: f64,
    prior: Option<&LinearPrior>,
) -> Option<(PolyFit, f64)> {
    let n = (2.0 * inlier_sum).max(1.0);
    let mut best: Option<(f64, PolyFit, f64)> = None;
    for order in 1..=max_order.max(1) {
        if order > 1 && n < obs_per_param * (2 * monomials(0, order).len()) as f64 {
            break;
        }
        let Some(fit) = fit_poly(crval0, crpix, pairs, order, prior) else {
            break;
        };
        let Some(c2) = chi2(&fit.wcs, pairs) else {
            break;
        };
        // Profile out the unknown noise variance (Schwarz 1978). A common rescaling
        // of the weights adds the same constant to every candidate BIC.
        let bic =
            n * (c2 / n).max(f64::MIN_POSITIVE).pln() + (2 * fit.cov.terms.len()) as f64 * n.pln();
        if best.as_ref().is_none_or(|b| bic < b.0) {
            best = Some((bic, fit, c2));
        }
    }
    best.map(|(_, f, c2)| (f, c2))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn grid(truth: &Wcs) -> Vec<Pair> {
        (0..60)
            .map(|k| {
                let px = [(k * 37 % 1000) as f64 + 3.3, (k * 91 % 1000) as f64 + 7.1];
                Pair {
                    px,
                    sky: truth.pixel_to_sky(PixelPos::new(px[0], px[1])),
                    w: 1.0,
                }
            })
            .collect()
    }

    #[test]
    fn recovers_a_known_tan() {
        let truth = TanWcs::new(
            SkyPos::from_degrees(83.8, -5.4).unwrap(),
            PixelPos::new(511.5, 511.5),
            [[-7.2e-4, 4.1e-4], [4.1e-4, 7.2e-4]],
        )
        .unwrap();
        let pairs = grid(&truth.into());
        let start = SkyPos::from_degrees(84.1, -5.2).unwrap();
        let fit = fit_tan(start, truth.crpix(), &pairs).unwrap();
        assert!(fayde_core::to_arcsec(fit.crval().separation(truth.crval())) < 1e-6);
        for i in 0..2 {
            for j in 0..2 {
                assert!((fit.cd()[i][j] - truth.cd()[i][j]).abs() < 1e-12);
            }
        }
        assert!(fit_tan(start, truth.crpix(), &pairs[..2]).is_none());
        assert!(fit_tan(start, truth.crpix(), &pairs[..3]).is_some());
    }

    #[test]
    fn recovers_sip_and_bic_picks_it() {
        let tan = TanWcs::new(
            SkyPos::from_degrees(150.0, 2.0).unwrap(),
            PixelPos::new(511.5, 511.5),
            [[-8.3e-4, 1e-5], [1e-5, 8.3e-4]],
        )
        .unwrap();
        let k = 2e-8;
        let truth = Wcs {
            tan,
            sip: Some(Sip {
                a: Poly {
                    order: 3,
                    terms: vec![(3, 0, k), (1, 2, k)],
                },
                b: Poly {
                    order: 3,
                    terms: vec![(2, 1, k), (0, 3, k)],
                },
                ap: Poly {
                    order: 0,
                    terms: vec![],
                },
                bp: Poly {
                    order: 0,
                    terms: vec![],
                },
            }),
        };
        // σ = 0.1 px → weight 100.
        let pairs: Vec<Pair> = grid(&truth)
            .into_iter()
            .map(|p| Pair { w: 100.0, ..p })
            .collect();
        let (fit, _) = fit_best_order(
            tan.crval(),
            tan.crpix(),
            &pairs,
            60.0,
            3,
            FINAL_OBS_PER_PARAM,
            None,
        )
        .unwrap();
        assert_eq!(fit.wcs.sip_order(), 3);
        for (x, y) in [(0.0, 0.0), (1023.0, 1023.0), (100.0, 900.0)] {
            let p = PixelPos::new(x, y);
            assert!(fit.wcs.pixel_to_sky(p).separation(truth.pixel_to_sky(p)) < 1e-12);
        }
        // Without distortion, BIC keeps TAN.
        let flat: Vec<Pair> = grid(&tan.into())
            .into_iter()
            .enumerate()
            .map(|(i, p)| Pair {
                px: [
                    p.px[0] + 0.1 * ((i * 7 % 11) as f64 / 5.0 - 1.0),
                    p.px[1] + 0.1 * ((i * 5 % 13) as f64 / 6.0 - 1.0),
                ],
                w: 100.0,
                ..p
            })
            .collect();
        let (fit, _) = fit_best_order(
            tan.crval(),
            tan.crpix(),
            &flat,
            60.0,
            3,
            FINAL_OBS_PER_PARAM,
            None,
        )
        .unwrap();
        assert_eq!(fit.wcs.sip_order(), 1);
    }

    #[test]
    fn interpolating_orders_are_not_bic_candidates() {
        // A plain TAN with 0.3 px noise and σ supplied. With 6 pairs order 2 interpolates, with
        // 10 pairs order 3 does: χ² → 0 and the profiled BIC would pick them every time.
        let tan = TanWcs::new(
            SkyPos::from_degrees(150.0, 2.0).unwrap(),
            PixelPos::new(511.5, 511.5),
            [[-8.3e-4, 1e-5], [1e-5, 8.3e-4]],
        )
        .unwrap();
        let noisy: Vec<Pair> = grid(&tan.into())
            .into_iter()
            .enumerate()
            .map(|(i, p)| Pair {
                px: [
                    p.px[0] + 0.3 * ((i * 7 % 11) as f64 / 5.0 - 1.0),
                    p.px[1] + 0.3 * ((i * 5 % 13) as f64 / 6.0 - 1.0),
                ],
                w: 1.0 / 0.09,
                ..p
            })
            .collect();
        // `n` pairs spread over the frame (the first few grid points are nearly collinear).
        let spread = |n: usize| -> Vec<Pair> { (0..n).map(|i| noisy[i * (60 / n)]).collect() };
        for n in [6usize, 7, 8, 10, 11] {
            let pairs = spread(n);
            let (fit, _) = fit_best_order(
                tan.crval(),
                tan.crpix(),
                &pairs,
                n as f64,
                3,
                FINAL_OBS_PER_PARAM,
                None,
            )
            .unwrap();
            assert_eq!(fit.wcs.sip_order(), 1, "{n} pairs");
            let corner = PixelPos::new(0.0, 1023.0);
            let err = fayde_core::to_arcsec(
                fit.wcs
                    .pixel_to_sky(corner)
                    .separation(tan.pixel_to_sky(corner)),
            );
            assert!(err < 15.0, "{n} pairs: corner off by {err}″");
        }
        // Order 2 becomes a candidate at 24 matches, order 3 at 40; noise alone still keeps TAN.
        for n in [12usize, 20, 30, 60] {
            let pairs = spread(n);
            let (fit, _) = fit_best_order(
                tan.crval(),
                tan.crpix(),
                &pairs,
                n as f64,
                3,
                FINAL_OBS_PER_PARAM,
                None,
            )
            .unwrap();
            assert_eq!(fit.wcs.sip_order(), 1, "{n} pairs");
        }
    }

    #[test]
    fn large_detector_order3_is_well_conditioned() {
        // 4096 px: raw u³ reaches ~1e10 and its square ~1e20 in the normal matrix; without
        // normalized offsets the cubic terms lose precision.
        let tan = TanWcs::new(
            SkyPos::from_degrees(30.0, -60.0).unwrap(),
            PixelPos::new(2047.5, 2047.5),
            [[-2.8e-4, 3e-6], [3e-6, 2.8e-4]],
        )
        .unwrap();
        let k = 1.2e-9; // ~0.9 px of barrel at the corners
        let truth = Wcs {
            tan,
            sip: Some(Sip {
                a: Poly {
                    order: 3,
                    terms: vec![(3, 0, k), (1, 2, k), (2, 0, 3e-7)],
                },
                b: Poly {
                    order: 3,
                    terms: vec![(2, 1, k), (0, 3, k), (1, 1, -2e-7)],
                },
                ap: Poly {
                    order: 0,
                    terms: vec![],
                },
                bp: Poly {
                    order: 0,
                    terms: vec![],
                },
            }),
        };
        let pairs: Vec<Pair> = (0..200)
            .map(|i| {
                let px = [
                    (i * 397 % 4096) as f64 + 0.25,
                    (i * 1009 % 4096) as f64 + 0.75,
                ];
                Pair {
                    px,
                    sky: truth.pixel_to_sky(PixelPos::new(px[0], px[1])),
                    w: 1.0,
                }
            })
            .collect();
        let fit = fit_poly(tan.crval(), tan.crpix(), &pairs, 3, None).unwrap();
        for (x, y) in [
            (0.0, 0.0),
            (4095.0, 4095.0),
            (0.0, 4095.0),
            (2047.5, 2047.5),
        ] {
            let p = PixelPos::new(x, y);
            let err =
                fayde_core::to_arcsec(fit.wcs.pixel_to_sky(p).separation(truth.pixel_to_sky(p)));
            assert!(err < 1e-6, "{err}″ at {p:?}");
        }
    }

    #[test]
    fn covariance_shrinks_with_more_stars_and_grows_off_centre() {
        let truth: Wcs = TanWcs::new(
            SkyPos::from_degrees(10.0, 10.0).unwrap(),
            PixelPos::new(511.5, 511.5),
            [[-8e-4, 0.0], [0.0, 8e-4]],
        )
        .unwrap()
        .into();
        let pairs: Vec<Pair> = grid(&truth)
            .into_iter()
            .map(|p| Pair { w: 1.0 / 0.09, ..p })
            .collect();
        let f60 = fit_poly(truth.crval(), truth.crpix(), &pairs, 1, None)
            .unwrap()
            .cov;
        let f15 = fit_poly(truth.crval(), truth.crpix(), &pairs[..15], 1, None)
            .unwrap()
            .cov;
        let c = PixelPos::new(511.5, 511.5);
        let corner = PixelPos::new(0.0, 0.0);
        assert!(f15.position_sigma_arcsec(c) > 1.5 * f60.position_sigma_arcsec(c));
        assert!(f60.position_sigma_arcsec(corner) > f60.position_sigma_arcsec(c));
        // Centre error ≈ σ/√N × scale: 0.3 px / √60 × 2.88″.
        let expect = 0.3 / 60f64.sqrt() * 2.88;
        let got = f60.position_sigma_arcsec(c);
        assert!((got / expect - 1.0).abs() < 0.3, "{got} vs {expect}");
    }
}
