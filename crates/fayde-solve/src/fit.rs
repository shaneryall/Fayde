use fayde_core::Portable;
use fayde_core::{PixelPos, SkyPos, UnitVector};
use fayde_wcs::{
    gnomonic_inverse,
    sip::{monomials, solve_dense, DenseFactor},
    Parity, Poly, Sip, TanWcs, TangentPlane, Wcs,
};

#[derive(Debug, Clone, Copy)]
pub(crate) struct Pair {
    pub px: [f64; 2],
    pub vec: [f64; 3],
    pub w: f64,
}

impl Pair {
    #[cfg(test)]
    pub(crate) fn at(px: [f64; 2], sky: SkyPos, w: f64) -> Self {
        Self {
            px,
            vec: sky.to_unit_vector().0,
            w,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct WcsCovariance {
    pub terms: Vec<(usize, usize)>,
    pub inv_normal: Vec<f64>,
    pub scale2: f64,
    pub rad_per_px: f64,
    pub crpix: PixelPos,
    pub norm: f64,
}

impl WcsCovariance {
    pub fn position_sigma_arcsec(&self, p: PixelPos) -> f64 {
        (self.basis_quadratic(p) * self.scale2).sqrt() * self.rad_per_px.to_degrees() * 3600.0
    }

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

    pub fn radial_sigma_arcsec(&self, p: PixelPos) -> f64 {
        std::f64::consts::SQRT_2 * self.position_sigma_arcsec(p)
    }
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct LinearPrior {
    pub scale: Option<(f64, f64)>,
    pub roll: Option<(f64, f64)>,
    pub parity: Option<Parity>,
    pub weight_scale2: f64,
}

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
    let s_ref = pr.scale.map_or(det.abs().sqrt(), |(s, _)| s * norm);
    let theta = pr.roll.map_or(c2.patan2(d2), |(r, _)| r);
    if !(s_ref > 0.0 && s_ref.is_finite() && theta.is_finite()) {
        return None;
    }
    let (sn, cs) = theta.psin_cos();
    let rad_per_px = s_ref / norm;
    let w0 = rad_per_px * rad_per_px * pr.weight_scale2;
    let mut rows: Vec<([f64; 4], f64, f64)> = Vec::new();
    if let Some((s0, sigma)) = pr.scale.filter(|&(_, sg)| sg > 0.0) {
        rows.push((
            [flip * cs / 2.0, sn / 2.0, -flip * sn / 2.0, cs / 2.0],
            s0 * norm,
            w0 / (sigma * norm).ppowi(2),
        ));
    }
    if let Some((_, sigma)) = pr.roll.filter(|&(_, sg)| sg > 0.0) {
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

pub(crate) struct PolyFit {
    pub wcs: Wcs,
    pub cov: WcsCovariance,
}

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
    let norm = pairs
        .iter()
        .filter(|p| p.w > 0.0)
        .map(|p| (p.px[0] - crpix.x).abs().max((p.px[1] - crpix.y).abs()))
        .fold(1.0f64, f64::max);
    let fitted: Vec<&Pair> = pairs.iter().filter(|p| p.w > 0.0).collect();
    let mut monomial = vec![0.0; fitted.len() * k];
    let mut pu = vec![0.0; order + 1];
    let mut pv = vec![0.0; order + 1];
    for (p, f) in fitted.iter().zip(monomial.chunks_exact_mut(k)) {
        let (u, v) = ((p.px[0] - crpix.x) / norm, (p.px[1] - crpix.y) / norm);
        for d in 0..=order {
            pu[d] = u.ppowi(d as i32);
            pv[d] = v.ppowi(d as i32);
        }
        for (fi, &(a, b)) in f.iter_mut().zip(&terms) {
            *fi = pu[a] * pv[b];
        }
    }
    let mut normal = vec![0.0; k * k];
    for (p, f) in fitted.iter().zip(monomial.chunks_exact(k)) {
        for i in 0..k {
            for j in 0..k {
                normal[i * k + j] += p.w * f[i] * f[j];
            }
        }
    }
    let (factor, prior) = match prior {
        Some(prior) => (DenseFactor::new(normal.clone(), k)?, Some((prior, normal))),
        None => (DenseFactor::new(normal, k)?, None),
    };
    let mut crval = crval0;
    let mut cx = vec![0.0; k];
    let mut cy = vec![0.0; k];
    for _ in 0..5 {
        let mut vx = vec![0.0; k];
        let mut vy = vec![0.0; k];
        let plane = TangentPlane::at(crval);
        for (p, f) in fitted.iter().zip(monomial.chunks_exact(k)) {
            let (xi, eta) = plane.project(UnitVector(p.vec))?;
            for i in 0..k {
                vx[i] += p.w * f[i] * xi;
                vy[i] += p.w * f[i] * eta;
            }
        }
        (cx, cy) = match &prior {
            Some((prior, normal)) => {
                let cx = factor.solve(vx.clone())?;
                let cy = factor.solve(vy.clone())?;
                solve_with_prior(normal, &vx, &vy, k, &cx, &cy, norm, prior).unwrap_or((cx, cy))
            }
            None => (factor.solve(vx)?, factor.solve(vy)?),
        };
        crval = gnomonic_inverse(crval, cx[0], cy[0]);
        if cx[0].phypot(cy[0]) < 1e-13 {
            cx[0] = 0.0;
            cy[0] = 0.0;
            break;
        }
    }
    for (i, &(a, b)) in terms.iter().enumerate() {
        let s = norm.ppowi((a + b) as i32);
        cx[i] /= s;
        cy[i] /= s;
    }
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
    let inv_normal = factor.inverse()?;
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

#[cfg(test)]
pub(crate) fn fit_tan(crval0: SkyPos, crpix: PixelPos, pairs: &[Pair]) -> Option<TanWcs> {
    fit_poly(crval0, crpix, pairs, 1, None).map(|f| f.wcs.tan)
}

pub(crate) fn chi2(wcs: &Wcs, pairs: &[Pair]) -> Option<f64> {
    let mut s = 0.0;
    for p in pairs.iter().filter(|p| p.w > 0.0) {
        let q = wcs.vector_to_pixel(UnitVector(p.vec))?;
        s += p.w * ((p.px[0] - q.x).ppowi(2) + (p.px[1] - q.y).ppowi(2));
    }
    Some(s)
}

pub(crate) const FINAL_OBS_PER_PARAM: f64 = 4.0;
pub(crate) const GROWTH_OBS_PER_PARAM: f64 = 2.0;

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
        let bic =
            n * (c2 / n).max(f64::MIN_POSITIVE).pln() + (2 * fit.cov.terms.len()) as f64 * n.pln();
        if best.as_ref().is_none_or(|b| bic < b.0) {
            best = Some((bic, fit, c2));
        }
    }
    best.map(|(_, f, c2)| (f, c2))
}

pub(crate) const MIN_DISTORTION_CELLS: usize = 4;

#[allow(clippy::too_many_arguments)]
pub(crate) fn held_out_error(
    crval0: SkyPos,
    crpix: PixelPos,
    pairs: &[Pair],
    order: usize,
    prior: Option<&LinearPrior>,
    cells: usize,
    size: [f64; 2],
    obs_per_param: f64,
) -> Option<f64> {
    let cell = |px: [f64; 2]| {
        let at = |v: f64, n: f64| (((v + 0.5) / n * cells as f64).max(0.0) as usize).min(cells - 1);
        at(px[1], size[1]) * cells + at(px[0], size[0])
    };
    let params = (2 * monomials(0, order).len()) as f64;
    let (mut sum, mut wsum) = (0.0f64, 0.0f64);
    for c in 0..cells * cells {
        let test: Vec<&Pair> = pairs
            .iter()
            .filter(|p| p.w > 0.0 && cell(p.px) == c)
            .collect();
        if test.is_empty() {
            continue;
        }
        let train: Vec<Pair> = pairs.iter().filter(|p| cell(p.px) != c).copied().collect();
        let n = 2.0 * train.iter().filter(|p| p.w > 0.0).count() as f64;
        if order > 1 && n < obs_per_param * params {
            return None;
        }
        let fit = fit_poly(crval0, crpix, &train, order, prior)?;
        for p in test {
            let q = fit.wcs.vector_to_pixel(UnitVector(p.vec))?;
            sum += p.w * ((p.px[0] - q.x).ppowi(2) + (p.px[1] - q.y).ppowi(2));
            wsum += p.w;
        }
    }
    (wsum > 0.0 && sum.is_finite()).then(|| sum / wsum)
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn held_out_order(
    crval0: SkyPos,
    crpix: PixelPos,
    pairs: &[Pair],
    order: usize,
    prior: Option<&LinearPrior>,
    cells: usize,
    size: [f64; 2],
    obs_per_param: f64,
    margin: f64,
) -> usize {
    let cell = |px: [f64; 2]| {
        let at = |v: f64, n: f64| (((v + 0.5) / n * cells as f64).max(0.0) as usize).min(cells - 1);
        at(px[1], size[1]) * cells + at(px[0], size[0])
    };
    let mut occupied = vec![false; cells * cells];
    for p in pairs.iter().filter(|p| p.w > 0.0) {
        occupied[cell(p.px)] = true;
    }
    if order > 1 && occupied.iter().filter(|&&o| o).count() < MIN_DISTORTION_CELLS {
        return 1;
    }
    let errs: Vec<Option<f64>> = (1..=order)
        .map(|o| held_out_error(crval0, crpix, pairs, o, prior, cells, size, obs_per_param))
        .collect();
    if errs[order - 1].is_none() {
        return order;
    }
    let Some(best) = errs.iter().flatten().copied().reduce(f64::min) else {
        return order;
    };
    errs.iter()
        .position(|e| e.is_some_and(|e| e <= margin * best))
        .map_or(order, |i| i + 1)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn grid(truth: &Wcs) -> Vec<Pair> {
        (0..60)
            .map(|k| {
                let px = [(k * 37 % 1000) as f64 + 3.3, (k * 91 % 1000) as f64 + 7.1];
                Pair::at(px, truth.pixel_to_sky(PixelPos::new(px[0], px[1])), 1.0)
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
    fn held_out_cells_reject_distortion_fitted_to_a_corner() {
        let tan = TanWcs::new(
            SkyPos::from_degrees(150.0, 2.0).unwrap(),
            PixelPos::new(511.5, 511.5),
            [[-8.3e-4, 1e-5], [1e-5, 8.3e-4]],
        )
        .unwrap();
        let truth: Wcs = tan.into();
        let pairs: Vec<Pair> = grid(&truth)
            .into_iter()
            .chain((0..40).map(|k| {
                let px = [700.0 + (k * 37 % 300) as f64, 700.0 + (k * 53 % 300) as f64];
                Pair::at(px, truth.pixel_to_sky(PixelPos::new(px[0], px[1])), 1.0)
            }))
            .map(|p| {
                let (u, v) = ((p.px[0] - 680.0).max(0.0), (p.px[1] - 680.0).max(0.0));
                let bend = 2e-7 * u * u * (u.min(v) / 300.0);
                Pair {
                    px: [p.px[0] + bend, p.px[1] + 0.5 * bend],
                    w: 100.0,
                    ..p
                }
            })
            .collect();
        let size = [1024.0, 1024.0];
        let (fit, _) = fit_best_order(
            tan.crval(),
            tan.crpix(),
            &pairs,
            100.0,
            3,
            FINAL_OBS_PER_PARAM,
            None,
        )
        .unwrap();
        assert!(fit.wcs.sip_order() > 1);
        let keep = held_out_order(
            tan.crval(),
            tan.crpix(),
            &pairs,
            fit.wcs.sip_order(),
            None,
            3,
            size,
            GROWTH_OBS_PER_PARAM,
            1.05,
        );
        assert_eq!(keep, 1);
        let corner: Vec<Pair> = pairs
            .iter()
            .filter(|p| p.px[0] > 340.0 && p.px[1] > 680.0)
            .copied()
            .collect();
        let keep = held_out_order(
            tan.crval(),
            tan.crpix(),
            &corner,
            3,
            None,
            3,
            size,
            GROWTH_OBS_PER_PARAM,
            1.05,
        );
        assert_eq!(keep, 1);
    }

    #[test]
    fn held_out_cells_keep_real_distortion() {
        let tan = TanWcs::new(
            SkyPos::from_degrees(150.0, 2.0).unwrap(),
            PixelPos::new(511.5, 511.5),
            [[-8.3e-4, 1e-5], [1e-5, 8.3e-4]],
        )
        .unwrap();
        let k = 2e-8;
        let poly = |terms| Poly { order: 3, terms };
        let truth = Wcs {
            tan,
            sip: Some(Sip {
                a: poly(vec![(3, 0, k), (1, 2, k)]),
                b: poly(vec![(2, 1, k), (0, 3, k)]),
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
        let pairs: Vec<Pair> = grid(&truth)
            .into_iter()
            .map(|p| Pair { w: 100.0, ..p })
            .collect();
        let size = [1024.0, 1024.0];
        let err = |o| {
            held_out_error(
                tan.crval(),
                tan.crpix(),
                &pairs,
                o,
                None,
                3,
                size,
                GROWTH_OBS_PER_PARAM,
            )
        };
        let (e1, e3) = (err(1).unwrap(), err(3).unwrap());
        assert!(e3.is_finite() && e3 * 1.05 < e1, "TAN {e1}, cubic {e3}");
        assert!(err(2).is_none_or(|e2| e3 * 1.05 < e2));
        let keep = held_out_order(
            tan.crval(),
            tan.crpix(),
            &pairs,
            3,
            None,
            3,
            size,
            GROWTH_OBS_PER_PARAM,
            1.05,
        );
        assert_eq!(keep, 3);
    }

    #[test]
    fn interpolating_orders_are_not_bic_candidates() {
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
        let tan = TanWcs::new(
            SkyPos::from_degrees(30.0, -60.0).unwrap(),
            PixelPos::new(2047.5, 2047.5),
            [[-2.8e-4, 3e-6], [3e-6, 2.8e-4]],
        )
        .unwrap();
        let k = 1.2e-9;
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
                Pair::at(px, truth.pixel_to_sky(PixelPos::new(px[0], px[1])), 1.0)
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
        let expect = 0.3 / 60f64.sqrt() * 2.88;
        let got = f60.position_sigma_arcsec(c);
        assert!((got / expect - 1.0).abs() < 0.3, "{got} vs {expect}");
    }
}
