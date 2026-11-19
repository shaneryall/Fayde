use fayde_core::Portable;
use fayde_core::{PixelPos, SkyPos, UnitVector};

use crate::{gnomonic, gnomonic_inverse, Parity, TanWcs, WcsError};

#[derive(Debug, Clone, PartialEq)]
pub struct Poly {
    pub order: usize,
    pub terms: Vec<(usize, usize, f64)>,
}

impl Poly {
    pub fn eval(&self, u: f64, v: f64) -> f64 {
        self.eval_at(&Powers::of(u, self.order), &Powers::of(v, self.order))
    }

    pub fn grad(&self, u: f64, v: f64) -> (f64, f64) {
        self.grad_at(&Powers::of(u, self.order), &Powers::of(v, self.order))
    }

    fn eval_at(&self, u: &Powers, v: &Powers) -> f64 {
        self.terms
            .iter()
            .map(|&(p, q, c)| c * u.get(p) * v.get(q))
            .sum()
    }

    fn grad_at(&self, u: &Powers, v: &Powers) -> (f64, f64) {
        self.terms.iter().fold((0.0, 0.0), |(du, dv), &(p, q, c)| {
            let gu = if p > 0 {
                c * p as f64 * u.get(p - 1) * v.get(q)
            } else {
                0.0
            };
            let gv = if q > 0 {
                c * q as f64 * u.get(p) * v.get(q - 1)
            } else {
                0.0
            };
            (du + gu, dv + gv)
        })
    }
}

const MAX_ORDER: usize = 9;

struct Powers {
    x: f64,
    table: [f64; MAX_ORDER + 1],
    known: usize,
}

impl Powers {
    fn of(x: f64, order: usize) -> Self {
        let known = order.min(MAX_ORDER) + 1;
        let mut table = [0.0; MAX_ORDER + 1];
        for (n, t) in table.iter_mut().enumerate().take(known) {
            *t = x.ppowi(n as i32);
        }
        Self { x, table, known }
    }

    fn get(&self, n: usize) -> f64 {
        if n < self.known {
            self.table[n]
        } else {
            self.x.ppowi(n as i32)
        }
    }
}

const SIP_INVERSE_TOL_PX: f64 = 1e-8;

pub fn monomials(min: usize, max: usize) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    for d in min..=max {
        for p in (0..=d).rev() {
            out.push((p, d - p));
        }
    }
    out
}

#[derive(Debug, Clone, PartialEq)]
pub struct Sip {
    pub a: Poly,
    pub b: Poly,
    pub ap: Poly,
    pub bp: Poly,
}

/// Tangent-plane (TAN) projection with optional SIP distortion. Pixel coordinates are zero-based.
#[derive(Debug, Clone, PartialEq)]
pub struct Wcs {
    pub tan: TanWcs,
    pub sip: Option<Sip>,
}

impl From<TanWcs> for Wcs {
    fn from(tan: TanWcs) -> Self {
        Self { tan, sip: None }
    }
}

impl Wcs {
    pub fn crval(&self) -> SkyPos {
        self.tan.crval()
    }
    pub fn crpix(&self) -> PixelPos {
        self.tan.crpix()
    }
    pub fn cd(&self) -> [[f64; 2]; 2] {
        self.tan.cd()
    }
    pub fn pixel_scale_arcsec(&self) -> f64 {
        self.tan.pixel_scale_arcsec()
    }
    pub fn parity(&self) -> Parity {
        self.tan.parity()
    }

    pub fn sip_order(&self) -> usize {
        self.sip.as_ref().map_or(1, |s| s.a.order.max(s.b.order))
    }

    fn distort(&self, u: f64, v: f64) -> (f64, f64) {
        match &self.sip {
            None => (u, v),
            Some(s) => (u + s.a.eval(u, v), v + s.b.eval(u, v)),
        }
    }

    pub fn pixel_to_sky(&self, p: PixelPos) -> SkyPos {
        let c = self.crpix();
        let (u, v) = self.distort(p.x - c.x, p.y - c.y);
        let cd = self.cd();
        let xi = (cd[0][0] * u + cd[0][1] * v).to_radians();
        let eta = (cd[1][0] * u + cd[1][1] * v).to_radians();
        gnomonic_inverse(self.crval(), xi, eta)
    }

    /// Returns `None` outside the TAN projection hemisphere or if the SIP inverse fails to converge.
    pub fn sky_to_pixel(&self, s: SkyPos) -> Option<PixelPos> {
        self.vector_to_pixel(s.to_unit_vector())
    }

    pub fn vector_to_pixel(&self, u: UnitVector) -> Option<PixelPos> {
        let lin = self.tan.vector_to_pixel(u)?;
        let Some(sip) = &self.sip else {
            return Some(lin);
        };
        let c = self.crpix();
        let (up, vp) = (lin.x - c.x, lin.y - c.y);
        let order = sip.a.order.max(sip.b.order);
        let powers = |u: f64, v: f64| (Powers::of(u, order), Powers::of(v, order));
        let f = |u: f64, v: f64| {
            let (pu, pv) = powers(u, v);
            (
                u + sip.a.eval_at(&pu, &pv) - up,
                v + sip.b.eval_at(&pu, &pv) - vp,
            )
        };
        let norm = |(a, b): (f64, f64)| a.phypot(b);
        let (mut u, mut v) = (up + sip.ap.eval(up, vp), vp + sip.bp.eval(up, vp));
        if !(u.is_finite() && v.is_finite()) {
            (u, v) = (up, vp);
        }
        let mut r = f(u, v);
        let mut size = norm(r);
        for _ in 0..50 {
            if size < 1e-12 {
                break;
            }
            let (pu, pv) = powers(u, v);
            let (au, av) = sip.a.grad_at(&pu, &pv);
            let (bu, bv) = sip.b.grad_at(&pu, &pv);
            let (j11, j12, j21, j22) = (1.0 + au, av, bu, 1.0 + bv);
            let det = j11 * j22 - j12 * j21;
            if !det.is_finite() || det.abs() < 1e-300 {
                return None;
            }
            let du = (j22 * r.0 - j12 * r.1) / det;
            let dv = (-j21 * r.0 + j11 * r.1) / det;
            let mut step = 1.0;
            loop {
                let (nu, nv) = (u - step * du, v - step * dv);
                let nr = f(nu, nv);
                let next = norm(nr);
                if next < size || step < 1e-6 {
                    (u, v, r, size) = (nu, nv, nr, next);
                    break;
                }
                step *= 0.5;
            }
        }
        (u.is_finite() && v.is_finite() && size <= SIP_INVERSE_TOL_PX)
            .then(|| PixelPos::new(u + c.x, v + c.y))
    }

    pub fn fit_inverse(&mut self, width: f64, height: f64, order: usize) {
        let Some(sip) = &self.sip else { return };
        let c = self.crpix();
        let order = order.min(MAX_ORDER);
        let terms = monomials(1, order);
        let s_norm = (0.5 * width.max(height)).max(1.0);
        let k = terms.len();
        let mut normal = vec![0.0; k * k];
        let mut ap = vec![0.0; k];
        let mut bp = vec![0.0; k];
        let mut basis = vec![0.0; k];
        let n = 25;
        for i in 0..n {
            for j in 0..n {
                let x = -0.5 + width * i as f64 / (n - 1) as f64;
                let y = -0.5 + height * j as f64 / (n - 1) as f64;
                let (u, v) = (x - c.x, y - c.y);
                let (up, vp) = (u + sip.a.eval(u, v), v + sip.b.eval(u, v));
                let (un, vn) = (up / s_norm, vp / s_norm);
                for (f, &(p, q)) in basis.iter_mut().zip(&terms) {
                    *f = un.ppowi(p as i32) * vn.ppowi(q as i32);
                }
                for a in 0..k {
                    for b in 0..k {
                        normal[a * k + b] += basis[a] * basis[b];
                    }
                    ap[a] += basis[a] * (u - up);
                    bp[a] += basis[a] * (v - vp);
                }
            }
        }
        let Some(factor) = DenseFactor::new(normal, k) else {
            return;
        };
        let (Some(ap), Some(bp)) = (factor.solve(ap), factor.solve(bp)) else {
            return;
        };
        let to_poly = |c: Vec<f64>| Poly {
            order,
            terms: terms
                .iter()
                .zip(c)
                .map(|(&(p, q), v)| (p, q, v / s_norm.ppowi((p + q) as i32)))
                .collect(),
        };
        if let Some(s) = &mut self.sip {
            s.ap = to_poly(ap);
            s.bp = to_poly(bp);
        }
    }

    /// Serialize the WCS as FITS header cards with one-based reference pixels.
    pub fn to_cards(&self) -> Vec<String> {
        let mut cards = self.tan.to_cards();
        let Some(sip) = &self.sip else { return cards };
        for c in &mut cards {
            if c.starts_with("CTYPE1") {
                *c = crate::header::str_card("CTYPE1", "RA---TAN-SIP", "TAN + SIP distortion");
            } else if c.starts_with("CTYPE2") {
                *c = crate::header::str_card("CTYPE2", "DEC--TAN-SIP", "TAN + SIP distortion");
            }
        }
        for (name, poly) in [
            ("A", &sip.a),
            ("B", &sip.b),
            ("AP", &sip.ap),
            ("BP", &sip.bp),
        ] {
            cards.push(crate::header::num_card(
                &format!("{name}_ORDER"),
                &poly.order.to_string(),
            ));
            for &(p, q, v) in &poly.terms {
                cards.push(crate::header::num_card(
                    &format!("{name}_{p}_{q}"),
                    &crate::header::format_real(v),
                ));
            }
        }
        cards
    }

    pub fn from_header(h: &crate::Header) -> Result<Self, WcsError> {
        let sip = h.str_trimmed("CTYPE1")? == Some("RA---TAN-SIP")
            && h.str_trimmed("CTYPE2")? == Some("DEC--TAN-SIP");
        if !sip {
            return TanWcs::from_header(h).map(Wcs::from);
        }
        let mut lin = h.clone();
        lin.set("CTYPE1", crate::Value::Str("RA---TAN".into()));
        lin.set("CTYPE2", crate::Value::Str("DEC--TAN".into()));
        let tan = TanWcs::from_header(&lin)?;
        let read = |name: &str, min: usize| -> Result<Poly, WcsError> {
            let order = match h.number(&format!("{name}_ORDER"))? {
                Some(o) if (0.0..=MAX_ORDER as f64).contains(&o) && o.fract() == 0.0 => o as usize,
                Some(o) => {
                    return Err(WcsError::InvalidValue {
                        key: format!("{name}_ORDER"),
                        reason: format!("{o}"),
                    })
                }
                None => 0,
            };
            let mut terms = Vec::new();
            for (p, q) in monomials(min, order) {
                let key = format!("{name}_{p}_{q}");
                if let Some(c) = h.number(&key)? {
                    terms.push((p, q, c));
                }
            }
            Ok(Poly { order, terms })
        };
        let (a, b) = (read("A", 2)?, read("B", 2)?);
        if a.order < 2 || b.order < 2 {
            return Err(WcsError::MissingKey("A_ORDER/B_ORDER".into()));
        }
        Ok(Wcs {
            tan,
            sip: Some(Sip {
                a,
                b,
                ap: read("AP", 0)?,
                bp: read("BP", 0)?,
            }),
        })
    }
}

pub struct DenseFactor {
    matrix: Vec<f64>,
    pivots: Vec<usize>,
}

impl DenseFactor {
    pub fn new(mut matrix: Vec<f64>, n: usize) -> Option<Self> {
        if matrix.len() != n.checked_mul(n)? {
            return None;
        }
        let mut pivots = Vec::with_capacity(n);
        for col in 0..n {
            let pivot = (col..n).max_by(|&a, &b| {
                matrix[a * n + col]
                    .abs()
                    .total_cmp(&matrix[b * n + col].abs())
            })?;
            if matrix[pivot * n + col].abs() < 1e-300 {
                return None;
            }
            pivots.push(pivot);
            if pivot != col {
                for k in col..n {
                    matrix.swap(col * n + k, pivot * n + k);
                }
            }
            for row in col + 1..n {
                let f = matrix[row * n + col] / matrix[col * n + col];
                matrix[row * n + col] = f;
                if f != 0.0 {
                    for k in col + 1..n {
                        matrix[row * n + k] -= f * matrix[col * n + k];
                    }
                }
            }
        }
        Some(Self { matrix, pivots })
    }

    pub fn solve(&self, mut rhs: Vec<f64>) -> Option<Vec<f64>> {
        self.solve_in_place(&mut rhs)?;
        Some(rhs)
    }

    fn solve_in_place(&self, rhs: &mut [f64]) -> Option<()> {
        let n = self.pivots.len();
        if rhs.len() != n {
            return None;
        }
        for (col, &pivot) in self.pivots.iter().enumerate() {
            rhs.swap(col, pivot);
            for row in col + 1..n {
                let f = self.matrix[row * n + col];
                if f != 0.0 {
                    rhs[row] -= f * rhs[col];
                }
            }
        }
        for row in (0..n).rev() {
            let s: f64 = (row + 1..n)
                .map(|k| self.matrix[row * n + k] * rhs[k])
                .sum();
            rhs[row] = (rhs[row] - s) / self.matrix[row * n + row];
        }
        rhs.iter().all(|v| v.is_finite()).then_some(())
    }

    pub fn inverse(&self) -> Option<Vec<f64>> {
        let n = self.pivots.len();
        let mut inv = vec![0.0; n * n];
        let mut rhs = vec![0.0; n];
        for c in 0..n {
            rhs.fill(0.0);
            rhs[c] = 1.0;
            self.solve_in_place(&mut rhs)?;
            for r in 0..n {
                inv[r * n + c] = rhs[r];
            }
        }
        Some(inv)
    }
}

pub fn solve_dense(m: Vec<f64>, v: Vec<f64>) -> Option<Vec<f64>> {
    DenseFactor::new(m, v.len())?.solve(v)
}

pub fn invert_dense(m: &[f64], n: usize) -> Option<Vec<f64>> {
    DenseFactor::new(m.to_vec(), n)?.inverse()
}

pub fn plane(crval: SkyPos, s: SkyPos) -> Option<(f64, f64)> {
    gnomonic(crval, s)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Header;

    #[test]
    fn dense_factor_replays_successive_pivots_for_each_rhs() {
        let matrix = vec![1.0, 8.0, 3.0, 2.0, 1.0, 5.0, 4.0, 2.0, 6.0];
        let factor = DenseFactor::new(matrix.clone(), 3).unwrap();
        assert_eq!(factor.pivots, [2, 2, 2]);
        assert_eq!(
            factor.solve(vec![26.0, 19.0, 26.0]).unwrap(),
            [1.0, 2.0, 3.0]
        );
        assert_eq!(
            factor.solve(vec![1.0, 1.0, -2.0]).unwrap(),
            [-2.0, 0.0, 1.0]
        );

        let inverse = factor.inverse().unwrap();
        for (i, row) in matrix.chunks_exact(3).enumerate() {
            for j in 0..3 {
                let dot: f64 = row
                    .iter()
                    .zip(inverse[j..].iter().step_by(3))
                    .map(|(a, b)| a * b)
                    .sum();
                assert!((dot - f64::from(i == j)).abs() < 1e-14);
            }
        }
        assert!(DenseFactor::new(vec![1.0, 2.0, 2.0, 4.0], 2).is_none());
    }

    #[test]
    fn malformed_sip_terms_and_inverse_order_are_errors() {
        let cards = [
            "CTYPE1  = 'RA---TAN-SIP'",
            "CTYPE2  = 'DEC--TAN-SIP'",
            "CRVAL1  = 10.0",
            "CRVAL2  = 20.0",
            "CRPIX1  = 1.0",
            "CRPIX2  = 1.0",
            "CD1_1   = -1E-3",
            "CD2_2   = 1E-3",
            "A_ORDER = 2",
            "B_ORDER = 2",
        ];
        for (card, bad_key) in [
            ("A_2_0   = broken", "A_2_0"),
            ("AP_ORDER= 'broken'", "AP_ORDER"),
        ] {
            let header = Header::parse_cards(cards.into_iter().chain([card]));
            let result = Wcs::from_header(&header);
            assert!(matches!(result, Err(WcsError::InvalidValue { key, .. }) if key == bad_key));
        }
    }

    #[test]
    fn sip_inverse_converges_where_fixed_point_iteration_fails() {
        let tan = TanWcs::new(
            SkyPos::from_degrees(10.0, 20.0).unwrap(),
            PixelPos::new(0.0, 0.0),
            [[-2e-4, 0.0], [0.0, 2e-4]],
        )
        .unwrap();
        let empty = Poly {
            order: 0,
            terms: vec![],
        };
        let w = Wcs {
            tan,
            sip: Some(Sip {
                a: Poly {
                    order: 3,
                    terms: vec![(3, 0, 1e-6)],
                },
                b: empty.clone(),
                ap: empty.clone(),
                bp: empty,
            }),
        };
        let p = PixelPos::new(600.0, 5.0);
        let q = w.sky_to_pixel(w.pixel_to_sky(p)).expect("converges");
        assert!(
            (q.x - p.x).abs() < 1e-6 && (q.y - p.y).abs() < 1e-6,
            "{q:?}"
        );
    }

    fn distorted() -> Wcs {
        let tan = TanWcs::new(
            SkyPos::from_degrees(150.0, 2.0).unwrap(),
            PixelPos::new(511.5, 511.5),
            [[-8.3e-4, 1e-5], [1e-5, 8.3e-4]],
        )
        .unwrap();
        let k = 2e-8;
        let a = Poly {
            order: 3,
            terms: vec![(3, 0, k), (1, 2, k)],
        };
        let b = Poly {
            order: 3,
            terms: vec![(2, 1, k), (0, 3, k)],
        };
        let mut w = Wcs {
            tan,
            sip: Some(Sip {
                a,
                b,
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
        w.fit_inverse(1024.0, 1024.0, 5);
        w
    }

    #[test]
    fn polynomial_from_shared_powers_equals_the_direct_sum() {
        let polys = [
            Poly {
                order: 3,
                terms: vec![
                    (3, 0, 2e-8),
                    (1, 2, -3e-8),
                    (2, 0, 1e-6),
                    (0, 3, 7e-9),
                    (1, 1, 4e-7),
                ],
            },
            Poly {
                order: 1,
                terms: vec![(2, 2, 1e-12), (0, 11, 1e-30), (12, 1, -1e-33), (0, 0, 0.25)],
            },
            Poly {
                order: 0,
                terms: vec![],
            },
        ];
        for poly in &polys {
            for (u, v) in [(0.0, 0.0), (311.25, -487.5), (-0.003, 912.0), (1e3, 1e3)] {
                let direct: f64 = poly
                    .terms
                    .iter()
                    .map(|&(p, q, c)| c * u.ppowi(p as i32) * v.ppowi(q as i32))
                    .sum();
                assert_eq!(poly.eval(u, v).to_bits(), direct.to_bits());
                let (gu, gv) = poly.terms.iter().fold((0.0, 0.0), |(du, dv), &(p, q, c)| {
                    let a = if p > 0 {
                        c * p as f64 * u.ppowi(p as i32 - 1) * v.ppowi(q as i32)
                    } else {
                        0.0
                    };
                    let b = if q > 0 {
                        c * q as f64 * u.ppowi(p as i32) * v.ppowi(q as i32 - 1)
                    } else {
                        0.0
                    };
                    (du + a, dv + b)
                });
                let got = poly.grad(u, v);
                assert_eq!(
                    (got.0.to_bits(), got.1.to_bits()),
                    (gu.to_bits(), gv.to_bits())
                );
            }
        }
    }

    #[test]
    fn vector_to_pixel_equals_sky_to_pixel() {
        let w = distorted();
        for (x, y) in [(0.0, 0.0), (1023.0, 17.0), (300.25, 900.5), (511.5, 511.5)] {
            let s = w.pixel_to_sky(PixelPos::new(x, y));
            assert_eq!(w.vector_to_pixel(s.to_unit_vector()), w.sky_to_pixel(s));
        }
        let behind = SkyPos::from_degrees(330.0, -2.0).unwrap();
        assert_eq!(w.vector_to_pixel(behind.to_unit_vector()), None);
    }

    #[test]
    fn round_trip_with_distortion() {
        let w = distorted();
        for (x, y) in [(0.0, 0.0), (1023.0, 17.0), (300.25, 900.5), (511.5, 511.5)] {
            let p = PixelPos::new(x, y);
            let q = w.sky_to_pixel(w.pixel_to_sky(p)).unwrap();
            assert!(
                (q.x - x).abs() < 1e-8 && (q.y - y).abs() < 1e-8,
                "{p:?} -> {q:?}"
            );
        }
        let lin = w
            .tan
            .sky_to_pixel(w.pixel_to_sky(PixelPos::new(0.0, 0.0)))
            .unwrap();
        assert!(lin.x.phypot(lin.y) > 1.0);
        let sip = w.sip.as_ref().unwrap();
        let s = w.pixel_to_sky(PixelPos::new(0.0, 0.0));
        let l = w.tan.sky_to_pixel(s).unwrap();
        let (up, vp) = (l.x - 511.5, l.y - 511.5);
        let (u, v) = (up + sip.ap.eval(up, vp), vp + sip.bp.eval(up, vp));
        assert!(
            (u + 511.5).abs() < 1e-2 && (v + 511.5).abs() < 1e-2,
            "{u} {v}"
        );
    }

    #[test]
    fn inverse_fit_on_a_large_detector() {
        let tan = TanWcs::new(
            SkyPos::from_degrees(30.0, -60.0).unwrap(),
            PixelPos::new(2047.5, 2047.5),
            [[-2.8e-4, 0.0], [0.0, 2.8e-4]],
        )
        .unwrap();
        let k = 1.2e-9;
        let mut w = Wcs {
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
        w.fit_inverse(4096.0, 4096.0, 5);
        let sip = w.sip.as_ref().unwrap();
        let mut worst = 0.0f64;
        for (x, y) in [
            (0.0, 0.0),
            (4095.0, 4095.0),
            (0.0, 4095.0),
            (1000.0, 3000.0),
        ] {
            let s = w.pixel_to_sky(PixelPos::new(x, y));
            let l = w.tan.sky_to_pixel(s).unwrap();
            let (up, vp) = (l.x - 2047.5, l.y - 2047.5);
            let (u, v) = (up + sip.ap.eval(up, vp), vp + sip.bp.eval(up, vp));
            worst = worst.max((u + 2047.5 - x).phypot(v + 2047.5 - y));
        }
        assert!(worst < 0.01, "AP/BP inversion error {worst} px");
    }

    #[test]
    fn cards_round_trip() {
        let w = distorted();
        let cards = w.to_cards();
        assert!(cards.iter().all(|c| c.len() == 80));
        assert!(cards
            .iter()
            .any(|c| c.starts_with("CTYPE1  = 'RA---TAN-SIP'")));
        let back =
            Wcs::from_header(&Header::parse_cards(cards.iter().map(String::as_str))).unwrap();
        for (x, y) in [(0.0, 0.0), (1023.0, 1023.0)] {
            let p = PixelPos::new(x, y);
            assert!(back.pixel_to_sky(p).separation(w.pixel_to_sky(p)) < 1e-14);
        }
        let tan_cards = w.tan.to_cards();
        let t =
            Wcs::from_header(&Header::parse_cards(tan_cards.iter().map(String::as_str))).unwrap();
        assert!(t.sip.is_none());
    }

    #[test]
    fn monomial_order() {
        assert_eq!(
            monomials(2, 3),
            vec![(2, 0), (1, 1), (0, 2), (3, 0), (2, 1), (1, 2), (0, 3)]
        );
    }
}
