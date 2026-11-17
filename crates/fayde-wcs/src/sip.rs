//! SIP polynomial distortion (Shupe et al. 2005, ASP Conf. Ser. 347, 491) on top of TAN.
//!
//! Forward (pixel → sky): with `u = x − CRPIX1`, `v = y − CRPIX2`,
//! `u' = u + Σ A_pq u^p v^q`, `v' = v + Σ B_pq u^p v^q` (2 ≤ p+q ≤ order), then
//! `(ξ, η) = CD · (u', v')` and the gnomonic deprojection. The inverse polynomials `AP`, `BP`
//! (1 ≤ p+q ≤ order) are fitted for FITS output; [`Wcs::sky_to_pixel`] refines them by
//! damped Newton iteration on the forward polynomial so round trips are exact to ~1e-10 px,
//! and reports no pixel when the iteration does not converge.

use fayde_core::Portable;
use fayde_core::{PixelPos, SkyPos, UnitVector};

use crate::{gnomonic, gnomonic_inverse, Parity, TanWcs, WcsError};

/// One polynomial: coefficients for monomials `u^p v^q`.
#[derive(Debug, Clone, PartialEq)]
pub struct Poly {
    /// Highest total degree.
    pub order: usize,
    /// `(p, q, coefficient)`.
    pub terms: Vec<(usize, usize, f64)>,
}

impl Poly {
    /// Evaluate at `(u, v)`.
    pub fn eval(&self, u: f64, v: f64) -> f64 {
        self.eval_at(&Powers::of(u, self.order), &Powers::of(v, self.order))
    }

    /// Partial derivatives `(∂/∂u, ∂/∂v)` at `(u, v)`.
    pub fn grad(&self, u: f64, v: f64) -> (f64, f64) {
        self.grad_at(&Powers::of(u, self.order), &Powers::of(v, self.order))
    }

    /// [`Poly::eval`] from the powers of `u` and `v`.
    fn eval_at(&self, u: &Powers, v: &Powers) -> f64 {
        self.terms
            .iter()
            .map(|&(p, q, c)| c * u.get(p) * v.get(q))
            .sum()
    }

    /// [`Poly::grad`] from the powers of `u` and `v`.
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

/// Highest SIP order: the convention's limit.
const MAX_ORDER: usize = 9;

/// The powers `x⁰ … xⁿ` of one number, each from [`Portable::ppowi`]. A polynomial of order n
/// has (n+1)(n+2)/2 terms, each a product of two powers, and the inverse iteration evaluates
/// two polynomials and their gradients at the same point: the same few powers, many times.
struct Powers {
    x: f64,
    /// `table[n]` is `x.ppowi(n)` for `n < known`.
    table: [f64; MAX_ORDER + 1],
    known: usize,
}

impl Powers {
    /// The powers of `x` up to `order` (those above the convention's limit are worked out when
    /// asked for).
    fn of(x: f64, order: usize) -> Self {
        let known = order.min(MAX_ORDER) + 1;
        let mut table = [0.0; MAX_ORDER + 1];
        for (n, t) in table.iter_mut().enumerate().take(known) {
            *t = x.ppowi(n as i32);
        }
        Self { x, table, known }
    }

    /// `x.ppowi(n)`.
    fn get(&self, n: usize) -> f64 {
        if n < self.known {
            self.table[n]
        } else {
            self.x.ppowi(n as i32)
        }
    }
}

/// Largest residual |forward(u, v) − target|, pixels, accepted from the SIP inversion.
const SIP_INVERSE_TOL_PX: f64 = 1e-8;

/// Monomials `(p, q)` with `min ≤ p+q ≤ max`, in FITS-conventional order (by degree, then p
/// descending).
pub fn monomials(min: usize, max: usize) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    for d in min..=max {
        for p in (0..=d).rev() {
            out.push((p, d - p));
        }
    }
    out
}

/// SIP distortion polynomials.
#[derive(Debug, Clone, PartialEq)]
pub struct Sip {
    /// Forward, x.
    pub a: Poly,
    /// Forward, y.
    pub b: Poly,
    /// Inverse, x (for FITS output and as a starting guess).
    pub ap: Poly,
    /// Inverse, y.
    pub bp: Poly,
}

/// A TAN WCS with optional SIP distortion.
#[derive(Debug, Clone, PartialEq)]
pub struct Wcs {
    /// Linear part.
    pub tan: TanWcs,
    /// Distortion, if any.
    pub sip: Option<Sip>,
}

impl From<TanWcs> for Wcs {
    fn from(tan: TanWcs) -> Self {
        Self { tan, sip: None }
    }
}

impl Wcs {
    /// Reference sky position.
    pub fn crval(&self) -> SkyPos {
        self.tan.crval()
    }
    /// Reference pixel (zero-based).
    pub fn crpix(&self) -> PixelPos {
        self.tan.crpix()
    }
    /// CD matrix, degrees per pixel.
    pub fn cd(&self) -> [[f64; 2]; 2] {
        self.tan.cd()
    }
    /// Pixel scale at the reference point, arcsec/px.
    pub fn pixel_scale_arcsec(&self) -> f64 {
        self.tan.pixel_scale_arcsec()
    }
    /// Parity of the linear part.
    pub fn parity(&self) -> Parity {
        self.tan.parity()
    }

    /// SIP order (1 when there is no distortion).
    pub fn sip_order(&self) -> usize {
        self.sip.as_ref().map_or(1, |s| s.a.order)
    }

    /// Distorted offsets `(u', v')` from raw offsets `(u, v)`.
    fn distort(&self, u: f64, v: f64) -> (f64, f64) {
        match &self.sip {
            None => (u, v),
            Some(s) => (u + s.a.eval(u, v), v + s.b.eval(u, v)),
        }
    }

    /// Sky position of a pixel.
    pub fn pixel_to_sky(&self, p: PixelPos) -> SkyPos {
        let c = self.crpix();
        let (u, v) = self.distort(p.x - c.x, p.y - c.y);
        let cd = self.cd();
        let xi = (cd[0][0] * u + cd[0][1] * v).to_radians();
        let eta = (cd[1][0] * u + cd[1][1] * v).to_radians();
        gnomonic_inverse(self.crval(), xi, eta)
    }

    /// Pixel of a sky position, or `None` behind the tangent plane.
    pub fn sky_to_pixel(&self, s: SkyPos) -> Option<PixelPos> {
        self.vector_to_pixel(s.to_unit_vector())
    }

    /// [`Wcs::sky_to_pixel`] for a direction given as a unit vector: the same result to the
    /// last bit for `s.to_unit_vector()` (see [`TanWcs::vector_to_pixel`]).
    pub fn vector_to_pixel(&self, u: UnitVector) -> Option<PixelPos> {
        let lin = self.tan.vector_to_pixel(u)?;
        let Some(sip) = &self.sip else {
            return Some(lin);
        };
        let c = self.crpix();
        let (up, vp) = (lin.x - c.x, lin.y - c.y);
        // Solve F(u, v) = (u + A − u', v + B − v') = 0 by Newton's method from the fitted
        // inverse, halving steps that do not reduce |F| (a fixed-point iteration converges
        // only while |∂A| < 1 and can settle on a wrong answer).
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
        // |r|, kept beside r: it is compared at every step of the line search.
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

    /// Fit `AP`/`BP` of order `order` on a grid covering `width × height` pixels, so FITS readers
    /// without iteration invert the forward polynomials accurately.
    pub fn fit_inverse(&mut self, width: f64, height: f64, order: usize) {
        let Some(sip) = &self.sip else { return };
        let c = self.crpix();
        let terms = monomials(1, order);
        let s_norm = (0.5 * width.max(height)).max(1.0);
        let mut rows = Vec::new();
        let n = 25;
        for i in 0..n {
            for j in 0..n {
                let x = -0.5 + width * i as f64 / (n - 1) as f64;
                let y = -0.5 + height * j as f64 / (n - 1) as f64;
                let (u, v) = (x - c.x, y - c.y);
                let (up, vp) = (u + sip.a.eval(u, v), v + sip.b.eval(u, v));
                // Normalized inputs for conditioning (rescaled after the solve).
                rows.push((up / s_norm, vp / s_norm, u - up, v - vp));
            }
        }
        let fit = |pick: fn(&(f64, f64, f64, f64)) -> f64| -> Option<Vec<f64>> {
            let k = terms.len();
            let mut m = vec![0.0; k * k];
            let mut rhs = vec![0.0; k];
            for r in &rows {
                let f: Vec<f64> = terms
                    .iter()
                    .map(|&(p, q)| r.0.ppowi(p as i32) * r.1.ppowi(q as i32))
                    .collect();
                for a in 0..k {
                    for b in 0..k {
                        m[a * k + b] += f[a] * f[b];
                    }
                    rhs[a] += f[a] * pick(r);
                }
            }
            solve_dense(m, rhs)
        };
        let (Some(ap), Some(bp)) = (fit(|r| r.2), fit(|r| r.3)) else {
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

    /// FITS cards: TAN cards, with `-SIP` and the `A`, `B`, `AP`, `BP` polynomials when present.
    pub fn to_cards(&self) -> Vec<String> {
        let mut cards = self.tan.to_cards();
        let Some(sip) = &self.sip else { return cards };
        for c in &mut cards {
            if c.starts_with("CTYPE1") {
                *c = crate::header::str_card_pub("CTYPE1", "RA---TAN-SIP", "TAN + SIP distortion");
            } else if c.starts_with("CTYPE2") {
                *c = crate::header::str_card_pub("CTYPE2", "DEC--TAN-SIP", "TAN + SIP distortion");
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

    /// Read TAN or TAN-SIP from header values.
    pub fn from_header(h: &crate::Header) -> Result<Self, WcsError> {
        let c1 = h
            .get("CTYPE1")
            .and_then(|v| v.as_str())
            .map(str::trim)
            .unwrap_or("");
        let c2 = h
            .get("CTYPE2")
            .and_then(|v| v.as_str())
            .map(str::trim)
            .unwrap_or("");
        let sip = c1 == "RA---TAN-SIP" && c2 == "DEC--TAN-SIP";
        if !sip {
            return TanWcs::from_header(h).map(Wcs::from);
        }
        let mut lin = h.clone();
        lin.set("CTYPE1", crate::Value::Str("RA---TAN".into()));
        lin.set("CTYPE2", crate::Value::Str("DEC--TAN".into()));
        let tan = TanWcs::from_header(&lin)?;
        let read = |name: &str, min: usize| -> Result<Poly, WcsError> {
            let order = match h.get(&format!("{name}_ORDER")).and_then(|v| v.as_f64()) {
                Some(o) if (0.0..=9.0).contains(&o) && o.fract() == 0.0 => o as usize,
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
                if let Some(v) = h.get(&key) {
                    // A coefficient that is not a finite number is an error, not an absent term.
                    let c = v.as_f64().ok_or_else(|| WcsError::InvalidValue {
                        key: key.clone(),
                        reason: format!("expected a number, got {v:?}"),
                    })?;
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

/// Solve a dense symmetric system (Gaussian elimination with partial pivoting).
pub fn solve_dense(mut m: Vec<f64>, mut v: Vec<f64>) -> Option<Vec<f64>> {
    let n = v.len();
    for col in 0..n {
        let piv =
            (col..n).max_by(|&a, &b| m[a * n + col].abs().total_cmp(&m[b * n + col].abs()))?;
        if m[piv * n + col].abs() < 1e-300 {
            return None;
        }
        if piv != col {
            for k in 0..n {
                m.swap(col * n + k, piv * n + k);
            }
            v.swap(col, piv);
        }
        for row in col + 1..n {
            let f = m[row * n + col] / m[col * n + col];
            if f == 0.0 {
                continue;
            }
            for k in col..n {
                m[row * n + k] -= f * m[col * n + k];
            }
            v[row] -= f * v[col];
        }
    }
    let mut x = vec![0.0; n];
    for row in (0..n).rev() {
        let s: f64 = (row + 1..n).map(|k| m[row * n + k] * x[k]).sum();
        x[row] = (v[row] - s) / m[row * n + row];
    }
    x.iter().all(|v| v.is_finite()).then_some(x)
}

/// Inverse of a dense matrix (n×n, row-major), column by column.
pub fn invert_dense(m: &[f64], n: usize) -> Option<Vec<f64>> {
    let mut inv = vec![0.0; n * n];
    for c in 0..n {
        let mut e = vec![0.0; n];
        e[c] = 1.0;
        let col = solve_dense(m.to_vec(), e)?;
        for r in 0..n {
            inv[r * n + c] = col[r];
        }
    }
    Some(inv)
}

/// Gnomonic projection re-exported for fitting code.
pub fn plane(crval: SkyPos, s: SkyPos) -> Option<(f64, f64)> {
    gnomonic(crval, s)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Header;

    #[test]
    fn sip_inverse_converges_where_fixed_point_iteration_fails() {
        // u' = u + 1e-6 u³: monotonic, but |∂A/∂u| = 3e-6 u² > 1 beyond u ≈ 577, where a
        // fixed-point iteration diverges or settles wrong.
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
        // Barrel-like: u' = u(1 + k r²) with k = 2e-8 → ~1.4 px at the corners.
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

    /// A polynomial evaluated from shared powers is the term-by-term sum, to the last bit:
    /// also when a term's degree is above the polynomial's stated order or the convention's.
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
            // Stated order below its terms, and a term beyond order 9.
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

    /// Projecting a unit vector through a distorted WCS is projecting its sky position.
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
        // Distortion is real: TAN alone disagrees by > 1 px at a corner.
        let lin = w
            .tan
            .sky_to_pixel(w.pixel_to_sky(PixelPos::new(0.0, 0.0)))
            .unwrap();
        assert!(lin.x.phypot(lin.y) > 1.0);
        // The fitted AP/BP alone (order = forward + 2) invert to < 0.01 px: what FITS readers
        // without iteration use.
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
        // Order-5 inverse on 4096 px: raw u^10 ≈ 1e33 in the normal matrix.
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
        // Plain TAN headers still read as TAN.
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
