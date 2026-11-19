use fayde_core::Portable;
use fayde_core::{PixelPos, SkyPos, UnitVector};

use crate::WcsError;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Parity {
    Negative,
    Positive,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TanWcs {
    crval: SkyPos,
    crpix: PixelPos,
    cd: [[f64; 2]; 2],
    cd_inv: [[f64; 2]; 2],
    basis: Basis,
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct Basis {
    centre: [f64; 3],
    east: [f64; 3],
    north: [f64; 3],
}

impl TanWcs {
    pub fn new(crval: SkyPos, crpix: PixelPos, cd: [[f64; 2]; 2]) -> Result<Self, WcsError> {
        if !crval.ra.is_finite() || !crval.dec.is_finite() {
            return Err(WcsError::NonFinite("CRVAL"));
        }
        let crval = SkyPos::new(crval.ra, crval.dec).ok_or_else(|| WcsError::InvalidValue {
            key: "CRVAL2".into(),
            reason: "declination outside [-90, 90] degrees".into(),
        })?;
        if !crpix.x.is_finite() || !crpix.y.is_finite() {
            return Err(WcsError::NonFinite("CRPIX"));
        }
        if cd.iter().flatten().any(|v| !v.is_finite()) {
            return Err(WcsError::NonFinite("CD"));
        }
        let det = cd[0][0] * cd[1][1] - cd[0][1] * cd[1][0];
        if det == 0.0 || !det.is_finite() {
            return Err(WcsError::SingularMatrix);
        }
        let cd_inv = [
            [cd[1][1] / det, -cd[0][1] / det],
            [-cd[1][0] / det, cd[0][0] / det],
        ];
        if cd_inv.iter().flatten().any(|v| !v.is_finite()) {
            return Err(WcsError::SingularMatrix);
        }
        Ok(Self {
            crval,
            crpix,
            cd,
            cd_inv,
            basis: Basis::at(crval),
        })
    }

    pub fn crval(&self) -> SkyPos {
        self.crval
    }

    pub fn crpix(&self) -> PixelPos {
        self.crpix
    }

    pub fn cd(&self) -> [[f64; 2]; 2] {
        self.cd
    }

    pub fn cd_det(&self) -> f64 {
        self.cd[0][0] * self.cd[1][1] - self.cd[0][1] * self.cd[1][0]
    }

    pub fn pixel_scale_arcsec(&self) -> f64 {
        self.cd_det().abs().sqrt() * 3600.0
    }

    pub fn parity(&self) -> Parity {
        if self.cd_det() < 0.0 {
            Parity::Negative
        } else {
            Parity::Positive
        }
    }

    pub fn pixel_to_plane(&self, p: PixelPos) -> (f64, f64) {
        let dx = p.x - self.crpix.x;
        let dy = p.y - self.crpix.y;
        let xi = self.cd[0][0] * dx + self.cd[0][1] * dy;
        let eta = self.cd[1][0] * dx + self.cd[1][1] * dy;
        (xi.to_radians(), eta.to_radians())
    }

    pub fn plane_to_pixel(&self, xi: f64, eta: f64) -> PixelPos {
        let (xi, eta) = (xi.to_degrees(), eta.to_degrees());
        PixelPos::new(
            self.crpix.x + self.cd_inv[0][0] * xi + self.cd_inv[0][1] * eta,
            self.crpix.y + self.cd_inv[1][0] * xi + self.cd_inv[1][1] * eta,
        )
    }

    pub fn pixel_to_sky(&self, p: PixelPos) -> SkyPos {
        let (xi, eta) = self.pixel_to_plane(p);
        self.basis.deproject(xi, eta).to_sky()
    }

    pub fn sky_to_pixel(&self, s: SkyPos) -> Option<PixelPos> {
        self.vector_to_pixel(s.to_unit_vector())
    }

    pub fn vector_to_pixel(&self, u: UnitVector) -> Option<PixelPos> {
        let (xi, eta) = self.basis.project(u)?;
        Some(self.plane_to_pixel(xi, eta))
    }
}

pub fn gnomonic(centre: SkyPos, s: SkyPos) -> Option<(f64, f64)> {
    TangentPlane::at(centre).project(s.to_unit_vector())
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TangentPlane(Basis);

impl TangentPlane {
    pub fn at(centre: SkyPos) -> Self {
        Self(Basis::at(centre))
    }

    pub fn project(&self, u: UnitVector) -> Option<(f64, f64)> {
        self.0.project(u)
    }
}

pub fn gnomonic_inverse(centre: SkyPos, xi: f64, eta: f64) -> SkyPos {
    Basis::at(centre).deproject(xi, eta).to_sky()
}

impl Basis {
    fn at(c: SkyPos) -> Self {
        let (sd, cd) = c.dec.psin_cos();
        let (sr, cr) = c.ra.psin_cos();
        Self {
            centre: [cd * cr, cd * sr, sd],
            east: [-sr, cr, 0.0],
            north: [-sd * cr, -sd * sr, cd],
        }
    }

    fn project(&self, u: UnitVector) -> Option<(f64, f64)> {
        let w = dot(u.0, self.centre);
        if w <= 0.0 {
            return None;
        }
        Some((dot(u.0, self.east) / w, dot(u.0, self.north) / w))
    }

    fn deproject(&self, xi: f64, eta: f64) -> UnitVector {
        const FAR: f64 = 1e150;
        let far = xi.abs().max(eta.abs());
        let (w, xi, eta) = if far > FAR {
            let unit = |v: f64| {
                if v.is_infinite() {
                    v.signum()
                } else if far.is_infinite() {
                    0.0
                } else {
                    v / far
                }
            };
            (0.0, unit(xi), unit(eta))
        } else {
            (1.0, xi, eta)
        };
        let v = [
            w * self.centre[0] + xi * self.east[0] + eta * self.north[0],
            w * self.centre[1] + xi * self.east[1] + eta * self.north[1],
            w * self.centre[2] + xi * self.east[2] + eta * self.north[2],
        ];
        UnitVector::normalize(v).unwrap_or(UnitVector(self.centre))
    }
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[cfg(test)]
mod tests {
    use super::*;

    fn simple() -> TanWcs {
        TanWcs::new(
            SkyPos::from_degrees(150.0, 2.0).unwrap(),
            PixelPos::new(511.5, 511.5),
            [[-8.3e-4, 0.0], [0.0, 8.3e-4]],
        )
        .unwrap()
    }

    #[test]
    fn reference_pixel_maps_to_crval() {
        let w = simple();
        assert!(w.pixel_to_sky(w.crpix()).separation(w.crval()) < 1e-15);
    }

    #[test]
    fn east_is_left_for_negative_parity() {
        let w = simple();
        assert_eq!(w.parity(), Parity::Negative);
        let right = w.pixel_to_sky(PixelPos::new(611.5, 511.5));
        assert!(right.ra < w.crval().ra);
        assert!((w.pixel_scale_arcsec() - 2.988).abs() < 1e-9);
    }

    #[test]
    fn round_trip_and_behind_plane() {
        let w = simple();
        let p = PixelPos::new(3.25, 1020.0);
        let q = w.sky_to_pixel(w.pixel_to_sky(p)).unwrap();
        assert!((p.x - q.x).abs() < 1e-9 && (p.y - q.y).abs() < 1e-9);
        let behind = SkyPos::from_degrees(330.0, -2.0).unwrap();
        assert!(w.sky_to_pixel(behind).is_none());
    }

    #[test]
    fn vector_forms_equal_the_sky_forms() {
        let w = simple();
        for (ra, dec) in [
            (150.0, 2.0),
            (150.4, 1.7),
            (149.1, 2.9),
            (10.0, -80.0),
            (330.0, -2.0),
        ] {
            let s = SkyPos::from_degrees(ra, dec).unwrap();
            assert_eq!(w.vector_to_pixel(s.to_unit_vector()), w.sky_to_pixel(s));
            let plane = TangentPlane::at(w.crval());
            assert_eq!(plane.project(s.to_unit_vector()), gnomonic(w.crval(), s));
        }
    }

    #[test]
    fn rejects_singular_and_non_finite() {
        let c = SkyPos::from_degrees(0.0, 0.0).unwrap();
        let p = PixelPos::new(0.0, 0.0);
        assert_eq!(
            TanWcs::new(
                SkyPos {
                    ra: f64::NAN,
                    dec: 0.0
                },
                p,
                [[1.0, 0.0], [0.0, 1.0]]
            ),
            Err(WcsError::NonFinite("CRVAL"))
        );
        assert!(TanWcs::new(SkyPos { ra: 0.0, dec: 2.0 }, p, [[1.0, 0.0], [0.0, 1.0]]).is_err());
        assert_eq!(
            TanWcs::new(c, p, [[1e-310, 0.0], [0.0, 1.0]]),
            Err(WcsError::SingularMatrix)
        );
        assert_eq!(
            TanWcs::new(c, PixelPos::new(0.0, 0.0), [[1.0, 2.0], [2.0, 4.0]]),
            Err(WcsError::SingularMatrix)
        );
        assert_eq!(
            TanWcs::new(c, PixelPos::new(f64::NAN, 0.0), [[1.0, 0.0], [0.0, 1.0]]),
            Err(WcsError::NonFinite("CRPIX"))
        );
    }
}
