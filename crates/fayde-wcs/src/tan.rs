//! Gnomonic (TAN) projection with the default `LONPOLE = 180°`.
//!
//! Source: Greisen & Calabretta 2002 (A&A 395, 1061) and Calabretta & Greisen 2002
//! (A&A 395, 1077), FITS WCS Papers I and II. See `docs/references.md`.
//!
//! The projection is computed in vector form: the tangent plane is spanned by the local east and
//! north unit vectors at the reference point. This stays exact near the celestial poles, where
//! the usual RA/Dec trigonometric formulas lose precision.

use fayde_core::Portable;
use fayde_core::{PixelPos, SkyPos, UnitVector};

use crate::WcsError;

/// Sign of the CD matrix determinant.
///
/// `Negative` (det < 0) is an un-mirrored view of the sky: with north up, east is to the left.
/// `Positive` is the mirror image, as produced by an odd number of reflections in the optics or
/// a flipped readout.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Parity {
    /// det(CD) < 0: sky as seen, north up / east left.
    Negative,
    /// det(CD) > 0: mirrored.
    Positive,
}

/// A TAN (gnomonic) world coordinate system without distortion.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TanWcs {
    crval: SkyPos,
    crpix: PixelPos,
    /// Degrees per pixel. `cd[i][j]` = `CD(i+1)_(j+1)`.
    cd: [[f64; 2]; 2],
    cd_inv: [[f64; 2]; 2],
    basis: Basis,
}

/// Reference direction and the east/north tangent-plane axes at it.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Basis {
    centre: [f64; 3],
    east: [f64; 3],
    north: [f64; 3],
}

impl TanWcs {
    /// Build a TAN WCS.
    ///
    /// - `crval`: sky position of the reference pixel.
    /// - `crpix`: reference pixel in Fayde's zero-based convention (FITS `CRPIX` − 1).
    /// - `cd`: linear transform in degrees per pixel, `cd[i][j]` = FITS `CD(i+1)_(j+1)`.
    pub fn new(crval: SkyPos, crpix: PixelPos, cd: [[f64; 2]; 2]) -> Result<Self, WcsError> {
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
        Ok(Self {
            crval,
            crpix,
            cd,
            cd_inv,
            basis: Basis::at(crval),
        })
    }

    /// Reference sky position (`CRVAL`).
    pub fn crval(&self) -> SkyPos {
        self.crval
    }

    /// Reference pixel, zero-based.
    pub fn crpix(&self) -> PixelPos {
        self.crpix
    }

    /// CD matrix, degrees per pixel.
    pub fn cd(&self) -> [[f64; 2]; 2] {
        self.cd
    }

    /// Determinant of the CD matrix, square degrees per pixel².
    pub fn cd_det(&self) -> f64 {
        self.cd[0][0] * self.cd[1][1] - self.cd[0][1] * self.cd[1][0]
    }

    /// Mean pixel scale at the reference point, arcseconds per pixel: `sqrt(|det CD|)`.
    pub fn pixel_scale_arcsec(&self) -> f64 {
        self.cd_det().abs().sqrt() * 3600.0
    }

    /// Image parity.
    pub fn parity(&self) -> Parity {
        if self.cd_det() < 0.0 {
            Parity::Negative
        } else {
            Parity::Positive
        }
    }

    /// Intermediate world coordinates (ξ, η) in **radians** for a pixel.
    pub fn pixel_to_plane(&self, p: PixelPos) -> (f64, f64) {
        let dx = p.x - self.crpix.x;
        let dy = p.y - self.crpix.y;
        let xi = self.cd[0][0] * dx + self.cd[0][1] * dy;
        let eta = self.cd[1][0] * dx + self.cd[1][1] * dy;
        (xi.to_radians(), eta.to_radians())
    }

    /// Pixel for intermediate world coordinates (ξ, η) in radians.
    pub fn plane_to_pixel(&self, xi: f64, eta: f64) -> PixelPos {
        let (xi, eta) = (xi.to_degrees(), eta.to_degrees());
        PixelPos::new(
            self.crpix.x + self.cd_inv[0][0] * xi + self.cd_inv[0][1] * eta,
            self.crpix.y + self.cd_inv[1][0] * xi + self.cd_inv[1][1] * eta,
        )
    }

    /// Sky position of a pixel. Never panics: a pixel so far off the image that its offset
    /// overflows maps to the point 90° from the reference point in that direction, and a NaN
    /// coordinate maps to the reference point.
    pub fn pixel_to_sky(&self, p: PixelPos) -> SkyPos {
        let (xi, eta) = self.pixel_to_plane(p);
        self.basis.deproject(xi, eta).to_sky()
    }

    /// Pixel of a sky position, or `None` if it is 90° or more from the reference point
    /// (behind the tangent plane).
    pub fn sky_to_pixel(&self, s: SkyPos) -> Option<PixelPos> {
        self.vector_to_pixel(s.to_unit_vector())
    }

    /// [`TanWcs::sky_to_pixel`] for a direction given as a unit vector. For
    /// `s.to_unit_vector()` the result is that of `sky_to_pixel(s)` to the last bit, without
    /// the sines and cosines: use it for positions projected many times.
    pub fn vector_to_pixel(&self, u: UnitVector) -> Option<PixelPos> {
        let (xi, eta) = self.basis.project(u)?;
        Some(self.plane_to_pixel(xi, eta))
    }
}

/// Gnomonic projection of `s` onto the tangent plane at `centre`: (ξ, η) in radians, ξ east,
/// η north. `None` if `s` is 90° or more from `centre`.
pub fn gnomonic(centre: SkyPos, s: SkyPos) -> Option<(f64, f64)> {
    TangentPlane::at(centre).project(s.to_unit_vector())
}

/// The tangent plane at one point, for projecting many directions onto it: [`gnomonic`] without
/// working out the plane's axes again for every one.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TangentPlane(Basis);

impl TangentPlane {
    /// The plane tangent to the sphere at `centre`.
    pub fn at(centre: SkyPos) -> Self {
        Self(Basis::at(centre))
    }

    /// (ξ, η) in radians of the direction `u`, ξ east, η north; `None` if it is 90° or more
    /// from the tangent point. For `s.to_unit_vector()` this is `gnomonic(centre, s)` to the
    /// last bit.
    pub fn project(&self, u: UnitVector) -> Option<(f64, f64)> {
        self.0.project(u)
    }
}

/// Inverse of [`gnomonic`].
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

    /// Unit vector → tangent-plane (ξ, η) in radians.
    fn project(&self, u: UnitVector) -> Option<(f64, f64)> {
        let w = dot(u.0, self.centre);
        if w <= 0.0 {
            return None;
        }
        Some((dot(u.0, self.east) / w, dot(u.0, self.north) / w))
    }

    /// Tangent-plane (ξ, η) in radians → unit vector. Defined for every input: an offset too
    /// large to square gives the direction it tends to, 90° from the centre, and a NaN offset
    /// gives the centre.
    fn deproject(&self, xi: f64, eta: f64) -> UnitVector {
        // Past this the squared norm overflows, and the centre term is long since negligible.
        const FAR: f64 = 1e150;
        let far = xi.abs().max(eta.abs());
        let (w, xi, eta) = if far > FAR {
            // Only the direction is left: an infinite offset keeps its sign, a finite one its
            // share of the larger.
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
        // For finite offsets |v| >= 1 (centre is orthogonal to east and north); only NaN fails.
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
        // Moving +x with CD1_1 < 0 decreases RA (west).
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

    /// Projecting a unit vector is projecting its sky position, to the last bit, and the
    /// tangent plane is the plane of `gnomonic`.
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
