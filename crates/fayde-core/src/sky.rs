//! Sky positions and unit direction vectors.

use crate::Portable;
use std::f64::consts::{PI, TAU};

/// A direction on the celestial sphere as ICRS right ascension and declination, in radians.
///
/// Right ascension is normalized to `[0, 2π)` and declination to `[-π/2, π/2]` by
/// [`SkyPos::new`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SkyPos {
    /// Right ascension, radians, `[0, 2π)`.
    pub ra: f64,
    /// Declination, radians, `[-π/2, π/2]`.
    pub dec: f64,
}

impl SkyPos {
    /// Build a position from radians, wrapping RA into `[0, 2π)`.
    ///
    /// Returns `None` for non-finite input or a declination outside `[-π/2, π/2]`.
    pub fn new(ra: f64, dec: f64) -> Option<Self> {
        if !ra.is_finite() || !dec.is_finite() || dec.abs() > PI / 2.0 {
            return None;
        }
        Some(Self {
            ra: wrap_two_pi(ra),
            dec,
        })
    }

    /// Build a position from degrees.
    pub fn from_degrees(ra_deg: f64, dec_deg: f64) -> Option<Self> {
        Self::new(ra_deg.to_radians(), dec_deg.to_radians())
    }

    /// `(ra, dec)` in degrees.
    pub fn to_degrees(self) -> (f64, f64) {
        (self.ra.to_degrees(), self.dec.to_degrees())
    }

    /// Unit vector in ICRS Cartesian axes.
    pub fn to_unit_vector(self) -> UnitVector {
        let (sd, cd) = self.dec.psin_cos();
        let (sr, cr) = self.ra.psin_cos();
        UnitVector([cd * cr, cd * sr, sd])
    }

    /// Angular separation in radians, accurate at all separations (Vincenty / atan2 form).
    pub fn separation(self, other: SkyPos) -> f64 {
        self.to_unit_vector().angle_to(other.to_unit_vector())
    }
}

/// A unit vector in ICRS Cartesian axes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UnitVector(pub [f64; 3]);

impl UnitVector {
    /// Normalize an arbitrary non-zero vector. Returns `None` for zero or non-finite input.
    pub fn normalize(v: [f64; 3]) -> Option<Self> {
        let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        if !n.is_finite() || n == 0.0 {
            return None;
        }
        Some(Self([v[0] / n, v[1] / n, v[2] / n]))
    }

    /// Dot product.
    pub fn dot(self, other: UnitVector) -> f64 {
        self.0[0] * other.0[0] + self.0[1] * other.0[1] + self.0[2] * other.0[2]
    }

    /// Cross product (not normalized).
    pub fn cross(self, other: UnitVector) -> [f64; 3] {
        let (a, b) = (self.0, other.0);
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    }

    /// Angle between two unit vectors, radians. Uses `atan2(|a×b|, a·b)` (the well-conditioned
    /// form, Vincenty 1975), which stays accurate for both tiny and near-antipodal separations.
    pub fn angle_to(self, other: UnitVector) -> f64 {
        let c = self.cross(other);
        (c[0] * c[0] + c[1] * c[1] + c[2] * c[2])
            .sqrt()
            .patan2(self.dot(other))
    }

    /// Back to RA/Dec.
    pub fn to_sky(self) -> SkyPos {
        let [x, y, z] = self.0;
        SkyPos {
            ra: wrap_two_pi(y.patan2(x)),
            dec: z.patan2((x * x + y * y).sqrt()),
        }
    }
}

/// Wrap an angle into `[0, 2π)`.
pub fn wrap_two_pi(angle: f64) -> f64 {
    let a = angle.rem_euclid(TAU);
    // rem_euclid can return exactly TAU for tiny negative inputs due to rounding.
    if a >= TAU {
        0.0
    } else {
        a
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ra_wraps_and_dec_is_validated() {
        let p = SkyPos::from_degrees(-10.0, 5.0).unwrap();
        assert!((p.to_degrees().0 - 350.0).abs() < 1e-12);
        assert!(SkyPos::from_degrees(0.0, 90.0001).is_none());
        assert!(SkyPos::new(f64::NAN, 0.0).is_none());
        assert_eq!(wrap_two_pi(-1e-20), 0.0);
    }

    #[test]
    fn unit_vector_round_trip() {
        for &(ra, dec) in &[(0.0, 0.0), (123.4, -56.7), (359.999, 89.9), (180.0, -90.0)] {
            let p = SkyPos::from_degrees(ra, dec).unwrap();
            let q = p.to_unit_vector().to_sky();
            assert!(p.separation(q) < 1e-14, "{ra} {dec}");
        }
    }

    #[test]
    fn separation_small_and_large() {
        let a = SkyPos::from_degrees(10.0, 20.0).unwrap();
        let b = SkyPos::from_degrees(10.0, 20.0 + 1.0 / 3600.0 / 1000.0).unwrap(); // 1 mas
        assert!((crate::to_arcsec(a.separation(b)) - 0.001).abs() < 1e-9);
        let c = SkyPos::from_degrees(190.0, -20.0).unwrap(); // antipode
        assert!((a.separation(c) - PI).abs() < 1e-12);
    }
}
