use crate::Portable;
use std::f64::consts::{PI, TAU};

/// Right ascension and declination in radians.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SkyPos {
    pub ra: f64,
    pub dec: f64,
}

impl SkyPos {
    pub fn new(ra: f64, dec: f64) -> Option<Self> {
        if !ra.is_finite() || !dec.is_finite() || dec.abs() > PI / 2.0 {
            return None;
        }
        Some(Self {
            ra: wrap_two_pi(ra),
            dec,
        })
    }

    pub fn from_degrees(ra_deg: f64, dec_deg: f64) -> Option<Self> {
        Self::new(ra_deg.to_radians(), dec_deg.to_radians())
    }

    pub fn to_degrees(self) -> (f64, f64) {
        (self.ra.to_degrees(), self.dec.to_degrees())
    }

    pub fn to_unit_vector(self) -> UnitVector {
        let (sd, cd) = self.dec.psin_cos();
        let (sr, cr) = self.ra.psin_cos();
        UnitVector([cd * cr, cd * sr, sd])
    }

    pub fn separation(self, other: SkyPos) -> f64 {
        self.to_unit_vector().angle_to(other.to_unit_vector())
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UnitVector(pub [f64; 3]);

impl UnitVector {
    pub fn normalize(v: [f64; 3]) -> Option<Self> {
        let squared = |v: [f64; 3]| v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
        let norm2 = squared(v);
        if norm2.is_normal() {
            let n = norm2.sqrt();
            return Some(Self(v.map(|x| x / n)));
        }
        if !v.iter().all(|x| x.is_finite()) {
            return None;
        }
        let scale = v.into_iter().map(f64::abs).fold(0.0, f64::max);
        if scale == 0.0 {
            return None;
        }
        let scaled = v.map(|x| x / scale);
        let n = squared(scaled).sqrt();
        Some(Self(scaled.map(|x| x / n)))
    }

    pub fn dot(self, other: UnitVector) -> f64 {
        self.0[0] * other.0[0] + self.0[1] * other.0[1] + self.0[2] * other.0[2]
    }

    pub fn cross(self, other: UnitVector) -> [f64; 3] {
        let (a, b) = (self.0, other.0);
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    }

    pub fn angle_to(self, other: UnitVector) -> f64 {
        let c = self.cross(other);
        (c[0] * c[0] + c[1] * c[1] + c[2] * c[2])
            .sqrt()
            .patan2(self.dot(other))
    }

    pub fn to_sky(self) -> SkyPos {
        let [x, y, z] = self.0;
        SkyPos {
            ra: wrap_two_pi(y.patan2(x)),
            dec: z.patan2((x * x + y * y).sqrt()),
        }
    }
}

pub fn wrap_two_pi(angle: f64) -> f64 {
    let a = angle.rem_euclid(TAU);
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
    fn finite_vectors_normalize_at_extreme_scales() {
        for value in [f64::MAX, f64::MIN_POSITIVE, f64::from_bits(1)] {
            let u = UnitVector::normalize([value, -value, value]).unwrap();
            assert!((u.dot(u) - 1.0).abs() < 1e-15);
            assert_eq!(u.0[0], -u.0[1]);
        }
        assert!(UnitVector::normalize([0.0; 3]).is_none());
        assert!(UnitVector::normalize([1.0, f64::NAN, 2.0]).is_none());
    }

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
        let b = SkyPos::from_degrees(10.0, 20.0 + 1.0 / 3600.0 / 1000.0).unwrap();
        assert!((crate::to_arcsec(a.separation(b)) - 0.001).abs() < 1e-9);
        let c = SkyPos::from_degrees(190.0, -20.0).unwrap();
        assert!((a.separation(c) - PI).abs() < 1e-12);
    }
}
