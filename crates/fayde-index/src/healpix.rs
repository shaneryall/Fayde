use fayde_core::Portable;
use std::f64::consts::{FRAC_PI_2, TAU};

use fayde_core::UnitVector;

const JRLL: [i64; 12] = [2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4];
const JPLL: [i64; 12] = [1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Healpix {
    nside: u32,
}

impl Healpix {
    pub fn new(nside: u32) -> Option<Self> {
        (nside.is_power_of_two() && nside <= 1 << 20).then_some(Self { nside })
    }

    pub fn for_pixel_side(max_side: f64) -> Self {
        let mut nside = 1u32;
        while nside < 1 << 20 && pixel_side(nside) > max_side {
            nside *= 2;
        }
        Self { nside }
    }

    pub fn nside(self) -> u32 {
        self.nside
    }

    pub fn npix(self) -> u64 {
        12 * u64::from(self.nside) * u64::from(self.nside)
    }

    pub fn pixel_side(self) -> f64 {
        pixel_side(self.nside)
    }

    pub fn pixel(self, v: UnitVector) -> u64 {
        let [x, y, z] = v.0;
        let phi = y.patan2(x);
        self.pixel_z_phi(z.clamp(-1.0, 1.0), phi)
    }

    fn pixel_z_phi(self, z: f64, phi: f64) -> u64 {
        let ns = i64::from(self.nside);
        let nsf = f64::from(self.nside);
        let za = z.abs();
        let tt = phi.rem_euclid(TAU) / FRAC_PI_2;
        let (face, ix, iy) = if za <= 2.0 / 3.0 {
            let t1 = nsf * (0.5 + tt);
            let t2 = nsf * z * 0.75;
            let jp = (t1 - t2) as i64;
            let jm = (t1 + t2) as i64;
            let ifp = jp / ns;
            let ifm = jm / ns;
            let face = if ifp == ifm {
                ifp | 4
            } else if ifp < ifm {
                ifp
            } else {
                ifm + 8
            };
            (face, jm & (ns - 1), ns - (jp & (ns - 1)) - 1)
        } else {
            let ntt = (tt as i64).min(3);
            let tp = tt - ntt as f64;
            let tmp = nsf * (3.0 * (1.0 - za)).sqrt();
            let jp = ((tp * tmp) as i64).min(ns - 1);
            let jm = (((1.0 - tp) * tmp) as i64).min(ns - 1);
            if z >= 0.0 {
                (ntt, ns - jm - 1, ns - jp - 1)
            } else {
                (ntt + 8, jp, jm)
            }
        };
        (face as u64) * (ns as u64) * (ns as u64) + spread(ix as u64) + (spread(iy as u64) << 1)
    }

    pub fn centre(self, pix: u64) -> UnitVector {
        let ns = i64::from(self.nside);
        let nsf = f64::from(self.nside);
        let npface = (ns * ns) as u64;
        let face = (pix / npface) as usize;
        let ipf = pix % npface;
        let ix = compress(ipf) as i64;
        let iy = compress(ipf >> 1) as i64;
        let jr = JRLL[face] * ns - ix - iy - 1;
        let (nr, z, kshift) = if jr < ns {
            let nr = jr;
            (nr, 1.0 - (nr * nr) as f64 / (3.0 * nsf * nsf), 0)
        } else if jr > 3 * ns {
            let nr = 4 * ns - jr;
            (nr, (nr * nr) as f64 / (3.0 * nsf * nsf) - 1.0, 0)
        } else {
            (ns, (2 * ns - jr) as f64 * 2.0 / (3.0 * nsf), (jr - ns) & 1)
        };
        let mut jp = (JPLL[face] * nr + ix - iy + 1 + kshift) / 2;
        if jp > 4 * ns {
            jp -= 4 * ns;
        }
        if jp < 1 {
            jp += 4 * ns;
        }
        let phi = (jp as f64 - (kshift + 1) as f64 * 0.5) * (FRAC_PI_2 / nr as f64);
        let s = (1.0 - z * z).max(0.0).sqrt();
        UnitVector([s * phi.pcos(), s * phi.psin(), z])
    }
}

fn pixel_side(nside: u32) -> f64 {
    (4.0 * std::f64::consts::PI / (12.0 * f64::from(nside) * f64::from(nside))).sqrt()
}

fn spread(mut v: u64) -> u64 {
    v &= 0xFFFF_FFFF;
    v = (v | (v << 16)) & 0x0000_FFFF_0000_FFFF;
    v = (v | (v << 8)) & 0x00FF_00FF_00FF_00FF;
    v = (v | (v << 4)) & 0x0F0F_0F0F_0F0F_0F0F;
    v = (v | (v << 2)) & 0x3333_3333_3333_3333;
    (v | (v << 1)) & 0x5555_5555_5555_5555
}

fn compress(mut v: u64) -> u64 {
    v &= 0x5555_5555_5555_5555;
    v = (v | (v >> 1)) & 0x3333_3333_3333_3333;
    v = (v | (v >> 2)) & 0x0F0F_0F0F_0F0F_0F0F;
    v = (v | (v >> 4)) & 0x00FF_00FF_00FF_00FF;
    v = (v | (v >> 8)) & 0x0000_FFFF_0000_FFFF;
    (v | (v >> 16)) & 0x0000_0000_FFFF_FFFF
}

#[cfg(test)]
mod tests {
    use super::*;
    use fayde_core::SkyPos;

    fn v(ra: f64, dec: f64) -> UnitVector {
        SkyPos::from_degrees(ra, dec).unwrap().to_unit_vector()
    }

    fn points(n: usize) -> impl Iterator<Item = UnitVector> {
        let golden = std::f64::consts::PI * (3.0 - 5f64.sqrt());
        (0..n).map(move |i| {
            let z = 1.0 - 2.0 * (i as f64 + 0.5) / n as f64;
            let s = (1.0 - z * z).sqrt();
            let phi = golden * i as f64;
            UnitVector([s * phi.pcos(), s * phi.psin(), z])
        })
    }

    #[test]
    fn base_faces() {
        let h = Healpix::new(1).unwrap();
        assert_eq!(h.pixel(v(10.0, 80.0)), 0);
        assert_eq!(h.pixel(v(100.0, 80.0)), 1);
        assert_eq!(h.pixel(v(0.0, 0.0)), 4);
        assert_eq!(h.pixel(v(90.0, 0.0)), 5);
        assert_eq!(h.pixel(v(10.0, -80.0)), 8);
        assert_eq!(h.pixel(v(280.0, -80.0)), 11);
        assert!(
            h.centre(4)
                .to_sky()
                .separation(SkyPos::from_degrees(0.0, 0.0).unwrap())
                < 1e-12
        );
    }

    #[test]
    fn nested_hierarchy() {
        for nside in [2u32, 16, 256, 4096] {
            let fine = Healpix::new(nside).unwrap();
            let coarse = Healpix::new(nside / 2).unwrap();
            for p in points(5000) {
                assert_eq!(fine.pixel(p) >> 2, coarse.pixel(p), "nside {nside}");
            }
        }
    }

    #[test]
    fn centres_map_to_their_pixel() {
        for nside in [1u32, 4, 64] {
            let h = Healpix::new(nside).unwrap();
            for pix in 0..h.npix() {
                assert_eq!(h.pixel(h.centre(pix)), pix, "nside {nside} pix {pix}");
            }
        }
        let h = Healpix::new(1 << 14).unwrap();
        for pix in (0..h.npix()).step_by(9_999_991) {
            assert_eq!(h.pixel(h.centre(pix)), pix);
        }
    }

    #[test]
    fn equal_area() {
        let h = Healpix::new(8).unwrap();
        let n = 768 * 400;
        let mut counts = vec![0u32; h.npix() as usize];
        for p in points(n) {
            counts[h.pixel(p) as usize] += 1;
        }
        assert!(
            counts.iter().all(|&c| (380..=420).contains(&c)),
            "{:?}",
            counts.iter().min()
        );
    }

    #[test]
    fn nside_for_side() {
        let h = Healpix::for_pixel_side(0.1f64.to_radians());
        assert!(h.pixel_side() <= 0.1f64.to_radians());
        assert!(Healpix::new(h.nside() / 2).unwrap().pixel_side() > 0.1f64.to_radians());
        assert!(Healpix::new(3).is_none());
    }
}
