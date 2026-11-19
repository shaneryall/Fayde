use fayde_core::Portable;
use fayde_core::{Epoch, SkyPos, UnitVector};

pub const MAS: f64 = std::f64::consts::PI / (180.0 * 3600.0 * 1000.0);

#[derive(Debug, Clone, Copy)]
pub(crate) struct LocalFrame {
    pub(crate) u0: [f64; 3],
    east: [f64; 3],
    north: [f64; 3],
}

impl LocalFrame {
    pub(crate) fn at(pos: SkyPos) -> Self {
        let (sd, cd) = pos.dec.psin_cos();
        let (sr, cr) = pos.ra.psin_cos();
        Self {
            u0: [cd * cr, cd * sr, sd],
            east: [-sr, cr, 0.0],
            north: [-sd * cr, -sd * sr, cd],
        }
    }
}

pub fn propagate(pos: SkyPos, pm_mas_yr: [f64; 2], from: Epoch, to: Epoch) -> SkyPos {
    let dt = to.years_since(from);
    if dt == 0.0 || (pm_mas_yr[0] == 0.0 && pm_mas_yr[1] == 0.0) {
        return pos;
    }
    propagate_in(&LocalFrame::at(pos), pos, pm_mas_yr, from, to).0
}

pub(crate) fn propagate_in(
    frame: &LocalFrame,
    pos: SkyPos,
    pm_mas_yr: [f64; 2],
    from: Epoch,
    to: Epoch,
) -> (SkyPos, [f64; 3]) {
    let dt = to.years_since(from);
    if dt == 0.0 || (pm_mas_yr[0] == 0.0 && pm_mas_yr[1] == 0.0) {
        return (pos, frame.u0);
    }
    let (u0, east, north) = (frame.u0, frame.east, frame.north);
    let (ma, md) = (pm_mas_yr[0] * MAS * dt, pm_mas_yr[1] * MAS * dt);
    let v = [
        u0[0] + ma * east[0] + md * north[0],
        u0[1] + ma * east[1] + md * north[1],
        u0[2] + ma * east[2] + md * north[2],
    ];
    match UnitVector::normalize(v) {
        Some(u) => (u.to_sky(), u.0),
        None => (pos, u0),
    }
}

pub fn grow_error(sigma_c: f64, central: Epoch, sigma_pm: f64, to: Epoch) -> f64 {
    let dt = to.years_since(central);
    (sigma_c * sigma_c + dt * dt * sigma_pm * sigma_pm).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use fayde_core::to_arcsec;

    fn e(y: f64) -> Epoch {
        Epoch::julian_year(y).unwrap()
    }

    #[test]
    fn moves_by_pm_times_dt() {
        let p = SkyPos::from_degrees(10.0, 0.0).unwrap();
        let q = propagate(p, [1000.0, 0.0], e(2000.0), e(2010.0));
        let expect = to_arcsec((10_000.0 * MAS).patan());
        let got = to_arcsec(p.separation(q));
        assert!((got - expect).abs() < 1e-9, "{got} vs {expect}");
        assert!(q.ra > p.ra, "positive μα* moves east");
        let back = propagate(q, [1000.0, 0.0], e(2010.0), e(2000.0));
        assert!(to_arcsec(p.separation(back)) < 1e-3);
    }

    #[test]
    fn exact_at_the_pole() {
        let p = SkyPos::from_degrees(0.0, 90.0).unwrap();
        let q = propagate(p, [0.0, -500.0], e(2000.0), e(2020.0));
        let expect = to_arcsec((10_000.0 * MAS).patan());
        let got = to_arcsec(p.separation(q));
        assert!((got - expect).abs() < 1e-9, "{got} vs {expect}");
    }

    #[test]
    fn error_growth() {
        assert!(
            (grow_error(20.0, e(1995.0), 2.0, e(2025.0)) - (400.0f64 + 3600.0).sqrt()).abs()
                < 1e-12
        );
        assert_eq!(grow_error(15.0, e(2000.0), 5.0, e(2000.0)), 15.0);
    }
}
