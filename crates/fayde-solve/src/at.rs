use fayde_catalog::{Catalog, Star, StarQuery};
use fayde_core::{to_arcsec, Epoch, PixelPos, Portable, SkyPos, UnitVector};
use fayde_wcs::Wcs;

use crate::SolveError;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AtPosition {
    Pixel(PixelPos, Option<f64>),
    Sky(SkyPos, Option<f64>),
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AtOptions {
    pub epoch: Option<Epoch>,
    pub default_sigma_px: f64,
    pub wcs_sigma_arcsec: f64,
    pub search_radius_arcsec: f64,
    pub max_stars: usize,
}

impl Default for AtOptions {
    fn default() -> Self {
        Self {
            epoch: None,
            default_sigma_px: 1.0,
            wcs_sigma_arcsec: 0.0,
            search_radius_arcsec: 30.0,
            max_stars: 3,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct AtStar {
    pub id: u64,
    pub pos: SkyPos,
    pub mag: Option<f32>,
    pub pixel: Option<PixelPos>,
    pub separation_arcsec: f64,
    pub sigma_arcsec: f64,
    pub p_consistent: f64,
    pub p_chance: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct AtResult {
    pub sky: SkyPos,
    pub pixel: Option<PixelPos>,
    pub position_sigma_arcsec: f64,
    pub stars: Vec<AtStar>,
    pub consistent: bool,
}

const SHARED_CONE_MAX: f64 = 2.0 * std::f64::consts::PI / 180.0;

pub const CONSISTENT_P: f64 = 0.0027;

pub fn stars_at(
    catalog: &dyn Catalog,
    wcs: &Wcs,
    positions: &[AtPosition],
    opt: &AtOptions,
) -> Result<Vec<AtResult>, SolveError> {
    let bad = |w: &str| Err(SolveError::InvalidInput(w.into()));
    let fin = |v: f64| v.is_finite();
    if !(fin(opt.default_sigma_px) && opt.default_sigma_px > 0.0) {
        return bad("default sigma must be > 0");
    }
    if !(fin(opt.wcs_sigma_arcsec) && opt.wcs_sigma_arcsec >= 0.0) {
        return bad("WCS sigma must be >= 0");
    }
    if !(fin(opt.search_radius_arcsec) && opt.search_radius_arcsec > 0.0) {
        return bad("search radius must be > 0");
    }
    if opt.epoch.is_some_and(|e| !e.is_plausible()) {
        return bad("epoch must be within 1800 to 2200 (Julian year)");
    }
    let mut placed = Vec::with_capacity(positions.len());
    for p in positions {
        placed.push(match *p {
            AtPosition::Pixel(px, s) => {
                if !(fin(px.x) && fin(px.y)) {
                    return bad("pixel positions must be finite");
                }
                let s = s.unwrap_or(opt.default_sigma_px);
                if !(fin(s) && s >= 0.0) {
                    return bad("position sigma must be >= 0");
                }
                (wcs.pixel_to_sky(px), s * local_scale_arcsec(wcs, px))
            }
            AtPosition::Sky(sky, s) => {
                let s = s.unwrap_or(0.0);
                if !(fin(s) && s >= 0.0) {
                    return bad("position sigma must be >= 0");
                }
                (sky, s)
            }
        });
    }
    if placed.is_empty() {
        return Ok(Vec::new());
    }
    let search = (opt.search_radius_arcsec / 3600.0).to_radians();
    if search >= std::f64::consts::FRAC_PI_2 {
        return bad("search radius must be under 90 degrees");
    }
    let query = StarQuery {
        epoch: opt.epoch,
        ..StarQuery::default()
    };
    let cone = |centre: SkyPos, r: f64| -> Result<(Vec<Star>, f64), SolveError> {
        let stars = catalog.cone(centre, r, &query)?;
        let area_arcsec2 =
            4.0 * std::f64::consts::PI * (0.5 * r).psin().ppowi(2) * to_arcsec(1.0).ppowi(2);
        let density = stars.len() as f64 / area_arcsec2.max(f64::MIN_POSITIVE);
        Ok((stars, density))
    };
    let mut c = [0.0f64; 3];
    for (s, _) in &placed {
        let v = s.to_unit_vector().0;
        c = [c[0] + v[0], c[1] + v[1], c[2] + v[2]];
    }
    let centre = UnitVector::normalize(c).map(UnitVector::to_sky);
    let spread = centre.map(|c| {
        placed
            .iter()
            .map(|(s, _)| c.separation(*s))
            .fold(0.0f64, f64::max)
    });
    let cones: Vec<(Vec<Star>, f64)> = match (centre, spread) {
        (Some(c), Some(sp)) if sp + search <= SHARED_CONE_MAX => vec![cone(c, sp + search)?],
        _ => placed
            .iter()
            .map(|(s, _)| cone(*s, search))
            .collect::<Result<_, _>>()?,
    };
    let shared = cones.len() == 1;
    Ok(placed
        .iter()
        .enumerate()
        .map(|(i, &(sky, pos_sigma))| {
            let (stars, density) = &cones[if shared { 0 } else { i }];
            let density = *density;
            let mut near: Vec<AtStar> = stars
                .iter()
                .filter_map(|st| {
                    let d = to_arcsec(sky.separation(st.pos));
                    (d <= opt.search_radius_arcsec).then(|| {
                        let [e1, e2] = st.pos_err_mas.map(f64::from);
                        let cat = (0.5 * (e1 * e1 + e2 * e2)).sqrt() / 1000.0;
                        let wcs_axis = opt.wcs_sigma_arcsec / std::f64::consts::SQRT_2;
                        let sigma =
                            (pos_sigma * pos_sigma + cat * cat + wcs_axis * wcs_axis).sqrt();
                        let z2 = if sigma > 0.0 {
                            (d / sigma).ppowi(2)
                        } else if d > 0.0 {
                            f64::INFINITY
                        } else {
                            0.0
                        };
                        AtStar {
                            id: st.id,
                            pos: st.pos,
                            mag: st.mag,
                            pixel: wcs.sky_to_pixel(st.pos),
                            separation_arcsec: d,
                            sigma_arcsec: sigma,
                            p_consistent: (-0.5 * z2).pexp(),
                            p_chance: 1.0 - (-density * std::f64::consts::PI * d * d).pexp(),
                        }
                    })
                })
                .collect();
            near.sort_by(|a, b| {
                a.separation_arcsec
                    .total_cmp(&b.separation_arcsec)
                    .then(a.id.cmp(&b.id))
            });
            let consistent = near.first().is_some_and(|s| s.p_consistent >= CONSISTENT_P);
            near.truncate(opt.max_stars);
            AtResult {
                sky,
                pixel: wcs.sky_to_pixel(sky),
                position_sigma_arcsec: pos_sigma,
                consistent,
                stars: near,
            }
        })
        .collect())
}

fn local_scale_arcsec(wcs: &Wcs, p: PixelPos) -> f64 {
    let c = wcs.pixel_to_sky(p);
    let dx = c.separation(wcs.pixel_to_sky(PixelPos::new(p.x + 1.0, p.y)));
    let dy = c.separation(wcs.pixel_to_sky(PixelPos::new(p.x, p.y + 1.0)));
    to_arcsec((dx * dy).sqrt())
}
