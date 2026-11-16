//! "Is there a catalogue star here?" (`fayde verify --at`).
//!
//! Given a WCS and positions (pixel or sky), report for each position the nearest catalogue
//! stars at the observation epoch: separation, magnitude, the combined positional error, and
//! two probabilities a caller can adjudicate with:
//!
//! - `p_consistent`: how often a source *at* the star would be found this far away or farther,
//!   given the combined 1σ error (isotropic 2-D Gaussian: `exp(−d²/2σ²)`, the Rayleigh tail);
//! - `p_chance`: how often an arbitrary position has some catalogue star this close, from the
//!   local star density (Poisson: `1 − exp(−ρ π d²)`).
//!
//! A position is `consistent` with its nearest star when `p_consistent ≥ 0.0027` (the 3σ
//! two-sided tail of a 1-D Gaussian, used here as the conventional cut). Combined error:
//! `σ² = σ_position² + σ_catalogue² + σ_wcs²`, with the pixel σ converted at the WCS's local
//! scale and the catalogue error propagated to the epoch.
//!
//! Sources: Rayleigh distribution of the radial offset of a 2-D isotropic Gaussian; Poisson
//! nearest-neighbour distance for a uniform density (e.g. Chandrasekhar 1943, Rev. Mod. Phys.
//! 15, 1, §III.7).

use fayde_catalog::{Catalog, Star, StarQuery};
use fayde_core::{to_arcsec, Epoch, PixelPos, Portable, SkyPos, UnitVector};
use fayde_wcs::Wcs;

use crate::SolveError;

/// One position to check.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AtPosition {
    /// Pixel position (internal 0-based convention) and its 1σ, px (`None`: the default).
    Pixel(PixelPos, Option<f64>),
    /// Sky position and its 1σ, arcsec (`None`: zero).
    Sky(SkyPos, Option<f64>),
}

/// Options of [`stars_at`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AtOptions {
    /// Catalogue epoch (the observation's); `None` is the catalogue epoch.
    pub epoch: Option<Epoch>,
    /// 1σ of a pixel position that has none of its own, px.
    pub default_sigma_px: f64,
    /// Radial 1σ of the WCS itself at the positions, arcsec: a solve's `worst_sigma_arcsec`
    /// (√2 × the per-axis σ). Converted to per-axis before it is combined.
    pub wcs_sigma_arcsec: f64,
    /// Stars farther than this are not reported, arcsec.
    pub search_radius_arcsec: f64,
    /// Stars reported per position, nearest first.
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

/// A catalogue star near a position.
#[derive(Debug, Clone, PartialEq)]
pub struct AtStar {
    /// Catalogue id (UCAC4: running number).
    pub id: u64,
    /// Position at the epoch.
    pub pos: SkyPos,
    /// Magnitude.
    pub mag: Option<f32>,
    /// Where the WCS puts the star, px (internal convention), if it projects.
    pub pixel: Option<PixelPos>,
    /// Separation from the position, arcsec.
    pub separation_arcsec: f64,
    /// Combined per-axis 1σ (position, catalogue, WCS), arcsec.
    pub sigma_arcsec: f64,
    /// Rayleigh tail: probability of this offset or more if the position is this star.
    pub p_consistent: f64,
    /// Probability that some catalogue star lies this close to an arbitrary position.
    pub p_chance: f64,
}

/// The answer for one position.
#[derive(Debug, Clone, PartialEq)]
pub struct AtResult {
    /// The position on the sky.
    pub sky: SkyPos,
    /// The position in pixels (internal convention), if the sky position projects.
    pub pixel: Option<PixelPos>,
    /// 1σ of the position itself, arcsec.
    pub position_sigma_arcsec: f64,
    /// Nearest stars within the search radius, nearest first.
    pub stars: Vec<AtStar>,
    /// The nearest star's `p_consistent ≥ 0.0027`.
    pub consistent: bool,
}

/// Positions closer together than this (plus the search radius) share one catalogue cone, rad.
const SHARED_CONE_MAX: f64 = 2.0 * std::f64::consts::PI / 180.0;

/// `p_consistent` cut: the two-sided 3σ tail of a 1-D Gaussian.
pub const CONSISTENT_P: f64 = 0.0027;

/// Check positions against the catalogue under `wcs`. Positions within about 2° share one
/// catalogue cone; others get one each.
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
    // Each position on the sky with its own 1σ (arcsec).
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
    // Stars and local density (per arcsec²) of a cone. The cap area is 4π sin²(r/2) sr, which
    // stays accurate for tiny radii (2π(1 − cos r) cancels to zero).
    let cone = |centre: SkyPos, r: f64| -> Result<(Vec<Star>, f64), SolveError> {
        let stars = catalog.cone(centre, r, &query)?;
        let area_arcsec2 =
            4.0 * std::f64::consts::PI * (0.5 * r).psin().ppowi(2) * to_arcsec(1.0).ppowi(2);
        let density = stars.len() as f64 / area_arcsec2.max(f64::MIN_POSITIVE);
        Ok((stars, density))
    };
    // Positions within a field share one cone around their mean direction; positions spread
    // wider (e.g. across the sky) get a cone each, so no query covers more sky than needed.
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
                        // Per-axis σ throughout; the WCS σ is given radial (as a solve's
                        // `worst_sigma_arcsec`), √2 × per-axis.
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
            near.truncate(opt.max_stars);
            AtResult {
                sky,
                pixel: wcs.sky_to_pixel(sky),
                position_sigma_arcsec: pos_sigma,
                consistent: near.first().is_some_and(|s| s.p_consistent >= CONSISTENT_P),
                stars: near,
            }
        })
        .collect())
}

/// Local pixel scale of `wcs` at `p`, arcsec/px: the geometric mean of the two axis steps.
fn local_scale_arcsec(wcs: &Wcs, p: PixelPos) -> f64 {
    let c = wcs.pixel_to_sky(p);
    let dx = c.separation(wcs.pixel_to_sky(PixelPos::new(p.x + 1.0, p.y)));
    let dy = c.separation(wcs.pixel_to_sky(PixelPos::new(p.x, p.y + 1.0)));
    to_arcsec((dx * dy).sqrt())
}
