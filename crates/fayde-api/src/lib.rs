#![forbid(unsafe_code)]

use std::time::Duration;

use fayde_catalog::{ucac4::Ucac4Catalog, Catalog};
use fayde_core::{Epoch, PixelConvention, PixelPos, Portable, SkyPos};
use fayde_detect::FrameMode;
use fayde_fits::{header_priors, parse_parity, InstrumentProfile};
use fayde_index::Index;
use fayde_solve::{
    refine_wcs, solve_auto, stars_at, verify_wcs, AtOptions, Frame, Match, Outcome, PointingPrior,
    Policy, Priors, RollPrior, ScalePrior, Solution, SolveRequest, Source, StageTrace, Trail,
    WcsCovariance,
};
use fayde_wcs::{Header, Wcs};
use serde_json::{json, Value};

mod calibrate;
mod image;

pub use calibrate::{calibrate_profile, CalibrateOptions, SolvedFrame, CALIBRATE_MIN_FRAMES};
pub use fayde_solve::AtPosition;
pub use image::{ImageRequestError, PreparedImage};
pub use serde_json;

/// CLI and Python request conventions. Angles are in degrees, pixel scale in arcseconds per pixel,
/// and epoch in Julian years.
#[derive(Debug, Clone, Default)]
pub struct RequestOptions {
    pub width: Option<u32>,
    pub height: Option<u32>,
    pub ra: Option<f64>,
    pub dec: Option<f64>,
    pub radius: Option<f64>,
    pub pointing_sigma: Option<f64>,
    pub scale: Option<f64>,
    /// Fractional one-sigma uncertainty in pixel scale (0.02 means 2%).
    pub scale_sigma: Option<f64>,
    pub roll: Option<f64>,
    pub roll_sigma: Option<f64>,
    pub roll_max: Option<f64>,
    pub parity: Option<String>,
    pub trail_length: Option<f64>,
    pub trail_angle: Option<f64>,
    pub epoch: Option<f64>,
    pub deadline: Option<Duration>,
    pub max_sip_order: Option<usize>,
    pub max_hypotheses: Option<u64>,
    pub pixel_error: Option<f64>,
    pub accept_log_odds: Option<f64>,
    pub no_priors: bool,
    pub no_pointing: bool,
    pub no_header_wcs: bool,
    pub search_outside_pointing: bool,
    pub relax_priors_on_failure: bool,
    pub header: Option<Header>,
    pub instrument: Option<InstrumentProfile>,
}

#[derive(Debug, Clone, Copy, Default)]
pub struct AtCheck {
    pub wcs_sigma_arcsec: Option<f64>,
    pub search_radius_arcsec: Option<f64>,
    pub max_stars: Option<usize>,
}

#[derive(Debug, Clone)]
pub struct Built {
    pub request: SolveRequest,
    pub priors: Value,
    pub warnings: Vec<String>,
}

pub fn build_request(
    sources: Vec<Source>,
    list_size: (Option<u32>, Option<u32>),
    o: &RequestOptions,
) -> Result<Built, String> {
    let header = o.header.as_ref();
    let no_header_wcs =
        o.no_header_wcs || o.instrument.as_ref().and_then(|p| p.header_wcs) == Some(false);
    let size = (o.width.or(list_size.0), o.height.or(list_size.1));
    let hp = header.map(|h| {
        if no_header_wcs {
            fayde_fits::header_priors_for_size(&h.without_wcs(), size)
        } else {
            fayde_fits::header_priors_for_size(h, size)
        }
    });
    let mut warnings: Vec<String> = hp.iter().flat_map(|h| h.warnings.iter().cloned()).collect();
    let ip = o
        .instrument
        .as_ref()
        .map(|p| p.priors(&header.cloned().unwrap_or_default()));
    let width = o
        .width
        .or(list_size.0)
        .or(hp.as_ref().and_then(|h| h.size.map(|s| s.0)));
    let height = o
        .height
        .or(list_size.1)
        .or(hp.as_ref().and_then(|h| h.size.map(|s| s.1)));
    let (Some(width), Some(height)) = (width, height) else {
        return Err("image size unknown: give width/height, IMAGEW/IMAGEH, or a header".into());
    };

    let no_priors = o.no_priors;
    let given = [
        ("scale", o.scale.is_some()),
        ("ra/dec", o.ra.is_some() || o.dec.is_some()),
        ("roll", o.roll.is_some()),
        ("parity", o.parity.is_some()),
    ];
    for (name, _) in given.iter().filter(|(_, g)| no_priors && *g) {
        warnings.push(format!("{name} was given but dropped by no-priors"));
    }
    if o.no_pointing && !no_priors && (o.ra.is_some() || o.dec.is_some()) {
        warnings.push("ra/dec was given but dropped by no-pointing".into());
    }
    let no_pointing = no_priors || o.no_pointing;
    let pointing = match (o.ra, o.dec) {
        (Some(ra), Some(dec)) => Some((
            SkyPos::from_degrees(ra, dec).ok_or("bad ra/dec")?,
            "command line",
        )),
        (None, None) => hp.as_ref().and_then(|h| h.pointing),
        _ => return Err("ra and dec go together".into()),
    }
    .filter(|_| !no_pointing);
    let profile = o.instrument.as_ref();
    let radius = o
        .radius
        .or(profile.and_then(|p| p.pointing_radius_deg))
        .unwrap_or(2.0);
    let pointing_sigma = o
        .pointing_sigma
        .or(profile.and_then(|p| p.pointing_sigma_deg))
        .unwrap_or(radius / 4.0);
    let cli_sigma = o.scale_sigma;
    let scale: Option<(f64, f64, &str)> = match o.scale {
        Some(s) => Some((s, s * cli_sigma.unwrap_or(0.02), "command line")),
        None => match ip.and_then(|p| p.scale) {
            Some((s, sig)) => Some((s, cli_sigma.map_or(sig, |f| f * s), "instrument profile")),
            None => hp
                .as_ref()
                .and_then(|h| h.scale_arcsec)
                .map(|(s, src)| (s, s * cli_sigma.unwrap_or(0.02), src)),
        },
    }
    .filter(|_| !no_priors);
    let roll_sigma = o.roll_sigma.map(f64::to_radians);
    let roll_max = o.roll_max.map(f64::to_radians);
    let roll: Option<(f64, f64, f64, &str)> = match o.roll {
        Some(r) => Some((
            r.to_radians(),
            0.5f64.to_radians(),
            3f64.to_radians(),
            "command line",
        )),
        None => ip
            .and_then(|p| p.roll)
            .map(|(a, s, m)| (a, s, m, "instrument profile")),
    }
    .map(|(a, s, m, src)| (a, roll_sigma.unwrap_or(s), roll_max.unwrap_or(m), src))
    .filter(|_| !no_priors);
    let parity = match &o.parity {
        Some(p) => Some((
            parse_parity(p).ok_or("parity must be negative or positive")?,
            "command line",
        )),
        None => ip.and_then(|p| p.parity).map(|p| (p, "instrument profile")),
    }
    .filter(|_| !no_priors);
    let epoch = match o.epoch {
        Some(y) => Some((
            Epoch::julian_year(y).ok_or("epoch must be a finite Julian year")?,
            "command line",
        )),
        None => hp.as_ref().and_then(|h| h.epoch),
    };
    if epoch.is_none() {
        warnings.push(
            "no observation time (JD-OBS, MJD-OBS, DATE-OBS): catalogue positions at J2000.0; \
             give --epoch for stars with large proper motions"
                .into(),
        );
    }
    let trail = match (o.trail_length, o.trail_angle) {
        (Some(l), Some(a)) => Some(Trail {
            length_px: l,
            angle: a.to_radians(),
        }),
        (None, None) => None,
        _ => return Err("trail length and trail angle go together".into()),
    };
    let defaults = Policy::default();
    let request = SolveRequest {
        frame: Frame {
            width,
            height,
            trail,
        },
        sources,
        priors: Priors {
            pointing: pointing.map(|(c, _)| PointingPrior {
                centre: c,
                sigma: pointing_sigma.to_radians(),
                max_radius: radius.to_radians(),
            }),
            scale: scale.map(|(s, sig, _)| ScalePrior {
                arcsec_per_px: s,
                sigma: sig,
            }),
            roll: roll.map(|(angle, sigma, max_dev, _)| RollPrior {
                angle,
                sigma,
                max_dev,
            }),
            parity: parity.map(|(p, _)| p),
        },
        policy: Policy {
            epoch: epoch.map(|(e, _)| e),
            max_sip_order: o.max_sip_order.unwrap_or(3),
            search_outside_pointing: o.search_outside_pointing,
            relax_priors_on_failure: o.relax_priors_on_failure,
            deadline: o.deadline,
            default_sigma_px: o.pixel_error.unwrap_or(defaults.default_sigma_px),
            accept_log_odds: o.accept_log_odds.unwrap_or(defaults.accept_log_odds),
            max_hypotheses: o.max_hypotheses.unwrap_or(defaults.max_hypotheses),
            ..defaults
        },
    };
    let priors = json!({
        "pointing": pointing.map(|(p, src)| json!({"ra_deg": p.to_degrees().0, "dec_deg": p.to_degrees().1, "radius_deg": radius, "sigma_deg": pointing_sigma, "origin": src})),
        "scale": scale.map(|(s, sig, src)| json!({"arcsec_per_px": s, "sigma_arcsec_per_px": sig, "origin": src})),
        "roll": roll.map(|(a, s, m, src)| json!({"deg": a.to_degrees(), "sigma_deg": s.to_degrees(), "max_dev_deg": m.to_degrees(), "origin": src})),
        "parity": parity.map(|(p, src)| json!({"parity": format!("{p:?}"), "origin": src})),
        "epoch": epoch.map(|(e, src)| json!({"julian_year": e.year(), "origin": src})),
        "instrument": o.instrument.as_ref().and_then(|p| p.name.clone()),
        "trail": trail.map(|t| json!({"length_px": t.length_px, "angle_deg": t.angle.to_degrees()})),
        "truncated_sources": request.sources.iter().filter(|s| s.truncated).count(),
        "pixel_error_px": request.policy.default_sigma_px,
    });
    Ok(Built {
        request,
        priors,
        warnings,
    })
}

pub struct Solver {
    catalog: Ucac4Catalog,
    indexes: Vec<Index>,
    index_paths: Vec<String>,
}

pub struct Report {
    pub json: Value,
    pub solution: Option<Box<Solution>>,
    pub error: bool,
}

pub enum Operation {
    Solve,
    Verify(Wcs),
    Refit(Wcs),
}

impl Solver {
    pub fn open(ucac4_dir: &str, index_paths: &[String]) -> Result<Self, String> {
        let absolute = |p: &str| {
            std::path::absolute(p)
                .map(|a| a.to_string_lossy().into_owned())
                .map_err(|e| format!("{p}: {e}"))
        };
        let catalog = Ucac4Catalog::open(absolute(ucac4_dir)?).map_err(|e| e.to_string())?;
        let index_paths: Vec<String> = index_paths
            .iter()
            .map(|p| absolute(p))
            .collect::<Result<_, _>>()?;
        let indexes = index_paths
            .iter()
            .map(|p| Index::read(p).map_err(|e| format!("{p}: {e}")))
            .collect::<Result<_, _>>()?;
        Ok(Self {
            catalog,
            indexes,
            index_paths,
        })
    }

    pub fn catalog(&self) -> &Ucac4Catalog {
        &self.catalog
    }

    pub fn run(&self, b: &Built, operation: &Operation) -> Report {
        match operation {
            Operation::Solve => self.solve(b),
            Operation::Verify(wcs) => self.verify(b, wcs),
            Operation::Refit(wcs) => self.refit(b, wcs),
        }
    }

    pub fn solve(&self, b: &Built) -> Report {
        let out = solve_auto(&b.request, &self.catalog, &self.indexes);
        self.outcome(out, b, true)
    }

    pub fn refit(&self, b: &Built, wcs: &Wcs) -> Report {
        let out = refine_wcs(&b.request, &self.catalog, wcs);
        self.outcome(out, b, false)
    }

    pub fn verify(&self, b: &Built, wcs: &Wcs) -> Report {
        let req = &b.request;
        let mut json = match verify_wcs(req, &self.catalog, wcs) {
            Ok(v) => json!({
                "accepted": v.accepted,
                "log_odds": v.log_odds,
                "threshold": req.policy.accept_log_odds,
                "matches": v.matches.len(),
                "confident_matches": v.matches.iter().filter(|m| m.inlier_prob > 0.5).count(),
                "rms_px": v.rms_px,
                "references": v.references,
                "sources": req.sources.len(),
                "pairs": pairs_json(&v.matches, req),
                "priors": b.priors,
            }),
            Err(e) => return self.failed(&e.to_string(), b, false),
        };
        json["identity"] = self.identity(req.policy.epoch, false);
        Report {
            json: with_warnings(json, b),
            solution: None,
            error: false,
        }
    }

    pub fn check_at(
        &self,
        wcs: &Wcs,
        positions: &[AtPosition],
        o: &RequestOptions,
        at: &AtCheck,
    ) -> Report {
        let epoch = match o.epoch {
            Some(y) => match Epoch::julian_year(y) {
                Some(e) => Some(e),
                None => return self.at_failed("epoch must be a finite Julian year"),
            },
            None => o
                .header
                .as_ref()
                .and_then(|h| header_priors(h).epoch)
                .map(|(e, _)| e),
        };
        let defaults = AtOptions::default();
        let opt = AtOptions {
            epoch,
            default_sigma_px: o.pixel_error.unwrap_or(defaults.default_sigma_px),
            wcs_sigma_arcsec: at.wcs_sigma_arcsec.unwrap_or(defaults.wcs_sigma_arcsec),
            search_radius_arcsec: at
                .search_radius_arcsec
                .unwrap_or(defaults.search_radius_arcsec),
            max_stars: at.max_stars.unwrap_or(defaults.max_stars),
        };
        let results = match stars_at(&self.catalog, wcs, positions, &opt) {
            Ok(r) => r,
            Err(e) => return self.at_failed(&e.to_string()),
        };
        let fits = |p: Option<PixelPos>| p.map(|p| PixelConvention::Fits.from_internal(p));
        let list: Vec<Value> = results
            .iter()
            .enumerate()
            .map(|(i, r)| {
                let (ra, dec) = r.sky.to_degrees();
                json!({
                    "index": i,
                    "ra_deg": ra,
                    "dec_deg": dec,
                    "pixel": fits(r.pixel).map(|(x, y)| [x, y]),
                    "position_sigma_arcsec": r.position_sigma_arcsec,
                    "consistent": r.consistent,
                    "stars": r.stars.iter().map(|s| {
                        let (ra, dec) = s.pos.to_degrees();
                        json!({
                            "star_id": s.id,
                            "ra_deg": ra,
                            "dec_deg": dec,
                            "mag": s.mag,
                            "pixel": fits(s.pixel).map(|(x, y)| [x, y]),
                            "separation_arcsec": s.separation_arcsec,
                            "sigma_arcsec": s.sigma_arcsec,
                            "p_consistent": s.p_consistent,
                            "p_chance": s.p_chance,
                        })
                    }).collect::<Vec<_>>(),
                })
            })
            .collect();
        let json = json!({
            "positions": list,
            "consistent": results.iter().filter(|r| r.consistent).count(),
            "consistent_p": fayde_solve::CONSISTENT_P,
            "options": {
                "epoch_julian_year": epoch.map(|e| e.year()),
                "default_sigma_px": opt.default_sigma_px,
                "wcs_sigma_arcsec": opt.wcs_sigma_arcsec,
                "search_radius_arcsec": opt.search_radius_arcsec,
                "max_stars": opt.max_stars,
            },
            "identity": self.identity(epoch, false),
        });
        Report {
            json,
            solution: None,
            error: false,
        }
    }

    fn at_failed(&self, reason: &str) -> Report {
        Report {
            json: json!({"error": true, "reason": reason, "identity": self.identity(None, false)}),
            solution: None,
            error: true,
        }
    }

    fn failed(&self, reason: &str, b: &Built, with_indexes: bool) -> Report {
        Report {
            json: with_warnings(
                json!({
                    "solved": false,
                    "error": true,
                    "reason": reason,
                    "priors": b.priors,
                    "identity": self.identity(b.request.policy.epoch, with_indexes),
                }),
                b,
            ),
            solution: None,
            error: true,
        }
    }

    fn outcome(
        &self,
        out: Result<Outcome, fayde_solve::SolveError>,
        b: &Built,
        with_indexes: bool,
    ) -> Report {
        let req = &b.request;
        let out = match out {
            Ok(o) => o,
            Err(e) => return self.failed(&e.to_string(), b, with_indexes),
        };
        let (mut json, solution) = match out {
            Outcome::Solved(sol) => {
                let (ra, dec) = sol.wcs.crval().to_degrees();
                let j = json!({
                    "solved": true,
                    "method": format!("{:?}", sol.method),
                    "crval_deg": [ra, dec],
                    "pixel_scale_arcsec": sol.wcs.pixel_scale_arcsec(),
                    "roll_deg": sol.wcs.cd()[0][1].patan2(sol.wcs.cd()[1][1]).to_degrees(),
                    "sip_order": sol.wcs.sip_order(),
                    "parity": format!("{:?}", sol.wcs.parity()),
                    "centre_sigma_arcsec": sol.covariance.radial_sigma_arcsec(sol.wcs.crpix()),
                    "worst_sigma_arcsec": sol.worst_sigma_arcsec,
                    "loo_rms_px": sol.loo_rms_px,
                    "log_odds": sol.log_odds,
                    "nfa": sol.nfa,
                    "hypotheses_tested": sol.hypotheses_tested,
                    "matches": sol.matches.len(),
                    "rms_px": sol.rms_px,
                    "error_scale": sol.error_scale,
                    "index_hash": sol.index_hash.map(|h| format!("{h:016x}")),
                    "seconds": sol.elapsed.as_secs_f64(),
                    "prior_violated": sol.prior_violated,
                    "violated_priors": sol.violated_priors,
                    "solved_by_stage": sol.search.last().map(|s| s.stage.as_str()),
                    "relaxed_priors": sol.search.last().is_some_and(|s| s.stage.starts_with("relaxed priors: ")),
                    "search": trace_json(&sol.search),
                    "priors": b.priors,
                    "wcs_cards": solution_cards(&sol, req),
                    "covariance": covariance_json(&sol.covariance),
                    "pairs": pairs_json(&sol.matches, req),
                });
                (j, Some(sol))
            }
            Outcome::Abstained {
                reason,
                best_log_odds,
                hypotheses_tested,
                search,
            } => (
                json!({
                    "solved": false,
                    "reason": reason,
                    "best_log_odds": best_log_odds,
                    "hypotheses_tested": hypotheses_tested,
                    "search": trace_json(&search),
                    "priors": b.priors,
                }),
                None,
            ),
        };
        json["identity"] = self.identity(req.policy.epoch, with_indexes);
        Report {
            json: with_warnings(json, b),
            solution,
            error: false,
        }
    }

    fn identity(&self, epoch: Option<Epoch>, with_indexes: bool) -> Value {
        let id = self.catalog.id();
        let n = if with_indexes { self.indexes.len() } else { 0 };
        json!({
            "fayde": {"version": env!("CARGO_PKG_VERSION"), "commit": commit()},
            "catalog": {
                "name": id.name,
                "release": id.release,
                "location": id.location,
                "reference_epoch": "J2000.0",
                "solve_epoch_julian_year": epoch.map(Epoch::year),
            },
            "indexes": self.index_paths.iter().zip(&self.indexes).take(n).map(|(path, i)| {
                let m = i.meta();
                json!({
                    "path": path,
                    "content_hash": format!("{:016x}", i.content_hash()),
                    "fov_deg": m.params.fov_deg,
                    "stars": m.star_count,
                    "patterns": m.pattern_count,
                    "epoch_year": m.params.epoch_year,
                })
            }).collect::<Vec<_>>(),
        })
    }
}

fn with_warnings(mut json: Value, b: &Built) -> Value {
    if !b.warnings.is_empty() {
        json["warnings"] = b.warnings.clone().into();
    }
    json
}

#[derive(Debug, Clone, Copy, Default)]
pub struct DetectOptions {
    pub fwhm: Option<f64>,
    pub threshold: Option<f64>,
    pub saturation: Option<f64>,
    pub gain: Option<f64>,
    pub max_sources: Option<usize>,
    pub trail_angle: Option<f64>,
    pub points_only: bool,
}

pub const DEFAULT_MAX_SOURCES: usize = 2000;

pub const FWHM_PER_SIGMA: f64 = 2.354_820_045_030_949;

/// Return zero-based solver sources and a detection report with one-based FITS positions.
pub fn detect_sources(
    image: &fayde_fits::Image,
    o: &DetectOptions,
) -> Result<(Vec<Source>, Value), String> {
    let result = detect_image(image, o, DetectionReport::Full)?;
    Ok((result.sources, result.report))
}

struct ImageDetection {
    sources: Vec<Source>,
    mode: FrameMode,
    trail: Option<(f64, f64)>,
    report: Value,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum DetectionReport {
    None,
    Summary,
    Full,
}

fn detect_image(
    image: &fayde_fits::Image,
    o: &DetectOptions,
    detail: DetectionReport,
) -> Result<ImageDetection, String> {
    use fayde_detect::Class;
    let hdr = |keys: &[&str]| {
        keys.iter().find_map(|k| {
            image
                .header
                .get(k)
                .and_then(|v| v.as_f64().or_else(|| v.as_str()?.trim().parse().ok()))
                .filter(|v: &f64| v.is_finite() && *v > 0.0)
        })
    };
    let positive = |v: Option<f64>| v.is_none_or(|x| x.is_finite() && x > 0.0);
    if o.fwhm.is_some_and(|f| !(0.5..=50.0).contains(&f)) {
        return Err("fwhm must be in [0.5, 50] px".into());
    }
    if o.threshold
        .is_some_and(|t| !(t.is_finite() && t >= fayde_detect::MIN_THRESHOLD))
    {
        return Err(format!(
            "the detection threshold must be at least {}",
            fayde_detect::MIN_THRESHOLD
        ));
    }
    if o.trail_angle.is_some_and(|a| !a.is_finite()) {
        return Err("the trail angle must be a number of degrees".into());
    }
    if !positive(o.saturation) || !positive(o.gain) {
        return Err("saturation and gain must be > 0".into());
    }
    let saturation = o.saturation.or_else(|| hdr(&["SATURATE"]));
    let gain = o.gain.or_else(|| hdr(&["GAIN", "EGAIN"]));
    let p = fayde_detect::DetectOptions {
        psf_sigma: o.fwhm.map(|f| f / FWHM_PER_SIGMA),
        threshold: o.threshold,
        saturation: saturation.map(|s| s as f32),
        gain,
        max_sources: o.max_sources.unwrap_or(DEFAULT_MAX_SOURCES),
        trail_angle: o.trail_angle.map(f64::to_radians),
        trails: !o.points_only,
        ..fayde_detect::DetectOptions::default()
    };
    let started = std::time::Instant::now();
    let r = fayde_detect::detect(&image.data, image.width, image.height, &p)
        .map_err(|e| e.to_string())?;
    let seconds = started.elapsed().as_secs_f64();
    let sources: Vec<Source> = r
        .stars()
        .map(|d| {
            let mut s = Source::at(d.x, d.y);
            s.flux = Some(d.flux);
            s.sigma_px = Some(d.sigma_px);
            s.p_star = Some(d.p_star);
            s.truncated = d.truncated;
            s
        })
        .collect();
    let mode = r.census.mode;
    let trail = r.trail.and_then(|t| {
        let angle = t.angle.to_degrees();
        (t.length.is_finite() && angle.is_finite()).then_some((t.length, angle))
    });
    if detail == DetectionReport::None {
        return Ok(ImageDetection {
            sources,
            mode,
            trail,
            report: Value::Null,
        });
    }
    let fits = |x: f64, y: f64| PixelConvention::Fits.from_internal(PixelPos::new(x, y));
    let finite = |v: f64| v.is_finite().then_some(v);
    let row = |d: &fayde_detect::Detected| -> Value {
        let (x, y) = fits(d.x, d.y);
        let major = d.sigma_major.min(10.0 * d.sigma_minor.max(0.1));
        let mean2 = 0.5 * (d.sigma_minor * d.sigma_minor + major * major);
        json!({
            "x": x, "y": y, "sigma_px": d.sigma_px,
            "sigma_along_px": d.sigma_along, "sigma_across_px": d.sigma_across,
            "flux": d.flux, "snr": d.snr,
            "p_star": d.p_star, "class": d.class.name(),
            "fwhm": FWHM_PER_SIGMA * mean2.sqrt(),
            "sigma_major_px": finite(d.sigma_major), "sigma_minor_px": d.sigma_minor,
            "theta_deg": d.theta.to_degrees(), "streak": d.streak,
            "truncated": d.truncated, "saturated": d.class == Class::Saturated,
        })
    };
    let count = |c: Class| r.sources.iter().filter(|s| s.class == c).count();
    let s = &r.stats;
    let mut warnings: Vec<String> = Vec::new();
    if !s.psf_measured && o.fwhm.is_none() {
        warnings.push("too few stars to measure their width: the filter's is assumed".into());
    }
    if let (Some(level), None) = (saturation, s.saturation) {
        warnings.push(format!(
            "saturation level {level} is at the background: not used (check SATURATE or the \
             saturation option)"
        ));
    }
    let mut report = json!({
        "image": {"width": image.width, "height": image.height, "hdu": image.hdu},
        "fwhm_px": s.psf_sigma * FWHM_PER_SIGMA,
        "fwhm_estimated": s.psf_measured,
        "background": s.background,
        "noise": s.noise,
        "threshold_sigma": s.threshold,
        "saturation": s.saturation,
        "gain": gain,
        "count": sources.len(),
        "classes": {
            "star": count(Class::Star),
            "saturated": count(Class::Saturated),
            "extended": count(Class::Extended),
            "sharp": count(Class::Sharp),
            "streak": count(Class::Streak),
            "point": count(Class::Point),
        },
        "mode": mode.name(),
        "census": {
            "reason": r.census.reason,
            "votes": r.census.votes,
            "mode_votes": r.census.mode_votes,
            "mode_share": r.census.mode_share,
            "second_share": r.census.second_share,
            "elongated_share": r.census.elongated_share,
            "shapes": r.census.shapes.map(|(a, l)| json!({"angle_deg": a.to_degrees(), "length_px": l})),
            "autocovariance": r.census.autocov.map(|(a, l)| json!({"angle_deg": a.to_degrees(), "length_px": l})),
            "regions": r.census.regions.0,
            "region_angle_deviation_deg": r.census.regions.1.to_degrees(),
            "region_length_deviation": r.census.regions.2,
        },
        "trail": r.trail.map(|t| json!({
            "length_px": t.length,
            "angle_deg": t.angle.to_degrees(),
            "width_sigma_px": t.width,
            "sigma_length_px": t.sigma_length,
            "sigma_angle_deg": t.sigma_angle.to_degrees(),
            "fitted": t.fitted,
        })),
        "kernel": r.kernel.map(|k| json!({
            "width_px": k.width,
            "height_px": k.height,
            "stacked": k.stacked,
            "path_length_px": k.path_length,
            "bend_px": k.bend,
            "regions": k.regions,
        })),
        "rejected_narrow": count(Class::Sharp),
        "bin": s.bin,
        "kernel_sigma_px": s.kernel_sigma,
        "passes": s.passes,
        "undefined_pixels": s.undefined,
        "seconds": seconds,
        "streaks": r.streaks.iter().map(|t| {
            let (x0, y0) = fits(t.x0, t.y0);
            let (x1, y1) = fits(t.x1, t.y1);
            json!({
                "x0": x0, "y0": y0, "x1": x1, "y1": y1,
                "length_px": t.length,
                "angle_deg": t.angle.to_degrees(),
                "width_sigma_px": t.width,
                "flux": t.flux,
                "snr": t.snr,
                "coverage": t.coverage,
                "touches_edge": t.touches_edge,
                "cosmic_ray": t.cosmic_ray,
                "photometry_masked": t.photometry_masked,
            })
        }).collect::<Vec<_>>(),
        "warnings": warnings,
    });
    if detail == DetectionReport::Full {
        report["sources"] = Value::Array(r.stars().map(row).collect());
        report["rejected"] = Value::Array(
            r.sources
                .iter()
                .filter(|d| !d.class.is_stellar())
                .map(row)
                .collect(),
        );
    }
    Ok(ImageDetection {
        sources,
        mode,
        trail,
        report,
    })
}

pub fn retry_as_points(detection: &Value) -> bool {
    FrameMode::POINT_RETRY_MODES
        .iter()
        .any(|m| Some(m.name()) == detection["mode"].as_str())
}

pub fn use_trail(opts: &mut RequestOptions, trail: Option<(f64, f64)>) {
    if opts.trail_length.is_none() && opts.trail_angle.is_none() {
        if let Some((l, a)) = trail {
            (opts.trail_length, opts.trail_angle) = (Some(l), Some(a));
        }
    }
}

pub fn detected_trail(report: &Value) -> Option<(f64, f64)> {
    let t = report.get("trail")?;
    Some((t.get("length_px")?.as_f64()?, t.get("angle_deg")?.as_f64()?))
}

pub fn header_from_json(m: &serde_json::Map<String, Value>) -> Result<Header, String> {
    use fayde_wcs::header::Value as H;
    let mut h = Header::default();
    for (k, v) in m {
        let key = k.trim().to_ascii_uppercase();
        let value = match v {
            Value::Bool(b) => H::Bool(*b),
            Value::Number(n) => match n.as_i64() {
                Some(i) => H::Int(i),
                None => H::Float(
                    n.as_f64()
                        .ok_or_else(|| format!("header {key}: bad number"))?,
                ),
            },
            Value::String(s) => H::Str(s.clone()),
            Value::Null => continue,
            _ => return Err(format!("header {key}: value must be a scalar")),
        };
        h.set(&key, value);
    }
    Ok(h)
}

pub fn commit() -> &'static str {
    env!("FAYDE_COMMIT")
}

pub fn solution_cards(sol: &Solution, req: &SolveRequest) -> Vec<String> {
    let mut cards = sol.wcs.to_cards();
    cards.push(format!(
        "{:<80}",
        format!("IMAGEW  = {:>20}", req.frame.width)
    ));
    cards.push(format!(
        "{:<80}",
        format!("IMAGEH  = {:>20}", req.frame.height)
    ));
    cards
}

fn covariance_json(c: &WcsCovariance) -> Value {
    let k = c.scale2 * c.rad_per_px * c.rad_per_px;
    let (cx, cy) = PixelConvention::Fits.from_internal(c.crpix);
    json!({
        "terms": c.terms,
        "matrix_rad2": c.inv_normal.iter().map(|v| v * k).collect::<Vec<_>>(),
        "norm_px": c.norm,
        "crpix_fits": [cx, cy],
    })
}

fn pairs_json(matches: &[Match], req: &SolveRequest) -> Value {
    matches
        .iter()
        .map(|m| {
            let s = &req.sources[m.source];
            let (x, y) = PixelConvention::Fits.from_internal(PixelPos::new(s.x, s.y));
            let (ra, dec) = m.star_pos.to_degrees();
            let (px, py) = PixelConvention::Fits.from_internal(PixelPos::new(
                s.x - m.residual_px[0],
                s.y - m.residual_px[1],
            ));
            json!({
                "source": m.source,
                "x": x,
                "y": y,
                "star_id": m.star_id,
                "ra_deg": ra,
                "dec_deg": dec,
                "mag": m.star_mag,
                "predicted_x": px,
                "predicted_y": py,
                "residual_px": m.residual_px,
                "inlier_prob": m.inlier_prob,
            })
        })
        .collect()
}

fn trace_json(trace: &[StageTrace]) -> Value {
    trace
        .iter()
        .map(|s| json!({"stage": s.stage, "hypotheses": s.hypotheses, "seconds": s.seconds, "outcome": s.outcome}))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn opts() -> RequestOptions {
        RequestOptions {
            width: Some(1024),
            height: Some(1024),
            ..RequestOptions::default()
        }
    }

    #[test]
    fn a_prior_dropped_by_a_switch_is_reported() {
        let b = build_request(
            vec![],
            (None, None),
            &RequestOptions {
                no_priors: true,
                scale: Some(3.0),
                ra: Some(10.0),
                dec: Some(2.0),
                ..opts()
            },
        )
        .unwrap();
        assert!(b.request.priors.scale.is_none() && b.request.priors.pointing.is_none());
        let dropped = |b: &Built| -> Vec<String> {
            b.warnings
                .iter()
                .filter(|w| !w.starts_with("no observation time"))
                .cloned()
                .collect()
        };
        assert_eq!(dropped(&b).len(), 2, "{:?}", b.warnings);
        assert!(dropped(&b)
            .iter()
            .all(|w| w.contains("dropped by no-priors")));
        let b = build_request(
            vec![],
            (None, None),
            &RequestOptions {
                no_pointing: true,
                scale: Some(3.0),
                ra: Some(10.0),
                dec: Some(2.0),
                ..opts()
            },
        )
        .unwrap();
        assert!(b.request.priors.scale.is_some() && b.request.priors.pointing.is_none());
        assert_eq!(dropped(&b), ["ra/dec was given but dropped by no-pointing"]);
        let quiet = RequestOptions {
            no_priors: true,
            ..opts()
        };
        assert!(dropped(&build_request(vec![], (None, None), &quiet).unwrap()).is_empty());
    }

    #[test]
    fn size_is_required() {
        let e = build_request(vec![], (None, None), &RequestOptions::default()).unwrap_err();
        assert!(e.contains("image size unknown"), "{e}");
        assert!(build_request(vec![], (Some(10), Some(20)), &RequestOptions::default()).is_ok());
    }

    #[test]
    fn pixel_error_and_threshold_reach_the_policy() {
        let b = build_request(
            vec![],
            (None, None),
            &RequestOptions {
                pixel_error: Some(3.0),
                accept_log_odds: Some(40.0),
                ..opts()
            },
        )
        .unwrap();
        assert_eq!(b.request.policy.default_sigma_px, 3.0);
        assert_eq!(b.request.policy.accept_log_odds, 40.0);
        assert_eq!(b.priors["pixel_error_px"], 3.0);
    }

    #[test]
    fn caller_scale_beats_profile_and_no_pointing_keeps_the_rest() {
        let profile = InstrumentProfile::from_json(
            r#"{"name": "cam", "scale_arcsec": 3.0, "scale_sigma_arcsec": 0.003, "parity": "positive"}"#,
        )
        .unwrap();
        let b = build_request(
            vec![],
            (None, None),
            &RequestOptions {
                ra: Some(10.0),
                dec: Some(20.0),
                scale: Some(2.0),
                no_pointing: true,
                instrument: Some(profile),
                ..opts()
            },
        )
        .unwrap();
        assert!(b.request.priors.pointing.is_none());
        assert_eq!(b.request.priors.scale.unwrap().arcsec_per_px, 2.0);
        assert!(b.request.priors.parity.is_some());
        assert_eq!(b.priors["instrument"], "cam");
    }

    #[test]
    fn profile_pointing_uncertainty_and_header_wcs_switch() {
        let header = Header::parse_cards([
            "NAXIS1  =                  512",
            "NAXIS2  =                  512",
            "CTYPE1  = 'RA---TAN'",
            "CTYPE2  = 'DEC--TAN'",
            "CRVAL1  =                 10.0",
            "CRVAL2  =                 20.0",
            "CD1_1   =               -0.001",
            "CD1_2   =                  0.0",
            "CD2_1   =                  0.0",
            "CD2_2   =                0.001",
            "RA      = '00 44 00.0'",
            "DEC     = '+21 30 00'",
        ]);
        let with = |o: RequestOptions| {
            build_request(
                vec![],
                (None, None),
                &RequestOptions {
                    header: Some(header.clone()),
                    ..o
                },
            )
            .unwrap()
        };
        let b = with(RequestOptions::default());
        let pt = b.request.priors.pointing.unwrap();
        assert!((pt.centre.to_degrees().0 - 10.0).abs() < 1e-9);
        assert!((pt.max_radius.to_degrees() - 2.0).abs() < 1e-12);
        assert!((pt.sigma.to_degrees() - 0.5).abs() < 1e-12);
        assert_eq!(b.priors["pointing"]["origin"], "CRVAL1/CRVAL2");
        assert!((b.request.priors.scale.unwrap().arcsec_per_px - 3.6).abs() < 1e-9);
        let b = with(RequestOptions {
            no_header_wcs: true,
            ..RequestOptions::default()
        });
        let pt = b.request.priors.pointing.unwrap();
        assert!((pt.centre.to_degrees().0 - 11.0).abs() < 1e-9);
        assert!((pt.centre.to_degrees().1 - 21.5).abs() < 1e-9);
        assert_eq!(b.priors["pointing"]["origin"], "RA/DEC");
        assert!(b.request.priors.scale.is_none());
        let profile = InstrumentProfile::from_json(
            r#"{"scale_arcsec": 3.0, "pointing_sigma_deg": 0.02, "pointing_radius_deg": 1.5,
                "header_wcs": false}"#,
        )
        .unwrap();
        let b = with(RequestOptions {
            instrument: Some(profile.clone()),
            ..RequestOptions::default()
        });
        let pt = b.request.priors.pointing.unwrap();
        assert_eq!(b.priors["pointing"]["origin"], "RA/DEC");
        assert!((pt.sigma.to_degrees() - 0.02).abs() < 1e-12);
        assert!((pt.max_radius.to_degrees() - 1.5).abs() < 1e-12);
        assert_eq!(b.priors["pointing"]["sigma_deg"], 0.02);
        let b = with(RequestOptions {
            instrument: Some(profile),
            radius: Some(3.0),
            pointing_sigma: Some(0.1),
            ..RequestOptions::default()
        });
        let pt = b.request.priors.pointing.unwrap();
        assert!((pt.sigma.to_degrees() - 0.1).abs() < 1e-12);
        assert!((pt.max_radius.to_degrees() - 3.0).abs() < 1e-12);
    }

    #[test]
    fn half_given_pairs_are_errors() {
        for o in [
            RequestOptions {
                ra: Some(1.0),
                ..opts()
            },
            RequestOptions {
                trail_length: Some(5.0),
                ..opts()
            },
            RequestOptions {
                parity: Some("sideways".into()),
                ..opts()
            },
        ] {
            assert!(build_request(vec![], (None, None), &o).is_err());
        }
    }
}
