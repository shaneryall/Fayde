//! The request-and-report layer shared by the `fayde` command line and the Python bindings.
//!
//! [`RequestOptions`] holds what a caller can say about one frame (the command-line request
//! options, typed). [`build_request`] combines it with the sources, the frame header and an
//! instrument profile into a [`SolveRequest`], applying the precedence command line (caller) >
//! instrument profile > header. [`Solver`] holds a catalogue and indexes opened once and turns
//! a request into the JSON report the CLI prints, so every front end reports the same fields.

#![forbid(unsafe_code)]

use std::time::Duration;

use fayde_catalog::{ucac4::Ucac4Catalog, Catalog};
use fayde_core::{Epoch, PixelConvention, PixelPos, Portable, SkyPos};
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

pub use calibrate::{calibrate_profile, CalibrateOptions, SolvedFrame, CALIBRATE_MIN_FRAMES};
pub use fayde_solve::AtPosition;
pub use serde_json;

/// What a caller says about one frame. `None` / `false` means "not given": the profile or the
/// header may then supply it. Angles in degrees, scale in arcsec/px.
#[derive(Debug, Clone, Default)]
pub struct RequestOptions {
    /// Image width, px. Else the source list's or header's.
    pub width: Option<u32>,
    /// Image height, px.
    pub height: Option<u32>,
    /// Pointing right ascension, degrees (with `dec`).
    pub ra: Option<f64>,
    /// Pointing declination, degrees.
    pub dec: Option<f64>,
    /// Hard pointing radius, degrees (default 2).
    pub radius: Option<f64>,
    /// 1σ of the pointing, degrees (default radius/4).
    pub pointing_sigma: Option<f64>,
    /// Pixel scale, arcsec/px.
    pub scale: Option<f64>,
    /// 1σ of the scale as a fraction (default 0.02, or the profile's).
    pub scale_sigma: Option<f64>,
    /// Roll, degrees (1σ 0.5°, limit 3° unless overridden).
    pub roll: Option<f64>,
    /// 1σ of the roll, degrees; overrides the profile's too.
    pub roll_sigma: Option<f64>,
    /// Hard roll limit, degrees; overrides the profile's too.
    pub roll_max: Option<f64>,
    /// `"negative"` or `"positive"`.
    pub parity: Option<String>,
    /// Trail length, px (with `trail_angle`).
    pub trail_length: Option<f64>,
    /// Trail angle, degrees.
    pub trail_angle: Option<f64>,
    /// Solve epoch, Julian year. Else the header's.
    pub epoch: Option<f64>,
    /// Wall-clock limit.
    pub deadline: Option<Duration>,
    /// Highest SIP order tried (default 3).
    pub max_sip_order: Option<usize>,
    /// Stop after this many hypotheses (default 5,000,000): bounds the work independently of
    /// the machine, unlike the deadline.
    pub max_hypotheses: Option<u64>,
    /// Default centroid 1σ, px, for sources without their own.
    pub pixel_error: Option<f64>,
    /// Log-odds needed to accept a match.
    pub accept_log_odds: Option<f64>,
    /// Drop every prior.
    pub no_priors: bool,
    /// Drop only the pointing.
    pub no_pointing: bool,
    /// Do not use a world coordinate system in the header for the pointing or the scale: read
    /// the pointing from the header's pointing keywords (`RA`/`DEC`, `OBJCTRA`/`OBJCTDEC`)
    /// only. For frames whose header WCS may be another solver's solution. The instrument
    /// profile can ask for the same (`header_wcs: false`).
    pub no_header_wcs: bool,
    /// Search the sky outside the pointing limit when the pointing finds nothing.
    pub search_outside_pointing: bool,
    /// Retry once with relaxed hard limits when the strict pass fails.
    pub relax_priors_on_failure: bool,
    /// The frame's primary header.
    pub header: Option<Header>,
    /// The instrument profile.
    pub instrument: Option<InstrumentProfile>,
}

/// Options of [`Solver::check_at`] beyond the request options. `None`: the default.
#[derive(Debug, Clone, Copy, Default)]
pub struct AtCheck {
    /// 1σ of the WCS at the positions, arcsec (e.g. a solve's `worst_sigma_arcsec`; default 0).
    pub wcs_sigma_arcsec: Option<f64>,
    /// Stars farther than this are not reported, arcsec (default 30).
    pub search_radius_arcsec: Option<f64>,
    /// Stars reported per position (default 3).
    pub max_stars: Option<usize>,
}

/// A request ready to solve, with the priors actually used (report field `priors`).
#[derive(Debug, Clone)]
pub struct Built {
    /// The request.
    pub request: SolveRequest,
    /// The priors used and where each came from.
    pub priors: Value,
    /// Warnings: header values that were not used (e.g. an unsupported reference frame, a
    /// clock error) and options dropped by `no_priors` or `no_pointing`.
    pub warnings: Vec<String>,
}

/// Build the request from the sources, the source list's image size (if it has one) and the
/// options.
pub fn build_request(
    sources: Vec<Source>,
    list_size: (Option<u32>, Option<u32>),
    o: &RequestOptions,
) -> Result<Built, String> {
    let header = o.header.as_ref();
    let no_header_wcs =
        o.no_header_wcs || o.instrument.as_ref().and_then(|p| p.header_wcs) == Some(false);
    let hp = header.map(|h| {
        if no_header_wcs {
            header_priors(&h.without_wcs())
        } else {
            header_priors(h)
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
    // A prior given together with the switch that drops it: the switch wins, as it always has
    // (a caller can keep its arguments and add the switch), but the report says so.
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
    // no_pointing drops the pointing only (header or caller); the instrument and caller scale,
    // roll and parity stay.
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
    // Hard limit and 1σ of the pointing: caller, else instrument profile, else 2° and a quarter
    // of the limit. The two are independent (a good hint with a long tail wants a small σ and a
    // wide limit).
    let profile = o.instrument.as_ref();
    let radius = o
        .radius
        .or(profile.and_then(|p| p.pointing_radius_deg))
        .unwrap_or(2.0);
    let pointing_sigma = o
        .pointing_sigma
        .or(profile.and_then(|p| p.pointing_sigma_deg))
        .unwrap_or(radius / 4.0);
    // Scale: caller, else instrument profile (with its σ), else header.
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
    // Roll: caller, else instrument rule. roll_sigma and roll_max override either source.
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

/// A catalogue and indexes, opened once and reused across frames. `Sync`: one `Solver` may serve
/// several threads.
pub struct Solver {
    catalog: Ucac4Catalog,
    indexes: Vec<Index>,
    index_paths: Vec<String>,
}

/// A report and, when solved or refitted, the solution behind it.
pub struct Report {
    /// The JSON report, as `fayde solve|verify|refit` prints it.
    pub json: Value,
    /// The solution, for writing a WCS file.
    pub solution: Option<Box<Solution>>,
    /// `true` when the solver failed (invalid request, catalogue error); `json` then has
    /// `"error": true` and `reason`.
    pub error: bool,
}

impl Solver {
    /// Open UCAC4 at `ucac4_dir` and read the index files.
    /// Relative paths are made absolute here, so a long-lived `Solver` does not depend on the
    /// process's working directory later.
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

    /// The catalogue.
    pub fn catalog(&self) -> &Ucac4Catalog {
        &self.catalog
    }

    /// Solve with the indexes (and the pointing / blind search the planner chooses).
    pub fn solve(&self, b: &Built) -> Report {
        let out = solve_auto(&b.request, &self.catalog, &self.indexes);
        self.outcome(out, b, true)
    }

    /// Re-fit a given WCS against the catalogue.
    pub fn refit(&self, b: &Built, wcs: &Wcs) -> Report {
        let out = refine_wcs(&b.request, &self.catalog, wcs);
        self.outcome(out, b, false)
    }

    /// Check a given WCS against the catalogue without changing it.
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
        json["identity"] = self.identity(req, false);
        Report {
            json: with_warnings(json, b),
            solution: None,
            error: false,
        }
    }

    /// "Is there a catalogue star here?" for each position under `wcs` (`fayde verify --at`).
    /// The epoch comes from `o` (caller, else header), the default pixel σ from
    /// `o.pixel_error` (default 1 px). Pixel fields in the report are FITS 1-based.
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
            "identity": self.identity_at(epoch),
        });
        Report {
            json,
            solution: None,
            error: false,
        }
    }

    fn at_failed(&self, reason: &str) -> Report {
        Report {
            json: json!({"error": true, "reason": reason, "identity": self.identity_at(None)}),
            solution: None,
            error: true,
        }
    }

    fn identity_at(&self, epoch: Option<Epoch>) -> Value {
        let mut req_for_id = SolveRequest {
            frame: Frame {
                width: 1,
                height: 1,
                trail: None,
            },
            sources: Vec::new(),
            priors: Priors::default(),
            policy: Policy::default(),
        };
        req_for_id.policy.epoch = epoch;
        self.identity(&req_for_id, false)
    }

    fn failed(&self, reason: &str, b: &Built, with_indexes: bool) -> Report {
        Report {
            json: with_warnings(
                json!({
                    "solved": false,
                    "error": true,
                    "reason": reason,
                    "priors": b.priors,
                    "identity": self.identity(&b.request, with_indexes),
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
                    // The stage that produced the solution, and whether it was the relaxed pass.
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
        json["identity"] = self.identity(req, with_indexes);
        Report {
            json: with_warnings(json, b),
            solution,
            error: false,
        }
    }

    /// What produced a report: Fayde version and commit, catalogue, epochs, indexes.
    fn identity(&self, req: &SolveRequest, with_indexes: bool) -> Value {
        let id = self.catalog.id();
        let n = if with_indexes { self.indexes.len() } else { 0 };
        json!({
            "fayde": {"version": env!("CARGO_PKG_VERSION"), "commit": commit()},
            "catalog": {
                "name": id.name,
                "release": id.release,
                "location": id.location,
                "reference_epoch": "J2000.0",
                "solve_epoch_julian_year": req.policy.epoch.map(|e| e.year()),
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

/// Warnings (header values not used, options dropped by a switch) go in the report as
/// `warnings`.
fn with_warnings(mut json: Value, b: &Built) -> Value {
    if !b.warnings.is_empty() {
        json["warnings"] = b.warnings.clone().into();
    }
    json
}

/// Detection settings a caller can give; `None`: the default, or the image header's.
#[derive(Debug, Clone, Copy, Default)]
pub struct DetectOptions {
    /// Stellar FWHM, px; default: measured from the image.
    pub fwhm: Option<f64>,
    /// Peak threshold of the filtered image, σ; default: the level above which one noise peak
    /// per frame is expected (about 5 for a 4-megapixel frame).
    pub threshold: Option<f64>,
    /// Saturation level, image units; default: the header's `SATURATE`, else the level of
    /// flat-topped stars, if any.
    pub saturation: Option<f64>,
    /// Gain, e⁻ per image unit; default: the header's `GAIN` or `EGAIN`, else none (the stars'
    /// own photon noise is then not in their position errors).
    pub gain: Option<f64>,
    /// Most sources kept, brightest first (default 2000).
    pub max_sources: Option<usize>,
}

/// Sources kept from an image when the caller sets no limit.
pub const DEFAULT_MAX_SOURCES: usize = 2000;

/// Full width at half maximum of a Gaussian over its σ: 2√(2 ln 2).
pub const FWHM_PER_SIGMA: f64 = 2.354_820_045_030_949;

/// Sources detected in `image`, ready for [`build_request`], and the detection report (JSON):
/// settings used, background and noise, counts by class, streaks, the sources (`sources`, in
/// the order handed to the solver) and what was set aside (`rejected`), with FITS 1-based
/// pixels.
///
/// The sources handed to the solver are the stars (classes `star`, `saturated`, `extended`),
/// each with its flux, 1σ position error and star probability. Hot pixels, cosmic-ray hits and
/// everything on a traced streak are listed in the report with their class and left out of the
/// sources (`docs/design/detection.md`).
pub fn detect_sources(
    image: &fayde_fits::Image,
    o: &DetectOptions,
) -> Result<(Vec<Source>, Value), String> {
    use fayde_detect::Class;
    // Numbers, or numbers written as strings (some cameras write GAIN = '1.10').
    let hdr = |keys: &[&str]| {
        keys.iter().find_map(|k| {
            image
                .header
                .get(k)
                .and_then(|v| v.as_f64().or_else(|| v.as_str()?.trim().parse().ok()))
                .filter(|v: &f64| v.is_finite() && *v > 0.0)
        })
    };
    // The caller's settings are checked here, under the names the caller knows.
    let positive = |v: Option<f64>| v.is_none_or(|x| x.is_finite() && x > 0.0);
    if o.fwhm.is_some_and(|f| !(0.5..=50.0).contains(&f)) {
        return Err("fwhm must be in [0.5, 50] px".into());
    }
    if !positive(o.threshold) {
        return Err("the detection threshold must be > 0".into());
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
    let fits = |x: f64, y: f64| PixelConvention::Fits.from_internal(PixelPos::new(x, y));
    // An infinite major axis (a ridge) has no JSON number.
    let finite = |v: f64| v.is_finite().then_some(v);
    let row = |d: &fayde_detect::Detected| -> Value {
        let (x, y) = fits(d.x, d.y);
        // Mean of the two axes; a ridge's endless major axis is capped.
        let major = d.sigma_major.min(10.0 * d.sigma_minor.max(0.1));
        let mean2 = 0.5 * (d.sigma_minor * d.sigma_minor + major * major);
        json!({
            "x": x, "y": y, "sigma_px": d.sigma_px, "flux": d.flux, "snr": d.snr,
            "p_star": d.p_star, "class": d.class.name(),
            "fwhm": FWHM_PER_SIGMA * mean2.sqrt(),
            "sigma_major_px": finite(d.sigma_major), "sigma_minor_px": d.sigma_minor,
            "theta_deg": d.theta.to_degrees(), "streak": d.streak,
            "truncated": d.truncated, "saturated": d.class == Class::Saturated,
        })
    };
    // `sources` is the list handed to the solver, in its order: `pairs[].source` of a solve
    // report indexes it. What was set aside is listed apart.
    let list: Vec<Value> = r.stars().map(row).collect();
    let rejected: Vec<Value> = r
        .sources
        .iter()
        .filter(|d| !d.class.is_stellar())
        .map(row)
        .collect();
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
    let report = json!({
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
        },
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
            })
        }).collect::<Vec<_>>(),
        "warnings": warnings,
        "sources": list,
        "rejected": rejected,
    });
    Ok((sources, report))
}

/// A FITS header from a JSON object of keyword → value (numbers, strings, booleans; `null`
/// skipped). Keys are upper-cased. Numbers given as strings (e.g. `"JD-OBS": "2459175.3"`) stay
/// strings; [`fayde_fits::header_priors`] parses those it needs.
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

/// The source commit this library was built from (`git describe`-style short hash).
pub fn commit() -> &'static str {
    env!("FAYDE_COMMIT")
}

/// WCS cards of a solution plus `IMAGEW`/`IMAGEH`.
pub fn solution_cards(sol: &Solution, req: &SolveRequest) -> Vec<String> {
    let mut cards = sol.wcs.to_cards();
    cards.push(format!("IMAGEW  = {:>20}", req.frame.width));
    cards.push(format!("IMAGEH  = {:>20}", req.frame.height));
    cards
}

/// The fit's parameter covariance: per tangent-plane axis (ξ and η share it), rad², over the
/// monomials `(u/norm)^p (v/norm)^q` with `(u, v)` the pixel offset from `crpix`.
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

/// Matched pairs, in input-table order, with source positions in the FITS convention (1-based).
fn pairs_json(matches: &[Match], req: &SolveRequest) -> Value {
    matches
        .iter()
        .map(|m| {
            let s = &req.sources[m.source];
            let (x, y) = PixelConvention::Fits.from_internal(PixelPos::new(s.x, s.y));
            let (ra, dec) = m.star_pos.to_degrees();
            // Where the WCS puts the star: the source position minus the residual.
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
        assert_eq!(b.warnings.len(), 2, "{:?}", b.warnings);
        assert!(b
            .warnings
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
        assert_eq!(b.warnings, ["ra/dec was given but dropped by no-pointing"]);
        // The switches alone say nothing.
        let quiet = RequestOptions {
            no_priors: true,
            ..opts()
        };
        assert!(build_request(vec![], (None, None), &quiet)
            .unwrap()
            .warnings
            .is_empty());
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

    /// The pointing's limit and 1σ come from the caller, else the profile, else the defaults;
    /// and a header WCS is left unread when the caller or the profile says so.
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
        // The header WCS is the pointing and the scale by default.
        let b = with(RequestOptions::default());
        let pt = b.request.priors.pointing.unwrap();
        assert!((pt.centre.to_degrees().0 - 10.0).abs() < 1e-9);
        assert!((pt.max_radius.to_degrees() - 2.0).abs() < 1e-12);
        assert!((pt.sigma.to_degrees() - 0.5).abs() < 1e-12);
        assert_eq!(b.priors["pointing"]["origin"], "CRVAL1/CRVAL2");
        assert!((b.request.priors.scale.unwrap().arcsec_per_px - 3.6).abs() < 1e-9);
        // Without it: the pointing keywords, and no scale from the header.
        let b = with(RequestOptions {
            no_header_wcs: true,
            ..RequestOptions::default()
        });
        let pt = b.request.priors.pointing.unwrap();
        assert!((pt.centre.to_degrees().0 - 11.0).abs() < 1e-9);
        assert!((pt.centre.to_degrees().1 - 21.5).abs() < 1e-9);
        assert_eq!(b.priors["pointing"]["origin"], "RA/DEC");
        assert!(b.request.priors.scale.is_none());
        // A profile can ask for the same, and gives the pointing's 1σ and limit.
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
        // The caller's values come before the profile's.
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
