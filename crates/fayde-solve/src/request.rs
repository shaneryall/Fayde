//! Solve inputs and outputs (a subset of the input contract, design doc §6).

use fayde_core::Portable;
use std::time::Duration;

use fayde_catalog::CatalogId;
use fayde_core::{Epoch, SkyPos};
use fayde_wcs::{Parity, Wcs};

use crate::fit::WcsCovariance;

/// One detected source (source-list row). Pixel convention: first pixel centre at `(0, 0)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Source {
    /// Column.
    pub x: f64,
    /// Row.
    pub y: f64,
    /// Integrated flux (any linear unit). Used for ordering.
    pub flux: Option<f64>,
    /// 1σ centroid error, pixels (isotropic).
    pub sigma_px: Option<f64>,
    /// Prior probability that this is a real star, 0–1.
    pub p_star: Option<f64>,
    /// The source is cut by the frame edge (or a mask), so its centroid is biased along the
    /// trail direction. Truncated sources still count as evidence for a match, with a gate
    /// stretched along [`Frame::trail`], but never enter the fit and never seed hypotheses.
    pub truncated: bool,
}

impl Source {
    /// A source with only a position.
    pub fn at(x: f64, y: f64) -> Self {
        Self {
            x,
            y,
            flux: None,
            sigma_px: None,
            p_star: None,
            truncated: false,
        }
    }
}

/// Image size, and the common star trail if the frame tracked something other than the stars.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Frame {
    /// Width, pixels.
    pub width: u32,
    /// Height, pixels.
    pub height: u32,
    /// Star trail shared by every source (track-rate imagery), if known.
    pub trail: Option<Trail>,
}

impl Frame {
    /// A frame of `width × height` pixels with point-like stars.
    pub fn new(width: u32, height: u32) -> Self {
        Self {
            width,
            height,
            trail: None,
        }
    }
}

/// The star trail common to a frame.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Trail {
    /// Full trail length, pixels.
    pub length_px: f64,
    /// Direction of the trail, radians from +x towards +y. Its sign does not matter.
    pub angle: f64,
}

/// Prior on the field centre.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PointingPrior {
    /// Expected centre.
    pub centre: SkyPos,
    /// 1σ, radians.
    pub sigma: f64,
    /// Hard limit, radians: no solution is accepted with its centre farther than this.
    pub max_radius: f64,
}

/// Prior on the pixel scale.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ScalePrior {
    /// Expected scale, arcsec/px.
    pub arcsec_per_px: f64,
    /// 1σ, arcsec/px.
    pub sigma: f64,
}

/// Prior on the roll: the position angle of the image +y axis, east of north, radians.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RollPrior {
    /// Expected roll.
    pub angle: f64,
    /// 1σ, radians (enters the hypothesis score).
    pub sigma: f64,
    /// Hard limit on |roll − angle|, radians.
    pub max_dev: f64,
}

/// All priors. Pair-seed (hinted) solving needs `pointing` and `scale`.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Priors {
    /// Pointing.
    pub pointing: Option<PointingPrior>,
    /// Scale.
    pub scale: Option<ScalePrior>,
    /// Roll.
    pub roll: Option<RollPrior>,
    /// Parity, if known.
    pub parity: Option<Parity>,
}

/// Solver policy.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Policy {
    /// Required log-odds after subtracting ln(hypotheses tested). Default ln(1e9).
    pub accept_log_odds: f64,
    /// Stop after this many hypotheses. Bounds the work of a search independently of the
    /// machine; a wide pointing radius costs hypotheses roughly in proportion to its area.
    pub max_hypotheses: u64,
    /// Wall-clock budget, from the start of the request (source ranking included). The search
    /// and a refit stop at it; the final fit of a hypothesis already accepted still completes.
    pub deadline: Option<Duration>,
    /// Number of top-ranked sources used to form seed pairs.
    pub seed_sources: usize,
    /// Centroid error for sources without `sigma_px`, pixels.
    pub default_sigma_px: f64,
    /// Prior fraction of sources that are not catalogue stars (distractors), 0–1.
    pub distractor: f64,
    /// Catalogue epoch (observation time). `None` uses catalogue-epoch positions.
    pub epoch: Option<Epoch>,
    /// Highest SIP order considered (1 = TAN only). The order used is chosen by BIC.
    pub max_sip_order: usize,
    /// Top-ranked sources used to form 4-star patterns for the pattern hash. Patterns are tried
    /// breadth-first (all from the 4 best, then those adding the 5th, …), so a larger value only
    /// costs time on frames the smaller one could not solve.
    pub pattern_sources: usize,
    /// After every stage allowed by the pointing prior has failed, also search the indexes
    /// outside the pointing limit. A solution found that way has `prior_violated = true`.
    /// Off by default: a hard limit is a hard limit unless the caller says otherwise.
    pub search_outside_pointing: bool,
    /// When every stage has failed, search once more with the hard limits relaxed: roll and
    /// parity dropped, the scale prior widened to at least 1%, and the pointing limit lifted (as
    /// `search_outside_pointing`). For priors that are usually right but sometimes stale or
    /// mis-calibrated: a wrong limit then costs time, not the solution. The second pass gets
    /// what is left of the deadline and hypothesis budget, the acceptance threshold counts the
    /// hypotheses of both passes, and a solution that breaks an original limit is reported in
    /// `Solution::violated_priors`. Off by default.
    pub relax_priors_on_failure: bool,
}

impl Default for Policy {
    fn default() -> Self {
        Self {
            accept_log_odds: 1e9f64.pln(),
            max_hypotheses: 5_000_000,
            deadline: None,
            seed_sources: 10,
            default_sigma_px: 1.0,
            distractor: 0.25,
            epoch: None,
            max_sip_order: 3,
            pattern_sources: 16,
            search_outside_pointing: false,
            relax_priors_on_failure: false,
        }
    }
}

/// A solve request.
#[derive(Debug, Clone, PartialEq)]
pub struct SolveRequest {
    /// Image size.
    pub frame: Frame,
    /// Detected sources (all of them; do not pre-filter).
    pub sources: Vec<Source>,
    /// Priors.
    pub priors: Priors,
    /// Policy.
    pub policy: Policy,
}

/// A source matched to a catalogue star in the accepted solution.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Match {
    /// Index into `SolveRequest::sources`.
    pub source: usize,
    /// Catalogue star identifier.
    pub star_id: u64,
    /// Catalogue star position at the solve epoch.
    pub star_pos: SkyPos,
    /// Catalogue magnitude of the star (the catalogue's default band), if any.
    pub star_mag: Option<f32>,
    /// Source − projected star, pixels.
    pub residual_px: [f64; 2],
    /// Posterior probability that this match is real.
    pub inlier_prob: f64,
}

/// How the accepted hypothesis was generated.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolveMethod {
    /// Star-pair seed in the pointing-prior cone.
    PairSeeds,
    /// 4-star pattern looked up in the index.
    PatternHash,
    /// Re-match and refit around an existing WCS ([`crate::refine_wcs`]).
    Refit,
}

/// One stage of the adaptive search plan, as executed.
#[derive(Debug, Clone, PartialEq)]
pub struct StageTrace {
    /// What the stage searched, e.g. `pair seeds r<2.00°` or `pattern hash 1.30° index r<0.50°`.
    pub stage: String,
    /// Hypotheses tested in this stage.
    pub hypotheses: u64,
    /// Wall-clock time of the stage.
    pub seconds: f64,
    /// How the stage ended: `"solved"`, `"exhausted"` (searched everything it covers),
    /// `"deadline"` or `"hypothesis limit"` (`Policy::max_hypotheses`).
    pub outcome: &'static str,
}

/// An accepted solution.
#[derive(Debug, Clone, PartialEq)]
pub struct Solution {
    /// Hypothesis generator that found it.
    pub method: SolveMethod,
    /// The search plan as executed, ending with the stage that solved.
    pub search: Vec<StageTrace>,
    /// True when the solution breaks a hard limit of the request's priors (only possible with
    /// `Policy::search_outside_pointing` or `Policy::relax_priors_on_failure`).
    pub prior_violated: bool,
    /// Which limits it breaks: any of `"pointing"`, `"scale"`, `"roll"`, `"parity"`.
    pub violated_priors: Vec<&'static str>,
    /// Content hash of the pattern index that found it (pattern-hash solves only).
    pub index_hash: Option<u64>,
    /// Fitted WCS (`CRPIX` at the frame centre); SIP distortion when BIC prefers it.
    pub wcs: Wcs,
    /// Parameter covariance (scaled by the error scale²).
    pub covariance: WcsCovariance,
    /// Log-odds of the accepted hypothesis (natural log).
    pub log_odds: f64,
    /// Number of false alarms: hypotheses tested × exp(−log-odds).
    pub nfa: f64,
    /// Hypotheses tested before acceptance.
    pub hypotheses_tested: u64,
    /// Matched sources.
    pub matches: Vec<Match>,
    /// Inlier-weighted RMS residual, pixels.
    pub rms_px: f64,
    /// Fitted scale factor on the supplied centroid errors (>1: errors under-reported).
    pub error_scale: f64,
    /// Worst predicted radial 1σ position error over a 5×5 grid covering the frame (corners
    /// included), arcsec: the accuracy where it is lowest, not at the centre.
    pub worst_sigma_arcsec: f64,
    /// Leave-one-out RMS residual of the confident matches, pixels: each match's residual from
    /// a fit without it. A fit-independent accuracy figure; `None` with too few matches.
    pub loo_rms_px: Option<f64>,
    /// Catalogue used.
    pub catalog: CatalogId,
    /// Wall-clock time.
    pub elapsed: Duration,
}

/// How well an existing WCS explains a source list ([`crate::verify_wcs`]).
#[derive(Debug, Clone, PartialEq)]
pub struct Verification {
    /// Data log-odds of the WCS against the null (no prior terms, nothing tested but this).
    pub log_odds: f64,
    /// Whether the log-odds clears `Policy::accept_log_odds`.
    pub accepted: bool,
    /// Matched sources.
    pub matches: Vec<Match>,
    /// Inlier-weighted RMS residual, pixels.
    pub rms_px: f64,
    /// Catalogue stars expected in the frame (the reference count used).
    pub references: usize,
}

/// Result of a solve.
#[derive(Debug, Clone, PartialEq)]
pub enum Outcome {
    /// A solution passed verification.
    Solved(Box<Solution>),
    /// No hypothesis passed. Not an error: this is a valid answer.
    Abstained {
        /// Why.
        reason: String,
        /// Best log-odds seen, if any hypothesis reached verification.
        best_log_odds: Option<f64>,
        /// Hypotheses tested.
        hypotheses_tested: u64,
        /// The search plan as executed.
        search: Vec<StageTrace>,
    },
}

/// Errors (bad input or catalogue failure, not "no solution").
#[derive(Debug)]
pub enum SolveError {
    /// The request is malformed.
    InvalidInput(String),
    /// The request needs a capability not built yet (e.g. blind solving).
    Unsupported(String),
    /// Catalogue failure.
    Catalog(fayde_catalog::CatalogError),
}

impl std::fmt::Display for SolveError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidInput(s) => write!(f, "invalid solve request: {s}"),
            Self::Unsupported(s) => write!(f, "unsupported: {s}"),
            Self::Catalog(e) => write!(f, "catalogue: {e}"),
        }
    }
}

impl std::error::Error for SolveError {}

impl From<fayde_catalog::CatalogError> for SolveError {
    fn from(e: fayde_catalog::CatalogError) -> Self {
        Self::Catalog(e)
    }
}
