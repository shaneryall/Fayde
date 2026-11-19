use fayde_core::Portable;
use std::time::Duration;

use fayde_catalog::CatalogId;
use fayde_core::{Epoch, SkyPos};
use fayde_wcs::{Parity, Wcs};

use crate::fit::WcsCovariance;

/// Source-list entry with a centroid in zero-based image pixels.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Source {
    pub x: f64,
    pub y: f64,
    pub flux: Option<f64>,
    pub sigma_px: Option<f64>,
    pub p_star: Option<f64>,
    /// Mark a clipped centroid for verification only. Seeding and fitting skip it.
    pub truncated: bool,
}

impl Source {
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

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Frame {
    pub width: u32,
    pub height: u32,
    pub trail: Option<Trail>,
}

impl Frame {
    pub fn new(width: u32, height: u32) -> Self {
        Self {
            width,
            height,
            trail: None,
        }
    }
}

/// Common stellar trail with orientation in radians from +x toward +y.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Trail {
    pub length_px: f64,
    pub angle: f64,
}

/// Prior on the field centre. Both the one-sigma uncertainty and the search radius are in radians.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PointingPrior {
    pub centre: SkyPos,
    pub sigma: f64,
    pub max_radius: f64,
}

/// Pixel scale and its one-sigma uncertainty, both in arcseconds per pixel.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ScalePrior {
    pub arcsec_per_px: f64,
    pub sigma: f64,
}

/// Position angle of the image +y axis, measured east of north. All fields are in radians.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RollPrior {
    pub angle: f64,
    pub sigma: f64,
    pub max_dev: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Priors {
    pub pointing: Option<PointingPrior>,
    pub scale: Option<ScalePrior>,
    pub roll: Option<RollPrior>,
    pub parity: Option<Parity>,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Policy {
    pub accept_log_odds: f64,
    pub max_hypotheses: u64,
    pub deadline: Option<Duration>,
    pub seed_sources: usize,
    pub default_sigma_px: f64,
    pub distractor: f64,
    pub epoch: Option<Epoch>,
    pub max_sip_order: usize,
    pub pattern_sources: usize,
    pub search_outside_pointing: bool,
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

#[derive(Debug, Clone, PartialEq)]
pub struct SolveRequest {
    pub frame: Frame,
    pub sources: Vec<Source>,
    pub priors: Priors,
    pub policy: Policy,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Match {
    /// Index into the original [`SolveRequest::sources`] list.
    pub source: usize,
    pub star_id: u64,
    pub star_pos: SkyPos,
    pub star_mag: Option<f32>,
    pub residual_px: [f64; 2],
    pub inlier_prob: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolveMethod {
    PairSeeds,
    PatternHash,
    Refit,
}

#[derive(Debug, Clone, PartialEq)]
pub struct StageTrace {
    pub stage: String,
    pub hypotheses: u64,
    pub seconds: f64,
    pub outcome: &'static str,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Solution {
    pub method: SolveMethod,
    pub search: Vec<StageTrace>,
    pub prior_violated: bool,
    pub violated_priors: Vec<&'static str>,
    pub index_hash: Option<u64>,
    pub wcs: Wcs,
    pub covariance: WcsCovariance,
    pub log_odds: f64,
    /// Model-based false-alarm score, `hypotheses_tested × exp(−log_odds)`.
    pub nfa: f64,
    pub hypotheses_tested: u64,
    pub matches: Vec<Match>,
    pub rms_px: f64,
    pub error_scale: f64,
    pub worst_sigma_arcsec: f64,
    pub loo_rms_px: Option<f64>,
    pub catalog: CatalogId,
    pub elapsed: Duration,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Verification {
    pub log_odds: f64,
    pub accepted: bool,
    pub matches: Vec<Match>,
    pub rms_px: f64,
    pub references: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub enum Outcome {
    Solved(Box<Solution>),
    Abstained {
        reason: String,
        best_log_odds: Option<f64>,
        hypotheses_tested: u64,
        search: Vec<StageTrace>,
    },
}

#[derive(Debug)]
pub enum SolveError {
    InvalidInput(String),
    Unsupported(String),
    Catalog(fayde_catalog::CatalogError),
}

impl std::fmt::Display for SolveError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidInput(s) => write!(f, "invalid solve request: {s}"),
            Self::Unsupported(s) => write!(f, "unsupported: {s}"),
            Self::Catalog(e) => write!(f, "{e}"),
        }
    }
}

impl std::error::Error for SolveError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Catalog(e) => Some(e),
            _ => None,
        }
    }
}

impl From<fayde_catalog::CatalogError> for SolveError {
    fn from(e: fayde_catalog::CatalogError) -> Self {
        Self::Catalog(e)
    }
}
