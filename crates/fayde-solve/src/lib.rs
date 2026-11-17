//! Plate solving for Fayde.
//!
//! [`solve_auto`] plans the search from the request's metadata (pair seeds, ringed region pattern
//! hash, scale-pruned or blind pattern hash across index bands), verifies hypotheses by log-odds
//! with prior terms and a false-alarm count, and fits TAN or TAN-SIP (order by BIC) with a
//! reported covariance. Sequence mode is not implemented yet. Algorithm sources are listed in
//! `docs/references.md`.

#![forbid(unsafe_code)]

mod at;
mod fit;
mod grid;
mod hash;
mod refset;
mod request;
mod solve;
mod verify;

pub use at::{stars_at, AtOptions, AtPosition, AtResult, AtStar, CONSISTENT_P};
pub use fit::WcsCovariance;
pub use request::{
    Frame, Match, Outcome, PointingPrior, Policy, Priors, RollPrior, ScalePrior, Solution,
    SolveError, SolveMethod, SolveRequest, Source, StageTrace, Trail, Verification,
};
pub use solve::{refine_wcs, solve, solve_auto, solve_with_index, solve_with_indexes, verify_wcs};
