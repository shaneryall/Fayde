//! Plate solving from image centroids and a star catalogue. Start with [`solve_auto`].

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
pub use solve::{
    refine_wcs, solve, solve_auto, solve_auto_after, solve_with_index, solve_with_indexes,
    verify_wcs,
};
