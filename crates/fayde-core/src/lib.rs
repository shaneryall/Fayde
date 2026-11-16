//! Core types shared by every Fayde crate.
//!
//! Conventions (see `docs/design/plate-solver.md` §2):
//! - Angles are stored in **radians** as `f64`. Use [`deg`] / [`arcsec`] to build them.
//! - Sky positions are ICRS right ascension / declination.
//! - Pixel positions are `(x, y)` = (column, row) with the **centre of the first pixel at
//!   `(0, 0)`** (NumPy / Astropy origin 0). Other conventions are converted only at I/O edges,
//!   via [`PixelConvention`].

#![forbid(unsafe_code)]

pub mod epoch;
/// Internal: shared by the Fayde crates for atomic file writes. Not part of the public API; it
/// may change or go away in any release.
#[doc(hidden)]
pub mod fsutil;
pub mod math;
pub mod pixel;
pub mod sky;

pub use epoch::Epoch;
pub use math::Portable;
pub use pixel::{PixelConvention, PixelPos};
pub use sky::{SkyPos, UnitVector};

/// Degrees → radians.
#[inline]
pub fn deg(value: f64) -> f64 {
    value.to_radians()
}

/// Arcseconds → radians.
#[inline]
pub fn arcsec(value: f64) -> f64 {
    (value / 3600.0).to_radians()
}

/// Radians → arcseconds.
#[inline]
pub fn to_arcsec(radians: f64) -> f64 {
    radians.to_degrees() * 3600.0
}
