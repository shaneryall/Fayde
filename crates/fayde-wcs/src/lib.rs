//! World coordinate systems for Fayde: the gnomonic (TAN) projection and FITS WCS header cards.
//!
//! Pixel coordinates follow the Fayde convention (first pixel centre at `(0, 0)`). FITS
//! `CRPIX` values are 1-based; [`TanWcs::from_header`] and [`TanWcs::to_cards`] convert at that
//! boundary and nowhere else.

#![forbid(unsafe_code)]

mod error;
pub mod header;
pub mod sip;
mod tan;

pub use error::WcsError;
pub use header::{Header, Value};
pub use sip::{Poly, Sip, Wcs};
pub use tan::{gnomonic, gnomonic_inverse, Parity, TanWcs};
