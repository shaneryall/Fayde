#![forbid(unsafe_code)]

pub mod epoch;
#[doc(hidden)]
pub mod fsutil;
pub mod math;
pub mod pixel;
pub mod sky;

pub use epoch::Epoch;
pub use math::Portable;
pub use pixel::{PixelConvention, PixelPos};
pub use sky::{SkyPos, UnitVector};

#[inline]
pub fn deg(value: f64) -> f64 {
    value.to_radians()
}

#[inline]
pub fn arcsec(value: f64) -> f64 {
    (value / 3600.0).to_radians()
}

#[inline]
pub fn to_arcsec(radians: f64) -> f64 {
    radians.to_degrees() * 3600.0
}
