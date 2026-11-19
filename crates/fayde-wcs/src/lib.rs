#![forbid(unsafe_code)]

mod error;
pub mod header;
pub mod sip;
mod tan;

pub use error::WcsError;
pub use header::{Header, Value};
pub use sip::{Poly, Sip, Wcs};
pub use tan::{gnomonic, gnomonic_inverse, Parity, TanWcs, TangentPlane};
