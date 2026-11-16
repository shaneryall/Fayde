//! Star catalogue access for Fayde.
//!
//! [`Catalog`] is the interface the index builder and the solver use. [`ucac4::Ucac4Catalog`]
//! reads the original UCAC4 distribution in place. All positions returned are propagated to the
//! requested epoch ([`propagate`]), with their uncertainty grown accordingly.

#![forbid(unsafe_code)]

mod error;
pub mod memory;
pub mod propagate;
pub mod ucac4;

use fayde_core::{Epoch, SkyPos};

pub use error::CatalogError;

/// Identifies a catalogue (and its local copy) in solution provenance.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CatalogId {
    /// Catalogue name, e.g. `"UCAC4"`.
    pub name: &'static str,
    /// Release, e.g. `"2012"`.
    pub release: &'static str,
    /// Where it was read from.
    pub location: String,
}

/// Photometric band to report as [`Star::mag`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Band {
    /// Catalogue default. For UCAC4: model magnitude (579–642 nm), falling back to aperture.
    #[default]
    Default,
    /// UCAC4 model magnitude only.
    UcacModel,
    /// UCAC4 aperture magnitude only.
    UcacAperture,
    /// APASS B.
    ApassB,
    /// APASS V.
    ApassV,
    /// APASS g′.
    ApassG,
    /// APASS r′.
    ApassR,
    /// APASS i′.
    ApassI,
    /// 2MASS J.
    TwoMassJ,
    /// 2MASS H.
    TwoMassH,
    /// 2MASS Ks.
    TwoMassK,
}

/// What to fetch.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StarQuery {
    /// Propagate positions to this epoch. `None` returns catalogue-epoch positions.
    pub epoch: Option<Epoch>,
    /// Keep only stars at least this bright in `band`. Stars without a magnitude in `band` are
    /// dropped when this is set.
    pub mag_max: Option<f32>,
    /// Band reported in [`Star::mag`] and used by `mag_max`.
    pub band: Band,
    /// Keep only stars with a measured proper motion in both components.
    pub require_proper_motion: bool,
    /// When `true`, a missing or unreadable catalogue part covering the cone is an error. When
    /// `false` it is skipped. Callers that read an empty result as "no star here" must use `true`.
    pub strict: bool,
}

impl Default for StarQuery {
    fn default() -> Self {
        Self {
            epoch: None,
            mag_max: None,
            band: Band::Default,
            require_proper_motion: false,
            strict: true,
        }
    }
}

/// One catalogue star at the query epoch.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Star {
    /// Catalogue-unique identifier (UCAC4: `rnm`).
    pub id: u64,
    /// Position at the query epoch (or catalogue epoch if none was given).
    pub pos: SkyPos,
    /// Magnitude in the requested band, if the catalogue has one.
    pub mag: Option<f32>,
    /// 1σ position error at the query epoch, mas, as `[RA·cos δ, Dec]`.
    pub pos_err_mas: [f32; 2],
    /// Proper motion, mas/yr, as `[RA·cos δ, Dec]`. A component is `None` when the catalogue
    /// has no measured value; the star then keeps its catalogue-epoch coordinate in that
    /// component and its error grows with a default proper-motion error.
    pub pm_mas_yr: [Option<f32>; 2],
}

/// Read access to a star catalogue.
pub trait Catalog {
    /// Identity for provenance.
    fn id(&self) -> CatalogId;

    /// All stars within `radius` (radians, true angular separation) of `centre`, at the query
    /// epoch, sorted by catalogue order (deterministic).
    fn cone(
        &self,
        centre: SkyPos,
        radius: f64,
        query: &StarQuery,
    ) -> Result<Vec<Star>, CatalogError>;

    /// Visit every star in the catalogue at the query epoch, in catalogue order. Used to build
    /// all-sky indexes without holding the whole catalogue in memory.
    fn scan(&self, query: &StarQuery, visit: &mut dyn FnMut(Star)) -> Result<(), CatalogError>;
}
