#![forbid(unsafe_code)]

mod error;
pub mod memory;
pub mod propagate;
pub mod ucac4;

use fayde_core::{Epoch, SkyPos};

pub use error::CatalogError;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CatalogId {
    pub name: &'static str,
    pub release: &'static str,
    pub location: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Band {
    #[default]
    Default,
    UcacModel,
    UcacAperture,
    ApassB,
    ApassV,
    ApassG,
    ApassR,
    ApassI,
    TwoMassJ,
    TwoMassH,
    TwoMassK,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StarQuery {
    pub epoch: Option<Epoch>,
    pub mag_max: Option<f32>,
    pub band: Band,
    pub require_proper_motion: bool,
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

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Star {
    pub id: u64,
    pub pos: SkyPos,
    pub mag: Option<f32>,
    pub pos_err_mas: [f32; 2],
    pub pm_mas_yr: [Option<f32>; 2],
}

pub trait Catalog {
    fn id(&self) -> CatalogId;

    /// Return stars within the angular radius (radians) at the query epoch.
    fn cone(
        &self,
        centre: SkyPos,
        radius: f64,
        query: &StarQuery,
    ) -> Result<Vec<Star>, CatalogError>;

    fn cone_brightest(
        &self,
        centre: SkyPos,
        radius: f64,
        query: &StarQuery,
        n: usize,
    ) -> Result<Vec<Star>, CatalogError> {
        Ok(brightest_first(self.cone(centre, radius, query)?, n))
    }

    fn scan(&self, query: &StarQuery, visit: &mut dyn FnMut(Star)) -> Result<(), CatalogError>;
}

pub fn brightest_first(stars: Vec<Star>, n: usize) -> Vec<Star> {
    let mut ranked: Vec<(f32, Star)> = stars
        .into_iter()
        .filter_map(|s| s.mag.filter(|m| m.is_finite()).map(|m| (m, s)))
        .collect();
    ranked.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.id.cmp(&b.1.id)));
    ranked.truncate(n);
    ranked.into_iter().map(|(_, s)| s).collect()
}
