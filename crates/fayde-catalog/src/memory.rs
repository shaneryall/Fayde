//! An in-memory catalogue: a list of stars with brute-force cone search.
//!
//! For tests, synthetic fields, and small user-supplied reference lists. Positions are treated
//! as already at the query epoch (no propagation).

use fayde_core::SkyPos;

use crate::{Catalog, CatalogError, CatalogId, Star, StarQuery};

/// Stars held in memory.
#[derive(Debug, Clone, Default)]
pub struct MemoryCatalog {
    stars: Vec<Star>,
}

impl MemoryCatalog {
    /// From a list of stars.
    pub fn new(stars: Vec<Star>) -> Self {
        Self { stars }
    }

    /// All stars.
    pub fn stars(&self) -> &[Star] {
        &self.stars
    }
}

impl Catalog for MemoryCatalog {
    fn id(&self) -> CatalogId {
        CatalogId {
            name: "memory",
            release: "",
            location: format!("{} stars", self.stars.len()),
        }
    }

    fn cone(
        &self,
        centre: SkyPos,
        radius: f64,
        query: &StarQuery,
    ) -> Result<Vec<Star>, CatalogError> {
        if !radius.is_finite() || radius < 0.0 {
            return Err(CatalogError::InvalidQuery(format!("radius {radius} rad")));
        }
        Ok(self
            .stars
            .iter()
            .filter(|s| centre.separation(s.pos) <= radius)
            .filter(|s| match query.mag_max {
                Some(lim) => s.mag.is_some_and(|m| m <= lim),
                None => true,
            })
            .filter(|s| {
                !query.require_proper_motion
                    || (s.pm_mas_yr[0].is_some() && s.pm_mas_yr[1].is_some())
            })
            .copied()
            .collect())
    }

    fn scan(&self, query: &StarQuery, visit: &mut dyn FnMut(Star)) -> Result<(), CatalogError> {
        for s in &self.stars {
            if query
                .mag_max
                .is_some_and(|lim| !s.mag.is_some_and(|m| m <= lim))
            {
                continue;
            }
            visit(*s);
        }
        Ok(())
    }
}
