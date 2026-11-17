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

#[cfg(test)]
mod tests {
    use super::*;

    /// The default `cone_brightest`: finite magnitudes only, brightest first, equal magnitudes
    /// by identifier and then in catalogue order, cut to `n`.
    #[test]
    fn brightest_are_ranked_by_magnitude_then_identifier() {
        let at = |ra: f64| SkyPos::from_degrees(ra, 0.0).unwrap();
        let star = |id: u64, ra: f64, mag: Option<f32>| Star {
            id,
            pos: at(ra),
            mag,
            pos_err_mas: [0.0; 2],
            pm_mas_yr: [None; 2],
        };
        let cat = MemoryCatalog::new(vec![
            star(5, 10.01, Some(9.0)),
            star(3, 10.02, Some(7.5)),
            star(9, 10.03, None),
            star(2, 10.04, Some(9.0)),
            star(7, 10.05, Some(f32::NAN)),
            star(2, 10.06, Some(9.0)), // the identifier of the fourth, listed after it
            star(1, 40.0, Some(1.0)),  // outside the cone
        ]);
        let got = |n| {
            cat.cone_brightest(at(10.0), 1f64.to_radians(), &StarQuery::default(), n)
                .unwrap()
                .iter()
                .map(|s| (s.id, s.pos.to_degrees().0))
                .collect::<Vec<_>>()
        };
        let all = got(usize::MAX);
        assert_eq!(all.iter().map(|s| s.0).collect::<Vec<_>>(), [3, 2, 2, 5]);
        assert!(all[1].1 < all[2].1, "equal stars keep catalogue order");
        assert_eq!(got(2), all[..2]);
        assert!(got(0).is_empty());
    }
}
