use fayde_core::SkyPos;

use crate::{Catalog, CatalogError, CatalogId, Star, StarQuery};

#[derive(Debug, Clone, Default)]
pub struct MemoryCatalog {
    stars: Vec<Star>,
}

impl MemoryCatalog {
    pub fn new(stars: Vec<Star>) -> Self {
        Self { stars }
    }

    pub fn stars(&self) -> &[Star] {
        &self.stars
    }
}

fn matches_query(star: &Star, query: &StarQuery) -> bool {
    query
        .mag_max
        .is_none_or(|limit| star.mag.is_some_and(|m| m <= limit))
        && (!query.require_proper_motion || star.pm_mas_yr.iter().all(Option::is_some))
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
            .filter(|s| matches_query(s, query) && centre.separation(s.pos) <= radius)
            .copied()
            .collect())
    }

    fn scan(&self, query: &StarQuery, visit: &mut dyn FnMut(Star)) -> Result<(), CatalogError> {
        for s in &self.stars {
            if matches_query(s, query) {
                visit(*s);
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scan_and_cone_require_both_proper_motion_components() {
        let pos = SkyPos::from_degrees(10.0, 0.0).unwrap();
        let star = Star {
            id: 1,
            pos,
            mag: Some(9.0),
            pos_err_mas: [0.0; 2],
            pm_mas_yr: [Some(0.0), Some(0.0)],
        };
        let cat = MemoryCatalog::new(vec![
            star,
            Star {
                id: 2,
                pm_mas_yr: [Some(0.0), None],
                ..star
            },
        ]);
        let query = StarQuery {
            require_proper_motion: true,
            ..StarQuery::default()
        };
        let mut scanned = Vec::new();
        cat.scan(&query, &mut |s| scanned.push(s)).unwrap();
        assert_eq!(scanned, [star]);
        assert_eq!(cat.cone(pos, 0.01, &query).unwrap(), scanned);
    }

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
            star(2, 10.06, Some(9.0)),
            star(1, 40.0, Some(1.0)),
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
