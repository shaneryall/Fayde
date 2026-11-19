#![forbid(unsafe_code)]

pub mod build;
mod file;
pub mod healpix;
pub mod kdtree;
pub mod pattern;

use serde::{Deserialize, Serialize};

pub use build::{
    build, build_from_catalog, build_from_catalog_with_progress, BuildError, BuildParams, BUILDER,
};
pub use file::IndexFileError;

use pattern::{
    angle, edge_shape, key_check, planar, shape_of, KeyWindow, PatternTable, Shape, Slot,
};

pub const FORMAT_VERSION: u32 = 2;

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct IndexMeta {
    pub format_version: u32,
    pub catalog_name: String,
    pub catalog_release: String,
    pub params: BuildParams,
    pub healpix_nside: u32,
    pub lattice_fields: usize,
    pub bins: u32,
    pub star_count: usize,
    pub pattern_count: usize,
    pub table_size: usize,
    pub builder: String,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IndexStar {
    pub vec: [f64; 3],
    pub mag: f32,
    pub id: u64,
    pub pattern: bool,
}

#[derive(Debug, Clone)]
pub struct Index {
    meta: IndexMeta,
    stars: Vec<IndexStar>,
    table: PatternTable,
    content_hash: u64,
    tree: std::sync::OnceLock<kdtree::KdTree<3>>,
}

impl PartialEq for Index {
    fn eq(&self, other: &Self) -> bool {
        self.content_hash == other.content_hash
            && self.meta == other.meta
            && self.stars == other.stars
            && self.table == other.table
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PatternMatch {
    pub stars: [u32; 4],
    pub ratio_error: f64,
    pub largest: f64,
}

impl Index {
    fn new(meta: IndexMeta, stars: Vec<IndexStar>, table: PatternTable) -> Self {
        let mut idx = Self {
            meta,
            stars,
            table,
            content_hash: 0,
            tree: std::sync::OnceLock::new(),
        };
        idx.content_hash = file::content_hash(&idx);
        idx
    }

    pub fn meta(&self) -> &IndexMeta {
        &self.meta
    }

    pub fn stars(&self) -> &[IndexStar] {
        &self.stars
    }

    pub fn star_tree(&self) -> &kdtree::KdTree<3> {
        self.tree
            .get_or_init(|| kdtree::KdTree::build(self.stars.iter().map(|s| s.vec).collect()))
    }

    pub fn table(&self) -> &PatternTable {
        &self.table
    }

    /// FNV-1a hash of the serialized index. A loaded index retains its original file hash.
    pub fn content_hash(&self) -> u64 {
        self.content_hash
    }

    pub fn lookup(
        &self,
        shape: &Shape,
        tol: f64,
        largest_range: Option<(f64, f64)>,
    ) -> Vec<PatternMatch> {
        self.lookup_with_cancel(shape, tol, largest_range, || false)
            .expect("lookup is never cancelled")
    }

    pub fn lookup_with_cancel(
        &self,
        shape: &Shape,
        tol: f64,
        largest_range: Option<(f64, f64)>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Option<Vec<PatternMatch>> {
        if cancelled() {
            return None;
        }
        if self.table.is_empty() {
            return Some(Vec::new());
        }
        let window = KeyWindow::new(&shape.ratios, tol, self.meta.bins);
        let scan = window.len() > self.table.slots().len() as u64;
        let mut work = 0usize;
        let mut stop = || {
            work = work.wrapping_add(1);
            work % 64 == 0 && cancelled()
        };
        let mut out = Vec::new();
        let mut seen = std::collections::HashSet::new();
        let mut consider = |slot: &Slot| {
            let largest = f64::from(slot.largest);
            if largest_range
                .is_some_and(|(lo, hi)| largest < lo * (1.0 - 1e-6) || largest > hi * (1.0 + 1e-6))
                || seen.contains(&slot.stars)
            {
                return;
            }
            let pts = slot.stars.map(|i| self.stars[i as usize].vec);
            let cat = edge_shape(&pts, angle);
            let err = (0..5)
                .map(|i| (cat.ratios[i] - shape.ratios[i]).abs())
                .fold(0.0, f64::max);
            if err <= tol && (!scan || window.contains(&cat.ratios)) {
                seen.insert(slot.stars);
                out.push(PatternMatch {
                    stars: slot.stars,
                    ratio_error: err,
                    largest: cat.largest,
                });
            }
        };
        if scan {
            for slot in self.table.slots() {
                if stop() {
                    return None;
                }
                if !slot.is_empty() {
                    consider(slot);
                }
            }
        } else {
            for key in window.keys() {
                if stop() {
                    return None;
                }
                let check = key_check(key);
                for slot in self.table.probe_chain(key) {
                    if stop() {
                        return None;
                    }
                    if slot.check == check {
                        consider(slot);
                    }
                }
            }
        }
        if cancelled() {
            return None;
        }
        out.sort_unstable_by(|a, b| {
            a.ratio_error
                .total_cmp(&b.ratio_error)
                .then(a.stars.cmp(&b.stars))
        });
        (!cancelled()).then_some(out)
    }
}

pub fn image_shape(points: &[[f64; 2]; 4]) -> (Shape, [usize; 4]) {
    shape_of(points, planar)
}

#[cfg(test)]
mod tests {
    use super::*;
    use fayde_core::SkyPos;

    fn sample_index(bins: u32) -> (Index, Shape) {
        let points = [(2.0, 3.0), (2.24, 3.08), (2.09, 3.29), (2.30, 3.25)]
            .map(|(ra, dec)| SkyPos::from_degrees(ra, dec).unwrap().to_unit_vector().0);
        let (shape, order) = shape_of(&points, angle);
        let indexed_points = order.map(|i| points[i]);
        let stars = indexed_points
            .into_iter()
            .chain(indexed_points)
            .enumerate()
            .map(|(i, vec)| IndexStar {
                vec,
                mag: 8.0,
                id: i as u64,
                pattern: true,
            })
            .collect::<Vec<_>>();
        let mut table = PatternTable::for_patterns(512);
        let key = pattern::key(&shape.ratios, bins);
        table.insert(key, [4, 5, 6, 7], shape.largest);
        table.insert(key, [0, 1, 2, 3], shape.largest);
        let meta = IndexMeta {
            format_version: FORMAT_VERSION,
            catalog_name: "synthetic".into(),
            catalog_release: String::new(),
            params: BuildParams::default(),
            healpix_nside: 1,
            lattice_fields: 1,
            bins,
            star_count: stars.len(),
            pattern_count: table.len(),
            table_size: table.slots().len(),
            builder: BUILDER.into(),
        };
        (Index::new(meta, stars, table), shape)
    }

    #[test]
    fn probing_and_scanning_return_the_same_ranked_matches() {
        let (coarse, shape) = sample_index(50);
        let (fine, _) = sample_index(1_000);
        let tol = 0.02;
        assert!(KeyWindow::new(&shape.ratios, tol, 50).len() < coarse.meta().table_size as u64);
        assert!(KeyWindow::new(&shape.ratios, tol, 1_000).len() > 100_000_000);
        let matches = coarse.lookup(&shape, tol, None);
        assert_eq!(matches.len(), 2);
        assert_eq!(matches[0].stars, [0, 1, 2, 3]);
        assert_eq!(fine.lookup(&shape, tol, None), matches);
        assert!(fine
            .lookup(&shape, tol, Some((0.0, shape.largest / 2.0)))
            .is_empty());
    }

    #[test]
    fn no_match_lookups_check_cancellation_in_both_paths() {
        for bins in [50, 1_000] {
            let (index, shape) = sample_index(bins);
            let range = Some((0.0, shape.largest / 2.0));
            let mut checks = 0;
            let result = index.lookup_with_cancel(&shape, 0.02, range, || {
                checks += 1;
                checks == 2
            });
            assert!(result.is_none());
            assert!(index.lookup(&shape, 0.02, range).is_empty());
        }
    }
}
