//! Pattern index for Fayde: uniform star selection, 4-star pattern hash table, and a versioned
//! on-disk format.
//!
//! See `docs/design/plate-solver.md` §5.2. Blind solving uses these pattern indexes in scale bands;
//! astrometry.net-style quad codes were considered and not adopted (`docs/references.md`).

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

use pattern::{angle, keys_within, planar, shape_of, PatternTable, Shape};

/// On-disk format version. Bump on any incompatible change.
pub const FORMAT_VERSION: u32 = 2;

/// Index metadata, stored as JSON in the file header.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct IndexMeta {
    /// [`FORMAT_VERSION`] at build time.
    pub format_version: u32,
    /// Source catalogue name.
    pub catalog_name: String,
    /// Source catalogue release.
    pub catalog_release: String,
    /// Build parameters (including epoch and region).
    pub params: BuildParams,
    /// HEALPix resolution used for star selection.
    pub healpix_nside: u32,
    /// Number of lattice fields over the whole sky.
    pub lattice_fields: usize,
    /// Bins per ratio in pattern keys.
    pub bins: u32,
    /// Stars in the index.
    pub star_count: usize,
    /// Patterns in the table.
    pub pattern_count: usize,
    /// Table slots.
    pub table_size: usize,
    /// Builder name and version.
    pub builder: String,
}

/// A star as stored in the index.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IndexStar {
    /// Unit vector at the index epoch.
    pub vec: [f64; 3],
    /// Magnitude in the build band.
    pub mag: f32,
    /// Catalogue identifier.
    pub id: u64,
    /// Member of the sparse pattern tier (patterns are built only from these).
    pub pattern: bool,
}

/// A loaded or freshly built index.
#[derive(Debug, Clone)]
pub struct Index {
    meta: IndexMeta,
    stars: Vec<IndexStar>,
    table: PatternTable,
    content_hash: u64,
    /// k-d tree over star vectors, built on first use and kept.
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

/// A catalogue pattern matching an image pattern.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PatternMatch {
    /// Index-star numbers in canonical order. Element `i` corresponds to the image star at
    /// canonical position `i` of the query.
    pub stars: [u32; 4],
    /// Largest |Δratio| between the catalogue and image shapes.
    pub ratio_error: f64,
    /// Catalogue largest edge, radians.
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

    /// Metadata.
    pub fn meta(&self) -> &IndexMeta {
        &self.meta
    }

    /// Stars, addressed by the numbers stored in patterns.
    pub fn stars(&self) -> &[IndexStar] {
        &self.stars
    }

    /// k-d tree over the star unit vectors (indices into [`Index::stars`]), built once.
    pub fn star_tree(&self) -> &kdtree::KdTree<3> {
        self.tree
            .get_or_init(|| kdtree::KdTree::build(self.stars.iter().map(|s| s.vec).collect()))
    }

    /// Pattern table.
    pub fn table(&self) -> &PatternTable {
        &self.table
    }

    /// 64-bit FNV-1a hash of the serialized content; goes into solution provenance.
    pub fn content_hash(&self) -> u64 {
        self.content_hash
    }

    /// Catalogue patterns whose shape matches `shape` within `tol` in every ratio, optionally
    /// restricted to a range of largest-edge sizes (radians). Results are ordered by ratio error,
    /// then by star numbers (deterministic).
    pub fn lookup(
        &self,
        shape: &Shape,
        tol: f64,
        largest_range: Option<(f64, f64)>,
    ) -> Vec<PatternMatch> {
        let mut out = Vec::new();
        let mut seen = std::collections::HashSet::new();
        for k in keys_within(&shape.ratios, tol, self.meta.bins) {
            for slot in self.table.probe(k) {
                if let Some((lo, hi)) = largest_range {
                    let l = f64::from(slot.largest);
                    if l < lo * (1.0 - 1e-6) || l > hi * (1.0 + 1e-6) {
                        continue;
                    }
                }
                if !seen.insert(slot.stars) {
                    continue;
                }
                let pts = slot.stars.map(|i| self.stars[i as usize].vec);
                let (cat, _) = shape_of(&pts, angle);
                let err = (0..5)
                    .map(|i| (cat.ratios[i] - shape.ratios[i]).abs())
                    .fold(0.0, f64::max);
                if err <= tol {
                    out.push(PatternMatch {
                        stars: slot.stars,
                        ratio_error: err,
                        largest: cat.largest,
                    });
                }
            }
        }
        out.sort_by(|a, b| {
            a.ratio_error
                .total_cmp(&b.ratio_error)
                .then(a.stars.cmp(&b.stars))
        });
        out
    }
}

/// Shape and canonical order of an image pattern from tangent-plane (or pixel) coordinates.
/// Distances in a ~1° tangent plane differ from angles by < 1e-4 relative, well inside pattern
/// tolerances.
pub fn image_shape(points: &[[f64; 2]; 4]) -> (Shape, [usize; 4]) {
    shape_of(points, planar)
}
