//! Building an index from catalogue stars.
//!
//! Sources: Lang et al. 2010 (AJ 139, 1782) for uniform selection; Górski et al. 2005 for
//! HEALPix; González 2010 (Math. Geosci. 42, 49) for the Fibonacci lattice; cedar-solve for
//! lattice fields and breadth-first patterns; Brown & Stubis 2017 (TETRA) / ESA tetra3 for
//! brightest-first minimum-separation pattern thinning. The proper-motion cut is Fayde's own. See
//! `docs/references.md`.
//!
//! 1. **Selection** (astrometry.net `uniformize-catalog` idea): stars are binned into HEALPix
//!    cells sized so that a field of view holds about `stars_per_fov` of them, and the
//!    `per_cell` brightest separated stars per cell are kept (including cross-cell checks).
//!    Stars too faint, without a magnitude, or moving fast
//!    enough to distort patterns within ±`motion_years` are dropped.
//!    The pattern tier is greedily thinned brightest-first with tetra3's minimum separation.
//! 2. **Fields** (cedar-solve "lattice fields"): field centres on a Fibonacci sphere, oversampled
//!    so fields overlap. Each field takes its `stars_per_field` brightest selected stars.
//! 3. **Patterns**: 4-combinations of those stars in breadth-first order (all patterns of the 4
//!    brightest, then those adding the 5th, …), at most `patterns_per_field` per field, skipping
//!    patterns with an edge shorter than `min_edge_deg` and duplicates across fields.

use fayde_core::Portable;
use std::collections::{HashMap, HashSet};

use fayde_catalog::{CatalogId, Star};
use fayde_core::{Epoch, SkyPos};
use serde::{Deserialize, Serialize};

use crate::{
    healpix::Healpix,
    kdtree::{chord, KdTree},
    pattern::{angle, bins_for, key, shape_of, PatternTable},
    Index, IndexMeta, IndexStar,
};

/// Build parameters. Defaults suit a ~1° field with 30–60 detected sources.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BuildParams {
    /// Field-of-view diameter the patterns are built for, degrees. Patterns never exceed it.
    pub fov_deg: f64,
    /// Shortest allowed pattern edge, degrees.
    pub min_edge_deg: f64,
    /// Target selected-star density: stars per field of view.
    pub stars_per_fov: f64,
    /// Stars kept per HEALPix cell.
    pub per_cell: usize,
    /// Minimum angular separation of dense-tier stars, arcseconds. Zero disables thinning.
    /// A brighter incoming star replaces nearby kept stars, including across cell boundaries.
    #[serde(default = "default_min_sep_arcsec")]
    pub min_sep_arcsec: f64,
    /// Stars per lattice field used for patterns (taken from the pattern tier, brightest first).
    pub stars_per_field: usize,
    /// Target density of the sparse pattern tier, stars per fov² area. Greedy brightest-first
    /// thinning uses tetra3's `r_min = 0.6 * fov / sqrt(target)`; achieved density depends on
    /// available dense stars. A frame with diagonal fov covers about half of fov², so a target
    /// of 20 gives roughly 10 pattern stars per frame. Only this tier supplies patterns.
    #[serde(default = "default_pattern_stars_per_fov")]
    pub pattern_stars_per_fov: f64,
    /// Patterns per lattice field.
    pub patterns_per_field: usize,
    /// Lattice oversampling: how many fields cover each point of the sky on average. cedar-solve
    /// uses 100; at 2× only a few image 4-sets fall inside one field's brightest stars.
    pub field_oversampling: f64,
    /// Ratio tolerance the table is binned for.
    pub pattern_max_error: f64,
    /// Faintest magnitude considered.
    pub mag_limit: Option<f32>,
    /// Epoch of the index star positions (the query epoch given to the catalogue).
    pub epoch_year: f64,
    /// Stars moving more than `max_motion_frac × pattern_max_error × fov/2` within
    /// ±`motion_years` of the epoch are left out of patterns. A shift δ of one star changes a
    /// pattern's ratios by about δ / (largest edge), and fov/2 is a typical largest edge, so this
    /// keeps motion-induced ratio error below `max_motion_frac` of the bin tolerance.
    pub max_motion_frac: f64,
    /// See `max_motion_frac`.
    pub motion_years: f64,
    /// Restrict to a region: centre RA/Dec (deg) and radius (deg). `None` = whole sky.
    pub region: Option<[f64; 3]>,
}

impl Default for BuildParams {
    fn default() -> Self {
        Self {
            fov_deg: 1.0,
            min_edge_deg: 0.05,
            stars_per_fov: 20.0,
            per_cell: 2,
            min_sep_arcsec: default_min_sep_arcsec(),
            stars_per_field: 12,
            pattern_stars_per_fov: 20.0,
            patterns_per_field: 50,
            field_oversampling: 100.0,
            pattern_max_error: 0.005,
            mag_limit: None,
            epoch_year: 2025.0,
            max_motion_frac: 0.5,
            motion_years: 25.0,
            region: None,
        }
    }
}

/// Errors building an index.
#[derive(Debug, Clone, PartialEq)]
pub enum BuildError {
    /// Parameters are inconsistent.
    InvalidParams(String),
    /// No stars survived selection.
    NoStars,
    /// Reading the catalogue failed.
    Catalog(String),
}

impl std::fmt::Display for BuildError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidParams(s) => write!(f, "invalid index parameters: {s}"),
            Self::NoStars => write!(f, "no stars survived selection"),
            Self::Catalog(e) => write!(f, "catalogue: {e}"),
        }
    }
}

impl std::error::Error for BuildError {}

/// Most lattice fields a build may ask for. The field centres are held in memory at 24 bytes
/// each, so this is a 24 GB ceiling; the recommended all-sky indexes use about 5 million.
pub const MAX_LATTICE_FIELDS: f64 = 1e9;

impl BuildParams {
    /// Lattice fields the parameters imply: `field_oversampling` fields over every sky point.
    fn lattice_fields(&self) -> f64 {
        let field_r = self.fov_deg.to_radians() / 2.0;
        let cap_area = 2.0 * std::f64::consts::PI * (1.0 - field_r.pcos());
        (self.field_oversampling * 4.0 * std::f64::consts::PI / cap_area)
            .ceil()
            .max(1.0)
    }

    /// Check the parameters. [`build`] does this itself; call it first when reading the
    /// catalogue is expensive, so a bad parameter is reported before the scan.
    pub fn validate(&self) -> Result<(), BuildError> {
        let bad = |s: &str| Err(BuildError::InvalidParams(s.into()));
        if !(self.fov_deg > 0.0 && self.fov_deg < 90.0) {
            return bad("fov_deg must be in (0, 90)");
        }
        if !(self.min_edge_deg >= 0.0 && self.min_edge_deg < self.fov_deg) {
            return bad("min_edge_deg must be in [0, fov_deg)");
        }
        if !(self.min_sep_arcsec.is_finite() && (0.0..=648_000.0).contains(&self.min_sep_arcsec)) {
            return bad("min_sep_arcsec must be finite and in [0, 648000]");
        }
        if !(self.pattern_stars_per_fov.is_finite() && self.pattern_stars_per_fov > 0.0) {
            return bad("pattern_stars_per_fov must be finite and positive");
        }
        if self.stars_per_field < 4 || self.per_cell == 0 || self.patterns_per_field == 0 {
            return bad("need stars_per_field >= 4, per_cell >= 1, patterns_per_field >= 1");
        }
        if !(self.pattern_max_error > 0.0 && self.pattern_max_error < 0.25) {
            return bad("pattern_max_error must be in (0, 0.25)");
        }
        if !(self.stars_per_fov > 0.0 && self.field_oversampling > 0.0) {
            return bad("stars_per_fov and field_oversampling must be positive");
        }
        if !(self.stars_per_fov.is_finite() && self.field_oversampling.is_finite()) {
            return bad("stars_per_fov and field_oversampling must be finite");
        }
        if self.lattice_fields() > MAX_LATTICE_FIELDS {
            return Err(BuildError::InvalidParams(format!(
                "fov_deg {} with field_oversampling {} needs {:.2e} lattice fields (limit {:.0e}): \
                 use a wider field of view or less oversampling",
                self.fov_deg,
                self.field_oversampling,
                self.lattice_fields(),
                MAX_LATTICE_FIELDS
            )));
        }
        if crate::pattern::bins_for(self.pattern_max_error) > crate::pattern::MAX_BINS {
            return bad("pattern_max_error is too small for the pattern key");
        }
        if self.mag_limit.is_some_and(|m| !m.is_finite()) {
            return bad("mag_limit must be finite");
        }
        // NaN and negative values would switch the proper-motion cut off without saying so.
        if !(self.max_motion_frac >= 0.0 && self.motion_years >= 0.0) {
            return bad("max_motion_frac and motion_years must be >= 0");
        }
        match Epoch::julian_year(self.epoch_year) {
            None => return bad("epoch_year must be finite"),
            Some(e) if !e.is_plausible() => return bad("epoch_year must be within 1800 to 2200"),
            Some(_) => {}
        }
        Ok(())
    }

    /// Epoch to request catalogue positions at.
    pub fn epoch(&self) -> Epoch {
        Epoch::julian_year(self.epoch_year).expect("validated")
    }
}

/// Default for [`BuildParams::pattern_stars_per_fov`] (also used for files written before it
/// existed).
pub fn default_pattern_stars_per_fov() -> f64 {
    20.0
}

/// Default dense-tier minimum separation, also used when reading older index metadata.
pub fn default_min_sep_arcsec() -> f64 {
    10.0
}

/// Select stars uniformly (step 1). Output is deterministic: by HEALPix cell, then brightness.
pub fn select_stars(stars: &[Star], p: &BuildParams) -> (Vec<Star>, Healpix) {
    let mut sel = Selector::new(p);
    for s in stars {
        sel.push(*s);
    }
    sel.finish()
}

/// A star with its sort key (magnitude, id).
type Ranked = (f32, u64, Star);

/// Streaming version of [`select_stars`]: keeps at most `per_cell` separated stars per cell,
/// so an all-sky catalogue can be selected without holding it in memory. Replacements use
/// magnitude, then ID. Given the same stream, selection is deterministic; rejected stars are
/// not reconsidered if a later replacement would free space for them.
pub struct Selector {
    hp: Healpix,
    per_cell: usize,
    mag_limit: Option<f32>,
    max_motion_mas: f64,
    motion_years: f64,
    cells: HashMap<u64, Vec<Ranked>>,
    separation: Option<SeparationGrid>,
    seen: u64,
}

impl Selector {
    /// A selector for `p`.
    pub fn new(p: &BuildParams) -> Self {
        let fov = p.fov_deg.to_radians();
        // A fov-sized square holds stars_per_fov / per_cell cells.
        let cells_per_fov = (p.stars_per_fov / p.per_cell as f64).max(1.0);
        Self {
            hp: Healpix::for_pixel_side(fov / cells_per_fov.sqrt()),
            per_cell: p.per_cell,
            mag_limit: p.mag_limit,
            max_motion_mas: p.max_motion_frac * p.pattern_max_error * (p.fov_deg / 2.0) * 3.6e6,
            motion_years: p.motion_years,
            cells: HashMap::new(),
            separation: (p.min_sep_arcsec > 0.0)
                .then(|| SeparationGrid::new((p.min_sep_arcsec / 3600.0).to_radians())),
            seen: 0,
        }
    }

    /// Offer one star.
    pub fn push(&mut self, s: Star) {
        self.seen += 1;
        let Some(m) = s.mag.filter(|m| m.is_finite()) else {
            return;
        };
        if self.mag_limit.is_some_and(|lim| m > lim) {
            return;
        }
        let pm = s.pm_mas_yr[0]
            .unwrap_or(0.0)
            .phypot(s.pm_mas_yr[1].unwrap_or(0.0));
        if f64::from(pm) * self.motion_years > self.max_motion_mas {
            return;
        }
        let v = s.pos.to_unit_vector();
        let pixel = self.hp.pixel(v);
        let rank = |c: &Ranked| c.0.total_cmp(&m).then(c.1.cmp(&s.id));
        // Most catalogue stars cannot enter a full cell. Reject them before spatial queries.
        if self.cells.get(&pixel).is_some_and(|cell| {
            cell.len() == self.per_cell && rank(cell.last().expect("nonempty cell")).is_lt()
        }) {
            return;
        }
        if let Some(grid) = &mut self.separation {
            let close: Vec<_> = grid.nearby(v.0).copied().collect();
            // Check every neighbour before changing anything: one better neighbour vetoes
            // replacement even when the candidate is brighter than other neighbours.
            if close
                .iter()
                .any(|c| !c.mag.total_cmp(&m).then(c.id.cmp(&s.id)).is_gt())
            {
                return;
            }
            for c in close {
                let cell = self
                    .cells
                    .get_mut(&self.hp.pixel(fayde_core::UnitVector(c.vec)))
                    .expect("grid contains only kept stars");
                cell.retain(|entry| entry.1 != c.id);
                grid.remove(c);
            }
        }
        let cell = self.cells.entry(pixel).or_default();
        let pos = cell.partition_point(|c| rank(c).is_lt());
        cell.insert(pos, (m, s.id, s));
        if let Some(grid) = &mut self.separation {
            if cell.len() > self.per_cell {
                let dropped = cell.pop().expect("overfull cell").2;
                grid.remove(index_star(&dropped));
            }
            grid.insert(index_star(&s));
        } else {
            cell.truncate(self.per_cell);
        }
    }

    /// Stars offered so far.
    pub fn seen(&self) -> u64 {
        self.seen
    }

    /// Selected stars, ordered by cell then brightness (deterministic).
    pub fn finish(self) -> (Vec<Star>, Healpix) {
        let mut cells: Vec<(u64, Vec<Ranked>)> = self.cells.into_iter().collect();
        cells.sort_by_key(|c| c.0);
        let out = cells
            .into_iter()
            .flat_map(|(_, v)| v.into_iter().map(|c| c.2))
            .collect();
        (out, self.hp)
    }
}

/// Builder identity stored in every index (and so in its content hash). It names the release
/// that introduced the current index-building code, not the running crate version: a release
/// that leaves the builder unchanged builds the same bytes, so published hashes stay valid.
/// Change it, to the version of the release, whenever a change alters the bytes an index build
/// produces (selection, patterns, table, serialization); the golden-bytes tests then fail until
/// re-pinned, and the release is breaking under the versioning rules (README).
pub const BUILDER: &str = "fayde-index 0.1.0";

/// Field centres on a Fibonacci sphere (deterministic).
pub fn lattice(n: usize) -> Vec<[f64; 3]> {
    // Golden angle as a fraction of a turn, (3 − √5)/2. The angle is reduced to one turn before
    // sin/cos: `golden_angle · i` reaches millions of radians, where argument reduction is slow
    // and loses the low bits.
    let golden_turn = (3.0 - 5f64.sqrt()) / 2.0;
    (0..n)
        .map(|i| {
            let z = 1.0 - 2.0 * (i as f64 + 0.5) / n as f64;
            let s = (1.0 - z * z).sqrt();
            let phi = std::f64::consts::TAU * (golden_turn * i as f64).fract();
            [s * phi.pcos(), s * phi.psin(), z]
        })
        .collect()
}

/// Build an index from catalogue stars already propagated to `params.epoch_year`.
pub fn build(
    stars: &[Star],
    catalog: CatalogId,
    params: &BuildParams,
) -> Result<Index, BuildError> {
    params.validate()?;
    let (selected, hp) = select_stars(stars, params);
    build_selected(selected, hp, catalog, params, &mut |_| {})
}

/// Build an index by streaming every star of `catalog` (at `params.epoch_year`) through the
/// selector. Memory scales with the selected stars, not the catalogue.
pub fn build_from_catalog(
    catalog: &dyn fayde_catalog::Catalog,
    params: &BuildParams,
) -> Result<Index, BuildError> {
    build_from_catalog_with_progress(catalog, params, &mut |_| {})
}

/// [`build_from_catalog`] reporting progress messages (stage changes and every ~1% of fields).
pub fn build_from_catalog_with_progress(
    catalog: &dyn fayde_catalog::Catalog,
    params: &BuildParams,
    progress: &mut dyn FnMut(&str),
) -> Result<Index, BuildError> {
    params.validate()?;
    progress("scanning catalogue");
    let query = fayde_catalog::StarQuery {
        epoch: Some(params.epoch()),
        mag_max: params.mag_limit,
        ..fayde_catalog::StarQuery::default()
    };
    let mut sel = Selector::new(params);
    match params.region {
        // A regional index reads only its cone.
        Some([ra, dec, r]) => {
            let c = SkyPos::from_degrees(ra, dec)
                .ok_or_else(|| BuildError::InvalidParams("region centre".into()))?;
            for s in catalog
                .cone(c, r.to_radians(), &query)
                .map_err(|e| BuildError::Catalog(e.to_string()))?
            {
                sel.push(s);
            }
        }
        None => catalog
            .scan(&query, &mut |s| sel.push(s))
            .map_err(|e| BuildError::Catalog(e.to_string()))?,
    }
    progress(&format!("scanned {} stars", sel.seen()));
    let (selected, hp) = sel.finish();
    progress(&format!(
        "selected {} stars (HEALPix nside {})",
        selected.len(),
        hp.nside()
    ));
    build_selected(selected, hp, catalog.id(), params, progress)
}

fn index_star(s: &Star) -> IndexStar {
    IndexStar {
        vec: s.pos.to_unit_vector().0,
        mag: s.mag.expect("selected stars have magnitudes"),
        id: s.id,
        pattern: false,
    }
}

// Incremental radius search in Cartesian unit-vector space. Buckets have side at least the
// chord radius, so any neighbour must be in one of the 27 touching buckets. Exact distances
// eliminate bucket-boundary artifacts, including at RA wrap and the poles. Only kept stars
// are stored; insertion/removal do not rebuild the static KdTree used later for field queries.
struct SeparationGrid {
    radius: f64,
    buckets: HashMap<[i64; 3], Vec<IndexStar>>,
}

impl SeparationGrid {
    fn new(angle: f64) -> Self {
        Self {
            radius: chord(angle.min(std::f64::consts::PI)),
            buckets: HashMap::new(),
        }
    }

    fn bucket(&self, v: [f64; 3]) -> [i64; 3] {
        // Keep bucket coordinates representable even for sub-nanoradian requested radii.
        v.map(|x| (x / self.radius.max(1e-12)).floor() as i64)
    }

    fn nearby(&self, v: [f64; 3]) -> impl Iterator<Item = &IndexStar> {
        let [x, y, z] = self.bucket(v);
        (-1..=1)
            .flat_map(move |dx| {
                (-1..=1).flat_map(move |dy| (-1..=1).map(move |dz| [x + dx, y + dy, z + dz]))
            })
            .filter_map(|key| self.buckets.get(&key))
            .flatten()
            .filter(move |s| {
                (0..3).map(|k| (v[k] - s.vec[k]).ppowi(2)).sum::<f64>() <= self.radius * self.radius
            })
    }

    fn insert(&mut self, s: IndexStar) {
        self.buckets.entry(self.bucket(s.vec)).or_default().push(s);
    }

    fn remove(&mut self, s: IndexStar) {
        let key = self.bucket(s.vec);
        let bucket = self.buckets.get_mut(&key).expect("kept star has a bucket");
        let i = bucket
            .iter()
            .position(|c| c.id == s.id)
            .expect("kept star is in bucket");
        bucket.swap_remove(i);
        if bucket.is_empty() {
            self.buckets.remove(&key);
        }
    }
}

fn pattern_min_sep(params: &BuildParams) -> f64 {
    // ESA tetra3 generate_database (Brown & Stubis 2017 / TETRA): the characteristic spacing
    // is fov/sqrt(N), scaled by 0.6 for greedy thinning. This is an empirical density relation,
    // not an exact count: sparse inputs undershoot and saturated inputs can overshoot N/fov².
    0.6 * params.fov_deg.to_radians() / params.pattern_stars_per_fov.sqrt()
}

// Kept separate from pattern enumeration so thinning can be measured on an all-sky input.
fn select_pattern_stars(index_stars: &mut [IndexStar], params: &BuildParams) -> Vec<usize> {
    let mut order: Vec<_> = (0..index_stars.len()).collect();
    order.sort_unstable_by(|&a, &b| {
        index_stars[a]
            .mag
            .total_cmp(&index_stars[b].mag)
            .then(index_stars[a].id.cmp(&index_stars[b].id))
    });
    let mut grid = SeparationGrid::new(pattern_min_sep(params));
    let mut selected = Vec::new();
    for i in order {
        let star = &mut index_stars[i];
        star.pattern = grid.nearby(star.vec).next().is_none();
        if star.pattern {
            grid.insert(*star);
            selected.push(i);
        }
    }
    // Preserve the index's stable order for field searches and pattern IDs.
    selected.sort_unstable();
    selected
}

/// Patterns in generation order, as their canonical star numbers only (16 bytes each). The hash
/// key and the largest edge are recomputed from the stars when the table is filled, and the
/// fixed-size chunks never reallocate: a single growing vector of (key, stars, edge) held twice
/// the bytes per pattern and doubled its allocation on the way up.
#[derive(Default)]
struct PatternLog {
    chunks: Vec<Vec<[u32; 4]>>,
    len: usize,
}

impl PatternLog {
    /// Patterns per chunk (16 MiB).
    const CHUNK: usize = 1 << 20;

    fn push(&mut self, stars: [u32; 4]) {
        match self.chunks.last_mut() {
            Some(chunk) if chunk.len() < Self::CHUNK => chunk.push(stars),
            _ => self.chunks.push(vec![stars]),
        }
        self.len += 1;
    }

    fn len(&self) -> usize {
        self.len
    }
}

fn build_selected(
    selected: Vec<Star>,
    hp: Healpix,
    catalog: CatalogId,
    params: &BuildParams,
    progress: &mut dyn FnMut(&str),
) -> Result<Index, BuildError> {
    let region = params
        .region
        .map(|[ra, dec, r]| {
            SkyPos::from_degrees(ra, dec)
                .map(|c| (c.to_unit_vector().0, r.to_radians()))
                .ok_or_else(|| BuildError::InvalidParams("region centre".into()))
        })
        .transpose()?;
    let in_region = |v: &[f64; 3], pad: f64| match region {
        None => true,
        Some((c, r)) => angle(v, &c) <= r + pad,
    };
    let selected: Vec<Star> = selected
        .into_iter()
        .filter(|s| in_region(&s.pos.to_unit_vector().0, 0.0))
        .collect();
    if selected.is_empty() {
        return Err(BuildError::NoStars);
    }
    let mut index_stars: Vec<IndexStar> = selected.iter().map(index_star).collect();
    let fov = params.fov_deg.to_radians();
    let pattern_stars = select_pattern_stars(&mut index_stars, params);
    let area = region.map_or(4.0 * std::f64::consts::PI, |(_, r)| {
        4.0 * std::f64::consts::PI * (r.min(std::f64::consts::PI) / 2.0).psin().ppowi(2)
    });
    let density = pattern_stars.len() as f64 * fov * fov / area;
    progress(&format!(
        "pattern tier: {} of {} index stars (min separation {:.3} deg, achieved {:.2} per fov² over {}, target {:.2})",
        pattern_stars.len(), index_stars.len(), pattern_min_sep(params).to_degrees(), density,
        if region.is_some() { "region" } else { "all sky" }, params.pattern_stars_per_fov
    ));
    let tree = KdTree::build(pattern_stars.iter().map(|&i| index_stars[i].vec).collect());

    let field_r = fov / 2.0;
    let n_fields = params.lattice_fields() as usize;
    let min_edge = params.min_edge_deg.to_radians();

    let mut seen: HashSet<[u32; 4]> = HashSet::new();
    let mut patterns = PatternLog::default();
    let bins = bins_for(params.pattern_max_error);
    let step = (n_fields / 100).max(1);
    for (fi, centre) in lattice(n_fields).into_iter().enumerate() {
        if fi % step == 0 {
            progress(&format!(
                "fields {fi}/{n_fields}, {} patterns",
                patterns.len()
            ));
        }
        if !in_region(&centre, field_r) {
            continue;
        }
        let mut members: Vec<usize> = tree
            .within(centre, chord(field_r))
            .into_iter()
            .map(|k| pattern_stars[k])
            .collect();
        members.sort_by(|&a, &b| {
            index_stars[a]
                .mag
                .total_cmp(&index_stars[b].mag)
                .then(a.cmp(&b))
        });
        members.truncate(params.stars_per_field);
        let mut made = 0;
        'field: for n in 4..=members.len() {
            // Breadth-first: every combination whose faintest member is members[n-1].
            let last = members[n - 1];
            for a in 0..n - 1 {
                for b in a + 1..n - 1 {
                    for c in b + 1..n - 1 {
                        let ids = [members[a], members[b], members[c], last];
                        let mut sorted = ids.map(|i| i as u32);
                        sorted.sort_unstable();
                        if seen.contains(&sorted) {
                            continue;
                        }
                        let pts = ids.map(|i| index_stars[i].vec);
                        let (shape, order) = shape_of(&pts, angle);
                        if shape.ratios[0] * shape.largest < min_edge || shape.largest > fov {
                            continue;
                        }
                        seen.insert(sorted);
                        patterns.push(order.map(|k| ids[k] as u32));
                        made += 1;
                        if made == params.patterns_per_field {
                            break 'field;
                        }
                    }
                }
            }
        }
    }

    // The duplicate check is over: free it before the table, the largest allocation, exists.
    drop(seen);
    let mut table = PatternTable::for_patterns(patterns.len());
    for chunk in patterns.chunks {
        for stars in chunk {
            // The six edges are the same numbers whatever the order of the four stars, so
            // the key and the largest edge are those found when the pattern was generated.
            let pts = stars.map(|i| index_stars[i as usize].vec);
            let (shape, _) = shape_of(&pts, angle);
            table.insert(key(&shape.ratios, bins), stars, shape.largest);
        }
        // `chunk` is freed here, so the log shrinks as the table fills.
    }
    let meta = IndexMeta {
        format_version: crate::FORMAT_VERSION,
        catalog_name: catalog.name.to_string(),
        catalog_release: catalog.release.to_string(),
        params: params.clone(),
        healpix_nside: hp.nside(),
        lattice_fields: n_fields,
        bins,
        star_count: index_stars.len(),
        pattern_count: table.len(),
        table_size: table.slots().len(),
        builder: BUILDER.to_string(),
    };
    Ok(Index::new(meta, index_stars, table))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn star(id: u64, ra: f64, dec: f64, mag: f32) -> Star {
        Star {
            id,
            pos: SkyPos::from_degrees(ra, dec).unwrap(),
            mag: Some(mag),
            pos_err_mas: [10.0, 10.0],
            pm_mas_yr: [Some(1.0), Some(1.0)],
        }
    }

    #[test]
    fn selection_caps_per_cell_and_drops_fast_and_faint() {
        let mut stars: Vec<Star> = (0..50)
            .map(|i| star(i, 10.0 + 0.0001 * i as f64, 5.0, 10.0 + i as f32 * 0.1))
            .collect();
        let mut fast = star(999, 10.0, 5.0, 1.0);
        fast.pm_mas_yr = [Some(5_000.0), Some(0.0)];
        stars.push(fast);
        let p = BuildParams {
            mag_limit: Some(12.0),
            min_sep_arcsec: 0.0, // Isolate the per-cell cap from separation thinning.
            ..BuildParams::default()
        };
        let (sel, _) = select_stars(&stars, &p);
        // All 50 are in one cell: only the 2 brightest survive; the fast star is dropped.
        assert_eq!(sel.iter().map(|s| s.id).collect::<Vec<_>>(), vec![0, 1]);
    }

    #[test]
    fn pattern_tier_has_no_coarse_cell_boundary_artifact() {
        // At this scale v3 uses nside 64 and keeps one star per coarse cell.
        // The two bright stars are 0.1*fov apart, both west of the RA=45° boundary.
        let p = BuildParams {
            fov_deg: 7.0,
            pattern_stars_per_fov: 64.0,
            ..BuildParams::default()
        };
        let stars = [
            star(1, 44.25, 0.0, 5.0),
            star(2, 44.95, 0.0, 6.0),
            star(3, 45.05, 0.0, 15.0),
        ];
        let hp = Healpix::new(64).unwrap();
        assert_eq!(
            hp.pixel(stars[0].pos.to_unit_vector()),
            hp.pixel(stars[1].pos.to_unit_vector())
        );
        assert_ne!(
            hp.pixel(stars[1].pos.to_unit_vector()),
            hp.pixel(stars[2].pos.to_unit_vector())
        );
        let mut index_stars: Vec<_> = stars
            .iter()
            .map(|s| IndexStar {
                vec: s.pos.to_unit_vector().0,
                mag: s.mag.unwrap(),
                id: s.id,
                pattern: false,
            })
            .collect();
        let selected = select_pattern_stars(&mut index_stars, &p);
        assert_eq!(
            selected
                .iter()
                .map(|&i| index_stars[i].id)
                .collect::<Vec<_>>(),
            [1, 2]
        );
        assert_eq!(
            index_stars.iter().map(|s| s.pattern).collect::<Vec<_>>(),
            [true, true, false]
        );
    }

    #[test]
    fn dense_tier_rejects_three_arcsec_double() {
        let bright = star(1, 10.0, 5.0, 5.0);
        let faint = star(2, 10.0 + 3.0 / 3600.0, 5.0, 6.0);
        let p = BuildParams::default();
        for stars in [[bright, faint], [faint, bright]] {
            let (selected, hp) = select_stars(&stars, &p);
            assert_eq!(
                hp.pixel(bright.pos.to_unit_vector()),
                hp.pixel(faint.pos.to_unit_vector())
            );
            assert_eq!(selected.iter().map(|s| s.id).collect::<Vec<_>>(), [1]);
        }
    }

    #[test]
    fn pattern_tier_uniform_sky_density_and_separation() {
        // Uniform all-sky positions and independent deterministic brightness ranks.
        // 1.6M / 41253 square degrees = 38.8 dense stars per 1° FOV squared.
        let p = BuildParams::default();
        let mut seed = 7u64;
        let mut rnd = || {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (seed >> 11) as f64 / (1u64 << 53) as f64
        };
        let mut stars: Vec<_> = (0..1_600_000)
            .map(|id| {
                let z = 2.0 * rnd() - 1.0;
                let phi = std::f64::consts::TAU * rnd();
                let r = (1.0 - z * z).sqrt();
                IndexStar {
                    vec: [r * phi.pcos(), r * phi.psin(), z],
                    mag: (20.0 * rnd()) as f32,
                    id,
                    pattern: false,
                }
            })
            .collect();
        let start = std::time::Instant::now();
        let selected = select_pattern_stars(&mut stars, &p);
        let elapsed = start.elapsed();
        let density =
            selected.len() as f64 * p.fov_deg.to_radians().ppowi(2) / (4.0 * std::f64::consts::PI);
        eprintln!(
            "pattern thinning: {} / {} stars, {density:.3} per fov² (target {}), {elapsed:.3?}",
            selected.len(),
            stars.len(),
            p.pattern_stars_per_fov
        );
        assert!((density / p.pattern_stars_per_fov - 1.0).abs() <= 0.3);
        // Independent static-tree oracle: every accepted point has only itself in its radius.
        let tree = KdTree::build(selected.iter().map(|&i| stars[i].vec).collect());
        let radius = chord(0.6 * p.fov_deg.to_radians() / p.pattern_stars_per_fov.sqrt());
        for i in 0..tree.len() {
            assert_eq!(
                tree.within(tree.point(i), radius),
                [i],
                "pattern neighbours at {i}"
            );
        }
    }

    #[test]
    fn dense_tier_checks_adjacent_cells_wrap_and_poles() {
        let p = BuildParams::default();
        for (a, b) in [
            ((45.0 - 1.5 / 3600.0, 0.0), (45.0 + 1.5 / 3600.0, 0.0)),
            ((359.9996, 0.0), (0.0004, 0.0)),
            ((0.0, 89.9996), (180.0, 89.9996)),
        ] {
            let bright = star(1, a.0, a.1, 5.0);
            // Equal magnitude tests the ID tie-break in both stream orders as well.
            for mag in [5.0, 6.0] {
                let faint = star(2, b.0, b.1, mag);
                for stars in [[bright, faint], [faint, bright]] {
                    let (selected, hp) = select_stars(&stars, &p);
                    assert_ne!(
                        hp.pixel(bright.pos.to_unit_vector()),
                        hp.pixel(faint.pos.to_unit_vector())
                    );
                    assert_eq!(selected.iter().map(|s| s.id).collect::<Vec<_>>(), [1]);
                }
            }
        }
    }

    #[test]
    fn dense_replacement_checks_all_neighbours_before_removing_any() {
        let p = BuildParams::default();
        let left = star(1, 45.0 - 8.0 / 3600.0, 0.0, 5.0);
        let right = star(2, 45.0 + 8.0 / 3600.0, 0.0, 7.0);
        for neighbours in [[left, right], [right, left]] {
            let mut sel = Selector::new(&p);
            for s in neighbours {
                sel.push(s);
            }
            sel.push(star(3, 45.0, 0.0, 6.0));
            let (kept, _) = sel.finish();
            let mut ids: Vec<_> = kept.iter().map(|s| s.id).collect();
            ids.sort_unstable();
            assert_eq!(ids, [1, 2]);
            let mut sel = Selector::new(&p);
            for s in neighbours {
                sel.push(s);
            }
            sel.push(star(4, 45.0, 0.0, 4.0));
            let (kept, _) = sel.finish();
            assert_eq!(kept.iter().map(|s| s.id).collect::<Vec<_>>(), [4]);
        }
    }

    #[test]
    fn dense_cell_eviction_removes_spatial_entry() {
        let p = BuildParams {
            per_cell: 1,
            ..BuildParams::default()
        };
        let a = star(1, 44.9996, 0.0, 5.0);
        let b = star(2, 44.9, 0.0, 4.0);
        let c = star(3, 45.0004, 0.0, 6.0);
        let (selected, hp) = select_stars(&[a, b, c], &p);
        assert_eq!(
            hp.pixel(a.pos.to_unit_vector()),
            hp.pixel(b.pos.to_unit_vector())
        );
        assert_ne!(
            hp.pixel(a.pos.to_unit_vector()),
            hp.pixel(c.pos.to_unit_vector())
        );
        let mut ids: Vec<_> = selected.iter().map(|s| s.id).collect();
        ids.sort_unstable();
        assert_eq!(ids, [2, 3]);
    }

    #[test]
    fn pattern_tier_matches_greedy_oracle_and_id_ties() {
        let p = BuildParams {
            fov_deg: 30.0,
            ..BuildParams::default()
        };
        let mut stars: Vec<_> = lattice(2000)
            .into_iter()
            .enumerate()
            .map(|(i, vec)| IndexStar {
                vec,
                mag: (i % 3) as f32,
                id: (2000 - i) as u64,
                pattern: false,
            })
            .collect();
        let mut order = stars.clone();
        order.sort_by(|a, b| a.mag.total_cmp(&b.mag).then(a.id.cmp(&b.id)));
        let mut expected: Vec<IndexStar> = Vec::new();
        for s in order {
            if expected
                .iter()
                .all(|a| angle(&a.vec, &s.vec) > pattern_min_sep(&p))
            {
                expected.push(s);
            }
        }
        let mut expected_ids: Vec<_> = expected.iter().map(|s| s.id).collect();
        expected_ids.sort_unstable();
        for _ in 0..2 {
            let selected = select_pattern_stars(&mut stars, &p);
            let mut ids: Vec<_> = selected.iter().map(|&i| stars[i].id).collect();
            ids.sort_unstable();
            assert_eq!(ids, expected_ids);
            stars.reverse();
        }
    }

    #[test]
    fn synthetic_build_is_deterministic_and_reports_achieved_density() {
        let stars: Vec<_> = (0..64)
            .map(|i| {
                star(
                    i,
                    44.6 + (i % 8) as f64 * 0.1,
                    -0.35 + (i / 8) as f64 * 0.1,
                    (i % 7) as f32,
                )
            })
            .collect();
        let catalog = CatalogId {
            name: "synthetic",
            release: "0",
            location: String::new(),
        };
        for region in [None, Some([45.0, 0.0, 1.0])] {
            let p = BuildParams {
                field_oversampling: 0.05,
                region,
                ..BuildParams::default()
            };
            let (selected, hp) = select_stars(&stars, &p);
            let mut messages = Vec::new();
            let first = build_selected(selected, hp, catalog.clone(), &p, &mut |s| {
                messages.push(s.to_string())
            })
            .unwrap();
            let second = build(&stars, catalog.clone(), &p).unwrap();
            assert_eq!(first.content_hash(), second.content_hash());
            assert_eq!(first.meta().format_version, 2);
            let n = first.stars().iter().filter(|s| s.pattern).count();
            let area = if region.is_some() {
                2.0 * std::f64::consts::PI * (1.0 - 1f64.to_radians().pcos())
            } else {
                4.0 * std::f64::consts::PI
            };
            let density = n as f64 * p.fov_deg.to_radians().ppowi(2) / area;
            assert!(messages
                .iter()
                .any(|s| s.contains(&format!("achieved {density:.2} per fov²"))));
        }
    }

    #[test]
    fn old_build_params_default_dense_separation() {
        let mut value = serde_json::to_value(BuildParams::default()).unwrap();
        value.as_object_mut().unwrap().remove("min_sep_arcsec");
        let params: BuildParams = serde_json::from_value(value).unwrap();
        assert_eq!(params.min_sep_arcsec, 10.0);
        assert_eq!(params, BuildParams::default());
        for bad in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(BuildParams {
                min_sep_arcsec: bad,
                ..params.clone()
            }
            .validate()
            .is_err());
            assert!(BuildParams {
                pattern_stars_per_fov: bad,
                ..params.clone()
            }
            .validate()
            .is_err());
        }
        assert!(BuildParams {
            min_sep_arcsec: 0.0,
            ..params
        }
        .validate()
        .is_ok());
    }

    #[test]
    fn rejects_bad_params() {
        let p = BuildParams {
            fov_deg: -1.0,
            ..BuildParams::default()
        };
        let id = CatalogId {
            name: "T",
            release: "0",
            location: String::new(),
        };
        assert!(matches!(
            build(&[], id, &p),
            Err(BuildError::InvalidParams(_))
        ));
    }
}
