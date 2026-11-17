//! The catalogue stars a solve checks its hypotheses against.
//!
//! A hypothesis is judged on the stars inside the frame it implies, which is a small part of
//! the sky a pointing prior allows: a 2° prior around a 0.6° frame covers 46 frame areas. A
//! [`RefSet`] over such a prior therefore reads the catalogue a frame-sized cone at a time, as
//! hypotheses ask for stars, and keeps what it has read. What a hypothesis sees is exactly what
//! it would see with the whole prior cone in memory: the same stars, in the same order.
//!
//! Reading in parts is Fayde's own arrangement; astrometry.net reaches the same end by storing
//! its reference stars in the index, in a k-d tree it searches per match (Lang et al. 2010, §2).

use std::cell::RefCell;

use fayde_catalog::{Catalog, Star, StarQuery};
use fayde_core::{SkyPos, UnitVector};
use fayde_index::{kdtree::chord, pattern::angle};

use crate::{hash::Map, request::SolveError};

/// Radius of a cone read for a hypothesis, as a multiple of the radius asked for (plus
/// [`LOCAL_PAD`]): the margin lets the refits that follow, which move the frame by a few
/// pixels, use the same cone.
const LOCAL_GROW: f64 = 1.25;
/// See [`LOCAL_GROW`]: keeps the margin above [`HOLD_SLACK`] for the smallest radii. Radians.
const LOCAL_PAD: f64 = 1e-6;
/// A cone is read on its own only while it is smaller than this fraction of the whole set's
/// radius; a larger one is most of the set, which is then read whole.
const LOCAL_MAX_FRACTION: f64 = 0.5;
/// Cones read on their own before the rest of the set is read in one go. Hypotheses that pass
/// the pre-filter by chance are scattered over the prior: past a handful of them, their cones
/// together cost more than the one read that serves all that follow.
const MAX_LOCAL_READS: usize = 6;
/// Slack, radians, in deciding that a cone already read holds the one asked for: far above
/// rounding, far below any star's position error.
const HOLD_SLACK: f64 = 1e-9;

/// Catalogue stars with a finite magnitude, addressed by a number that stays the same for the
/// life of the set.
pub(crate) struct RefSet<'a> {
    state: RefCell<State>,
    /// Where a set that grows on demand reads from. `None`: the set was read in one go.
    source: Option<Source<'a>>,
}

/// The cone a growing set stands for, and the catalogue it is read from.
struct Source<'a> {
    catalog: &'a dyn Catalog,
    query: StarQuery,
    centre: SkyPos,
    centre_vec: [f64; 3],
    radius: f64,
}

#[derive(Default)]
struct State {
    stars: Vec<Star>,
    /// `stars[k].pos.to_unit_vector()`.
    vecs: Vec<[f64; 3]>,
    /// Of a set that grows on demand: the identity of each star, by number.
    keys: Vec<StarKey>,
    /// The number of each star, so that a star met again in a later read keeps it. Filled from
    /// `keys` when a second read happens (`keys[..indexed]` are in it): most solves read once.
    known: Map<StarKey, u32>,
    indexed: usize,
    /// Cones read from the catalogue: on their own (at most [`MAX_LOCAL_READS`]) or the whole
    /// set.
    cones: Vec<Cone>,
    /// Which of `cones` is the whole set, once it has been read.
    whole: Option<usize>,
    /// The part of the whole set around the latest place asked for.
    part: Option<Cone>,
    /// The first catalogue failure: the set is then incomplete and says so ([`RefSet::check`]).
    error: Option<SolveError>,
}

/// A star's identity across reads: its identifier, position and magnitude, and its place among
/// the stars of one read that share all three (a catalogue may list a star twice).
type StarKey = (u64, u64, u64, u32, u32);

/// Which stars [`RefSet::within`] looks through.
#[derive(Clone, Copy)]
enum Held {
    /// One of `State::cones`.
    Cone(usize),
    /// `State::part`.
    Part,
}

/// The stars of one cone, brightest first (equal magnitudes by identifier, then in catalogue
/// order).
struct Cone {
    centre: [f64; 3],
    radius: f64,
    members: Vec<u32>,
}

impl Cone {
    /// Whether every star within `need` radians of `q` is one of this cone's.
    fn holds(&self, q: &[f64; 3], need: f64) -> bool {
        angle(&self.centre, q) + need + HOLD_SLACK <= self.radius
    }
}

/// Squared distance between two unit vectors, in the form a k-d tree search uses.
fn chord2(p: &[f64; 3], q: &[f64; 3]) -> f64 {
    (0..3).map(|k| (p[k] - q[k]) * (p[k] - q[k])).sum()
}

impl State {
    /// Add the stars of one read and return their numbers, brightest first. `keep` drops stars
    /// outside the set. With `again` the set grows on demand: a star already in it keeps its
    /// number.
    fn add(&mut self, stars: Vec<Star>, keep: impl Fn(&[f64; 3]) -> bool, again: bool) -> Vec<u32> {
        let mut ranked: Vec<(f32, Star, [f64; 3])> = stars
            .into_iter()
            .filter_map(|s| {
                let mag = s.mag.filter(|m| m.is_finite())?;
                let v = s.pos.to_unit_vector().0;
                keep(&v).then_some((mag, s, v))
            })
            .collect();
        ranked.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.id.cmp(&b.1.id)));
        // Nothing read before: every star is new, and no star needs looking up.
        let first = self.stars.is_empty();
        if again && !first {
            for k in self.indexed..self.keys.len() {
                self.known.insert(self.keys[k], k as u32);
            }
        }
        let place = |s: &Star| (s.pos.ra.to_bits(), s.pos.dec.to_bits());
        let mut members = Vec::with_capacity(ranked.len());
        // Stars of equal magnitude and identifier are together in `ranked`: `group` is where
        // the current run of them starts.
        let mut group = 0;
        for (i, &(mag, star, v)) in ranked.iter().enumerate() {
            let number = if again {
                if (ranked[group].0.to_bits(), ranked[group].1.id) != (mag.to_bits(), star.id) {
                    group = i;
                }
                let nth = ranked[group..i]
                    .iter()
                    .filter(|r| place(&r.1) == place(&star))
                    .count() as u32;
                let (ra, dec) = place(&star);
                let key = (star.id, ra, dec, mag.to_bits(), nth);
                match (!first).then(|| self.known.get(&key).copied()).flatten() {
                    Some(k) => k,
                    None => {
                        self.stars.push(star);
                        self.vecs.push(v);
                        self.keys.push(key);
                        let k = (self.stars.len() - 1) as u32;
                        if !first {
                            self.known.insert(key, k);
                        }
                        k
                    }
                }
            } else {
                self.stars.push(star);
                self.vecs.push(v);
                (self.stars.len() - 1) as u32
            };
            members.push(number);
        }
        if again && !first {
            self.indexed = self.keys.len();
        }
        members
    }
}

impl<'a> RefSet<'a> {
    /// The stars of a cone, read now.
    pub(crate) fn from_cone(
        catalog: &dyn Catalog,
        centre: SkyPos,
        radius: f64,
        query: &StarQuery,
    ) -> Result<Self, SolveError> {
        // A cone wider than the sky is the whole sky.
        let radius = radius.min(std::f64::consts::PI - 1e-9);
        let mut state = State::default();
        let members = state.add(catalog.cone(centre, radius, query)?, |_| true, false);
        state.cones.push(Cone {
            centre: centre.to_unit_vector().0,
            radius,
            members,
        });
        state.whole = Some(0);
        Ok(Self {
            state: RefCell::new(state),
            source: None,
        })
    }

    /// The stars of a cone, read from `catalog` as [`RefSet::within`] is asked for them.
    pub(crate) fn on_demand(
        catalog: &'a dyn Catalog,
        centre: SkyPos,
        radius: f64,
        query: StarQuery,
    ) -> Self {
        Self {
            state: RefCell::new(State::default()),
            source: Some(Source {
                catalog,
                query,
                centre,
                centre_vec: centre.to_unit_vector().0,
                radius: radius.min(std::f64::consts::PI - 1e-9),
            }),
        }
    }

    /// Star number `k`.
    pub(crate) fn star(&self, k: usize) -> Star {
        self.state.borrow().stars[k]
    }

    /// Unit vector of star number `k`: `star(k).pos.to_unit_vector()`.
    pub(crate) fn vec(&self, k: usize) -> [f64; 3] {
        self.state.borrow().vecs[k]
    }

    /// The failure of a catalogue read, if one happened since the last call. The stars handed
    /// out after it were incomplete, so whatever was concluded from them does not stand.
    pub(crate) fn check(&self) -> Result<(), SolveError> {
        match self.state.borrow_mut().error.take() {
            Some(e) => Err(e),
            None => Ok(()),
        }
    }

    /// Hand `visit` the set's stars within `radius` radians of `centre`, brightest first: the
    /// star's number, the star, and its unit vector. `visit` returns `false` to stop.
    ///
    /// "Within" is by chord length between unit vectors, the test a k-d tree over the whole set
    /// would make.
    pub(crate) fn within(
        &self,
        centre: SkyPos,
        radius: f64,
        mut visit: impl FnMut(usize, &Star, &[f64; 3]) -> bool,
    ) {
        let q = centre.to_unit_vector().0;
        let Some(held) = self.hold(centre, &q, radius) else {
            return;
        };
        let reach = chord(radius);
        let reach2 = reach * reach;
        let state = self.state.borrow();
        let cone = match held {
            Held::Cone(i) => &state.cones[i],
            Held::Part => state.part.as_ref().expect("made by hold"),
        };
        for &k in &cone.members {
            let k = k as usize;
            let v = &state.vecs[k];
            if chord2(v, &q) <= reach2 && !visit(k, &state.stars[k], v) {
                break;
            }
        }
    }

    /// Stars in memory that include every star of the set within `need` of `q`, read from the
    /// catalogue if they are not there yet. `None` after a catalogue failure.
    fn hold(&self, centre: SkyPos, q: &[f64; 3], need: f64) -> Option<Held> {
        let mut state = self.state.borrow_mut();
        let Some(source) = &self.source else {
            return state.whole.map(Held::Cone);
        };
        if state.error.is_some() {
            return None;
        }
        // The latest first: the refits of one hypothesis ask for nearly the same place.
        if state.part.as_ref().is_some_and(|p| p.holds(q, need)) {
            return Some(Held::Part);
        }
        if let Some(i) = (0..state.cones.len())
            .rev()
            .find(|&i| state.whole != Some(i) && state.cones[i].holds(q, need))
        {
            return Some(Held::Cone(i));
        }
        let local = need * LOCAL_GROW + LOCAL_PAD;
        // Not for a radius that is no radius (negative, not a number): the whole set answers
        // that as a search of all of it would.
        let small = need >= 0.0 && local < LOCAL_MAX_FRACTION * source.radius;
        if state.whole.is_none() {
            let own = small && state.cones.len() < MAX_LOCAL_READS;
            let (at, radius) = if own {
                (centre, local)
            } else {
                (source.centre, source.radius)
            };
            let stars = match source.catalog.cone(at, radius, &source.query) {
                Ok(stars) => stars,
                Err(e) => {
                    state.error = Some(e.into());
                    return None;
                }
            };
            // A cone reaching past the set's edge brings stars that are not the set's.
            let inside =
                !own || angle(&source.centre_vec, q) + radius + HOLD_SLACK <= source.radius;
            let edge = UnitVector(source.centre_vec);
            let members = state.add(
                stars,
                |v| {
                    // The catalogue's own cone test: "not farther than the radius".
                    let outside = edge.angle_to(UnitVector(*v)) > source.radius;
                    inside || !outside
                },
                true,
            );
            state.cones.push(Cone {
                centre: at.to_unit_vector().0,
                radius,
                members,
            });
            let read = state.cones.len() - 1;
            if own {
                return Some(Held::Cone(read));
            }
            state.whole = Some(read);
        }
        let whole = state.whole?;
        if !small {
            return Some(Held::Cone(whole));
        }
        // The whole set is in memory: take this place's part of it, so that a hypothesis's
        // refits look through a frame's worth of stars and not the prior's.
        let reach = chord(local);
        let reach2 = reach * reach;
        let members: Vec<u32> = state.cones[whole]
            .members
            .iter()
            .copied()
            .filter(|&k| chord2(&state.vecs[k as usize], q) <= reach2)
            .collect();
        state.part = Some(Cone {
            centre: *q,
            radius: local,
            members,
        });
        Some(Held::Part)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fayde_catalog::memory::MemoryCatalog;

    /// A grid of stars around (150°, 2°), every third one listed twice.
    fn catalog() -> MemoryCatalog {
        let mut stars = Vec::new();
        for i in 0..60i64 {
            for j in 0..60i64 {
                let star = Star {
                    id: (i * 60 + j) as u64,
                    pos: SkyPos::from_degrees(
                        150.0 + (i - 30) as f64 * 0.11,
                        2.0 + (j - 30) as f64 * 0.11,
                    )
                    .unwrap(),
                    mag: (!(i == 7 && j == 9))
                        .then_some(8.0 + ((i * 37 + j * 91) % 800) as f32 / 100.0),
                    pos_err_mas: [20.0, 20.0],
                    pm_mas_yr: [None, None],
                };
                stars.push(star);
                if (i + j) % 3 == 0 {
                    stars.push(star);
                }
            }
        }
        MemoryCatalog::new(stars)
    }

    fn seen(set: &RefSet, centre: SkyPos, radius: f64, most: usize) -> Vec<(u64, u64)> {
        let mut out = Vec::new();
        set.within(centre, radius, |_, star, v| {
            assert_eq!(*v, star.pos.to_unit_vector().0);
            out.push((star.id, star.pos.ra.to_bits()));
            out.len() < most
        });
        out
    }

    /// A set read as it is asked for hands out the stars of the set read whole, in the same
    /// order, wherever it is asked: inside, across its edge, after many reads, and with a
    /// radius of any size.
    #[test]
    fn on_demand_equals_whole() {
        let cat = catalog();
        let centre = SkyPos::from_degrees(150.0, 2.0).unwrap();
        let radius = 2.5f64.to_radians();
        let query = StarQuery::default();
        let whole = RefSet::from_cone(&cat, centre, radius, &query).unwrap();
        let lazy = RefSet::on_demand(&cat, centre, radius, query);
        let mut state = 12345u64;
        let mut rnd = || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (state >> 11) as f64 / (1u64 << 53) as f64
        };
        let mut numbers: std::collections::HashMap<(u64, u64), Vec<usize>> = Default::default();
        for round in 0..400 {
            // Mostly small cones, some reaching past the set's edge, some larger than the set.
            // The first rounds are small: they are answered by reads of their own.
            let at = SkyPos::from_degrees(150.0 + 6.4 * (rnd() - 0.5), 2.0 + 6.4 * (rnd() - 0.5))
                .unwrap();
            let r = match round % 9 {
                8 => 4.0 * rnd(),
                7 => 1e-7,
                _ => 0.05 + 0.5 * rnd(),
            }
            .to_radians();
            if round == 3 {
                let state = lazy.state.borrow();
                assert!(state.whole.is_none() && !state.cones.is_empty());
            }
            let most = if round % 4 == 0 { 7 } else { usize::MAX };
            assert_eq!(
                seen(&lazy, at, r, most),
                seen(&whole, at, r, most),
                "round {round}"
            );
            // A star keeps its number however often it is read.
            lazy.within(at, r, |k, star, _| {
                assert_eq!(lazy.star(k), *star);
                let ks = numbers.entry((star.id, star.pos.ra.to_bits())).or_default();
                if !ks.contains(&k) {
                    ks.push(k);
                }
                true
            });
        }
        // Stars listed once have one number, stars listed twice have two.
        assert!(numbers.values().all(|ks| ks.len() <= 2));
        assert!(numbers.values().any(|ks| ks.len() == 2));
        assert!(lazy.check().is_ok());
        let state = lazy.state.borrow();
        assert_eq!(state.cones.len(), MAX_LOCAL_READS + 1);
        assert!(state.whole.is_some() && state.part.is_some());
        // The star without a magnitude is in neither set.
        assert!(state.stars.iter().all(|s| s.mag.is_some()));
    }

    /// A radius that is no radius (not a number, negative, wider than the sky) is answered as
    /// a search of the whole set answers it, without a failed catalogue read.
    #[test]
    fn odd_radii_are_answered_as_by_the_whole_set() {
        let cat = catalog();
        let centre = SkyPos::from_degrees(150.0, 2.0).unwrap();
        let lazy = RefSet::on_demand(&cat, centre, 1f64.to_radians(), StarQuery::default());
        let whole =
            RefSet::from_cone(&cat, centre, 1f64.to_radians(), &StarQuery::default()).unwrap();
        for radius in [f64::NAN, f64::INFINITY, -1.0, 0.0, 7.0] {
            assert_eq!(
                seen(&lazy, centre, radius, usize::MAX),
                seen(&whole, centre, radius, usize::MAX),
                "radius {radius}"
            );
        }
        assert!(seen(&lazy, centre, f64::NAN, usize::MAX).is_empty());
        assert!(lazy.check().is_ok());
    }
}
