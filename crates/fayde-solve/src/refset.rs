use std::cell::RefCell;

use fayde_catalog::{Catalog, Star, StarQuery};
use fayde_core::{SkyPos, UnitVector};
use fayde_index::{kdtree::chord, pattern::angle};

use crate::{hash::Map, request::SolveError};

const LOCAL_GROW: f64 = 1.25;
const LOCAL_PAD: f64 = 1e-6;
const LOCAL_MAX_FRACTION: f64 = 0.5;
const MAX_LOCAL_READS: usize = 6;
const HOLD_SLACK: f64 = 1e-9;

pub(crate) struct RefSet<'a> {
    state: RefCell<State>,
    source: Option<Source<'a>>,
}

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
    vecs: Vec<[f64; 3]>,
    keys: Vec<StarKey>,
    known: Map<StarKey, u32>,
    indexed: usize,
    cones: Vec<Cone>,
    whole: Option<usize>,
    part: Option<Cone>,
    error: Option<SolveError>,
}

type StarKey = (u64, u64, u64, u32, u32);

#[derive(Clone, Copy)]
enum Held {
    Cone(usize),
    Part,
}

struct Cone {
    centre: [f64; 3],
    radius: f64,
    members: Vec<u32>,
}

impl Cone {
    fn holds(&self, q: &[f64; 3], need: f64) -> bool {
        angle(&self.centre, q) + need + HOLD_SLACK <= self.radius
    }
}

fn chord2(p: &[f64; 3], q: &[f64; 3]) -> f64 {
    (0..3).map(|k| (p[k] - q[k]) * (p[k] - q[k])).sum()
}

impl State {
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
        let first = self.stars.is_empty();
        if again && !first {
            for k in self.indexed..self.keys.len() {
                self.known.insert(self.keys[k], k as u32);
            }
        }
        let place = |s: &Star| (s.pos.ra.to_bits(), s.pos.dec.to_bits());
        let mut members = Vec::with_capacity(ranked.len());
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
    pub(crate) fn from_cone(
        catalog: &dyn Catalog,
        centre: SkyPos,
        radius: f64,
        query: &StarQuery,
    ) -> Result<Self, SolveError> {
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

    pub(crate) fn star(&self, k: usize) -> Star {
        self.state.borrow().stars[k]
    }

    pub(crate) fn vec(&self, k: usize) -> [f64; 3] {
        self.state.borrow().vecs[k]
    }

    pub(crate) fn check(&self) -> Result<(), SolveError> {
        match self.state.borrow_mut().error.take() {
            Some(e) => Err(e),
            None => Ok(()),
        }
    }

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

    fn hold(&self, centre: SkyPos, q: &[f64; 3], need: f64) -> Option<Held> {
        let mut state = self.state.borrow_mut();
        let Some(source) = &self.source else {
            return state.whole.map(Held::Cone);
        };
        if state.error.is_some() {
            return None;
        }
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
            let inside =
                !own || angle(&source.centre_vec, q) + radius + HOLD_SLACK <= source.radius;
            let edge = UnitVector(source.centre_vec);
            let members = state.add(
                stars,
                |v| {
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
            lazy.within(at, r, |k, star, _| {
                assert_eq!(lazy.star(k), *star);
                let ks = numbers.entry((star.id, star.pos.ra.to_bits())).or_default();
                if !ks.contains(&k) {
                    ks.push(k);
                }
                true
            });
        }
        assert!(numbers.values().all(|ks| ks.len() <= 2));
        assert!(numbers.values().any(|ks| ks.len() == 2));
        assert!(lazy.check().is_ok());
        let state = lazy.state.borrow();
        assert_eq!(state.cones.len(), MAX_LOCAL_READS + 1);
        assert!(state.whole.is_some() && state.part.is_some());
        assert!(state.stars.iter().all(|s| s.mag.is_some()));
    }

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
