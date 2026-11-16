//! Four-star patterns: shape, canonical star order, hash keys, and the pattern hash table.
//!
//! Sources: Brown & Stubis 2017 (Tetra) and ESA tetra3 for the pattern shape, key and canonical
//! order; cedar-solve for the table layout (linear probing, key check, largest-edge prune);
//! Knuth, TAOCP vol. 3 §6.4, for multiplicative hashing. See `docs/references.md`.
//!
//! Following tetra3 / cedar-solve: a pattern's shape is its 6 pairwise angular distances,
//! sorted; dividing the 5 smaller ones by the largest gives 5 ratios in (0, 1] that do not
//! change with rotation, translation or scale. The ratios are binned into an integer key, and
//! the key is hashed into an open-addressed table with linear probing (cedar-solve). Each slot
//! also stores a 16-bit check derived from the key, so most probe collisions are rejected
//! without touching the star table, and the largest edge, so a scale prior can prune by size.

use fayde_core::Portable;
/// Shape of a 4-star pattern.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Shape {
    /// The 5 smaller edges divided by the largest, ascending.
    pub ratios: [f64; 5],
    /// Largest edge (radians on the sky, or tangent-plane radians for image patterns).
    pub largest: f64,
}

/// Shape of a pattern and its stars in canonical order (ascending distance from the pattern
/// centroid, ties by input index). Canonical order lets a catalogue pattern and an image pattern
/// be matched star to star.
pub fn shape_of<const D: usize>(
    points: &[[f64; D]; 4],
    dist: impl Fn(&[f64; D], &[f64; D]) -> f64,
) -> (Shape, [usize; 4]) {
    let mut edges = [0.0f64; 6];
    let mut k = 0;
    for i in 0..4 {
        for j in i + 1..4 {
            edges[k] = dist(&points[i], &points[j]);
            k += 1;
        }
    }
    edges.sort_by(f64::total_cmp);
    let largest = edges[5];
    let mut ratios = [0.0; 5];
    for i in 0..5 {
        ratios[i] = if largest > 0.0 {
            edges[i] / largest
        } else {
            0.0
        };
    }
    let mut centroid = [0.0; D];
    for p in points {
        for d in 0..D {
            centroid[d] += p[d] / 4.0;
        }
    }
    let mut order = [0usize, 1, 2, 3];
    let r = |i: usize| -> f64 { (0..D).map(|d| (points[i][d] - centroid[d]).ppowi(2)).sum() };
    order.sort_by(|&a, &b| r(a).total_cmp(&r(b)).then(a.cmp(&b)));
    (Shape { ratios, largest }, order)
}

/// Angle between unit vectors (accurate at all separations).
pub fn angle(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    let c = [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ];
    (c[0] * c[0] + c[1] * c[1] + c[2] * c[2])
        .sqrt()
        .patan2(a[0] * b[0] + a[1] * b[1] + a[2] * b[2])
}

/// Euclidean distance in a plane.
pub fn planar(a: &[f64; 2], b: &[f64; 2]) -> f64 {
    (a[0] - b[0]).phypot(a[1] - b[1])
}

/// Largest bin count whose 5-digit keys fit in a `u64` (`bins^5 < 2^64`).
pub const MAX_BINS: u32 = 7131;

/// Number of bins per ratio for a ratio tolerance `max_error` (tetra3: `1 / (4·max_error)`).
pub fn bins_for(max_error: f64) -> u32 {
    ((0.25 / max_error).round() as u32).max(1)
}

fn bin(r: f64, bins: u32) -> u64 {
    ((r * f64::from(bins)).floor().max(0.0) as u64).min(u64::from(bins) - 1)
}

/// Key of a shape: the binned ratios as digits base `bins`.
pub fn key(ratios: &[f64; 5], bins: u32) -> u64 {
    ratios
        .iter()
        .rev()
        .fold(0u64, |acc, &r| acc * u64::from(bins) + bin(r, bins))
}

/// All keys whose bins intersect `ratios ± tol`, in ascending order of distance from the
/// central key (nearest bins first, as in cedar-solve), ties by key.
pub fn keys_within(ratios: &[f64; 5], tol: f64, bins: u32) -> Vec<u64> {
    let ranges: Vec<(u64, u64)> = ratios
        .iter()
        .map(|&r| (bin(r - tol, bins), bin(r + tol, bins)))
        .collect();
    let centre: Vec<u64> = ratios.iter().map(|&r| bin(r, bins)).collect();
    let mut out: Vec<(u64, u64)> = vec![(0, 0)];
    for i in (0..5).rev() {
        let (lo, hi) = ranges[i];
        let mut next = Vec::with_capacity(out.len() * (hi - lo + 1) as usize);
        for &(k, d) in &out {
            for b in lo..=hi {
                next.push((k * u64::from(bins) + b, d + b.abs_diff(centre[i])));
            }
        }
        out = next;
    }
    out.sort_by_key(|&(k, d)| (d, k));
    out.into_iter().map(|(k, _)| k).collect()
}

const MAGIC_MULT: u128 = 2_654_435_761;

fn home(key: u64, size: usize) -> usize {
    ((u128::from(key) * MAGIC_MULT) % size as u128) as usize
}

/// 16-bit check stored with each slot.
pub fn key_check(key: u64) -> u16 {
    (key.wrapping_mul(0x9E37_79B9_7F4A_7C15) >> 48) as u16
}

/// One table slot.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Slot {
    /// Index-star numbers in canonical order. `stars[0] == u32::MAX` marks an empty slot.
    pub stars: [u32; 4],
    /// [`key_check`] of the pattern's key.
    pub check: u16,
    /// Largest edge, radians.
    pub largest: f32,
}

impl Slot {
    /// An empty slot.
    pub const EMPTY: Slot = Slot {
        stars: [u32::MAX; 4],
        check: 0,
        largest: 0.0,
    };

    /// Whether this slot is empty.
    pub fn is_empty(&self) -> bool {
        self.stars[0] == u32::MAX
    }
}

/// Open-addressed pattern table.
#[derive(Debug, Clone, PartialEq)]
pub struct PatternTable {
    slots: Vec<Slot>,
    len: usize,
}

impl PatternTable {
    /// A table sized for `n` patterns: the next prime above `3n` (load factor ≤ 1/3).
    pub fn for_patterns(n: usize) -> Self {
        Self {
            slots: vec![Slot::EMPTY; next_prime(3 * n.max(1) + 1)],
            len: 0,
        }
    }

    /// From stored slots (index file).
    pub fn from_slots(slots: Vec<Slot>) -> Self {
        let len = slots.iter().filter(|s| !s.is_empty()).count();
        Self { slots, len }
    }

    /// Slots.
    pub fn slots(&self) -> &[Slot] {
        &self.slots
    }

    /// Number of stored patterns.
    pub fn len(&self) -> usize {
        self.len
    }

    /// Whether no patterns are stored.
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Insert a pattern under `key`. Panics if the table is full (cannot happen when sized with
    /// [`PatternTable::for_patterns`]).
    pub fn insert(&mut self, key: u64, stars: [u32; 4], largest: f64) {
        assert!(self.len < self.slots.len(), "pattern table full");
        let size = self.slots.len();
        let mut i = home(key, size);
        while !self.slots[i].is_empty() {
            i = (i + 1) % size;
        }
        self.slots[i] = Slot {
            stars,
            check: key_check(key),
            largest: largest as f32,
        };
        self.len += 1;
    }

    /// Slots stored under `key` (plus rare 16-bit check collisions; callers verify the shape).
    pub fn probe(&self, key: u64) -> impl Iterator<Item = &Slot> + '_ {
        let size = self.slots.len();
        let check = key_check(key);
        // An empty table has no home slot (and `0..size` below yields nothing).
        let start = if size == 0 { 0 } else { home(key, size) };
        (0..size)
            .map(move |k| &self.slots[(start + k) % size])
            .take_while(|s| !s.is_empty())
            .filter(move |s| s.check == check)
    }
}

fn next_prime(mut n: usize) -> usize {
    let is_prime = |n: usize| {
        if n < 2 {
            return false;
        }
        let mut d = 2;
        while d * d <= n {
            if n % d == 0 {
                return false;
            }
            d += 1;
        }
        true
    };
    while !is_prime(n) {
        n += 1;
    }
    n
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn an_empty_table_probes_without_panicking() {
        let t = PatternTable::from_slots(Vec::new());
        assert_eq!(t.probe(12345).count(), 0);
    }

    fn rot(v: [f64; 3], axis: [f64; 3], a: f64) -> [f64; 3] {
        // Rodrigues rotation.
        let (s, c) = a.psin_cos();
        let d = axis[0] * v[0] + axis[1] * v[1] + axis[2] * v[2];
        let x = [
            axis[1] * v[2] - axis[2] * v[1],
            axis[2] * v[0] - axis[0] * v[2],
            axis[0] * v[1] - axis[1] * v[0],
        ];
        [0, 1, 2].map(|i| v[i] * c + x[i] * s + axis[i] * d * (1.0 - c))
    }

    fn unit(ra: f64, dec: f64) -> [f64; 3] {
        let (ra, dec) = (ra.to_radians(), dec.to_radians());
        [dec.pcos() * ra.pcos(), dec.pcos() * ra.psin(), dec.psin()]
    }

    #[test]
    fn shape_is_rotation_invariant_and_order_follows_stars() {
        let p = [
            unit(10.0, 5.0),
            unit(10.3, 5.1),
            unit(10.1, 5.4),
            unit(10.25, 4.85),
        ];
        let (s1, o1) = shape_of(&p, angle);
        let axis = {
            let a = [0.3f64, -0.5, 0.8];
            let n = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
            a.map(|x| x / n)
        };
        let q = p.map(|v| rot(v, axis, 1.234));
        let (s2, o2) = shape_of(&q, angle);
        for i in 0..5 {
            assert!((s1.ratios[i] - s2.ratios[i]).abs() < 1e-12);
        }
        assert!((s1.largest - s2.largest).abs() < 1e-14);
        assert_eq!(o1, o2);
        // Shuffle input order: canonical order still names the same stars.
        let shuffled = [p[2], p[0], p[3], p[1]];
        let (_, o3) = shape_of(&shuffled, angle);
        let names = [2, 0, 3, 1];
        assert_eq!(o3.map(|i| names[i]), o1);
    }

    #[test]
    fn keys_within_contains_true_key_nearest_first() {
        let bins = bins_for(0.005);
        assert_eq!(bins, 50);
        let r = [0.21, 0.33, 0.501, 0.62, 0.9];
        let k = key(&r, bins);
        let ks = keys_within(&r, 0.004, bins);
        assert_eq!(ks[0], k);
        assert!(ks.len() >= 2 && ks.len() <= 32);
        let mut sorted = ks.clone();
        sorted.dedup();
        assert_eq!(sorted.len(), ks.len());
    }

    #[test]
    fn table_insert_probe() {
        let mut t = PatternTable::for_patterns(1000);
        for k in 0..1000u64 {
            t.insert(k * 7919, [k as u32, 1, 2, 3], 0.01);
        }
        assert_eq!(t.len(), 1000);
        for k in 0..1000u64 {
            let hits: Vec<_> = t.probe(k * 7919).collect();
            assert!(hits.iter().any(|s| s.stars[0] == k as u32), "key {k}");
        }
        assert_eq!(next_prime(3001), 3001);
        assert_eq!(next_prime(3002), 3011);
    }
}
