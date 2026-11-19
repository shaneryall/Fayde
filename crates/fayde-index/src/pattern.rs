use fayde_core::Portable;
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Shape {
    pub ratios: [f64; 5],
    pub largest: f64,
}

pub fn shape_of<const D: usize>(
    points: &[[f64; D]; 4],
    dist: impl Fn(&[f64; D], &[f64; D]) -> f64,
) -> (Shape, [usize; 4]) {
    let shape = edge_shape(points, dist);
    let mut centroid = [0.0; D];
    for p in points {
        for d in 0..D {
            centroid[d] += p[d] / 4.0;
        }
    }
    let mut order = [0usize, 1, 2, 3];
    let r = |i: usize| -> f64 { (0..D).map(|d| (points[i][d] - centroid[d]).ppowi(2)).sum() };
    order.sort_by(|&a, &b| r(a).total_cmp(&r(b)).then(a.cmp(&b)));
    (shape, order)
}

pub(crate) fn edge_shape<const D: usize>(
    points: &[[f64; D]; 4],
    dist: impl Fn(&[f64; D], &[f64; D]) -> f64,
) -> Shape {
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
    Shape { ratios, largest }
}

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

pub fn planar(a: &[f64; 2], b: &[f64; 2]) -> f64 {
    (a[0] - b[0]).phypot(a[1] - b[1])
}

pub const MAX_BINS: u32 = 7131;

pub fn bins_for(max_error: f64) -> u32 {
    ((0.25 / max_error).round() as u32).max(1)
}

fn bin(r: f64, bins: u32) -> u64 {
    ((r * f64::from(bins)).floor().max(0.0) as u64).min(u64::from(bins) - 1)
}

pub fn key(ratios: &[f64; 5], bins: u32) -> u64 {
    pack_key(ratios.map(|r| bin(r, bins)), bins)
}

fn pack_key(digits: [u64; 5], bins: u32) -> u64 {
    digits
        .into_iter()
        .rev()
        .fold(0, |key, digit| key * u64::from(bins) + digit)
}

pub(crate) struct KeyWindow {
    bounds: [(u64, u64); 5],
    bins: u32,
}

impl KeyWindow {
    pub(crate) fn new(ratios: &[f64; 5], tol: f64, bins: u32) -> Self {
        Self {
            bounds: ratios.map(|r| (bin(r - tol, bins), bin(r + tol, bins))),
            bins,
        }
    }

    pub(crate) fn len(&self) -> u64 {
        self.bounds
            .iter()
            .map(|&(lo, hi)| hi.checked_sub(lo).map_or(0, |d| d + 1))
            .product()
    }

    pub(crate) fn contains(&self, ratios: &[f64; 5]) -> bool {
        self.bounds
            .iter()
            .zip(ratios)
            .all(|(&(lo, hi), &ratio)| (lo..=hi).contains(&bin(ratio, self.bins)))
    }

    pub(crate) fn keys(&self) -> impl Iterator<Item = u64> + '_ {
        self.digits().map(|digits| pack_key(digits, self.bins))
    }

    fn digits(&self) -> impl Iterator<Item = [u64; 5]> + '_ {
        let mut next = self
            .bounds
            .iter()
            .all(|&(lo, hi)| lo <= hi)
            .then(|| self.bounds.map(|(lo, _)| lo));
        std::iter::from_fn(move || {
            let digits = next.take()?;
            let mut following = digits;
            for (digit, &(lo, hi)) in following.iter_mut().zip(&self.bounds) {
                if *digit < hi {
                    *digit += 1;
                    next = Some(following);
                    break;
                }
                *digit = lo;
            }
            Some(digits)
        })
    }
}

pub fn keys_within(ratios: &[f64; 5], tol: f64, bins: u32) -> Vec<u64> {
    let window = KeyWindow::new(ratios, tol, bins);
    let centre = ratios.map(|r| bin(r, bins));
    let mut ranked: Vec<_> = window
        .digits()
        .map(|digits| {
            let distance: u64 = digits
                .into_iter()
                .zip(centre)
                .map(|(digit, at)| digit.abs_diff(at))
                .sum();
            (distance, pack_key(digits, bins))
        })
        .collect();
    ranked.sort_unstable();
    ranked.into_iter().map(|(_, key)| key).collect()
}

const MAGIC_MULT: u128 = 2_654_435_761;

fn home(key: u64, size: usize) -> usize {
    ((u128::from(key) * MAGIC_MULT) % size as u128) as usize
}

pub fn key_check(key: u64) -> u16 {
    (key.wrapping_mul(0x9E37_79B9_7F4A_7C15) >> 48) as u16
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Slot {
    pub stars: [u32; 4],
    pub check: u16,
    pub largest: f32,
}

impl Slot {
    pub const EMPTY: Slot = Slot {
        stars: [u32::MAX; 4],
        check: 0,
        largest: 0.0,
    };

    pub fn is_empty(&self) -> bool {
        self.stars[0] == u32::MAX
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct PatternTable {
    slots: Vec<Slot>,
    len: usize,
}

impl PatternTable {
    pub fn for_patterns(n: usize) -> Self {
        Self {
            slots: vec![Slot::EMPTY; next_prime(3 * n.max(1) + 1)],
            len: 0,
        }
    }

    pub fn from_slots(slots: Vec<Slot>) -> Self {
        let len = slots.iter().filter(|s| !s.is_empty()).count();
        Self { slots, len }
    }

    pub fn slots(&self) -> &[Slot] {
        &self.slots
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

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

    pub fn probe(&self, key: u64) -> impl Iterator<Item = &Slot> + '_ {
        let check = key_check(key);
        self.probe_chain(key).filter(move |s| s.check == check)
    }

    pub(crate) fn probe_chain(&self, key: u64) -> impl Iterator<Item = &Slot> + '_ {
        let size = self.slots.len();
        let start = if size == 0 { 0 } else { home(key, size) };
        (0..size)
            .map(move |k| &self.slots[(start + k) % size])
            .take_while(|s| !s.is_empty())
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
