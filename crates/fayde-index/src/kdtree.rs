use fayde_core::Portable;
#[derive(Debug, Clone)]
pub struct KdTree<const D: usize> {
    points: Vec<[f64; D]>,
    order: Vec<u32>,
}

impl<const D: usize> KdTree<D> {
    pub fn build(points: Vec<[f64; D]>) -> Self {
        assert!(points.len() < u32::MAX as usize, "too many points");
        let mut order: Vec<u32> = (0..points.len() as u32).collect();
        build_rec(&points, &mut order, 0);
        Self { points, order }
    }

    pub fn len(&self) -> usize {
        self.points.len()
    }

    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }

    pub fn point(&self, i: usize) -> [f64; D] {
        self.points[i]
    }

    pub fn within(&self, q: [f64; D], r: f64) -> Vec<usize> {
        let mut out = Vec::new();
        self.within_rec(0, self.order.len(), 0, q, r * r, r, &mut out);
        out.sort_unstable();
        out
    }

    pub fn any_within(&self, q: [f64; D], r: f64) -> bool {
        self.any_rec(0, self.order.len(), 0, q, r * r, r)
    }

    fn any_rec(&self, lo: usize, hi: usize, depth: usize, q: [f64; D], r2: f64, r: f64) -> bool {
        if lo >= hi {
            return false;
        }
        let mid = lo + (hi - lo) / 2;
        let p = self.points[self.order[mid] as usize];
        let d2: f64 = (0..D).map(|k| (p[k] - q[k]) * (p[k] - q[k])).sum();
        if d2 <= r2 {
            return true;
        }
        let diff = q[depth % D] - p[depth % D];
        (diff - r <= 0.0 && self.any_rec(lo, mid, depth + 1, q, r2, r))
            || (diff + r >= 0.0 && self.any_rec(mid + 1, hi, depth + 1, q, r2, r))
    }

    #[allow(clippy::too_many_arguments)]
    fn within_rec(
        &self,
        lo: usize,
        hi: usize,
        depth: usize,
        q: [f64; D],
        r2: f64,
        r: f64,
        out: &mut Vec<usize>,
    ) {
        if lo >= hi {
            return;
        }
        let mid = lo + (hi - lo) / 2;
        let idx = self.order[mid] as usize;
        let p = self.points[idx];
        let d2: f64 = (0..D).map(|k| (p[k] - q[k]) * (p[k] - q[k])).sum();
        if d2 <= r2 {
            out.push(idx);
        }
        let axis = depth % D;
        let diff = q[axis] - p[axis];
        if diff - r <= 0.0 {
            self.within_rec(lo, mid, depth + 1, q, r2, r, out);
        }
        if diff + r >= 0.0 {
            self.within_rec(mid + 1, hi, depth + 1, q, r2, r, out);
        }
    }
}

fn build_rec<const D: usize>(points: &[[f64; D]], order: &mut [u32], depth: usize) {
    if order.len() <= 1 {
        return;
    }
    let axis = depth % D;
    let mid = order.len() / 2;
    order.select_nth_unstable_by(mid, |&a, &b| {
        points[a as usize][axis]
            .total_cmp(&points[b as usize][axis])
            .then(a.cmp(&b))
    });
    let (left, right) = order.split_at_mut(mid);
    build_rec(points, left, depth + 1);
    build_rec(points, &mut right[1..], depth + 1);
}

pub fn chord(angle: f64) -> f64 {
    2.0 * (angle / 2.0).psin()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_brute_force() {
        let mut s = 12345u64;
        let mut rnd = || {
            s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (s >> 11) as f64 / (1u64 << 53) as f64
        };
        let pts: Vec<[f64; 3]> = (0..3000).map(|_| [rnd(), rnd(), rnd()]).collect();
        let tree = KdTree::build(pts.clone());
        for _ in 0..200 {
            let q = [rnd(), rnd(), rnd()];
            let r = 0.15 * rnd();
            let brute: Vec<usize> = (0..pts.len())
                .filter(|&i| (0..3).map(|k| (pts[i][k] - q[k]).ppowi(2)).sum::<f64>() <= r * r)
                .collect();
            assert_eq!(tree.within(q, r), brute);
            assert_eq!(tree.any_within(q, r), !brute.is_empty());
        }
        for r in [0.0, 1e-9, f64::NAN] {
            let q = [rnd(), rnd(), rnd()];
            assert_eq!(tree.any_within(q, r), !tree.within(q, r).is_empty());
        }
        assert!(KdTree::<3>::build(vec![]).within([0.0; 3], 1.0).is_empty());
    }
}
