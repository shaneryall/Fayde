//! Radius search among points in a frame, on a uniform grid of cells (the cell method for
//! fixed-radius near neighbours; Bentley & Friedman 1979, ACM Comput. Surv. 11, 397).
//!
//! Matching a source list to the reference stars of one hypothesis searches a new set of a few
//! hundred points a few hundred times, and a hypothesis is refitted some dozens of times. A
//! k-d tree answers the same question, but building one costs a sort per level, where a grid
//! is filled in one pass; and for points spread evenly over a frame a cell holds the few
//! points a tree would reach after descending its whole depth.
//!
//! [`Grid::within`] returns what `KdTree::<2>::within` returns for the same points: the same
//! indices, in the same order.

/// Points per cell aimed for.
const PER_CELL: f64 = 2.0;
/// Most cells along one axis.
const MAX_CELLS: usize = 1024;

/// A uniform grid over points in a plane.
pub(crate) struct Grid {
    points: Vec<[f64; 2]>,
    /// Lowest coordinates of the grid.
    origin: [f64; 2],
    cell: f64,
    /// Cells along x and y.
    cells: [usize; 2],
    /// `items[start[c]..start[c + 1]]` are the points of cell `c`, ascending.
    start: Vec<u32>,
    items: Vec<u32>,
}

impl Grid {
    /// Build from points. Indices returned by [`Grid::within`] refer to this order.
    pub(crate) fn build(points: Vec<[f64; 2]>) -> Self {
        assert!(points.len() < u32::MAX as usize, "too many points");
        // The extent of the finite coordinates (`min` and `max` skip NaN).
        let (mut lo, mut hi) = ([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]);
        for p in &points {
            for k in 0..2 {
                if p[k].is_finite() {
                    lo[k] = lo[k].min(p[k]);
                    hi[k] = hi[k].max(p[k]);
                }
            }
        }
        let size = [0, 1].map(|k| if hi[k] > lo[k] { hi[k] - lo[k] } else { 0.0 });
        let origin = [0, 1].map(|k| if lo[k].is_finite() { lo[k] } else { 0.0 });
        // Square cells holding PER_CELL points each when the points are spread evenly, no
        // smaller than a MAX_CELLS-th of the longer side.
        let even = (size[0] * size[1] * PER_CELL / points.len().max(1) as f64).sqrt();
        let cell = even.max(size[0].max(size[1]) / MAX_CELLS as f64);
        let cell = if cell > 0.0 && cell.is_finite() {
            cell
        } else {
            1.0
        };
        let cells = size.map(|s| ((s / cell) as usize + 1).min(MAX_CELLS));
        let mut grid = Self {
            points,
            origin,
            cell,
            cells,
            start: vec![0; cells[0] * cells[1] + 1],
            items: Vec::new(),
        };
        // Counting sort by cell: within a cell the points stay in ascending order.
        let of: Vec<usize> = grid.points.iter().map(|p| grid.cell_of(p)).collect();
        for &c in &of {
            grid.start[c + 1] += 1;
        }
        for c in 0..cells[0] * cells[1] {
            grid.start[c + 1] += grid.start[c];
        }
        let mut next = grid.start.clone();
        grid.items = vec![0; grid.points.len()];
        for (i, &c) in of.iter().enumerate() {
            grid.items[next[c] as usize] = i as u32;
            next[c] += 1;
        }
        grid
    }

    /// Cell along axis `k` of coordinate `v`, clamped to the grid (a NaN lands in cell 0).
    fn along(&self, k: usize, v: f64) -> usize {
        // `as usize` saturates: below the grid is 0, NaN is 0, beyond it is clamped.
        (((v - self.origin[k]) / self.cell) as usize).min(self.cells[k] - 1)
    }

    fn cell_of(&self, p: &[f64; 2]) -> usize {
        self.along(1, p[1]) * self.cells[0] + self.along(0, p[0])
    }

    /// Indices of all points within Euclidean distance `r` (not negative) of `q`, ascending,
    /// into `out` (cleared first).
    pub(crate) fn within(&self, q: [f64; 2], r: f64, out: &mut Vec<usize>) {
        out.clear();
        let r2 = r * r;
        // The cells the circle reaches, taken a hair wide: the bounds below and the distance
        // test round differently, and a point exactly on the circle must not fall between.
        let reach = |k: usize| r + 1e-9 * (q[k].abs() + r + self.cell);
        let (x0, x1) = (
            self.along(0, q[0] - reach(0)),
            self.along(0, q[0] + reach(0)),
        );
        let (y0, y1) = (
            self.along(1, q[1] - reach(1)),
            self.along(1, q[1] + reach(1)),
        );
        for y in y0..=y1 {
            let row = y * self.cells[0];
            let (from, to) = (self.start[row + x0], self.start[row + x1 + 1]);
            for &i in &self.items[from as usize..to as usize] {
                let p = &self.points[i as usize];
                // The sum a k-d tree search makes, term for term.
                let d2: f64 = (0..2).map(|k| (p[k] - q[k]) * (p[k] - q[k])).sum();
                if d2 <= r2 {
                    out.push(i as usize);
                }
            }
        }
        out.sort_unstable();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fayde_index::kdtree::KdTree;

    fn rng(seed: u64) -> impl FnMut() -> f64 {
        let mut s = seed;
        move || {
            s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (s >> 11) as f64 / (1u64 << 53) as f64
        }
    }

    /// The grid answers every search as the k-d tree does, for points spread over a frame,
    /// bunched, on a line, repeated, and for radii from nothing to everything.
    #[test]
    fn answers_as_the_kd_tree_does() {
        let mut rnd = rng(99);
        let mut sets: Vec<Vec<[f64; 2]>> = vec![
            Vec::new(),
            vec![[3.0, 4.0]],
            vec![[3.0, 4.0]; 9],
            (0..300).map(|i| [i as f64 * 3.3, 500.0]).collect(),
            (0..300).map(|i| [7.0, i as f64 * 0.01]).collect(),
        ];
        for n in [2usize, 17, 300, 2000] {
            sets.push(
                (0..n)
                    .map(|_| [1024.0 * rnd() - 0.5, 768.0 * rnd() - 0.5])
                    .collect(),
            );
            // Most of the points in one corner.
            sets.push(
                (0..n)
                    .map(|i| {
                        let s = if i % 5 == 0 { 1000.0 } else { 20.0 };
                        [s * rnd(), s * rnd()]
                    })
                    .collect(),
            );
        }
        let mut found = Vec::new();
        for (which, points) in sets.into_iter().enumerate() {
            let tree = KdTree::build(points.clone());
            let grid = Grid::build(points.clone());
            for round in 0..300 {
                let q = match round % 4 {
                    0 if !points.is_empty() => points[(rnd() * points.len() as f64) as usize],
                    1 => [3000.0 * rnd() - 1000.0, 3000.0 * rnd() - 1000.0],
                    _ => [1024.0 * rnd(), 768.0 * rnd()],
                };
                let r = match round % 7 {
                    0 => 0.0,
                    1 => 1e-9,
                    2 => 5000.0,
                    3 => f64::INFINITY,
                    _ => 60.0 * rnd() * rnd(),
                };
                grid.within(q, r, &mut found);
                assert_eq!(found, tree.within(q, r), "set {which}, q {q:?}, r {r}");
            }
        }
    }

    /// A radius or a point that is not a number matches nothing, as in the tree.
    #[test]
    fn nan_matches_nothing() {
        let points = vec![[1.0, 2.0], [f64::NAN, 2.0], [5.0, f64::NAN], [4.0, 4.0]];
        let tree = KdTree::build(points.clone());
        let grid = Grid::build(points);
        let mut found = Vec::new();
        for (q, r) in [
            ([1.0, 2.0], f64::NAN),
            ([f64::NAN, 2.0], 10.0),
            ([4.0, 4.0], 100.0),
        ] {
            grid.within(q, r, &mut found);
            assert_eq!(found, tree.within(q, r), "q {q:?}, r {r}");
        }
    }
}
