use fayde_core::Portable;

use crate::grid::Grid;

#[derive(Debug, Clone, Copy)]
pub(crate) struct TestSource {
    pub px: [f64; 2],
    pub sigma: f64,
    pub sigma_defaulted: bool,
    pub distractor: f64,
    pub along: f64,
    pub angle: f64,
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct RefStar {
    pub px: [f64; 2],
    pub star: usize,
    pub sigma_px: f64,
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct Hit {
    pub source: usize,
    pub reference: usize,
    pub inlier_prob: f64,
    pub sigma: f64,
}

#[derive(Debug, Clone)]
pub(crate) struct Assessment {
    pub log_odds: f64,
    pub hits: Vec<Hit>,
}

const SIGMA_FLOOR_PX: f64 = 1e-3;

pub(crate) fn assess(
    sources: &[TestSource],
    refs: &[RefStar],
    area: f64,
    sigma_scale: f64,
    seed: &[usize],
) -> Assessment {
    let n_ref = refs.len().max(1) as f64;
    let grid = Grid::build(refs.iter().map(|r| r.px).collect());
    let mut cands = Vec::new();
    let mut ranked: Vec<(f64, usize)> = Vec::new();
    let mut used = vec![false; refs.len()];
    let mut log_odds = 0.0;
    let mut best = f64::NEG_INFINITY;
    let mut hits = Vec::new();
    let max_ref_sigma = refs.iter().map(|r| r.sigma_px).fold(0.0, f64::max);
    for (i, s) in sources.iter().enumerate() {
        let sigma = (s.sigma
            * if s.sigma_defaulted {
                sigma_scale
            } else {
                sigma_scale.max(1.0)
            })
        .max(SIGMA_FLOOR_PX);
        let sigma_along = sigma.phypot(s.along);
        let (ca, sa) = (s.angle.pcos(), s.angle.psin());
        let split = |k: usize| {
            let (dx, dy) = (refs[k].px[0] - s.px[0], refs[k].px[1] - s.px[1]);
            (dx * ca + dy * sa, -dx * sa + dy * ca)
        };
        let comb = |k: usize| {
            (
                sigma_along.phypot(refs[k].sigma_px),
                sigma.phypot(refs[k].sigma_px),
            )
        };
        let m2 = |k: usize| {
            let (da, dc) = split(k);
            let (a, c) = comb(k);
            (da / a).ppowi(2) + (dc / c).ppowi(2)
        };
        grid.within(s.px, 5.0 * sigma_along.phypot(max_ref_sigma), &mut cands);
        ranked.clear();
        ranked.extend(
            cands
                .iter()
                .map(|&k| (m2(k), k))
                .filter(|&(m, _)| m <= 25.0),
        );
        ranked.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
        match ranked.iter().map(|&(_, k)| k).find(|&k| !used[k]) {
            Some(k) => {
                let (sa_k, sc_k) = comb(k);
                let (da, dc) = split(k);
                let g = (-0.5 * ((da / sa_k).ppowi(2) + (dc / sc_k).ppowi(2))).pexp()
                    / (2.0 * std::f64::consts::PI * sa_k * sc_k);
                let sigma = sc_k;
                let fg = (1.0 - s.distractor) * (area / n_ref) * g;
                let ratio = fg + s.distractor;
                if !seed.contains(&i) {
                    log_odds += ratio.pln();
                }
                used[k] = true;
                hits.push(Hit {
                    source: i,
                    reference: k,
                    inlier_prob: fg / ratio,
                    sigma,
                });
            }
            None if seed.contains(&i) => {}
            None => log_odds += s.distractor.max(1e-12).pln(),
        }
        best = best.max(log_odds);
    }
    Assessment {
        log_odds: if sources.is_empty() { 0.0 } else { best },
        hits,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gates_include_catalogue_error() {
        let refs = [RefStar {
            px: [50.1, 50.0],
            star: 0,
            sigma_px: 0.2,
        }];
        let sources = [TestSource {
            px: [50.0, 50.0],
            sigma: 0.01,
            sigma_defaulted: false,
            distractor: 0.25,
            along: 0.0,
            angle: 0.0,
        }];
        let a = assess(&sources, &refs, 1e4, 1.0, &[]);
        assert_eq!(a.hits.len(), 1);
        assert!(a.hits[0].inlier_prob > 0.9);
    }

    #[test]
    fn only_default_sigmas_shrink_in_mixed_sources() {
        let refs = [
            RefStar {
                px: [10.0, 10.0],
                star: 0,
                sigma_px: 0.2,
            },
            RefStar {
                px: [90.0, 90.0],
                star: 1,
                sigma_px: 0.2,
            },
        ];
        let sources = [
            TestSource {
                px: [10.0, 10.0],
                sigma: 1.0,
                sigma_defaulted: true,
                distractor: 0.25,
                along: 0.0,
                angle: 0.0,
            },
            TestSource {
                px: [90.0, 90.0],
                sigma: 1.0,
                sigma_defaulted: false,
                distractor: 0.25,
                along: 0.0,
                angle: 0.0,
            },
        ];
        let score = assess(&sources, &refs, 10_000.0, 0.1, &[]).log_odds;
        let variance_default = 0.1f64.ppowi(2) + 0.2f64.ppowi(2);
        let variance_supplied = 1.0 + 0.2f64.ppowi(2);
        let odds = |v: f64| (0.75 * 5000.0 / (2.0 * std::f64::consts::PI * v) + 0.25).pln();
        assert!((score - odds(variance_default) - odds(variance_supplied)).abs() < 1e-12);
    }

    #[test]
    fn real_field_scores_high_and_random_scores_low() {
        let refs: Vec<RefStar> = (0..40)
            .map(|k| RefStar {
                px: [(k * 97 % 1000) as f64, (k * 61 % 1000) as f64],
                star: k,
                sigma_px: 0.0,
            })
            .collect();
        let src = |p: [f64; 2]| TestSource {
            px: p,
            sigma: 1.0,
            sigma_defaulted: false,
            distractor: 0.25,
            along: 0.0,
            angle: 0.0,
        };
        let good: Vec<TestSource> = refs
            .iter()
            .map(|r| src([r.px[0] + 0.5, r.px[1] - 0.3]))
            .collect();
        let a = assess(&good, &refs, 1e6, 1.0, &[]);
        assert!(a.log_odds > 1e9f64.pln() * 5.0, "{}", a.log_odds);
        assert_eq!(a.hits.len(), 40);
        let bad: Vec<TestSource> = (0..40)
            .map(|k| src([(k * 13 % 1000) as f64 + 0.5, (k * 29 % 1000) as f64 + 0.5]))
            .collect();
        assert!(assess(&bad, &refs, 1e6, 1.0, &[]).log_odds < 1e9f64.pln() / 2.0);
    }

    #[test]
    fn faint_junk_does_not_outvote_the_stars_and_seed_stars_score_nothing() {
        let refs: Vec<RefStar> = (0..20)
            .map(|k| RefStar {
                px: [(k * 97 % 1000) as f64, (k * 61 % 1000) as f64],
                star: k,
                sigma_px: 0.0,
            })
            .collect();
        let src = |p: [f64; 2]| TestSource {
            px: p,
            sigma: 1.0,
            sigma_defaulted: false,
            distractor: 0.25,
            along: 0.0,
            angle: 0.0,
        };
        let stars: Vec<TestSource> = refs.iter().map(|r| src([r.px[0] + 0.4, r.px[1]])).collect();
        let clean = assess(&stars, &refs, 1e6, 1.0, &[]).log_odds;
        let mut junk = stars.clone();
        junk.extend((0..300).map(|k| src([(k * 13 % 1000) as f64 + 0.5, 990.0])));
        let scored = assess(&junk, &refs, 1e6, 1.0, &[]);
        assert!(
            (scored.log_odds - clean).abs() < 1e-9,
            "{} vs {clean}",
            scored.log_odds
        );
        assert_eq!(scored.hits.len(), 20);
        let seeded = assess(&stars, &refs, 1e6, 1.0, &[0, 1, 2, 3]);
        assert_eq!(seeded.hits.len(), 20);
        assert!((seeded.log_odds - clean * 16.0 / 20.0).abs() < 1e-9);
    }
}
