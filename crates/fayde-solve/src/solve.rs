//! Solving: hypothesis generation, verification and fitting (design doc §5.3–§5.5).
//!
//! Two hypothesis generators share one pipeline:
//!
//! - **Pair seeds** (tight scale prior and a pointing prior): the brightest catalogue stars in
//!   the prior cone, every pair closer than the frame diagonal, sorted by separation; a range
//!   lookup in this list plays the role of Mortari's k-vector (Mortari 1997; Mortari & Neta
//!   2000). An image pair plus a catalogue pair, an assignment and a parity fix a similarity
//!   transform.
//! - **Pattern hash** (loose or missing priors; needs an [`Index`]): 4-star image patterns, shape
//!   lookup in the index (Brown & Stubis 2017; tetra3; cedar-solve), edge-compatible permutations
//!   fitted by four-parameter complex least squares (planar Umeyama 1991), both parities.
//!
//! Shared pipeline: sources ranked by `p_star × flux` (PROSAC; Chum & Matas 2005) → hypothesis →
//! hard prior limits → count hypothesis → cheap pre-filter → local refit (LO-RANSAC; Chum,
//! Matas & Kittler 2003; Cauchy IRLS, Holland & Welsch 1977), with a distortion tune-up when a
//! linear model is not conclusive → log-odds verification (Lang et al. 2010) that leaves out
//! the seed stars, accepted when `log-odds − ln(hypotheses tested)` clears the threshold (NFA)
//! → IRLS final fit to convergence with an error-scale factor; the fitted WCS must pass the
//! same threshold.

use fayde_core::Portable;
use std::{collections::HashMap, f64::consts::PI, time::Instant};

use fayde_catalog::{Catalog, Star, StarQuery};
use fayde_core::{PixelPos, SkyPos, UnitVector};
use fayde_index::{
    image_shape,
    kdtree::{chord, KdTree},
    pattern::angle,
    Index,
};
use fayde_wcs::{gnomonic, gnomonic_inverse, Parity, TanWcs, Wcs};

use crate::{
    fit::{
        fit_best_order, LinearPrior, Pair, WcsCovariance, FINAL_OBS_PER_PARAM, GROWTH_OBS_PER_PARAM,
    },
    request::{
        Match, Outcome, Priors, ScalePrior, Solution, SolveError, SolveMethod, SolveRequest,
        StageTrace, Verification,
    },
    verify::{assess, Assessment, RefStar, TestSource},
};

/// Catalogue stars per frame area kept for seeding, per seed source.
const SEED_STARS_PER_SOURCE: f64 = 3.0;
/// Shortest image pair used as a seed, pixels (short pairs give poor rotation).
const MIN_SEED_PAIR_PX: f64 = 20.0;
/// Sources projected by the pre-filter.
const PREFILTER_SOURCES: usize = 12;
/// A scale prior is "tight" (pair seeds first) when σ/scale is at most this.
const TIGHT_SCALE: f64 = 0.03;
/// Sources closer than this to a better-ranked source are treated as duplicates, pixels.
const DUPLICATE_PX: f64 = 3.0;
/// Shortest edge of an image pattern, pixels.
const MIN_PATTERN_EDGE_PX: f64 = 10.0;
/// Local refits with fewer matches than this fit a similarity (4 parameters) instead of an
/// affine transform (6): a handful of stars cannot constrain shear and anisotropy.
const AFFINE_MIN_MATCHES: usize = 8;
/// Match radius of the distortion tune-up's wide stage, pixels.
const TUNE_RADIUS_PX: f64 = 12.0;
/// Fewest wide-radius matches for which a distortion tune-up is attempted (order 2 has 12
/// parameters; chance hypotheses collect about half this many matches at that radius).
const TUNE_MIN_MATCHES: usize = 12;
/// A linear hypothesis is suspected of unmodelled distortion when the wide radius matches more
/// than this factor (plus 3) times the sources its own match gates do: stars sit just outside
/// the gates. Undistorted frames gain about one chance match at the wide radius.
const TUNE_SUSPECT_FACTOR: f64 = 1.25;
/// Most re-matching passes of the tune-up's wide stage.
const TUNE_WIDE_PASSES: usize = 8;
/// A distortion polynomial is kept only when its confident matches cover the sources: the frame
/// is cut into `SIP_CELLS` × `SIP_CELLS` cells, and every cell holding at least
/// `SIP_CELL_MIN_SOURCES` ranked sources must have a matched fraction of at least
/// `SIP_CELL_MIN_FRACTION` × the frame's overall matched fraction. A polynomial constrained by
/// one part of the frame is wrong by tens of pixels in the rest. Measured on 240 synthetic
/// solves: good fits had a lowest cell ratio of 0.47, fits locked onto part of the frame 0.13 or
/// less.
const SIP_CELLS: usize = 3;
/// See [`SIP_CELLS`].
const SIP_CELL_MIN_SOURCES: usize = 6;
/// See [`SIP_CELLS`].
const SIP_CELL_MIN_FRACTION: f64 = 0.25;
/// Most IRLS iterations of the final fit (it stops earlier once the match set repeats).
const FINAL_FIT_MAX_ITERATIONS: usize = 20;
/// Gate widening for a truncated source when the frame's trail is unknown.
const TRUNCATED_SIGMA_FACTOR: f64 = 4.0;
/// Most passes of a similarity fit while its tangent point moves to the fitted `CRVAL`.
const SIMILARITY_TANGENT_PASSES: usize = 8;

struct Ranked {
    /// Index into the request's sources.
    input: usize,
    px: [f64; 2],
    sigma: f64,
    distractor: f64,
    /// Cut by the frame edge: evidence for a match, never fitted, never a seed.
    truncated: bool,
}

/// Catalogue stars sorted by magnitude, with a tree over their unit vectors.
struct RefSet {
    stars: Vec<Star>,
    vecs: Vec<[f64; 3]>,
    tree: KdTree<3>,
}

impl RefSet {
    fn from_cone(
        catalog: &dyn Catalog,
        centre: SkyPos,
        radius: f64,
        query: &StarQuery,
    ) -> Result<Self, SolveError> {
        // A cone wider than the sky is the whole sky.
        let radius = radius.min(std::f64::consts::PI - 1e-9);
        let mut stars: Vec<Star> = catalog
            .cone(centre, radius, query)?
            .into_iter()
            .filter(|s| s.mag.is_some_and(f32::is_finite))
            .collect();
        stars.sort_by(|a, b| {
            a.mag
                .unwrap()
                .total_cmp(&b.mag.unwrap())
                .then(a.id.cmp(&b.id))
        });
        let vecs: Vec<[f64; 3]> = stars.iter().map(|s| s.pos.to_unit_vector().0).collect();
        let tree = KdTree::build(vecs.clone());
        Ok(Self { stars, vecs, tree })
    }
}

/// A final fit iterated to convergence.
struct Converged {
    wcs: Wcs,
    error_scale: f64,
    covariance: WcsCovariance,
    /// References and their assessment under `wcs` at `error_scale`.
    refs: Vec<RefStar>,
    assessment: Assessment,
    /// Posterior log-odds of `wcs`: data log-odds plus the prior terms.
    log_odds: f64,
}

/// Why the hypothesis loop stopped.
enum Stop {
    Solved(Box<Solution>),
    Budget,
}

/// State shared by every hypothesis of one solve.
struct Ctx<'a> {
    req: &'a SolveRequest,
    catalog: &'a dyn Catalog,
    query: StarQuery,
    ranked: Vec<Ranked>,
    tests: Vec<TestSource>,
    w: f64,
    h: f64,
    crpix: PixelPos,
    n_ref: usize,
    /// Hard scale limits, arcsec/px.
    s_lo: f64,
    s_hi: f64,
    started: Instant,
    method: SolveMethod,
    index_hash: Option<u64>,
    deadline: Option<Instant>,
    tested: u64,
    best: Option<f64>,
    /// Ignore the pointing prior's hard limit (out-of-prior fallback stage).
    outside_pointing: bool,
    /// Catalogue cone covering the whole pointing prior, built on first use.
    prior_cone: Option<std::rc::Rc<RefSet>>,
    /// Stages executed so far.
    trace: Vec<StageTrace>,
    /// Seeds skipped, and refined fits rejected, for lying outside the hard limits: reported
    /// when the search abstains, so a limit that blocks a solution is visible.
    limit_seeds: std::cell::Cell<u64>,
    limit_fits: std::cell::Cell<u64>,
}

impl Ctx<'_> {
    fn half_diag_px(&self) -> f64 {
        0.5 * self.w.phypot(self.h)
    }

    fn out_of_budget(&self) -> bool {
        self.tested >= self.req.policy.max_hypotheses || self.past_deadline()
    }

    fn past_deadline(&self) -> bool {
        self.deadline.is_some_and(|d| Instant::now() > d)
    }

    /// A refit stopped by its deadline (a refit has no hypothesis limit): traced and reported.
    fn refit_deadline(&mut self) -> Outcome {
        let stage = "refit around the given WCS";
        self.trace.push(StageTrace {
            stage: stage.into(),
            hypotheses: 1,
            seconds: self.started.elapsed().as_secs_f64(),
            outcome: "deadline",
        });
        self.abstain(&format!(
            "deadline reached during {stage}, {:.2} s after the start",
            self.started.elapsed().as_secs_f64()
        ))
    }

    /// Which budget stopped the search: the trace outcome and the abstention reason.
    fn budget_stop(&self, stage: &str) -> (&'static str, String) {
        if self.tested >= self.req.policy.max_hypotheses {
            (
                "hypothesis limit",
                format!(
                    "hypothesis limit ({}) reached during {stage}",
                    self.req.policy.max_hypotheses
                ),
            )
        } else {
            (
                "deadline",
                format!(
                    "deadline reached during {stage}, {:.2} s after the start, {} hypotheses tested",
                    self.started.elapsed().as_secs_f64(),
                    self.tested
                ),
            )
        }
    }

    /// The cone holding the frame under `wcs`: the sky position of the frame centre and its
    /// largest separation from the corners and edge midpoints. Exact for any CRPIX and for
    /// anisotropic or distorted WCSs, unlike a radius from the mean scale around CRVAL.
    fn footprint(&self, wcs: &Wcs) -> (SkyPos, f64) {
        let (x0, y0, x1, y1) = (-0.5, -0.5, self.w - 0.5, self.h - 0.5);
        let (xm, ym) = (0.5 * (x0 + x1), 0.5 * (y0 + y1));
        let centre = wcs.pixel_to_sky(PixelPos::new(xm, ym));
        let r = [
            (x0, y0),
            (x1, y0),
            (x0, y1),
            (x1, y1),
            (xm, y0),
            (xm, y1),
            (x0, ym),
            (x1, ym),
        ]
        .into_iter()
        .map(|(x, y)| centre.separation(wcs.pixel_to_sky(PixelPos::new(x, y))))
        .fold(0.0, f64::max);
        (centre, r)
    }

    /// The brightest `n_ref` catalogue stars inside the frame under `wcs`.
    fn refs(&self, set: &RefSet, wcs: &Wcs) -> Vec<RefStar> {
        let scale = wcs.pixel_scale_arcsec();
        let (centre, r) = self.footprint(wcs);
        let mut idx = set.tree.within(centre.to_unit_vector().0, chord(r * 1.05));
        idx.sort_unstable(); // RefSet is in magnitude order
        let mut out = Vec::new();
        for k in idx {
            if let Some(px) = wcs.sky_to_pixel(set.stars[k].pos) {
                if px.x > -0.5 && px.y > -0.5 && px.x < self.w - 0.5 && px.y < self.h - 0.5 {
                    let e = set.stars[k].pos_err_mas;
                    // Per-axis 1σ: RMS of the two axis errors, mas → px.
                    let sigma_mas =
                        (0.5 * (f64::from(e[0]).ppowi(2) + f64::from(e[1]).ppowi(2))).sqrt();
                    out.push(RefStar {
                        px: [px.x, px.y],
                        star: k,
                        sigma_px: sigma_mas / 1000.0 / scale,
                    });
                    if out.len() == self.n_ref {
                        break;
                    }
                }
            }
        }
        out
    }

    /// Whether the ranked sources `matched` cover the source list well enough to constrain a
    /// distortion polynomial across the frame ([`SIP_CELLS`]).
    fn covers_sources(&self, matched: impl Iterator<Item = usize>) -> bool {
        let cell = |px: [f64; 2]| {
            let at = |v: f64, size: f64| {
                (((v + 0.5) / size * SIP_CELLS as f64) as usize).min(SIP_CELLS - 1)
            };
            at(px[1], self.h) * SIP_CELLS + at(px[0], self.w)
        };
        let mut sources = [0usize; SIP_CELLS * SIP_CELLS];
        let mut hits = [0usize; SIP_CELLS * SIP_CELLS];
        for r in &self.ranked {
            sources[cell(r.px)] += 1;
        }
        let mut total = 0usize;
        for i in matched {
            hits[cell(self.ranked[i].px)] += 1;
            total += 1;
        }
        let overall = total as f64 / self.ranked.len().max(1) as f64;
        (0..SIP_CELLS * SIP_CELLS).all(|c| {
            sources[c] < SIP_CELL_MIN_SOURCES
                || hits[c] as f64 >= SIP_CELL_MIN_FRACTION * overall * sources[c] as f64
        })
    }

    /// The scale and roll priors as terms of the final fit ([`LinearPrior`]), with the weight
    /// convention of [`Ctx::converge`]'s pairs at `error_scale`.
    fn linear_prior(&self, error_scale: f64) -> LinearPrior {
        let pr = &self.req.priors;
        LinearPrior {
            scale: pr.scale.filter(|s| s.sigma > 0.0).map(|s| {
                (
                    (s.arcsec_per_px / 3600.0).to_radians(),
                    (s.sigma / 3600.0).to_radians(),
                )
            }),
            roll: pr
                .roll
                .filter(|r| r.sigma > 0.0)
                .map(|r| (r.angle, r.sigma)),
            parity: pr.parity,
            weight_scale2: error_scale * error_scale,
        }
    }

    /// Hard limits from the priors.
    fn within_priors(&self, wcs: &TanWcs) -> bool {
        let s = wcs.pixel_scale_arcsec();
        if s < self.s_lo || s > self.s_hi {
            return false;
        }
        let pr = &self.req.priors;
        if !self.outside_pointing
            && pr
                .pointing
                .is_some_and(|pt| wcs.crval().separation(pt.centre) > pt.max_radius)
        {
            return false;
        }
        if pr.parity.is_some_and(|par| wcs.parity() != par) {
            return false;
        }
        if let Some(roll) = pr.roll {
            let cd = wcs.cd();
            let r = cd[0][1].patan2(cd[1][1]);
            let d = (r - roll.angle + PI).rem_euclid(2.0 * PI) - PI;
            if d.abs() > roll.max_dev {
                return false;
            }
        }
        true
    }

    /// At least `need` of the next-ranked sources (excluding `used`) land on stars of `tree`,
    /// within a tolerance that grows with distance from `anchor` (seed transform error).
    fn prefilter(
        &self,
        wcs: &TanWcs,
        used: &[usize],
        anchor: [f64; 2],
        tree: &KdTree<3>,
        need: usize,
    ) -> bool {
        let scale = (wcs.pixel_scale_arcsec() / 3600.0).to_radians();
        let mut hits = 0;
        for (k, r) in self
            .ranked
            .iter()
            .enumerate()
            .take(PREFILTER_SOURCES + used.len())
        {
            if used.contains(&k) {
                continue;
            }
            let dist = (r.px[0] - anchor[0]).phypot(r.px[1] - anchor[1]);
            let tol_px = 3.0 * r.sigma + 2.0 + 0.01 * dist;
            let sky = wcs.pixel_to_sky(PixelPos::new(r.px[0], r.px[1]));
            if !tree
                .within(sky.to_unit_vector().0, chord(tol_px * scale))
                .is_empty()
            {
                hits += 1;
                if hits >= need {
                    return true;
                }
            }
        }
        false
    }

    /// Greedy unique matching of ranked sources to references within `radius(source)` pixels.
    fn match_sources(
        &self,
        refs: &[RefStar],
        radius: impl Fn(&Ranked) -> f64,
    ) -> Vec<(usize, usize)> {
        let tree = KdTree::build(refs.iter().map(|r| r.px).collect());
        let mut used = vec![false; refs.len()];
        let mut out = Vec::new();
        for (i, r) in self.ranked.iter().enumerate() {
            let mut c = tree.within(r.px, radius(r));
            let d2 =
                |k: usize| (refs[k].px[0] - r.px[0]).ppowi(2) + (refs[k].px[1] - r.px[1]).ppowi(2);
            c.sort_by(|&x, &y| d2(x).total_cmp(&d2(y)).then(x.cmp(&y)));
            if let Some(k) = c.into_iter().find(|&k| !used[k]) {
                used[k] = true;
                out.push((i, k));
            }
        }
        out
    }

    /// Match ranked sources to the references within `radius_px` under `wcs` and refit: a
    /// similarity for a handful of matches, otherwise a polynomial of order ≤ `max_order` chosen
    /// by BIC among the orders with `obs_per_param` observations per parameter (order 1 is an
    /// affine TAN). Returns the match count and the fit.
    fn rematch_fit(
        &self,
        set: &RefSet,
        wcs: &Wcs,
        radius_px: f64,
        max_order: usize,
        obs_per_param: f64,
    ) -> Option<(usize, Wcs)> {
        let refs = self.refs(set, wcs);
        // Truncated sources gate matches but never enter a fit.
        let m: Vec<(usize, usize)> = self
            .match_sources(&refs, |r| radius_px.max(3.0 * r.sigma))
            .into_iter()
            .filter(|&(i, _)| !self.ranked[i].truncated)
            .collect();
        if m.len() < 4 {
            return None;
        }
        if m.len() < AFFINE_MIN_MATCHES {
            let pairs: Vec<([f64; 2], SkyPos)> = m
                .iter()
                .map(|&(i, k)| (self.ranked[i].px, set.stars[refs[k].star].pos))
                .collect();
            let (fitted, _) = similarity_seed(&pairs, wcs.parity(), self.crpix)?;
            return Some((m.len(), fitted.into()));
        }
        let pairs: Vec<Pair> = m
            .iter()
            .map(|&(i, k)| {
                let px = self.ranked[i].px;
                let variance = self.ranked[i].sigma.ppowi(2) + refs[k].sigma_px.ppowi(2);
                let r2 = (px[0] - refs[k].px[0]).ppowi(2) + (px[1] - refs[k].px[1]).ppowi(2);
                // Cauchy IRLS (Holland & Welsch 1977), with Fayde's r0 = 3σ_tot.
                Pair {
                    px,
                    sky: set.stars[refs[k].star].pos,
                    w: 1.0 / (variance * (1.0 + r2 / (9.0 * variance))),
                }
            })
            .collect();
        let (fit, _) = fit_best_order(
            wcs.crval(),
            self.crpix,
            &pairs,
            m.len() as f64,
            max_order,
            obs_per_param,
            None,
        )?;
        Some((m.len(), fit.wcs))
    }

    /// Re-match with shrinking radii and refit a linear WCS (LO-RANSAC).
    fn local_refit(&self, set: &RefSet, seed: TanWcs) -> Option<TanWcs> {
        let mut wcs = Wcs::from(seed);
        for radius_px in [12.0f64, 6.0, 3.0] {
            wcs = self
                .rematch_fit(set, &wcs, radius_px, 1, FINAL_OBS_PER_PARAM)?
                .1;
        }
        Some(wcs.tan)
    }

    /// Distortion tune-up of a converged linear hypothesis, for frames whose lens distortion
    /// leaves the outer stars outside the match gate of any linear model. Fayde's own;
    /// astrometry.net also refits a promising match before deciding (its "tune-up"), without
    /// distortion at that point.
    ///
    /// At the wide radius the order is chosen by BIC and the stars are re-matched until no more
    /// join (a distortion fit brings outer stars into range); then the radius shrinks as in the
    /// linear refit. `None` when there are too few matches to fit distortion or BIC keeps the
    /// linear model.
    fn tune(&self, set: &RefSet, linear: &TanWcs) -> Option<Wcs> {
        let max = self.req.policy.max_sip_order;
        let mut wcs = Wcs::from(*linear);
        let mut prev = 0usize;
        for _ in 0..TUNE_WIDE_PASSES {
            let (n, fitted) =
                self.rematch_fit(set, &wcs, TUNE_RADIUS_PX, max, GROWTH_OBS_PER_PARAM)?;
            if n < TUNE_MIN_MATCHES {
                return None;
            }
            wcs = fitted;
            if n <= prev {
                break;
            }
            prev = n;
        }
        let order = wcs.sip.as_ref()?.a.order;
        for radius_px in [6.0f64, 3.0] {
            wcs = self
                .rematch_fit(set, &wcs, radius_px, order, GROWTH_OBS_PER_PARAM)?
                .1;
        }
        Some(wcs)
    }

    /// Log of the prior density of `wcs` relative to a uniform density over the hard limits,
    /// for each prior given (design doc §6.1: σ ranks hypotheses, the limit bounds the search).
    /// Small-angle forms: pointing is a 2-D Gaussian against a uniform cap of the hard radius,
    /// scale and roll are 1-D Gaussians against uniform ranges.
    fn prior_log_ratio(&self, wcs: &Wcs) -> f64 {
        let pr = &self.req.priors;
        let mut lr = 0.0;
        if let (Some(pt), false) = (pr.pointing, self.outside_pointing) {
            if pt.sigma > 0.0 && pt.max_radius > 0.0 {
                let d = wcs.crval().separation(pt.centre);
                lr += (pt.max_radius.ppowi(2) / (2.0 * pt.sigma.ppowi(2))).pln()
                    - d * d / (2.0 * pt.sigma.ppowi(2));
            }
        }
        if let Some(sc) = pr.scale {
            if sc.sigma > 0.0 && self.s_hi.is_finite() {
                let width = self.s_hi - self.s_lo;
                let z = (wcs.pixel_scale_arcsec() - sc.arcsec_per_px) / sc.sigma;
                lr += (width / ((2.0 * PI).sqrt() * sc.sigma)).pln() - 0.5 * z * z;
            }
        }
        if let Some(roll) = pr.roll {
            if roll.sigma > 0.0 && roll.max_dev > 0.0 {
                let cd = wcs.cd();
                let r = cd[0][1].patan2(cd[1][1]);
                let d = (r - roll.angle + PI).rem_euclid(2.0 * PI) - PI;
                lr += (2.0 * roll.max_dev / ((2.0 * PI).sqrt() * roll.sigma)).pln()
                    - d * d / (2.0 * roll.sigma.ppowi(2));
            }
        }
        lr
    }

    /// Refit and verify an already-counted hypothesis; `Some` once accepted. `seed_sources`
    /// are the ranked sources the hypothesis was built from: they are not scored.
    fn try_hypothesis(
        &mut self,
        set: &RefSet,
        seed: TanWcs,
        seed_sources: &[usize],
    ) -> Option<Box<Solution>> {
        let linear = self.local_refit(set, seed)?;
        if !self.within_priors(&linear) {
            self.limit_fits.set(self.limit_fits.get() + 1);
            return None;
        }
        let area = self.w * self.h;
        let accept = self.req.policy.accept_log_odds;
        let refined = Wcs::from(linear);
        let refs = self.refs(set, &refined);
        let gated = assess(&self.tests, &refs, area, 1.0, seed_sources);
        // Posterior odds: data log-odds plus the prior terms.
        let mut score = gated.log_odds + self.prior_log_ratio(&refined);
        let mut distorted = None;
        // Distortion is suspected when the linear model is not conclusive, or when it leaves
        // many stars just outside its match gates.
        let suspect = self.req.policy.max_sip_order >= 2
            && (!passes(score, self.tested, accept) || {
                let wide = self.match_sources(&refs, |_| TUNE_RADIUS_PX).len();
                wide >= TUNE_MIN_MATCHES
                    && wide as f64 > TUNE_SUSPECT_FACTOR * gated.hits.len() as f64 + 3.0
            });
        if suspect {
            // Try the distortion tune-up, as one more hypothesis, and keep whichever model
            // explains the whole source list better.
            if let Some(tuned) = self
                .tune(set, &linear)
                .filter(|t| self.within_priors(&t.tan))
            {
                self.tested += 1;
                let tuned_score = assess(
                    &self.tests,
                    &self.refs(set, &tuned),
                    area,
                    1.0,
                    seed_sources,
                )
                .log_odds
                    + self.prior_log_ratio(&tuned);
                if tuned_score > score {
                    (score, distorted) = (tuned_score, Some(tuned));
                }
            }
        }
        self.best = Some(self.best.map_or(score, |b| b.max(score)));
        if !passes(score, self.tested, accept) {
            return None;
        }
        self.final_fit(set, distorted.unwrap_or(refined), seed_sources)
            .map(Box::new)
    }

    /// IRLS from `start` until the match set and the error scale repeat. A start that carries
    /// distortion (from the tune-up) keeps it from the first iteration; a linear start settles
    /// its match set with TAN first and then lets BIC choose the order.
    fn converge(&self, set: &RefSet, start: Wcs, seed_sources: &[usize]) -> Option<Converged> {
        let area = self.w * self.h;
        let tuned = start.sip.is_some();
        let mut cap = self.req.policy.max_sip_order;
        loop {
            let mut wcs = start.clone();
            let mut error_scale = 1.0f64;
            let mut covariance = None;
            let mut last: Option<(Vec<(usize, usize)>, f64)> = None;
            for iteration in 0..FINAL_FIT_MAX_ITERATIONS {
                let refs = self.refs(set, &wcs);
                let a = assess(&self.tests, &refs, area, error_scale, seed_sources);
                let matched: Vec<(usize, usize)> = a
                    .hits
                    .iter()
                    .map(|h| (h.source, refs[h.reference].star))
                    .collect();
                if iteration >= 4
                    && last
                        .as_ref()
                        .is_some_and(|(m, s)| *m == matched && (s - error_scale).abs() < 1e-3)
                {
                    break;
                }
                last = Some((matched, error_scale));
                let pairs: Vec<Pair> = a
                    .hits
                    .iter()
                    .map(|h| Pair {
                        px: self.ranked[h.source].px,
                        sky: set.stars[refs[h.reference].star].pos,
                        // Preserve the effective relative weights, including each source's
                        // sigma floor and the unscaled catalogue error. Factor out the common
                        // scale for BIC/Birge estimation; covariance restores it below.
                        w: if self.ranked[h.source].truncated {
                            0.0
                        } else {
                            h.inlier_prob / h.sigma.ppowi(2) * error_scale.ppowi(2)
                        },
                    })
                    .collect();
                // Observations that constrain the fit: truncated sources carry zero weight, so
                // they count for neither the BIC sample size nor the error-scale dof.
                let inlier_sum: f64 = a
                    .hits
                    .iter()
                    .filter(|h| !self.ranked[h.source].truncated)
                    .map(|h| h.inlier_prob)
                    .sum();
                let max_order = if iteration < 2 && !tuned { 1 } else { cap };
                let prior = self.linear_prior(error_scale);
                let (fit, chi2) = fit_best_order(
                    wcs.crval(),
                    self.crpix,
                    &pairs,
                    inlier_sum,
                    max_order,
                    FINAL_OBS_PER_PARAM,
                    Some(&prior),
                )?;
                // Birge ratio: default centroid errors may shrink. Producer-supplied
                // errors retain their floor in assess, including mixed input lists.
                let floor = if a.hits.iter().all(|h| !self.tests[h.source].sigma_defaulted) {
                    1.0
                } else {
                    f64::EPSILON
                };
                let k = fit.cov.terms.len() as f64;
                error_scale = (chi2 / (2.0 * inlier_sum - 2.0 * k).max(1.0))
                    .sqrt()
                    .max(floor);
                wcs = fit.wcs;
                let mut cov = fit.cov;
                cov.scale2 = error_scale * error_scale;
                covariance = Some(cov);
            }
            let refs = self.refs(set, &wcs);
            let assessment = assess(&self.tests, &refs, area, error_scale, seed_sources);
            // Distortion needs support across the sources. Without it a tuned start has no
            // solution, and a linear start is refitted as TAN.
            let confident = assessment
                .hits
                .iter()
                .filter(|h| h.inlier_prob > 0.5)
                .map(|h| h.source);
            if wcs.sip.is_some() && !self.covers_sources(confident) {
                if tuned {
                    return None;
                }
                cap = 1;
                continue;
            }
            let log_odds = assessment.log_odds + self.prior_log_ratio(&wcs);
            return Some(Converged {
                wcs,
                error_scale,
                covariance: covariance?,
                refs,
                assessment,
                log_odds,
            });
        }
    }

    /// Final fit of an accepted hypothesis, from its refined model (linear, or the tune-up's).
    fn final_fit(&self, set: &RefSet, start: Wcs, seed_sources: &[usize]) -> Option<Solution> {
        let accept = self.req.policy.accept_log_odds;
        let best = self.converge(set, start, seed_sources);
        let Converged {
            mut wcs,
            error_scale,
            covariance,
            refs,
            assessment: a,
            log_odds,
        } = best?;
        // The WCS that is returned must itself pass, not only the hypothesis it grew from, and
        // must lie within the hard limits: the final fit can move it.
        if !passes(log_odds, self.tested, accept) {
            return None;
        }
        if !self.within_priors(&wcs.tan) {
            self.limit_fits.set(self.limit_fits.get() + 1);
            return None;
        }
        if wcs.sip.is_some() {
            let order = wcs.sip_order() + 2;
            wcs.fit_inverse(self.w, self.h, order);
        }
        let mut matches = Vec::new();
        let (mut ss, mut wsum) = (0.0f64, 0.0f64);
        for h in &a.hits {
            let r = &self.ranked[h.source];
            let star = &set.stars[refs[h.reference].star];
            let px = wcs.sky_to_pixel(star.pos)?;
            let res = [r.px[0] - px.x, r.px[1] - px.y];
            ss += h.inlier_prob * (res[0] * res[0] + res[1] * res[1]);
            wsum += h.inlier_prob;
            matches.push(Match {
                source: r.input,
                star_id: star.id,
                star_pos: star.pos,
                star_mag: star.mag,
                residual_px: res,
                inlier_prob: h.inlier_prob,
            });
        }
        matches.sort_by_key(|m| m.source);
        // Leave-one-out residuals of the confident, fitted matches: r / (1 − h), with leverage
        // h = w · fᵀ(AᵀWA)⁻¹f from the final fit (standard least-squares identity).
        let mut loo = Vec::new();
        for h in a.hits.iter().filter(|h| h.inlier_prob > 0.5) {
            let r = &self.ranked[h.source];
            if r.truncated {
                continue;
            }
            let star = &set.stars[refs[h.reference].star];
            let Some(px) = wcs.sky_to_pixel(star.pos) else {
                continue;
            };
            let w_fit = h.inlier_prob / h.sigma.ppowi(2) * error_scale.ppowi(2);
            let lev = w_fit * covariance.basis_quadratic(PixelPos::new(r.px[0], r.px[1]));
            if lev < 0.999 {
                let res2 = (r.px[0] - px.x).ppowi(2) + (r.px[1] - px.y).ppowi(2);
                loo.push(res2 / (1.0 - lev).ppowi(2));
            }
        }
        let loo_rms_px = (loo.len() > covariance.terms.len() + 1)
            .then(|| (loo.iter().sum::<f64>() / loo.len() as f64).sqrt());
        let worst_sigma_arcsec = covariance.worst_radial_sigma_arcsec(self.w, self.h);
        let violated_priors = violated_priors(&self.req.priors, &wcs.tan);
        Some(Solution {
            worst_sigma_arcsec,
            loo_rms_px,
            method: self.method,
            search: self.trace.clone(),
            prior_violated: !violated_priors.is_empty(),
            violated_priors,
            index_hash: self.index_hash,
            wcs,
            covariance,
            log_odds,
            nfa: self.tested as f64 * (-log_odds).pexp(),
            hypotheses_tested: self.tested,
            matches,
            rms_px: (ss / wsum.max(1e-12)).sqrt(),
            error_scale,
            catalog: self.catalog.id(),
            elapsed: self.started.elapsed(),
        })
    }

    fn abstain(&self, reason: &str) -> Outcome {
        let (seeds, fits) = (self.limit_seeds.get(), self.limit_fits.get());
        let reason = if seeds + fits > 0 {
            format!(
                "{reason}; outside the prior hard limits: {fits} refined fits rejected, \
                 {seeds} seeds skipped"
            )
        } else {
            reason.to_string()
        };
        Outcome::Abstained {
            reason,
            best_log_odds: self.best,
            hypotheses_tested: self.tested,
            search: self.trace.clone(),
        }
    }

    /// The catalogue cone covering the whole pointing prior (pointing + scale priors needed).
    fn prior_cone(&mut self) -> Result<Option<std::rc::Rc<RefSet>>, SolveError> {
        let (Some(pt), true) = (self.req.priors.pointing, self.s_hi.is_finite()) else {
            return Ok(None);
        };
        if self.prior_cone.is_none() {
            let r = pt.max_radius + self.half_diag_px() * (self.s_hi / 3600.0).to_radians();
            self.prior_cone = Some(std::rc::Rc::new(RefSet::from_cone(
                self.catalog,
                pt.centre,
                r,
                &self.query,
            )?));
        }
        Ok(self.prior_cone.clone())
    }
}

fn new_ctx<'a>(req: &'a SolveRequest, catalog: &'a dyn Catalog) -> Result<Ctx<'a>, SolveError> {
    // The deadline counts from here: source ranking and duplicate detection are part of it.
    let started = Instant::now();
    let p = &req.policy;
    if req.frame.width == 0 || req.frame.height == 0 {
        return Err(SolveError::InvalidInput("frame has zero size".into()));
    }
    if !(0.0..1.0).contains(&p.distractor) || p.default_sigma_px <= 0.0 {
        return Err(SolveError::InvalidInput(
            "distractor must be in [0,1), sigma > 0".into(),
        ));
    }
    let (w, h) = (f64::from(req.frame.width), f64::from(req.frame.height));
    let (s_lo, s_hi) = match req.priors.scale {
        Some(sc) => {
            if !(sc.arcsec_per_px > 0.0 && sc.sigma >= 0.0) {
                return Err(SolveError::InvalidInput(
                    "scale prior must be positive".into(),
                ));
            }
            // A frame wider than a hemisphere cannot be imaged through a gnomonic projection:
            // such a "prior" is a placeholder value (e.g. CDELT = 1), not information.
            let half_diag = 0.5 * w.phypot(h) * (sc.arcsec_per_px / 3600.0).to_radians();
            if half_diag >= std::f64::consts::FRAC_PI_2 {
                return Err(SolveError::InvalidInput(format!(
                    "scale prior {} arcsec/px implies a field wider than a hemisphere",
                    sc.arcsec_per_px
                )));
            }
            scale_limits(sc)
        }
        None => (0.0, f64::INFINITY),
    };
    validate(req)?;

    // Ranking key: (group, score). Sources with a measured flux (group 1) rank by p_star × flux;
    // sources without one (group 0) follow, in input order weighted by p_star.
    let n = req.sources.len();
    let mut ranked: Vec<((u8, f64), f64, Ranked)> = req
        .sources
        .iter()
        .enumerate()
        .filter(|(_, s)| {
            s.x.is_finite()
                && s.y.is_finite()
                && s.x > -0.5
                && s.y > -0.5
                && s.x < w - 0.5
                && s.y < h - 0.5
        })
        .map(|(i, s)| {
            let p_star = s.p_star.unwrap_or(1.0).clamp(0.0, 1.0);
            let key = match s.flux.filter(|f| f.is_finite()) {
                Some(f) => (1u8, p_star * f.max(0.0)),
                None => (0u8, p_star * (n - i) as f64),
            };
            let sigma = s
                .sigma_px
                .filter(|v| *v > 0.0)
                .unwrap_or(p.default_sigma_px);
            (
                key,
                p_star,
                Ranked {
                    input: i,
                    px: [s.x, s.y],
                    sigma,
                    distractor: 0.0,
                    truncated: s.truncated,
                },
            )
        })
        .collect();
    let order = |a: &((u8, f64), f64, Ranked), b: &((u8, f64), f64, Ranked)| {
        b.0 .0
            .cmp(&a.0 .0)
            .then(b.0 .1.total_cmp(&a.0 .1))
            .then(a.2.input.cmp(&b.2.input))
    };
    ranked.sort_by(order);
    // Duplicates (design doc §5.0): a source within DUPLICATE_PX of a better-ranked one is most
    // likely a fragment of it. Its p_star drops by 10× (neutral in verification) and it moves
    // behind the others so it cannot seed pairs or patterns.
    // A grid of DUPLICATE_PX cells holds the better-ranked sources seen so far: a duplicate can
    // only be in the 3×3 cells around a source, so this is linear in the number of sources
    // (the all-pairs check took seconds for tens of thousands).
    let mut demoted = Vec::new();
    let mut grid: HashMap<(i64, i64), Vec<[f64; 2]>> = HashMap::new();
    let cell = |v: f64| (v / DUPLICATE_PX).floor() as i64;
    for (i, r) in ranked.iter_mut().enumerate() {
        let (x, y) = (r.2.px[0], r.2.px[1]);
        let (cx, cy) = (cell(x), cell(y));
        let dup = (cx - 1..=cx + 1)
            .flat_map(|gx| (cy - 1..=cy + 1).map(move |gy| (gx, gy)))
            .filter_map(|k| grid.get(&k))
            .flatten()
            .any(|p| (p[0] - x).phypot(p[1] - y) < DUPLICATE_PX);
        grid.entry((cx, cy)).or_default().push([x, y]);
        if dup {
            r.1 *= 0.1;
            demoted.push(i);
        }
    }
    for &i in &demoted {
        // Below every non-duplicate (scores are ≥ 0); duplicates keep input order.
        ranked[i].0 = (0, -1.0);
    }
    // Truncated sources are real stars with biased centroids: after every whole source (so they
    // never seed a hypothesis), before duplicates.
    for r in ranked
        .iter_mut()
        .filter(|r| r.2.truncated && r.0 != (0, -1.0))
    {
        r.0 = (0, -0.5);
    }
    ranked.sort_by(order);
    for r in &mut ranked {
        r.2.distractor = (1.0 - (1.0 - p.distractor) * r.1).clamp(1e-6, 1.0 - 1e-6);
    }
    // Effective real-source count, including duplicate demotion (Fayde choice).
    let effective_sources: f64 = ranked.iter().map(|r| r.1).sum();
    let ranked: Vec<Ranked> = ranked.into_iter().map(|(_, _, r)| r).collect();
    let tests = ranked
        .iter()
        .map(|r| {
            // A truncated trail's midpoint is off by up to half the trail length along the
            // trail; without a known trail, widen its gate in every direction instead.
            let (sigma, along, angle) = match (r.truncated, req.frame.trail) {
                (false, _) => (r.sigma, 0.0, 0.0),
                (true, Some(t)) => (r.sigma, 0.5 * t.length_px, t.angle),
                (true, None) => (TRUNCATED_SIGMA_FACTOR * r.sigma, 0.0, 0.0),
            };
            TestSource {
                px: r.px,
                sigma,
                sigma_defaulted: req.sources[r.input].sigma_px.is_none(),
                distractor: r.distractor,
                along,
                angle,
            }
        })
        .collect();
    Ok(Ctx {
        req,
        catalog,
        query: StarQuery {
            epoch: p.epoch,
            ..StarQuery::default()
        },
        n_ref: ((1.5 * effective_sources).round() as usize).clamp(10, 2000),
        ranked,
        tests,
        w,
        h,
        crpix: PixelPos::new((w - 1.0) / 2.0, (h - 1.0) / 2.0),
        s_lo,
        s_hi,
        started,
        method: SolveMethod::PairSeeds,
        index_hash: None,
        // A deadline past what `Instant` can represent is no deadline (no overflow panic).
        deadline: p.deadline.and_then(|d| started.checked_add(d)),
        tested: 0,
        best: None,
        outside_pointing: false,
        prior_cone: None,
        limit_seeds: std::cell::Cell::new(0),
        limit_fits: std::cell::Cell::new(0),
        trace: Vec::new(),
    })
}

/// One stage of an adaptive search plan.
#[derive(Debug, Clone)]
enum Stage {
    /// Star-pair seeds in the pointing-prior cone.
    PairSeeds,
    /// Pattern hash in one index, optionally restricted to a sky region (centre, radius of the
    /// field centre), optionally ignoring the pointing limit.
    Hash {
        index: usize,
        region: Option<(SkyPos, f64)>,
        outside: bool,
    },
}

/// Plan the search from the metadata the request carries (design doc §5.3):
///
/// | Known                     | Stages                                                        |
/// |---------------------------|---------------------------------------------------------------|
/// | pointing + tight scale    | pair seeds → hash inner ring (2σ) → hash to the hard radius     |
/// | pointing (loose/no scale) | hash inner ring → hash to the hard radius                       |
/// | scale only                | hash, whole index, scale-pruned                                 |
/// | nothing                   | hash, whole index, every band (blind)                           |
/// | + `search_outside_pointing` | then hash outside the pointing limit (flagged)               |
///
/// Within a stage, indexes are ordered by how well their field size matches the frame (scale
/// prior) or largest first; indexes whose region cannot contain the field are skipped.
fn plan(ctx: &Ctx<'_>, indexes: &[Index]) -> Vec<Stage> {
    let pr = &ctx.req.priors;
    let mut stages = Vec::new();
    let tight = matches!((pr.pointing, pr.scale), (Some(_), Some(sc)) if sc.sigma <= TIGHT_SCALE * sc.arcsec_per_px);
    if tight {
        stages.push(Stage::PairSeeds);
    }
    let mut order: Vec<usize> = (0..indexes.len()).collect();
    match pr.scale {
        Some(sc) => {
            let diag_deg = 2.0 * ctx.half_diag_px() * sc.arcsec_per_px / 3600.0;
            let d = |i: &usize| (indexes[*i].meta().params.fov_deg / diag_deg).pln().abs();
            order.sort_by(|a, b| d(a).total_cmp(&d(b)).then(a.cmp(b)));
        }
        None => order.sort_by(|a, b| {
            indexes[*b]
                .meta()
                .params
                .fov_deg
                .total_cmp(&indexes[*a].meta().params.fov_deg)
                .then(a.cmp(b))
        }),
    }
    // Can index `i` hold a field centred within `radius` of `centre`?
    let covers = |i: usize, centre: SkyPos, radius: f64| match indexes[i].meta().params.region {
        None => true,
        Some([ra, dec, r]) => SkyPos::from_degrees(ra, dec)
            .is_some_and(|c| c.separation(centre) <= r.to_radians() + radius),
    };
    match pr.pointing {
        Some(pt) => {
            let inner = 2.0 * pt.sigma;
            let rings: Vec<f64> = if inner < 0.8 * pt.max_radius {
                vec![inner, pt.max_radius]
            } else {
                vec![pt.max_radius]
            };
            for r in rings {
                for &i in order.iter().filter(|&&i| covers(i, pt.centre, r)) {
                    stages.push(Stage::Hash {
                        index: i,
                        region: Some((pt.centre, r)),
                        outside: false,
                    });
                }
            }
            if ctx.req.policy.search_outside_pointing {
                for &i in &order {
                    stages.push(Stage::Hash {
                        index: i,
                        region: None,
                        outside: true,
                    });
                }
            }
        }
        None => {
            for &i in &order {
                stages.push(Stage::Hash {
                    index: i,
                    region: None,
                    outside: false,
                });
            }
        }
    }
    stages
}

/// Whether `log_odds` clears `accept` after `tested` hypotheses. NaN never passes.
fn passes(log_odds: f64, tested: u64, accept: f64) -> bool {
    log_odds - (tested.max(1) as f64).pln() >= accept
}

/// Reject non-finite or out-of-range numbers in the priors and policy: one NaN would otherwise
/// propagate into a score and through comparisons that are false for NaN.
fn validate(req: &SolveRequest) -> Result<(), SolveError> {
    let bad = |what: &str| Err(SolveError::InvalidInput(what.into()));
    let fin = |v: f64| v.is_finite();
    let pr = &req.priors;
    if let Some(pt) = pr.pointing {
        let (ra, dec) = pt.centre.to_degrees();
        let ok = fin(ra)
            && fin(dec)
            && fin(pt.sigma)
            && pt.sigma >= 0.0
            && fin(pt.max_radius)
            && pt.max_radius >= 0.0;
        if !ok {
            return bad("pointing prior needs a finite centre, sigma >= 0 and max_radius >= 0");
        }
    }
    if let Some(sc) = pr.scale {
        if !(fin(sc.arcsec_per_px) && fin(sc.sigma)) {
            return bad("scale prior must be finite");
        }
    }
    if let Some(r) = pr.roll {
        if !(fin(r.angle) && fin(r.sigma) && r.sigma >= 0.0 && fin(r.max_dev) && r.max_dev >= 0.0) {
            return bad("roll prior needs a finite angle, sigma >= 0 and max_dev >= 0");
        }
    }
    if let Some(t) = req.frame.trail {
        if !(fin(t.length_px) && t.length_px >= 0.0 && fin(t.angle)) {
            return bad("trail needs a finite length >= 0 and angle");
        }
    }
    let p = &req.policy;
    if !fin(p.accept_log_odds) || !fin(p.default_sigma_px) || !fin(p.distractor) {
        return bad("policy values must be finite");
    }
    if p.default_sigma_px <= 0.0 {
        return bad("default_sigma_px (pixel error) must be > 0");
    }
    if p.epoch.is_some_and(|e| !e.is_plausible()) {
        return bad("epoch must be within 1800 to 2200 (Julian year)");
    }
    if p.max_sip_order > 9 {
        return bad("max_sip_order must be at most 9, the limit of the SIP convention");
    }
    Ok(())
}

/// Hard scale limits of a scale prior: ±3σ, never below a tenth of the scale.
fn scale_limits(sc: ScalePrior) -> (f64, f64) {
    (
        (sc.arcsec_per_px - 3.0 * sc.sigma).max(sc.arcsec_per_px * 0.1),
        sc.arcsec_per_px + 3.0 * sc.sigma,
    )
}

/// The hard limits of `priors` that `wcs` breaks.
fn violated_priors(priors: &Priors, wcs: &TanWcs) -> Vec<&'static str> {
    let mut out = Vec::new();
    if priors
        .pointing
        .is_some_and(|pt| wcs.crval().separation(pt.centre) > pt.max_radius)
    {
        out.push("pointing");
    }
    if let Some(sc) = priors.scale {
        let (lo, hi) = scale_limits(sc);
        let s = wcs.pixel_scale_arcsec();
        if s < lo || s > hi {
            out.push("scale");
        }
    }
    if let Some(roll) = priors.roll {
        let cd = wcs.cd();
        let d = (cd[0][1].patan2(cd[1][1]) - roll.angle + PI).rem_euclid(2.0 * PI) - PI;
        if d.abs() > roll.max_dev {
            out.push("roll");
        }
    }
    if priors.parity.is_some_and(|p| wcs.parity() != p) {
        out.push("parity");
    }
    out
}

/// Least relative 1σ of the scale prior in the relaxed pass of
/// [`crate::Policy::relax_priors_on_failure`].
const RELAXED_SCALE_SIGMA: f64 = 0.01;

/// The request for the relaxed second pass, or `None` when relaxing changes nothing.
fn relaxed_request(req: &SolveRequest) -> Option<SolveRequest> {
    let mut r = req.clone();
    r.priors.roll = None;
    r.priors.parity = None;
    if let Some(sc) = r.priors.scale.as_mut() {
        sc.sigma = sc.sigma.max(RELAXED_SCALE_SIGMA * sc.arcsec_per_px);
    }
    r.policy.search_outside_pointing |= r.priors.pointing.is_some();
    r.policy.relax_priors_on_failure = false;
    (r.priors != req.priors
        || r.policy.search_outside_pointing != req.policy.search_outside_pointing)
        .then_some(r)
}

/// Solve adaptively: plan stages from the request's metadata and run them in order until one
/// produces a verified solution. Hypotheses are counted across all stages, so the false-alarm
/// accounting covers the whole search; the executed plan is returned in the outcome.
///
/// With [`crate::Policy::relax_priors_on_failure`], a failed search is retried once with the
/// hard limits relaxed.
pub fn solve_auto(
    req: &SolveRequest,
    catalog: &dyn Catalog,
    indexes: &[Index],
) -> Result<Outcome, SolveError> {
    let started = Instant::now();
    let first = solve_auto_once(req, catalog, indexes, 0)?;
    let Outcome::Abstained {
        reason: reason0,
        best_log_odds: best0,
        hypotheses_tested: tested0,
        search: search0,
    } = &first
    else {
        return Ok(first);
    };
    let Some(mut relaxed) = req
        .policy
        .relax_priors_on_failure
        .then(|| relaxed_request(req))
        .flatten()
    else {
        return Ok(first);
    };
    if let Some(d) = req.policy.deadline {
        let left = d.saturating_sub(started.elapsed());
        if left.is_zero() {
            return Ok(first);
        }
        relaxed.policy.deadline = Some(left);
    }
    if *tested0 >= req.policy.max_hypotheses {
        return Ok(first);
    }
    // The second pass continues the first pass's hypothesis count, so every candidate it tests
    // is judged against the false-alarm budget of both passes.
    let second = solve_auto_once(&relaxed, catalog, indexes, *tested0)?;
    let tag = |t: Vec<StageTrace>| {
        t.into_iter().map(|mut s| {
            s.stage = format!("relaxed priors: {}", s.stage);
            s
        })
    };
    let mut search = search0.clone();
    Ok(match second {
        Outcome::Solved(mut s) => {
            search.extend(tag(std::mem::take(&mut s.search)));
            s.search = search;
            s.violated_priors = violated_priors(&req.priors, &s.wcs.tan);
            s.prior_violated = !s.violated_priors.is_empty();
            Outcome::Solved(s)
        }
        Outcome::Abstained {
            reason,
            best_log_odds,
            hypotheses_tested,
            search: s2,
        } => {
            search.extend(tag(s2));
            Outcome::Abstained {
                reason: format!("{reason0}; relaxed priors: {reason}"),
                best_log_odds: match (*best0, best_log_odds) {
                    (Some(a), Some(b)) => Some(a.max(b)),
                    (a, b) => a.or(b),
                },
                hypotheses_tested,
                search,
            }
        }
    })
}

fn solve_auto_once(
    req: &SolveRequest,
    catalog: &dyn Catalog,
    indexes: &[Index],
    already_tested: u64,
) -> Result<Outcome, SolveError> {
    let mut ctx = new_ctx(req, catalog)?;
    ctx.tested = already_tested;
    let stages = plan(&ctx, indexes);
    if stages.is_empty() && !indexes.is_empty() {
        return Ok(ctx.abstain("no index covers the pointing prior"));
    }
    if stages.is_empty() {
        // A valid request this configuration cannot search (e.g. a scale prior dropped for
        // unknown binning, and no index): an answer, not an error.
        return Ok(ctx.abstain(
            "nothing to search: without both a pointing and a tight scale prior, a pattern index is needed",
        ));
    }
    if ctx.ranked.len() < 4 {
        return Ok(ctx.abstain("fewer than 4 usable sources"));
    }
    for stage in stages {
        let (t0, h0) = (Instant::now(), ctx.tested);
        let (label, result) = match stage {
            Stage::PairSeeds => {
                let r = ctx
                    .req
                    .priors
                    .pointing
                    .map_or(0.0, |p| p.max_radius)
                    .to_degrees();
                (format!("pair seeds r<{r:.2}°"), pair_seeds(&mut ctx)?)
            }
            Stage::Hash {
                index,
                region,
                outside,
            } => {
                ctx.outside_pointing = outside;
                let fov = indexes[index].meta().params.fov_deg;
                let label = match (region, outside) {
                    (Some((_, r)), _) => {
                        format!("pattern hash {fov:.2}° index r<{:.2}°", r.to_degrees())
                    }
                    (None, true) => format!("pattern hash {fov:.2}° index outside pointing limit"),
                    (None, false) => format!("pattern hash {fov:.2}° index all sky"),
                };
                (label, pattern_hash(&mut ctx, &indexes[index], region)?)
            }
        };
        let (outcome, budget_reason) = match &result {
            Some(Stop::Solved(_)) => ("solved", None),
            Some(Stop::Budget) => {
                let (o, r) = ctx.budget_stop(&label);
                (o, Some(r))
            }
            None => ("exhausted", None),
        };
        ctx.trace.push(StageTrace {
            stage: label,
            hypotheses: ctx.tested - h0,
            seconds: t0.elapsed().as_secs_f64(),
            outcome,
        });
        match result {
            Some(Stop::Solved(mut s)) => {
                s.search = ctx.trace.clone();
                return Ok(Outcome::Solved(s));
            }
            Some(Stop::Budget) => {
                return Ok(ctx.abstain(&budget_reason.expect("budget stop has a reason")))
            }
            None => ctx.outside_pointing = false,
        }
    }
    Ok(ctx.abstain("no hypothesis passed verification"))
}

/// Catalogue stars around the field of an existing WCS.
fn cone_around(ctx: &Ctx<'_>, wcs: &Wcs) -> Result<RefSet, SolveError> {
    let (centre, r) = ctx.footprint(wcs);
    RefSet::from_cone(ctx.catalog, centre, r * 1.2, &ctx.query)
}

/// Score an existing WCS (from any solver) against a source list with Fayde's verifier: the data
/// log-odds of the WCS as a single tested hypothesis, its matches and their residuals. Priors in
/// the request are ignored; `policy.epoch` sets the catalogue epoch.
pub fn verify_wcs(
    req: &SolveRequest,
    catalog: &dyn Catalog,
    wcs: &Wcs,
) -> Result<Verification, SolveError> {
    let ctx = new_ctx(req, catalog)?;
    let set = cone_around(&ctx, wcs)?;
    let refs = ctx.refs(&set, wcs);
    let a = assess(&ctx.tests, &refs, ctx.w * ctx.h, 1.0, &[]);
    let mut matches = Vec::new();
    let (mut ss, mut wsum) = (0.0f64, 0.0f64);
    for h in &a.hits {
        let r = &ctx.ranked[h.source];
        let star = &set.stars[refs[h.reference].star];
        let Some(px) = wcs.sky_to_pixel(star.pos) else {
            continue;
        };
        let res = [r.px[0] - px.x, r.px[1] - px.y];
        ss += h.inlier_prob * (res[0] * res[0] + res[1] * res[1]);
        wsum += h.inlier_prob;
        matches.push(Match {
            source: r.input,
            star_id: star.id,
            star_pos: star.pos,
            star_mag: star.mag,
            residual_px: res,
            inlier_prob: h.inlier_prob,
        });
    }
    matches.sort_by_key(|m| m.source);
    Ok(Verification {
        log_odds: a.log_odds,
        accepted: a.log_odds >= req.policy.accept_log_odds,
        matches,
        rms_px: (ss / wsum.max(1e-12)).sqrt(),
        references: refs.len(),
    })
}

/// Refit an existing WCS: re-match the source list around it with shrinking radii (12 → 6 →
/// 3 px, linear), then the same converged IRLS fit as a solve, with SIP order chosen by BIC up to
/// `policy.max_sip_order`, covariance and accuracy figures. Counts as one tested hypothesis; the
/// result must clear `policy.accept_log_odds` like any solution. For repairing stored solutions
/// (e.g. a per-frame SIP fitted to too few stars) without searching again.
pub fn refine_wcs(
    req: &SolveRequest,
    catalog: &dyn Catalog,
    wcs: &Wcs,
) -> Result<Outcome, SolveError> {
    let mut ctx = new_ctx(req, catalog)?;
    ctx.method = SolveMethod::Refit;
    ctx.tested = 1;
    if ctx.ranked.len() < 4 {
        return Ok(ctx.abstain("fewer than 4 usable sources"));
    }
    let set = cone_around(&ctx, wcs)?;
    let mut start = wcs.clone();
    for radius_px in [12.0f64, 6.0, 3.0] {
        if ctx.past_deadline() {
            return Ok(ctx.refit_deadline());
        }
        match ctx.rematch_fit(&set, &start, radius_px, 1, FINAL_OBS_PER_PARAM) {
            Some((_, fitted)) => start = fitted,
            None => return Ok(ctx.abstain("too few matches around the given WCS")),
        }
    }
    if ctx.past_deadline() {
        return Ok(ctx.refit_deadline());
    }
    let t0 = Instant::now();
    let solution = ctx.final_fit(&set, start, &[]);
    ctx.trace.push(StageTrace {
        stage: "refit around the given WCS".into(),
        hypotheses: 1,
        seconds: t0.elapsed().as_secs_f64(),
        outcome: if solution.is_some() {
            "solved"
        } else {
            "exhausted"
        },
    });
    Ok(match solution {
        Some(mut s) => {
            s.search = ctx.trace.clone();
            Outcome::Solved(Box::new(s))
        }
        None => ctx.abstain("the refitted WCS did not pass verification"),
    })
}

/// Pair seeds only (needs pointing and a tight scale prior); see [`solve_auto`].
pub fn solve(req: &SolveRequest, catalog: &dyn Catalog) -> Result<Outcome, SolveError> {
    let tight = matches!((req.priors.pointing, req.priors.scale), (Some(_), Some(sc)) if sc.sigma <= TIGHT_SCALE * sc.arcsec_per_px);
    if req.priors.pointing.is_none() || req.priors.scale.is_none() || !tight {
        return Err(SolveError::Unsupported(
            "pair seeds need a pointing and a tight (σ ≤ 3%) scale prior; use solve_auto with an index".into(),
        ));
    }
    solve_auto(req, catalog, &[])
}

/// [`solve_auto`] with one index.
pub fn solve_with_index(
    req: &SolveRequest,
    catalog: &dyn Catalog,
    index: &Index,
) -> Result<Outcome, SolveError> {
    solve_auto(req, catalog, std::slice::from_ref(index))
}

/// [`solve_auto`] with a set of indexes (scale bands).
pub fn solve_with_indexes(
    req: &SolveRequest,
    catalog: &dyn Catalog,
    indexes: &[Index],
) -> Result<Outcome, SolveError> {
    solve_auto(req, catalog, indexes)
}

/// Pair-seed hypotheses. `Ok(None)` when all seeds are exhausted without a solution.
fn pair_seeds(ctx: &mut Ctx<'_>) -> Result<Option<Stop>, SolveError> {
    let pointing = ctx.req.priors.pointing.expect("checked by caller");
    let scale = ctx.req.priors.scale.expect("checked by caller");
    let as2rad = |a: f64| (a / 3600.0).to_radians();
    let field_r_max = ctx.half_diag_px() * as2rad(ctx.s_hi);
    let cone_r = pointing.max_radius + field_r_max;
    let set = RefSet::from_cone(ctx.catalog, pointing.centre, cone_r, &ctx.query)?;
    if set.stars.len() < 4 {
        return Ok(None);
    }

    let field_area = ctx.w * ctx.h * as2rad(scale.arcsec_per_px).ppowi(2);
    let cone_area = 2.0 * PI * (1.0 - cone_r.pcos());
    let n_seed_src = ctx.req.policy.seed_sources.clamp(3, ctx.ranked.len());
    let n_seed_cat = ((SEED_STARS_PER_SOURCE * n_seed_src as f64 * cone_area / field_area).ceil()
        as usize)
        .clamp(n_seed_src.min(set.stars.len()), set.stars.len());
    let seed_tree = KdTree::build(set.vecs[..n_seed_cat].to_vec());
    let mut cat_pairs: Vec<(f64, u32, u32)> = Vec::new();
    for i in 0..n_seed_cat {
        for j in seed_tree.within(set.vecs[i], chord(2.0 * field_r_max)) {
            if j > i {
                cat_pairs.push((angle(&set.vecs[i], &set.vecs[j]), i as u32, j as u32));
            }
        }
    }
    cat_pairs.sort_by(|a, b| a.0.total_cmp(&b.0).then((a.1, a.2).cmp(&(b.1, b.2))));
    let parities: Vec<Parity> = match ctx.req.priors.parity {
        Some(par) => vec![par],
        None => vec![Parity::Negative, Parity::Positive],
    };

    for b in 1..n_seed_src {
        for a in 0..b {
            if ctx.ranked[a].truncated || ctx.ranked[b].truncated {
                continue;
            }
            let (pa, pb) = (ctx.ranked[a].px, ctx.ranked[b].px);
            let d_px = (pb[0] - pa[0]).phypot(pb[1] - pa[1]);
            if d_px < MIN_SEED_PAIR_PX {
                continue;
            }
            let tol = 3.0 * (ctx.ranked[a].sigma + ctx.ranked[b].sigma) * as2rad(ctx.s_hi);
            let lo = d_px * as2rad(ctx.s_lo) - tol;
            let hi = d_px * as2rad(ctx.s_hi) + tol;
            let start = cat_pairs.partition_point(|c| c.0 < lo);
            let end = cat_pairs.partition_point(|c| c.0 <= hi);
            let candidates = &cat_pairs[start..end];
            let mean = d_px * as2rad(scale.arcsec_per_px);
            for &(_, ci, cj) in pairs_from_mean(candidates, mean) {
                for (c1, c2) in [(ci, cj), (cj, ci)] {
                    for &parity in &parities {
                        if ctx.out_of_budget() {
                            return Ok(Some(Stop::Budget));
                        }
                        let Some(seed) = seed_wcs(
                            &set.stars,
                            c1 as usize,
                            c2 as usize,
                            pa,
                            pb,
                            parity,
                            ctx.crpix,
                        ) else {
                            continue;
                        };
                        if !ctx.within_priors(&seed) {
                            ctx.limit_seeds.set(ctx.limit_seeds.get() + 1);
                            continue;
                        }
                        ctx.tested += 1;
                        let anchor = [(pa[0] + pb[0]) / 2.0, (pa[1] + pb[1]) / 2.0];
                        if !ctx.prefilter(&seed, &[a, b], anchor, &seed_tree, 2) {
                            continue;
                        }
                        if let Some(sol) = ctx.try_hypothesis(&set, seed, &[a, b]) {
                            return Ok(Some(Stop::Solved(sol)));
                        }
                    }
                }
            }
        }
    }
    Ok(None)
}

/// Start at the prior mean and alternate outward through the sorted separation list.
/// Prior-guided ordering is Fayde's extension of PROSAC (Chum & Matas 2005).
fn pairs_from_mean(pairs: &[(f64, u32, u32)], mean: f64) -> impl Iterator<Item = &(f64, u32, u32)> {
    let mut right = pairs.partition_point(|p| p.0 < mean);
    let mut left = right;
    let mut take_left =
        left > 0 && (right == pairs.len() || mean - pairs[left - 1].0 <= pairs[right].0 - mean);
    std::iter::from_fn(move || {
        let index = if left > 0 && (take_left || right == pairs.len()) {
            left -= 1;
            take_left = false;
            left
        } else if right < pairs.len() {
            let i = right;
            right += 1;
            take_left = true;
            i
        } else {
            return None;
        };
        Some(&pairs[index])
    })
}

/// Pattern-hash hypotheses. `Ok(None)` when all patterns are exhausted without a solution.
fn pattern_hash(
    ctx: &mut Ctx<'_>,
    index: &Index,
    region: Option<(SkyPos, f64)>,
) -> Result<Option<Stop>, SolveError> {
    ctx.method = SolveMethod::PatternHash;
    ctx.index_hash = Some(index.content_hash());
    let istars = index.stars();
    let index_tree = index.star_tree();
    let index_fov = index.meta().params.fov_deg.to_radians();
    let as2rad = |a: f64| (a / 3600.0).to_radians();
    // Within the pointing prior one catalogue cone serves every hypothesis.
    let shared = if ctx.outside_pointing {
        None
    } else {
        ctx.prior_cone()?
    };
    let region = region.map(|(c, r)| (c.to_unit_vector().0, r));
    let m = ctx.req.policy.pattern_sources.max(4).min(ctx.ranked.len());
    // Without a pointing prior each hypothesis needs its own catalogue cone. Cache cones by
    // (HEALPix cell of the field centre, field-radius bucket): a cone around the cell centre,
    // padded by the cell size, serves every hypothesis whose centre falls in that cell.
    let mut cones: std::collections::HashMap<(u64, i32), std::rc::Rc<RefSet>> =
        std::collections::HashMap::new();
    for n in 4..=m {
        let last = n - 1;
        for i in 0..last {
            for j in i + 1..last {
                for k in j + 1..last {
                    let ids = [i, j, k, last];
                    if ids.iter().any(|&q| ctx.ranked[q].truncated) {
                        continue;
                    }
                    let pts = ids.map(|q| ctx.ranked[q].px);
                    let (shape, _) = image_shape(&pts);
                    if shape.ratios[0] * shape.largest < MIN_PATTERN_EDGE_PX {
                        continue;
                    }
                    // Ratio tolerance from centroid noise: δr ≈ 2√2·σ / L, ×3.
                    let sig = ids.iter().map(|&q| ctx.ranked[q].sigma).fold(0.0, f64::max);
                    let tol = (8.5 * sig / shape.largest).clamp(0.002, 0.02);
                    let range = ctx.s_hi.is_finite().then(|| {
                        (
                            shape.largest * as2rad(ctx.s_lo),
                            shape.largest * as2rad(ctx.s_hi),
                        )
                    });
                    if range.is_some_and(|(lo, _)| lo > index_fov) {
                        continue;
                    }
                    // Lookup ratios are scale-free; feed the index plane coordinates in pixels.
                    for found in index.lookup(&shape, tol, range) {
                        if ctx.out_of_budget() {
                            return Ok(Some(Stop::Budget));
                        }
                        // Implied scale (rad/px). Index patterns come from fields no larger than
                        // the index field size, so a hypothesis implying a frame more than twice
                        // that size cannot be this index's sky: skip it before any catalogue work.
                        let implied = found.largest / shape.largest;
                        if implied.is_nan()
                            || implied <= 0.0
                            || ctx.half_diag_px() * implied > index_fov
                        {
                            continue;
                        }
                        // Region: every pattern star must lie within the field of a frame
                        // centred inside the region.
                        if let Some((c, r)) = region {
                            let reach = r + ctx.half_diag_px() * implied;
                            if found
                                .stars
                                .iter()
                                .any(|&s| angle(&istars[s as usize].vec, &c) > reach)
                            {
                                continue;
                            }
                        }
                        let skies = found
                            .stars
                            .map(|s| UnitVector(istars[s as usize].vec).to_sky());
                        // Canonical centroid ranks can swap under noise. Test all edge-compatible
                        // assignments and both parities, best four-point residual first.
                        for (_, seed) in pattern_seeds(&pts, &skies, tol, sig, ctx.crpix) {
                            if ctx.out_of_budget() {
                                return Ok(Some(Stop::Budget));
                            }
                            if !ctx.within_priors(&seed) {
                                ctx.limit_seeds.set(ctx.limit_seeds.get() + 1);
                                continue;
                            }
                            ctx.tested += 1;
                            let used: Vec<usize> = ids.to_vec();
                            let anchor = [
                                pts.iter().map(|p| p[0]).sum::<f64>() / 4.0,
                                pts.iter().map(|p| p[1]).sum::<f64>() / 4.0,
                            ];
                            if !ctx.prefilter(&seed, &used, anchor, index_tree, 1) {
                                continue;
                            }
                            let local: std::rc::Rc<RefSet>;
                            let set = match &shared {
                                Some(s) => s,
                                None => {
                                    let field_r =
                                        ctx.half_diag_px() * as2rad(seed.pixel_scale_arcsec());
                                    // Radius bucket: next power of √2 above the field radius.
                                    let bucket = (field_r.pln() / std::f64::consts::SQRT_2.pln())
                                        .ceil()
                                        as i32;
                                    let r_bucket = std::f64::consts::SQRT_2.ppowi(bucket);
                                    let hp = fayde_index::healpix::Healpix::for_pixel_side(
                                        r_bucket / 4.0,
                                    );
                                    let cell = hp.pixel(seed.crval().to_unit_vector());
                                    local = match cones.get(&(cell, bucket)) {
                                        Some(c) => c.clone(),
                                        None => {
                                            let centre = hp.centre(cell).to_sky();
                                            // Any centre in the cell is within ~0.75 side of the cell centre.
                                            let r = r_bucket * 1.05 + hp.pixel_side() * 0.75;
                                            // A cone the catalogue rejects (e.g. larger than the
                                            // sky) means an implausible hypothesis, not a failure.
                                            let c = match RefSet::from_cone(
                                                ctx.catalog,
                                                centre,
                                                r,
                                                &ctx.query,
                                            ) {
                                                Ok(set) => std::rc::Rc::new(set),
                                                Err(SolveError::Catalog(
                                                    fayde_catalog::CatalogError::InvalidQuery(_),
                                                )) => continue,
                                                Err(e) => return Err(e),
                                            };
                                            cones.insert((cell, bucket), c.clone());
                                            c
                                        }
                                    };
                                    &local
                                }
                            };
                            if let Some(sol) = ctx.try_hypothesis(set, seed, &ids) {
                                return Ok(Some(Stop::Solved(sol)));
                            }
                        }
                    }
                }
            }
        }
    }
    Ok(None)
}

/// Similarity transform mapping image pixels `pa`, `pb` onto catalogue stars `c1`, `c2`.
fn seed_wcs(
    stars: &[Star],
    c1: usize,
    c2: usize,
    pa: [f64; 2],
    pb: [f64; 2],
    parity: Parity,
    crpix: PixelPos,
) -> Option<TanWcs> {
    similarity_seed(&[(pa, stars[c1].pos), (pb, stars[c2].pos)], parity, crpix)
        .map(|(seed, _)| seed)
}

/// Four-parameter similarity by closed-form complex least squares, the planar equivalent
/// of Umeyama (1991). Centre both sets, then α = Σ(z_c conj(z_p)) / Σ|z_p|².
/// Negative parity conjugates the pixels first.
///
/// The fit starts in the tangent plane of the catalogue centroid and is repeated with the
/// tangent point moved to the fitted `CRVAL` until it stops moving: a transform found in another
/// point's tangent plane is rotated by (east offset × tan δ) when reused at `CRVAL` (convergence
/// of the meridians), which is degrees near the poles.
/// Returns the seed and radial RMS in pixels in the final tangent plane.
fn similarity_seed(
    pairs: &[([f64; 2], SkyPos)],
    parity: Parity,
    crpix: PixelPos,
) -> Option<(TanWcs, f64)> {
    let mut tangent = UnitVector::normalize(pairs.iter().fold([0.0; 3], |mut sum, (_, sky)| {
        let v = sky.to_unit_vector().0;
        for k in 0..3 {
            sum[k] += v[k];
        }
        sum
    }))?
    .to_sky();
    let flip = if matches!(parity, Parity::Negative) {
        -1.0
    } else {
        1.0
    };
    let n = pairs.len() as f64;
    for pass in 0..SIMILARITY_TANGENT_PASSES {
        let projected: Vec<_> = pairs
            .iter()
            .map(|(p, sky)| {
                gnomonic(tangent, *sky)
                    .map(|(x, y)| ([p[0] - crpix.x, flip * (p[1] - crpix.y)], [x, y]))
            })
            .collect::<Option<_>>()?;
        let mut mp = [0.0; 2];
        let mut mc = [0.0; 2];
        for (p, c) in &projected {
            for k in 0..2 {
                mp[k] += p[k] / n;
                mc[k] += c[k] / n;
            }
        }
        let (mut den, mut ar, mut ai) = (0.0, 0.0, 0.0);
        for (p, c) in &projected {
            let (px, py) = (p[0] - mp[0], p[1] - mp[1]);
            let (cx, cy) = (c[0] - mc[0], c[1] - mc[1]);
            den += px * px + py * py;
            ar += cx * px + cy * py;
            ai += cy * px - cx * py;
        }
        if den <= 0.0 {
            return None;
        }
        ar /= den;
        ai /= den;
        let scale2 = ar * ar + ai * ai;
        if scale2 <= 0.0 {
            return None;
        }
        // Tangent-plane position of CRPIX (z_p = 0).
        let bx = mc[0] - (ar * mp[0] - ai * mp[1]);
        let by = mc[1] - (ai * mp[0] + ar * mp[1]);
        let crval = gnomonic_inverse(tangent, bx, by);
        if bx.phypot(by) > 1e-10 && pass + 1 < SIMILARITY_TANGENT_PASSES {
            tangent = crval;
            continue;
        }
        let ss: f64 = projected
            .iter()
            .map(|(p, c)| {
                (c[0] - (ar * p[0] - ai * p[1] + bx)).ppowi(2)
                    + (c[1] - (ai * p[0] + ar * p[1] + by)).ppowi(2)
            })
            .sum();
        let cd = [[ar, -flip * ai], [ai, flip * ar]];
        let wcs = TanWcs::new(crval, crpix, cd.map(|r| r.map(f64::to_degrees))).ok()?;
        return Some((wcs, (ss / (n * scale2)).sqrt()));
    }
    None
}

/// Enumerate at most 24 assignments. Compare all six labelled edges, not just their
/// sorted multiset, to avoid relying on unstable canonical ranks (Fayde choice).
fn pattern_seeds(
    pts: &[[f64; 2]; 4],
    skies: &[SkyPos; 4],
    tolerance: f64,
    sigma: f64,
    crpix: PixelPos,
) -> Vec<(f64, TanWcs)> {
    let mut image_edges = [[0.0f64; 4]; 4];
    let mut cat_edges = [[0.0f64; 4]; 4];
    let (mut image_max, mut cat_max) = (0.0f64, 0.0f64);
    for i in 0..4 {
        for j in i + 1..4 {
            let dp = (pts[i][0] - pts[j][0]).phypot(pts[i][1] - pts[j][1]);
            let dc = skies[i].separation(skies[j]);
            image_edges[i][j] = dp;
            image_edges[j][i] = dp;
            cat_edges[i][j] = dc;
            cat_edges[j][i] = dc;
            image_max = image_max.max(dp);
            cat_max = cat_max.max(dc);
        }
    }
    if image_max <= 0.0 || cat_max <= 0.0 {
        return Vec::new();
    }
    let mut seeds = Vec::new();
    for a in 0..4 {
        for b in 0..4 {
            for c in 0..4 {
                if a == b || a == c || b == c {
                    continue;
                }
                let order = [a, b, c, 6 - a - b - c];
                if (0..4).any(|i| {
                    (i + 1..4).any(|j| {
                        (image_edges[order[i]][order[j]] / image_max - cat_edges[i][j] / cat_max)
                            .abs()
                            > tolerance
                    })
                }) {
                    continue;
                }
                let pairs = std::array::from_fn::<_, 4, _>(|i| (pts[order[i]], skies[i]));
                for parity in [Parity::Negative, Parity::Positive] {
                    if let Some((seed, rms)) = similarity_seed(&pairs, parity, crpix) {
                        if rms <= 3.0 * sigma {
                            seeds.push((rms, seed));
                        }
                    }
                }
            }
        }
    }
    seeds.sort_by(|a, b| a.0.total_cmp(&b.0));
    seeds
}

#[cfg(test)]
mod tests;
