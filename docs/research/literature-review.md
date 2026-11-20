# Fayde plate solver: literature review

Companion to `docs/design/plate-solver.md` (the "design doc"; section
numbers like §5.3 refer to it). The baseline is taken as given: a tetra3/cedar-solve pattern hash for
hinted solving, astrometry.net quads for blind solving, astrometry.net log-odds verification, a
weighted TAN+SIP fit, and UCAC4.

**Scope and method.** Web search, then reading of papers, abstracts and project pages. Each entry
says how much was read:

- **[full]**: the full text was read, or the relevant sections of it.
- **[abstract]**: the abstract or publisher summary only.
- **[secondary]**: known only through citing papers, surveys or indexes.

Most IEEE/AIAA papers are paywalled and were assessed at abstract level. Verdicts assume the regime
Fayde is built for: about 1° FOV, about 3″/px, 30–60 sources, many false detections, stars often
trailed, and usually a pointing prior.

Already covered by the design doc and not repeated here: Lang et al. 2010 (astrometry.net:
geometric hashing with quads plus Bayesian verification), tetra3 (ESA's rewrite of Tetra, Brown & Stubis
2017: 4-star edge-ratio hash), cedar-solve (a uniform-database fork of tetra3).

---

## Summary of verdicts

| # | Item | Verdict | Design section affected |
|---|---|---|---|
| 1.1 | Pyramid + k-vector (Mortari et al. 2004; Mortari 1997/2000) | **ADOPT (idea)**: 2-star seeds through a sorted pair list (k-vector) inside the hint cone | §5.3 |
| 1.2 | Grid / polar-grid (Padgett & Kreutz-Delgado 1997; Clouse & Padgett 2000) | SKIP as a matcher; CONSIDER Clouse & Padgett's Bayesian thresholding as support for §5.4 | §5.4 |
| 1.3 | Geometric voting (Kolomenkin et al. 2008) | CONSIDER as a hinted fallback when the hash fails on false bright sources | §5.3 |
| 1.4 | Liebe triangle (1993) | SKIP (fragile with false stars) | — |
| 1.5 | SP-Search / recursive mode (Samaan, Mortari & Junkins 2005) | Already covered by §5.6; no change | §5.6 |
| 1.6 | Multi-Poles Algorithm (Schiattarella et al. 2017, 2020) | CONSIDER: evidence that magnitude filtering helps under heavy false-object load | §5.4 |
| 1.7 | Non-dimensional star-ID (Leake, Arnas & Mortari 2020) | SKIP (focal-length drift is not a problem here; the scale prior handles it) | — |
| 1.8 | Reverse-attitude statistics (Dong et al. 2024) | CONSIDER: attitude-space voting as the hinted fallback (same family as 1.3) | §5.3 |
| 1.9 | Surveys and benchmarks (Spratling & Mortari 2009; Rijlaarsdam et al. 2020; ESA Kelvins 2017) | CONSIDER the Kelvins dataset and Rijlaarsdam's protocol as public regression tests | §9 |
| 2.1 | Houtz & Frueh 2021 (polar-grid plate solving of streaked stars) | **ADOPT (idea)**: common trail geometry as a star/non-star classifier and as a prior | §5.0, §6.3a, §6.4 |
| 2.2 | Sease, Koglin & Flewelling 2013 (attitude plus rate from streak endpoints) | **ADOPT**: confirms the time-dependent WCS plan | §5.5 |
| 2.3 | AURORAS "dynamic astrometric calibration" (Krantz et al. 2025) | ADOPT (confirms §5.5); CONSIDER a quadratic time term | §5.5 |
| 2.4 | Vereš et al. 2012 trail fitting | ADOPT as producer guidance (covariance and endpoints) | §6.3, §6.9 |
| 2.5 | Rolling-shutter streak ID (Schiattarella et al. 2020) | SKIP (global shutter) | — |
| 2.6 | Astrometry24.NET (Sybilska et al. 2019) | SKIP (closed service; no algorithmic detail to take) | — |
| 3.1 | RANSAC minimal-sample argument (Fischler & Bolles 1981) | **ADOPT**: the reason to prefer 2-star seeds in hinted mode | §5.3 |
| 3.2 | PROSAC (Chum & Matas 2005) | **ADOPT**: order hypotheses by `p_star` × brightness | §5.3, §6.3 |
| 3.3 | LO-RANSAC (Chum, Matas & Kittler 2003) | **ADOPT**: refit the inliers before the final log-odds decision | §5.4, §5.5 |
| 3.4 | MAGSAC / MAGSAC++ (Barath et al. 2019/2020) | SKIP the algorithm; CONSIDER a fitted global error-scale factor | §5.4, §5.5 |
| 3.5 | A contrario / NFA framing | CONSIDER: report the number of false alarms (NFA) next to the log-odds | §5.4, §7 |
| 3.6 | Mixture-model fitting (Hogg, Bovy & Lang 2010) | **ADOPT**: one inlier/outlier model shared by the verifier and the fitter | §5.4, §5.5 |
| 3.7 | Photometric terms (Clouse & Padgett 2000; MPA) | ADOPT (already a "Should"), with heavy-tailed magnitude residuals | §5.4 |
| 4.1 | Wahba solvers (QUEST, ESOQ/ESOQ2, SVD; Markley & Mortari 2000) | SKIP QUEST/ESOQ; keep SVD (Kabsch/Umeyama); a closed form for 2-star seeds | §5.5 |
| 4.2 | SIP vs TPV vs Zernike (Shupe et al. 2005; Calabretta & Greisen 2002; SCAMP) | Keep SIP; SKIP TPV and Zernike | §5.5 |
| 4.3 | Persistent instrument distortion plus per-exposure terms (Bernstein et al. 2017; SCAMP) | **ADOPT**: fit distortion across many frames, not per frame | §5.5, §6.5, §9 |
| 4.4 | Model selection (F-test vs BIC vs cross-validation) | ADOPT BIC, or leave-one-out cross-validation, in place of the F-test | §5.5 |
| 5.1 | HEALPix uniformization (Górski et al. 2005; astrometry.net) | Keep | §5.2 |
| 5.2 | Fibonacci lattice (González 2010) / spherical spiral (Bauer 2000) | Keep for pattern fields (via cedar-solve) | §5.2 |
| 5.3 | Rigorous epoch propagation (Butkevich & Lindegren 2014) | **ADOPT** as the reference formulae, including covariance propagation | §5.1 |
| 5.4 | UCAC4 error budget (Zacharias et al. 2013) | ADOPT: fixed-epoch index, query-time propagation | §5.1, §5.2 |
| 5.5 | Relativistic aberration (Klioner 2003; SOFA/ERFA `ab`) | ADOPT the SOFA/ERFA formulation; aberration moves the field centre and changes the scale, nothing else | §6.6 |
| 6.1 | Learned star-ID (SIGNN, WACV 2025; ICPRAM 2023; RPNet) | **SKIP** (not credible for production verification) | — |

---

## 1. Star identification / lost-in-space (star-tracker literature)

**Context that drives every verdict below.** Star-tracker algorithms are tuned for wide FOVs
(8–20°), bright stars, tiny onboard catalogues and no prior. Fayde's hinted path is different: it
has a pointing prior, an often tight scale prior, sometimes a roll prior, and a narrow field crowded
with false detections. What carries over is mostly about *robustness to false stars*, and the main
lesson is simple: **the fewer true stars a hypothesis needs, the more robust the matcher.**

### 1.1 Pyramid and the k-vector

- **Citations.** Mortari, Samaan, Bruccoleri & Junkins, "The Pyramid star identification
  technique", *Navigation* 51(3), 2004. Mortari, "Search-less algorithm for star pattern
  recognition", *J. Astronautical Sci.* 45(2), 1997
  ([Springer](https://link.springer.com/doi/10.1007/BF03546375)). Mortari & Neta, "k-vector range
  searching techniques", AAS 00-128, 2000. **[secondary]** (described consistently in
  [Leake et al. 2020](https://arxiv.org/abs/2003.13736) and in surveys; originals not read.)
- **Core idea.** The k-vector is a precomputed index over the *sorted* angular distances of all
  catalogue star pairs. A range query ("every pair with separation in [d−ε, d+ε]") then costs O(1)
  plus the number of results, with no comparisons. Pyramid takes image star triples, looks each
  edge up in the k-vector, intersects the candidate sets to identify a triangle, then confirms it
  with a fourth "pyramid" star. A triangle that the fourth star cannot confirm is rejected, which
  is how false stars get discarded.
- **Evidence.** Widely called the nominal lost-in-space algorithm. In a six-algorithm comparison
  (cited in the search results; not read in full) Pyramid had the highest average accuracy with
  false stars and was the fastest. Its known weakness is a need for accurate focal length
  ([Leake et al. 2020](https://arxiv.org/abs/2003.13736)). Fayde's scale prior covers that. The
  winner of the ESA Kelvins star-ID challenge ("Super k-ID") was built on k-vector techniques
  ([results](https://kelvins.esa.int/star-trackers-first-contact/results/), [full]).
- **Verdict: ADOPT the idea, not the algorithm (§5.3).** In hinted mode, build a sorted pair list
  (a k-vector) of the catalogue stars inside the hint cone, at query time. That is a few hundred
  stars and about 10⁵ pairs, built in well under a millisecond. Then seed hypotheses from **image
  star pairs**. With a scale prior and known parity, two stars fix a full similarity transform. See
  §3.1 for why 2-star seeds beat 4-star patterns when many sources are false. Keep the tetra3 hash
  as the first, cheapest attempt.

### 1.2 Grid algorithm, and its Bayesian and polar variants

- **Citations.** Padgett & Kreutz-Delgado, "A grid algorithm for autonomous star identification",
  *IEEE TAES* 33(1), 1997. Clouse & Padgett, "Small field-of-view star identification using
  Bayesian decision theory", *IEEE TAES* 36(3), 773–783, 2000
  (doi:10.1109/7.869495). Lee et al. 2007 (polar grid).
  **[abstract]**
- **Core idea.** Rotate the field so that a pivot star's nearest neighbour lies on the x-axis,
  rasterize the neighbours onto a coarse grid, and compare bit patterns. Clouse & Padgett add a
  Bayesian, adaptive acceptance threshold aimed at a 2° FOV.
- **Evidence.** Clouse & Padgett report 96.0% correct with <0.3% false positives at 0.5 px noise
  and 0.8 mag noise. Grid methods depend on the axis-defining neighbour being real and detected,
  which false or missing stars break. Houtz & Frueh (§2.1) say this explicitly for streaked images.
- **Verdict: SKIP as a matcher.** CONSIDER Clouse & Padgett only as further support for
  likelihood-ratio acceptance with magnitude information (§3.7), which Fayde already plans.

### 1.3 Geometric voting

- **Citation.** Kolomenkin, Pollak, Shimshoni & Lindenbaum, "Geometric voting algorithm for star
  trackers", *IEEE TAES* 44(2), 441–456, 2008, doi:10.1109/TAES.2008.4560198. **[abstract]**
- **Core idea.** Every catalogue pair whose angular distance matches an image pair votes for both
  identities. Each image star takes the identity with the most votes, and a second voting round
  checks the result. The method uses geometry only and draws on the whole image rather than a
  chosen subset.
- **Evidence.** The authors report robustness to non-star objects. The LOST open-source tracker
  ([docs](https://uwcubesat.github.io/lost/classlost_1_1GeometricVotingStarIdAlgorithm.html))
  notes "no statistical guarantees" on the wrong-answer rate.
- **Verdict: CONSIDER as a hinted fallback (§5.3).** Inside a hint cone with a tight scale prior,
  pair voting is cheap and degrades gracefully when the brightest sources are false. Fayde's
  log-odds verifier provides the statistical guarantee that voting lacks. Run it only after the
  hash and the 2-star seeds have failed.

### 1.4 Liebe triangle

- **Citation.** Liebe, "Pattern recognition of star constellations for spacecraft applications",
  *IEEE AES Magazine* 8(1), 31–39, 1993, doi:10.1109/62.180383. **[secondary]**
- **Core idea.** A star plus its two nearest neighbours gives a triangle feature (two distances and
  an angle).
- **Evidence.** The ICPRAM 2023 related-work section and others: fast, but it "suffered limitations
  with the introduction of false stars". A single false nearest neighbour breaks it.
- **Verdict: SKIP.**

### 1.5 SP-Search and recursive mode

- **Citation.** Samaan, Mortari & Junkins, "Recursive mode star identification algorithms",
  *IEEE TAES* 41(4), 2005. **[secondary]**
- **Core idea.** Given the previous attitude, predict which catalogue stars fall in the FOV and
  match locally, with no pattern search.
- **Verdict:** this is §5.6 sequence mode. No change.

### 1.6 Multi-Poles Algorithm (MPA)

- **Citations.** Schiattarella, Spiller & Curti, "A novel star identification technique robust to
  high presence of false objects: the Multi-Poles Algorithm", *Adv. Space Res.* 59(8), 2133–2147,
  2017. Schiattarella, Spiller & Curti, "Star identification robust to angular rates and false
  objects with rolling shutter compensation", *Acta Astronautica*, 2020
  ([preprint](https://iris.uniroma1.it/bitstream/11573/1349415/4/Schiattarella_preprint_star-identification_2020.pdf),
  **[full]** for the 2020 paper, which summarizes the 2017 one).
- **Core idea.** Several "pole" stars are each tried as the centre of a pair-distance pattern, in
  three phases: acceptance, verification, confirmation.
- **Evidence.** 100% success with up to 185 false objects against 33 catalogued stars, falling to
  66% at 486 false objects against 31 stars (stable attitude). Second in the ESA Kelvins challenge.
  The ESA-contest version's main improvement was **magnitude filtering**, which made the
  confirmation phase unnecessary.
- **Verdict: CONSIDER, as evidence.** The magnitude-filtering result is independent support for
  a photometric term in the verifier (§3.7). The algorithm itself adds nothing beyond 1.1 and 1.3.

### 1.7 Non-dimensional star-ID

- **Citation.** Leake, Arnas & Mortari, "Non-dimensional star-identification", *Sensors* 20(9):2697,
  2020 ([arXiv:2003.13736](https://arxiv.org/abs/2003.13736)). **[abstract]**
- **Core idea.** Use dihedral angles of spherical star triangles, which do not depend on focal
  length, so identification survives a miscalibrated focal length or optical axis.
- **Verdict: SKIP.** Fayde estimates scale itself, and the astrometry.net quad codes are already
  scale-invariant.

### 1.8 Reverse attitude statistics

- **Citation.** Dong, Wang, Wang & Wang, "Reverse attitude statistics based star map identification
  method", arXiv:2410.23758, 2024 (submitted to *IEEE Sensors J.*)
  ([arXiv](https://arxiv.org/abs/2410.23758)). **[abstract]**
- **Core idea.** Build candidate attitudes from many matched pairs (with spatial-hash indexing),
  then take the mode of the attitude distribution. Matching and attitude solving happen together.
  Aimed at missing stars, false stars and *star trailing* in near-space conditions.
- **Evidence.** The authors claim a >14.3% better identification rate and >28.5% faster solving
  than the state of the art, with field and on-orbit tests. The baselines were not visible in the
  abstract.
- **Verdict: CONSIDER**, as the same family as 1.3: attitude-space voting. If Fayde implements the
  hinted fallback, voting in the (roll, offset) space of 2-star seeds is the natural form.

### 1.9 Surveys and benchmarks

- **Citations.** Spratling & Mortari, "A survey on star identification algorithms", *Algorithms*
  2(1), 2009. Rijlaarsdam, Yous, Byrne, Oddenino, Furano & Moloney, "A survey of lost-in-space star
  identification algorithms since 2009", *Sensors* 20(9):2579, 2020 (**[abstract]**; the MDPI and
  PMC full text was blocked to automated fetch). Toloei et al., "A comparative analysis of star
  identification algorithms", *Astrophys. Space Sci.*, 2020
  ([Springer](https://link.springer.com/article/10.1007/s10509-020-03775-9), **[abstract]**). ESA
  Kelvins "Star Trackers: First Contact" (2016–2017)
  ([data](https://kelvins.esa.int/star-trackers-first-contact/data/),
  [results](https://kelvins.esa.int/star-trackers-first-contact/results/), **[full]**).
- **Findings.**
  - Rijlaarsdam et al. find that published comparisons "can lead to inconsistent conclusions" and
    give considerations for simulation-based evaluation.
  - Kelvins is the only shared public benchmark found: 10,000 scenes, 10° FOV, 1920×1440, up to
    **80 false artifacts per scene**, centroid noise under 1 px, magnitude noise, Hipparcos.
  - Final Kelvins accuracy scores out of a perfect 9921: Super k-ID 9708, MPA 9703, **"Tetra"
    (the original hash-table approach behind tetra3) 9672**, then others from 9605 down.
  - So a hash-table matcher is within about 0.4% of the best on a false-star-heavy benchmark.
- **Verdict: CONSIDER.** The Kelvins scenes are synthetic, public and false-star heavy, which
  suits a reproducible regression test of the pattern-hash and seed matchers. LOST
  ([GitHub](https://github.com/uwcubesat/lost)) implements Pyramid, geometric voting and
  tetra-style matchers, an independent baseline.

---

## 2. Matching trailed / streaked stars

**Overall finding.** The literature is thin. No paper found does blind or hinted plate solving of
*trailed* stars with per-source covariance, endpoints and a time-dependent WCS together. Fayde's
§5.5/§6.3 design is already ahead of the published SSA pipelines, which mostly either stack the
frame so stars become static or feed streak centroids to astrometry.net. The useful contributions
are narrower; they are listed below.

### 2.1 Houtz & Frueh: streak detection and plate solving

- **Citation.** Houtz & Frueh, "Efficient streak detection and plate-solving techniques", Proc. 8th
  European Conf. on Space Debris, ESA SDO, 2021
  ([PDF](https://conference.sdo.esoc.esa.int/proceedings/sdc8/paper/76/SDC8-paper76.pdf),
  **[full]**).
- **Core idea.**
  - Streaks are detected with cumulative-sum line integration.
  - For plate solving, only streaks whose length and orientation match the frame's most common
    values are kept as stars.
  - A polar-grid pattern is built around the brightest star, with the second-brightest streak as
    the axis, using spherical-triangle radii and angles. It is matched by rotational correlation
    against catalogue patterns, centred on Tycho-2 stars near the a-priori pointing.
  - All stars are then identified from the radius/angle residuals, and QUEST gives the attitude at
    mid-exposure.
- **Evidence.** 15 frames from the POGS 35 cm telescope: on average 11.7 streaked stars identified
  per frame, with star-vector errors under 1″. A small sample with no false-solve statistics.
- **Verdict: ADOPT the idea, not the matcher.** This is the closest published analogue to Fayde's
  case. It confirms that (a) trail length and direction consistency is the strongest
  star/non-star cue (§5.0 already does this, via `p_star`), and (b) a hinted, pointing-centred
  search is the practical approach. Fayde can go one step further. In target-tracking mode the
  common star trail vector (`trail_model`) is the sky's apparent motion relative to the tracked
  target, projected onto the detector. Given the commanded tracking rate (from the TLE and observer
  state) and the exposure time, the trail's **direction gives roll** (with its sign, **parity**) and
  its **length gives scale**. This is Fayde's own inference, not a published method; Fayde does
  not derive priors from it.

### 2.2 Sease et al.: attitude and rate from streak endpoints

- **Citation.** Sease, Koglin & Flewelling, "Long integration star tracker image processing for
  combined attitude–attitude rate estimation", *Proc. SPIE* 8739, 2013, doi:10.1117/12.2020004
  ([UCF record](https://stars.library.ucf.edu/scopus2010/6089)). **[abstract]**, plus a
  catalogue-free follow-up and a UCF thesis ([core.ac.uk](https://core.ac.uk/download/236295411.pdf),
  **[secondary]**).
- **Core idea.** Treat streak endpoints as corner features. Group the endpoints into "start" and
  "end" sets using a cross-product test of the rotation sense, then run ordinary star-ID and
  attitude determination on each set. Two attitudes from one image give the attitude rate.
- **Evidence.** Monte Carlo: endpoint accuracy better than 0.1 px when more than 20° from the
  rotation axis, for streaks up to 40% of the FOV. The thesis reports rates within 2% in
  simulation.
- **Verdict: ADOPT. It confirms §5.5.** Fayde's joint fit of the mid-exposure WCS plus pointing
  and roll rates is the least-squares generalization of this, and it is better: one fit, with
  covariances, instead of two independent solves. Note the asymmetry: endpoint error *along* the
  trail is usually much larger than across it (§2.4), so endpoint covariances must be anisotropic.

### 2.3 AURORAS: time-dependent astrometric calibration

- **Citation.** Krantz, Bloch, Serda, Flechtner & McMahon-Crabtree, "AURORAS: orbit determination
  with just one look", AMOS 2025
  ([PDF](https://amostech.com/TechnicalPapers/2025/Astrodynamics/Krantz.pdf), **[full]**, relevant
  sections).
- **Core idea.** Event-camera or fast-frame data are shift-stacked so the stars become static and
  can be plate solved against Gaia. Then a "dynamic astrometric calibration" fits an
  (x, y, t) → (RA, Dec, t) transform. Residuals before the time-dependent fit show offsets and
  curvature from the moving tangent plane plus mount jitter; after it they are flat.
- **Evidence.** Qualitative residual plots; no numbers on astrometric accuracy were given.
- **Verdict: ADOPT, as confirmation of §5.5.** CONSIDER an optional **quadratic time term**
  (pointing acceleration) when trail `profile` or curvature data show non-linear trails, as in
  long exposures tracking fast LEO targets.

### 2.4 Trail fitting for astrometry (producer side)

- **Citation.** Vereš, Jedicke, Denneau, Wainscoat, Holman & Lin, "Improved asteroid astrometry
  and photometry with trail fitting", *PASP* 124, 1197, 2012
  ([arXiv:1209.6106](https://arxiv.org/abs/1209.6106)).
  **[abstract]**
- **Core idea.** Model a trail as a Gaussian PSF convolved with a line segment, and fit its
  centre, length, angle and flux.
- **Evidence.** About 3× better astrometry and up to 2 mag better photometry than PSF or aperture
  methods for trails longer than about 10 px (Pan-STARRS1).
- **Verdict: ADOPT, as guidance in §6.3/§6.9.** Recommend that producers (streak detectors
  included) report centroids, endpoints and covariances from a trail-model fit. The
  resulting covariance (narrow across the trail, wide along it) is exactly what §5.4/§5.5 consume.
  No Fayde code change.

### 2.5 Rolling-shutter streaks

- **Citation.** Schiattarella et al. 2020 (§1.6), **[full]**. Curti et al. (rate estimation from
  streaks up to 10–25°/s; [Sapienza](https://iris.uniroma1.it/handle/11573/852340),
  **[abstract]**).
- **Finding.** A rolling shutter distorts streak length by about 25% (Curti et al.). A first-order
  compensation gives reliable identification up to 5°/s.
- **Verdict: SKIP.** Fayde's main target uses a global-shutter CCD.

### 2.6 Astrometry24.NET

- **Citation.** Sybilska et al., "Astrometry24.NET – precise astrometry for SST and NEO", Proc. 1st
  NEO and Debris Detection Conf., ESA, 2019
  ([PDF](https://conference.sdo.esoc.esa.int/proceedings/neosst1/paper/147/NEOSST1-paper147.pdf),
  **[full]** for the introduction and architecture).
- **Finding.** A commercial cloud service that fits streak profiles and solves frames tracked at
  non-sidereal rates. The paper describes the architecture and accuracy campaigns but gives little
  detail on the matching algorithm.
- **Verdict: SKIP.** Nothing to port. It could serve as an external accuracy cross-check if
  real frames were ever run through it.

---

## 3. Verification, false-match probability and robust sampling

### 3.1 The minimal-sample argument (RANSAC)

- **Citation.** Fischler & Bolles, "Random sample consensus", *Comm. ACM* 24(6), 1981.
  **[secondary; standard result]**
- **Core point.** If a fraction *w* of sources are real catalogue stars, a hypothesis built from *k*
  sources is all-true with probability about *wᵏ*. At *w* = 0.5, a 4-star pattern (tetra3 or
  astrometry.net) succeeds 6.3% of the time, a 2-star seed 25%, and a 1-star seed 50%. At
  *w* = 0.3 the figures are 0.8%, 9% and 30%. For a 30–60 source xylist with many false
  detections, this ratio matters more than any detail of hash design.
- **Verdict: ADOPT (§5.3).** When priors fix scale (and parity), the hinted path should seed from
  2-star pairs (§1.1). When roll is also known, for example from spacecraft attitude or the trail
  direction (§2.1), seed from single stars: each image–catalogue pairing proposes a translation,
  and the translations vote. The 4-star hash remains the blind and weak-prior path.

### 3.2 PROSAC

- **Citation.** Chum & Matas, "Matching with PROSAC – progressive sample consensus", CVPR 2005,
  doi:10.1109/CVPR.2005.221
  ([PDF](https://dspace.cvut.cz/bitstream/10467/9496/1/2005-Matching-with-PROSAC-progressive-sample-consensus.pdf)).
  **[abstract]**, plus reviews.
- **Core idea.** Sample hypotheses from a progressively growing set of the top-ranked
  correspondences rather than uniformly, which is still correct in the limit. The gain depends
  entirely on how well the ranking predicts correctness.
- **Evidence.** Up to more than 100× fewer samples than RANSAC on CV benchmarks (in one example,
  9 samples against about 10⁵); no advantage when the ordering is random.
- **Verdict: ADOPT (§5.3, §6.3).** tetra3/cedar's brightest-first, breadth-first enumeration is
  already a PROSAC-like ordering, but it ranks by brightness alone. Fayde has `p_star` and the
  §5.0 trail-consistency adjustments, which predict "real star" much better than flux alone (the
  brightest sources in SSA frames are often the target, cosmic rays or blends). Rank sources by a
  score such as `p_star` × (normalized brightness rank), and feed all matchers in that order. The
  ordering must be deterministic (requirement 6).

### 3.3 LO-RANSAC (local optimization)

- **Citation.** Chum, Matas & Kittler, "Locally optimized RANSAC", DAGM 2003, LNCS 2781.
  **[secondary; standard]**
- **Core idea.** When a hypothesis looks promising, refit the model on its inliers (and iterate)
  before scoring it. Minimal-sample models are noisy, so scoring the raw model underestimates
  support.
- **Verdict: ADOPT (§5.4/§5.5).** This matters most for trailed stars and for 1- or 2-star seeds.
  A similarity transform from two stars, ignoring trail rate and distortion, can leave a true match
  just below threshold. Fayde's flow (§5.4):
  1. A cheap pre-filter.
  2. If it passes, a weighted refit on its matches (similarity → TAN).
  3. Re-match, then compute the log-odds on the **refined** model.

  astrometry.net does a version of this ("tweak" during verification).

### 3.4 MAGSAC / MAGSAC++

- **Citation.** Barath, Matas & Noskova, "MAGSAC: marginalizing sample consensus", CVPR 2019;
  Barath, Noskova, Ivashechkin & Matas, "MAGSAC++, a fast, reliable and accurate robust estimator",
  CVPR 2020 ([arXiv:1912.05909](https://arxiv.org/abs/1912.05909)). **[abstract]**
- **Core idea.** Remove the hard inlier threshold by marginalizing the model score over the noise
  scale σ, up to a loose upper bound, and use weighted least squares (IRLS) in place of a hard
  inlier set.
- **Verdict: SKIP the algorithm.** Fayde has per-source covariances, so the log-odds already uses
  an explicit noise model, which is what MAGSAC approximates when one is missing. **CONSIDER the
  underlying idea:** producers' covariances will often be miscalibrated, so fit a single global
  error-scale factor *s* (covariance × *s*², with a prior on *s* near 1), or marginalize over it
  in the verifier. Report *s* in the diagnostics, so a detector that under-reports its errors is
  visible.

### 3.5 A contrario / number of false alarms (NFA)

- **Citations.** The a-contrario detection framework (Desolneux, Moisan & Morel; Abergel & Moisan
  on trajectories). **[secondary]** No paper was found that applies NFA to star identification.
- **Core idea.** Accept a detection when (number of tests) × P(an event at least this good under
  the null hypothesis) < ε. The threshold is calibrated against the *search size*, not set per
  hypothesis.
- **Verdict: CONSIDER (§5.4, §7).** astrometry.net's 10⁹ odds is a fixed per-hypothesis threshold.
  Blind solves test many more hypotheses than hinted ones, so the family-wise false-accept rate
  differs between them. Cheap improvement: count the hypotheses tested and report
  `NFA ≈ n_tested × exp(−log_odds)` (an upper-bound style figure) next to the log-odds. Optionally
  require NFA < ε in blind mode. This also gives a calling pipeline a single number to record.

### 3.6 Mixture-model fitting

- **Citation.** Hogg, Bovy & Lang, "Data analysis recipes: fitting a model to data",
  arXiv:1008.4686, 2010. **[secondary; well-known]**
- **Core idea.** Model each data point as a mixture: inlier (Gaussian with its covariance) or
  outlier (broad background), and fit or marginalize rather than hard-clip at kσ.
- **Verdict: ADOPT (§5.4/§5.5).** This is the *same* likelihood the astrometry.net verifier uses (a
  Gaussian match against a uniform distractor). Use it in the fitter as well, with per-source
  outlier priors from `p_star`. The match loop in §5.5 then becomes IRLS on the mixture weights,
  replacing "re-match within kσ". One noise model, used consistently, is easier to reason about
  and gives better-calibrated covariances.

### 3.7 Photometric consistency in verification

- **Evidence.** Clouse & Padgett (1.2) use 0.8 mag noise inside a Bayesian decision. The
  MPA authors (1.6) found magnitude filtering was the main improvement for their contest entry. A
  brightness-referenced method ([Sensors 2014](https://dx.doi.org/10.3390/s141018498),
  **[abstract]**) reports speed gains from brightness ordering.
- **Verdict: ADOPT (already a "Should" in §2 and §5.4).** Specifics from the literature and from
  the regime:
  - Fit the zero point robustly.
  - Use a **heavy-tailed** (Student-t) magnitude residual, because UCAC4 `magm` against an
    unfiltered or wide band, trail photometry, saturation and blends all produce outliers.
  - Weight the magnitude term well below the positional term until it has been calibrated on a
    corpus of real frames.

  The magnitude term is a second discriminator, not a primary one.

---

## 4. WCS fitting

### 4.1 Wahba's problem: QUEST, ESOQ, SVD

- **Citations.** Shuster & Oh, "Three-axis attitude determination from vector observations",
  *J. Guid. Control* 4(1), 1981 (QUEST). Mortari, "ESOQ: a closed-form solution to the Wahba
  problem", *J. Astronaut. Sci.* 45(2), 1997. Markley & Mortari, "Quaternion attitude estimation
  using vector observations", *J. Astronaut. Sci.* 48(2–3), 359, 2000, doi:10.1007/BF03546284 (preprint, under the title "New
  developments in quaternion estimation from vector observations":
  [NTRS](https://ntrs.nasa.gov/citations/20000034107)). **[abstract]** Umeyama,
  "Least-squares estimation of transformation parameters between two point patterns", *IEEE
  TPAMI* 13(4), 1991. **[secondary; standard]**
- **Core finding.** Davenport's q-method and SVD are the most robust; QUEST, ESOQ and ESOQ2 are
  faster and adequate for star sensors. Robustness only matters when measurement accuracies differ
  widely (Markley & Mortari).
- **Verdict.**
  - **SKIP QUEST/ESOQ.** At 30–60 stars, speed is irrelevant. Fayde's measurement accuracies *do*
    differ widely and are anisotropic (trails), and none of these solvers handles anisotropic
    weights anyway.
  - **Keep SVD (Kabsch, or Umeyama for a similarity with scale)** only for the initial rotation,
    then go straight to the weighted TAN least-squares fit, which handles full 2×2 covariances.
  - For 2-star seeds (§3.1) use the closed-form two-point similarity in the tangent plane.
  - The design's parity check from the SVD determinant is sound at 1° FOV, because star spread
    (~10⁻² rad) is far larger than noise (~10⁻⁵ rad). For tiny or degenerate configurations (all
    stars nearly collinear), fit both parities and compare the log-odds instead.

### 4.2 Distortion representation: SIP vs TPV vs Zernike

- **Citations.** Shupe et al., "The SIP convention for representing distortion in FITS image
  headers", ASP Conf. Ser. 347, 2005. Calabretta & Greisen, "Representations of celestial
  coordinates in FITS", *A&A* 395, 2002 (TPV is a registered convention, not part of the
  standard). Bertin, "Automatic astrometric and photometric calibration with SCAMP", ASP Conf. Ser.
  351, 2006. Shupe et al., SPIE 2012 (`sip2pv`/`pv2sip` conversion). **[secondary]**, from survey
  pipeline papers ([PTF](https://arxiv.org/pdf/1404.1953), [SAGES](https://arxiv.org/pdf/2306.15611)).
- **Finding.** SIP (pixel-space polynomials) and TPV (polynomials in intermediate world
  coordinates) are interconvertible to good precision. SCAMP writes TPV, and pipelines routinely
  convert it to SIP because most tools read SIP. Zernike is useful only for wide-field optics with
  a circular pupil, and has no FITS convention.
- **Verdict: keep SIP; SKIP TPV and Zernike (§5.5, "Won't" confirmed).**

### 4.3 Persistent instrument distortion with per-exposure terms

- **Citations.** Bernstein et al., "Astrometric calibration and performance of the Dark Energy
  Camera", *PASP* 129:074503, 2017 ([arXiv:1703.01679](https://arxiv.org/abs/1703.01679),
  **[abstract]**). SCAMP (Bertin 2006), which fits distortion across many exposures at once.
- **Core idea.** Split the astrometric model into a *persistent* instrument component (optics,
  detector), fitted once from many exposures, and a *per-exposure* component, kept low-order.
  DECam reaches 3–6 mas per exposure this way.
- **Verdict: ADOPT (§5.5, §6.5, §9).** With 30–60 stars per frame, many of them trailed with large
  along-trail errors, a per-frame SIP order 3 (14 free coefficients on top of 6 linear ones) is at
  serious risk of overfitting. Fayde already has `instrument.distortion` in §6.5; make it the
  normal case:
  1. Add `fayde refit --calibrate`: a joint fit over N solved frames that estimates the persistent
     SIP plus a per-frame affine part.
  2. Store the result in the `InstrumentProfile`.
  3. Per-frame fits then solve only the linear terms, plus at most a low-order correction chosen
     by model selection (4.4).

  This also supplies the `instrument.distortion` input listed in §6.9.

### 4.4 Model selection and covariance

- **Finding (standard statistics).** The F-test assumes nested models with Gaussian, correctly
  scaled errors. With miscalibrated covariances it over-selects. BIC penalizes parameters by
  log *N*. Leave-one-out cross-validation on the matched stars is cheap at N ≤ 60 and robust to
  error-scale problems.
- **Verdict: ADOPT (§5.5).**
  - Choose the SIP order by BIC, or by leave-one-out residuals, instead of an F-test.
  - Report both the formal covariance, (JᵀWJ)⁻¹, and the covariance scaled by the reduced χ² (or
    by the fitted *s*² from 3.4).
  - Include the catalogue position error at the epoch in W (the design already does this), so the
    covariance is not optimistic when UCAC4 dominates the error.

---

## 5. Catalogue and index design

### 5.1 HEALPix uniformization

- **Citations.** Górski et al., "HEALPix", *ApJ* 622, 2005. Lang et al. 2010
  (`uniformize-catalog`).
- **Verdict: keep (§5.2).** Nothing better found.

### 5.2 Fibonacci lattice and spherical-spiral catalogues

- **Citations.** González, "Measurement of areas on a sphere using Fibonacci and
  latitude–longitude lattices", *Math. Geosci.* 42, 2010
  ([arXiv:0912.4540](https://arxiv.org/abs/0912.4540), **[abstract]**). Bauer, "Distribution of
  points on a sphere with application to star catalogs", *J. Guid. Control Dyn.* 23(1), 2000
  (**[secondary]**; the spherical spiral method). Variants: charged particles, k-means clustering
  (GWKM: Somayehee, Nikkhah & Roshanian 2019, *J. Navig.* 72(4), 948,
  doi:10.1017/S0373463318001029, **[abstract]**).
- **Finding.** Fibonacci lattices give near-equal-area cells for any point count, with at least
  40% lower RMS error in area estimates than latitude–longitude grids. Spherical-spiral catalogues
  are the star-tracker standard for uniform onboard catalogues.
- **Verdict: keep (§5.2).** cedar-solve's Fibonacci lattice fields for pattern generation plus
  HEALPix for star selection is a sound combination. No change.

### 5.3 Rigorous epoch propagation

- **Citation.** Butkevich & Lindegren, "Rigorous treatment of barycentric stellar motion",
  *A&A* 570, A62, 2014 ([arXiv:1407.4664](https://arxiv.org/abs/1407.4664), **[abstract]**). The
  same formulae appear in the
  [Gaia DR3 documentation](https://gea.esac.esa.int/archive/documentation/GDR3/Data_processing/chap_cu3ast/sec_cu3ast_intro/ssec_cu3ast_intro_tansforms.html).
- **Core idea.** Closed-form propagation of all six astrometric parameters between epochs,
  reversible, with a full **covariance** transformation.
- **Verdict: ADOPT as the reference (§5.1).** The design's "full 3-D space-motion vector" should
  point to these equations. With UCAC4 (no parallax or radial velocity) they reduce to the
  zero-parallax case. With Gaia DR3 they are exact. Use the covariance transformation to produce
  `pos_err` at the epoch, which feeds the fit weights in §5.5.

### 5.4 UCAC4 error budget and epoch policy

- **Citation.** Zacharias et al., "The fourth US Naval Observatory CCD Astrograph Catalog
  (UCAC4)", *AJ* 145:44, 2013 ([arXiv:1212.6182](https://arxiv.org/abs/1212.6182)). **[abstract]**
- **Finding.**
  - Formal position errors are 15–100 mas at the mean epoch, plus 10–20 mas systematic.
  - Proper-motion systematics are about 1–4 mas/yr; random errors are larger for faint stars.
  - Propagated to 2026 (about 26 yr from 2000, somewhat less from UCAC4's central epochs), typical
    position errors are therefore about 0.05–0.3″. At 3″/px that is **0.02–0.1 px**: irrelevant
    for matching and small for the fit.
  - Gaia DR3 is about 100× better and would only matter for sub-0.1 px accuracy claims.
- **Verdict: ADOPT a two-level policy (§5.1, §5.2).**
  1. **Build the index at a fixed epoch** (for example 2025.0). For 99.9% or more of stars,
     |µ| < 100 mas/yr, so ±10 yr moves them less than 1″ (0.3 px), well within the code
     tolerance.
  2. **Exclude high proper-motion stars** (for example |µ| × Δt_max > 0.3 px) from *patterns*, but
     keep them in the star table.
  3. **Propagate at query time**, to the mid-exposure epoch, for verification and fitting.

  Rebuilding the index for each campaign is then unnecessary. Gaia DR3 matters only for
  accuracy requirements below about 0.1″.

### 5.5 Aberration for spacecraft observers

- **Citations.** Klioner, "A practical relativistic model for microarcsecond astrometry in space",
  *AJ* 125, 2003. **[secondary]** SOFA/ERFA `iauAb`/`eraAb` (fully relativistic to about 1 µas).
  Gravity Probe B context: orbital aberration about ±5.2″, annual about ±20.5″ (via search,
  **[secondary]**).
- **Finding (computed).** For a telescope in low Earth orbit, β = |v_annual + v_orbital|/c is at
  most about 1.25×10⁻⁴.
  That shifts the field centre by up to about 26″ (the design says about 23″, which is consistent).
  The *differential* effect across a 1° field is a scale change of order β, at most about 0.2″ at
  the field edge, about 0.07 px.
- **Verdict: ADOPT ERFA's `ab` formulation in `fayde-ephem` (§6.6).** Note in the design that, at
  this FOV, aberration affects only CRVAL and the plate scale. Matching and verification can safely
  ignore it, and it matters only for the `output_frame` choice. That makes `observer` purely an
  output-accuracy input, never a solve-rate input. Also say this in the §6.6 "Effect" column so
  producers know what they gain.

### 5.6 Index compression

- **Finding.** No result found that improves on cedar-solve's layout (16-bit key check plus
  float16 edge in each slot) for Fayde's scale. Fayde's 1° full-sky band is modest anyway: about
  10⁵–10⁶ stars and 50 patterns per lattice field. No change.

---

## 6. Learned / ML star identification

- **Citations.**
  - Hepburn-Dickins, Jones, Edwards, Morgan & Bell, "SIGNN – star identification using graph neural networks",
    WACV 2025
    ([CVF page](https://openaccess.thecvf.com/content/WACV2025/html/Hepburn-Dickins_SIGNN_-_Star_Identification_using_Graph_Neural_Networks_WACV_2025_paper.html),
    **[abstract]**). A graph-attention network over a Hipparcos star graph with edges below 3°.
  - Hepburn-Dickins & Edwards, "Noise robustness of data-driven star classification", ICPRAM 2023
    ([PDF](https://scitepress.org/publishedPapers/2023/118040/pdf/index.html), **[full]**).
  - RPNet (Xu et al. 2019), 1-D CNN classifiers (Yang et al. 2022), and Wang et al., *Sensors*
    2021 (**[secondary]**).
- **Evidence.**
  - ICPRAM 2023: 98.5% accuracy on augmented samples, but **13.5% of stars misidentified** on
    simulated Stellarium images.
  - A 0.999 confidence cut removed 99.85% of wrong identifications but discarded 39% of correct
    ones.
  - No head-to-head comparison against Pyramid or tetra3.
  - Most learned methods are trained for wide-FOV Hipparcos star trackers. None addresses trailed
    stars, sub-degree FOVs with UCAC4 density, or calibrated false-accept rates.
- **Verdict: SKIP.** Not credible for production where a wrong WCS is worse than no WCS. Fayde's
  correctness rests on a calibrated verifier, which ML methods do not provide. If anything, a
  learned *source classifier* belongs on the producer side, feeding `p_star`, which the contract
  already supports.

---

## 7. Source and streak detection (2026-11-16)

Survey made when detection was brought into Fayde for frames with point-like stars
(`docs/design/detection.md`). Reference implementations were read at the head of their default
branches; **[full]**, **[abstract]** and **[secondary]** as above, **[code]** for source code.

### 7.1 Star detection in use today

| Tool | Background / noise | Filter | Threshold | Shape tests | Centroid | Licence |
|---|---|---|---|---|---|---|
| astrometry.net `simplexy` (Lang et al. 2010) **[code]** | Median in 201 px boxes on a grid; noise from pixel pairs 5 px apart | Gaussian σ = 1 px, fixed | 8σ (4σ for 8-bit) | none; peaks split at a 5σ saddle | Parabolas on the 3 × 3 of the smoothed image | BSD-3 own code, GPLv3 as distributed |
| SExtractor / SEP (Bertin & Arnouts 1996; Barbary 2016) **[code]** | 64 px mesh, clipped mode `2.5·median − 1.5·mean`, 3 × 3 median, spline | 3 × 3 kernel "FWHM 2 px", fixed | 1.5σ on 5 connected pixels | elongation from second moments; multi-threshold deblending | Isophotal or Gaussian-windowed moments (3 to 5 iterations) | GPLv3 / LGPLv3 |
| DAOFIND, photutils (Stetson 1987) **[code]** | Zero-sum kernel fits a local sky | Truncated Gaussian of a **user-given** FWHM | 4 to 5σ | sharpness 0.2 to 1.0, roundness ±1 | Marginal Gaussian fits | BSD-3 (photutils) |
| tetra3 (ESA) **[code]** | 25 px box mean; global RMS (not robust) | none | 2σ, binary opening, area 5 to 100 px | optional axis ratio | Centre of gravity of the region | Apache-2.0 |
| cedar-detect (smroid) **[code]** | Darkest of 15 strips | none; 7-pixel row gate, then rings around a blob; manual 2 to 8 × binning | 8σ typical | hot pixel if its four neighbours hold under half its excess; perimeter must be uniform | Centre of gravity, core + 2 px | FSL-1.1 (source-available, not open source); 8-bit only |
| ASTAP (H. Kleijn) **[code]** | Histogram mode, global | none | ≥ 3.5σ pixel, aperture SNR > 10 | HFD limits, fill fraction, single-pixel reject | Centre of gravity above 3σ | MPL-2.0 |
| Siril, PixInsight **[code]** | Global median + noise; Gaussian high-pass | Gaussian σ = 2 px fixed; none | 5σ; 3σ then SNR | roundness; kurtosis, coverage | Gaussian or Moffat fit (Levenberg–Marquardt) | GPLv3; PCL licence |

Published speed: SExtractor "up to about 50 Mpixel/s"; cedar-detect "less than 10 ms" per
megapixel on a Raspberry Pi 4B. None for the others.

**What none of them does**: measure the star width and adapt the filter to it; report a per-star
position error to the matcher; model cosmic rays or hot pixels other than by minimum area or
sharpness; find the centre of a saturated star with a bleed trail; trace streaks before solving
(astrometry.net removes sources in over-dense rows and columns only).

### 7.2 Centroid accuracy

- Cramér–Rao bound, Gaussian star: `√2·σ / SNR` per axis background-limited (0.60 FWHM / SNR),
  `σ / √N` source-limited (King 1983 **[full]**; Mendez, Silva & Lobos 2013 **[full]**; Mortensen
  et al. 2010 **[full]**, whose eq. 6 gives `(16/9)·σ²/N + 8π·σ⁴·b/N²` for a Gaussian-weighted
  centroid).
- A matched filter followed by a 3 × 3 quadratic fit reaches the bound for FWHM 2.8 to 5.6 px at
  SNR 5 to 100; a fixed kernel does not at the extremes; a 7 × 7 centre of light is several times
  worse (Vakili & Hogg 2016 **[full]**).
- A parabola through three samples of a Gaussian peak is biased towards the pixel centre; the
  same fit on the logarithms has no bias (Naidu & Fisher 1991 **[full]**).
- Undersampled stars (σ 0.66 px): moment methods have pixel-phase floors of 0.02 to 0.04 px; only
  a pixel-integrated profile fit goes lower (Fialho 2025, arXiv:2507.17560 **[full]**, preprint).

**Verdict: ADOPT** the matched filter with a log-quadratic 3 × 3 fit; derive the error per source.

### 7.3 Statistics of noise peaks

- Height distribution and density of local maxima of smoothed Gaussian noise: Cheng & Schwartzman
  2015, 2017, 2018 **[full]**; Bond & Efstathiou 1987 **[full]**. For a Gaussian kernel of σ = s,
  maxima above `u` per unit area tend to `u·exp(−u²/2) / (2^{5/2}·π^{3/2}·s²)`, the figure LSST
  uses for false-peak counts (DMTN-006 **[full]**).
- On pixels, with a kernel of about 1 px, the continuous formula is 25 to 35% conservative (Lin,
  Schwartzman & Davenport 2025, arXiv:2501.13239 **[full]**, preprint).
- False discovery rate for astronomical sources: Miller et al. 2001 **[full]**, Hopkins et al.
  2002 **[full]** (pixel-level); local rate: Efron 2008 **[full]**. No published use of the
  peak-height test for optical star detection was found (it is used on CMB and radio maps).

**Verdict: ADOPT** the peak density as a threshold and null model, bounded by the pixel count.

### 7.4 Streaks

| Method | Cost | Reaches | Note |
|---|---|---|---|
| Fast Radon transform (Nir, Zackay & Ofek 2018 **[full]**) | 0.9 s per megapixel | 0.1 × the pixel noise for a long streak | Stars must be removed first; GPL-3 |
| Median Radon (Stark et al. 2022 **[full]**), Fourier-slice Radon (Hickson 2018 **[full]**) | 3 to 10 s per megapixel | 0.13 to 0.4 × the pixel noise | |
| Canny + Hough (`acstools.satdet` **[code]**); kernel-based Hough (Rubin `maskStreaks`, DMTN-197 **[full]**) | 4 to 22 s per image | the per-pixel detection threshold | Needs long streaks |
| Blobs, moments and a line fit (ZTF, PTF **[full]**; ASTRiDE **[code]**) | 0.6 to 25 s per image | 1 to 3σ per pixel | Many false candidates |
| Grouping point detections into lines (Rubin `FindGlintTrails` **[code]**; Storkey et al. 2004 **[full]**) | milliseconds | detected fragments | |
| Hessian line points (Steger 1998 **[secondary]**; Lowe 2004 **[full]** for the curvature ratio) | per pixel | per-pixel contrast | The local test Fayde uses at peaks |
| Neural networks (DeepStreaks, MaxiMask, ASTA **[abstract]**) | 1 to 25 s per frame on a CPU | | Outside a few ms per megapixel |

No published pipeline combines a per-peak ridge test with a trace from that peak, and no plate
solver runs a streak detector before solving.

**Verdict: ADOPT** ridge-like peaks as seeds and trace from them (Fayde's own tracer); **SKIP**
Radon and Hough (cost; they find streaks too faint to disturb the star list).

## Recommended changes to the design doc (ranked)

1. **§5.3: add a minimal-sample hinted matcher in front of (or beside) the 4-star hash.** When the
   scale (and parity) prior is tight:
   - build a sorted pair list (k-vector) of the catalogue stars in the hint cone at query time;
   - seed similarity hypotheses from image-star **pairs**, and from single stars when roll is
     known;
   - verify each with the existing verifier.

   With half the sources false, a 2-star seed is all-true 4× as often as a 4-star pattern, and a
   1-star seed 8× as often (§3.1). Pair or attitude voting (§1.3, §1.8) is the fallback. The 4-star
   hash stays the primary path for weak priors, and the quad kd-tree stays the blind path.
   Basis: Fischler & Bolles 1981; Mortari et al. 2004; Kolomenkin et al. 2008.
2. **§5.3, §6.3: PROSAC-style ordering.** Order sources for every matcher by a deterministic score
   built from `p_star` (after the §5.0 adjustments) and brightness, not by brightness alone. In SSA
   frames the brightest detections are disproportionately non-stars. Basis: Chum & Matas 2005.
3. **§6.4, §6.3a: derive roll, scale and parity priors from the trail model.** In target-tracking
   mode, the common star-trail vector together with the predicted apparent target rate (TLE plus
   observer) and the exposure time fixes roll (with sign, hence parity) and scale. Treat this as a
   derived prior with σ, recorded in the provenance with its origin. It turns many "hinted" solves
   into the 1-star-seed case of item 1. This is Fayde's own inference: the literature (Houtz &
   Frueh 2021) uses common trail geometry only as a filter.
4. **§5.4: local optimization before acceptance, and an NFA report.** Score hypotheses on a
   weighted refit (similarity → TAN, plus rate terms when there are endpoints), not on the raw
   seed transform. Count the hypotheses tested and report an NFA-style figure alongside the
   log-odds, optionally requiring NFA < ε in blind mode. Basis: Chum et al. 2003; a-contrario
   framework.
5. **§5.5: one mixture likelihood shared by the verifier and the fitter.** Replace "re-match within
   kσ" with IRLS over the same Gaussian-plus-distractor mixture the verifier uses, with per-source
   outlier priors from `p_star`. Fit and report a global error-scale factor *s*. Choose the SIP
   order by BIC or leave-one-out instead of an F-test. Basis: Hogg, Bovy & Lang 2010; MAGSAC (the
   idea only).
6. **§5.5, §6.5, §9: persistent instrument distortion.** Add `fayde refit --calibrate`, a joint fit
   over many solved frames for a persistent SIP stored in the `InstrumentProfile`. Per-frame fits
   then solve only the linear terms plus at most a low-order correction. This avoids overfitting
   SIP order 3 on 30–60 partly trailed stars per frame. Basis: Bernstein et al. 2017; SCAMP.
7. **§5.1, §5.2, §11 question 2: epoch policy.** Build the index at a fixed epoch and exclude
   high proper-motion stars from patterns. Propagate at query time with the Butkevich & Lindegren
   2014 formulae, including covariance, for verification and fitting. UCAC4 errors at 2026 are
   0.02–0.1 px at 3″/px, so Gaia DR3 is needed only for accuracy goals below about 0.1″. In §6.6,
   note that aberration (ERFA `ab`) affects only CRVAL and scale at 1° FOV, never solve rate.
8. **§9: public benchmarks.** Add the ESA Kelvins "Star Trackers: First Contact" scenes (10k
   scenes, up to 80 false objects each; a Tetra entry scored within 0.4% of the winner) as a public
   regression test for the matchers, after checking their terms. Follow Rijlaarsdam et al. 2020's
   evaluation considerations, and use LOST's Pyramid and voting implementations as independent
   baselines on the synthetic trailed-star fixtures.

**Explicitly not recommended:**

- learned star-ID;
- QUEST/ESOQ (SVD plus weighted TAN is enough);
- TPV or Zernike output;
- MAGSAC as an algorithm;
- grid and Liebe-triangle matchers;
- non-dimensional star-ID;
- rolling-shutter models in v1.
