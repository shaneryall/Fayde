# Fayde plate solver: literature review

Companion to `docs/design/plate-solver.md` (the "design doc"; section numbers like §5.3 refer to
it, and "entry 3.1" to this survey's own entries). Fayde's design: pair seeds when a pointing and a
tight scale are given, a 5-ratio pattern hash otherwise, across indexes for different field widths,
astrometry.net-style log-odds verification, a weighted TAN+SIP fit, and UCAC4 or Tycho-2.

Sources were found by web search, then papers, abstracts and project pages were read. Each entry
says how much was read: **[full]**, the full text or the relevant sections of it; **[abstract]**,
the abstract or publisher summary only; **[secondary]**, known only through citing papers, surveys
or indexes. In sections 1 to 6, most entries give the citation with this marker, then the idea,
the evidence and a paragraph starting "Verdict:", in that order.

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
| 1.5 | SP-Search / recursive mode (Samaan, Mortari & Junkins 2005) | §5.6 tries a neighbor's predicted WCS, last, when the search finds nothing | §5.6 |
| 1.6 | Multi-Poles Algorithm (Schiattarella et al. 2017, 2020) | CONSIDER: evidence that magnitude filtering helps under heavy false-object load | §5.4 |
| 1.7 | Non-dimensional star-ID (Leake, Arnas & Mortari 2020) | SKIP (focal-length drift is not a problem here; the scale prior handles it) | — |
| 1.8 | Reverse-attitude statistics (Dong et al. 2024) | CONSIDER: attitude-space voting as the hinted fallback (same family as 1.3) | §5.3 |
| 1.9 | Surveys and benchmarks (Spratling & Mortari 2009; Rijlaarsdam et al. 2020; ESA Kelvins 2017) | CONSIDER the Kelvins dataset and Rijlaarsdam's protocol as public regression tests | §9 |
| 2.1 | Houtz & Frueh 2021 (polar-grid plate solving of trailed stars) | **ADOPT (idea)**: common trail geometry as a star/non-star cue (in detection); CONSIDER it as a roll and scale prior (not built) | streak-detection §4.1, §4.4 |
| 2.2 | Sease, Koglin & Flewelling 2013 (attitude plus rate from streak endpoints) | CONSIDER: a WCS with pointing and roll rates (Fayde fits a static WCS) | §5.5 |
| 2.3 | AURORAS "dynamic astrometric calibration" (Krantz et al. 2025) | CONSIDER a time-dependent WCS, with a quadratic term (Fayde fits a static WCS) | §5.5 |
| 2.4 | Vereš et al. 2012 trail fitting | ADOPT the trail fit, in Fayde's detector; the solver takes one isotropic error per source | streak-detection §4.3 |
| 2.5 | Rolling-shutter streak ID (Schiattarella et al. 2020) | SKIP (global shutter) | — |
| 2.6 | Astrometry24.NET (Sybilska et al. 2019) | SKIP (closed service; no algorithmic detail to take) | — |
| 3.1 | RANSAC minimal-sample argument (Fischler & Bolles 1981) | **ADOPT**: the reason to prefer 2-star seeds in hinted mode | §5.3 |
| 3.2 | PROSAC (Chum & Matas 2005) | **ADOPT**: order hypotheses by `p_star` × brightness | §5.3, §6.3 |
| 3.3 | LO-RANSAC (Chum, Matas & Kittler 2003) | **ADOPT**: refit the inliers before the final log-odds decision | §5.4, §5.5 |
| 3.4 | MAGSAC / MAGSAC++ (Barath et al. 2019/2020) | SKIP the algorithm; Fayde fits and reports a global error-scale factor | §5.5 |
| 3.5 | A contrario / NFA framing | Done: Fayde reports the number of false alarms (NFA) next to the log-odds | §5.4 |
| 3.6 | Mixture-model fitting (Hogg, Bovy & Lang 2010) | **ADOPT**: one inlier/outlier model shared by the verifier and the fitter | §5.4, §5.5 |
| 3.7 | Photometric terms (Clouse & Padgett 2000; MPA) | CONSIDER (not built), with heavy-tailed magnitude residuals | §5.4 |
| 4.1 | Wahba solvers (QUEST, ESOQ/ESOQ2, SVD; Markley & Mortari 2000) | SKIP QUEST/ESOQ and SVD; a closed-form similarity for 2-star seeds | §5.5 |
| 4.2 | SIP vs TPV vs Zernike (Shupe et al. 2005; Calabretta & Greisen 2002; SCAMP) | Keep SIP; SKIP TPV and Zernike | §5.5 |
| 4.3 | Persistent instrument distortion plus per-exposure terms (Bernstein et al. 2017; SCAMP) | CONSIDER: fit distortion across many frames (Fayde fits it per frame) | §5.5 |
| 4.4 | Model selection (F-test vs BIC vs cross-validation) | ADOPT BIC, checked on held-out cells of the frame | §5.5 |
| 5.1 | HEALPix uniformization (Górski et al. 2005; astrometry.net) | Keep | §5.2 |
| 5.2 | Fibonacci lattice (González 2010) / spherical spiral (Bauer 2000) | Keep for pattern fields (via cedar-solve) | §5.2 |
| 5.3 | Rigorous epoch propagation (Butkevich & Lindegren 2014) | **ADOPT** the zero-parallax case; errors grow per axis, with no covariance | §5.1 |
| 5.4 | UCAC4 error budget (Zacharias et al. 2013) | ADOPT: fixed-epoch index, query-time propagation | §5.1, §5.2 |
| 5.5 | Relativistic aberration (Klioner 2003; SOFA/ERFA `ab`) | CONSIDER the SOFA/ERFA formulation (not built); aberration moves the field center and changes the scale, nothing else | §6.6 |
| 6.1 | Learned star-ID (SIGNN, WACV 2025; ICPRAM 2023; RPNet) | **SKIP** (not credible for production verification) | — |

---

## 1. Star identification / lost-in-space (star-tracker literature)

Star-tracker algorithms are tuned for wide FOVs (8–20°), bright stars, tiny onboard catalogs and no
prior. Fayde's hinted path is different: it has a pointing prior, an often tight scale prior,
sometimes a roll prior, and a narrow field crowded with false detections. What carries over is
mostly about *robustness to false stars*, and the main lesson is simple: **the fewer true stars a
hypothesis needs, the more robust the matcher.**

### 1.1 Pyramid and the k-vector

Mortari, Samaan, Bruccoleri & Junkins, "The Pyramid star identification technique", *Navigation*
51(3), 2004. Mortari, "Search-less algorithm for star pattern recognition", *J. Astronautical Sci.*
45(2), 1997 ([Springer](https://link.springer.com/doi/10.1007/BF03546375)). Mortari & Neta,
"k-vector range searching techniques", AAS 00-128, 2000. **[secondary]** (described consistently
in [Leake et al. 2020](https://arxiv.org/abs/2003.13736) and in surveys; originals not read.)

The k-vector is a precomputed index over the *sorted* angular distances of all catalog star pairs.
A range query ("every pair with separation in [d−ε, d+ε]") then costs O(1) plus the number of
results, with no comparisons. Pyramid takes image star triples, looks each edge up in the k-vector,
intersects the candidate sets to identify a triangle, then confirms it with a fourth "pyramid"
star. A triangle that the fourth star cannot confirm is rejected, which is how false stars get
discarded.

Pyramid is widely called the nominal lost-in-space algorithm. In a six-algorithm comparison (cited
in the search results; not read in full) Pyramid had the highest average accuracy with false stars
and was the fastest. Its known weakness is a need for accurate focal length
([Leake et al. 2020](https://arxiv.org/abs/2003.13736)). Fayde's scale prior covers that. The
winner of the ESA Kelvins star-ID challenge ("Super k-ID") was built on k-vector techniques
([results](https://kelvins.esa.int/star-trackers-first-contact/results/), **[full]**).

Verdict: ADOPT the idea, not the algorithm (§5.3). In hinted mode, build a sorted pair list (a
k-vector) of the catalog stars inside the hint cone, at query time. That is a few hundred stars
and about 10⁵ pairs, built in well under a millisecond. Then seed hypotheses from **image star
pairs**. With a scale prior and known parity, two stars fix a full similarity transform. See entry
3.1 for why 2-star seeds beat 4-star patterns when many sources are false. Fayde tries the pair
seeds before its pattern hash (§5.2a).

### 1.2 Grid algorithm, and its Bayesian and polar variants

Padgett & Kreutz-Delgado, "A grid algorithm for autonomous star identification", *IEEE TAES*
33(1), 1997. Clouse & Padgett, "Small field-of-view star identification using Bayesian decision
theory", *IEEE TAES* 36(3), 773–783, 2000 (doi:10.1109/7.869495). Lee et al. 2007 (polar grid).
**[abstract]**

The grid algorithm rotates the field so that a pivot star's nearest neighbor lies on the x-axis,
rasterizes the neighbors onto a coarse grid, and compares bit patterns. Clouse & Padgett add a
Bayesian, adaptive acceptance threshold aimed at a 2° FOV.

Clouse & Padgett report 96.0% correct with <0.3% false positives at 0.5 px noise and 0.8 mag
noise. Grid methods depend on the axis-defining neighbor being real and detected, which false or
missing stars break. Houtz & Frueh (entry 2.1) say this explicitly for images of trailed stars.

Verdict: SKIP as a matcher. CONSIDER Clouse & Padgett only as further support for likelihood-ratio
acceptance with magnitude information (entry 3.7), which Fayde does not use (§5.4).

### 1.3 Geometric voting

Kolomenkin, Pollak, Shimshoni & Lindenbaum, "Geometric voting algorithm for star trackers",
*IEEE TAES* 44(2), 441–456, 2008, doi:10.1109/TAES.2008.4560198. **[abstract]**

Every catalog pair whose angular distance matches an image pair votes for both identities. Each
image star takes the identity with the most votes, and a second voting round checks the result.
The method uses geometry only and draws on the whole image rather than a chosen subset.

The authors report robustness to non-star objects. The LOST open-source tracker
([docs](https://uwcubesat.github.io/lost/classlost_1_1GeometricVotingStarIdAlgorithm.html))
notes "no statistical guarantees" on the wrong-answer rate.

Verdict: CONSIDER as a hinted fallback (§5.3). Inside a hint cone with a tight scale prior, pair
voting is cheap and degrades gracefully when the brightest sources are false. Fayde has no voting
stage (§5.3); its log-odds test, with the hypothesis-count correction, would gate the votes.

### 1.4 Liebe triangle

Liebe, "Pattern recognition of star constellations for spacecraft applications", *IEEE AES
Magazine* 8(1), 31–39, 1993, doi:10.1109/62.180383. **[secondary]**

A star plus its two nearest neighbors gives a triangle feature (two distances and an angle).

The ICPRAM 2023 related-work section and others describe it as fast, but it "suffered limitations
with the introduction of false stars". A single false nearest neighbor breaks it.

Verdict: SKIP.

### 1.5 SP-Search and recursive mode

Samaan, Mortari & Junkins, "Recursive mode star identification algorithms",
*IEEE TAES* 41(4), 2005. **[secondary]**

Given the previous attitude, recursive mode predicts which catalog stars fall in the FOV and
matches them locally, with no pattern search.

Verdict: §5.6 sequence mode is a narrower form: one WCS predicted from a solved neighbor, tried
last, when the search finds nothing.

### 1.6 Multi-Poles Algorithm (MPA)

Schiattarella, Spiller & Curti, "A novel star identification technique robust to high presence of
false objects: the Multi-Poles Algorithm", *Adv. Space Res.* 59(8), 2133–2147, 2017. Schiattarella,
Spiller & Curti, "Star identification robust to angular rates and false objects with rolling shutter
compensation", *Acta Astronautica*, 2020
([preprint](https://iris.uniroma1.it/bitstream/11573/1349415/4/Schiattarella_preprint_star-identification_2020.pdf),
**[full]** for the 2020 paper, which summarizes the 2017 one).

Several "pole" stars are each tried as the center of a pair-distance pattern, in three phases:
acceptance, verification, confirmation.

The authors report 100% success with up to 185 false objects against 33 cataloged stars, falling
to 66% at 486 false objects against 31 stars (stable attitude). MPA came second in the ESA Kelvins
challenge. The ESA-contest version's main improvement was **magnitude filtering**, which made the
confirmation phase unnecessary.

Verdict: CONSIDER, as evidence. The magnitude-filtering result is independent support for a
photometric term in the verifier (entry 3.7). The algorithm itself adds nothing beyond 1.1 and 1.3.

### 1.7 Non-dimensional star-ID

Leake, Arnas & Mortari, "Non-dimensional star-identification", *Sensors* 20(9):2697, 2020
([arXiv:2003.13736](https://arxiv.org/abs/2003.13736)). **[abstract]**

The method uses dihedral angles of spherical star triangles, which do not depend on focal length,
so identification survives a miscalibrated focal length or optical axis.

Verdict: SKIP. Fayde estimates scale itself, and its 5-ratio edge patterns are already
scale-invariant (edges over the longest).

### 1.8 Reverse attitude statistics

Dong, Wang, Wang & Wang, "Reverse attitude statistics based star map identification method",
arXiv:2410.23758, 2024 (submitted to *IEEE Sensors J.*)
([arXiv](https://arxiv.org/abs/2410.23758)). **[abstract]**

The method builds candidate attitudes from many matched pairs (with spatial-hash indexing), then
takes the mode of the attitude distribution. Matching and attitude solving happen together. It is
aimed at missing stars, false stars and *star trailing* in near-space conditions.

The authors claim a >14.3% better identification rate and >28.5% faster solving than the best
existing methods, with field and on-orbit tests. The baselines were not visible in the abstract.

Verdict: CONSIDER, as the same family as 1.3: attitude-space voting. If Fayde implements the hinted
fallback, voting in the (roll, offset) space of 2-star seeds is the natural form.

### 1.9 Surveys and benchmarks

Spratling & Mortari, "A survey on star identification algorithms", *Algorithms* 2(1), 2009.
Rijlaarsdam, Yous, Byrne, Oddenino, Furano & Moloney, "A survey of lost-in-space star
identification algorithms since 2009", *Sensors* 20(9):2579, 2020 (**[abstract]**; the MDPI and
PMC full text was blocked to automated fetch). Toloei et al., "A comparative analysis of star
identification algorithms", *Astrophys. Space Sci.*, 2020
([Springer](https://link.springer.com/article/10.1007/s10509-020-03775-9), **[abstract]**). ESA
Kelvins "Star Trackers: First Contact" (2016–2017)
([data](https://kelvins.esa.int/star-trackers-first-contact/data/),
[results](https://kelvins.esa.int/star-trackers-first-contact/results/), **[full]**).

Rijlaarsdam et al. find that published comparisons "can lead to inconsistent conclusions" and
give considerations for simulation-based evaluation. Kelvins is the only shared public benchmark
found: 10,000 scenes, 10° FOV, 1920×1440, up to **80 false artifacts per scene**, centroid noise
under 1 px, magnitude noise, Hipparcos. Its final accuracy scores, out of a perfect 9921, are Super
k-ID 9708, MPA 9703, **"Tetra" (the original hash-table approach behind tetra3) 9672**, then
others from 9605 down. So a hash-table matcher is within about 0.4% of the best on a
false-star-heavy benchmark.

Verdict: CONSIDER. The Kelvins scenes are synthetic, public and false-star heavy, which suits a
reproducible regression test of the pattern-hash and seed matchers. LOST
([GitHub](https://github.com/uwcubesat/lost)) implements Pyramid, geometric voting and tetra-style
matchers, an independent baseline.

---

## 2. Matching trailed stars

The literature is thin. No paper found does blind or hinted plate solving of *trailed* stars with
per-source covariance, endpoints and a time-dependent WCS together. The published SSA pipelines
mostly either stack the frame so stars become static or feed trail centroids to astrometry.net.
For each straight star trail, Fayde's solver receives a center and the larger of its two position
errors (streak-detection §7); its input has no trail endpoints, per-source covariance matrix or
time term. The useful contributions are narrower.

### 2.1 Houtz & Frueh: streak detection and plate solving

Houtz & Frueh, "Efficient streak detection and plate-solving techniques", Proc. 8th European Conf.
on Space Debris, ESA SDO, 2021
([PDF](https://conference.sdo.esoc.esa.int/proceedings/sdc8/paper/76/SDC8-paper76.pdf),
**[full]**).

Streaks are detected with cumulative-sum line integration. For plate solving, only streaks whose
length and orientation match the frame's most common values are kept as stars. A polar-grid
pattern is built around the brightest star, with the second-brightest streak as the axis, using
spherical-triangle radii and angles. It is matched by rotational correlation against catalog
patterns, centered on Tycho-2 stars near the a-priori pointing. All stars are then identified from
the radius/angle residuals, and QUEST gives the attitude at mid-exposure.

On 15 frames from the POGS 35 cm telescope, on average 11.7 trailed stars were identified per
frame, with star-vector errors under 1″. The sample is small and has no false-solve statistics.

Verdict: ADOPT the idea, not the matcher. This is the closest published analog to Fayde's case. It
confirms that (a) trail length and direction consistency is the strongest star/non-star cue
(Fayde's detector uses it: the census and the trail-against-point comparison set `p_star`,
streak-detection §4.1 and §4.4), and (b) a hinted, pointing-centered search is the practical
approach. One step further, not built: in target-tracking mode the common star trail is the sky's
apparent motion relative to the tracked target, projected onto the detector. Given the commanded
tracking rate and the exposure time, the trail's **direction could constrain roll** and its
**length the pixel scale**. This is an inference of this survey, not a published method.

### 2.2 Sease et al.: attitude and rate from streak endpoints

Sease, Koglin & Flewelling, "Long integration star tracker image processing for combined
attitude–attitude rate estimation", *Proc. SPIE* 8739, 2013, doi:10.1117/12.2020004
([UCF record](https://stars.library.ucf.edu/scopus2010/6089)). **[abstract]**, plus a
catalog-free follow-up and a UCF thesis
([core.ac.uk](https://core.ac.uk/download/236295411.pdf), **[secondary]**).

The method treats streak endpoints as corner features. It groups the endpoints into "start" and
"end" sets using a cross-product test of the rotation sense, then runs ordinary star-ID and
attitude determination on each set. Two attitudes from one image give the attitude rate.

In Monte Carlo runs, endpoint accuracy is better than 0.1 px when more than 20° from the rotation
axis, for streaks up to 40% of the FOV. The thesis reports rates within 2% in simulation.

Verdict: CONSIDER. A joint fit of the mid-exposure WCS plus pointing and roll rates would be the
least-squares generalization of this: one fit, with covariances, instead of two independent
solves. Fayde fits a static WCS (§5.5). Note the asymmetry: endpoint error *along* the trail is
usually much larger than across it (entry 2.4), so endpoint covariances must be anisotropic.

### 2.3 AURORAS: time-dependent astrometric calibration

Krantz, Bloch, Serda, Flechtner & McMahon-Crabtree, "AURORAS: orbit determination with just one
look", AMOS 2025
([PDF](https://amostech.com/TechnicalPapers/2025/Astrodynamics/Krantz.pdf), **[full]**, relevant
sections).

Event-camera or fast-frame data are shift-stacked so the stars become static and can be plate
solved against Gaia. Then a "dynamic astrometric calibration" fits an (x, y, t) → (RA, Dec, t)
transform. Residuals before the time-dependent fit show offsets and curvature from the moving
tangent plane plus mount jitter; after it they are flat.

The paper shows qualitative residual plots and gives no numbers on astrometric accuracy.

Verdict: CONSIDER a time-dependent WCS, with an optional **quadratic time term** (pointing
acceleration) for non-linear trails, as in long exposures tracking fast LEO targets. Fayde fits a
static WCS (§5.5); curved star paths are handled in detection (streak-detection §8).

### 2.4 Trail fitting for astrometry (producer side)

Vereš, Jedicke, Denneau, Wainscoat, Holman & Lin, "Improved asteroid astrometry and photometry
with trail fitting", *PASP* 124, 1197, 2012 ([arXiv:1209.6106](https://arxiv.org/abs/1209.6106)).
**[abstract]**

A trail is modeled as a Gaussian PSF convolved with a line segment, and its center, length, angle
and flux are fitted.

The fit gives about 3× better astrometry and up to 2 mag better photometry than PSF or aperture
methods for trails longer than about 10 px (Pan-STARRS1).

Verdict: ADOPT the trail-model fit. Fayde's detector measures a shared trail length, angle and
width, fits each trail with them held fixed, and reports errors along and across it
(streak-detection §4.3). The solver takes one isotropic error per source (§6.3); errors that are
narrow across the trail and wide along it, and endpoints, would need a change to the input contract
and the verifier.

### 2.5 Rolling-shutter streaks

Schiattarella et al. 2020 (entry 1.6), **[full]**. Curti et al. (rate estimation from streaks up
to 10–25°/s; [Sapienza](https://iris.uniroma1.it/handle/11573/852340), **[abstract]**).

A rolling shutter distorts streak length by about 25% (Curti et al.). A first-order compensation
gives reliable identification up to 5°/s.

Verdict: SKIP. The frames this survey targets come from global-shutter sensors.

### 2.6 Astrometry24.NET

Sybilska et al., "Astrometry24.NET – precise astrometry for SST and NEO", Proc. 1st NEO and Debris
Detection Conf., ESA, 2019
([PDF](https://conference.sdo.esoc.esa.int/proceedings/neosst1/paper/147/NEOSST1-paper147.pdf),
**[full]** for the introduction and architecture).

Astrometry24.NET is a commercial cloud service that fits streak profiles and solves frames tracked
at non-sidereal rates. The paper describes the architecture and accuracy campaigns but gives
little detail on the matching algorithm.

Verdict: SKIP. Nothing to port. It could serve as an external accuracy cross-check if real frames
were ever run through it.

---

## 3. Verification, false-match probability and robust sampling

### 3.1 The minimal-sample argument (RANSAC)

Fischler & Bolles, "Random sample consensus", *Comm. ACM* 24(6), 1981.
**[secondary; standard result]**

If a fraction *w* of sources are real catalog stars, a hypothesis built from *k* sources is
all-true with probability about *wᵏ*. At *w* = 0.5, a 4-star pattern (tetra3 or astrometry.net)
succeeds 6.3% of the time, a 2-star seed 25%, and a 1-star seed 50%. At *w* = 0.3 the figures are
0.8%, 9% and 30%. For a source list of 30–60 entries with many false detections, this ratio
matters more than any detail of hash design.

Verdict: ADOPT (§5.3). When priors fix scale (and parity), the hinted path should seed from 2-star
pairs (entry 1.1). Single-star seeds, when roll is also known, would be the next step: each
image–catalog pairing proposes a translation, and the translations vote (not built). The pattern
hash remains the blind and weak-prior path.

### 3.2 PROSAC

Chum & Matas, "Matching with PROSAC – progressive sample consensus", CVPR 2005,
doi:10.1109/CVPR.2005.221
([PDF](https://dspace.cvut.cz/bitstream/10467/9496/1/2005-Matching-with-PROSAC-progressive-sample-consensus.pdf)).
**[abstract]**, plus reviews.

PROSAC samples hypotheses from a progressively growing set of the top-ranked correspondences
rather than uniformly, which is still correct in the limit. The gain depends entirely on how well
the ranking predicts correctness.

It needs up to more than 100× fewer samples than RANSAC on CV benchmarks (in one example, 9
samples against about 10⁵), and gives no advantage when the ordering is random.

Verdict: ADOPT (§5.3, §6.3). tetra3/cedar's brightest-first, breadth-first enumeration is already
a PROSAC-like ordering, but it ranks by brightness alone. Fayde ranks by `p_star` × flux (§5.0),
and its detector's `p_star` carries the trail consistency, which predicts "real star" much better
than flux alone (the brightest sources in SSA frames are often the target, cosmic rays or blends).
Fayde feeds all matchers in that order (`p_star` × input-order rank when there is no flux). The
ordering must be deterministic (requirement 6).

### 3.3 LO-RANSAC (local optimization)

Chum, Matas & Kittler, "Locally optimized RANSAC", DAGM 2003, LNCS 2781.
**[secondary; standard]**

When a hypothesis looks promising, LO-RANSAC refits the model on its inliers (and iterates)
before scoring it. Minimal-sample models are noisy, so scoring the raw model underestimates
support.

Verdict: ADOPT (§5.4/§5.5). This matters most for trailed stars and for 2-star seeds. A similarity
transform from two stars, ignoring trail rate and distortion, can leave a true match just below
threshold. Fayde's flow (§5.4):

1. A cheap pre-filter.
2. If it passes, a weighted refit on its matches (similarity → TAN).
3. Re-match, then compute the log-odds on the **refined** model.

astrometry.net does a version of this ("tweak" during verification).

### 3.4 MAGSAC / MAGSAC++

Barath, Matas & Noskova, "MAGSAC: marginalizing sample consensus", CVPR 2019; Barath, Noskova,
Ivashechkin & Matas, "MAGSAC++, a fast, reliable and accurate robust estimator", CVPR 2020
([arXiv:1912.05909](https://arxiv.org/abs/1912.05909)). **[abstract]**

MAGSAC removes the hard inlier threshold by marginalizing the model score over the noise scale σ,
up to a loose upper bound. MAGSAC++ fits with weighted least squares (IRLS) in place of a hard
inlier set.

Verdict: SKIP the algorithm. Fayde has per-source position errors, so the log-odds already uses an
explicit noise model, which is what MAGSAC approximates when one is missing. The underlying idea
is built: producers' errors are often miscalibrated, so Fayde fits a single global error-scale
factor *s* (covariance × *s*²) and reports it (§5.5), so a detector that under-reports its errors
is visible. It has no prior on *s* and does not marginalize over it.

### 3.5 A contrario / number of false alarms (NFA)

The a-contrario detection framework (Desolneux, Moisan & Morel; Abergel & Moisan on
trajectories). **[secondary]** No paper was found that applies NFA to star identification.

A detection is accepted when (number of tests) × P(an event at least this good under the null
hypothesis) < ε. The threshold is calibrated against the *search size*, not set per hypothesis.

Verdict: done (§5.4). astrometry.net's 10⁹ odds is a fixed per-hypothesis threshold. Blind solves
test many more hypotheses than hinted ones, so the family-wise false-accept rate differs between
them. Fayde counts the hypotheses tested, reports `NFA = n_tested × exp(−log_odds)` next to the
log-odds, and in every search mode compares `log_odds − ln(n_tested)` with the log-odds threshold,
`ln(10⁹)` by default. This also gives a calling pipeline a single number to record.

### 3.6 Mixture-model fitting

Hogg, Bovy & Lang, "Data analysis recipes: fitting a model to data", arXiv:1008.4686, 2010.
**[secondary; well-known]**

Each data point is modeled as a mixture: inlier (Gaussian with its covariance) or outlier (broad
background). The mixture is fitted or marginalized in place of hard clipping at kσ.

Verdict: ADOPT (§5.4/§5.5). This is the *same* likelihood the astrometry.net verifier uses (a
Gaussian match against a uniform distractor). Fayde's final fit uses it too, as IRLS on the
mixture weights with per-source outlier priors from `p_star`; the local refit before scoring still
re-matches within shrinking radii. One noise model, used consistently, is easier to reason about
and gives better-calibrated covariances.

### 3.7 Photometric consistency in verification

Clouse & Padgett (1.2) use 0.8 mag noise inside a Bayesian decision. The MPA authors (1.6) found
magnitude filtering was the main improvement for their contest entry. A brightness-referenced
method ([Sensors 2014](https://dx.doi.org/10.3390/s141018498), **[abstract]**) reports speed gains
from brightness ordering.

Verdict: CONSIDER (not built: the verifier does not compare brightnesses, §5.4). The literature and
this regime suggest the following:

- Fit the zero point robustly.
- Use a **heavy-tailed** (Student-t) magnitude residual, because UCAC4 `magm` against an
  unfiltered or wide band, trail photometry, saturation and blends all produce outliers.
- Weight the magnitude term well below the positional term until it has been calibrated on a
  corpus of real frames.

The magnitude term is a second discriminator, not a primary one.

---

## 4. WCS fitting

### 4.1 Wahba's problem: QUEST, ESOQ, SVD

Shuster & Oh, "Three-axis attitude determination from vector observations", *J. Guid. Control*
4(1), 1981 (QUEST). Mortari, "ESOQ: a closed-form solution to the Wahba problem",
*J. Astronaut. Sci.* 45(2), 1997. Markley & Mortari, "Quaternion attitude estimation using vector
observations", *J. Astronaut. Sci.* 48(2–3), 359, 2000, doi:10.1007/BF03546284 (preprint, under
the title "New developments in quaternion estimation from vector observations":
[NTRS](https://ntrs.nasa.gov/citations/20000034107)). **[abstract]** Umeyama, "Least-squares
estimation of transformation parameters between two point patterns", *IEEE TPAMI* 13(4), 1991.
**[secondary; standard]**

Davenport's q-method and SVD are the most robust; QUEST, ESOQ and ESOQ2 are faster and adequate
for star sensors. Robustness only matters when measurement accuracies differ widely (Markley &
Mortari).

Verdict: SKIP QUEST/ESOQ and SVD. At 30–60 stars, speed is irrelevant, and none of these solvers
handles anisotropic weights. Fayde's pair seeds use the closed-form two-point similarity in the
tangent plane (entry 3.1), once for each allowed parity; pattern seeds fit a similarity to all four
stars. The one that verifies is kept. The weighted TAN least-squares fit follows, with one weight
per matched pair.

### 4.2 Distortion representation: SIP vs TPV vs Zernike

Shupe et al., "The SIP convention for representing distortion in FITS image headers", ASP Conf.
Ser. 347, 2005. Calabretta & Greisen, "Representations of celestial coordinates in FITS", *A&A*
395, 2002 (TPV is a registered convention, not part of the standard). Bertin, "Automatic
astrometric and photometric calibration with SCAMP", ASP Conf. Ser. 351, 2006. Shupe et al., SPIE
2012 (`sip2pv`/`pv2sip` conversion). **[secondary]**, from survey pipeline papers
([PTF](https://arxiv.org/pdf/1404.1953), [SAGES](https://arxiv.org/pdf/2306.15611)).

SIP (pixel-space polynomials) and TPV (polynomials in intermediate world coordinates) are
interconvertible to good precision. SCAMP writes TPV, and pipelines routinely convert it to SIP
because most tools read SIP. Zernike is useful only for wide-field optics with a circular pupil,
and has no FITS convention.

Verdict: keep SIP; SKIP TPV and Zernike (§5.5, "Won't" confirmed).

### 4.3 Persistent instrument distortion with per-exposure terms

Bernstein et al., "Astrometric calibration and performance of the Dark Energy Camera", *PASP*
129:074503, 2017 ([arXiv:1703.01679](https://arxiv.org/abs/1703.01679), **[abstract]**). SCAMP
(Bertin 2006), which fits distortion across many exposures at once.

The astrometric model is split into a *persistent* instrument component (optics, detector), fitted
once from many exposures, and a *per-exposure* component, kept low-order. DECam reaches 3–6 mas
per exposure this way.

Verdict: CONSIDER. With 30–60 stars per frame, many of them trailed with large along-trail errors, a
per-frame SIP order 3 (14 free coefficients on top of 6 linear ones) is at serious risk of
overfitting. A persistent SIP fitted across solved frames could live in the instrument profile,
leaving each frame only the linear terms and a low-order correction. Fayde fits SIP per frame and
guards the order with held-out cells (§5.5); its profile holds no distortion.

### 4.4 Model selection and covariance

These are standard statistical results. The F-test assumes nested models with Gaussian, correctly
scaled errors. With miscalibrated covariances it over-selects. BIC penalizes parameters by
log *N*. Leave-one-out cross-validation on the matched stars is cheap at N ≤ 60 and robust to
error-scale problems.

Verdict: ADOPT (§5.5). Choose the SIP order by BIC and check it on held-out cells (Fayde's rule,
§5.5). Report the covariance scaled by the fitted *s*² (entry 3.4), with *s* alongside, as Fayde
does. Include the catalog position error at the epoch in W (the design already does this), so the
covariance is not optimistic when UCAC4 dominates the error.

---

## 5. Catalog and index design

### 5.1 HEALPix uniformization

Górski et al., "HEALPix", *ApJ* 622, 2005. Lang et al. 2010 (`uniformize-catalog`).

Verdict: keep (§5.2). Nothing better found.

### 5.2 Fibonacci lattice and spherical-spiral catalogs

González, "Measurement of areas on a sphere using Fibonacci and latitude–longitude lattices",
*Math. Geosci.* 42, 2010 ([arXiv:0912.4540](https://arxiv.org/abs/0912.4540), **[abstract]**).
Bauer, "Distribution of points on a sphere with application to star catalogs", *J. Guid. Control
Dyn.* 23(1), 2000 (**[secondary]**; the spherical spiral method). Variants: charged particles,
k-means clustering (GWKM: Somayehee, Nikkhah & Roshanian 2019, *J. Navig.* 72(4), 948,
doi:10.1017/S0373463318001029, **[abstract]**).

Fibonacci lattices give near-equal-area cells for any point count, with at least 40% lower RMS
error in area estimates than latitude–longitude grids. Spherical-spiral catalogs are the
star-tracker standard for uniform onboard catalogs.

Verdict: keep (§5.2). cedar-solve's Fibonacci lattice fields for pattern generation plus HEALPix
for star selection is a sound combination. No change.

### 5.3 Rigorous epoch propagation

Butkevich & Lindegren, "Rigorous treatment of barycentric stellar motion", *A&A* 570, A62, 2014
([arXiv:1407.4664](https://arxiv.org/abs/1407.4664), **[abstract]**). The same formulae appear in
the [Gaia DR3
documentation](https://gea.esac.esa.int/archive/documentation/GDR3/Data_processing/chap_cu3ast/sec_cu3ast_intro/ssec_cu3ast_intro_tansforms.html).

The paper gives closed-form propagation of all six astrometric parameters between epochs,
reversible, with a full **covariance** transformation.

Verdict: ADOPT as the reference (§5.1). With UCAC4 (no parallax or radial velocity) the equations
reduce to the zero-parallax case, which Fayde's file-backed catalogs use. Fayde grows each axis's
position error with the proper-motion error, with no covariance transformation, and the result
feeds the fit weights in §5.5.

### 5.4 UCAC4 error budget and epoch policy

Zacharias et al., "The fourth US Naval Observatory CCD Astrograph Catalog (UCAC4)", *AJ* 145:44,
2013 ([arXiv:1212.6182](https://arxiv.org/abs/1212.6182)). **[abstract]**

Formal position errors are 15–100 mas at the mean epoch, plus 10–20 mas systematic. Proper-motion
systematics are about 1–4 mas/yr; random errors are larger for faint stars. Propagated to 2026
(about 26 yr from 2000, somewhat less from UCAC4's central epochs), typical position errors are
therefore about 0.05–0.3″. At 3″/px that is **0.02–0.1 px**: irrelevant for matching and small
for the fit. Gaia DR3 is about 100× better and would only matter for accuracy claims finer than
about 0.1 px.

Verdict: ADOPT a two-level policy (§5.1, §5.2).

1. Build the index at a fixed epoch (for example 2025.0). For 99.9% or more of stars,
   |µ| < 100 mas/yr, so ±10 yr moves them less than 1″ (0.3 px), well within the pattern
   tolerance.
2. Exclude high proper-motion stars from the index (Fayde's cut is about 180 mas/yr per degree
   of field, §5.2). Query-time catalog reads still see them.
3. Propagate at query time, to the mid-exposure epoch, for verification and fitting.

Rebuilding the index for each campaign is then unnecessary. Gaia DR3 matters only for accuracy
requirements finer than about 0.1 px.

### 5.5 Aberration for spacecraft observers

Klioner, "A practical relativistic model for microarcsecond astrometry in space", *AJ* 125, 2003.
**[secondary]** SOFA/ERFA `iauAb`/`eraAb` (fully relativistic to about 1 µas). Gravity Probe B
context: orbital aberration about ±5.2″, annual about ±20.5″ (via search, **[secondary]**).

For a telescope in low Earth orbit, β = |v_annual + v_orbital|/c is at most about 1.25×10⁻⁴
(computed for this survey). That shifts the field center by up to about 26″. The *differential*
effect across a 1° field is a scale change of order β, at most about 0.2″ at the field edge,
about 0.07 px.

Verdict: CONSIDER ERFA's `ab` formulation (§6.6). At this FOV, aberration affects only CRVAL and the
plate scale, so matching and verification can safely ignore it, and observer state would be an
output-accuracy input, never a solve-rate input. Fayde takes no observer state and applies no
aberration.

### 5.6 Index compression

No result was found that improves on cedar-solve's layout (a 16-bit key check plus the longest
edge in each slot; Fayde stores the edge as `f32`) for Fayde's scale. Index size is covered in
§5.2.

---

## 6. Learned / ML star identification

Hepburn-Dickins, Jones, Edwards, Morgan & Bell, "SIGNN – star identification using graph neural
networks", WACV 2025 ([CVF
page](https://openaccess.thecvf.com/content/WACV2025/html/Hepburn-Dickins_SIGNN_-_Star_Identification_using_Graph_Neural_Networks_WACV_2025_paper.html),
**[abstract]**). Hepburn-Dickins & Edwards, "Noise robustness of data-driven star classification",
ICPRAM 2023 ([PDF](https://scitepress.org/publishedPapers/2023/118040/pdf/index.html), **[full]**).
RPNet (Xu et al. 2019), 1-D CNN classifiers (Yang et al. 2022), and Wang et al., *Sensors* 2021
(**[secondary]**).

SIGNN is a graph-attention network over a Hipparcos star graph with edges below 3°.

ICPRAM 2023 reports 98.5% accuracy on augmented samples, but **13.5% of stars misidentified** on
simulated Stellarium images. A 0.999 confidence cut removed 99.85% of wrong identifications but
discarded 39% of correct ones. There is no head-to-head comparison against Pyramid or tetra3. Most
learned methods are trained for wide-FOV Hipparcos star trackers. None addresses trailed stars,
sub-degree FOVs with UCAC4 density, or calibrated false-accept rates.

Verdict: SKIP. Not credible for production where a wrong WCS is worse than no WCS. Fayde's
correctness rests on a likelihood-ratio verifier with a hypothesis-count correction, which ML
methods do not provide. If anything, a learned *source classifier* belongs on the producer side,
feeding `p_star`, which the contract already supports.

---

## 7. Source and streak detection (2026-11-16)

Survey for detection in frames with point-like stars (`docs/design/detection.md`). Reference
implementations were read at the head of their default branches; **[full]**, **[abstract]** and
**[secondary]** as above, **[code]** for source code.

### 7.1 Star detectors in use

| Tool | Background / noise | Filter | Threshold | Shape tests | Centroid | License |
|---|---|---|---|---|---|---|
| astrometry.net `simplexy` (Lang et al. 2010) **[code]** | Median in 201 px boxes on a grid; noise from pixel pairs 5 px apart | Gaussian σ = 1 px, fixed | 8σ (4σ for 8-bit) | none; peaks split at a 5σ saddle | Parabolas on the 3 × 3 of the smoothed image | BSD-3 own code, GPLv3 as distributed |
| SExtractor / SEP (Bertin & Arnouts 1996; Barbary 2016) **[code]** | 64 px mesh, clipped mode `2.5·median − 1.5·mean`, 3 × 3 median, spline | 3 × 3 kernel "FWHM 2 px", fixed | 1.5σ on 5 connected pixels | elongation from second moments; multi-threshold deblending | Isophotal or Gaussian-windowed moments (3 to 5 iterations) | GPLv3 / LGPLv3 |
| DAOFIND, photutils (Stetson 1987) **[code]** | Zero-sum kernel fits a local sky | Truncated Gaussian of a **user-given** FWHM | 4 to 5σ | sharpness 0.2 to 1.0, roundness ±1 | Marginal Gaussian fits | BSD-3 (photutils) |
| tetra3 (ESA) **[code]** | 25 px box mean; global RMS (not robust) | none | 2σ, binary opening, area 5 to 100 px | optional axis ratio | Center of gravity of the region | Apache-2.0 |
| cedar-detect (smroid) **[code]** | Darkest of 15 strips | none; 7-pixel row gate, then rings around a blob; manual 2 to 8 × binning | 8σ typical | hot pixel if its four neighbors hold under half its excess; perimeter must be uniform | Center of gravity, core + 2 px | FSL-1.1 (source-available, not open source); 8-bit only |
| ASTAP (H. Kleijn) **[code]** | Histogram mode, global | none | ≥ 3.5σ pixel, aperture SNR > 10 | HFD limits, fill fraction, single-pixel reject | Center of gravity above 3σ | MPL-2.0 |
| Siril, PixInsight **[code]** | Global median + noise; Gaussian high-pass | Gaussian σ = 2 px fixed; none | 5σ; 3σ then SNR | roundness; kurtosis, coverage | Gaussian or Moffat fit (Levenberg–Marquardt) | GPLv3; PCL license |

Published speed: SExtractor "up to about 50 Mpixel/s"; cedar-detect "less than 10 ms" per
megapixel on a Raspberry Pi 4B. None for the others.

None of them measures the star width and adapts the filter to it, reports a per-star position
error to the matcher, models cosmic rays or hot pixels other than by minimum area or sharpness,
finds the center of a saturated star with a bleed trail, or traces streaks before solving
(astrometry.net removes sources in over-dense rows and columns only).

### 7.2 Centroid accuracy

- Cramér–Rao bound, Gaussian star: `√2·σ / SNR` per axis background-limited (0.60 FWHM / SNR),
  `σ / √N` source-limited (King 1983 **[full]**; Mendez, Silva & Lobos 2013 **[full]**; Mortensen
  et al. 2010 **[full]**, whose eq. 6 gives `(16/9)·σ²/N + 8π·σ⁴·b/N²` for a Gaussian-weighted
  centroid).
- A matched filter followed by a 3 × 3 quadratic fit reaches the bound for FWHM 2.8 to 5.6 px at
  SNR 5 to 100; a fixed kernel does not at the extremes; a 7 × 7 center of light is several times
  worse (Vakili & Hogg 2016 **[full]**).
- A parabola through three samples of a Gaussian peak is biased toward the pixel center; the
  same fit on the logarithms has no bias (Naidu & Fisher 1991 **[full]**).
- Undersampled stars (σ 0.66 px): moment methods have pixel-phase floors of 0.02 to 0.04 px; only
  a pixel-integrated profile fit goes lower (Fialho 2025, arXiv:2507.17560 **[full]**, preprint).

Verdict: ADOPT the matched filter with a log-quadratic 3 × 3 fit; derive the error per source.

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

Verdict: ADOPT the peak density as a threshold and null model, bounded by the pixel count.

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

Steger 1998 links Hessian line points into lines, and in trailed frames Houtz & Frueh 2021 (entry
2.1) and Sharma et al. 2023 (entry 8.4) detect trails before plate solving. For frames with
point-like stars, no surveyed plate solver traces streaks before solving.

Verdict: ADOPT ridge-like peaks as seeds and trace from them (Fayde's own tracer); SKIP Radon and
Hough (cost; they find streaks too faint to disturb the star list).

## 8. Trailed stars and faint streaks (2026-11-17)

Survey for one detector that handles trailed stars, streaks and faint streaks
(`docs/design/streak-detection.md`). Sources already in entries 2 and 7.4 are not repeated. Many
publisher PDFs could not be read; each entry says what was. Figures marked *derived* were worked
out during the survey and are not from a source.

### 8.1 The frame's common trail

| Source | Idea | Numbers | Verdict |
|---|---|---|---|
| Danescu & Turcu 2026, *Sensors* 26(5), 1628, doi:10.3390/s26051628 **[full]** | LoG filter, top 5% of pixels, connected components, second moments; components with eccentricity over 0.65 vote; the frame's trail is the median direction and median major-axis length; template = 1 px line ⊛ Gaussian | astrometry.net solved 39 of 40 frames after filtering with the template, 15 of 40 before; under 1 s for 1536 × 1024; fails when trails reach about 12% of the frame width | ADOPT: cheap, deterministic, the median ignores a satellite or cosmic rays |
| Lévesque 2011, AMOS poster **[abstract]** | Track-rate frames: trail direction from the Radon transform of the Fourier modulus; length next; an iterative matched filter; star trails measured before satellites are sought | none obtained | CONSIDER as the segmentation-free check |
| Li, Zhang, Yao & Sun 2013, *Opt. Precision Eng.* 21(10), 2656, doi:10.3788/OPE.20132110.2656 **[abstract]** | Direction from PCA of the binarized cepstrum; length from the minima of the spectrum's Radon profile along it | 0.14° and 0.27 px on natural images | CONSIDER for refining long trails |
| Cannon 1976, *IEEE TASSP* 24(1), 58; Moghaddam & Jamzad 2007, *EURASIP JASP* 068985; Oliveira, Figueiredo & Bioucas-Dias 2014, *IEEE TIP* 23(1), 466 **[abstract]** | Blur kernel from the cepstrum or the Radon transform of the spectrum | natural images only | SKIP as primary: noise-sensitive, not tested on star fields |
| Sun, Xing, You & Wei 2013, *Opt. Express* 21, 20096 **[abstract]**; Hou et al. 2021, *Math. Probl. Eng.* 8819277 **[abstract]** | Star-tracker smear direction and length from image derivatives, bispectrum or the spectrum's Radon profile | not obtained | CONSIDER |

*Derived*: the power spectrum of a field of identical trails is `|PSF̂|²·sinc²` times a speckle
factor from the star positions, so in an ideal frame the sinc's zeros survive any star layout.
Cropping, noise and trails that vary across the frame fill them in. For trails under about
2 FWHM the first zero falls where the PSF has already removed the power: the spectrum can only
help with long trails.

### 8.2 Aimed filters and the trailing loss

| Source | Idea | Numbers | Verdict |
|---|---|---|---|
| Nir, Zackay & Ofek 2018 (entry 7.4) **[full]** | Optimal streak detector: the PSF-broadened line; `(S/N) = ξ·√(L / (2√π·σ·B))` for flux ξ per unit length | | ADOPT the algebra |
| Jones, "Trailing losses for moving objects", Rubin SMTN-003 **[full]**; Jones et al. 2018, *Icarus* 303, 181, arXiv:1711.10621 **[full]** | Loss fits `Δm = 1.25 log₁₀(1 + a·x² / (1 + b·x))`, `x` = length / FWHM; point filter `a = 0.420, b = 0.003`, trail filter `a = 0.761, b = 1.162` (SMTN-003; Jones et al. 2018 give `a = 0.67, b = 1.16` for the trail filter) | at `x = 5`: 1.32 against 0.72 mag | ADOPT as an independent check of the gain |
| Soille, Breen & Jones 1996, *IEEE TPAMI* 18(5), 562, doi:10.1109/34.494646 **[abstract]** | Recursive filters along translates of a Bresenham line at any angle: a running sum costs a few operations per pixel whatever the length | | ADOPT for the box along the trail |
| Geusebroek, Smeulders & van de Weijer 2003, *IEEE TIP* 12(8), 938 **[abstract]** | Anisotropic Gaussian at any angle: one pass along x and one along a skewed axis, recursive | 512² in 40 ms on 2003 hardware | CONSIDER; an isotropic Gaussian after the box is already the matched filter |
| Cvrček & Šára 2021, 8th European Conf. on Space Debris, paper 80 **[abstract]** | Steerable basis filters for streaks of any angle and length | over 100 × faster than spatial convolution | SKIP: short streaks, which the tracer covers |

*Derived*, checked numerically: for long trails the aimed filter's response with an angle error
`Δθ` is `(√π / 2z)·erf(z)` of the matched one, `z = L·|Δθ| / 4σ`: 2% lost at `σ/L`, 7.7% at
`2σ/L`, 25% at `4σ/L`. A length off by a factor `k` keeps `1/√k`: 29% lost at `k = 2`, 16% at
`k = √2` (the worst case on a grid of factor-2 steps).

### 8.3 Fitting a trail and its error

| Source | Idea | Numbers | Verdict |
|---|---|---|---|
| Bouquillon et al. 2017, *A&A* 606, A27, arXiv:1707.01447 **[full]**, main sections | Erf model (Vereš et al. 2012 with a sign corrected); Cramér–Rao bounds along (U) and across (V) the trail | background-limited, `L ≥ 2 FWHM`: `Var_U ∝ L̄²`, `Var_V ∝ L̄² / (√π·L̄ − 1)`, `L̄ = L / 2σ`; source-limited: `Var_U ≈ 1.11·(σ²/F)·L̄`, `Var_V = σ²/F`; unweighted least squares about 30% above the bound across | ADOPT as the error model and the benchmark |
| Blake et al. 2026, "DebrisWatch II", *J. Astronaut. Sci.*, doi:10.1007/s40295-026-00602-1, arXiv:2605.27516 **[full]**, pipeline section | Gaussian across, a difference of arctangents with tilt terms along; candidates pre-filtered by expected angle and length; a trail is a blend when its fitted centroid departs from the detection's | median offset from about 2″ to under 0.3 px | CONSIDER the tilt and the blend test (not built; Fayde refits neighbors instead) |
| Kouprianov 2008, *Adv. Space Res.* 41, 1029 **[abstract]** | Trailed reference stars limit GEO astrometry; jitter and extinction along the trail are common to the frame, so one empirical profile is fitted to every star | | CONSIDER |
| Fraser et al. 2016, "TRIPPy", *AJ* 151, 158, arXiv:1604.00031 **[abstract]** | Moffat ⊛ line PSF with a look-up table; "pill" aperture | aperture corrections to about 10 mmag | CONSIDER the pill aperture for fluxes |
| Yan, Jiang & Zhang 2016, *Opt. Express* 24, 5961 **[abstract]**; Wei, Tan, Li & Zhang 2014, *Sensors* 14, 4914 **[full]**, partial | Star-tracker smeared-spot models; the peak pixel holds about `0.264/l` of the flux for trails over 9 px | | CONSIDER; shows why peak thresholds fail on trails |
| Zhang, Quan & Guo 2012, *Sensors* 12, 6712 **[full]**; Mu et al. 2019, *Sensors* 19, 4127; Wang et al. 2018, *Sensors* 18, 2662 | Blur kernel from a gyro, then deconvolution | | SKIP: need the rate as input |
| Kakarala & Hero 1992, *IEEE TPAMI* 14(7), 777 **[abstract]** | Cramér–Rao bounds for edge position; Canny's estimator about twice the bound at high signal | | CONSIDER: simple end estimators come within a factor of two |

*Derived* (Gaussian PSF of width `s`, pixel noise σ, flux λ per pixel of length; leading order for
long, well-sampled trails in white background noise): an end is located to
`σ_end ≥ 2√π·σ·s / λ`, a trail's offset across to `2π^{1/4}·σ·s^{3/2} / (λ√L)`, its angle to
`√12·σ_offset / L`. A whole trail's center, located from both ends, has `σ_end / √2`, consistent
with Bouquillon et al. 2017, App. A.3. Their scalings, the factor 1.11 included, hold at fixed total
flux in the limiting regimes stated, for the model's own parameters; a fit with more nuisance
parameters does worse.

### 8.4 Pipelines for trailed frames

| Source | What it does | Numbers | Verdict |
|---|---|---|---|
| Sharma et al. 2023, "Astreaks", *MNRAS*, arXiv:2306.16519 **[full]** | Template of points 1 px apart ⊛ Gaussian, rotated; angle checked against the tracking; streaks within `L/2` of the edge discarded; stars re-rendered as points for astrometry.net | 371 images, all solved, O − C 0.52″; minutes per image; target found by hand | SKIP the discard rule; Fayde flags edge trails and fits their visible pixels (streak-detection §4.6) |
| Danescu & Turcu 2026 (entry 8.1) | Target as a round component in a median stack of 3 to 7 frames; RANSAC on position against time | | multi-frame; noted for sequence mode |
| Spiller & Curti, *Acta Astronautica*, doi:10.1016/j.actaastro.2020.11.043 **[abstract]** | Under rotation, star trails are arcs; rate from chords | rate within 1% at 1 to 5°/s | SKIP unless trails curve |
| Samadzadegan & Alidoost 2013, *ISPRS Archives* XL-1/W3, 25 **[abstract]** | FFT-based streak features | under 5″ in simulation | SKIP |

### 8.5 Faint streaks: line integrals

| Source | Idea | Numbers | Verdict |
|---|---|---|---|
| Ofir, Galun, Alpert, Brandt, Nadler & Basri, "On detection of faint edges in noisy images", *IEEE TPAMI* 42(4), 894, 2020, arXiv:1706.07717 **[full]** | Bottom-up hierarchy: a length-2L line integral from two length-L ones, or for a new angle by interpolating four; positions along the line subsampled at `L/2`, so about 8N responses per length and O(N log N) in all; threshold for the *mean* of `w·L` pixels `σ·√(2 ln(K/δ) / (w·L))` per length (for a sum, `σ·√(2·w·L·ln(K/δ))`), conservative, checked on pure noise | about 3 s for 257² in research code | CONSIDER the hierarchy and the threshold, with a ridge profile in place of a step (not built) |
| Brandt & Dym 1999, *SIAM J. Sci. Comput.* 20(4), 1417 **[abstract]** | Fast approximation of all line integrals at all lengths, angles and positions within a set error | about n log n | CONSIDER (abstract only) |
| Arias-Castro, Donoho & Huo 2005, *IEEE Trans. Inf. Theory* 51(7), 2402 **[abstract]** | Multiscale (beamlet) detectors reach near-optimal thresholds where naive searches do not | | the theory behind the hierarchy |
| Brady 1998, *SIAM J. Comput.* 27(1), 107 **[abstract]**; Press 2006, *PNAS* 103(51), 19249 **[full]**; Götz & Druckmüller 1996, *Pattern Recognition* 29, 711 **[secondary]** | Dyadic approximate Radon transform, additions only, about `1.05 to 1.17 · n² log₂ n`; Press's median variant recovers a 0.5σ line among bright clutter | | basis of Nir et al. 2018 |
| Karpenko & Ershov 2021, *Probl. Inf. Transm.* 57(3) **[abstract]**, cited in Kazimirov et al. 2024, arXiv:2411.07351 **[full]** | Dyadic lines stray from true ones by up to `log₂n / 6` px, `n` the image side | 1.7 px at n = 1024 | CONSIDER: a band across the line wide enough for this and the PSF; the bound is for dyadic paths, not directly for Ofir's interpolated hierarchy |
| Horev, Nadler, Arias-Castro, Galun & Basri 2015, *SIAM J. Imaging Sci.* 8(1), 458 **[abstract]** | Long edges found by sampling a fraction of the pixels | | CONSIDER |
| Averbuch et al. 2008, *SIAM J. Sci. Comput.* 30(2), 764 and 785 **[abstract]** | Fast slant stack via the pseudopolar FFT | | SKIP: no locality for short segments |
| Šára & Cvrček 2017, 7th European Conf. on Space Debris, paper 403 **[abstract]**; Cvrček & Šára 2019, VISAPP, doi:10.5220/0007399804980509 **[abstract]** | Shear instead of rotation, per-angle background calibration, dynamic programming, branch-and-bound over angle, model selection | over 90% certified for streaks over 500 px at −10 dB | CONSIDER the per-image calibration |
| Dawson, Schneider & Kamath 2016, arXiv:1609.07158 **[full]** | Maximum likelihood reduces to the PSF-broadened line filter | 10 × below the pixel noise; unmasked stellar wings make false alarms | SKIP the search; heed the wings |
| Vananti, Schild & Schildknecht 2020, *Adv. Space Res.* 65, 364 **[abstract]** | Bank of streak filters | 2 to 3 h per 4k × 4k (as quoted elsewhere) | SKIP |
| Virtanen et al. 2016 (StreakDet), *Adv. Space Res.* 57, 1607 **[abstract]**; Pöntinen et al. 2020, *A&A* 644, A35, arXiv:2011.05297 **[full]** | Difference of mean filters, morphology, classification, PSF fit | 90% found at S/N over 1 for streaks over 100 px; about 1 min per 4k × 4k; 30 false per true before post-processing | SKIP |
| Han, Wen, Li & Qiao 2025, *Remote Sens.* 17(4), 631 **[abstract]** | Coarse-to-fine streak filter with a perpendicular cross-check | | CONSIDER the perpendicular veto |
| Bektešević & Vinković 2017, *MNRAS* 471, 2626 **[abstract]** | Known objects removed, then Hough; satellite trails narrower than meteor trails | | ADOPT the width veto |
| Stark et al., arXiv:2602.17816 (2026 posting of ACS ISR 2022-08) **[full]** | Median Radon | 10 s per 2048 × 1024 on 24 cores; 0.13 × pixel noise at 4580 px | SKIP |
| Tagawa et al. 2016, *Adv. Space Res.* 57(4), 929; Reed, Gagliardi & Shao 1983, *IEEE TAES* 19(6), 898; Blake et al. 2026 stacking **[abstract]** | Multi-frame | | SKIP for one frame |

### 8.6 False alarms and model choice

| Source | Idea | Verdict |
|---|---|---|
| von Gioi, Jakubowicz, Morel & Randall 2010, *IEEE TPAMI* 32(4), 722; IPOL 2 (2012), 35 **[abstract]**, and OpenCV's `lsd.cpp` **[code]** | Number of false alarms = number of tests × binomial tail of aligned pixels; accept at NFA ≤ 1; about `(WH)^{5/2}` tests | CONSIDER the bookkeeping, not the detector: gradient alignment is noise for a streak under the pixel noise |
| Grosjean & Moisan 2009, *J. Math. Imaging Vis.* 33(3), 313 **[secondary]** | NFA of a linear filter on colored Gaussian noise = tests × Gaussian tail | CONSIDER: the form of the line-sum test |
| Lezama, Morel, Randall & von Gioi 2015 (entry 7.4) **[abstract]** | A-contrario alignment of points, aware of local density | CONSIDER for dotted (tumbling) tracks |
| Killick, Fearnhead & Eckley 2012, *JASA* 107(500), 1590, arXiv:1101.1438 **[abstract]** | Exact penalized changepoint search; linear time under conditions on the changepoints, not in the worst case | CONSIDER for cutting a line's profile into on and off |
| Adler & Taylor 2007, *Random Fields and Geometry*, Springer **[secondary]** | Density of maxima of a smooth stationary Gaussian field above a high level `u`: `(√det Λ / 2π)·u·φ(u)`, `Λ` the covariance of the gradient | ADOPT for the aimed filter's threshold, with `Λ` from the model profiles |
| Self & Liang 1987, *JASA* 82, 605; Chernoff 1954; Davies 1987, *Biometrika* 74(1), 33 **[secondary]** | Likelihood ratios on a parameter boundary are a mixture of χ²; with a parameter undefined under the null, no standard distribution, and the usual BIC approximation is unjustified | CONSIDER the consequence: calibrate point against segment by simulation (Fayde's ratio is not calibrated) |
| Toldo & Fusiello 2008 (J-linkage); Magri & Fusiello 2014 (T-linkage); Barath & Matas 2019 (Progressive-X) **[abstract]** | Multi-model robust fitting | SKIP: randomized, or quadratic in points |

### 8.7 Line-segment detectors from computer vision

| Source | Cost | Verdict |
|---|---|---|
| LSD (von Gioi et al. 2010) | 36.5 ms for 640 × 480 | SKIP as a detector (entry 8.6) |
| EDLines, Akinlar & Topal 2011, *Pattern Recognit. Lett.* 32(13), 1633 **[abstract]** | 7.6 ms for 640 × 480 | SKIP: gradient thresholds |
| ELSED, Suárez, Buenaposada & Baumela 2022, *Pattern Recognit.* 127, 108619, arXiv:2108.03144 **[full]** | 5.4 ms for 640 × 480, about 17 ms per megapixel | SKIP as a detector; its gap jumping is worth copying into the tracer |
| MCMLSD, Almazàn, Tal, Qian & Elder 2017, CVPR, 2031 **[abstract]** | about 4.7 s for 640 × 480 | CONSIDER the structure: find the line globally, then cut it in one dimension |
| CannyLines, Lu, Yao, Li & Li 2015, ICIP, 507 **[abstract]** | | SKIP: edges in man-made scenes |

Verdict: ADOPT an aimed segment filter for trailed frames, with the frame's trail from moments
and checked by a segmentation-free estimate, and the erf fit, whose errors come from its inverse
normal matrix, with Bouquillon et al.'s bounds as a comparison (both built: streak-detection §4).
CONSIDER multiscale line sums with an a-contrario threshold for faint streaks (not built). The
parts are known; no surveyed pipeline combines them in one frame. Detection takes a median 40 ms on
Fayde's frames in trails mode (streak-detection §7).
