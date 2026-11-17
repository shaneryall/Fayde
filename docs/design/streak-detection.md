# Trails and streaks in one detector

This design extends `fayde-detect` ([`detection.md`](detection.md)) to frames whose stars are
trails: a straight-trail pass (§4, measured in §7) and empirical kernels for curved paths (§8).
Literature: [`docs/research/literature-review.md`](../research/literature-review.md) §7.4 and
§8. Sources: [`docs/references.md`](../references.md).

## 1. Goal

One detector, `fayde-detect`, that takes any frame from an imager pointed at the sky and returns
what the solver needs, whatever moved during the exposure:

| Frame | Stars | Other things in it | Fayde |
|---|---|---|---|
| Sidereal tracking, or a spacecraft holding its attitude | points | satellite or debris streaks, cosmic-ray tracks | handled; streaks only above about 5σ per PSF area |
| Tracking a satellite, a slew, a drift | trails, all the same length and direction | the tracked object (a point), other satellites (streaks), cosmic rays | handled for straight trails at least 7 FWHM long (§7) |
| Either, with a faint satellite | as above | a streak below the per-pixel noise | not found |

In trails mode the output adds the frame's star trail with its evidence, point-like sources
(among them the tracked object), and per-source errors along and across the trail. The
constraints are the point pass's (`detection.md` §1): no parameter the user must know, the same
bits on every machine, no `unsafe`, no new dependency. Its cost is in §7.

**Supported geometry.** The straight-trail detector fits a common length and direction across
the frame. The empirical kernel can measure a shared curved path or field rotation when it
learns a reliable template (§8). When Fayde detects trails but cannot establish a consistent
model, the report records the reason and the detector uses its point measurements.

## 2. What exists

**The point pass.** `fayde-detect`'s pipeline for point sources applies a Gaussian matched filter of
the measured width, estimates the background and noise of the filtered image on a mesh, finds peaks,
fits each with a log-quadratic (position, size, axis ratio, direction), and classes them (star,
saturated, extended, sharp), then gives each a `p_star` from its class and its local false discovery
rate. Streaks are traced from ridge-like peaks (axis ratio at least twice the stars') along a
running least-squares line and reported with their ends; every peak on one is kept out of the star
list. A drift that smears every star a little is handled by judging shape against the frame's own
stars (`detection.md` §3.3). 109 of 109 real point-source frames solve, and point-pass detection
takes a median 2.2 ms per frame (§7).

**An earlier detector, `fayde-trails`.** About 3,100 lines, not merged:

- a line matched filter computed directly (cost proportional to pixels × trail length), with an
  AVX2 path and a thread pool;
- a clip loop, `min(image, 2 · line mean)`, iterated so that compact bright sources do not leak
  into the line response;
- the frame's trail from the flux-weighted second moments of connected components:
  `L = √(12 · (λmajor − λminor))`, direction from the major axis, with a verdict (`Trails`,
  `Compact`, `Ambiguous`, `PriorMismatch`, `NoGeometry`) and an optional prior that selects the
  population.

On 74 real frames tracking satellites, measured against a reference WCS per
frame: with a header prior its source lists solved 74 of 74; without one it refused 2 frames as
ambiguous and the other 72 solved. Median centre error 0.31 px against 0.50 px for
`fayde-detect`'s streak table used as a source list; 25 catalogue matches against 16. Those
figures predate its change to portable maths.

**The solver.** `Frame::trail` (length, angle) is accepted. A source flagged `truncated` counts as
evidence with a gate stretched along the trail but is never fitted or used as a seed. The
verifier can add a variance along one direction to a gate (`verify.rs`, an extra term in
quadrature). Everything else is isotropic: `Source::sigma_px` is one number, a match carries one
σ, a fitted pair one weight, and the WCS covariance one normal matrix shared by both coordinates
(`verify.rs`, `fit.rs`). Sources whose position falls outside the image are dropped
(`solve.rs`).

**Why a second pass.** The point pass alone solves 9 of 84 frames tracking satellites (§7), and
the costs of `fayde-trails` (`unsafe`, a thread pool, cost growing with trail length) conflict
with the workspace's rules.

## 3. What the literature gives

Detail and reading levels in the literature review, §8. The facts this design rests on:

1. **A trail needs its own filter.** The matched filter for a trail is the PSF convolved with the
   segment. Its signal-to-noise ratio over a point filter's tends to `√(L / (2√π·σ))` for
   `L ≫ σ` (Nir, Zackay & Ofek 2018). The LSST trailing-loss fits (Jones, SMTN-003) agree:

   | Trail length / FWHM | 2 | 5 | 10 | 20 |
   |---|---|---|---|---|
   | Gain, exact for a Gaussian PSF in white noise (computed) | 0.12 mag | 0.60 | 1.00 | 1.39 |
   | Gain, long-trail limit `√(L / 1.505 FWHM)` | 0.15 mag | 0.65 | 1.03 | 1.40 |
   | Gain, from SMTN-003's fits (point `a = 0.420, b = 0.003`; trail `a = 0.761, b = 1.162`) | 0.18 mag | 0.60 | 0.97 | 1.34 |

   Below about 2 FWHM the point filter gives up little; that is where the point pass stays.
2. **The filter must be aimed.** For long trails the matched response falls with an angle error
   `Δθ` as `(√π / 2z)·erf(z)`, `z = L·|Δθ| / 4σ`: 2% lost at `Δθ = σ/L`, 7.7% at `2σ/L`, 25% at
   `4σ/L`. A length off by a factor `k` keeps `1/√k` (a factor of 2 loses 29%, in the long-trail
   limit; 26% computed at `L = 40σ`), so lengths are measured, not taken from a coarse grid.
   (Derived; checked numerically. No published closed form was found.)
3. **One trail per frame is cheap to find.** The median orientation and major-axis length of
   bright elongated components gives it, and the median ignores a satellite streak or cosmic rays
   (Danescu & Turcu 2026: astrometry.net went from 15 to 39 of 40 frames solved once the image was
   filtered with that trail). The second moments give the length only for whole, evenly lit
   trails and a round PSF; components cut by the threshold read short.
4. **An aimed box costs O(1) per pixel, approximately.** A running sum along translates of one
   digital line (Soille, Breen & Jones 1996) costs two additions per pixel whatever the length.
   It is not translation invariant: the pixels a window covers change with its phase along the
   digital line, so the effective length and centring depend on the angle. The Gaussian applied
   after it hides most of that; how much is left has to be measured against a direct filter.
5. **The error of a trail is not round.** For a Gaussian PSF convolved with a segment (Vereš et al.
   2012; erf model as corrected by Bouquillon et al. 2017), when the background dominates, the
   error along the trail grows as `L` and across it as `√L`; when the star dominates, the error
   across is a static star's and the error along grows as `√L` (Bouquillon et al. 2017, Cramér–Rao
   bounds, for the model's own parameters; a fit with more nuisance parameters does worse). A
   brightness slope along the trail needs a tilt term (Blake et al. 2026); tracking jitter is
   common to every star in a frame, so an empirical profile stacked from bright trails beats the
   ideal model (Kouprianov 2008).
6. **Faint streaks are found by multiscale line sums.** Line integrals at a subsampled set of
   positions (spacing `L/2` along the line), angles and dyadic lengths, built bottom-up (Brandt &
   Dym 1999; Ofir et al. 2020), about 8N responses per length and O(N log N) in all. For the
   *mean* of `w·L` independent pixels the noise threshold from `K` tests is
   `σ·√(2 ln(K/δ) / (w·L))` (Ofir et al. 2020). Dyadic lines stray from true lines by up to
   `log₂ n / 6` px for an image side `n` (Karpenko & Ershov 2021): 1.7 px at `n = 1024`.
7. **One acceptance rule for every line, with an honest count.** The number of false alarms,
   `K · P(noise ≥ observed)` (von Gioi et al. 2010 and the a-contrario school; Grosjean & Moisan
   2009 for linear filters in Gaussian noise), is a union bound: it needs no independence between
   tests, but it needs valid null probabilities and every search counted. Clipping, masking,
   fitted-star subtraction, data-chosen geometry and endpoint optimisation all change the null,
   so a whole pipeline has to be calibrated on simulated frames, not only the formula.
8. **Point against segment is not a regular test.** Length 0 is on the boundary of the parameter
   space and the angle is undefined under it (Self & Liang 1987; Davies 1987). Wilks's χ² and the
   usual BIC approximation do not apply; penalties calibrated by simulation do.
9. **What others found useful as vetoes.** A streak as wide as the PSF (narrower is likely a
   cosmic ray: Bektešević & Vinković 2017); the response across against along (Han et al. 2025);
   a centroid that moves when the model is fitted means a blend (Blake et al. 2026); a trail
   within `L/2` of the edge is cut (Sharma et al. 2023 discard it; §4.6 uses it).

## 4. Method

```
image
  │ A. point pass (unchanged): filter, mesh, peaks, fits, classes, bright streaks
  │ B. trail census: one straight star trail common to the frame? → points | trails | refused
  │ C. segment pass (trails mode): masks, aimed segment filter, seeds
  │ D. one fit per trail on raw pixels: shared length, angle, width; position, flux, background free
  │ E. classes from a trail model and a point model: star, point, sharp
  ▼
stars (x, y, σ, σ along, σ across, flux, p_star) + streaks + point-like sources + frame trail
```

**Regression boundary.** The trail pass runs unless it is turned off (`--no-trails`). When the point
pass binned the image, B uses a full-resolution point pass, saved before binning or run again. A
frame that ends as `points` or `refused` keeps the point pass's star list bit for bit (the test
`point_frames_are_untouched` in `crates/fayde-detect/tests/trails.rs`). An empirical kernel (§8) can
still take a frame over: a field-rotation kernel after a `points` census, a shared curved template
after a refusal, or over an accepted straight trail under the stronger conditions in §8. The frame's
mode is then `kernel`.

### 4.1 Trail census (B)

Two estimates that fail in different ways; trails mode only when they agree.

1. **From shapes.** The filtered image is cut into 8-connected components of pixels above 4σ. A
   component votes when its peak is not masked and is at least 6σ, its variance ratio is at
   least 2.25 and it is at least 3 px long. Its vote is its direction and
   `L = √(12·(λmajor − λminor))`, weighted by its peak signal-to-noise ratio capped at 30, so one
   bright satellite cannot outvote fifty stars. The votes go into a circular histogram of the
   axis (an axis has no sign) in 36 bins of 5°, smoothed with weights 1, 2, 1. The largest mode
   gives `θ₁`, the weighted mean axis of the votes within 12.5° of it. A second mode more than
   15° away, with at least 5 votes and half the first mode's weight, is a second population. `L₁`
   is the weighted median length of the mode's votes whose peaks reach 18σ with no saturated
   pixel, or of all the mode's votes when fewer than 3 do. Taking the median length of the
   components follows Danescu & Turcu 2026 and `fayde-trails`; the weights and the 18σ choice are
   Fayde's. Components within 2 px of the edge do not vote for the length.
2. **From pixel pairs.** The autocovariance of the filtered image. Unmasked pixels above 3σ serve as
   anchors; at least 50 are needed, and large sets are thinned to about 5,000. Each partner counts
   as the ratio of its signal-to-noise ratio to its anchor's, at most 2, and a masked partner counts
   as 0. Masked are the frame's mask, a zone around traced streaks that are not star trails, and
   pixels at 0.8 of saturation or more, grown by 3 px. The direction `θ₂` is found coarse to fine:
   36 directions at a lag of four trail widths (less for short trails), refined with a parabola,
   then the lag doubled up to `L₁/2` with 17 angles a quarter of a width apart at that lag. A field
   of identical trails gives an elongated central peak; along `θ₂` its profile is a triangle of
   half-width `L` blurred by a Gaussian of `√2·σₜ` (σₜ the trail width with the filter's), plus a
   near-flat floor from the pairs of different stars. `L₂` is the `L` of the best fit of that model
   (amplitude and floor fitted linearly for each trial `L` from `L₁/2` to `2L₁` in steps of 0.25
   px), not where the profile reaches the floor.

The decision (`census.reason` in the report says which rule decided):

| Outcome | Rule |
|---|---|
| `points` | No vote; fewer than 5 votes in the mode; under half the weight in the mode (without a prior); elongated components carrying under 5% of the capped weight of all bright components; or `L₁` under 7 FWHM across the trail (real frames meant to show points drift by up to 6.8 FWHM, §7) |
| `trails` | `θ₁` and `θ₂` within `max(3σₜ/L₁, 0.08 rad)`, lengths within 20%, the straightness check passes, and the refit below accepts the geometry |
| `refused` | Two populations; the estimates disagree; the trails change across the frame; trails longer than the frame (no length votes, or fewer than 3 whole trails while at least 75% of the votes touch the edge); more than four round sources of the stars' width per trail vote; fewer than 50 anchors, or no fit of the second estimate with a positive amplitude; or the refit rejects the geometry |

**Straightness check.** The frame is cut into 3 × 3 regions. A region with at least 4 votes gives an
angle, and with at least 5 whole, unsaturated votes reaching 18σ a length. The frame is refused when
a region's angle is off by more than twice the angle tolerance, or its length by more than 50%. With
fewer than 2 regions measured, Fayde refuses the frame if fewer than 3 trails fit, or if the spread
of their fitted angles exceeds `max(2w/L, 0.05)` rad.

**Refit.** The frame's length, angle and width are the medians of up to 15 of the brightest isolated
whole trails, each fitted with all three free (at most 30 steps, χ²/dof at most 4); the census's
values stand when none settles. When at least 3 fits, and over half of all attempts, come out more
than 20% shorter than both the census length and their components' lengths, the geometry is
rejected: such fits follow part of a curve.

**Cleaned retry.** When radiation hits crowd the frame, or no trail emerges while at least 5
bright components are long and narrow (axis ratio 4 or more), the census runs once more on a
cleaned image, with only its anchors and votes changed (§7 lists the cleaning steps).

**Prior.** An angle (`--trail-angle`, Python `trail_angle`) selects a mode within 15° of it and then
waives the share and second-population rules; the geometry is still measured from the image, as
in `fayde-trails`. Detection takes no length, and no prior from the instrument profile or the
header. In `refused`, the point-pass star list is returned. The report's `census` holds both
estimates, the vote counts and shares, the elongated share, and the number of regions with their
worst angle and length deviations; `trail` holds the refined geometry with its errors.

### 4.2 Segment pass (C)

1. **Masks.** A point convolved with a segment kernel is a segment-shaped response, so every hot
   pixel and cosmic-ray hit would become a faint false trail. The filter's input is set to zero on
   the point pass's `sharp` detections and cosmic-ray tracks (when it ran at full resolution),
   pixels that are sharp at full resolution (`sharp_in_raw`) and undefined pixels. Saturated pixels
   stay in the filter and are left out of the fits. Traced streaks that are not star trails are
   subtracted with an erf model first. The raw image is kept for the fit (D).
2. **Aimed filter.** A running sum along lines at `θ`, taken on sheared rows and interpolated
   between rows, with fractional ends weighted, so that its length is the physical `L` at every
   angle (the running box of Soille, Breen & Jones 1996 uses digital lines instead). Then a
   Gaussian of the trail's width (σ between 0.6 and 4 px). A unit test checks the box against
   sums along the line.
3. **Threshold.** For a smooth stationary Gaussian field of unit variance, the density of maxima
   above a high level `u` tends to `(√det Λ / 2π)·u·φ(u)`, where `Λ` is the covariance of the
   field's gradient (Adler & Taylor 2007). The point pass's threshold is the smaller of that bound
   for a Gaussian kernel and the pixel tail. For the aimed filter, `Λ` comes from the model profiles
   along and across the trail, sampled every 0.25 px; along a long trail it scales as `1/(√π·σ·L)`.
   On frames whose sky is clipped at a floor the automatic threshold is at least 5σ. The bound is
   checked on empty frames for the point pass only.
4. **Seeds.** Filtering a trail with a box of its own length gives a blurred triangular response
   along the trail with its maximum at the trail's centre, so each maximum is a seed near the
   centre. A maximum within `L/2` along the trail and 1.5 filter widths across of a brighter one
   is dropped; at most 5,000 seeds are kept.
5. **Noise.** The mesh measures noise from differences of pixel pairs set apart (`mesh.rs`). In
   trails mode the pairs are set more than a trail apart (`L + 3σk` for the filtered image,
   `L + 6w` for the fit's weights), so that two pixels of one trail are not paired and the noise
   is not underestimated.

### 4.3 One fit per trail (D)

1. **Model.** An erf along the trail (Bouquillon et al. 2017) times a Gaussian across, of variance
   `w² + 1/12` for the pixel's width, evaluated at pixel centres on the raw image. The frame's `L`,
   `θ` and width are fixed; position, flux and a background plane (a gradient, for stray light) are
   free. Weights are `1 / (σ² + s / gain)`, with σ the mesh's sky noise and `s` the model's signal
   above the background plus any subtracted neighbour or streak light (the photon term only when the
   gain is known), recomputed at each step. On frames whose sky is clipped at a floor, a prior holds
   the background.
2. **Start** from the centre of light along the line (halfway there on clipped-sky frames) and the
   seed's position across it. A fit that
   moves more than 1 px is started once more from where it ended.
3. **Iterate** at most 8 bounded Gauss–Newton steps (24 on clipped-sky frames). A step is capped at
   `max(0.1·L, 1)` px, and a step that reverses is halved. The fit converges when its step is under
   a tenth of the position error (between 0.001 and 0.05 px) and counts when it is under half. A fit
   that leaves its stamp, cannot be solved or does not settle is dropped, and so is a trail known no
   better than half a trail along (one FWHM on clipped-sky frames) or two widths across, or a point
   known no better than two widths either way. Fixed limits keep the work and the bits fixed.
4. **Errors** from the inverse normal matrix at the parameters returned, projected along and
   across the trail (`σ_along`, `σ_across`) and floored. An optional test
   (`error_claims_cover_the_truth`) checks them on simulated frames: for matched stars away from
   the edge, 58 to 78% of the residuals within 1σ and at least 86% within 2σ, on each axis.
5. **Neighbours.** Fits that converge within 1.5 widths of each other are merged when both fit
   well (χ²/dof at most 2). The others are fitted again with their neighbours' models subtracted
   where their light reaches the stamp; the first fit stays if the second fails.

### 4.4 Classes (E)

Each seed is fitted as a trail. It is also fitted as a point (the PSF, `L` = 0) when 40% of
its light along the line lies within `2w + 1` px of its brightest column, or its trail fit has
χ²/dof above 2. The two are compared on a common stamp centred on the fitted trail:

| Model | Free | Class (morphology) | To the solver |
|---|---|---|---|
| Star trail: the frame's `L`, `θ`, width | position, flux, background | `star` | yes |
| Point: the PSF | position, flux, background | `point` | no; reported |
| Point narrower than a star can be | | `sharp` | no |

With `d = ½(χ²point − χ²trail)` on that stamp, the trail's probability against a point is
`p_trail = 1 / (1 + e^(−d))` (equal priors): a seed with `p_trail` of at least 0.5 is a `star`,
otherwise a `point`. A point whose brightest pixel, less 3σ, holds more than 2.5 times the share
of light a star puts in its brightest pixel is `sharp` (a hot pixel or a cosmic ray). `p_star` is
(1 − local false discovery rate) × `p_trail`. Likelihood ratios between these models do not
follow χ² (§3, item 8), so `p_trail` is a relative weight, not a calibrated probability.

The classes describe shape, not what the object is: a `point` in a trailed frame is often the
tracked object but can be anything moving with it. A cosmic-ray track that lies along the trail
can pass as a star.

### 4.5 Streaks

The point pass traces streaks before the census, with no trail direction: its seeds are peaks
elongated at least twice as much as the stars. After the census, a traced streak is removed with
its sources when it looks like a star trail (within `3·max(2w/L, 0.05)` rad of the trail, at most
3 trails long, not a cosmic ray) or a chain of star trails (within 12° and uneven along its
length); the aimed filter's seeds then measure that light. Other streaks are subtracted before the
trail pass, and a seed on one is kept only if its fit has χ²/dof of at most 1.5.

### 4.6 Trails cut by the edge

A trail whose outer end lies beyond the edge has a centroid biased by up to `L/2` along it when it
is measured as a blob. The fit holds the frame's `L` and `θ` fixed and fits the visible pixels, so
it places a cut trail's centre where the whole trail's would be, which can fall outside the image
(on clipped-sky frames such sources are dropped). A trail is flagged `truncated` when its centre
is closer to an edge than half the trail's length projected towards that edge plus two widths.
The solver scores a truncated source with a gate stretched along the trail but neither seeds nor
fits with it, and drops a source whose centre is outside the image (§2).

On a square frame whose trails are a tenth of its width, about one star in eight is cut (12.4%
averaged over directions, 10 to 13.6% depending on the direction).

## 5. Interfaces

| Where | Change | Compatibility |
|---|---|---|
| Detection report | Top-level `mode` (`points`, `trails`, `kernel`, `refused`), `census` with its evidence, `trail` for the straight-trail model and `kernel` for the empirical model; per source `sigma_along_px`, `sigma_across_px`; class `point` | adds fields |
| Source list (FITS) | Header keywords `TRAILLEN` (pixels) and `TRAILANG` (degrees from +x towards +y), written by `fayde detect -o` and read by `fayde solve --sources`; `SIGMA` is the larger of the two errors | adds keywords |
| Command line, Python | No new required option; `--trail-angle` (Python `trail_angle`) chooses among trail directions in detection (the length is not needed); Python `Solver.solve_detection` | adds options. New public fields in `fayde-detect`, `fayde-api` and `fayde-fits` break Rust code that builds those structs with a literal, and an exhaustive match on `Class` must handle `Point` |

## 6. Integration choices, and what may be new

Most of this is known technique put together: a segment-shaped PSF fitted with its neighbours'
models subtracted is PSF photometry (Stetson 1987) with a longer PSF, and Houtz & Frueh 2021 and
Sharma et al. 2023 already detect trails before plate solving. One part has no antecedent that the
survey (§8 of the literature review) found, which does not make it new: the census's refusal rule, a
segmentation-based estimate checked against a segmentation-free one, with regional straightness, and
no trails mode unless they agree.

## 7. As built

**Choices made in testing.** Each row is a first design that the synthetic or real frames changed;
§4 describes the result.

| First design | Built | Why |
|---|---|---|
| Points below 2 FWHM | Points below 7 FWHM across the trail | Real frames meant to show points drift by up to 6.8 FWHM; the point pass solves them, and trails mode solved two of them worse |
| Autocovariance weighted by brightness | Each anchor counts once: partner over anchor, at most 2 | Weighted by brightness, saturated stars' bleed trails set the direction; clipped, bright trails read 20 to 30% long |
| Bright elongated components vote by capped signal-to-noise ratio | When short, narrow, differently oriented tracks are common, mask sharp charge deposits before smoothing and exclude narrow tracks from the new votes. Keep traced streaks supported by the cleaned pixels. | Radiation tracks can outvote stars and bias the point pass’s width and masks. |
| Use the point pass’s background and widths | On a retry, subtract a finer local background from the raw pixels, rebuild the filtered image and estimate the width from whole, bright components. Use the cleaned pixels for trail measurements if the census passes. | Smooth structure left by the coarse mesh contributes to both estimates and hides faint trails. |
| Remove known periodic patterns from an instrument profile | On a retry, also find strong short-period waves in the frame and subtract them with the existing tiled fit. | A strong readout pattern can set the autocovariance direction when no profile is supplied. |
| Use every filtered pixel above 3σ as an autocovariance anchor | On a cleaned frame, anchor on bright elongated components from all directions. Sample partners anywhere and search the full 180°. | Background pixels can overwhelm the peak. Sampling all eligible directions keeps the second estimate independent of the winning direction. |
| Check regions for agreement, or check individual fits when too few regions have enough votes | A cleaned census also requires three isolated sources to fit the straight-trail model. | A curve repeated across the frame can pass the regional check. |
| Refine the geometry with free fits to isolated trails | Reject the geometry if a majority of attempts, with at least three fits, produce segments more than 20% shorter than both the census length and their component lengths. Keep an accepted original census if cleaning reduces its regional coverage. | A fit can shorten a segment to follow part of a curve, even when both census estimates agree. Cleaning can also reduce the number of regions with usable trails. |
| Census on the point pass's image | At full resolution when the point pass binned | Trails read as wide stars make the point pass bin; binned, short trails look round |
| Seeds start the fit | The fit starts from the centre of light along the line, and is repeated once from where it ended | A seed lies anywhere on the filter's plateau within half a trail of a point source |
| One fit per trail | Fitted again with the neighbours' models subtracted where their light reaches the stamp | Trails in line or side by side share pixels |
| Streaks kept as streaks | Traced streaks that are not star trails are subtracted (as traced) before the trail pass; a seed on one is kept only if the better of the two models fits it | Seeds along a satellite near the trail direction made false stars; stars it crosses are kept |
| Hot pixels from the point pass | Also by their own pixels at the trail pass's resolution (`sharp_in_raw`) | A binned point pass misses them at full resolution |
| Point model for every seed | Only for seeds whose light is concentrated along the line, or whose trail fit is poor | Time |
| `erf` and `exp` per pixel | Tables 1/512 apart, interpolated (error under 1e-6) | Time |
| Fits as written | Steps that reverse are halved; a fit counts when its last step is under half its position error; χ² and covariance at the parameters returned; positions known no better than half a trail along or two widths across are dropped; trail and point compared on the same pixels | Found in testing: fits on faint long trails and trails cut by the edge wander along the trail; a fit with no light has an unbounded error |
| Satellite light subtracted | Its photon noise kept in the weights; saturation judged on the image before subtraction | Found in testing |
| One sheared copy of the image | Two sheared rows at a time, each only where it crosses the image | A frame far wider than tall needed its width squared |
| Census for mixed frames | Refused when round sources of the stars' width outnumber the trails four to one | Point stars with aligned satellite segments could pass as trailed; radiation hits on real trailed frames reach 3.5 to one |
| Trails mode final | A frame measured as trails that gives no solution is solved again with its stars measured as points (`fayde solve --image`, Python `solve_image`) | A slew ending in a dwell leaves curved trails with a bright point at the end: trails mode failed where the points solved |

**Synthetic frames** (`crates/fayde-detect/tests/trails.rs`; 768 × 768, PSF σ 1.3 px, 120 stars,
12 frames per length, errors with no floor):

| Trail | Stars | Median error along / across | Within 1σ along / across | Within 2σ along / across |
|---|---|---|---|---|
| 24 px | 1,218 | −0.014 / −0.001 px | 0.69 / 0.69 | 0.95 / 0.95 |
| 36 px | 1,049 | −0.001 / 0.000 px | 0.71 / 0.68 | 0.92 / 0.94 |
| 60 px | 624 | −0.006 / 0.012 px | 0.66 / 0.67 | 0.88 / 0.92 |

The length and direction come out within 0.2 px and 0.05°. The tracked object is found as a
`point`; hot pixels, cosmic rays and a satellite crossing the frame make no false star. Faint
long trails have heavier tails along the trail than a Gaussian (0.88 within 2σ at 60 px).

**Real frames** (one laptop core, an instrument profile for the camera, UCAC4 with its index, the
frame's header pointing):

| | Point pass only | With the trail pass |
|---|---|---|
| 109 frames whose stars are points: star lists identical to the point pass's | | 109 of 109 |
| Detection time, median, those frames | 2.2 ms | 3.0 ms |
| 84 frames tracking satellites, trails 9 to 85 px (header rate × exposure): solved | 9 | 66 |
| Of those, in trails mode | | 58, all solved |
| Detection time, median, trails mode | | 40 ms |
| Median residual of the solutions, trailed frames solved | 1.6 px (9 frames) | 1.64 px (66 frames) |

No frame that the point pass alone solves is lost.

**Limits.**

- Cost: trails mode takes a median 40 ms per frame, against 2.2 ms for the point pass alone and
  3.0 ms with the census on point-source frames. The fits are most of it.
- 26 of the 84 frames are not in trails mode (14 called points, most with trails under 20 px; 12
  refused, mostly because the two estimates disagree on real trails, whose brightness is not uniform
  along them), and 18 stay unsolved. The census also misses faint long trails: a 160 px trail
  spreads a star's light so thin that it may not stand out in the point filter the census votes
  from.
- The census refuses two trail populations of comparable weight unless an angle prior selects
  one, and refuses the geometry when too few trails have both ends in the frame.
- For straight trails, the solver gets one isotropic error per source, the larger of the errors
  along and across the trail, so the error across the trail is overstated.
- A cut trail's error is the fit's own: the error of the frame's shared length and angle is not
  carried as one term common to every cut trail. On real frames, 411 cut trails had a spread
  along the trail of 1.65 times the claimed error, against 1.25 for whole trails.
- Trails that nearly coincide are not split by a joint fit. Where fitted neighbours overlap, a
  seed is fitted again with their models subtracted; the first fit stays if the second fails.
- The straight-trail fit has no term for brightness that changes along the trail.
- Streaks fainter than about 5σ per PSF area are not found.

## 8. Empirical kernels for curved star paths

Tracking, rolling, or wobbling during an exposure can spread each star along a curve.
A straight-trail filter can lock onto fragments of the path, giving the solver too few
stars and positions displaced from their exposure-averaged coordinates.

### How it works

1. **Select training stars.** Reuse the census components and its defect mask. For a
   shared path, keep up to 64 bright, isolated, unsaturated components away from the
   frame edge, with lengths and orientations in the dominant mode. Subtract the
   background with the mesh, normalize by flux, and register the cutouts by their
   centres of light. Align them twice to the current mean using subpixel correlation
   peaks. Compute a clipped, noise-weighted mean and recenter it on its own light.
   The empirical PSF follows the approach of DAOPHOT and PSFEx. The origin
   represents the exposure-averaged position of the star.
2. **Choose the mode.** When the straight-trail census accepts or refuses a dominant
   trail, measure the template's centre of light in 2 px slices along its major axis.
   Use the empirical template when the path's RMS deviation from the best straight
   line is at least 0.5 px. When the census has already accepted one straight trail,
   the template overrides it only if it stacks at least 15 stars and bends by at least
   1 px. On real satellite-tracking frames, templates of 6 to 11 stars bent by 0.6 to
   3 px moved solutions by up to 1.7 px from the straight-trail fits, which agreed
   with Gaia to 0.05–0.1 px. Every template must stack at least 10 stars: its centre of
   light is every position's origin, so its error moves the whole solution (a
   template of 5 stars on a sparse simulated frame put the solution 2.6 px off). The
   sum of distances between slices estimates the path length. Evidence of field
   rotation can also select `kernel` mode when the census finds no common direction.
3. **Correlate and find peaks.** Correlate the background-subtracted full-resolution
   image with the taps that hold 99% of the template's light. Mask undefined pixels,
   hot pixels, and cosmic-ray hits using the trail pass's mask. Use the mesh and the
   existing peak finder with an equivalent Gaussian scale from the template's
   gradient covariance. Visit the peaks in descending correlation amplitude.
4. **Fit and suppress sidelobes.** Start with a quadratic peak interpolation and
   refine the position and flux against the template in the image. Require the
   fitted signal to clear the noise threshold and agree with the template across its
   support. Subtract each accepted source in image space. Testing the correlation
   on the residual image applies the predicted autocorrelation of every brighter
   source with its actual subpixel shift and mask. A second correlation pass
   recovers peaks that were hidden by bright neighbours. This is a bounded,
   CLEAN-style subtraction.
5. **Report the fit.** Use the matched-filter amplitude for flux and the template's
   gradient Fisher matrix for position errors. The weights include photon noise
   when the gain is available and an allowance for template mismatch. Compute
   `p_star` from the local false discovery rate. Report template size, training
   count, path length, bend, and region count. An unsuccessful image solve retries
   with point detections.

### Field rotation and limits

For small field rotations, elongated components trace tangents to circles about a
common centre. Fit the centre from their directions and require their lengths to scale
with radius. Stack samples across the middle of isolated trails to learn a common
transverse profile. Use it with a Gaussian-smoothed longitudinal segment along a
quadratic approximation to the local arc. The grid spacing limits the change in
endpoint displacement within a region to about one pixel, or about 2.4 px at fast
rotation, where the spacing reaches its 32 px minimum. For each fit, evaluate
the template at the candidate's position.

The rotation model assumes a single centre in the frame and an approximately uniform angular sweep.
The shared template models one stable shape. Compact wobble that the census reads as point-like
follows the point path. Large background variations and overlapping trails can still hide faint
stars. The dense correlation is about 2–4 × 10¹⁰ multiply-adds on a 4096 × 4096 frame with a
large template (several seconds), and rotation mode builds a local template for each candidate
that reaches fitting. Next to a very bright star, residuals of a few percent of template mismatch
could leave peaks along its path; such peaks must pass the fit-quality, signal-to-noise and
separation tests, and none has been seen on the test frames.

### Measured results

The simulations use 1024 × 1024 images at 2 arcsec/px. Each group has five
frames. The scorer matches sources to exposure-averaged truth positions within
3 px. Completeness counts matches among detectable stars. Purity counts matched
detections. Centroid p90 is the 90th-percentile error of the matched detections.
Each cell shows the range from the five frames before and after this change,
with the detector before this change as the baseline.

| Frame group | Completeness, before → after | Purity, before → after | Centroid p90 in px, before → after |
|---|---|---|---|
| `base-L0-curve` | 0.182–0.265 → 0.901–0.971 | 0.285–0.367 → 0.980–0.994 | 2.567–2.794 → 0.384–0.503 |
| `base-L0-roll` | 0.282–0.429 → 0.902–0.954 | 0.394–0.608 → 0.948–0.983 | 2.235–2.543 → 0.485–0.740 |
| `base-L0-still` | 0.955–0.994 → 0.955–0.994 | 0.988–1.000 → 0.988–1.000 | 0.099–0.112 → 0.099–0.112 |
| `base-L0-trail` | 0.731–0.897 → 0.731–0.897 | 0.793–0.940 → 0.793–0.940 | 0.948–1.459 → 0.948–1.459 |
| `cosmics-L2-curve` | 0.017–0.259 → 0.864–0.931 | 0.011–0.330 → 0.933–1.000 | 2.652–2.931 → 0.380–0.479 |
| `hot_pixels-L2-curve` | 0.049–0.272 → 0.922–0.960 | 0.038–0.327 → 0.934–0.989 | 2.820–2.916 → 0.409–0.465 |
| `smooth_field-L2-curve` | 0.000–0.098 → 0.562–0.723 | 0.000–1.000 → 0.827–0.978 | 1.730–2.652 → 0.216–0.330 |
| `wobble-L3-still` | 0.957–0.977 → 0.957–0.977 | 0.977–1.000 → 0.977–1.000 | 0.162–0.626 → 0.162–0.626 |
| `wobble-L3-trail` | 0.780–0.892 → 0.922–0.977 | 0.796–0.898 → 0.943–0.977 | 1.252–1.686 → 0.449–0.608 |
| `wobble-L4-still` | 0.927–0.982 → 0.927–0.982 | 0.524–1.000 → 0.524–1.000 | 0.498–2.037 → 0.498–2.037 |
| `wobble-L4-trail` | 0.000–0.723 → 0.927–0.988 | 0.000–0.904 → 0.944–0.982 | 1.394–2.561 → 0.352–0.444 |

One baseline frame each in `smooth_field-L2-curve` and `wobble-L4-trail` had no matches
and is omitted from the corresponding p90 range. All curved, rotating, and trailed
wobble frames used `kernel` mode after the change. All ten point and straight-trail
control FITS source lists compared byte-for-byte with `cmp`.

On arm64 macOS, median detection time for `base-L0-curve-f0` was 56.7 ms before and
76.6 ms after. These are the detector's reported times for 11 interleaved runs per
binary after one warm-up each. A separate synthetic curved frame solved from its
image detections within 0.3 px at the centre and all four corners.
