# Star and streak detection — design

Crate `fayde-detect`. Survey of existing methods:
[`docs/research/literature-review.md`](../research/literature-review.md) §7. Sources of every
algorithm: [`docs/references.md`](../references.md).

## 1. Scope

Input: one image whose stars are points (a frame tracked at the sidereal rate, or a spacecraft
holding its attitude), with possibly a satellite streak or two, hot pixels, cosmic-ray hits,
saturated stars, stray light. Output: what the solver's input contract asks of a source list
(`plate-solver.md` §6.3), filled honestly:

| Per source | How it is obtained |
|---|---|
| `x`, `y` | Closed-form fit to the peak of the matched-filtered image (§3.4) |
| `sigma_px` | Noise theory of that estimator, per source (§4) |
| `flux` | Integral of the fitted Gaussian |
| `p_star` | Local false discovery rate × shape class (§3.8) |
| `truncated` | The fit window reaches the frame edge |

Streaks are traced and reported (ends, angle, width, flux), and every detection on one is kept
out of the star list. Frames whose stars are trailed (tracking a satellite) are out of scope:
their star list still comes from the host's own detector, as before.

Design targets, in order: no parameter the user must know (star width, noise and threshold
are measured or derived; gain is read from the header when it is there); robust to cosmic rays
by the hundred and to stray light; position errors that the solver's verification and fit can
trust; speed; and the same bits on every machine (`fayde_core::Portable`, no platform maths).

## 2. Why not an existing detector

From the survey (literature review §7), what is in use today:

- None measures the star width and adapts its filter. astrometry.net's `simplexy` smooths with
  a fixed σ = 1 px; DAOFIND and photutils need the FWHM from the user; SExtractor a fixed 3 × 3
  kernel; cedar-detect a fixed 7-pixel gate with manual binning.
- None reports a per-star position error to the matcher.
- None runs a streak detector before solving. Streak fragments are rejected only as a side
  effect of size or shape gates, and astrometry.net removes only horizontal and vertical lines.
- Hot pixels and cosmic rays are handled by minimum-area or sharpness heuristics that conflict
  with undersampled stars.
- Licences: SExtractor GPL, SEP LGPL, astrometry.net GPL as distributed, cedar-detect
  source-available but not open source. Everything here is written from the papers, so no other
  project's code, and none of its licence terms, enters Fayde.

## 3. Method

```
image
  │ 1. matched filter (Gaussian; width measured, image binned if stars are wide)
  │ 2. background and noise mesh of the filtered image
  │ 3. peaks: local maxima above the false-alarm threshold
  │ 4. one closed-form fit per peak → position, size, shape, amplitude
  │ 5. local test: each peak against the ring around it
  │ 6. classes: star, saturated, extended, sharp
  │ 7. streaks: traced from ridge-like peaks
  │ 8. probabilities and errors
  ▼
stars (x, y, σ, flux, p_star) + streaks
```

### 3.1 Matched filter and scale

The image is convolved with a Gaussian of the stars' width, the filter of highest
signal-to-noise ratio for a point source in white noise. The width is not asked for. A first
pass at σk = 1.2 px (within 7% of the best ratio for stars of FWHM 1.9 to 4.1 px) measures it
(§3.3); a second pass is run only when the ratio `2·σk·σp / (σk² + σp²)` falls below 0.93. Stars
wider than σ = 2 px are binned first (block means), so wide-star frames cost less than
narrow-star ones. At most three passes.

### 3.2 Background and noise

Tiles of 64 px (the same sky area on a binned image, down to 24 binned px), measured on the
*filtered* image so that the noise is that of the thresholded quantity whatever the input's
noise correlation.

- Background: median of a sample of the tile, clipped once at median + 3σ.
- Noise: `1.4826 × MAD / √2` of the differences of pixel pairs separated by more than the
  kernel's reach, clipped up to three times. A difference cancels the background and its
  gradient, which a spread of values would count as noise.
- A tile standing above what its neighbours predict (a bright star, a galaxy) is replaced by the
  prediction: the largest, over the four lines through the tile, of the mean of its two
  neighbours on that line. SExtractor's 3 × 3 median flattens the corner of a frame lit from
  outside; this filter does not (the glow of stray light curves up along the direction to its
  source and down across it, and the line nearest that direction predicts at least the tile).
- Between tile centres the values are bilinear; beyond the outermost centres the background is
  extrapolated in a straight line, so a gradient running off the frame is followed to the edge.

**Clipped sky.** If at least a quarter of the valid raw pixels equal the frame minimum, Fayde
checks the distribution of those pixels in each tile. A tile is padding if at least 99% of its
valid pixels equal the floor, or if transitions between floor and above-floor pixels occur at
most a quarter as often as in a random mix. The detector checks the transition rate in tiles
with at least eight floor pixels.

Fayde uses the clipped-sky model when more than half of the sky tiles with at least 24 valid
pixels yield a fit. To fit a tile, at least a quarter of its valid pixels must be at the floor
and at least 16 must be above it. The singly left-censored normal fit uses the floor fraction
and an empirical quantile near the middle of the surviving values. Tiles without a fit use the
frame's typical fit. When the model applies, floor pixels in padding tiles are undefined for
detection. Otherwise, Fayde uses the standard filtered mesh.

The fit estimates the sky before clipping. Its clipped moments give the background and noise of
the recorded pixels. Fayde propagates these through the point filter and accounts for the skew
of the filtered noise when setting the peak threshold. On a simulated frame of 100 stars with
90% of the sky clipped at zero, 95 of 579 detections were stars before; with the model, 100 of
100, and every star is found. On 35 real frames with 25–36% of their pixels at zero, every frame
still solves and the fit residual fell on 31 of them (median 0.40 to 0.35 px).

### 3.3 Star width

For a peak of the filtered image, the variance of the star's profile is the fitted variance
(§3.4) less the kernel's. The frame's star width is taken from bright, round peaks, but they
are not all stars: on a warm or irradiated detector hot pixels and cosmic-ray hits outnumber
stars ten to one or more. Nothing that came through the optics is sharper than the optics allow,
and stars share one width where artefacts have many. So the stars are the **widest distinct
mode** of the size distribution: in a histogram of ln(size), the widest local maximum holding at
least 5% of the tallest one's count and separated from narrower sizes by a dip to half its
height. Blends and galaxies form a tail, not a mode.

The stars' *shape* is measured with their width: the median narrow-axis variance and the median
axis ratio of the peaks in that mode. A drift during the exposure smears every star alike. Each
source's narrow axis is then well under the stars' mean width, so narrowness (§3.6) and
elongation (§3.6, §3.7) are judged against what the frame's stars show, not against a round
star of the mean width.

### 3.4 One fit per peak

A quadratic surface is fitted by least squares to the *logarithm* of the 3 × 3 filtered pixels
around the maximum. For an elliptical Gaussian this is exact, and it gives at once:

- the position (no pixel-phase bias, unlike a quadratic fitted to the values themselves, whose
  bias reaches 0.03 px for these widths);
- the variances along both axes and the orientation (eigen-decomposition of the fitted Hessian
  `H`: variance along an axis = −1/eigenvalue);
- the amplitude, and so the flux `2π·A·σmajor·σminor`.

Cost: nine logarithms and about 40 multiplications per peak.

### 3.5 Local test

The mesh follows the frame on the scale of its tiles. Stray light and halos change faster, and
what the mesh misses is left as broad bumps. A star is narrow: three widths from its centre
nothing of it is left, so its height above the median of a ring of that radius (eight samples,
robust to a neighbour or a streak crossing the ring) is its height; a bump two or more
star-widths across falls short. Where the ring agrees with the mesh, the mesh (less noisy)
stands; elsewhere the ring is the background and the peak must exceed it by the threshold less
one σ. Saturated stars are exempt.

### 3.6 Classes

| Class | Test | Goes to the solver |
|---|---|---|
| `sharp` | Narrowest axis under 0.4 of the stars' narrow-axis variance, beyond the fit's noise; or, in the unfiltered pixels, the second faintest of the four neighbours of the brightest pixel under half of `exp(−1/σ²)` of it (light through the optics reaches all four; charge deposited in the detector does not) | no |
| `saturated` | Top at the saturation level (given, from the header, or the level of flat tops at the frame's brightest) | yes, `p` × 0.9 |
| `extended` | Elongated over 1.5 times the stars' own axis ratio, or over twice the stars' variance, beyond the fit's noise | yes, `p` × 0.5 |
| `streak` | On a traced streak | no |
| `star` | The rest | yes |

A saturated star's centre is the mean position of the *core* of its saturated region: the
pixels whose runs of saturated pixels along the row and along the column are each at least half
the longest. A bleed trail, a few pixels wide, has short runs across it and does not move the
centre. A saturated region that is a line of even width is a streak, not a star.

The work is bounded whatever the frame and the settings: each saturated pixel is flooded once
per frame, a region is followed for at most 200,000 pixels, a saturation level within 20σ of
the background is ignored, at most 200,000 peaks are measured and 1,000 streaks traced, and
binning never leaves an image under 16 pixels on a side.

### 3.7 Streaks

A peak on a streak is a ridge: the fitted surface curves down across the streak and hardly at
all along it. Each ridge-like peak (axis ratio at least twice the stars') seeds a trace: the ridge is followed one
pixel at a time in both directions along a straight line refitted by least squares as points are
added, re-centring on the brightest of five samples across the track, bridging gaps in
proportion to the length traced. Holding the trace to a line keeps it on course past a bright
star beside the streak or another streak crossing it. Then:

- isolated points at the ends, short groups beyond a gap, and a bright star touching an end are
  removed;
- the line is a total-least-squares fit; each end is put where the ridge falls to half its level
  (the true end of a uniform streak blurred by the PSF);
- the streak is traced once more from the middle of the fitted line, whose direction is better
  than the seed's;
- it is kept if it is at least 10 px and 8 cross-section σ long, half covered, with a total
  signal-to-noise ratio of 10, and if at least half of the bright peaks on it are ridges (a row
  of stars is not a streak);
- a streak narrower than the stars is flagged `cosmic_ray`.

Only streaks bright enough to raise a peak above the detection threshold are found (about 5σ per PSF
area). Fainter ones do not disturb the star list. Fayde does not search for them: an integration
along the line (Radon transform; Nir, Zackay & Ofek 2018) costs 10 to 100 times everything here.

### 3.8 Threshold and probabilities

The threshold is set from the number of noise peaks the caller accepts per frame (default 1):
the expected number of local maxima of filtered noise above `u` per pixel is bounded by
`min(P(Z > u), u·φ(u) / (4π·σk²))`, the pixel count and the peak density of a continuous
Gaussian field. On empty synthetic frames the bound is within 45% of the count for σk from 0.7
to 2 px. For a 4-megapixel frame the threshold is 5.0σ.

`p_star` = (1 − local false discovery rate) × class factor. The local rate is the expected
number of noise peaks within a quarter σ of a source's height over the number of sources seen
there. A source at the threshold in an empty region of height gets a low probability; the same
height among many real faint stars gets a high one.

## 4. Position error

The position is the maximum of the filtered image. Noise moves a maximum by (noise gradient) /
(signal curvature). For white noise σn filtered by a Gaussian of σk, and a Gaussian star whose
filtered profile has variance σt² = σp² + σk² and peak A:

```
Var(∂n/∂x) = N² / (2σk²),  N = filtered noise;   curvature = A / σt²
σx = σt² / (√2 · σk · SNR),  SNR = A / N
```

For a matched kernel (σk = σp) this is `√2·σp / SNR`, the Cramér–Rao bound for a Gaussian star
in background noise (0.60 FWHM / SNR). Everything in it is measured per source.

The star's own photon noise, with gain `g` given (N electrons): the same calculation with noise
variance proportional to the star's profile gives

```
Var = σp² · σt⁸ / (N · σk⁴ · (σk² + 2σp²)²)     = (16/9) · σp² / N  for a matched kernel
```

the known result for a Gaussian-weighted centroid (Mortensen et al. 2010, eq. 6). Without a gain
it is not modelled, and a floor of 0.05 px stands in. Saturated stars: 0.1 px + 3% of the radius
of the saturated core (measured on synthetic frames).

## 5. Measured

Synthetic frames (`fayde_detect::synth`: pixel-integrated Gaussian or Moffat stars, sky and read
noise, photon noise, streaks, hot pixels, cosmic rays, stray light, saturation). One laptop
core.

| Quantity | Result |
|---|---|
| Speed, 4 Mpx frame, stars of σ 1.5 px | 10 ms (2.4 ms/Mpx), one pass |
| Speed, stars of σ 0.6 / 2.5 / 4 / 8 px | 4 to 7 ms/Mpx (two or three passes; binned from σ 2.5 px up) |
| Star width | within 1% of the truth, σ 0.6 to 8 px |
| Noise | within 1% on sparse frames; reads up to 12% high when stars cover much of the frame |
| Completeness | every star above SNR 8; 21 of 28 between SNR 5 and 8 |
| Error / claimed error (robust spread; gain given) | 0.8 to 1.25 for σ 0.9, 1.5, 2.6 px |
| Error / Cramér–Rao bound, faint stars | 0.8 to 1.3 |
| False stars on empty frames | 3 in six 4 Mpx frames at the default (1 per frame asked) |
| Hot pixels (60 per frame) flagged | 60 of 60 at σ 0.7 and 1.5 px; 55 of 60 at σ 3 px (binned) |
| Cosmic rays, 15 to 30 hits per star | 0.3 to 0.6% of hits pass as faint stars at σ 1.5 px, 0.8% at σ 1.0 px; star width unaffected |
| Stray light, 20 × sky from a corner | 1 false star; without the local test, 33 |
| Streak ends | within 2 px (0.2 px typical); faint streaks (5σ ridge) within 10 px |
| Saturated stars, 3 to 300 × saturation, with bleed trails | centre within 0.6 px; one detection each |

End to end (`crates/fayde-api/tests/image_solve.rs`, a synthetic sky through a known WCS, two
streaks, hot pixels, saturated stars): worst sky error over the frame 0.11″ at 3″/px for
well-sampled stars (0.035 px), 0.08″ undersampled, 0.49″ for seeing-limited stars 7 px wide. The
solver's fitted error scale is 1.01 on the well-sampled frame: the claimed errors are the real
ones.

Real frames: 109 frames whose stars are points (two pixel scales; 3 to 10 s exposures, stars of
FWHM 2.2 to 3.3 px, 4 to 540 hot pixels and radiation hits per frame, up to 4 streaks or
cosmic-ray tracks in most frames), solved from the image file alone against UCAC4 with the
header pointing (2° radius) and a scale prior:

| | |
|---|---|
| Solved | 109 of 109 |
| Detection time | 2.0 ms per frame (median) |
| Matched stars | 181 (median) |
| Residual | 0.40 px RMS (median) |
| Error scale fitted by the solver | 1.58 (median): on these frames the claimed errors are too small by that factor, with or without the gain |

Two of those frames shaped the method:

- One, taken with a warm detector, has 300 radiation hits to 36 stars among its bright compact
  sources. It is what the widest-mode rule (§3.3) was made for.
- One is a crowded field whose stars are smeared two to one by a drift. With
  narrowness judged against the stars' mean width, its 36 brightest stars were classed as
  cosmic rays (117 "sharp" sources; 7 once judged against the stars' narrow axis), and with the
  brightest stars missing the solver's search, which starts from them, found nothing although
  417 of the remaining sources matched the catalogue. It now solves with 451 matched stars.

The detector this one replaced (after SExtractor: clipped-mode mesh, pixel cleaning, matched
filter, connected pixels, windowed centroids) was run on the same 109 frames with the same
solver and priors: 105 solved, median residual 0.79 px, error scale 2.65, about 20 ms per
frame. On the 105 frames both solve, the field centres agree to 0.19″ (median); on the one
where they differ by 10″, this detector's solution has a residual of 1.5 px and the other's 5 px.

## 6. Limits

1. Faint streaks (under about 5σ per PSF area) are not found (§3.7).
2. With stars of FWHM under about 2.5 px, multi-pixel cosmic-ray clumps of signal-to-noise ratio
   under about 30 cannot be told from faint stars by a 3 × 3 shape: in the worst synthetic case
   (3000 hits, 200 stars of FWHM 2.4 px, stray light and streaks together) 8% of the hits pass
   as faint stars. They are distractors
   the solver is built for.
3. The claimed position errors assume a Gaussian core and white noise. On real frames the
   solver's error scale shows them 1.5 times too small; the cause (pixel response, PSF shape,
   pointing jitter) has not been separated.
4. Undersampled stars (σ under 0.7 px): errors about 1.3 times the claim.
5. One thread. The filter and the scan split into strips without changing the result.
6. FITS images only (8 to 64-bit integers and floats); tile-compressed images are not read.
7. Trailed stars (tracking a satellite) are out of scope here. The streak tracer does find star
   trails, and turning them into sources with `Frame::trail` is a possible extension.
