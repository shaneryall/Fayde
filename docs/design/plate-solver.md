# Fayde plate solver — design

This document describes the solver design. For the crate layout and integration contract, see
§4 and §6. For commands and report fields, see the [CLI guide](../usage.md) and
[Python guide](../python.md). Source references in the
prior-art survey (§3) point to the upstream repositories.

## 1. Purpose

Fayde is a public, open-source astrometry tool written in Rust. It goes from a source list to a
verified WCS with measured accuracy. It must fully meet the needs of a host pipeline, so that the
pipeline can retire `solve-field` and the solving-related code around it.

**Input is an xylist, as with `solve-field`, or an image.** The user can run whatever source
detector suits their sensor and give Fayde the resulting list of sources (§6) along with any
header or metadata they have. Fayde also takes the image itself when its stars
are points (a frame tracked at the sidereal rate, a spacecraft holding its attitude), detects the
stars, keeps streaks and artefacts out of the list, and fills the same contract:
[`detection.md`](detection.md). It does the same when every star is a trail of
one length and direction (a telescope tracking a satellite):
[`streak-detection.md`](streak-detection.md).

It is not a wrapper. It ports the strongest algorithms from three solvers, and it takes over
capabilities that pipelines otherwise build around a solver, generalized for public use:

| Concern | Taken from | Reason |
|---|---|---|
| Fast hinted match | tetra3 / cedar-solve; Mortari's pair lists | Pair seeds when the scale is tight, otherwise a 4-star edge-ratio pattern hash; millisecond solves |
| Blind / wide fallback | astrometry.net, adapted | Its series of indexes for different field widths, searched with Fayde's 5-ratio pattern hash (quad codes considered, not adopted) |
| Match acceptance | astrometry.net | Bayesian log-odds verification against the catalogue |
| Final fit | astrometry.net, extended | TAN + SIP fit, upgraded to weighted iterative refit with outlier rejection |
| Source detection | new, for point-like and trailed stars; otherwise the user supplies the xylist | One fixed extractor cannot fit every sensor; astrometry.net's built-in `simplexy` solved about 12% of one host pipeline's frames. Fayde defines a rich xylist contract (§6) that any detector can fill, and fills it itself from an image of point-like stars with a detector that measures the star width, noise and position errors instead of asking for them ([`detection.md`](detection.md)). |
| Catalogue | new | UCAC4 and Tycho-2 read natively behind a catalogue trait |
| Header priors, neighbour hints | new | Scale and pointing resolution and a prediction from a solved neighbour, which a pipeline otherwise builds around the solver |

### Standalone first, host pipelines first-class

Fayde is a general-purpose public tool. A host pipeline is a user, not a dependency. Fayde takes
on capabilities that pipelines otherwise build around a solver when they belong in a public
astrometry tool.

Detection of point-like and trailed stars is one of them. A host may still prefer its own
detector and give the xylist.

**In Fayde** (so that a host pipeline need not build them):

| Capability | Fayde home |
|---|---|
| Star and streak detection from the image, points or trails | `fayde-detect` |
| Source ranking and demotion of near-duplicate centroids | `fayde-solve` (§5.0) |
| Scale prior from focal length and pixel size, or from CD or PC+CDELT with unit checks; rejection of mixed or invalid WCS | `fayde-fits` |
| Pointing prior from the header, never treating a commanded WCS as solved | `fayde-fits` |
| Exposure time resolution (DATE-OBS, JD-OBS, TIME-OBS, EXPTIME and variants) | `fayde-fits` |
| UCAC4 cone search; "is there a catalogue star here?" check | `fayde-catalog` (cone search) + `fayde verify --at` |
| Verifying or repairing existing WCS headers | `fayde verify`, `fayde refit` |
| Writing the solution atomically as a separate FITS WCS file, leaving the image untouched | `fayde-fits` (WCS file), `fayde-cli` (provenance cards) |

**Stays in the host pipeline** (pipeline policy, not astrometry):

- WCS cache reuse policy, versioned header cards and policy fingerprints.
- Run-summary and evidence schemas.
- Satellite identification, correlation, stacking, and drift.
- A detector of the host's own, where it prefers one to Fayde's.

A host's solving path is then: its detector → xylist → Fayde → adapter → its own cache and
records.

General rules:

- Core crates know nothing about any host pipeline: no host field names, config keys, cache
  records or FITS history cards.
- Fayde's own inputs and outputs are general: a source list, hints, and a solution with full
  evidence. They are designed so that everything a host pipeline needs maps onto them with no
  loss (§6).
- A host adapter is thin and lives at the edge: `fayde-py` plus a backend in the host's own
  repository. Host-specific provenance (cache records, run summaries) is built there from
  Fayde's generic provenance.
- Other users get the same features through the CLI and the Rust API. These include ground
  telescopes, star trackers, and any `solve-field` user who is willing to switch.

### Why not keep astrometry.net

- Maintenance is sporadic; build and packaging are C/Python/SWIG and hard to vendor.
- Every solve is a subprocess with scratch files, a 180 s budget and text-file I/O.
- It does not expose solution covariance, so the accuracy of an accepted solution is not
  quantified.
- Index files are fixed 2MASS/Tycho-2/Gaia-DR2 builds (series 4200/5200); epoch handling is not
  declared, and we cannot build from UCAC4 with our own epoch and density policy.
- No native notion of an image sequence; a pipeline has to implement neighbour hints and
  borrowing around it.

## 2. Requirements

### Must

1. Define an input contract optimized for solution quality (§6). Only the sources and image size
   are required. Every other field (uncertainties, flux, `p_star`, exposure time,
   priors with σ) improves the result when present. **Internal pixel convention:** `(x, y)` = (column, row), with the
   centre of the first pixel at `(0, 0)` (NumPy/Astropy). Conversions happen only at the edges:
   +1 for FITS/astrometry.net xylists, +0.5 and an `(x, y)`→`(y, x)` swap for tetra3/cedar.
2. Use priors when present: pointing, roll, scale and parity, each given as a σ plus a hard limit.
   The limit bounds the search; the σ ranks hypotheses (§6.4).
3. Solve blind (no position hint) within a bounded time when the hint is missing or wrong.
4. Build indexes from **UCAC4** zone files (`z001`–`z900`, 78-byte records), with proper motion
   propagated to a chosen epoch and per-cell density limiting. Catalogue access goes through a
   trait so Tycho-2, Gaia DR3 and others can be added.
5. Output a FITS-compatible TAN or TAN-SIP WCS, plus: matched pairs, per-star residuals, RMS,
   parameter covariance, the log-odds score, index/catalogue identity and epoch.
6. Deterministic: same input, same index → bit-identical output. Required by a host pipeline's
   WCS cache fingerprints.
7. Library first (no subprocess), with a CLI and Python bindings for host pipelines.
8. Read xylists in Fayde's native format (§6.8) and in astrometry.net's FITS xylist format. A
   `solve-field` user should be able to point Fayde at their existing `.xyls` files and get a
   solution with no other changes.
9. Build priors (scale, pointing, exposure time) from standard FITS headers with strict validation.
10. Public tool hygiene: permissive licence, documented CLI, published crates and wheels, and
    reproducible index builds anyone can run from public catalogues.

### Should

1. Sequence mode: when the search finds no solution, verify a WCS predicted from a solved frame
   of the same pass.
2. Per-source centroid uncertainty: each source may carry a position error that the verifier
   and fitter use (§6.3). Sources are treated as point centroids; producers with streaked stars
   supply the streak centroid.
3. Solve the target case (sub-degree FOV, ~3″/px, 30–60 sources) in < 100 ms hinted and
   < 5 s blind on a single core.

### Won't

- Lens-distortion models beyond SIP (TPV, lookup tables).

## 3. Prior art

### 3.1 astrometry.net (`dstndstn/astrometry.net`)

- **Input:** an image or a FITS xylist (`X`, `Y`, optional `FLUX`/`BACKGROUND`). Coordinates are
  1-based, so the centre of the first pixel is (1, 1) (`util/image2xy.c:102-105`). The list is
  sorted brightest first, and `--depth` limits how many objects are tried. The built-in extractor
  (simplexy) does median background subtraction, 8σ peaks and a 3×3 centroid
  (`include/astrometry/simplexy.h:10-20`).
- **Hints:** scale in degwidth, arcminwidth, arcsec/px or focal length
  (`solver/augment-xylist.c:83-89`), an RA/Dec centre with a radius (`solver/engine.c:496-576`),
  parity, and a CPU limit.
- **Pattern:** a quad of 4 stars. A and B are the widest pair. C and D are expressed in the frame
  where A=(0,0) and B=(1,1), giving a 4-D code that does not change under rotation, scaling or
  translation (`solver/quad-utils.c:15-55`). Requiring cx ≤ dx and flipping for parity removes
  duplicate orderings (`quad-utils.c:102-150`). Code tolerance is 0.01 (`solver.h:30`). Field stars
  are added one at a time and every quad containing the new star is tried
  (`solver/solver.c:898-990`), using a kd-tree range search in code space.
- **Index:** FITS files, each holding a star kd-tree (unit xyz vectors), a quad table and a code
  kd-tree. One file covers one band of quad sizes, in steps of √2. The 4100 series (Tycho-2)
  covers 22′–2000′. The 5200 series covers 2′–22′ and is split into 12 or 48 HEALPix tiles. Small
  scales need tens of GB. Star density is evened out by `uniformize-catalog` (N brightest per
  HEALPix cell, in sweeps) and `hpquads`.
- **Verification:** Bayesian log-odds (`solver/verify.c:448-620`). Each test star is scored either
  as a true match (a Gaussian) or as a spurious "distractor" star (uniform over the field), with
  conflicts resolved. Defaults: distractor ratio 0.25, verify_pix 1 px, accept at odds 1e9
  (`engine.c:672-673`).
- **Output:** a full TAN + SIP WCS from `tweak2` (`solver.c:126-220`), plus `.corr`, `.match` and
  `.rdls` files.
- **Strengths and weaknesses:** the most robust of the three and needs no FOV. It is slow to start
  (seconds), its indexes are large, and the C/FITS stack is heavy.
- **Catalogues:** Tycho-2, 2MASS, USNO-B, NOMAD, UCAC3/4/5 (`catalogs/`). UCAC4 is supported
  (`catalogs/ucac4.c:51-125`) through `ucac4tofits` → `hpsplit` → `build-astrometry-index`.
  **It does not propagate proper motion**; positions stay at J2000.

### 3.2 tetra3 (`esa/tetra3`)

- **Input:** `solve_from_image` uses its own blob extractor (`tetra3.py:1796`).
  `solve_from_centroids` takes **(y, x)**, with the origin at the top-left pixel edge (so the first
  pixel centre is 0.5), sorted brightest first, plus `(height, width)`.
- **Pattern:** 4 stars give 6 edges. They are sorted and divided by the longest edge to give 5
  ratios, binned at `1/(4·max_err)` (`tetra3.py:1060-1071`). The key is a polynomial over the bins,
  multiplied by 2654435761, in an open-addressed table with quadratic probing
  (`tetra3.py:110-153`). It tries every 4-combination of the brightest 8 centroids
  (`tetra3.py:1381`).
- **Database:** an `.npz` with a star table (float32) and a pattern table of star indices (2× table
  oversize), about 99 MB. Built per FOV scale (step 1.5×) with a minimum star separation of
  `0.6·fov/√stars_per_fov` (`tetra3.py:925-1013`).
- **Verification:** rotation by SVD (Kabsch) **with no determinant check, so it can accept a
  mirror-image (wrong parity) solution** (`tetra3.py:236-244`). Projected catalogue stars are
  matched by nearest neighbour and accepted using a binomial false-positive probability
  (`tetra3.py:1589`).
- **Output:** RA, Dec, roll, FOV, one radial distortion term k (`tetra3.py:1641`), RMSE, match count
  and probability. No WCS or SIP.
- **Strengths and weaknesses:** 10–500 ms in Python. Pattern density is uneven across the sky, and
  limiting the search to the brightest 8 stars fails when bright detections are spurious.
- **Catalogues:** BSC5, Hipparcos, Tycho (`tetra3.py:702-870`). Proper motion is applied linearly
  from 1991.25 to the requested epoch, and stars without proper motion are dropped. No UCAC4 or
  Gaia support.

### 3.3 cedar-solve (`smroid/cedar-solve`)

A tetra3 fork with the same API and pixel convention. What changed:

- **Uniform database:** patterns are generated in "lattice fields" on a Fibonacci sphere, with a
  fixed 50 patterns per field, brightest first, in breadth-first combination order
  (`tetra3.py:1182-1308`). The table size is the next prime above 3× the pattern count, with linear
  probing (`tetra3.py:1316-1320`). Each slot stores a 16-bit key hash, which rejects probe
  collisions without loading stars (`tetra3.py:1382`, `2268`), plus a float16 longest edge. The
  magnitude limit is set automatically from the star density needed (`tetra3.py:1118-1143`). The
  database is about 38 MB (about 42k stars).
- **Solve:** tries every image combination breadth-first, so there is no fixed limit of 8 stars
  (`tetra3.py:1819`). Hash bins are visited nearest first (`~1855-1866`). Match tolerance is 0.002,
  the default timeout is 5 s, and solves can be cancelled. It adds P90E/MAXE errors and the rotation
  matrix. The distortion range search was removed.
- **Catalogues:** same as tetra3.
- **Rust:** there is no Rust solver. smroid's **cedar-detect** is a Rust star detector:
  - a row gate followed by a 2-D gate, noise estimation, rejection of hot pixels and trails, and
    2/4/8× binning;
  - about 10 ms per Mpx on a Raspberry Pi 4;
  - takes 8-bit images only and returns centroids with the first pixel centre at 0.5
    (`src/algorithm.rs:509-528`, `958`).

  It could become the CLI's extractor if it were extended to 16-bit images.

### 3.4 Comparison

| | astrometry.net | tetra3 | cedar-solve |
|---|---|---|---|
| Input | image or FITS xylist, 1-based (x, y) | image or (y, x), first pixel centre at 0.5 | same as tetra3 |
| Hash | 4-D quad code, kd-tree | 5 binned edge ratios, open addressing | same + 16-bit key check, linear probing |
| Index | FITS per scale band and HEALPix tile, GBs | npz, ~99 MB | npz, ~38 MB, uniform |
| Verify | Bayesian log-odds, odds 1e9 | binomial false-positive probability | binomial false-positive probability |
| Output | TAN + SIP WCS | RA/Dec/roll/FOV/k | + P90E, MAXE, rotation matrix |
| Speed | 0.1–10 s | 10–500 ms | 10–100 ms |
| Catalogues | Tycho-2, 2MASS, Gaia, USNO-B, UCAC4/5 | BSC5, HIP, TYC | BSC5, HIP, TYC |
| Proper motion | none | linear, drops stars without PM | linear, drops stars without PM |

### 3.5 What Fayde takes from each

The wider literature (star-tracker identification, robust estimation, distortion fitting,
benchmarks) is reviewed in [`docs/research/literature-review.md`](../research/literature-review.md).
Sections 5, 6 and 9 describe what Fayde implements.

1. **Database layout from cedar-solve:**
   - uniform lattice fields with a fixed budget of patterns per field;
   - brightest-first breadth-first combinations;
   - a prime-sized table with linear probing;
   - a 16-bit key check and the longest edge (as `f32`) in each slot.

      A lookup probes every key in the tolerance window (or scans the table when the window holds
      more keys than the table has slots) and ranks the hits by their exact shape error. The index
      is a versioned flat binary, read into memory, rather than npz.
2. **Verification from astrometry.net** (Bayesian log-odds with distractors), behind a cheap
   pre-filter that counts sources landing near catalogue stars (§5.4).
3. **A real WCS output** (TAN, with SIP of order 2–3), plus RA, Dec, pixel scale and roll. Fayde
   fits each seed with an explicit parity and, when the parity is unknown, tries both and verifies
   each, so the parity is right, which tetra3 does not check.
4. **Several indexes for different field widths, like astrometry.net's series**, tried in order
   of closeness to the frame's diagonal (the widest first when the scale is unknown). A pointing
   hint rejects a matched pattern with a star farther from the pointing than the search radius
   plus the frame's half-diagonal.
5. **Pluggable catalogues with proper motion applied to every star.** Each position moves along
   its proper motion on the sphere, which holds at the poles, and stars without proper motion are
   kept with a larger error: all three take shortcuts near the poles or drop such stars.

## 4. Architecture

Cargo workspace:

```
fayde/
  crates/
    fayde-core      sky, pixel and epoch types; portable numerical primitives
    fayde-catalog   Catalog trait, UCAC4/Tycho-2/star-table/.fcat readers, proper motion
    fayde-index     HEALPix selection, pattern hash, streamed index file reads and writes
    fayde-solve     Source, SolveRequest, Outcome and Solution; search, verification, fits and covariance
    fayde-wcs       TAN/SIP projections, WCS header cards, WCS prediction from pointing and roll
    fayde-fits      FITS images and source lists, header priors, profiles, WCS file output
    fayde-detect    star, trail and streak detection in an image
    fayde-api       request building and the JSON report, shared by the command line and Python
    fayde-cli       `fayde index build`, `fayde index info`, `fayde catalog`, `fayde solve`,
                    `fayde verify`, `fayde refit`, `fayde detect`, `fayde calibrate`, `fayde bench`,
                    `fayde compare`
    fayde-py        PyO3 bindings (general Python API)
```

The solver uses `fayde-core`, `fayde-catalog`, `fayde-index` and `fayde-wcs`. FITS, detection
and request orchestration sit outside the solver. `fayde-api` combines those layers for the
CLI and Python. The index reader verifies a checksum in one streamed pass, then parses the
index into memory in a second pass.

Dependencies kept small: `serde` and `serde_json` (reports, profiles), `libm` (portable maths),
`pyo3` (Python bindings), `flate2` (catalogue downloads, in pure Rust). FITS reading and writing,
the k-d tree and the linear algebra are Fayde's own. No C toolchain at runtime.

## 5. Solve pipeline

```
sources + hints
   ├─ validate geometry and options, filter and rank usable centroids (§5.0)
   ├─ plan pair-seed and pattern-hash stages from the hints and indexes (§5.2a)
   ├─ try hypotheses, refit locally, verify, and select the final TAN/SIP fit (§5.3–5.5)
   ├─ try the predicted WCS if the search found no solution and a prediction was supplied (§5.6)
   └─ report the solution, abstention, or error
```

Search stages share the request's deadline and hypothesis limit. The relaxed-prior pass, when
enabled, uses the remaining budget and continues the hypothesis count.

### 5.0 Xylist validation (`fayde-solve`)

A source list can come from any detector, so Fayde checks it before solving:

- Validate the frame dimensions, search priors, and solver policy. The Python adapter also
  checks that source arrays have equal lengths.
- Ignore non-finite centroids and centroids outside the frame. Rank the remaining sources by
  `p_star` and flux when supplied. With no flux, use `p_star` and input order.
- Demote candidates near a better-ranked centroid. Keep truncated sources as match evidence
  and exclude them from seeding and fitting.

### 5.1 Catalogue layer (`fayde-catalog`)

The [`Catalog` trait](../../crates/fayde-catalog/src/lib.rs) exposes fallible cone search and
scans. Queries specify an optional epoch and magnitude limit plus a magnitude band. Catalogue
stars carry a sky position, identifier, optional magnitude, and position errors.

UCAC4 reader:

- Reads the original zone files in place, seeking to the record ranges for each cone. There
  are 900 zones (`z001`–`z900`), each a 0.2° Dec band starting at −90° and sorted by RA, about 113M
  stars and 8.9 GB in total.
- Each record is 78 bytes, little-endian, decoded with explicit little-endian conversion (the field
  layout matches `catalogs/ucac4.c:51-119` in astrometry.net):
  - `ra` u32 (mas) and `spd` u32 (mas from the south pole, so Dec = spd − 90°);
  - `magm`/`maga` i16 (mmag), `sigra`/`sigdc` i8 (+128 mas);
  - `pmrac`/`pmdc` i16 (0.1 mas/yr, where pmRA already includes cos δ), `sigpmr`/`sigpmd` i8 (+128);
  - 2MASS J/H/K, APASS B/V/g/r/i, and `rnm` u32 (the unique star ID).
- **Sentinels:** a magnitude of 20000 means missing. A proper motion of 32767 means the real value
  is in `u4i/u4hpm.dat`, which is loaded once.
- `u4i/u4index.unf` holds, for each of the 900 zones × 1440 RA bins (0.25° each), a start offset and
  a count, so a cone search seeks straight to the rows it needs. The file is optional: a zone
  file is in RA order, so without it the same rows are found by binary search of the zone file.
- **Epoch:** positions are ICRS mean positions at epoch 2000.0. `cepra`/`cepdc` are central epochs,
  not the reference epoch: positions do not move from them, but position errors grow from them.
  Proper motion moves each position to the query epoch as a unit vector on the sphere, the
  zero-parallax, zero-radial-velocity case of Butkevich & Lindegren (2014), which holds at the
  poles. The error on each axis grows as `√(σ₀² + (t − t_c)²·σμ²)` from the central epoch
  `t_c`. Stars without valid
  proper motion keep their 2000.0 position with an inflated position error, rather than being
  dropped as tetra3 does.
- **Magnitude band:** `magm` (579–642 nm, close to unfiltered silicon) by default, falling back to
  `maga`. This is the same default rule as the
  independent reader behind the test oracle, so results can be compared row for row.
- **Bright stars:** UCAC4 is saturated or incomplete brighter than about mag 7–8. A Tycho-2
  index covers the bright stars for wide fields.

Test vectors come from an independent UCAC4 reader
(`crates/fayde-catalog/tests/reader_oracle.rs`).

### 5.2 Index (`fayde-index`)

Each scale band is built offline from a `Catalog` and stored in one `.fidx` file:

- **Star selection:** HEALPix cells at the index epoch. Keep the K brightest stars per cell, in
  a streaming pass, with a minimum separation between stars (default 10″,
  `min_sep_arcsec`). A nearby kept star, including across HEALPix boundaries, rejects the
  candidate unless the candidate is brighter (ties by catalogue ID); a better candidate
  replaces all its fainter neighbours. Rejected stars are not buffered for reconsideration,
  so results are deterministic for a given catalogue stream. This follows astrometry.net's
  `uniformize-catalog` approach. It keeps the density even, so the galactic plane does not swamp
  the hash table and the poles are not empty. The magnitude limit is optional
  (`--mag-limit`); stars without a magnitude are left out.
- **Field widths:** each index is built for one field width (`--fov`). A solve can load several,
  tried in order of closeness to the frame's diagonal (the widest first when the scale is
  unknown, §5.3); [`index-set.md`](../index-set.md)
  lists a recommended set.
- **Pattern hash table** (cedar-solve design):
  - lattice fields with a fixed number of patterns per field, chosen brightest first in
    breadth-first combination order;
  - the key is 5 sorted edge ratios, binned;
  - the table size is a prime above 3× the pattern count, with linear probing;
  - each slot holds a 16-bit key check, a 32-bit floating-point longest edge, and four star indices.
  - lookup visits neighboring keys or checks stored patterns directly when there are fewer
    table slots than keys to visit. Candidates are ranked by exact ratio error.
- **Coverage vs size (measured on a regional UCAC4 index, 0.4° FOV):** with 2× lattice
  oversampling 0 of 70 image 4-sets were in the index; 30× found 11, 100× (cedar-solve's value)
  found 35, and the true match always ranked first. Patterns grow roughly linearly with
  oversampling (431 → 11,645 for the region), so an all-sky index for a sub-degree FOV at 100× is
  large (rough estimate: tens of millions of patterns, gigabytes). This is a reason hinted solves
  should start from pair seeds built at solve time (§5.3), with the stored pattern index reserved
  for loose priors and blind solving.
- **Two tiers:** the first all-sky benchmark solved 89/100 frames blind with
  no wrong solutions, but every failure was a dense field whose true patterns were absent from
  the table (0 of 126–1001 true 4-star sets of the frame's top index stars). Cause: about 40 index
  stars per field of view, but patterns built from each field's 10 brightest, so a frame's
  brightest stars rarely formed a stored pattern. The fix is tetra3/cedar-solve's two densities:
  the dense tier (20 stars per field of view by default) feeds the pattern-seed pre-filter, and
  patterns come only from a sparse pattern tier.
  First attempt (v2): one star per cell of side fov/√10 rounded *down* to a HEALPix level put
  both tiers on nside 256, so half the index became pattern stars and fields were truncated
  again (6 failures in 149 frames, from 8). The v3 selection kept k brightest stars per
  coarse HEALPix cell, creating boundary artifacts: a bright star could lose its slot while a
  nearby, much fainter star in the next cell was promoted. Current (v4): visit dense-tier stars
  brightest-first (ties by catalogue ID), accepting each only if no accepted pattern star is
  within `r_min = 0.6 × fov / sqrt(N)`, with N = `pattern_stars_per_fov` (default 20).
  This is ESA tetra3's minimum-separation thinning (Brown & Stubis 2017 / TETRA): the
  characteristic spacing fov/√N sets an approximate density N/fov², or roughly 10 pattern
  stars in a square frame whose diagonal is fov. The 0.6 factor is empirical; achieved density
  depends on the available dense stars, and sparse catalogues can undershoot the target.
  An incremental Cartesian unit-vector grid checks 27 neighbouring buckets with exact chord
  distances, avoiding cell-boundary and polar artifacts; the existing static k-d tree then
  supplies lattice-field neighbours. Progress reports the number of pattern stars per fov² over
  the requested region (or all sky). Pattern membership remains stored per star, with index
  **file format v2 unchanged**; v4 names the selection method, not the file format.
  Validation: a synthetic uniform all-sky input of 1.6 million dense stars at
  fov = 1° retained 751,083 pattern stars, **18.207 per fov²** (−9.0% from target 20), in
  about **1.8–1.9 s** for sorting and thinning in release mode. An independent k-d-tree check
  verifies separation for every retained star. The coarse-boundary and 3″-double regressions
  both failed before the fixes and pass with v4. Regional UCAC4 validation retains 67 dense
  stars and 9,411 patterns; 5/70 image patterns are present, all ranked first. These v4 results
  are separate from the earlier oversampling measurements above.
- **Quad kd-tree** (astrometry.net design) — *considered, not adopted* (§5.3): 4-D quad codes for blind solving, code tolerance about
  0.01.
- **Epoch policy** (see §10): the index is built at a fixed epoch (2025.0 by default). Stars whose
  proper motion would move them more than half the pattern tolerance at half the field width in 25
  years (about 180 mas/yr per degree of field) are left out of the index, so it stays valid for
  decades; a star with no proper motion counts as still. Pair seeds, verification and fitting do not
  use index positions: they read catalogue stars propagated to the observation epoch at query time
  (§5.1), so they still see those stars. UCAC4 propagation errors at 2026 are about 0.06–0.3″
  (0.02–0.1 px at 3″/px); Gaia matters only for accuracy goals finer than that.
- The header records catalogue ID, epoch, magnitude limit, cell level and build parameters, and
  a trailer holds a content hash. A pattern-hash solution names its index's hash, and the report
  lists the hash of every index loaded.

### 5.2a Adaptive search plan

Fayde's default entry point, `solve_auto`, plans the search from the metadata the request carries
instead of defaulting to a blind search. Stages run in order and share one hypothesis count, so
the false-alarm accounting covers the whole search. Each solution and each abstention returns the
executed plan (`search`: stage, hypotheses, seconds).

| Metadata | Stages |
|---|---|
| Pointing + tight scale (σ/scale ≤ 3%) | Pair seeds within 3σ of the pointing, when 3σ is at most half the hard radius → pair seeds in the whole prior cone → pattern hash in an inner ring (2σ) → pattern hash to the hard radius |
| Pointing, loose or no scale | Pattern hash, inner ring → hard radius (scale-pruned when a scale prior exists) |
| Scale only | Pattern hash over the whole index, scale-pruned |
| Nothing | Blind: pattern hash over every index band |
| `search_outside_pointing` (opt-in) | After the above, pattern hash outside the pointing limit; the solution is flagged `prior_violated` |

The pair seeds look close to the pointing first when its σ is small. The work of a pair-seed
stage grows with the area of its cone, and a pipeline's pointing is often far better than the
hard limit it is willing to set: measured on 146 real
source lists, the pointing was within 2.7′ while the limit was 2°, and a solve cost 8.7 ms at
the 2° limit against 1.1 ms at 0.25°. The first look is skipped when 3σ is more than half the
hard radius, which includes the default σ of a quarter of it: a first look nearly as wide as
the limit costs almost as much as the search it would spare. The σ comes from the request
(`--pointing-sigma`) or the instrument profile, which `fayde calibrate` fills in from solved
frames.

With a pointing and a tight scale the pair seeds go first even when an index is loaded. The
reverse was tried: one hash stage ahead of the pair seeds halves the median solve
on real source lists, but with an all-sky index a stage that finds nothing costs 37 to 395 ms,
against 3 ms for the pair-seed solve it delays, and the mean went from 4.5 ms to 24 ms
(`docs/references.md`, "Considered and not adopted").

The inner ring is skipped when 2σ ≥ 0.8 × hard radius (the second pass would add nothing).
Indexes whose region cannot contain the field are skipped, and implausible priors are rejected
rather than searched: a scale prior that implies a field wider than a hemisphere (header
placeholders such as `CDELT = 1`) is an input error, and `fayde-fits` drops such header scales
with a warning. A sequence stage (predict from the previous frame, §5.6) will go first once
implemented.

### 5.3 Matchers (`fayde-solve`)

- **Source order (PROSAC; Chum & Matas 2005):** image sources are tried in order of
  `p_star` × brightness, not brightness alone, widening the set progressively. In
  frames tracking satellites the brightest detections are often not stars.
- **Hinted, tight scale prior → pair seeds** (k-vector; Mortari): fetch the brightest catalogue
  stars in the hint cone at solve time and build a list of their pair separations. The fainter
  stars that a hypothesis is verified against are read around the hypotheses that pass the
  pre-filter, not for the whole cone, which is tens of frame areas. Each image pair from
  the top of the source order is looked up by separation (tolerance from the scale σ and
  position errors). One matched pair and a choice of parity fix rotation, scale and translation.
  The roll prior can rule out a seed. With half the sources false, a pair is all-real 25%
  of the time and a 4-star pattern only 6%.
- **Hinted, loose scale prior → pattern hash:** build patterns from the ordered image sources and
  look up only hash entries whose stars lie in the cone. The scale range prunes by edge length
  before lookup.
- **No geometric-voting fallback.** Fayde has no voting stage (Kolomenkin et al. 2008), in
  which image pairs vote for catalogue pairs of matching separation. A source list dominated by
  false detections can use up the search budget before a seed finds real stars.
- **Blind:** pattern-hash lookups across the indexes loaded, each built for one field width (as
  astrometry.net's index series), the closest to the frame's diagonal first (the widest when the
  scale is unknown), hypotheses counted across all of them. astrometry.net's 4-D quad codes were
  considered and not adopted: the 5 edge ratios describe the same four-point shape, and the hash
  lookup already gives continuous-tolerance matching (`docs/references.md`, "Considered and not
  adopted").
- Each stage tries candidate transforms in deterministic order. The search returns a hypothesis
  once it passes verification and the final fit.

### 5.4 Verifier

The pre-filter projects the first 12 ranked sources outside the seed and counts those that land
within `3σ + 2 px + 0.01·d` of a star, `d` being the distance in pixels from the seed's centre: at
least 2 among the catalogue stars read for pair seeds, at least 1 among the index's stars for
pattern seeds. Hypotheses that pass go on to the full test.

**Local refit before deciding** (LO-RANSAC; Chum, Matas & Kittler 2003): a seed transform from 2–4
stars is rough. The hypothesis is refitted on its provisional inliers (weighted, §5.5), re-matched,
and only then scored. This rescues correct hypotheses whose seed transform was slightly off.

Bayesian log-odds as in astrometry.net (`solver/verify.c:448-620`, distractor ratio 0.25 by
default): for each image source in rank order (§5.0), the probability
it matches a projected catalogue star (Gaussian with the combined position error) versus being a
distractor (uniform over the field). The default acceptance threshold is ln(10⁹) after the
hypothesis-count correction. Fayde combines the centroid error (`sigma_px` or the policy default)
with catalogue uncertainty, uses `p_star` to adjust the distractor probability, and stretches
the match gate along a known trail for truncated sources. The solve also scores its pointing,
roll and scale priors. The verifier treats each position error as isotropic and doesn't compare
brightnesses.

**Number of false alarms (NFA):** each solution reports
`hypotheses_tested × exp(−log_odds)`. Acceptance uses the log-odds minus the natural logarithm
of the hypothesis count. The fixed-model likelihood argument assumes independent uniform null sources. Refitting on provisional matches uses some of the scored sources to select the model. End-to-end false-accept rates therefore need empirical calibration for each search regime (§9).

### 5.5 Fit (`fayde-solve`)

- Seed transforms use the supplied parity or try both image parities.
- The weighted least-squares TAN fit uses the centroid and catalogue position error at the
  observation epoch.
- **Robust fit on the verifier's likelihood:** iteratively reweighted least squares (IRLS) over
  the same Gaussian-plus-distractor mixture the verifier uses. Each match gets a soft weight from
  its probability of being real; there is no hard kσ clipping. The fit and the verifier therefore
  agree on what an outlier is.
- **Error-scale factor:** a global factor s on the supplied covariances is fitted and reported.
  s ≫ 1 means the producer under-reports its centroid errors, which is useful feedback to
  detector authors.
- **Distortion:** Fayde fits SIP per frame, with the order chosen by BIC. Fayde doesn't fit one
  distortion model across many frames, as a SCAMP-style global solution does (Bernstein et al.
  2017). Among orders up to BIC's choice, the lowest whose weighted mean squared error on held-out
  cells of the frame is within 5% of the best is kept (with matches in fewer than four of the nine
  cells, no distortion; when the order cannot be scored, BIC's choice stands). SIP is also dropped,
  or the hypothesis rejected, when a cell of a 3 × 3 grid with at least 6 sources has confident
  matches at under a quarter of the frame's rate. When the request has a scale or roll prior, a
  solution with distortion is converged again as TAN with scale and roll free of the priors: a TAN
  outside the hint's limits that explains the sources at least as well replaces it and is then
  refused, and a TAN inside them replaces it only with higher posterior odds. Distortion must not
  stand in for a scale or roll the hint held back.
- The solution carries CRVAL, CRPIX at image centre, CD, optional SIP, parameter covariance,
  per-star residuals, weighted RMS, and a leave-one-out RMS when there are enough matches.

### 5.6 Sequence mode

The solver keeps state for a sequence ID: last WCS, timestamps, fitted pointing drift. For frame N
it predicts the WCS, runs §5.4 only, and accepts if it passes. This replaces a pipeline-side
neighbour-hint ladder for the common case and stays separate from WCS *borrowing*: Fayde never
returns a WCS that did not pass verification on the frame itself.

## 6. Input contract

The [Rust request types](../../crates/fayde-solve/src/request.rs) define what the solver
consumes. The [shared API](../../crates/fayde-api/src/lib.rs) resolves header values, instrument
profiles and caller overrides before building that request. See the [FITS source-list
format](../source-list.md) for the file contract.

### 6.1 Principles

Give Fayde the frame dimensions, source centroids, and any measured centroid errors, fluxes,
and star probabilities. Supply uncertainty with a pointing, scale, or roll estimate so the
search can rank hypotheses and apply the intended limits. Fayde records where each prior
came from in the CLI and Python reports.

### 6.2 Frame

| Field | Input |
|---|---|
| `width`, `height` | Image dimensions in pixels. Supply them directly, through the image, or through the source list and its header. |
| `trail` | Optional common star trail. The Rust request takes full length in pixels and angle in radians from +x toward +y. The CLI, Python options and FITS source list take degrees. |
| Observation epoch | Set `Policy::epoch` in Rust or use an epoch override or observation header in the CLI and Python. Without one, file-backed catalogues use J2000.0, and the CLI and Python warn. |

### 6.3 Sources

| `Source` field | Meaning |
|---|---|
| `x`, `y` | Centroid in pixels. Rust and Python default to zero-based pixel centres. FITS source-list positions are one-based. |
| `flux` | Optional integrated flux in linear units, used for ranking. |
| `sigma_px` | Optional isotropic 1σ centroid error in pixels. The policy supplies the default. |
| `p_star` | Optional probability that a detection is a star, from 0 to 1. |
| `truncated` | Whether the frame edge or a mask cuts the source. Fayde uses truncated sources to score a match and excludes them from seeding and fitting. A known trail widens the match gate along its direction. |

The solver ranks centroids itself. Keep uncertain candidates with a low `p_star` when the
upstream detector can estimate one. The verifier takes one isotropic position error per source,
and it doesn't compare source brightness with catalogue magnitudes.

### 6.3a Frame-level detector outputs

The detection report contains its mode, background and noise estimates, and per-source
measurements. In straight-trail mode, the common length and angle accompany the source list.
For a source measured along and across a trail, the solver receives the larger position error
as `sigma_px`. See [Read the detection report](../usage.md#read-the-detection-report).

### 6.4 Priors

| Input | Effect |
|---|---|
| Pointing centre, 1σ error and maximum radius | Set the search region and contribute to the score. |
| Scale in arcsec/px and 1σ error | Select index bands and prune hypotheses. A tight scale and pointing enable pair-seed search. |
| Roll, 1σ error and maximum deviation | Constrain image orientation. |
| Parity | Choose the image handedness. |

### 6.5 Instrument model

The JSON instrument profile can supply a scale and its uncertainty, binning, parity, a roll
rule based on a header keyword, pointing uncertainty and radius, and whether to use the
frame's WCS as a hint. Caller overrides take precedence, followed by the profile and then the
image header. See the [instrument profile example](../../profiles/example.json) and the
[calibration guide](../usage.md#calibrate-an-instrument-profile).

### 6.6 Observer

Fayde doesn't model the observer. It applies no stellar aberration from the observer's velocity,
no light-time correction and no refraction, and it computes no observer state from a TLE or a
ground site. Fayde solves each frame on its own: a solved neighbour's WCS can seed the next frame
(§5.6), but Fayde keeps no state across a sequence and never interpolates a WCS between frames.

### 6.7 Solve policy

| Field | Default | Effect |
|---|---|---|
| `deadline` | None | Starts when the solver receives the request. An accepted hypothesis can finish its final fit. |
| `max_hypotheses` | 5,000,000 | Limits scored hypotheses across the search stages. |
| `accept_log_odds` | ln(10⁹) | Threshold after subtracting the log of the hypothesis count. |
| `default_sigma_px` | 1 px | Position error for a source without `sigma_px`. |
| `max_sip_order` | 3 | Upper limit for SIP model selection. |
| `epoch` | None: J2000.0 for file-backed catalogues | Target epoch for catalogue proper motion. |

### 6.8 Encodings

- **Rust:** construct `SolveRequest` with `frame`, `sources`, `priors`, `policy`, and `previous`.
  Pass the catalogue and any indexes to the solver. `fayde_api::build_request` can resolve
  the header and profile first.
- **Python:** use `Solver.solve` with coordinate sequences, `Solver.solve_image` with a FITS
  path, NumPy array or Astropy HDU, or `Solver.solve_detection` with a detection report.
  Header cards or a keyword dictionary provide frame metadata.
- **CLI:** use `fayde solve --image` with a FITS image or `fayde solve --sources` with a FITS
  binary table. Add `--header`, `--instrument`, or explicit options as needed.

### 6.9 Integrating with a host pipeline

1. Send the detector's source list and the frame header, or pass the image to Fayde's detector.
2. Supply a catalogue path and any indexes needed for the search. In Python, reuse one `Solver`
   across frames so it opens them once.
3. Read the solution or abstention before using the WCS. A solved report includes WCS cards,
   matched pairs and residuals, uncertainty estimates, and the catalogue and index identity.
   The Rust solution also exposes the parameter covariance.
4. Compare solutions against the host's current backend on the same frames. Measure solve
   rate, residuals, and latency before switching the default.

## 7. Error handling

- An exhausted search returns `Outcome::Abstained` with its reason, best log-odds if available,
  hypothesis count, and search trace.
- Errors are for bad input (non-finite hints, zero image size) and unreadable index/catalogue.
- A hypothesis limit or deadline that runs out also returns `Outcome::Abstained`, with that
  reason, the best log-odds seen (no WCS), the hypothesis count and the trace.

## 8. Licensing

Fayde is GPL-3.0-or-later from 0.3.3 (`LICENSE`, `NOTICE`); 0.3.0 to 0.3.2 were Apache-2.0.

- tetra3 and cedar-solve (Apache-2.0): code may be ported with notices kept.
- astrometry.net: mostly BSD-3, which is compatible with notices kept. Some components are GPL.
  Fayde takes no code from them: code under another author's copyleft licence would stop the
  copyright holder from relicensing Fayde. Rule: port from papers (Lang et al. 2010) and BSD files
  only; check the licence header of every file before porting from it.
- Catalogue data: Fayde's releases contain no catalogue data and no built indexes.

## 9. Validation

1. **Unit:** projection round trips, edge-ratio codes, UCAC4 record decode against an
   independent reader on the same zones.
2. **Synthetic:** `fayde bench` solves source lists projected from the catalogue with a known
   WCS, with position noise, false sources and missing stars, and measures the recovery error.
   Every other list holds only random positions, to count false solutions.
3. **Real frames:** `fayde compare` checks each solution against a stored reference WCS on a
   5 × 5 grid of image points.
4. **False-positive guard:** `fayde bench` counts and reports the random lists that solve.
   Whether the false-accept rate stays within what the reported NFA predicts (§5.4) isn't
   measured.

Fayde's results haven't been compared with other solvers on shared public frames.

## 10. Decisions

1. **Catalogue epoch.** The index is built at a fixed epoch, high-proper-motion stars are
   excluded from patterns, and positions are propagated rigorously at query time (§5.2).
2. **Licence.** GPL-3.0-or-later from 0.3.3, with every dependency under a permissive license, so
   the copyright holder can relicense. Versions 0.3.0 to 0.3.2 were Apache-2.0.
