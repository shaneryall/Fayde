# Fayde plate solver — design

Status: the design that versions 0.x implement (written 2026-11-15); where it and the code
differ, the code and `CHANGELOG.md` are right. Prior art (§3) comes from a source-level read of astrometry.net,
tetra3, cedar-solve and cedar-detect. File:line references there are to their upstream repos.

## 1. Purpose

Fayde is a public, open-source astrometry tool written in Rust. It goes from a source list to a
verified WCS with measured accuracy. It must fully meet the needs of a host pipeline, so that the
pipeline can retire `solve-field` and the solving-related code around it.

**Input is an xylist, as with `solve-field`, or an image of point-like stars.** The user runs
whatever source detector suits their sensor, then gives Fayde the resulting list of sources (§6)
along with any header or metadata they have. That remains the contract for frames whose stars
are trailed. Since 2026-11-16 Fayde also takes the image itself when its stars are points (a
frame tracked at the sidereal rate, a spacecraft in fine pointing), detects the stars, keeps
streaks and artefacts out of the list, and fills the same contract:
[`detection.md`](detection.md).

It is not a wrapper. It ports the strongest algorithms from three solvers, and it takes over
capabilities that pipelines otherwise build around a solver, generalized for public use:

| Concern | Taken from | Reason |
|---|---|---|
| Fast hinted match | tetra3 / cedar-solve | 4-star edge-ratio pattern hash, O(1) lookup, millisecond solves |
| Blind / wide fallback | astrometry.net, adapted | Its √2 scale-band index series, searched with Fayde's 5-ratio pattern hash (quad codes considered, not adopted) |
| Match acceptance | astrometry.net | Bayesian log-odds verification against the catalogue |
| Final fit | astrometry.net, extended | TAN + SIP fit, upgraded to weighted iterative refit with outlier rejection |
| Source detection | new, for point-like stars; otherwise the user supplies the xylist | One fixed extractor cannot fit every sensor; astrometry.net's built-in `simplexy` solved about 12% of one host pipeline's frames. Fayde defines a rich xylist contract (§6) that any detector can fill, and fills it itself from an image of point-like stars with a detector that measures the star width, noise and position errors instead of asking for them ([`detection.md`](detection.md)). |
| Catalogue | new | UCAC4 read natively behind a catalogue trait; Gaia later |
| Header priors, sequences, observer state | new | Scale and pointing resolution, neighbour hints and astrometric place, which a pipeline otherwise builds around the solver |

### Standalone first, host pipelines first-class

Fayde is a general-purpose public tool. A host pipeline is a user, not a dependency. Fayde takes
on capabilities that pipelines otherwise build around a solver when they belong in a public
astrometry tool.

**Detection of trailed stars is not one of them.** The host's detector for such frames stays in
the host and writes the xylist. Frames with point-like stars can be given as images.

**Implemented in Fayde** (so that a host pipeline need not build them):

| Capability | Fayde home |
|---|---|
| Source sanitation: dedupe within 3 px, giant-component and linearity filters, brightness ranking | `fayde-solve` xylist validation (§5.0), expressed as `p_star` adjustments rather than hard drops |
| Scale prior from focal length and pixel size, or from CD or PC+CDELT with unit checks; rejection of mixed or invalid WCS | `fayde-fits` |
| Pointing prior from the header, never treating a commanded WCS as solved | `fayde-fits` |
| Exposure time resolution (DATE-OBS, JD-OBS, TIME-OBS, EXPTIME and variants) | `fayde-fits` |
| Annual and orbital aberration, light time | `fayde-ephem` |
| Observer state from a TLE (SGP4) | `fayde-ephem` |
| Neighbour hints, time-ordered sequence solving | `fayde-seq` |
| WCS interpolation across a sequence (borrowing) | `fayde-seq`, always labelled `Interpolated` and never as solved |
| UCAC4 cone search; "is there a catalogue star here?" check | `fayde-catalog` (cone search) + `fayde verify --at` |
| Verifying or repairing existing WCS headers | `fayde verify`, `fayde refit` |
| Writing a FITS WCS atomically while keeping the original header cards | `fayde-fits` (generic provenance cards) |

**Stays in the host pipeline** (pipeline policy, not astrometry):

- WCS cache reuse policy, versioned header cards and policy fingerprints.
- Run-summary and evidence schemas.
- Satellite identification, correlation, stacking, and drift.
- The detector for frames with trailed stars, which writes the xylist Fayde reads.

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
   are required. Every other field (uncertainties, flux, `p_star`, exposure time, observer state,
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

1. Sequence mode: predict frame N from N−1, verify the prediction, fall back to a search only on
   failure.
2. Per-source centroid uncertainty: each source may carry a 2×2 covariance that the verifier and
   fitter use (§6.3). Sources are treated as point centroids; producers with streaked stars
   supply the streak centroid.
3. Photometric consistency: fit a zero point from the matches and use magnitude agreement as a
   second discriminator in the verifier.
4. Aberration-aware output, given observer state (§6.6).
5. Solve the target case (sub-degree FOV, ~3″/px, 30–60 sources) in < 100 ms hinted and
   < 5 s blind on a single core.

### Won't (v1)

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
Its adopted recommendations are folded into §5, §6 and §9.

1. **Database layout from cedar-solve:**
   - uniform lattice fields with a fixed budget of patterns per field;
   - brightest-first breadth-first combinations;
   - a prime-sized table with linear probing;
   - a 16-bit key hash and float16 longest edge in each slot;
   - nearest-bin-first lookup.

   Stored as a versioned, memory-mappable flat binary rather than npz.
2. **Verification from astrometry.net** (Bayesian log-odds with distractors). The binomial test is
   kept only as a cheap pre-filter.
3. **A real WCS output** (TAN, with SIP of order 2–3), plus RA/Dec/roll/FOV. The SVD determinant is
   checked so the parity is right, which tetra3 does not do.
4. **Scale bands like astrometry.net's** (√2 steps) for the blind path, with RA/Dec hints used to
   prune lattice cells before probing the hash.
5. **Pluggable catalogues with proper motion done properly.** Applying proper motion with the full
   3-D space-motion vector fixes the shortcuts all three take near the poles or by dropping stars.

## 4. Architecture

Cargo workspace:

```
fayde/
  crates/
    fayde-core      types: Source, Hint, SkyPos, Wcs, Solution, errors
    fayde-catalog   Catalog trait; UCAC4 reader; epoch propagation; HEALPix cell queries
    fayde-index     index builder + on-disk format (mmap); pattern generation
    fayde-solve     adaptive planner, matchers (pair seeds, pattern hash), verifier, fit
    fayde-wcs       TAN/SIP projection, weighted fit, covariance
    fayde-fits      xylist read (native + astrometry.net), header→priors, WCS + provenance write
    fayde-ephem     time scales, observer state (SGP4/TLE, ground site), aberration, refraction
    fayde-seq       sequence solving, prediction, labelled interpolation
    fayde-cli       `fayde index build`, `fayde solve`, `fayde verify`, `fayde refit`
    fayde-py        PyO3 bindings (general Python API)
```

Dependency direction: `core` ← `catalog`, `wcs`, `ephem` ← `index` ← `solve` ← `seq` ← `fits`,
`cli`, `py`. `fayde-solve` depends on neither `fits` nor `seq`, so the library stays small for
embedded users such as star trackers.

Dependencies kept small: `memmap2`, `rayon` (index build only), `nalgebra` (fits),
`fitsio` or a minimal FITS writer, `pyo3`. No C toolchain at runtime.

## 5. Solve pipeline

```
sources + hints
   │
   ├─ 0. validate      contract + source checks (§5.0); order by p_star × brightness
   ├─ 1. predict       (sequence mode) previous WCS + motion model → candidate WCS
   │                    └─ verify → accept, or continue
   ├─ 2. hinted match  catalogue cone at the hint: pair seeds (tight scale) or
   │                    pattern hash (loose scale); geometric voting as fallback
   ├─ 3. blind match   pattern hash across index scale bands (no priors)
   ├─ 4. verify        local refit on inliers → log-odds + NFA; threshold
   ├─ 5. refine        IRLS fit on the verifier's mixture likelihood; SIP only
   │                    by BIC, on top of the instrument distortion model
   └─ 6. report        WCS + residuals + covariance + error scale + provenance
```

Each stage has its own time budget. The overall budget is a single per-frame deadline, matching
a host pipeline's per-frame solver budget.

### 5.0 Xylist validation (`fayde-solve`)

Fayde never sees the image, so it checks the xylist before solving:

- **Contract checks:** finite values, valid covariances, positions inside the frame, a consistent
  number of rows in every column. A violation is an input error that names the row and field.
- **Source checks:**
  - duplicates within a few pixels;
  - outliers in component size or width (hot pixels and cosmic rays are too narrow, blends too
    large).

  These adjust `p_star` rather than dropping rows, so a poor detector shows up in the diagnostics
  instead of silently degrading the solve.
- **Diagnostics in the solution:** source counts by check, the detection limit Fayde estimated
  when none was given, and which optional fields were missing and what they would have added.
  This tells detector authors what to improve.

### 5.1 Catalogue layer (`fayde-catalog`)

```rust
pub trait Catalog {
    fn id(&self) -> CatalogId;               // name, version, source files hash
    fn cone(&self, centre: SkyPos, radius: Angle, q: &StarQuery) -> Vec<Star>;
}
pub struct StarQuery { pub epoch: Epoch, pub mag_max: f32, pub band: Band }
pub struct Star { pub pos: SkyPos, pub mag: f32, pub pos_err: Option<[f32;2]>, pub id: StarId }
```

UCAC4 reader:

- Reads the original zone files in place, one memory map per zone, with no conversion step. There
  are 900 zones (`z001`–`z900`), each a 0.2° Dec band starting at −90° and sorted by RA, about 113M
  stars and 8.9 GB in total.
- Each record is 78 bytes, little-endian, decoded with explicit little-endian conversion (the field
  layout matches `catalogs/ucac4.c:51-119` in astrometry.net):
  - `ra` u32 (mas) and `spd` u32 (mas from the south pole, so Dec = spd − 90°);
  - `magm`/`maga` i16 (mmag), `sigra`/`sigdc` i8 (+128 mas);
  - `pmrac`/`pmdc` i16 (0.1 mas/yr, where pmRA already includes cos δ), `sigpmr`/`sigpmd` u8;
  - 2MASS J/H/K, APASS B/V/g/r/i, and `rnm` u32 (the unique star ID).
- **Sentinels:** a magnitude of 20000 means missing. A proper motion of 32767 means the real value
  is in `u4i/u4hpm.dat`, which is loaded once.
- `u4i/u4index.unf` holds, for each of the 900 zones × 1440 RA bins (0.25° each), a start offset and
  a count, so a cone search seeks straight to the rows it needs. The file is optional: a zone
  file is in RA order, so without it the same rows are found by binary search of the zone file.
- **Epoch:** positions are ICRS mean positions at epoch 2000.0. `cepra`/`cepdc` are central epochs,
  not the reference epoch, and are not used for propagation. Proper motion is applied to the query
  epoch with the rigorous formulae of Butkevich & Lindegren (2014), including propagation of the
  position covariance; this is correct at the poles. Stars without valid
  proper motion keep their 2000.0 position with an inflated position error, rather than being
  dropped as tetra3 does.
- **Magnitude band:** `magm` (579–642 nm, close to unfiltered silicon) by default, falling back to
  `maga`; APASS V or r can be selected per index. This is the same default rule as the
  independent reader behind the test oracle, so results can be compared row for row.
- **Bright stars:** UCAC4 is saturated or incomplete brighter than about mag 7–8. Indexes merge
  Tycho-2 (or Hipparcos) for bright stars and remove duplicates within about 1″.

Test vectors come from an independent UCAC4 reader
(`crates/fayde-catalog/tests/reader_oracle.rs`).

### 5.2 Index (`fayde-index`)

Built offline from a `Catalog`, stored as one mmap-able file per scale band:

- **Star selection:** HEALPix cells at the index epoch. Keep the K brightest stars per cell, in
  a streaming pass, with a minimum separation between stars (default 10″,
  `min_sep_arcsec`). A nearby kept star, including across HEALPix boundaries, rejects the
  candidate unless the candidate is brighter (ties by catalogue ID); a better candidate
  replaces all its fainter neighbours. Rejected stars are not buffered for reconsideration,
  so results are deterministic for a given catalogue stream. This follows astrometry.net's
  `uniformize-catalog` approach. It keeps the density even, so the galactic plane does not swamp
  the hash table and the poles are not empty. The magnitude limit can be set automatically from the
  density needed, as cedar-solve does.
- **Scale bands:** √2 steps in pattern size, as astrometry.net uses. Bands are built only for the
  instrument classes we actually run, not for every scale over the whole sky.
- **Pattern hash table** (cedar-solve design):
  - lattice fields with a fixed number of patterns per field, chosen brightest first in
    breadth-first combination order;
  - the key is 5 sorted edge ratios, binned;
  - the table size is a prime above 3× the pattern count, with linear probing;
  - each slot holds a 16-bit key check and a float16 longest edge, plus 4 star IDs.
- **Coverage vs size (measured 2026-11-15, regional UCAC4 index, 0.4° FOV):** with 2× lattice
  oversampling 0 of 70 image 4-sets were in the index; 30× found 11, 100× (cedar-solve's value)
  found 35, and the true match always ranked first. Patterns grow roughly linearly with
  oversampling (431 → 11,645 for the region), so an all-sky index for a sub-degree FOV at 100× is
  large (rough estimate: tens of millions of patterns, gigabytes). This is a reason hinted solves
  should start from pair seeds built at solve time (§5.3), with the stored pattern index reserved
  for loose priors and blind solving.
- **Two tiers (fixed 2026-11-15):** the first all-sky benchmark solved 89/100 frames blind with
  no wrong solutions, but every failure was a dense field whose true patterns were absent from
  the table (0 of 126–1001 true 4-star sets of the frame's top index stars). Cause: about 40 index
  stars per field of view, but patterns built from each field's 10 brightest, so a frame's
  brightest stars rarely formed a stored pattern. The fix is tetra3/cedar-solve's two densities:
  the dense tier (~40 per field of view) stays for verification and fitting, and patterns come
  only from a sparse pattern tier.
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
  supplies lattice-field neighbours. Progress reports the actual pattern count per fov² over
  the requested region (or all sky). Pattern membership remains stored per star, with index
  **file format v2 unchanged**; v4 names the selection method, not the file format.
  Validation (2026-11-15): a synthetic uniform all-sky input of 1.6 million dense stars at
  fov = 1° retained 751,083 pattern stars, **18.207 per fov²** (−9.0% from target 20), in
  about **1.8–1.9 s** for sorting and thinning in release mode. An independent k-d-tree check
  verifies separation for every retained star. The coarse-boundary and 3″-double regressions
  both failed before the fixes and pass with v4. Regional UCAC4 validation retains 67 dense
  stars and 9,411 patterns; 5/70 image patterns are present, all ranked first. These v4 results
  are separate from the earlier oversampling measurements above.
- **Quad kd-tree** (astrometry.net design) — *considered, not adopted* (§5.3): 4-D quad codes for blind solving, code tolerance about
  0.01.
- **Epoch policy** (decided, see §11): the index is built at a fixed epoch (e.g. 2025.0). Stars
  whose proper motion would move them more than a fraction of the pattern tolerance over ±25 years
  are excluded from patterns, so the index stays valid for decades. Verification and fitting do
  not use index positions: they fetch catalogue stars propagated to the observation epoch at
  query time (§5.1). UCAC4 propagation errors at 2026 are about 0.02–0.1 px at 3″/px, so Gaia is
  needed only for accuracy goals below about 0.1″.
- Header records catalogue ID, epoch, magnitude limit, cell level, build parameters and a content
  hash. The hash goes into every solution's provenance.

### 5.2a Adaptive search plan (implemented 2026-11-15)

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

The pair seeds look close to the pointing first when its σ is small (added 2026-11-16). The
work of a pair-seed stage grows with the area of its cone, and a spacecraft's pointing is
usually far better than the hard limit a pipeline is willing to set: measured on 146 real
source lists, the pointing was within 2.7′ while the limit was 2°, and a solve cost 8.7 ms at
the 2° limit against 1.1 ms at 0.25°. The first look is skipped when 3σ is more than half the
hard radius, which includes the default σ of a quarter of it: a first look nearly as wide as
the limit costs almost as much as the search it would spare. The σ comes from the request
(`--pointing-sigma`) or the instrument profile, which `fayde calibrate` fills in from solved
frames.

With a pointing and a tight scale the pair seeds go first even when an index is loaded. The
reverse was tried (2026-11-16): one hash stage ahead of the pair seeds halves the median solve
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
  satellite-tracking frames the brightest detections are often not stars.
- **Hinted, tight scale prior → pair seeds** (k-vector; Mortari): fetch the brightest catalogue
  stars in the hint cone at solve time and build a list of their pair separations. The fainter
  stars that a hypothesis is verified against are read around the hypotheses that pass the
  pre-filter, not for the whole cone, which is tens of frame areas. Each image pair from
  the top of the source order is looked up by separation (tolerance from the scale σ and
  position errors). One matched pair plus the parity prior fixes rotation, scale and translation,
  so it is a complete hypothesis. With known roll, single-star seeds suffice. Fewer stars per seed
  is what matters with many false detections: with half the sources false, a pair is all-real 25%
  of the time and a 4-star pattern only 6%.
- **Hinted, loose scale prior → pattern hash:** build patterns from the ordered image sources and
  look up only hash entries whose stars lie in the cone. The scale range prunes by edge length
  before lookup.
- **Hinted fallback → geometric voting** (Kolomenkin et al. 2008): each image pair votes for
  catalogue pairs of matching separation; the stars with the most consistent votes seed
  hypotheses. Used when seeds and hash fail, e.g. many bright false sources.
- **Blind (implemented 2026-11-15):** pattern-hash lookups across a series of indexes in √2
  field-size bands (astrometry.net's index series), best-matched band first, hypotheses counted
  across all bands. astrometry.net's 4-D quad codes were considered and not adopted: the 5-ratio
  edge pattern carries more shape information and the hash lookup already gives
  continuous-tolerance matching (`docs/references.md`, "Considered and not adopted").
- All matchers emit hypotheses `(image stars ↔ catalogue stars, similarity transform)` into a
  verifier queue, ordered by prior probability and seed quality. The first hypothesis over the
  acceptance threshold wins; ties are broken deterministically.

### 5.4 Verifier

Cheap pre-filter first: the binomial false-positive probability that tetra3 uses, computed from a
nearest-neighbour match count. Hypotheses that pass go on to the full test.

**Local refit before deciding** (LO-RANSAC; Chum, Matas & Kittler 2003): a seed transform from 2–4
stars is rough. The hypothesis is refitted on its provisional inliers (weighted, §5.5), re-matched,
and only then scored. This rescues correct hypotheses whose seed transform was slightly off.

Bayesian log-odds as in astrometry.net (`solver/verify.c:448-620`, distractor ratio 0.25 by
default): for each image source in brightness order, the probability
it matches a projected catalogue star (Gaussian with the source's covariance) versus being a
distractor (uniform over the field). Accept at log-odds ≥ ln(10⁹), reject at ≤ −ln(10⁹)
(tunable). Upgrades over the original, all driven by the input contract (§6):

- per-source covariance (`pos_cov`, or the default from `fwhm` and `snr`);
- per-source distractor prior from `p_star`;
- catalogue stars counted only inside `valid_region` and brighter than `detection_limit`;
- a magnitude-consistency term from a fitted zero point;
- a pointing, roll and scale prior term, so the result is a posterior and not just a likelihood.

**Number of false alarms (NFA):** alongside the log-odds, each solution reports
NFA = (hypotheses tested) × P(a false hypothesis scores this well). The acceptance threshold is
applied to the log-odds minus ln(hypotheses tested), so searching harder does not raise the
false-accept rate. The false-accept rate is measured on random source lists (§9).

### 5.5 Fit (`fayde-wcs`)

- The first rotation comes from an SVD (Kabsch) fit. Its determinant sign gives the parity, which
  is checked against the parity hint. tetra3 skips this check and can return a mirror-image
  solution.
- Weighted least-squares TAN fit on the matched pairs. Weights come from each source's covariance
  and the catalogue position error at the observation epoch.
- **Robust fit on the verifier's likelihood:** iteratively reweighted least squares (IRLS) over
  the same Gaussian-plus-distractor mixture the verifier uses. Each match gets a soft weight from
  its probability of being real; there is no hard kσ clipping. The fit and the verifier therefore
  agree on what an outlier is.
- **Error-scale factor:** a global factor s on the supplied covariances is fitted and reported.
  s ≫ 1 means the producer under-reports its centroid errors, which is useful feedback to
  detector authors.
- **Distortion:** fitted once per instrument across many frames (`fayde calibrate`, SCAMP-style
  global solution; Bernstein et al. 2017) and supplied as the instrument distortion prior (§6.5).
  Per-frame SIP is only a small correction on top of that prior, with the order chosen by BIC
  (or leave-one-out residuals), never by default. Fitting order-3 SIP per frame on 30–60 stars
  overfits.
- Output: CRVAL, CRPIX (at image centre, matching `--crpix-center`), CD, SIP A/B (+AP/BP),
  parameter covariance, per-star residuals, RMS, median and 95th-percentile residual.

### 5.6 Sequence mode

The solver keeps state for a sequence ID: last WCS, timestamps, fitted pointing drift. For frame N
it predicts the WCS, runs §5.4 only, and accepts if it passes. This replaces a pipeline-side
neighbour-hint ladder for the common case and stays separate from WCS *borrowing*: Fayde never
returns a WCS that did not pass verification on the frame itself.

## 6. Input contract

The input contract is designed for **solution quality first**. Producers adapt to it, rather
than the other way round. Every field below exists because it improves the
solve rate, rejects false solutions, or makes the solution more accurate. The "Effect" column says
how.

### 6.1 Principles

1. **Every measurement carries its uncertainty.** A position without an error forces the solver
   to guess a weight. Uncertainties drive the verifier's likelihood (§5.4) and the fit weights
   (§5.5).
2. **Prior knowledge comes as a probability distribution plus a hard limit, never as a bare
   bound.** For example, pointing is given as "centre, 1σ = 0.3°, never beyond 8°". The σ
   ranks hypotheses and enters the score. The limit bounds the search. A wrong prior can cost time
   but cannot force a wrong answer.
3. **Only the source list and image size are required.** Each optional field the producer leaves
   out is either estimated by Fayde (for example the zero point or detection limit) or treated as
   unknown. Whichever happened is recorded in the solution's provenance.
4. **Give raw facts, not pre-digested ones.** Send exposure times rather than an epoch, and
   measured fluxes rather than magnitudes. Fayde derives the rest the same way for every producer.
5. **The core crates contain nothing producer-specific.** Each field can carry a free-text
   `origin` tag (for example `"tle_target"` or `"header_commanded"`), which is passed through to
   the provenance.

### 6.2 Frame

| Field | Level | Effect |
|---|---|---|
| `width`, `height` (px) | **required** | Field geometry, edge handling |
| `valid_region`: bitmask (any downsampling) or polygons | recommended | The verifier does not count catalogue stars as missing inside masked areas (bad columns, edges, vignetting, the target itself). Gives sharper log-odds and fewer false rejections. |
| `detection_limit`: flux threshold and/or estimated limiting magnitude, optional completeness map | recommended | The verifier expects only catalogue stars the detector could have seen. This is the largest single improvement to log-odds discrimination on sparse fields. Fayde estimates it if absent. |
| `exposure`: start time with time scale (UTC/TAI/TT), duration | recommended | Proper-motion epoch and aberration. Without it the solve uses epoch 2000.0 and a position error inflated by proper motion. |
| `band`: filter name or effective wavelength range | optional | Chooses the catalogue magnitude band (`magm`, APASS V/r, Gaia G…) so the magnitude-consistency check works |

### 6.3 Sources

Fayde accepts the xylist as **external input**: it never runs or calls the detector. But the
xylist is not limited to `x, y`. It carries everything a good detector measures. The more of it
a producer fills in, the better the solve. Columnar: a struct of arrays in Rust, NumPy arrays in
Python, Arrow/Parquet or a FITS binary table on disk.

**Position**

| Field | Level | Effect |
|---|---|---|
| `x`, `y` (px, centre of the first pixel = 0) | **required** | Centroid. Every source is treated as a point. A producer with streaked stars supplies the streak centroid, which should be its mid-exposure position. |
| `truncated` (bool; xylist column `TRUNCATED`) (**implemented**) | recommended (track-rate) | The source is cut by the frame edge, so its centroid is biased along the trail. It counts as evidence with a gate stretched along `Frame::trail` (1σ = half the trail length; ×4 isotropic when no trail is given), never enters a fit and never seeds. Synthetic test: 17 edge trails biased along a 90 px trail gave 17.96″ grid RMS unflagged, 0.35″ flagged, all 17 still matched. |
| `Frame::trail` (length px, angle) (**implemented**; CLI `--trail-length/--trail-angle`) | recommended (track-rate) | Direction and length of the gate for truncated sources. |
| `pos_cov` (σx², σxy, σy²) | recommended | Fit weights and the verifier's match likelihood. Defaults to an isotropic error from `fwhm` and `snr` when absent. |

**Shape**

| Field | Level | Effect |
|---|---|---|
| `fwhm` (px) | recommended | Default centroid error; rejects cosmic rays and hot pixels, which are too narrow |
| `npix` (pixels in the component) | optional | Rejects giant blended components |
| `bbox` (row/column min/max) | optional | Diagnostics only |

**Photometry**

| Field | Level | Effect |
|---|---|---|
| `flux`, `flux_err` (integrated, ADU) | recommended | Brightness order; the fitted zero point and the magnitude-consistency check in the verifier, which is a second, independent test against false matches |
| `peak` (ADU), `mf_peak` and `mf_snr` (matched-filter response and its signal-to-noise) | recommended | Fallback for brightness order and weighting when integrated flux is missing. |
| `background`, `noise` (local, ADU) | optional | Weights the flux error; marks sources on bright background |
| `saturated_pixels` (count) | optional | Left out of the zero-point fit, with a larger centroid error |

**Classification and provenance**

| Field | Level | Effect |
|---|---|---|
| `p_star` (0–1) | optional | Prior probability the source is a real star. Sets a separate distractor weight per source, so known targets and doubtful detections can be included without risk. |
| `class` (enum) with `p_target`, `p_cosmic`, `p_hot` | optional | A known satellite target becomes masked rather than counted as a missing star |
| `blend` (bool) or `neighbors` (count within 2·fwhm) | optional | Inflated position error |
| `stage` (string: detector pass, e.g. `primary`, `low_count_fallback`, `long_streak`) | optional | Lets Fayde weight passes differently and report solve rate per pass |
| `rank` (brightness rank from the detector) | optional | Tie-breaking and traceability |
| `id` (u64) | optional | Passed through to the matched-pairs output |
| `extra.*` (any further columns) | optional | Carried through unchanged into the diagnostics, so producer-specific fields are never lost |

Fayde takes **every** source. It does its own staging in `p_star` × brightness order (§5.3; the
30/45/60 tiers become internal), so producers should not cap or pre-filter the list. A doubtful source should be
sent with a low `p_star` rather than left out.

### 6.3a Frame-level detector outputs

The detector also learns things about the whole frame. Fayde wants those too:

| Field | Level | Effect |
|---|---|---|
| `psf_fwhm`, `background`, `noise` (frame-wide) | recommended | Defaults for sources that lack them; feeds the estimate of the detection limit |
| `detection_threshold` (SNR and/or ADU) with completeness estimate | recommended | Becomes `detection_limit` (§6.2) |
| `masked_regions`: masks for alarms or targets, bad columns, detector edges | recommended | Becomes `valid_region` (§6.2) |
| `detector`: name, version, config hash | recommended | Recorded in the solution's provenance, so results can be traced to the exact detector run |

### 6.4 Priors

| Field | Level | Effect |
|---|---|---|
| `pointing`: RA/Dec centre, 1σ, hard radius, origin | recommended | Restricts the search to a cone and ranks hypotheses by prior probability. Turns most solves into hinted solves (§5.3). |
| `roll`: angle, 1σ, hard limit, origin (**implemented**: `RollPrior`; CLI `--roll/--roll-sigma/--roll-max`; or an instrument roll rule, §6.5) | optional | Prunes seeds before they are counted and adds a prior term to the score. On sparse synthetic frames (4 real stars + 3 false), roll + parity + a tight scale took 8/10 → 10/10 solved and 493 → 10 hypotheses (`tests/priors_and_accuracy.rs`). |
| `scale`: arcsec/px, 1σ, hard range | recommended | Selects index bands and prunes by pattern edge length |
| `parity`: negative (sky as seen), positive (mirrored), or unknown (**implemented**; CLI `--parity`; instrument profile) | recommended | Halves the hypothesis space and guards against mirror-image solutions |
| `previous`: an earlier solution (WCS + covariance + exposure time) and an optional pointing rate | optional | Sequence mode (§5.6): predict, verify, done |

### 6.5 Instrument model

Reusable across frames. Stored as a separate `InstrumentProfile` (JSON, `fayde-fits::instrument`;
CLI `--instrument`) so standalone users set it up once. Precedence for every prior: command line,
then instrument profile, then frame header. **Implemented 2026-11-15**:
scale with σ, parity, and a roll rule `roll = offset + sign × header[keyword]`. An example
profile is in `profiles/example.json`. On a real frame whose header carries the roll angle, a
calibrated profile cut pair-seed hypotheses from 172 to 2 with an identical solution.

| Field | Level | Effect |
|---|---|---|
| `pixel_pitch` (µm), `focal_length` (mm) with σ, `binning` | recommended | Scale prior when `scale` is not given |
| `distortion`: SIP coefficients from `fayde calibrate` (§5.5), with the calibration date | optional | Fits only small corrections on top of a known distortion. Gives accurate solutions with few stars, where fitting SIP from scratch would overfit. |
| `parity` | optional | Default for the per-frame parity prior |
| `scale_arcsec`, `scale_sigma_arcsec` | recommended | Instrument scale prior; overrides header-derived scales |
| `roll`: `{keyword, offset_deg, sign, sigma_deg, max_dev_deg}` | optional | Per-frame roll prior from a header keyword |

### 6.6 Observer

| Field | Level | Effect |
|---|---|---|
| `observer`: GCRS position (km) and velocity (km/s) at mid-exposure, **or** a ground site (lon, lat, height) | recommended | Fayde computes annual plus orbital aberration. It reports the WCS in either astrometric (catalogue) or apparent coordinates, without folding aberration into an apparent plate-scale error. For a telescope in low Earth orbit the aberration is up to about 23″ (about 7 px at 3″/px). |
| `weather`: temperature, pressure, humidity (ground only) | optional | Refraction, including its variation across wide fields |

### 6.7 Solve policy

| Field | Default | Effect |
|---|---|---|
| `deadline` | none | Wall-clock budget for the whole frame |
| `accept_log_odds` | ln(1e9) | Acceptance threshold (§5.4) |
| `max_sip_order` | 3 | Upper limit; the order actually used is chosen by model selection |
| `output_frame` | `icrs_astrometric` | Or `apparent` (needs `observer`) |
| `output_epoch` | mid-exposure | The epoch catalogue positions are propagated to |
| `indexes` | all loaded | Restrict to named indexes or catalogues |

### 6.8 Encodings

- **Rust:** `SolveRequest { frame, sources, priors, instrument, observer, policy }`, built with a
  builder that checks invariants (finite values, positive σ, a covariance with no negative
  variance).
- **Python (`fayde-py`):** NumPy arrays for sources and dicts or dataclasses for the rest. Zero
  copy where possible.
- **Files (CLI):** `request.toml` plus `sources.parquet` (or Arrow IPC). Importers read
  astrometry.net FITS xylists (1-based, converted on read) and tetra3 `(y, x)` arrays, so existing
  users can switch without changing their tools.

### 6.9 Integrating with a host pipeline

A host pipeline runs its own detector, writes the xylist and calls Fayde. What the detector and
the pipeline should export to give Fayde the best possible input:

| Fayde field | What the host supplies |
|---|---|
| `x`, `y` | Zero-based centroids from the detector's source table |
| `fwhm`, `npix`, `bbox`, `rank` | PSF width, pixel count, bounding box and brightness rank per source |
| `mf_peak`, `mf_snr` | Matched-filter response and its signal-to-noise |
| `flux`, `flux_err` | Integrated flux and its error, for the zero-point fit |
| `stage` | Tag each row with the detector pass that found it |
| `p_star`, `class` | Candidate scores and star/target vetoes, as probabilities instead of filtering upstream |
| `saturated_pixels` | Count per source |
| `masked_regions` → `valid_region` | The mask used during detection |
| `detection_threshold` → `detection_limit` | The threshold actually used |
| `detector` | Name, version and config hash |
| `exposure` | The frame header: `fayde-fits` resolves DATE-OBS / JD-OBS / TIME-OBS and EXPTIME |
| `observer` | For a spacecraft, its own GCRS state at mid-exposure |
| `pointing` | Header CRVAL, a TLE-derived target position or a solved neighbour, each as a σ plus a hard radius according to how far it is trusted |
| `scale` | A σ, not a padded range. With an instrument profile (§6.5): the profile's scale times √(`XBINNING`·`YBINNING`); the profile drops the scale when the header lacks those keywords, so give `--scale` then |
| `roll`, `parity` | An instrument profile (§6.5): a roll rule on a header keyword holding the roll angle, and the instrument's parity |
| `previous` | A solved neighbouring frame's solution, not only its centre |
| `truncated`, `Frame::trail` | Flag sources cut by the frame edge and export the frame trail length/angle |
| `instrument.distortion` | Build one from a calibration run (§9) |
| star/not-star check | `fayde verify --at` / `Solver.stars_at`: nearest stars, combined σ, `p_consistent`, `p_chance` |

What the host consumes back: the TAN/TAN-SIP header, matched pairs with residuals, and index and
catalogue identity with epoch. Fayde returns all of these plus covariance and full provenance.
The host's adapter builds its own cache records and run summaries from them.

Integration path:

1. `fayde-py` exposes `fayde.Solver(ucac4, indexes)` with `solve`, `verify`, `refit` and `stars_at`,
   each returning the CLI's JSON report as a dict (`crates/fayde-py`).
2. The host adds a Fayde backend next to its existing solver backend, behind a config switch.
   Inputs from the table are added one by one; each should give a measurable improvement in the
   corpus A/B (§9).
3. Run both on the existing corpus and compare solve rate, residuals and latency before changing
   the default.

## 7. Error handling

- No solution is a valid result, returned as `Solution::Abstained { reason, best_log_odds, timings }`,
  never an error.
- Errors are for bad input (non-finite hints, zero image size) and unreadable index/catalogue.
- Budget expiry returns the best hypothesis seen and its score, marked not accepted.

## 8. Licensing

Fayde is planned as Apache-2.0 or MIT/Apache dual.

- tetra3 and cedar-solve (Apache-2.0): code may be ported with notices kept.
- astrometry.net: mostly BSD-3, which is compatible with notices kept. Some components are GPL;
  porting from those would make Fayde GPL. Rule: port from papers (Lang et al. 2010) and BSD files
  only; check the licence header of every file before porting from it.
- UCAC4 data: check redistribution terms before shipping built indexes publicly.

## 9. Validation

1. **Unit:** projection round trips, edge-ratio codes, UCAC4 record decode against an
   independent reader on the same zones.
2. **Synthetic:** catalogue-projected fields with known WCS, added noise, false sources and
   missing stars; measure recovery error and false-accept rate.
3. **Corpus A/B:** run Fayde and `solve-field` on the same real frames; compare solve rate,
   residual RMS against UCAC4, latency, and agreement between the two WCS.
4. **False-positive guard:** shuffled or random source lists must never be accepted, and the
   measured false-accept rate must not exceed what the reported NFA predicts (§5.4).
5. **Public benchmarks:** the ESA Kelvins star-tracker scenes (10,000 scenes, up to 80 false
   objects each), scored with the Rijlaarsdam et al. (2020) protocol, as a reproducible regression
   test anyone can run. The open-source LOST star tracker serves as an independent baseline.

## 10. Milestones

| # | Deliverable | Rough size |
|---|---|---|
| 1 | Workspace, core types, TAN projection, FITS WCS writer | 1–2 days |
| 2 | UCAC4 reader + epoch propagation, verified against an independent reader | 2–3 days |
| 3 | Index builder (HEALPix density, pattern hash) | 3–5 days |
| 4 | Hinted matcher (pair seeds + pattern hash), verifier with local refit and NFA, IRLS TAN fit; solves synthetic fields | 1–2 weeks |
| 5 | Python bindings + host-pipeline backend (sources in); corpus A/B | 3–5 days |
| 6 | Blind multi-band search ✅, SIP + BIC + covariance ✅, geometric-voting fallback, `fayde calibrate`, sequence mode | 2–3 weeks |
| 7 | `fayde-fits` header priors + `fayde-ephem` (time, SGP4, aberration), checked against an independent implementation | 1 week |
| 8 | Xylist validation and diagnostics (§5.0); astrometry.net `.xyls` and tetra3 importers; a host pipeline writes the native xylist | 3–5 days |
| 9 | `fayde verify` / `refit`, public packaging (crates, wheels, docs, index build recipes) | 1 week |

## 11. Open questions

1. Instrument classes to index for: sub-degree space-telescope fields only, or also ground
   sensors?
2. ~~Catalogue epoch policy~~ **Decided:** fixed-epoch index, high-proper-motion stars excluded
   from patterns, rigorous propagation at query time (§5.2).
3. Licence choice: Apache-2.0 only, or MIT/Apache dual?
4. Zero-point fit: one free zero point per frame, or a colour term as well when the catalogue
   provides two bands?
5. Should `p_star` priors be calibrated per producer from corpus statistics, or left to the
   producer?
6. Native xylist on disk: Parquet + `request.toml`, or a single FITS file (a binary table plus a
   header with the frame, prior and observer fields) that stays close to astrometry.net's
   `.xyls`?
