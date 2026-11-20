# Changelog

Versioning rules: `docs/versioning.md`.

Versions up to 0.2.1 were development versions and were not published; their notes are kept
because the index files and their hashes date from 0.1.0.

## Unreleased

- Trailed stars from an image. When every star is a straight trail of one length and direction
  (the telescope tracked a satellite or slewed), `fayde detect`, `--image` and
  `fayde.detect`/`solve_image` measure the stars as trails: a trail census decides `points`,
  `trails` or `refused` from two independent estimates and a regional check, an aimed filter
  finds the trails, and one fit per trail gives its position with errors along and across it.
  The tracked object comes out as class `point`. `fayde solve --image` passes the measured
  trail to the solver. Report: `mode`, `census`, `trail`, `sigma_along_px`, `sigma_across_px`,
  class `point`. Trails shorter than 7 FWHM are left to the point detector. On 109 real frames
  whose stars are points the star lists are unchanged; on 84 real frames tracking satellites, 66
  solve (9 before). `--trail-angle` (Python `detect(trail_angle=)`) chooses among trail
  directions; `fayde detect -o` writes the trail into the source list (`TRAILLEN`, `TRAILANG`)
  and `fayde solve --sources` reads it back. Design and measurements:
  `docs/design/streak-detection.md`. For Rust users of `fayde-detect`, `fayde-api` and
  `fayde-fits` the change is breaking: `DetectOptions`, `Detected`, `Detection` and `SourceList`
  gain fields and `Class` gains `Point`.
- Detect stars smeared along a curved path using an empirical template. Field rotation uses
  templates that vary across the image. The detection report adds `kernel` mode and a
  `kernel` field with template dimensions, star count, path length, bend, and region count.
  Image solves in this mode retry with point detections if they find no solution. For Rust
  users, `Detection` gains a field and `FrameMode` gains `Kernel`.
- `--no-trails` (`fayde detect`, `fayde solve`; Python `no_trails`) measures every star as a
  point, skipping the trail census.
- A trailed frame that gives no solution is solved again with its stars measured as points
  (`fayde solve --image`, Python `solve_image`); the report's `retried_as_points` says so.
- The trail census holds up against cosmic rays and background structure. On frames crowded
  with radiation hits, sharp charge deposits are masked and narrow tracks are kept out of the
  direction vote; the second estimate takes its anchors from bright sources, after a finer local
  background and any strong periodic pattern are removed. On simulated trailed frames with
  80–1,000 hits per cm² per second, or a smooth background field of 5–20 ADU on 5 ADU noise,
  the census went from points or refused to trails on all 30, with the direction within
  0.09°; frames that passed before give the same source lists byte for byte.
- The clipped-sky background estimate also applies to frames with trailed stars.
- On frames whose stars all trail, radiation hits no longer set the stars' width. When the round
  peaks are as narrow as hits and the elongated peaks along the frame's trail direction are more
  than twice as wide, detection takes the trails' width across the trail. On an archive of
  46,215 real trailed frames, 511 more frames solved and 41 fewer.
- Trailed frames: a tracked object is reported at its position with its flux; radiation hits
  are flagged by comparing their sharpness with the trails' width across the trail; a hit is no
  longer reported as a tracked object; and roughly collinear trails are no longer traced as a
  satellite streak, whose light is even along its length. On simulated track-rate frames, hits
  flagged went from 0–4 to 120 of 120, false streaks from 5 to 1, and trails found from 212 to
  224. Frames whose stars are points give the same results to the bit.
- Detection of trailed stars is about a third faster, with the same results to the bit.
- Detection does less redundant work, with identical output (every report the same on 46,219
  real frames and 2,540 simulated ones): the periodic-pattern search takes its spectrum's median
  by selection instead of a full sort, kernel templates interpolate without the slopes they do
  not use, and a straight-trail frame skips training a template that cannot override it. 1 to 6%
  faster by frame type.
- `Solver.solve_detection(report)` in Python solves a `fayde.detect` report directly. Its
  positions are FITS 1-based, which `Solver.solve` takes only with `origin=1`; passing them
  without it moved the solution by a pixel.
- `FAYDE_UCAC4` names the catalogue directory when `--ucac4` is left out, and for
  `fayde.Solver()` with no directory.
- Catalogue errors say to download the zone again.
- Image solves report point retry errors.

## 0.3.3 (2026-11-21)

### Added

- The README and `CITATION.cff` give the DOI of the archived releases,
  [10.5281/zenodo.24379728](https://doi.org/10.5281/zenodo.24379728).
- Each release file has a build attestation from GitHub. To check that a download was built
  from the tagged source by Fayde's release workflow, run `gh attestation verify` (see
  `docs/install.md`).

### Changed

- Fayde 0.3.3 and later are licensed under the GNU General Public License, version 3 or later
  (`GPL-3.0-or-later`). Versions 0.3.0 to 0.3.2 stay under the Apache License 2.0. Using Fayde
  is unchanged; distributing a modified Fayde, or software that includes it, now requires the
  GPL's terms. Contributions are accepted under the terms in `CONTRIBUTING.md`.

### Fixed

- `fayde calibrate` takes each frame's header pointing at the centre of the image it solved,
  using the image's size when the header does not give it, as image and source-list solves
  have done since 0.3.2.

## 0.3.2 (2026-11-19)

- Header coordinates whose last field reads exactly 60, such as `RA = '20:29:60.0'` (59.96 s
  rounded up by the writer), are read as the next minute instead of being rejected; an RA of
  exactly 24 hours is 0. Before, such a frame lost its header pointing and, without an index, was
  not searched at all: 80 of 46,219 frames in a real archive. A decimal RA or DEC written with
  padding, such as `RA = ' 12.5'`, is read as degrees; before, the padding made it read as
  hours.
- Distortion is fitted only where the stars support it. A SIP order is kept only if no lower
  order predicts held-out cells of the frame as well. When the hint has a scale or roll, a
  solution with distortion is fitted again as TAN with its scale and roll free; a TAN outside
  the hint's limits that fits the sources at least as well means no solution, and one inside
  them is kept if its posterior odds are higher. Before, a hint whose scale was wrong by more
  than its limits, or matches with errors in one region, could give a solution bent by up to
  146 px in the corners.
- Image and source-list solves use the WCS position of the delivered frame centre, including
  for cropped images. `--header` and `fayde calibrate` read the metadata of image extensions.
  North-pole TAN headers use the FITS default orientation. Eighth- and ninth-order SIP
  solutions write inverse polynomials that Fayde can read back.
- Streak photometry masks bright stars that cross a streak. The report marks affected streaks
  with `photometry_masked`; their flux includes an estimate for the covered ridge, and their SNR
  uses only unmasked samples. For Rust users, `fayde_detect::Streak` gains the field
  `photometry_masked`.
- Distortion tune-ups honour the hypothesis cap, an expired pair-seed search stops before
  reading the catalogue, and relaxed-prior timings include both attempts. Legacy indexes keep
  their verified file hash in the solution provenance. `fayde verify` ignores the request's
  priors, as documented.
- `-o` may not name any input: `fayde detect` refuses its `--image`, and `fayde calibrate`
  its `--solved` list and every frame and WCS the list names. `fayde calibrate -o` writes the
  profile atomically.
- A detection threshold below 1.5σ is refused: it kept noise peaks with no light, reported as
  stars with negative flux.
- When the stars' FWHM is at most one pixel at full resolution, Fayde keeps bright sources
  that are too sharp for the usual shape checks; hot pixels and compact cosmic-ray hits can then
  reach the source list. On simulated frames of FWHM 0.45 and 0.6 px, completeness rose from
  0.23–0.84 to 0.98–1.00 with purity 0.98–1.00. Source lists for 2 px FWHM frames are
  unchanged.
- Detection estimates the sky background of frames that background subtraction clipped at a
  floor, such as zero, by fitting each tile as a left-censored normal. On a simulated frame with
  90% of the sky at the floor, 95 of 579 detections were stars before and 100 of 100 after.
  On 35 real frames with 25–36% of their pixels at zero, all still solve; the fit residual
  fell on 31 (median 0.40 to 0.35 px) and solutions moved by at most 0.14″. Constant padding
  alone leaves the standard background estimate in use.
- The stars' width is no longer taken from a crowd of cosmic-ray hits: in the histogram of
  sizes, a mode at least three quarters as tall as the tallest counts as the stars' when it is
  the widest, dip or not.
- A warning when no observation time is found (catalogue positions at J2000.0). The
  `IMAGEW`/`IMAGEH` cards of `wcs_cards` are whole 80-character cards, so the cards join into a
  header.
- Python: an image passed with `header=None` uses the header of its HDU, and arrays are
  converted to 32-bit floats once instead of through a 64-bit copy.
- Comments that restated the code are removed from the source, scripts and workflows.
- All 109 real point-source frames of the test set solve with 0.3.1 and with 0.3.2.

## 0.3.1 (2026-11-17)

- `pip install fayde` also installs the `fayde` command. It runs the same code as the release
  binary (`fayde-cli` is now a library as well, with `fayde_cli::run`), so a solve finds the same
  solution either way; only the timings in the report differ.
- Solving is faster, with the same results to the last bit. Measured on 109 real frames, one
  thread: a solve with a pointing and a scale takes 6.0 ms on average, from 22.5 ms; with a
  pattern index and a pointing, 2.3 ms from 18.8 ms; with no priors, 1.9 ms from 3.7 ms.
  - The catalogue stars of a pointing prior are read as hypotheses need them. Only the brightest
    stars of the prior's cone are read for seeding; the fainter ones are read around each
    hypothesis that passes the pre-filter. A 2° prior around a 0.6° frame is 46 frame areas, and
    reading all of it was most of a solve. A consequence: a missing zone file is an error when a
    seed or a hypothesis needs it, no longer whenever it lies inside the prior.
  - The UCAC4 cone search works out each star's direction once, not three times, and decides
    the cone from it; the exact separation is computed only for a star at the cone's edge.
  - The pair-seed table is sorted only where a lookup reaches it; sources are matched to
    reference stars on a grid of cells; a fit projects each star from its unit vector.
- A pointing known better than its hard limit is searched close in first. When the pointing's
  1σ (`--pointing-sigma`, or the instrument profile) is at most a sixth of the limit, pair seeds
  look within 3σ of the pointing and widen to the limit only when that finds nothing. Measured
  on 146 real source lists: 8.7 ms per list at a 2° limit, 1.1 ms at 0.25°. With the default σ
  (a quarter of the limit) nothing changes.
- `fayde calibrate` writes an instrument profile from frames Fayde has solved: scale, parity,
  the roll rule for a header keyword, and the pointing's 1σ. A profile needs no other solver:
  solve some frames with their header pointing and a nominal `--scale`, then calibrate.
- Instrument profiles can carry `pointing_sigma_deg`, `pointing_radius_deg` (the hard limit, in
  place of `--radius`) and `header_wcs` (`false`: do not read the pointing or scale from a WCS
  in the header).
- `--no-header-wcs` (Python: `no_header_wcs=True`) reads the pointing from `RA`/`DEC` or
  `OBJCTRA`/`OBJCTDEC` only, for headers whose WCS may be another solver's solution.
- `fayde compare` builds its requests as `fayde solve` does and takes `--instrument`,
  `--no-header-wcs` and `--hinted-only`. It reports how far apart Fayde's solutions with priors
  and blind are; another solver's WCS in the frame file is a second opinion.
- A UCAC4 copy without `u4i/u4index.unf` is searched by RA in each zone file, which is in RA
  order, instead of read whole: a solve that took 350 ms takes 6 ms, as with the index. Keep
  `u4hpm.dat`, which holds the proper motions of the fastest stars.
- New in the library: `Catalog::cone_brightest` (the n brightest stars of a cone; a catalogue
  can override it to rank by magnitude before working out positions), `brightest_first`,
  `Ucac4Catalog::cone_hits_ra_search`, `TanWcs::vector_to_pixel`, `Wcs::vector_to_pixel`,
  `TangentPlane`, `KdTree::any_within`, `Header::without_wcs`, `calibrate_profile` and
  `RequestOptions::no_header_wcs`.

## 0.3.0 (2026-11-16)

The first public release: on GitHub, PyPI and crates.io.

- Giving neither `--sources` nor `--image` (or both) is reported before anything else runs, so
  like every command-line mistake it leaves an existing `-o` file in place. `--help` now says
  which errors come with a JSON report (from the moment the inputs are read) and which are on
  stderr only.

- Source lists are Fayde's own format, named for what they are: `--sources` (and
  `--at-sources`) replace `--xyls` (and `--at-xyls`), which are still accepted. The file is
  unchanged in substance, a FITS binary table with `X`, `Y` and optional columns, so earlier
  files still read; by convention it is named after its frame, `frame.sources.fits`. Lists that
  Fayde writes are labelled `EXTNAME = 'SOURCES'`, and giving an image where a source list is
  expected (or the reverse) is reported as such. In the library, `read_xylist`/`write_xylist`/
  `Xylist` are `read_sources`/`write_sources`/`SourceList`.

- Fayde finds the stars itself when given a FITS image instead of a source list:
  `fayde solve|verify|refit --image frame.fits` (the image's header is used for the priors
  unless `--header` is given), `fayde detect --image frame.fits [-o frame.sources.fits]`, and
  `fayde.detect(...)` / `Solver.solve_image(...)` in Python (a path, a 2-D array or an Astropy
  HDU). For point sources (sidereal tracking or stare). FITS images of every `BITPIX` are read
  (not fpack-compressed ones). Detection options: `--fwhm`, `--detect-threshold`,
  `--saturation`, `--gain`, `--max-sources`; `SATURATE` and `GAIN`/`EGAIN` are read from the
  header.

- The detector (crate `fayde-detect`) is Fayde's own design, `docs/design/detection.md`. It
  needs no settings: it measures the stars' width and matches its filter to it (wide stars are
  binned), and sets its threshold from one noise peak per frame. Each star gets a position, a
  1σ position error and a star probability from one closed-form fit at its peak. Hot pixels
  and cosmic-ray hits are set aside, saturated stars are centred from their clipped region, and
  streaks are traced and reported with their ends, angle, width and flux; nothing on a streak
  is handed to the solver. It is built for difficult frames: the star width is the widest
  distinct size among compact sources (radiation hits can outnumber stars), each peak is tested
  against a ring around it (stray light), and narrowness is judged against the stars' own
  narrow axis (frames smeared by a drift). 2 to 7 ms per megapixel on one core, the same bits
  on every machine. Its work is bounded on hostile input (a width far larger than the frame,
  a saturation level in the background, a frame saturated all over). The detection report gains `classes`, `streaks`, `rejected` (what was set
  aside) and, per source, `p_star` and `class`; `sources` is the list handed to the solver.

- `fayde verify` no longer touches the path given with `-o`. It names the option as unused, as
  before; 0.2.0 and 0.2.1 then went on to delete a file at that path when the check failed or
  produced no solution.
- Solved `.wcs` files name the commit as well as the version in their comment cards.
- Releases are built by GitHub Actions: the `fayde` command and Python wheels for Linux, macOS
  and Windows, every file run on its platform before it is published, with `SHA256SUMS`.
- `THIRD-PARTY-LICENSES.md` holds the licences of the Rust crates compiled into binaries and
  wheels; both carry it.
- Metadata for crates.io and PyPI. The Fayde crates require exactly their own version of each
  other, and each carries `LICENSE`, `NOTICE` and a README.

## 0.2.1 (2026-11-16)

- `-o`: any failed run (an unreadable input, a wrong catalogue directory, a solver error, no
  solution, a failed write) leaves no file at the output path. An error before the solve used to
  leave an earlier run's WCS there, where it would pass for this frame's.
- `-o` naming one of the run's inputs (`--wcs`, `--xyls`, `--header`, `--instrument`,
  `--ucac4`, an `--index`) is refused before anything runs. `fayde refit --wcs a.wcs -o a.wcs`
  with no solution used to delete the WCS it was given.
- `fayde_core::fsutil` is internal (hidden from the documentation, not covered by the
  versioning rules).

## 0.2.0 (2026-11-16)

Fixes from a review with hostile inputs. Several change what Fayde does with input it used to
accept, so this is a minor release (see the versioning rules). Index bytes and content hashes are
unchanged from 0.1.0. On a real 146-frame comparison set the header-prior solve rate and
agreement are identical to 0.1.1.

- No panic on extreme or non-finite WCS and pixel values (`fayde verify --at 1e308,1e308`, a
  `NaN` SIP coefficient): `pixel_to_sky` is defined for every input, and a non-finite header
  number is an error naming the keyword.
- A damaged `u4i/u4index.unf` is reported as such instead of crashing a cone search.
- `fayde index build` checks its parameters and its output path before the catalogue scan.
  Rejected now: a lattice of more than 1e9 fields (`--oversampling 1e12`), a non-finite
  `--mag-limit`, an epoch outside 1800 to 2200.
- A closed stdout pipe (`fayde solve ... | head -1`) ends quietly instead of panicking.
- `CROTA2` (with `CDELT`, no CD or PC matrix) is applied; such a WCS was read unrotated.
- Epochs outside 1800 to 2200 are clock errors: from a header (`MJD-OBS` holding a Julian Date,
  an unset `DATE-OBS`) they are skipped with a warning; from `--epoch` they are an error.
- `-o`: the output path is checked before solving, the report is printed before the file is
  written, a failed write leaves no temporary file, and a run that finds no solution removes a
  file already at that path instead of leaving an earlier WCS there.
- Header pointing: coordinates the header declares as FK4 or of another equinox are not used
  (with a warning); a right ascension out of range is rejected instead of wrapped; the first
  card of a repeated keyword is used, as in cfitsio and Astropy (it was the last).
- Options dropped by `--no-priors` or `--no-pointing`, and options `fayde verify` and
  `fayde refit` do not use, are named in `warnings` or on stderr instead of ignored in silence.
- `--max-sip-order` above 9 is an error (the SIP convention's limit).
- Error messages name the file they are about.
- Temporary files (the `-o` write check, `.wcs` and index writes) are created with
  `O_CREAT|O_EXCL` under unpredictable names: an existing file or symlink at that name is never
  followed or truncated, and a failed write removes only the file Fayde created.
- The example instrument profile is `profiles/example.json`; the catalogue cross-check reads
  `ucac4_oracle.json`.

## 0.1.1 (2026-11-16)

- Licence: Apache-2.0 (was `MIT OR Apache-2.0`, provisional). `LICENSE`, `NOTICE` and
  `CITATION.cff` added; wheels and binary archives carry them.
- Index set built and its content hashes published (`docs/index-set.md`):
  1.0° `1ddf14fa50ad59dd`, 1.3° `4ad74722223054e1`. Index bytes unchanged from 0.1.0: the
  builder identity stored in an index (`fayde_index::BUILDER`) now changes only when the
  index-building code does, not with every crate version, so a patch release keeps every hash.
- Python: `py.typed`, typed reports (`fayde.Report`, `fayde.AtReport`), checked with
  `mypy --strict` in CI.
- `fayde solve|verify|refit --help` and `fayde index build --help` list options with defaults.
- The commit id is recorded in builds from `git archive` source archives (`.git_archival.txt`).
- `fayde compare` reads headers as `fayde solve --header` does.
- Solved `.wcs` files ask for a citation in their comment cards.

## 0.1.0 (2026-11-16)

First release: solve, verify (including `--at`), refit; Python bindings; identical index bytes
on macOS, Linux and Windows; Linux binaries and wheels.
