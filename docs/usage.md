# Use the command line

To install the `fayde` command, see [Install Fayde](install.md). To download the catalog for
solving, see [Get the UCAC4 catalog](catalog.md).

## Solve a frame

To solve an image with the local catalog and check the result, follow these steps:

1. In a terminal, pass the FITS image and the catalog to `fayde solve`:

   ```bash
   fayde solve --image FRAME --ucac4 UCAC4_DIR -o OUTPUT
   ```

   Replace the following:

   - `FRAME`: The path to the FITS image, for example `frame.fits`
   - `UCAC4_DIR`: The path to the directory that contains `u4b`
   - `OUTPUT`: The path to the WCS file to write, for example `frame.wcs`

   The directory of `OUTPUT` must exist, and `OUTPUT` can't be one of the input files. If a
   path contains spaces, enclose it in quotation marks.

   Fayde finds the stars in the image and reads the pointing, pixel scale, observation time,
   and image size from its header. Fayde then prints a JSON report to standard output.

2. Check the report's `solved` field:

   - When `solved` is `true`, Fayde writes the solution to `OUTPUT` as FITS WCS header cards.
   - When `solved` is `false`, Fayde writes no solution, and the report's `reason` field says
     why.

   Fayde returns an exit status of 0 when the run completes, with or without a solution. If the
   run fails, Fayde returns an exit status of 1 and prints the error to standard error. If the
   solver fails after Fayde reads the inputs, for example because a zone file that the solve
   needs is missing, Fayde also prints a JSON report to standard output. That report has
   `error` set to `true` and the error in `reason`.

   For more information about the report fields, see [Read the report](#read-the-report).

When a run writes no solution, Fayde removes any file that already exists at `OUTPUT`. Fayde
keeps the file if it rejects the command before the run starts, if its initial check finds
that it can't write to `OUTPUT`, or if `OUTPUT` is one of the input files.

### Set the catalog directory once

To leave out `--ucac4`, set the `FAYDE_UCAC4` environment variable to the path to the catalog
directory. Each `fayde` command that has a `--ucac4` option uses the variable when you leave out
the option. A `--ucac4` option on the command line takes precedence. In Python, `fayde.Solver()`
without a catalog directory also uses the variable.

On Linux or macOS, set the variable for the current terminal session:

```bash
export FAYDE_UCAC4="UCAC4_DIR"
```

On Windows, set the variable for the current PowerShell session:

```powershell
$env:FAYDE_UCAC4 = "UCAC4_DIR"
```

Replace `UCAC4_DIR` with the path to the directory that contains `u4b`.

### Use your own source list

If you detect the stars yourself, pass a source list with `--sources` in place of the image. A
source list is a FITS binary table with `X` and `Y` columns. For the columns that Fayde reads,
see the [source list format](source-list.md).

To solve a source list, pass it together with the FITS file that has the frame's header:

```bash
fayde solve --sources SOURCES --header FRAME --ucac4 UCAC4_DIR -o OUTPUT
```

Replace `SOURCES` with the path to the source list, for example `frame.sources.fits`. Fayde uses
only the header from the file at `FRAME`.

### Detect stars without solving

To find the stars in an image and write them as a source list, use `fayde detect`:

```bash
fayde detect --image FRAME -o SOURCES
```

Fayde prints a JSON report with the background, the noise, the width of the stars, and every source
that it found. The report also lists the streaks that cross the image, such as satellite streaks,
and the hot pixels and cosmic-ray hits that Fayde set aside. With `-o`, Fayde also writes the
sources to `SOURCES`, in the format that `--sources` reads. For more information about the report
fields, see [Read the detection report](#read-the-detection-report).

When Fayde measures stars as points and their full width at half maximum (FWHM) is about
1 pixel or less at full resolution, it keeps sources that are too sharp for the usual shape
checks. A star this narrow can put nearly all its light in one pixel, so the source list can
also include hot pixels and compact cosmic-ray hits.

If your instrument profile has a `charge_transfer` field, add `--instrument PROFILE` so that
the report and `SOURCES` hold corrected positions. For more information, see
[Correct charge-transfer smear](#correct-charge-transfer-smear).

To list the detection options and their default values, run `fayde detect --help`. The same
options work with `fayde solve --image`. The `--detect-threshold` option sets the peak threshold
in σ of the filtered image, and its value must be at least 1.5.

If background subtraction clipped many sky pixels to the same minimum value, such as
zero, Fayde accounts for the clipping when estimating the background and setting detection
thresholds. The clipped-sky model requires floor pixels mixed with sky noise across most
of the image's sky tiles. Frames with a constant border and ordinary sky use the standard
background estimate.

Tracking a satellite, slewing, or rolling a camera can smear every star during an exposure.
Fayde checks each image for this motion and reports its decision in `mode`, which has one of
the following values:

- `points`: The stars are points, or trails too short to measure as trails.
- `trails`: Every star is a straight trail of the same length and direction. The `trail` field
  gives the length and direction, and each source has a position error along the trail
  (`sigma_along_px`) and across it (`sigma_across_px`). Point-like sources, such as the tracked
  object, have the class `point` and aren't used to solve the frame.
- `kernel`: Fayde measures the stars with a template learned from the image. A shared
  curved path uses one template. Field rotation uses templates that vary with position in
  the image. Each reported position is the center of the star's light over the exposure.
- `refused`: Fayde finds evidence of trails but can't establish a consistent model, so it
  measures the stars as points. The `census` field gives the reason.

To measure every star as a point, add `--no-trails` to `fayde detect` or `fayde solve`.

If an image measured in `trails` or `kernel` mode gives no solution before the time and hypothesis
limits are reached, `fayde solve` retries with point detections. The `retried_as_points` field
records the retry outcome in `solved` and the original census in `trails_census`. If the retry fails
to read the catalog, the report carries the error and its reason, plus the first attempt's search
trace and hypothesis count.

The `fayde solve --image` command passes a measured trail to the solver. With `-o`,
`fayde detect` records the measured trail in the source list, and `fayde solve --sources`
passes that trail to the solver. For more information about the keywords that hold the trail,
see [Trail keywords](source-list.md#trail-keywords).

If an image shows star trails in more than one direction, give the expected direction with
`--trail-angle`, in degrees from +x toward +y. Fayde uses the angle to choose a direction and
measures the trail from the image. The `fayde detect` command accepts `--trail-angle` alone. It
reports an error if you give `--trail-length` without `--trail-angle`.

With `fayde solve`, `fayde verify`, or `fayde refit`, give `--trail-length` in pixels and
`--trail-angle` together, or neither. With `--image`, the angle also guides detection. The two
values replace the measured or recorded trail in the solve.

### Detect streaks without a catalog

To find streaks in an image and save the detection report to a file, use `fayde detect`:

```bash
fayde detect --image frame.fits > frame.detections.json
```

Read the JSON report's `streaks` array for straight streak candidates, longest first. This
command works without a star catalog or WCS. The optional `-o` FITS output contains the
stars used for solving, so read streak measurements from the JSON report.

Each candidate has these fields.

| Fields | Meaning |
|---|---|
| `x0`, `y0`, `x1`, `y1` | Estimated visible endpoints in FITS 1-based pixel coordinates |
| `length_px`, `angle_deg` | Visible length in pixels and line orientation from +x toward +y, modulo 180° |
| `width_sigma_px`, `flux` | Cross-streak Gaussian width in pixels and integrated signal above the background |
| `photometry_masked` | A recognized crossing star was masked. Flux estimates the missing samples at the median ridge height. SNR uses the unmasked samples only. |
| `snr`, `coverage` | Integrated signal-to-noise ratio of the unmasked ridge and fraction of the visible length above the tracing threshold |
| `touches_edge`, `cosmic_ray` | An endpoint reaches the image edge, or the candidate is narrow enough to be flagged as a possible detector track, respectively |

The tracer starts from a detectable elongated peak, so it can miss a faint streak that becomes
visible only when integrated along its length. In trailed-star frames, Fayde
uses the shared stellar trails as astrometric sources and reports the remaining streak
candidates in `streaks`. Treat `cosmic_ray` as a shape flag and the endpoint order as a
geometric convention. Object identity and direction of motion require other evidence.

### Specify the pointing and scale

Command-line values take precedence over the instrument profile, which takes precedence over
the frame header. To solve an image whose header lacks the pointing or the pixel scale, specify
them:

```bash
fayde solve --image FRAME --ra RA --dec DEC --scale SCALE --ucac4 UCAC4_DIR -o OUTPUT
```

Replace the following:

- `RA` and `DEC`: The right ascension and declination of the image center in degrees
- `SCALE`: The pixel scale in arcseconds per pixel

When you solve a source list without a frame header, Fayde also needs the image size. Add
`--width` and `--height`, in pixels, unless the source list has the `IMAGEW` and `IMAGEH`
keywords.

To decide whether you need these options, check the frame header. Fayde reads keywords such as
the following:

| Information | Header keywords |
|---|---|
| Pointing | Image center projected through a supported header WCS, otherwise its reference coordinates (`CRVAL1`, `CRVAL2`), `RA` and `DEC`, or `OBJCTRA` and `OBJCTDEC` |
| Pixel scale | `XPIXSZ` with `FOCALLEN`, the `CD` matrix, `CDELT1` and `CDELT2`, `PIXSCALE`, or `SECPIX` |
| Observation time | `JD-OBS`, `MJD-OBS`, or `DATE-OBS`, with `EXPTIME` or `EXPOSURE` |
| Image size | `NAXIS1` and `NAXIS2` |

If the header has no observation time and you don't set the epoch, Fayde uses the catalog
positions at epoch J2000.0 and adds a warning that starts with `no observation time` to the
report. To set the epoch, use `--epoch`.

Fayde searches within 2° of the pointing. To change the limit, use `--radius`. To specify the
pointing accuracy, use `--pointing-sigma` with the 1σ error in degrees. When that error is at
most a sixth of the limit and the pixel scale is known to within 3%, Fayde first searches
within 3σ of the pointing. It widens the search to the limit only when it finds no solution
there.

To ignore the pointing and pixel scale from a WCS in the frame header, for example a WCS that
another plate solver wrote, add `--no-header-wcs`. Fayde then reads the header pointing only
from `RA` and `DEC`, or from `OBJCTRA` and `OBJCTDEC`.

To list every option and its default value, run `fayde solve --help`.

### Use an instrument profile

An instrument profile is a JSON file with settings that frames from one instrument share:
the pixel scale, parity, a rule for deriving the roll from a header keyword, and the pointing
accuracy in the frame headers. To use a profile, add `--instrument`:

```bash
fayde solve --image FRAME --instrument PROFILE --ucac4 UCAC4_DIR -o OUTPUT
```

Replace `PROFILE` with the path to the profile file.

All of the following profile fields are optional:

| Field | Meaning |
|---|---|
| `name` | A name for the profile. The report shows it in `priors.instrument`. |
| `scale_arcsec` | The pixel scale of an unbinned frame, in arcseconds per pixel |
| `scale_sigma_arcsec` | The 1σ uncertainty of `scale_arcsec`, in arcseconds per pixel. The default is 0.5% of `scale_arcsec`. |
| `parity` | `negative` for an image that shows the sky as seen, with east to the left when north is up, or `positive` for a mirrored image |
| `roll` | The rule that gives the roll from a header keyword, with these fields: `keyword`, the header keyword that holds the instrument's roll in degrees; `offset_deg` and `sign`, which give the roll as `offset_deg` + `sign` × the keyword's value; `sigma_deg`, the 1σ error of that roll; and `max_dev_deg`, the hard limit on its deviation. `sign` is 1 or -1, and the angles are in degrees. |
| `binning` | The header keywords that hold the binning factors, with these fields: `x_keyword` along x and `y_keyword` along y. Fayde multiplies the scale and its uncertainty by the square root of the product of the two factors. Without `y_keyword`, Fayde uses the x factor for both axes. |
| `pointing_sigma_deg` | The 1σ error of the pointing in the frame headers, for each axis, in degrees |
| `pointing_radius_deg` | The search limit around the pointing, in degrees, in place of the default of 2° |
| `header_wcs` | When `false`, has the same effect as `--no-header-wcs` on every frame |
| `charge_transfer` | The CCD's charge-transfer smear, which Fayde removes from the positions it detects. For more information, see [Correct charge-transfer smear](#correct-charge-transfer-smear). |
| `periodic_pattern` | Stripes from electronic pickup during readout, which Fayde subtracts before it detects stars. For more information, see [Remove a periodic background pattern](#remove-a-periodic-background-pattern). |
| `source` | A description of where the values come from. Fayde doesn't read this field. |

The roll is the position angle of the image +y axis, measured east of north:
`atan2(CD1_2, CD2_2)`. The report's `roll_deg` field uses the same convention.

The repository contains an example profile, `profiles/example.json`. A release archive has the
same file in the extracted archive's `profiles` directory. The example's values are
illustrative. To create a profile for your instrument, follow these steps:

1. Copy `profiles/example.json`.
2. Edit the copy with your instrument's values. The `source` field of the example explains the
   scale and roll values.
3. If your frame headers lack the keywords that `binning` names, delete `binning` from the copy.

   Without those keywords, Fayde doesn't use the profile's pixel scale for the frame.

4. If your frame headers lack the keyword that `roll` names, delete `roll` from the copy.

   Without that keyword, Fayde doesn't use the profile's roll for the frame.

If you don't know your instrument's values, use Fayde to measure them from its solved frames.

### Calibrate an instrument profile

You don't need another plate solver to calibrate a profile. For the initial solves, you need
the pointing in the frame headers and an approximate pixel scale. To create a profile from
these solutions, follow these steps:

1. In a terminal, solve at least five frames with an approximate pixel scale:

   ```bash
   fayde solve --image FRAME --scale SCALE --ucac4 UCAC4_DIR -o WCS
   ```

   Replace the following:

   - `SCALE`: The pixel scale in arcseconds per pixel, correct to within 5%
   - `WCS`: The path to the WCS file to write, for example `frame.wcs`. Use a different path
     for each frame.

   For each frame whose report has `solved` set to `true`, Fayde writes the solution to `WCS`.
   More frames give a better profile.

2. In a text file, list only the frames whose report has `solved` set to `true`, with one frame
   per line. Each line contains the path to the FITS image, a tab character, and the path to
   its WCS file, in that order.

3. In a terminal, create the profile from the list with `fayde calibrate`:

   ```bash
   fayde calibrate --solved LIST --roll-keyword ROLL_KEYWORD --binning XBINNING,YBINNING -o PROFILE
   ```

   Replace the following:

   - `LIST`: The path to the text file from the previous step
   - `ROLL_KEYWORD`: The header keyword that holds the instrument's roll angle in degrees
   - `PROFILE`: The path to the profile to write, for example `camera.json`

   If the headers have no roll keyword, omit `--roll-keyword`. If the headers don't record the
   binning in `XBINNING` and `YBINNING`, omit `--binning` or name the keywords that they use.

   Fayde prints the profile and writes it to `PROFILE`.

A calibrated profile has the following fields:

| Field | Meaning |
|---|---|
| `scale_arcsec`, `scale_sigma_arcsec` | The median pixel scale of the solutions and its 1σ uncertainty, in arcseconds per pixel. With `--binning`, this is the scale of an unbinned frame. |
| `parity` | The parity of the solutions. Fayde omits this field when the solutions disagree. |
| `roll` | With `--roll-keyword`, the rule that gives the roll from the roll keyword, with its 1σ error and hard limit in degrees |
| `binning` | With `--binning`, the binning keywords that you named |
| `pointing_sigma_deg` | The 1σ error of the pointing in the headers, in degrees. Fayde omits this field when fewer than five headers have a pointing. |
| `charge_transfer` | With `--charge-transfer`, the CCD's charge-transfer smear. Fayde omits this field when the frames show no smear. |
| `calibration` | The numbers behind these values: frames used, spreads, and largest deviations |

You can add the `pointing_radius_deg` and `header_wcs` fields manually. To include
`header_wcs` in the profile during calibration, add `--no-header-wcs` to the `fayde calibrate`
command.

### Correct charge-transfer smear

A CCD moves each pixel's charge across the chip to its readout. Traps in the silicon hold back
part of the charge and release it later, behind the rest. A star's measured position therefore
moves away from the readout edge. The shift is larger for stars farther from that edge, and
larger for faint stars than for bright ones. The WCS fit can't absorb a shift that depends on
brightness, so it remains in the residuals. Radiation damage adds traps over time, so the
smear can grow over the years.

To measure the smear, add `--charge-transfer` to the `fayde calibrate` command in
[Calibrate an instrument profile](#calibrate-an-instrument-profile):

```bash
fayde calibrate --solved LIST --charge-transfer --ucac4 UCAC4_DIR -o PROFILE
```

Each line of `LIST` must name the FITS image, not only a header, because Fayde detects the
stars again. Fayde matches them to the catalog under each frame's solution, and fits how their
residuals depend on brightness and on distance from each edge of the frame. Fayde keeps the
result only if the shift is positive and at least five times its uncertainty. If it isn't,
Fayde prints a warning and omits the field. Calibration needs at least 200 matched stars.

The `charge_transfer` field has these fields:

| Field | Meaning |
|---|---|
| `axis` | `y` if the charge moves along the columns to a row of the frame, or `x` if it moves along the rows to a column |
| `readout_edge` | `low` if the charge moves toward row or column 1, or `high` if it moves toward the last row or column |
| `shift_px` | The shift in pixels per 1,000 transfers per factor of e (2.718) in flux below `flux_ref` |
| `flux_ref` | The flux, in image units, at and above which Fayde doesn't move a star |
| `per_year`, `epoch` | The change of `shift_px` per year, and the Julian year at which `shift_px` holds. Fayde fits these when at least five dated frames, each with at least six matched stars, span at least 1 year. Frames without an observation time then don't take part. |

For a frame observed in Julian year `t`, the shift is
`s` = `shift_px` + `per_year` × (`t` − `epoch`), and never less than 0. Without an observation
time, Fayde uses `shift_px`. A star `n` rows or columns from the readout edge, with flux `F`,
moves away from that edge by `s` × (`n` / 1000) × ln(`flux_ref` / `F`). A star with a flux of at
least `flux_ref` doesn't move. With `readout_edge` set to `high`, Fayde counts `n` from the last
row or column, and the star moves toward row or column 1. Fayde removes that shift from every
star that it detects: with `fayde solve`, `fayde verify`, or `fayde refit` and `--image`, and
with `fayde detect --instrument PROFILE`. Fayde doesn't change the positions in a source list
that you give with `--sources`, because it calibrates the term against the flux that its own
detection measures.

### Remove a periodic background pattern

Some cameras pick up electronic interference while they read out a frame. The interference adds
a wave of fixed frequency to every frame, which appears as fine stripes. Stripes a fraction of
the noise high raise the measured noise, so faint stars go undetected. They also line up noise
peaks along the stripes, which can make Fayde take a frame of trailed stars for something else.

To find the stripes, add `--periodic-pattern` to the `fayde calibrate` command in
[Calibrate an instrument profile](#calibrate-an-instrument-profile):

```bash
fayde calibrate --solved LIST --periodic-pattern -o PROFILE
```

Each line of `LIST` must name the FITS image. Fayde averages the power spectra of the frames'
backgrounds and keeps up to four peaks that stand far above the rest. If it finds none, it
prints a warning and omits the field. Otherwise the profile gets a `periodic_pattern` field with
a `frequencies` list. Each frequency is a pair of numbers, in cycles per pixel along x and y.

With such a profile, Fayde fits each wave's amplitude and phase in 32 × 32 pixel tiles of every
frame, away from the stars, and subtracts the waves before it detects stars. The detection
report's `periodic_pattern` field gives the frequencies, the RMS of what Fayde subtracted
(`rms`), and each wave's mean amplitude (`amplitude`), in image units. Without such a profile,
`periodic_pattern` is `null`. If you give `--periodic-pattern` and `--charge-transfer`
together, Fayde subtracts the pattern before it measures the charge-transfer smear.

A frame that is binned differently can have the stripes at another frequency. At a frequency
with no stripes, the fitted amplitude is near zero and the frame barely changes.

### Solve without a pointing or scale

When you know the pointing and have a precise pixel scale, Fayde needs only the catalog. In
that case, omit `--index`. Fayde searches a pattern index only if it finds no solution around
the pointing. This extra search delays the report when a frame has no solution. When the
pointing or scale is uncertain or unknown, Fayde also needs a pattern index, which you build
once from the catalog. To build the recommended indexes, see
[Build pattern indexes](index-set.md).

To solve with a pattern index and ignore all priors, use `--index` and `--no-priors`:

```bash
fayde solve --image FRAME --ucac4 UCAC4_DIR --index INDEX --no-priors -o OUTPUT
```

Replace `INDEX` with the path to a pattern index, for example `allsky-1.0.fidx`. To use
several indexes, add `--index` once for each file.

### Limit the search time

Search time increases with the area within the pointing radius. The following measurements
show exhausted searches of real frames on a 32-core x86_64 Linux computer, with the pixel
scale known to 0.1% and the roll unconstrained:

| Pointing radius | Hypotheses until exhausted | Time |
|---|---|---|
| 0.5° | 250–760 | 0.01–0.02 s |
| 1° | 830–2,500 | 0.02–0.05 s |
| 2° | 2,900–10,500 | 0.07–0.16 s |
| 4° | 11,000–40,000 | 0.24–0.62 s |
| 8° | 47,000–157,000 | 0.9–2.3 s |

A search with a 1° radius takes about 50 times less work than one with an 8° radius. Frames
with hundreds of sources took more than 18 s at an 8° radius. To limit a search, use these
options:

- Use `--max-hypotheses` to limit the work independently of the computer's speed.
- Use `--deadline` to limit elapsed search time. For image input the clock starts after the
  initial source detection. A retry with point detections shares the remaining time and
  hypothesis budget. Fayde checks the deadline again after measuring the points.

Each `search` stage in the report ends with an `outcome`: `solved`, `exhausted`, `deadline`, or
`hypothesis limit`. When Fayde finds no solution, the report names the stage that a limit stopped.

## Read the report

The `fayde solve` and `fayde refit` commands use the same report fields. The fields present
depend on the outcome:

- With a solution, the report has `solved` set to `true` and the fields in the following table,
  except `reason`, `error`, `best_log_odds`, `warnings`, and `detection`.
- When a search finishes without a solution, the report has `solved` set to `false`, `reason`,
  `best_log_odds`, `hypotheses_tested`, `search`, `priors`, and `identity`.
- After a solver error, the report has `solved` set to `false`, `error` set to `true`, `reason`,
  `priors`, and `identity`.

Each report also has `warnings` when there are warnings, and `detection` when you give an image
with `--image`. Pixel positions in the report are 1-based, as in FITS.

To state the uncertainty of a solution, quote `worst_sigma_arcsec`. It's the largest predicted
position error at the sampled positions across the frame.

| Field | Meaning |
|---|---|
| `solved` | `true` when Fayde found a solution that passed verification |
| `reason` | Why Fayde found no solution, or the error |
| `error` | `true` when the solver failed |
| `crval_deg` | The right ascension and declination of the image center, in degrees. The WCS has its reference pixel at the image center. |
| `pixel_scale_arcsec` | The pixel scale at the image center, in arcseconds per pixel |
| `roll_deg` | The position angle of the image +y axis, measured east of north, in degrees from −180° to 180°: `atan2(CD1_2, CD2_2)` |
| `parity` | `Negative` for an image that shows the sky as seen, with east to the left when north is up, or `Positive` for a mirrored image |
| `centre_sigma_arcsec` | The predicted radial 1σ position error at the image center, in arcseconds |
| `worst_sigma_arcsec` | The largest predicted radial 1σ position error on a 5 × 5 grid that covers the frame, corners included, in arcseconds |
| `rms_px` | The root mean square (RMS) residual of the matched sources, weighted by their inlier probability, in pixels |
| `loo_rms_px` | The leave-one-out RMS residual of the confident matches, in pixels: each match's residual from a fit without it. The value is `null` when there are too few matches. |
| `error_scale` | The factor that the fit applies to the source position errors. A value greater than 1 means that the errors were too small. |
| `matches` | The number of matched sources |
| `log_odds` | The natural logarithm of the odds of the solution against sources at random positions |
| `best_log_odds` | Without a solution, the highest `log_odds` that a hypothesis reached, or `null` |
| `nfa` | Model-based false-alarm score, `hypotheses_tested × exp(−log_odds)`. Calibration of the adaptive search requires null trials. |
| `method` | How Fayde found the solution: `PairSeeds` from star pairs around the pointing, `PatternHash` from a pattern index, or `Refit` from the WCS that you gave |
| `index_hash` | The content hash of the pattern index that found the solution, or `null` |
| `solved_by_stage` | The name of the search stage that found the solution |
| `sip_order` | The order of the SIP distortion polynomial, or 1 when the solution has no distortion terms |
| `priors` | The priors that Fayde used, and in each prior's `origin` field, where the prior came from: the command line, the instrument profile, or a header keyword. It also gives the instrument profile's name, the trail, the number of truncated sources, and the default position error. |
| `warnings` | Header values that Fayde didn't use, and options that `--no-priors` or `--no-pointing` dropped. When the header has no observation time and you don't set the epoch, a warning starts with `no observation time`, and Fayde uses the catalog positions at epoch J2000.0. |
| `prior_violated`, `violated_priors` | Whether the solution breaks a hard limit of the priors, and which: `pointing`, `scale`, `roll`, or `parity`. This can happen only with `--search-outside-pointing` or `--relax-priors-on-failure`. |
| `relaxed_priors` | `true` when the solution comes from the retry with wider limits that `--relax-priors-on-failure` allows |
| `hypotheses_tested` | The number of hypotheses tested |
| `search` | The search stages in order. Each stage has `stage`, `hypotheses`, `seconds`, and `outcome`. |
| `seconds` | The wall-clock time of the solve, in seconds |
| `wcs_cards` | The solution as 80-character FITS header cards, including `IMAGEW` and `IMAGEH` |
| `pairs` | The matched sources. Each pair has the source's index in the source list (`source`) and position (`x`, `y`), the catalog star's UCAC4 running number (`star_id`), position at the solve epoch (`ra_deg`, `dec_deg`), and magnitude (`mag`), the star's predicted position (`predicted_x`, `predicted_y`), the source position minus the predicted position in pixels (`residual_px`), and the probability that the match is real (`inlier_prob`). |
| `covariance` | The covariance of the fitted coefficients, the same for each tangent-plane axis. `matrix_rad2` is the matrix row by row, in square radians. `terms` gives the powers of the pixel offsets from `crpix_fits`, divided by `norm_px`, that the coefficients multiply. |
| `detection` | The detection report without its source lists or its separate `identity` field. For more information about its fields, see [Read the detection report](#read-the-detection-report). |
| `identity` | The Fayde version and commit, the catalog and the solve epoch, and the pattern indexes with their content hashes |

In a `fayde verify` report, `accepted` is `true` when `log_odds` reaches `threshold`. The report
also has `matches`, `confident_matches`, `rms_px`, `references`,
`sources`, `pairs`, `priors`, and `identity`.

### Read the detection report

The `fayde detect` command prints a detection report. The `detection` field in a
`fayde solve --image` report omits the `sources` and `rejected` lists and the detection
report's separate `identity` field.

The `mode` field is `points`, `trails`, `kernel`, or `refused`. For more information about these values,
see [Detect stars without solving](#detect-stars-without-solving). The `census` field holds
the evidence for the decision. Its `reason` field has one of the following values:

| `reason` | `mode` | Meaning |
|---|---|---|
| `no elongated sources` | `points` | No bright source is elongated. |
| `no common trail direction` | `points` | Too few elongated sources share one direction. |
| `trails shorter than 7 FWHM` | `points` | The trails are shorter than seven times the stars' full width at half maximum (FWHM) across the trail. |
| `one straight trail` | `trails` | Every star is a straight trail of one length and direction. |
| `shared curved star path` | `kernel` | The template learned from the stars has a measurable bend. |
| `star paths vary with field rotation` | `kernel` | The star trails are tangential to a common center and grow longer with distance from it. |
| `point stars as well as trails (...)` | `refused` | The image has round stars as well as trails. The reason gives the number of each. |
| `two trail populations` | `refused` | The trails have two directions. To choose one, use `--trail-angle`. |
| `trails longer than the frame` | `refused` | Most trails touch the frame edge, so Fayde can't measure their length. |
| `too little light for the second estimate` | `refused` | The image has too little light for a second, independent estimate of the trail. |
| `the two estimates disagree` | `refused` | Two independent estimates of the trail disagree. |
| `the trails change across the frame` | `refused` | The trail direction or length differs between regions of the frame. |
| `too few trails to check that they are straight` | `refused` | Fayde fitted too few trails to check that they're straight. |

The `trail` field is `null` except in `trails` mode, when it contains the following fields:

| Field | Meaning |
|---|---|
| `length_px` | The full length of the trail, in pixels |
| `angle_deg` | The direction of the trail, in degrees from +x toward +y, from −90° to 90° |
| `width_sigma_px` | The Gaussian σ of the stars across the trail, in pixels |
| `sigma_length_px` | The 1σ error of `length_px`, in pixels |
| `sigma_angle_deg` | The 1σ error of `angle_deg`, in degrees |
| `fitted` | The number of trails that Fayde fitted to measure the trail. When the value is 0, Fayde uses the census values. |

The `kernel` field contains the following measurements in `kernel` mode. Its value is
`null` in the other modes.

| Field | Meaning |
|---|---|
| `width_px`, `height_px` | Bounding dimensions of the correlation taps. For field rotation, these are the maximum dimensions over the regions. |
| `stacked` | Number of stars that contributed to the template or common transverse profile. |
| `path_length_px` | Estimated path length. For field rotation, this is the median over the regions. |
| `bend_px` | RMS deviation of the path from its best straight line. For field rotation, this is the median over the regions. |
| `regions` | Number of regions with distinct templates. |

The `sources` field lists the stars that Fayde passes to the solver, and the `rejected` field
lists the detections that it set aside. Each detection has a 1-based position (`x`, `y`), a
1σ position error (`sigma_px`), and 1σ errors along the trail (`sigma_along_px`) and across it
(`sigma_across_px`). In `trails` mode, a star's `sigma_px` is the larger of `sigma_along_px`
and `sigma_across_px`. In `kernel` mode, the errors are along and across the template's
major axis, and `sigma_px` is the largest 1σ error over all directions. When the stars are
points, the three errors are equal. The `class` field describes each detection.

- In `sources`: `star`, `saturated` for a star with saturated pixels, or `extended` for a
  source that is wider than the stars
- In `rejected`: `sharp` for a hot pixel or cosmic-ray hit, `streak` for a detection on a traced
  streak, or `point` for a point-like source in an image whose stars are trails, such as the
  tracked object

With a profile that has a `charge_transfer` field, the positions in `sources` are corrected.
Each star's `charge_transfer_px` field holds the shift that Fayde removed, in pixels. The
report's `charge_transfer` field holds the term that Fayde used, the shift at the frame's
epoch, the number of stars moved (`moved`), and the largest shift (`max_px`). Without such a
profile, `charge_transfer` is `null`.

## Check or refit an existing WCS

To check a WCS from any solver against the catalog, use `fayde verify`:

```bash
fayde verify --image FRAME --ucac4 UCAC4_DIR --wcs WCS
```

To fit the WCS again from the matched stars and write the result, use `fayde refit`:

```bash
fayde refit --image FRAME --ucac4 UCAC4_DIR --wcs WCS -o OUTPUT
```

Replace the following:

- `WCS`: The path to the WCS file to check, for example `frame.wcs`
- `OUTPUT`: The path to the WCS file to write, which must differ from `WCS`

Both commands also accept a source list: use `--sources SOURCES --header FRAME` in place of
`--image FRAME`.

## Find catalog stars at a position

To find catalog stars at positions under a given WCS, use `fayde verify`. Specify pixel
positions with `--at`, sky positions with `--at-sky`, or a source list with `--at-sources`.

To check positions from a source list, use `--at-sources`:

```bash
fayde verify --ucac4 UCAC4_DIR --wcs WCS --header FRAME --at-sources SOURCES
```

To check a pixel position and a sky position at a specified epoch, use `--at` and `--at-sky`:

```bash
fayde verify --ucac4 UCAC4_DIR --wcs WCS --epoch 2020.89 --at 412.3,88.1 --at-sky 172.5,1.6
```

Pixel positions are 1-based, as in FITS: the center of the first pixel is at (1, 1). Sky
positions are the right ascension and declination in degrees.

For each position, the report lists the following:

- The nearest stars at the observation epoch and their separations
- The combined 1σ uncertainty of the position, catalog, and WCS (`--wcs-sigma`)
- The probability of an offset this large or larger if the position is the star (`p_consistent`)
- The probability that some star lies this close to an arbitrary position (`p_chance`)
