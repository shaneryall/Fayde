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
   run fails, Fayde returns an exit status of 1 and prints the error to standard error.

When a run writes no solution, Fayde removes any file that already exists at `OUTPUT`. Fayde
keeps the file if it rejects the command before the run starts, if it can't write to `OUTPUT`, or if
`OUTPUT` is one of the input files.

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

Fayde prints a JSON report with the background, the noise, the width of the stars, and every
source that it found. The report also lists the streaks that cross the image, such as satellite
trails, and the hot pixels and cosmic-ray hits that Fayde set aside. With `-o`, Fayde also
writes the sources to `SOURCES`, in the format that `--sources` reads. To
list the detection options and their default values, run `fayde detect --help`. The same
options work with `fayde solve --image`.

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
| Pointing | `CRVAL1` and `CRVAL2` of a header WCS, `RA` and `DEC`, or `OBJCTRA` and `OBJCTDEC` |
| Pixel scale | `XPIXSZ` with `FOCALLEN`, the `CD` matrix, `CDELT1` and `CDELT2`, `PIXSCALE`, or `SECPIX` |
| Observation time | `JD-OBS`, `MJD-OBS`, or `DATE-OBS`, with `EXPTIME` or `EXPOSURE` |
| Image size | `NAXIS1` and `NAXIS2` |

If the header has no observation time, Fayde uses the epoch 2000. To set the epoch, use
`--epoch`.

Fayde searches within 2° of the pointing. To change the limit, use `--radius`. To specify the
pointing accuracy, use `--pointing-sigma` with the 1-sigma error in degrees. When
that error is at most a sixth of the limit and the pixel scale is known to within 3%, Fayde
first searches within 3 sigma of the pointing. It widens the search to the limit only when it
finds no solution there.

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

The repository contains an example profile, `profiles/example.json`. A release archive has the
same file in the extracted archive's `profiles` directory. The example's values are
illustrative. To create a profile for your instrument, follow these steps:

1. Copy `profiles/example.json`.
2. Edit the copy with your instrument's values. The `source` field of the example explains the
   scale and roll values.

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

2. List only the frames whose report has `solved` set to `true` in a text file, with one frame
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
| `scale_arcsec`, `scale_sigma_arcsec` | The median pixel scale of the solutions and its 1-sigma uncertainty, in arcseconds per pixel. With `--binning`, this is the scale of an unbinned frame. |
| `parity` | The parity of the solutions. Fayde omits this field when the solutions disagree. |
| `roll` | With `--roll-keyword`, the rule that gives the roll from the roll keyword, with its 1-sigma error and hard limit in degrees |
| `binning` | With `--binning`, the binning keywords that you named |
| `pointing_sigma_deg` | The 1-sigma error of the pointing in the headers, in degrees. Fayde omits this field when fewer than five headers have a pointing. |
| `calibration` | The numbers behind these values: frames used, spreads, and largest deviations |

You can add two optional fields manually. The `pointing_radius_deg` field sets the search limit
around the pointing in degrees, in place of `--radius`. The `header_wcs` field, when set to
`false`, has the same effect as `--no-header-wcs` on every frame. To include `header_wcs` in the
profile during calibration, add `--no-header-wcs` to the `fayde calibrate` command.

### Solve without a pointing or scale

When you know the pointing and have a precise pixel scale, Fayde needs only the catalog. In
that case, omit `--index`. Fayde searches a pattern index only if it finds no solution around
the pointing. This extra search delays the report when a frame has no solution. When
the pointing or scale is uncertain or unknown, Fayde also needs a pattern index, which you
build once from the catalog. To build the recommended indexes, see the
[index guide](index-set.md).

To solve with a pattern index and ignore all priors, use `--index` and `--no-priors`:

```bash
fayde solve --image FRAME --ucac4 UCAC4_DIR --index INDEX --no-priors -o OUTPUT
```

Replace `INDEX` with the path to a pattern index, for example `allsky-1.0.fidx`. To use
several indexes, add `--index` once for each file.

### Limit the search time

Search time increases with the area within the pointing radius. The following measurements
show exhausted searches on a 32-core x86_64 Linux computer. The frames had a pixel scale of
3 arcseconds per pixel, a scale uncertainty of 0.1%, and unconstrained roll:

| Pointing radius | Hypotheses until exhausted | Time |
|---|---|---|
| 0.5° | 250–760 | 0.01–0.02 s |
| 1° | 830–2,500 | 0.02–0.05 s |
| 2° | 2,900–10,500 | 0.07–0.16 s |
| 4° | 11,000–40,000 | 0.24–0.62 s |
| 8° | 47,000–157,000 | 0.9–2.3 s |

Frames with hundreds of sources took more than 18 s at an 8° radius. To limit a search, use
these options:

- Use `--max-hypotheses` to limit the work independently of the computer's speed.
- Use `--deadline` to limit wall-clock time from the start of the request.

Each `search` stage in the report ends with an `outcome`: `solved`, `exhausted`, `deadline`, or
`hypothesis limit`. When Fayde finds no solution, the report names the stage that a limit
stopped. A search with a 1° radius takes about 50 times less work than one with an 8° radius.

## Check or refit an existing WCS

Replace `WCS` with the path to the WCS file to check, for example `frame.wcs`. `OUTPUT` must be
a different file from `WCS`.

To check a WCS from any solver against the catalog, use `fayde verify`:

```bash
fayde verify --image FRAME --ucac4 UCAC4_DIR --wcs WCS
```

To fit the WCS again from the matched stars and write the result, use `fayde refit`:

```bash
fayde refit --image FRAME --ucac4 UCAC4_DIR --wcs WCS -o OUTPUT
```

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

For each position, the report lists the following:

- The nearest stars at the observation epoch and their separations
- The combined 1σ uncertainty of the position, catalog, and WCS (`--wcs-sigma`)
- The probability of an offset this large or larger if the position is the star (`p_consistent`)
- The probability that some star lies this close to an arbitrary position (`p_chance`)
