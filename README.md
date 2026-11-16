# Fayde

Fayde is an extremely fast and accurate astrometric plate solver. You provide a FITS image of a star field
and information about the frame, such as its pointing, pixel scale, and observation time.
Fayde finds the stars in the image, matches them against the UCAC4 star catalog, and returns a
verified world coordinate system (WCS) solution with its measured accuracy.

If you have your own star detector, you can give Fayde a list of sources instead of an image.

You can use Fayde in two ways:

- From the command line, with the `fayde` command
- From Python, with the `fayde` package

Fayde is at version 0.x. The command-line options, the JSON report, and the index file format
change only as described in [Versioning and releases](#versioning-and-releases). For the changes
in each version, see the [changelog](CHANGELOG.md).

## Prerequisites

To install Fayde, use a computer that runs one of these systems:

- Linux on x86_64 or 64-bit Arm, with glibc 2.17 or later
- macOS 11 or later, on Apple silicon or an Intel processor
- Windows 10 or later, on x86_64

You also need the following:

- About 9 GB of free disk space for the UCAC4 catalog
- For the Python package, CPython 3.10 or later

To solve a frame, you also need a FITS image from your own observation. The image must meet
these requirements:

- The stars are point sources, as in an image taken with sidereal tracking. Fayde doesn't yet
  detect trailed stars.
- The image isn't compressed with `fpack`. To decompress an image, use the `funpack` program
  from CFITSIO.
- The header gives the pointing and the pixel scale. If it doesn't, you can
  [specify them yourself](#specify-the-pointing-and-scale).

Install the command-line tool, the Python package, or both. Then
[get the UCAC4 catalog](#get-the-ucac4-catalog).

## Install the command-line tool

Each release includes a prebuilt `fayde` binary. You don't need Rust to install it. After you
install the tool, continue to [Get the UCAC4 catalog](#get-the-ucac4-catalog).

### Install on Linux or macOS

To install the latest release, follow these steps:

1. In a terminal, identify your operating system and processor architecture:

   ```bash
   uname -sm
   ```

2. Set the `platform` variable to the value that matches the output:

   | Output of `uname -sm` | Value |
   |---|---|
   | `Linux x86_64` | `linux-x86_64` |
   | `Linux aarch64` | `linux-aarch64` |
   | `Darwin arm64` | `macos-arm64` |
   | `Darwin x86_64` | `macos-x86_64` |

   For example, for the output `Darwin arm64`, set the variable to `macos-arm64`:

   ```bash
   platform=macos-arm64
   ```

3. Download the archive:

   ```bash
   curl -LO "https://github.com/shaner306/Fayde/releases/latest/download/fayde-$platform.tar.gz"
   ```

   Optional: Before you extract the archive, check it as described in
   [Verify a download](#verify-a-download).

4. Extract the archive:

   ```bash
   tar -xzf "fayde-$platform.tar.gz"
   ```

   Extraction creates a directory named `fayde-` followed by your platform value. The
   directory contains the `fayde` binary, the license files, and an example instrument profile.

5. Install the binary in `/usr/local/bin`, which is on your `PATH`:

   ```bash
   sudo mkdir -p /usr/local/bin && sudo install -m 755 "fayde-$platform/fayde" /usr/local/bin/fayde
   ```

   If `sudo` prompts for your password, enter it.

6. If macOS blocks a binary that you downloaded with a web browser instead of `curl`, remove
   the quarantine attribute to allow it to run:

   ```bash
   sudo xattr -d com.apple.quarantine /usr/local/bin/fayde
   ```

7. Check the installation:

   ```bash
   fayde --version
   ```

   The output is similar to the following:

   ```text
   fayde 0.2.2 (0123456789ab)
   ```

### Install on Windows

To install the latest release, follow these steps:

1. In PowerShell, download the archive:

   ```powershell
   curl.exe -LO https://github.com/shaner306/Fayde/releases/latest/download/fayde-windows-x86_64.zip
   ```

   Optional: Before you extract the archive, check it as described in
   [Verify a download](#verify-a-download).

2. Extract the archive:

   ```powershell
   Expand-Archive fayde-windows-x86_64.zip -DestinationPath .
   ```

   Extraction creates a directory named `fayde-windows-x86_64`. The directory contains
   `fayde.exe`, the license files, and an example instrument profile.

3. In the current PowerShell session, add the extracted directory to `Path`:

   ```powershell
   $env:Path += ";$PWD\fayde-windows-x86_64"
   ```

   The `Path` change lasts until you close PowerShell.

4. Check the installation:

   ```powershell
   fayde --version
   ```

   The output is similar to the following:

   ```text
   fayde 0.2.2 (0123456789ab)
   ```

Optional: To keep the `Path` change after you close PowerShell, follow these steps:

1. In the Windows search box, enter `environment variables`.
2. Select **Edit environment variables for your account**.
3. Select **Path**.
4. Click **Edit**.
5. Click **New**.
6. Enter the full path of the extracted directory.
7. Click **OK**.

### Verify a download

Each release includes a `SHA256SUMS` file that lists the SHA-256 checksum of every other file in
the release.

To check an archive on Linux or macOS, follow these steps:

1. In the directory that contains the archive, download the checksum file:

   ```bash
   curl -LO https://github.com/shaner306/Fayde/releases/latest/download/SHA256SUMS
   ```

2. Compare the checksums. On Linux, use `sha256sum`:

   ```bash
   sha256sum --check --ignore-missing SHA256SUMS
   ```

   On macOS, use `shasum`:

   ```bash
   shasum -a 256 --check --ignore-missing SHA256SUMS
   ```

   The output names the archive, followed by `OK`.

To check the archive on Windows, follow these steps:

1. In PowerShell, in the directory that contains the archive, download the checksum file:

   ```powershell
   curl.exe -LO https://github.com/shaner306/Fayde/releases/latest/download/SHA256SUMS
   ```

2. Print the checksum of the archive:

   ```powershell
   (Get-FileHash fayde-windows-x86_64.zip -Algorithm SHA256).Hash
   ```

3. Print the expected checksum:

   ```powershell
   Select-String fayde-windows-x86_64.zip SHA256SUMS
   ```

4. Compare the checksums from the two outputs. The archive is intact when they match.
   Uppercase and lowercase letters are equivalent.

### Install with Cargo

You need [Rust](https://www.rust-lang.org/tools/install) 1.85 or later to build the `fayde`
command from crates.io. To build and install the command, follow these steps:

1. In a terminal, install the `fayde-cli` crate:

   ```bash
   cargo install --locked fayde-cli
   ```

   Cargo installs `fayde` in its binary directory, which is `~/.cargo/bin` unless you changed it.
   If your shell can't find `fayde`, add that directory to your `PATH`.

2. Check the installation:

   ```bash
   fayde --version
   ```

   A build from crates.io shows `unknown` in place of the commit. The output is similar to the
   following:

   ```text
   fayde 0.2.2 (unknown)
   ```

## Install the Python package

The `fayde` package on PyPI includes a prebuilt wheel for every supported platform. Each wheel
supports CPython 3.10 and every later version, and you don't need Rust to install it. The
package contains the Python module only, not the `fayde` command.

If PyPI has no wheel for your platform, `pip` builds the package from its source distribution.
You need Rust 1.85 or later for that build.

To install the package, follow these steps:

1. On Linux or macOS, in a terminal, create a virtual environment:

   ```bash
   python3 -m venv .venv
   ```

   On Debian and Ubuntu, this command needs the `python3-venv` package.

   On Windows, in PowerShell, create the virtual environment:

   ```powershell
   py -m venv .venv
   ```

2. On Linux or macOS, activate the virtual environment:

   ```bash
   source .venv/bin/activate
   ```

   On Windows, activate the virtual environment:

   ```powershell
   .venv\Scripts\Activate.ps1
   ```

   If PowerShell reports that running scripts is disabled, follow these steps:

   1. Allow scripts for the current session:

      ```powershell
      Set-ExecutionPolicy -Scope Process -ExecutionPolicy RemoteSigned
      ```

   2. Activate the environment again.

3. Install the package:

   ```bash
   python -m pip install fayde
   ```

4. Check the installation:

   ```bash
   python -c "import fayde; print(fayde.__version__)"
   ```

   The output is the version number, for example `0.2.2`.


To solve a frame from Python, see [Use Fayde from Python](#use-fayde-from-python).

## Get the UCAC4 catalog

Fayde reads the UCAC4 catalog (Zacharias et al. 2013) from a local copy of its original binary
files. The catalog isn't included with Fayde. The Strasbourg astronomical Data Center (CDS)
publishes it as [catalog I/322A](https://vizier.cfa.harvard.edu/ftp/cats/I/322A/). The download is about 8.9 GB.

To download the catalog, use a terminal on Linux or macOS, or PowerShell on Windows. In the
download commands, use `curl.exe` in place of `curl` on Windows:

1. On Linux or macOS, create the catalog directories:

   ```bash
   mkdir -p UCAC4/u4b UCAC4/u4i
   ```

   On Windows, create the directories:

   ```powershell
   mkdir UCAC4\u4b, UCAC4\u4i
   ```

2. Go to the `u4b` directory:

   ```bash
   cd UCAC4/u4b
   ```

3. Download the 900 zone files, `z001` to `z900`:

   ```bash
   curl -fO "https://vizier.cfa.harvard.edu/ftp/cats/I/322A/UCAC4/u4b/z[001-900]"
   ```

   If the download stops, start it again from the file that was downloading. For example, if it
   stopped at `z412`, change the range to `z[412-900]`.

4. Go to the `u4i` directory:

   ```bash
   cd ../u4i
   ```

5. Download the two index files:

   ```bash
   curl -fO "https://vizier.cfa.harvard.edu/ftp/cats/I/322A/UCAC4/u4i/{u4index.unf,u4hpm.dat}"
   ```

   Fayde uses `u4index.unf` to seek directly to the stars that it needs and `u4hpm.dat` for the
   proper motions of fast-moving stars.

6. Return to the directory where you started:

   ```bash
   cd ../..
   ```

7. On Linux or macOS, count the zone files:

   ```bash
   ls UCAC4/u4b | wc -l
   ```

   On Windows, count the zone files:

   ```powershell
   (Get-ChildItem UCAC4\u4b).Count
   ```

   Check that the output is `900`. If it's lower, download the missing files again.

The catalog directory has this layout:

```text
UCAC4/
├── u4b/
│   ├── z001
│   ├── ...
│   └── z900
└── u4i/
    ├── u4hpm.dat
    └── u4index.unf
```

When a Fayde command has a `--ucac4` option, give it the path to the `UCAC4` directory, which
is the directory that contains `u4b`. Download all 900 zone files, because Fayde reports an
error when a zone that it needs is missing.

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
see the [source list format](docs/source-list.md).

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

Fayde searches within 2° of the pointing. To change the limit, use `--radius`. To list every
option and its default value, run `fayde solve --help`.

### Use an instrument profile

An instrument profile is a JSON file with settings that frames from one instrument share:
the pixel scale, parity, and a rule for deriving the roll from a header keyword. To use a
profile, add `--instrument`:

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

### Solve without a pointing or scale

When you know the pointing and have a precise pixel scale, Fayde needs only the catalog. When
the pointing or scale is uncertain or unknown, Fayde also needs a pattern index, which you
build once from the catalog. To build the recommended indexes, see the
[index guide](docs/index-set.md).

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

## Use Fayde from Python

To open the catalog and load the indexes only once, reuse one solver for every frame. To
find the stars in an image and solve it, create the solver and call `solve_image`:

```python
import fayde

solver = fayde.Solver("/path/to/UCAC4")
report = solver.solve_image("frame.fits")
if report["solved"]:
    cards = report["wcs_cards"]
```

Replace `/path/to/UCAC4` with the path to your catalog directory and `frame.fits` with the path
to your FITS image.

The `image` argument accepts the path to a FITS file, a 2-D NumPy array, or an Astropy HDU. For
an array, pass the frame's header with `header`. You can use Fayde without NumPy or Astropy.

If you detect the stars yourself, pass their positions to `solve`. In this example, `x`, `y`,
and `flux` are sequences of the source positions and fluxes from your detector, and `hdr` is
the image header. The `indexes` and `instrument` arguments are optional. To solve the source
list, check the WCS, and find catalog stars at the source positions, use these methods:

```python
solver = fayde.Solver("/path/to/UCAC4", indexes=["allsky-1.0.fidx"])
report = solver.solve(x, y, flux=flux, header=hdr, instrument="profiles/example.json",
                      deadline=5.0)
if report["solved"]:
    cards = report["wcs_cards"]
    check = solver.verify(x, y, cards, header=hdr)
    stars = solver.stars_at(cards, x, y, header=hdr)
```

To find the stars in an image without solving it, use `fayde.detect`:

```python
detection = fayde.detect("frame.fits")
sources = detection["sources"]
```

The Python interface follows the command line:

- The `solve`, `solve_image`, and `verify` methods take the command-line options under the same
  names with underscores, such as `pointing_sigma` and `relax_priors_on_failure`. In Python, run
  `help(fayde)` to list them. The `detect` and `stars_at` functions have their own options. To
  list those options, run `help(fayde.detect)` and `help(fayde.Solver.stars_at)`.
- Input coordinates `x` and `y` are 0-based unless you pass `origin=1`. Pixel fields in the
  report are 1-based, as in FITS and in the command-line report.
- The `header` argument accepts a dictionary, a list of header cards, or an
  `astropy.io.fits.Header`.
- The `instrument` argument accepts a profile dictionary or the path to a profile file.
- The report is a dictionary with the same fields as the command-line JSON report.
- An invalid option raises `ValueError`. A solver failure raises `fayde.FaydeError`, which
  stores the report in its `report` attribute.
- Fayde releases the global interpreter lock (GIL) while it solves, so threads can share one
  `Solver`.

## Develop Fayde

Fayde is a Cargo workspace. All crates share one release version. The workspace contains
the following crates:

| Crate | Contents |
|---|---|
| `fayde-core` | Sky positions, unit vectors, pixel conventions |
| `fayde-wcs` | TAN and SIP projection, and FITS WCS header cards |
| `fayde-catalog` | The `Catalog` trait, an in-memory catalog, the UCAC4 reader, and proper-motion propagation |
| `fayde-index` | HEALPix, k-d tree, uniform star selection, the 4-star pattern hash table, and the index file |
| `fayde-solve` | Hypothesis generation, verification, and the WCS fit |
| `fayde-fits` | FITS reader, images, source lists, priors from a header, and the `.wcs` writer |
| `fayde-detect` | Star and streak detection for FITS images |
| `fayde-api` | Request building and the JSON report, shared by the command line and Python |
| `fayde-cli` | The `fayde` command |
| `fayde-py` | The Python package |

To build Fayde from a clone of the repository, you need Git and Rust 1.85 or later. Follow
these steps:

1. In a terminal, clone the repository:

   ```bash
   git clone https://github.com/shaner306/Fayde
   ```

2. Go to the repository directory:

   ```bash
   cd Fayde
   ```

3. Build and install the `fayde` command:

   ```bash
   cargo install --locked --path crates/fayde-cli
   ```

4. Optional: To build and install the Python package, activate a virtual environment.

5. Optional: In the active virtual environment, build and install the Python package with `pip`:

   ```bash
   python -m pip install ./crates/fayde-py
   ```

Continuous integration runs code, documentation, packaging, and platform checks. To check a
change locally, run these commands:

```bash
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test
python3 scripts/check_release.py
python3 scripts/gen_third_party_licenses.py --check
```

After a change to `Cargo.lock`, run `python3 scripts/gen_third_party_licenses.py` without
`--check` to update the third-party license file.

The test vectors in `tests/fixtures/` are synthetic. The scripts in `scripts/` generate them
with Astropy.

Some tests compare Fayde with the real catalog. These tests use `#[ignore]` because the catalog
isn't in the repository. To include them, follow these steps:

1. Put a copy of UCAC4 at `data/ucac4/`, or set the `FAYDE_UCAC4_DIR` environment variable to
   the path to a copy.

2. Put `ucac4_oracle.json` in the same directory. This file holds cone-search results from an
   independent UCAC4 reader. The top of `crates/fayde-catalog/tests/reader_oracle.rs` describes
   its format.

3. Run every test:

   ```bash
   cargo test -- --include-ignored
   ```

The Python tests in `crates/fayde-py/tests` run with `pytest` against an installed wheel. They
skip the tests that need the catalog when it's missing.

For more information about the algorithms, see the [references](docs/references.md).
For more information about the architecture, see the [design document](docs/design/plate-solver.md).

## Versioning and releases

Fayde follows semantic versioning. While the version is 0.x, a breaking change increases the
minor number. The following changes are breaking:

- Removing a command-line or Python option, or changing its meaning
- Removing or renaming a JSON report field
- Changing the index file format or the bytes of an index, including its `content_hash`

Adding a report field isn't a breaking change. If a release changes indexes, you must rebuild
your indexes as directed in the release notes.

A fix that leaves reports and indexes unchanged is a patch release.

A tag `vX.Y.Z` that matches the workspace version starts the release pipeline. The pipeline
builds the `fayde` command and the Python wheels for Linux, macOS, and Windows, runs each one on
its own platform, and then publishes them as a GitHub release with SHA-256 checksums. Publish
the same version to PyPI and crates.io. For release instructions, see
[Release Fayde](docs/releasing.md).

## Cite Fayde

If you use Fayde in work that leads to a publication, cite it as described in
[`CITATION.cff`](CITATION.cff): Ryall 2026, Fayde, with the version that you used. The comment
cards of a solved `.wcs` file record the version and the commit.

## License

Fayde uses the Apache License 2.0. For the license terms, see [`LICENSE`](LICENSE).
For the copyright notice (2026 Shane Ryall), see [`NOTICE`](NOTICE). If you redistribute Fayde,
include both files.

The binaries and wheels contain compiled code from other Rust crates. The included
[`THIRD-PARTY-LICENSES.md`](THIRD-PARTY-LICENSES.md) file lists the licenses of those crates.
