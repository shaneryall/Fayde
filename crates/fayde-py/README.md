# Fayde

The `fayde` Python package provides the [Fayde](https://github.com/shaneryall/Fayde) astrometric
plate solver. You provide a FITS image of a star field and information about the frame, such
as its pointing, pixel scale, and observation time. Fayde finds the stars in the image,
matches them against the UCAC4 star catalog, and returns a verified world coordinate system
(WCS) solution with its measured accuracy.

## Prerequisites

To use the package, you need the following:

- CPython 3.10 or later. One wheel works with every later version, and you don't need Rust
  to install it.
- A local copy of the UCAC4 catalog. The catalog isn't included with this package.
  For download instructions, see
  [Get the UCAC4 catalog](https://github.com/shaneryall/Fayde/blob/main/docs/catalog.md).

To install the package from PyPI, follow these steps:

1. Activate a virtual environment.

2. Install the package with `pip`:

   ```bash
   python -m pip install fayde
   ```

   The package also installs the `fayde` command. It runs the same code as the release
   binaries.

For more information about installation, see
[Install the Python package](https://github.com/shaneryall/Fayde/blob/main/docs/install.md#install-the-python-package).

## Solve a frame

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
to your FITS image. If the `FAYDE_UCAC4` environment variable contains the path to the catalog
directory, you can omit the catalog argument: `fayde.Solver()`.

The `image` argument accepts the path to a FITS file, a 2-D NumPy array, or an Astropy HDU. For
an array, pass the frame's header with `header`. You can use Fayde without NumPy or Astropy.

If you detect the stars yourself, pass their positions to `solve`. In this example, `x`, `y`,
and `flux` are sequences of the source positions and fluxes from your detector, and `hdr` is
the image header. The `indexes` and `instrument` arguments are optional. To solve the source
list, check the WCS, fit the WCS again, and find catalog stars at the source positions, use
these methods:

```python
solver = fayde.Solver("/path/to/UCAC4", indexes=["allsky-1.0.fidx"])
report = solver.solve(x, y, flux=flux, header=hdr, instrument="profiles/example.json",
                      deadline=5.0)
if report["solved"]:
    cards = report["wcs_cards"]
    check = solver.verify(x, y, cards, header=hdr)
    refit = solver.refit(x, y, cards, header=hdr)
    stars = solver.stars_at(cards, x, y, header=hdr)
```

The `verify` method checks the WCS against the catalog without changing it. The `refit` method
fits the WCS again from the matched stars. It returns the same report fields as `solve`,
including `wcs_cards` when it finds a solution.

The `solve`, `verify`, and `refit` methods take `x` and `y` as 0-based positions, as NumPy and
Astropy index pixels. If your positions are 1-based, as in FITS, pass `origin=1`. These methods
also accept the following optional sequences, with one value for each source:

- `flux`: The integrated flux in any linear unit
- `sigma_px`: The 1σ position error in pixels
- `p_star`: The prior probability that the source is a real star, from 0 to 1
- `truncated`: Whether the frame edge or a mask cuts the source, as booleans or as 0 and 1

Fayde ignores `sigma_px` values that aren't positive and `p_star` values outside that range.
In `flux`, `sigma_px`, and `p_star`, NaN means that the source has no value.

To find the stars in an image without solving it, use `fayde.detect`:

```python
detection = fayde.detect("frame.fits")
sources = detection["sources"]
```

The `x` and `y` positions in the detection report are 1-based, as in FITS. To solve these
sources, pass the report to `solve_detection`:

```python
report = solver.solve_detection(detection, header=hdr)
```

The `solve_detection` method passes the image size and each source's position, flux,
`sigma_px`, `p_star`, and `truncated` value to the solver. If the report identifies star
trails, the method also passes the measured trail unless you give `trail_length` and
`trail_angle`. It doesn't read the image header, so pass the header with `header`. If you pass
the detection report's positions to `solve` without `origin=1`, Fayde interprets each coordinate
as one pixel larger than its actual value.

If an image shows star trails in more than one direction, give the expected direction to
`fayde.detect` with `trail_angle`, in degrees from +x toward +y. Fayde still measures the
trail from the image. With `solve_image`, give `trail_length` and `trail_angle` together. The
angle guides detection, and the two values replace the measured trail in the solve.

The Python interface follows the `fayde` command line:

- The `solve`, `solve_detection`, `solve_image`, `verify`, and `refit` methods take the
  command-line options under the same names with underscores, such as `pointing_sigma` and
  `relax_priors_on_failure`. In Python, run `help(fayde)` to list them. The `detect` function
  and the `stars_at` method have their own options. To list those options, run
  `help(fayde.detect)` and `help(fayde.Solver.stars_at)`.
- The `solve_image` method also takes the options of `detect`, such as `detect_threshold`,
  which must be at least 1.5.
- Input coordinates `x` and `y` of `solve`, `verify`, `refit`, and `stars_at` are 0-based unless
  you pass `origin=1`. Pixel fields in the reports, including the source positions from
  `detect`, are 1-based, as in FITS.
- The `header` argument accepts a dictionary, a list of header cards, or an
  `astropy.io.fits.Header`. The WCS argument of `verify`, `refit`, and `stars_at` accepts the
  same types.
- The `instrument` argument accepts a profile dictionary or the path to a profile file.
- The report is a dictionary with the same fields as the command-line JSON report. For more
  information about these fields, see
  [Read the report](https://github.com/shaneryall/Fayde/blob/main/docs/usage.md#read-the-report).
- An invalid option raises `ValueError`. A solver failure raises `fayde.FaydeError`, which
  stores the report in its `report` attribute.
- Fayde releases the global interpreter lock (GIL) while it solves, so threads can share one
  `Solver`.

The package includes type information: the `py.typed` marker and the `fayde.Report`,
`fayde.AtReport`, and `fayde.DetectReport` types.

For more information, see the following pages:

- Priors and the report:
  [Use the command line](https://github.com/shaneryall/Fayde/blob/main/docs/usage.md)
- Pattern indexes:
  [Pattern indexes for sub-degree fields](https://github.com/shaneryall/Fayde/blob/main/docs/index-set.md)
- Versioning:
  [Versioning](https://github.com/shaneryall/Fayde/blob/main/docs/versioning.md)
- Changes in each version: [Changelog](https://github.com/shaneryall/Fayde/blob/main/CHANGELOG.md)

## License

Fayde uses the GNU General Public License, version 3 or later. The wheel includes `LICENSE`,
`NOTICE`, and
`THIRD-PARTY-LICENSES.md`. The latter lists the licenses of the Rust crates compiled into the
wheel. If you use Fayde in a publication, cite it as described in
[`CITATION.cff`](https://github.com/shaneryall/Fayde/blob/main/CITATION.cff).
