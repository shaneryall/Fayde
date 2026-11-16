# Fayde

The `fayde` Python package provides the [Fayde](https://github.com/shaner306/Fayde) astrometric
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
  [Get the UCAC4 catalog](https://github.com/shaner306/Fayde#get-the-ucac4-catalog).

To install the package from PyPI, follow these steps:

1. Activate a virtual environment.

2. Install the package with `pip`:

   ```bash
   python -m pip install fayde
   ```

For more information about installation, see
[Install the Python package](https://github.com/shaner306/Fayde#install-the-python-package).

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

The Python interface follows the `fayde` command line:

- The `solve`, `solve_image`, and `verify` methods take the command-line options under the same
  names with underscores, such as `pointing_sigma` and `relax_priors_on_failure`. In Python, run
  `help(fayde)` to list them. The `detect` and `stars_at` functions have their own options. To
  list those options, run `help(fayde.detect)` and `help(fayde.Solver.stars_at)`.
- Input coordinates `x` and `y` are 0-based unless you pass `origin=1`. Pixel fields in the
  report are 1-based, as in FITS.
- The `header` argument accepts a dictionary, a list of header cards, or an
  `astropy.io.fits.Header`.
- The `instrument` argument accepts a profile dictionary or the path to a profile file.
- The report is a dictionary with the same fields as the command-line JSON report.
- An invalid option raises `ValueError`. A solver failure raises `fayde.FaydeError`, which
  stores the report in its `report` attribute.
- Fayde releases the global interpreter lock (GIL) while it solves, so threads can share one
  `Solver`.

The package includes type information: the `py.typed` marker and the `fayde.Report`,
`fayde.AtReport`, and `fayde.DetectReport` types.

For more information about priors, pattern indexes, and versioning, see the
[Fayde README](https://github.com/shaner306/Fayde#readme). For the changes in each version, see the
[changelog](https://github.com/shaner306/Fayde/blob/main/CHANGELOG.md).

## License

Fayde uses the Apache License 2.0. The wheel includes `LICENSE`, `NOTICE`, and
`THIRD-PARTY-LICENSES.md`. The latter lists the licenses of the Rust crates compiled into the
wheel. If you use Fayde in a publication, cite it as described in
[`CITATION.cff`](https://github.com/shaner306/Fayde/blob/main/CITATION.cff).
