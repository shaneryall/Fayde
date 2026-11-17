# Use Fayde from Python

To install the package, see [Install the Python package](install.md#install-the-python-package).

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
