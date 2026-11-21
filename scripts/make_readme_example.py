# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy", "astropy"]
# ///
r"""The example report in the README's quick start: `fayde solve` on a frame drawn from catalog stars.

    uv run scripts/make_readme_example.py --fayde FAYDE --catalog CATALOG -o OUT

FAYDE is the `fayde` command to run, for example `target/release/fayde`. CATALOG is the catalog
that the quick start downloads; give `--fov 0.7` for UCAC4. OUT is a directory for the frame,
the solution and the full report.

The script asks FAYDE for the catalog stars in the field (`fayde verify --at-sources` on a grid
of positions), at the epoch of the frame. It draws them into a 1024 × 1024 frame centred at
RA 150°, Dec +20°, rotated 12°, with Gaussian stars (σ 1.5 px), 200 counts of sky, Poisson noise
and 6 counts of read noise from random seed 7. The header gives a pointing 0.05° off in each
axis, the pixel scale to two decimals, and the time. The script then runs the quick start's
`fayde solve --image` on the frame, prints the report with only the fields that the README
shows, and prints how far the solution is from the WCS that drew the frame.
"""

import argparse
import json
import math
import subprocess
from pathlib import Path

import numpy as np
from astropy.io import fits

SIDE = 1024
CENTRE = (150.0, 20.0)  # RA, Dec of the frame centre, degrees
ROLL = 12.0  # position angle of the image +y axis, east of north, degrees
OFFSET = 0.05  # error of the header pointing in RA and Dec, degrees
DATE = "2024-03-15T03:00:00"
EPOCH = 2024.2  # Julian year of DATE, to two decimals
SIGMA = 1.5  # Gaussian star width, px
SKY, READ = 200.0, 6.0  # counts
SHOWN = {"crval_deg", "log_odds", "matches", "pixel_scale_arcsec", "rms_px", "roll_deg", "solved",
         "worst_sigma_arcsec"}


def cd_matrix(scale_deg):
    """The CD matrix of a frame that shows the sky as seen (east left when north is up)."""
    s, c = math.sin(math.radians(ROLL)), math.cos(math.radians(ROLL))
    return np.array([[-scale_deg * c, scale_deg * s], [scale_deg * s, scale_deg * c]])


def wcs_header(scale_deg):
    h = fits.Header()
    h["CTYPE1"], h["CTYPE2"] = "RA---TAN", "DEC--TAN"
    h["CRVAL1"], h["CRVAL2"] = CENTRE
    h["CRPIX1"] = h["CRPIX2"] = (SIDE + 1) / 2
    cd = cd_matrix(scale_deg)
    h["CD1_1"], h["CD1_2"], h["CD2_1"], h["CD2_2"] = cd[0, 0], cd[0, 1], cd[1, 0], cd[1, 1]
    h["IMAGEW"] = h["IMAGEH"] = SIDE
    return h


def to_pixels(ra, dec, scale_deg):
    """Gnomonic projection about CENTRE, then the inverse CD matrix: FITS 1-based pixels."""
    a0, d0 = np.radians(CENTRE[0]), np.radians(CENTRE[1])
    a, d = np.radians(ra), np.radians(dec)
    cosc = np.sin(d0) * np.sin(d) + np.cos(d0) * np.cos(d) * np.cos(a - a0)
    xi = np.degrees(np.cos(d) * np.sin(a - a0) / cosc)
    eta = np.degrees((np.cos(d0) * np.sin(d) - np.sin(d0) * np.cos(d) * np.cos(a - a0)) / cosc)
    dx, dy = np.linalg.solve(cd_matrix(scale_deg), np.vstack([xi, eta]))
    return dx + (SIDE + 1) / 2, dy + (SIDE + 1) / 2


def catalog_stars(fayde, catalog_args, out, scale_deg):
    """Every catalog star near the frame, from `fayde verify --at-sources` on a grid."""
    step, most = 64, 1000
    grid = np.arange(step / 2, SIDE + step / 2, step) + 0.5
    x, y = (g.ravel() for g in np.meshgrid(grid, grid))
    table = fits.BinTableHDU.from_columns([fits.Column(name="X", format="D", array=x),
                                           fits.Column(name="Y", format="D", array=y)])
    table.header["IMAGEW"] = table.header["IMAGEH"] = SIDE
    fits.HDUList([fits.PrimaryHDU(), table]).writeto(out / "grid.fits", overwrite=True)
    fits.PrimaryHDU(header=wcs_header(scale_deg)).writeto(out / "truth.wcs", overwrite=True)
    radius = 0.75 * step * scale_deg * 3600  # covers each grid cell, corners included
    run = subprocess.run(
        [fayde, "verify", *catalog_args, "--wcs", str(out / "truth.wcs"), "--epoch", str(EPOCH),
         "--at-sources", str(out / "grid.fits"), "--search-radius", f"{radius:.1f}",
         "--max-stars", str(most)],
        capture_output=True, text=True, check=True)
    report = json.loads(run.stdout)
    stars = {}
    for position in report["positions"]:
        if len(position["stars"]) >= most:  # verify returns only the nearest stars
            raise SystemExit("too many stars near one grid position: use a smaller step")
        for star in position["stars"]:
            stars[star["star_id"]] = (star["ra_deg"], star["dec_deg"], star["mag"])
    return np.array(sorted(stars.values()))


def draw(stars, scale_deg):
    rng = np.random.default_rng(7)
    x, y = to_pixels(stars[:, 0], stars[:, 1], scale_deg)
    mag = stars[:, 2]
    faint = np.nanpercentile(mag, 95)
    flux = 1600.0 * 10 ** (-0.4 * (mag - faint))  # the faintest stars at about 20σ
    image = np.full((SIDE, SIDE), SKY)
    yy, xx = np.mgrid[-10:11, -10:11]
    for xs, ys, f in zip(x - 1, y - 1, flux):  # 0-based pixel centres
        if not (np.isfinite(f) and -10 <= xs < SIDE + 10 and -10 <= ys < SIDE + 10):
            continue
        cx, cy = int(round(xs)), int(round(ys))
        stamp = np.exp(-((xx + cx - xs) ** 2 + (yy + cy - ys) ** 2) / (2 * SIGMA**2))
        stamp *= f / (2 * math.pi * SIGMA**2)
        y0, x0 = cy - 10, cx - 10
        sy = slice(max(y0, 0), min(y0 + 21, SIDE))
        sx = slice(max(x0, 0), min(x0 + 21, SIDE))
        image[sy, sx] += stamp[sy.start - y0:sy.stop - y0, sx.start - x0:sx.stop - x0]
    return (rng.poisson(image) + rng.normal(0, READ, image.shape)).astype(np.float32)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--fayde", required=True)
    which = ap.add_mutually_exclusive_group(required=True)
    which.add_argument("--catalog", help="the catalog, as `fayde solve --catalog` takes it")
    which.add_argument("--ucac4", help="a UCAC4 directory, for versions without --catalog")
    ap.add_argument("--fov", type=float, default=2.0, help="frame width in degrees (default 2)")
    ap.add_argument("-o", "--out", required=True, type=Path)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    catalog_args = ["--catalog", a.catalog] if a.catalog else ["--ucac4", a.ucac4]
    scale_deg = a.fov / SIDE

    stars = catalog_stars(a.fayde, catalog_args, a.out, scale_deg)
    h = fits.Header()
    h["RA"], h["DEC"] = CENTRE[0] + OFFSET, CENTRE[1] - OFFSET
    h["PIXSCALE"] = (round(scale_deg * 3600, 2), "arcsec per pixel")
    h["DATE-OBS"], h["EXPTIME"] = DATE, 30.0
    fits.PrimaryHDU(draw(stars, scale_deg), header=h).writeto(a.out / "frame.fits", overwrite=True)

    run = subprocess.run(
        [a.fayde, "solve", "--image", str(a.out / "frame.fits"), *catalog_args,
         "-o", str(a.out / "frame.wcs")],
        capture_output=True, text=True)
    (a.out / "report.json").write_text(run.stdout)
    report = json.loads(run.stdout)
    print(json.dumps({k: v for k, v in report.items() if k in SHOWN}, indent=2))  # in report order
    if report.get("solved"):
        ra, dec = report["crval_deg"]
        cos_sep = (math.sin(math.radians(dec)) * math.sin(math.radians(CENTRE[1]))
                   + math.cos(math.radians(dec)) * math.cos(math.radians(CENTRE[1]))
                   * math.cos(math.radians(ra - CENTRE[0])))
        error = math.degrees(math.acos(min(1.0, cos_sep))) * 3600
        x, y = to_pixels(stars[:, 0], stars[:, 1], scale_deg)
        inside = np.sum((x >= 0.5) & (x <= SIDE + 0.5) & (y >= 0.5) & (y <= SIDE + 0.5))
        print(f"{inside} catalog stars in the frame; centre {error:.3f} arcsec from the drawn WCS")


if __name__ == "__main__":
    main()
