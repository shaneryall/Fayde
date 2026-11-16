"""Generate synthetic TAN WCS test vectors with astropy (zero-based pixel origin)."""
import json
import sys

import numpy as np
from astropy.wcs import WCS

CASES = [
    # name, crval (deg), crpix (FITS 1-based), cd (deg/px), size
    ("equator_north_up", (150.0, 2.0), (512.5, 512.5),
     [[-8.3e-4, 0.0], [0.0, 8.3e-4]], (1024, 1024)),
    ("rotated_30deg", (83.82, -5.39), (1024.5, 768.5),
     [[-8.3e-4 * np.cos(np.radians(30)), 8.3e-4 * np.sin(np.radians(30))],
      [8.3e-4 * np.sin(np.radians(30)), 8.3e-4 * np.cos(np.radians(30))]], (2048, 1536)),
    ("near_pole", (12.0, 88.5), (512.5, 512.5),
     [[-1.0e-3, 2.0e-5], [3.0e-5, 1.0e-3]], (1024, 1024)),
    ("ra_wrap", (359.95, -30.0), (512.5, 512.5),
     [[-8.3e-4, 0.0], [0.0, 8.3e-4]], (1024, 1024)),
    ("flipped_parity", (210.0, 45.0), (300.5, 200.5),
     [[8.3e-4, 1.0e-5], [-1.0e-5, 8.3e-4]], (600, 400)),
    ("wide_field_10deg", (270.0, -60.0), (1000.5, 1000.5),
     [[-5.0e-3, 0.0], [0.0, 5.0e-3]], (2000, 2000)),
]

rng = np.random.default_rng(20261115)
out = []
for name, crval, crpix, cd, (nx, ny) in CASES:
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = crval
    w.wcs.crpix = crpix
    w.wcs.cd = np.array(cd)
    w.wcs.radesys = "ICRS"
    w.wcs.set()
    px = np.vstack([
        [[0.0, 0.0], [nx - 1.0, ny - 1.0], [nx - 1.0, 0.0], [0.0, ny - 1.0],
         [crpix[0] - 1.0, crpix[1] - 1.0]],
        rng.uniform([-0.5, -0.5], [nx - 0.5, ny - 0.5], size=(15, 2)),
    ])
    world = w.wcs_pix2world(px, 0)
    header = w.to_header()
    header["NAXIS1"] = nx
    header["NAXIS2"] = ny
    out.append({
        "name": name,
        "crval_deg": list(crval),
        "crpix_fits": list(crpix),
        "cd": [list(map(float, r)) for r in cd],
        "naxis": [nx, ny],
        "pixels_zero_based": px.tolist(),
        "world_deg": world.tolist(),
        "pixel_scale_arcsec": float(np.sqrt(abs(np.linalg.det(np.array(cd)))) * 3600.0),
        "header_cards": [card.image for card in header.cards],
    })

json.dump({"generator": "astropy " + __import__("astropy").__version__,
           "pixel_origin": 0, "cases": out}, sys.stdout, indent=1)
