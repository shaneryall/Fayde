"""Synthetic FITS fixtures for fayde-fits (written with Astropy, as other tools would).

Writes into the directory given as argv[1]:
- frame.sources.fits: BINTABLE with X, Y (float32, FITS 1-based), FLUX (float64), IMAGEW/IMAGEH;
- frame.fits:  small primary image whose header carries the usual pointing/scale/time keys;
- frame_pc.fits: same with PC+CDELT instead of CD, and numeric RA/DEC;
- frame_mixed.fits: CD and PC both present (must be rejected as ambiguous).
"""
import os
import sys

import numpy as np
from astropy.io import fits

out = sys.argv[1]
os.makedirs(out, exist_ok=True)

cols = [
    fits.Column(name="X", format="E", array=np.array([10.5, 200.25, 767.0, 1.0], dtype=np.float32)),
    fits.Column(name="Y", format="E", array=np.array([20.5, 300.75, 1.0, 768.0], dtype=np.float32)),
    fits.Column(name="FLUX", format="D", array=np.array([1000.0, 500.5, 250.0, 125.0])),
]
t = fits.BinTableHDU.from_columns(cols)
t.header["IMAGEW"] = 768
t.header["IMAGEH"] = 768
fits.HDUList([fits.PrimaryHDU(), t]).writeto(os.path.join(out, "frame.sources.fits"), overwrite=True)

img = np.zeros((4, 6), dtype=np.int16)
h = fits.Header()
h["CTYPE1"] = "RA---TAN"
h["CTYPE2"] = "DEC--TAN"
h["CRVAL1"] = 80.9212084123
h["CRVAL2"] = 0.121863492816
h["CRPIX1"] = 384.5
h["CRPIX2"] = 384.5
h["CD1_1"] = 0.000799901108658
h["CD1_2"] = 0.000225
h["CD2_1"] = -0.000225
h["CD2_2"] = 0.000799901108658
h["RA"] = "05 22 57.30"
h["DEC"] = "+00 11 24.3"
h["DATE-OBS"] = "2020-11-21T16:46:40.074"
h["EXPTIME"] = 3.001
fits.PrimaryHDU(img, header=h).writeto(os.path.join(out, "frame.fits"), overwrite=True)

h2 = fits.Header()
h2["CTYPE1"] = "RA---TAN"
h2["CTYPE2"] = "DEC--TAN"
h2["CRVAL1"] = 210.5
h2["CRVAL2"] = -45.25
h2["CRPIX1"] = 3.5
h2["CRPIX2"] = 2.5
h2["CDELT1"] = -0.0008
h2["CDELT2"] = 0.0008
h2["PC1_1"] = 1.0
h2["PC2_2"] = 1.0
h2["MJD-OBS"] = 59174.69907
fits.PrimaryHDU(img, header=h2).writeto(os.path.join(out, "frame_pc.fits"), overwrite=True)

h3 = h.copy()
h3["PC1_1"] = 1.0
fits.PrimaryHDU(img, header=h3).writeto(os.path.join(out, "frame_mixed.fits"), overwrite=True)
print("ok")
