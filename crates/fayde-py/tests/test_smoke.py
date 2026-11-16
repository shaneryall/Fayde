"""Checks that need no catalogue data: run everywhere, including CI on every OS.

The catalogue is a few empty UCAC4 zone files, so every solve sees no stars.
"""

import math

import pytest

import fayde


@pytest.fixture()
def solver(tmp_path):
    (tmp_path / "u4b").mkdir()
    for z in range(455, 468):  # dec about +0.8 to +3.4: every zone the cones below touch
        (tmp_path / "u4b" / f"z{z}").write_bytes(b"")
    return fayde.Solver(tmp_path)


ARGS = dict(width=100, height=100, ra=150.0, dec=2.0, radius=0.1, scale=2.0, scale_sigma=0.001)
X, Y = [10.0, 50.0, 80.0, 30.0], [20.0, 60.0, 15.0, 90.0]


def test_version_and_commit():
    assert fayde.__version__
    assert fayde.commit


def test_missing_catalogue_is_a_value_error(tmp_path):
    with pytest.raises(ValueError, match="catalogue not found"):
        fayde.Solver(tmp_path)


def test_no_stars_abstains_with_the_report(solver):
    r = solver.solve(X, Y, **ARGS)
    assert r["solved"] is False
    for k in ("reason", "best_log_odds", "hypotheses_tested", "search", "priors", "identity"):
        assert k in r, k
    assert r["identity"]["fayde"]["version"] == fayde.__version__


def test_bad_options_are_value_errors(solver):
    with pytest.raises(ValueError, match="raduis"):
        solver.solve(X, Y, raduis=1.0, **ARGS)
    with pytest.raises(ValueError, match="y has"):
        solver.solve(X, Y[:-1], **ARGS)
    with pytest.raises(ValueError, match="truncated"):
        solver.solve(X, Y, truncated=[math.nan] * 4, **ARGS)
    with pytest.raises(ValueError, match="needs wcs"):
        solver._run("verify", X, Y, dict(ARGS))


def test_solver_errors_raise_faydeerror(solver):
    with pytest.raises(fayde.FaydeError) as e:
        solver.solve(X, Y, pixel_error=-1.0, **ARGS)
    assert e.value.report["error"] is True


def test_huge_deadline_does_not_crash(solver):
    assert solver.solve(X, Y, deadline=1e19, **ARGS)["solved"] is False


def test_stars_at_with_no_stars(solver):
    cards = ["CTYPE1  = 'RA---TAN'", "CTYPE2  = 'DEC--TAN'", "CRVAL1  = 150.0", "CRVAL2  = 2.0",
             "CRPIX1  = 50.5", "CRPIX2  = 50.5", "CD1_1   = -0.000555", "CD1_2   = 0.0",
             "CD2_1   = 0.0", "CD2_2   = 0.000555"]
    r = solver.stars_at(cards, [10.0, 20.0], [30.0, 40.0], ra=[150.0], dec=[2.0])
    assert len(r["positions"]) == 3
    assert r["consistent"] == 0 and all(p["stars"] == [] for p in r["positions"])
    with pytest.raises(ValueError, match="x with y"):
        solver.stars_at(cards, [1.0])


def _star_image(stars, size=200, fwhm=3.0, noise=5.0, seed=1):
    np = pytest.importorskip("numpy")
    rng = np.random.default_rng(seed)
    img = 500.0 + noise * rng.standard_normal((size, size))
    yy, xx = np.mgrid[0:size, 0:size]
    sig = fwhm / 2.3548
    for x, y, flux in stars:
        img += flux / (2 * math.pi * sig**2) * np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / (2 * sig**2))
    return img


def test_detect_finds_stars_in_an_array():
    np = pytest.importorskip("numpy")
    stars = [(70.3, 50.6, 20000.0), (150.8, 120.2, 8000.0), (30.0, 170.5, 5000.0)]
    img = _star_image(stars)
    img[10, 10] += 5000.0  # a hot pixel
    r = fayde.detect(img)
    assert r["count"] == 3, r["count"]
    found = sorted((s["x"], s["y"]) for s in r["sources"])
    # Array row 0 is FITS row 1: a star at column 70.3, row 50.6 is FITS (71.3, 51.6).
    want = sorted((x + 1, y + 1) for x, y, _ in stars)
    for (fx, fy), (wx, wy) in zip(found, want):
        assert abs(fx - wx) < 0.05 and abs(fy - wy) < 0.05
    # The hot pixel is set aside as "sharp", at FITS (11, 11).
    assert r["classes"]["sharp"] >= 1 and r["rejected_narrow"] == r["classes"]["sharp"]
    assert any(
        s["class"] == "sharp" and abs(s["x"] - 11) < 1 and abs(s["y"] - 11) < 1 for s in r["rejected"]
    )
    assert all(s["class"] == "star" and s["p_star"] > 0.9 for s in r["sources"])
    assert r["streaks"] == []
    # The stars' width is measured from the three stars.
    assert r["fwhm_estimated"] is True and abs(r["fwhm_px"] - 3.0) < 0.3
    # A float32 array and a masked array work too.
    assert fayde.detect(img.astype(np.float32))["count"] == 3
    masked = np.ma.masked_array(img, mask=np.zeros_like(img, bool))
    masked.mask[0:5, :] = True
    assert fayde.detect(masked)["count"] == 3


def test_detect_errors():
    np = pytest.importorskip("numpy")
    with pytest.raises(ValueError, match="2-D"):
        fayde.detect(np.zeros((4, 4, 4)))
    with pytest.raises(ValueError, match="missing.fits"):
        fayde.detect("missing.fits")
    with pytest.raises(ValueError, match="fwhm"):
        fayde.detect(np.zeros((64, 64)), fwhm=-1.0)
    with pytest.raises(ValueError, match="does not take pixel_error"):
        fayde.detect(np.zeros((64, 64)), pixel_error=1.0)


def test_detection_options_need_an_image(solver):
    with pytest.raises(ValueError, match="apply to images"):
        solver.solve(X, Y, fwhm=3.0, **ARGS)
