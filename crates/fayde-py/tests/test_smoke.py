"""Checks that need no catalogue data: run everywhere, including CI on every OS.

The catalogue is a few empty UCAC4 zone files, so every solve sees no stars.
"""

import math

import pytest

import fayde


@pytest.fixture()
def solver(tmp_path):
    (tmp_path / "u4b").mkdir()
    for z in range(455, 468):
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
    img[10, 10] += 5000.0
    r = fayde.detect(img)
    assert r["count"] == 3, r["count"]
    found = sorted((s["x"], s["y"]) for s in r["sources"])
    want = sorted((x + 1, y + 1) for x, y, _ in stars)
    for (fx, fy), (wx, wy) in zip(found, want):
        assert abs(fx - wx) < 0.05 and abs(fy - wy) < 0.05
    assert r["classes"]["sharp"] >= 1 and r["rejected_narrow"] == r["classes"]["sharp"]
    assert any(
        s["class"] == "sharp" and abs(s["x"] - 11) < 1 and abs(s["y"] - 11) < 1 for s in r["rejected"]
    )
    assert all(s["class"] == "star" and s["p_star"] > 0.9 for s in r["sources"])
    assert r["streaks"] == []
    assert r["fwhm_estimated"] is True and abs(r["fwhm_px"] - 3.0) < 0.3
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


@pytest.mark.parametrize("limit", [{"deadline": 0}, {"max_hypotheses": 0}])
def test_exhausted_pair_search_does_not_read_missing_catalogue_zones(solver, limit):
    report = solver.solve(X, Y, **dict(ARGS, dec=50), **limit)
    assert not report["solved"] and report["hypotheses_tested"] == 0
    assert report["search"][-1]["outcome"] in ("deadline", "hypothesis limit")


def test_bright_crossing_star_is_masked_in_streak_photometry():
    np = pytest.importorskip("numpy")
    y, x = np.mgrid[:256, :256]
    sigma = 1.3
    a = np.deg2rad(28)
    u = (x - 128) * np.cos(a) + (y - 128) * np.sin(a)
    v = -(x - 128) * np.sin(a) + (y - 128) * np.cos(a)
    erf = np.vectorize(math.erf)
    ends = (erf((u + 90) / (sigma * 2**0.5)) - erf((u - 90) / (sigma * 2**0.5))) / 2
    signal = 30 * np.exp(-0.5 * (v / sigma)**2) * ends
    sky = 100 + np.random.default_rng(0).normal(size=x.shape) + signal
    def trace(image):
        (streak,) = fayde.detect(image, fwhm=2.35482 * sigma)["streaks"]
        return streak
    clean = trace(sky)
    blended = trace(sky + 3000 * np.exp(-((x - 128)**2 + (y - 128)**2) / (2 * sigma**2)))
    assert not clean["photometry_masked"] and blended["photometry_masked"]
    assert blended["flux"] == pytest.approx(signal.sum(), rel=0.04)
    assert 0.85 < blended["snr"] / clean["snr"] <= 1
    for key in ("x0", "y0", "x1", "y1"):
        assert blended[key] == pytest.approx(clean[key], abs=1)


def test_cutout_prior_uses_delivered_size_with_header_only_wcs(solver):
    header = dict(CTYPE1="RA---TAN", CTYPE2="DEC--TAN", CRVAL1=150.4, CRVAL2=2.0,
                  CRPIX1=50.5 - 960, CRPIX2=50.5, CDELT1=-1.5/3600, CDELT2=1.5/3600)
    options = dict(header=header, radius=0.1, max_hypotheses=0)
    image = _star_image([(x, y, 12000) for x, y in zip(X, Y)], size=100)
    for r in (solver.solve(X, Y, width=100, height=100, **options), solver.solve_image(image, **options)):
        pointing = r["priors"]["pointing"]
        assert pointing["origin"] == "WCS image centre"
        assert pointing["ra_deg"] == pytest.approx(149.99976, abs=1e-5)
        assert pointing["dec_deg"] == pytest.approx(1.99995, abs=1e-5)
