"""Tests of the Python bindings.

Need the UCAC4 zones and oracle the Rust data tests use (``FAYDE_UCAC4_DIR``, default
``<repo>/data/ucac4``); skipped without them. The frame is synthesized from the oracle's stars
through a known TAN WCS, so the solution can be checked against the truth.
"""

import json
import math
import os
import pathlib
import threading

import pytest

import fayde

REPO = pathlib.Path(__file__).resolve().parents[3]
UCAC4 = pathlib.Path(os.environ.get("FAYDE_UCAC4_DIR", REPO / "data" / "ucac4"))
ORACLE = UCAC4 / "ucac4_oracle.json"

pytestmark = pytest.mark.skipif(
    not ORACLE.is_file(), reason=f"needs UCAC4 data and oracle at {UCAC4}"
)

CRVAL = (150.0, 2.0)
SCALE = 1.5  # arcsec/px
ROLL = 20.0  # deg
SIZE = 700


def tan_project(ra, dec):
    """Pixel (0-based) of (ra, dec) for the true WCS: CRPIX at the frame centre, parity
    negative (north up, east left before roll)."""
    a0, d0 = map(math.radians, CRVAL)
    a, d = math.radians(ra), math.radians(dec)
    cosc = math.sin(d0) * math.sin(d) + math.cos(d0) * math.cos(d) * math.cos(a - a0)
    xi = math.cos(d) * math.sin(a - a0) / cosc
    eta = (math.cos(d0) * math.sin(d) - math.sin(d0) * math.cos(d) * math.cos(a - a0)) / cosc
    xi, eta = math.degrees(xi) * 3600 / SCALE, math.degrees(eta) * 3600 / SCALE
    r = math.radians(ROLL)
    u = -xi * math.cos(r) + eta * math.sin(r)
    v = xi * math.sin(r) + eta * math.cos(r)
    c = (SIZE - 1) / 2
    return c + u, c + v


@pytest.fixture(scope="module")
def frame():
    cones = json.loads(ORACLE.read_text())["cones"]
    stars = next(c for c in cones if c["name"] == "z461_j2000")["stars"]
    stars = sorted(stars, key=lambda s: s["mag"] if s["mag"] is not None else 99)
    xs, ys, flux = [], [], []
    for s in stars:
        x, y = tan_project(s["ra_deg"], s["dec_deg"])
        if 5 < x < SIZE - 5 and 5 < y < SIZE - 5:
            xs.append(x)
            ys.append(y)
            flux.append(10 ** (-0.4 * (s["mag"] or 16)))
        if len(xs) == 60:
            break
    assert len(xs) >= 30
    return xs, ys, flux


@pytest.fixture(scope="module")
def solver():
    return fayde.Solver(UCAC4)


PRIORS = dict(
    width=SIZE,
    height=SIZE,
    # Local data holds zones z458-z463 only (dec +1.4 to +2.6): keep the search inside them.
    ra=CRVAL[0] + 0.03,
    dec=CRVAL[1] - 0.02,
    radius=0.1,
    scale=SCALE,
    epoch=2000.0,
    deadline=20.0,
)


def test_solves_and_reports_like_the_cli(solver, frame):
    x, y, flux = frame
    r = solver.solve(x, y, flux=flux, **PRIORS)
    assert r["solved"], r.get("reason")
    ra, dec = r["crval_deg"]
    assert abs(r["pixel_scale_arcsec"] - SCALE) / SCALE < 1e-3
    # The CLI report's fields are all there.
    for k in ("log_odds", "nfa", "search", "violated_priors", "solved_by_stage", "identity",
              "pairs", "wcs_cards", "covariance", "priors"):
        assert k in r, k
    assert r["identity"]["catalog"]["name"].startswith("UCAC4")
    assert r["matches"] >= 30
    p = r["pairs"][0]
    # Pairs are FITS 1-based: the input was 0-based.
    assert abs(p["x"] - 1 - x[p["source"]]) < 1e-9
    # CRPIX is the frame centre, so CRVAL is close to the true centre.
    assert math.hypot((ra - CRVAL[0]) * math.cos(math.radians(dec)), dec - CRVAL[1]) < 0.3


def test_verify_and_refit_take_cards_or_a_dict(solver, frame):
    x, y, flux = frame
    r = solver.solve(x, y, flux=flux, **PRIORS)
    cards = r["wcs_cards"]
    v = solver.verify(x, y, cards, width=SIZE, height=SIZE, epoch=2000.0)
    assert v["accepted"] and v["confident_matches"] >= 30
    as_dict = {}
    for card in cards:
        key, _, rest = card.partition("=")
        val = rest.split("/")[0].strip().strip("'").strip()
        try:
            as_dict[key.strip()] = float(val) if "." in val or "E" in val else int(val)
        except ValueError:
            as_dict[key.strip()] = val
    f = solver.refit(x, y, as_dict, width=SIZE, height=SIZE, epoch=2000.0)
    assert f["solved"]
    assert abs(f["pixel_scale_arcsec"] - r["pixel_scale_arcsec"]) < 1e-6


def test_numpy_arrays_and_origin(solver, frame):
    np = pytest.importorskip("numpy")
    x, y, flux = frame
    a = solver.solve(np.asarray(x), np.asarray(y), flux=np.asarray(flux), **PRIORS)
    b = solver.solve(np.asarray(x) + 1, np.asarray(y) + 1, flux=flux, origin=1, **PRIORS)
    assert a["solved"] and b["solved"]
    assert a["crval_deg"] == pytest.approx(b["crval_deg"], abs=1e-9)


def test_threads_share_one_solver(solver, frame):
    x, y, flux = frame
    out = []
    ts = [threading.Thread(target=lambda: out.append(solver.solve(x, y, flux=flux, **PRIORS)))
          for _ in range(4)]
    for t in ts:
        t.start()
    for t in ts:
        t.join()
    assert len(out) == 4 and all(r["solved"] for r in out)
    assert len({tuple(r["crval_deg"]) for r in out}) == 1


def test_abstains_on_noise(solver):
    import random

    rng = random.Random(1)
    x = [rng.uniform(0, SIZE) for _ in range(40)]
    y = [rng.uniform(0, SIZE) for _ in range(40)]
    r = solver.solve(x, y, **dict(PRIORS, deadline=2.0))
    assert not r["solved"]
    assert "reason" in r and "hypotheses_tested" in r and "search" in r


def test_bad_options_are_errors(solver, frame):
    x, y, _ = frame
    with pytest.raises(ValueError, match="raduis"):
        solver.solve(x, y, raduis=1.0, **PRIORS)
    with pytest.raises(ValueError, match="y has"):
        solver.solve(x, y[:-1], **PRIORS)
    with pytest.raises(ValueError, match="flux has"):
        solver.solve(x, y, flux=[1.0], **PRIORS)
    with pytest.raises(ValueError, match="ra and dec"):
        solver.solve(x, y, **dict(PRIORS, dec=None))
    with pytest.raises(ValueError, match="image size"):
        solver.solve(x, y, ra=150.0, dec=2.0)
    with pytest.raises(ValueError):
        solver.solve(x, y, **dict(PRIORS, deadline=float("nan")))


def test_solver_errors_raise_with_the_report(solver, frame):
    x, y, _ = frame
    with pytest.raises(fayde.FaydeError) as e:
        solver.solve(x, y, **dict(PRIORS, pixel_error=-1.0))
    assert e.value.report["error"] is True
    assert "identity" in e.value.report


def test_header_and_instrument_profile(solver, frame):
    x, y, flux = frame
    header = {"CTYPE1": "RA---TAN", "CTYPE2": "DEC--TAN", "CRVAL1": CRVAL[0],
              "CRVAL2": CRVAL[1], "IMAGEW": SIZE, "IMAGEH": SIZE,
              "EPOCH": 2000.0}
    profile = {"name": "test camera", "scale_arcsec": SCALE, "scale_sigma_arcsec": 0.003,
               "parity": "negative"}
    r = solver.solve(x, y, flux=flux, header=header, instrument=profile, epoch=2000.0,
                     radius=0.1, deadline=20.0)
    assert r["solved"], r.get("reason")
    assert r["priors"]["instrument"] == "test camera"
    assert r["priors"]["pointing"]["origin"].startswith("CRVAL")
    assert r["priors"]["scale"]["origin"] == "instrument profile"


def test_optional_columns_follow_the_source_list_rules(solver, frame):
    np = pytest.importorskip("numpy")
    x, y, flux = frame
    base = solver.solve(x, y, flux=flux, **PRIORS)
    n = len(x)
    # Out-of-range values are absent, as in a source list: same result as not giving them.
    for kw in (dict(sigma_px=[0.0] * n), dict(sigma_px=[math.inf] * n),
               dict(p_star=[-1.0] * n), dict(p_star=[math.nan] * n)):
        r = solver.solve(x, y, flux=flux, **kw, **PRIORS)
        assert r["log_odds"] == base["log_odds"], kw
    # Masked entries are missing values.
    masked = np.ma.masked_array(np.full(n, 0.5), mask=np.ones(n, bool))
    r = solver.solve(x, y, flux=flux, sigma_px=masked, **PRIORS)
    assert r["log_odds"] == base["log_odds"]


def test_bad_columns_are_errors(solver, frame):
    np = pytest.importorskip("numpy")
    x, y, _ = frame
    n = len(x)
    with pytest.raises(ValueError, match="truncated"):
        solver.solve(x, y, truncated=[math.nan] * n, **PRIORS)
    with pytest.raises(ValueError, match="one-dimensional"):
        solver.solve(x, y, truncated=np.zeros((n, 1), bool), **PRIORS)
    with pytest.raises(ValueError, match="one-dimensional"):
        solver.solve(np.asarray(x)[:, None], y, **PRIORS)
    r = solver.solve(x, y, truncated=np.zeros(n, int), **PRIORS)
    assert r["priors"]["truncated_sources"] == 0


def test_long_header_keywords(solver, frame):
    x, y, flux = frame
    profile = {"scale_arcsec": SCALE, "scale_sigma_arcsec": 0.003, "parity": "negative",
               "roll": {"keyword": "DETECTOR_ROLL", "offset_deg": 0.0, "sign": 1.0,
                        "sigma_deg": 0.5, "max_dev_deg": 3.0}}
    r = solver.solve(x, y, flux=flux, instrument=profile, header={"DETECTOR_ROLL": 20.0},
                     **PRIORS)
    assert r["priors"]["roll"] is not None


def test_solver_keeps_working_after_chdir(frame, tmp_path, monkeypatch):
    x, y, flux = frame
    monkeypatch.chdir(UCAC4.parent)
    s = fayde.Solver(UCAC4.name)
    monkeypatch.chdir(tmp_path)
    assert s.solve(x, y, flux=flux, **PRIORS)["solved"]


def test_stars_at(solver, frame):
    x, y, flux = frame
    r = solver.solve(x, y, flux=flux, **PRIORS)
    cards = r["wcs_cards"]
    # Every frame source was made from a catalogue star: all consistent.
    at = solver.stars_at(cards, x[:10], y[:10], epoch=2000.0, pixel_error=0.5,
                         wcs_sigma_arcsec=r["worst_sigma_arcsec"])
    assert len(at["positions"]) == 10
    assert at["consistent"] == 10, [p["stars"][:1] for p in at["positions"]]
    p0 = at["positions"][0]
    assert abs(p0["pixel"][0] - 1 - x[0]) < 1e-6  # FITS 1-based out, 0-based in
    assert p0["stars"][0]["p_consistent"] > 0.0027
    # Sky positions: 40" north of the first source. A star may lie nearby, but not within
    # the position's error.
    ra, dec = p0["ra_deg"], p0["dec_deg"] + 40 / 3600
    off = solver.stars_at(cards, ra=[ra], dec=[dec], epoch=2000.0, search_radius_arcsec=120)
    q = off["positions"][0]
    assert not q["consistent"]
    assert q["stars"] and q["stars"][0]["p_consistent"] < 0.0027
    with pytest.raises(ValueError, match="stars_at options"):
        solver.solve(x, y, max_stars=2, **PRIORS)


def test_solve_image(solver, frame):
    np = pytest.importorskip("numpy")
    x, y, flux = frame
    rng = np.random.default_rng(3)
    img = 500.0 + 5.0 * rng.standard_normal((SIZE, SIZE))
    yy, xx = np.mgrid[0:SIZE, 0:SIZE]
    sig = 2.5 / 2.3548
    for sx, sy, f in zip(x, y, flux):
        # Oracle magnitudes as relative fluxes, scaled to S/N of tens to thousands.
        amp = 4e5 * f / max(flux)
        img += amp / (2 * math.pi * sig**2) * np.exp(-((xx - sx) ** 2 + (yy - sy) ** 2) / (2 * sig**2))
    r = solver.solve_image(img, **PRIORS)
    assert r["solved"], r.get("reason")
    assert abs(r["pixel_scale_arcsec"] - SCALE) / SCALE < 1e-3
    assert r["detection"]["count"] >= 20
    assert "sources" not in r["detection"] and "rejected" not in r["detection"]
    assert r["detection"]["streaks"] == []
