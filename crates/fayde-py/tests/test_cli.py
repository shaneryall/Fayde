"""Tests of the ``fayde`` command that the package installs.

They run without the catalog: the command's version, its usage error and its exit statuses.
"""

import json
import os
import pathlib
import subprocess
import sys
import sysconfig

import pytest

import fayde


def script() -> str:
    """The ``fayde`` console script of the environment that runs the tests."""
    name = "fayde.exe" if os.name == "nt" else "fayde"
    path = pathlib.Path(sysconfig.get_path("scripts")) / name
    assert path.exists(), f"{path}: the package did not install the fayde command"
    return str(path)


def run(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run([script(), *args], capture_output=True, text=True, timeout=60)


def test_version_matches_the_package() -> None:
    r = run("--version")
    assert r.returncode == 0, r.stderr
    assert r.stdout.strip() == f"fayde {fayde.__version__} ({fayde.commit})"


def test_usage_error_exits_1() -> None:
    r = run("no-such-command")
    assert r.returncode == 1
    assert r.stderr.startswith("fayde: usage:")


def test_help_exits_0() -> None:
    r = run("solve", "--help")
    assert r.returncode == 0, r.stderr
    assert "--ucac4" in r.stdout


def test_module_entry_point() -> None:
    r = subprocess.run(
        [sys.executable, "-m", "fayde._cli", "--version"], capture_output=True, text=True, timeout=60
    )
    assert r.returncode == 0, r.stderr
    assert r.stdout.startswith("fayde ")


def test_image_retry_report_matches_python(tmp_path: pathlib.Path) -> None:
    np = pytest.importorskip("numpy")
    (tmp_path / "u4b").mkdir()
    (tmp_path / "u4b" / "z001").touch()
    solver = fayde.Solver(tmp_path)

    size = 320
    yy, xx = np.mgrid[:size, :size]
    rng = np.random.default_rng(7)
    image = 500.0 + 5.0 * rng.standard_normal((size, size))
    for y in range(35, size - 20, 48):
        for x in range(35, size - 20, 60):
            sx, sy = x + rng.uniform(-8, 8), y + rng.uniform(-8, 8)
            along = np.maximum(np.abs(xx - sx) - 12, 0)
            image += rng.uniform(220, 550) * np.exp(-((yy - sy)**2 + along**2) / (2 * 1.3**2))

    cards = ["SIMPLE  =                    T", "BITPIX  =                  -32",
             "NAXIS   =                    2", f"NAXIS1  = {size:>20}", f"NAXIS2  = {size:>20}",
             "END"]
    header = "".join(f"{card:80}" for card in cards).encode("ascii").ljust(2880)
    pixels = image.astype(">f4").tobytes()
    path = tmp_path / "trailed.fits"
    path.write_bytes(header + pixels + bytes(-len(pixels) % 2880))

    for points_only in (False, True):
        flags = ["--no-trails"] if points_only else []
        cli = run("solve", "--image", str(path), "--ucac4", str(tmp_path), "--no-priors", *flags)
        assert cli.returncode == 0, cli.stderr
        native = solver.solve_image(path, no_priors=True, no_trails=points_only)
        report = json.loads(cli.stdout)
        detection = report["detection"]
        assert report["solved"] is False
        assert "sources" not in detection and "rejected" not in detection
        if points_only:
            assert detection["mode"] == "points"
            assert "retried_as_points" not in report
        else:
            assert detection["mode"] == "trails"
            assert report["retried_as_points"] == {"solved": False, "trails_census": detection["census"]}
        del detection["seconds"], native["detection"]["seconds"]
        assert report == native

    expired = solver.solve_image(path, no_priors=True, deadline=0.0)
    assert expired["solved"] is False
    assert "retried_as_points" not in expired
