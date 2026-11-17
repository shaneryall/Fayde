"""Tests of the ``fayde`` command that the package installs.

They run without the catalog: the command's version, its usage error and its exit statuses.
"""

import os
import pathlib
import subprocess
import sys
import sysconfig

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
