"""Fayde astrometric plate solver.

Open a :class:`Solver` once (catalogue and pattern indexes), then call :meth:`Solver.solve`,
:meth:`Solver.verify` or :meth:`Solver.refit` per frame. Each returns the report ``fayde solve``,
``fayde verify`` or ``fayde refit`` prints, as a dict.

Conventions:

- Input ``x``, ``y`` are 0-based (NumPy / Astropy) unless ``origin=1`` (FITS).
- Report pixel fields (``pairs[].x``, ``predicted_x``, ``covariance.crpix_fits``, ``wcs_cards``)
  are FITS 1-based, as in the command-line report.
- Angles are degrees, the scale arcsec/px, the deadline seconds.

Options (keyword arguments), the same as the command line's:

``width``, ``height``, ``ra``, ``dec``, ``radius``, ``pointing_sigma``, ``scale``,
``scale_sigma`` (fraction), ``roll``, ``roll_sigma``, ``roll_max``, ``parity``
(``"negative"``/``"positive"``), ``trail_length``, ``trail_angle``, ``epoch`` (Julian year),
``deadline``, ``max_hypotheses``, ``max_sip_order``, ``pixel_error``, ``accept_log_odds``, ``no_priors``,
``no_pointing``, ``search_outside_pointing``, ``relax_priors_on_failure``, ``header``,
``instrument``, ``origin``.

``header`` is a dict, a list of 80-character cards, or an ``astropy.io.fits.Header``.
``instrument`` is a profile dict or the path of a profile JSON file.
"""

from __future__ import annotations

import json
import math
import os
from typing import Any, Dict, Iterable, List, Optional, Sequence, cast

from . import _native
from ._types import AtReport, DetectReport, Report

__all__ = [
    "Solver",
    "FaydeError",
    "detect",
    "Report",
    "AtReport",
    "DetectReport",
    "__version__",
    "commit",
]

__version__: str = _native.__version__
#: Source commit of the compiled extension.
commit: str = _native.commit


class FaydeError(RuntimeError):
    """The solver failed (invalid request, catalogue error). ``report`` holds the JSON report,
    with ``"error": True`` and ``reason``, as the command line prints it."""

    def __init__(self, report: Dict[str, Any]):
        super().__init__(report.get("reason", "fayde error"))
        self.report = report


def _column(name: str, values: Any, missing: Any) -> List[Any]:
    """A 1-D column as a list. Masked entries (numpy.ma) become ``missing``."""
    if hasattr(values, "filled"):  # numpy masked array
        values = values.filled(missing) if hasattr(values, "mask") else values
    if getattr(values, "ndim", 1) != 1:
        raise ValueError(f"{name} must be one-dimensional, not shape {values.shape}")
    if hasattr(values, "tolist"):  # NumPy arrays, astropy columns
        values = values.tolist()
    out: List[Any] = list(values)
    if any(isinstance(v, (list, tuple)) for v in out):
        raise ValueError(f"{name} must be one-dimensional")
    return out


def _floats(name: str, values: Any) -> Optional[List[float]]:
    """Floats; for the optional columns, NaN (and a masked entry) means "no value"."""
    if values is None:
        return None
    if hasattr(values, "astype") and hasattr(values, "mask"):
        values = values.astype(float)
    out: List[float] = []
    for v in _column(name, values, math.nan):
        if v is None:
            v = math.nan
        out.append(float(v))
    return out


def _required(name: str, values: Any) -> List[float]:
    """A required float column; ``None`` means an empty one."""
    return _floats(name, values) or []


def _bools(name: str, values: Any) -> Optional[List[bool]]:
    """Flags: booleans, or the integers 0 and 1. Anything else (NaN, 2, "yes") is an error, not
    silently true."""
    if values is None:
        return None
    out: List[bool] = []
    for v in _column(name, values, False):
        if isinstance(v, bool):
            out.append(v)
        elif isinstance(v, int) and v in (0, 1):
            out.append(bool(v))
        else:
            raise ValueError(f"{name} must hold booleans (or 0/1), found {v!r}")
    return out


# Header keywords with no value to pass on.
_COMMENTARY = {"", "COMMENT", "HISTORY", "CONTINUE", "END"}


def _header(h: Any) -> Any:
    """A JSON-ready header: a dict, or a list of 80-character cards.

    An ``astropy.io.fits.Header`` is passed as a dict of its values: Astropy has already joined
    ``CONTINUE`` strings and unwrapped ``HIERARCH`` keywords, which card images would lose.
    """
    if h is None or isinstance(h, dict):
        return h
    if hasattr(h, "cards") and hasattr(h, "items"):  # astropy.io.fits.Header
        out = {}
        for k, v in h.items():
            if k.upper() in _COMMENTARY:
                continue
            if hasattr(v, "item"):
                v = v.item()
            if v is None or isinstance(v, (bool, int, float, str)):
                out[k] = v
        return out
    if isinstance(h, str):  # a header block
        return h.splitlines() if "\n" in h else [h[i : i + 80] for i in range(0, len(h), 80)]
    return [str(c) for c in h]


def _json_scalar(v: Any) -> Any:
    if hasattr(v, "item"):  # NumPy scalar
        v = v.item()
    if isinstance(v, float) and not math.isfinite(v):
        raise ValueError(f"non-finite option value {v!r}")
    return v


# Options of source detection (images only).
_DETECT_OPTIONS = ("fwhm", "detect_threshold", "saturation", "gain", "max_sources")


def _image_args(image: Any) -> "tuple[Optional[str], Optional[bytes], int, int, Any]":
    """(path, float32 bytes, width, height, header) for the native calls.

    ``image`` is a FITS path, a 2-D array (NumPy, masked entries blank), or an object with
    ``data`` and ``header`` such as an ``astropy.io.fits`` HDU.
    """
    if isinstance(image, (str, os.PathLike)):
        return os.fspath(image), None, 0, 0, None
    header = None
    if hasattr(image, "data") and hasattr(image, "header") and not hasattr(image, "ndim"):
        header, image = image.header, image.data
    import numpy as np  # only needed for arrays

    a = np.ma.filled(np.ma.asarray(image, dtype=np.float64), np.nan) if np.ma.isMaskedArray(
        image
    ) else np.asarray(image, dtype=np.float64)
    if a.ndim != 2:
        raise ValueError(f"image must be 2-D, not shape {a.shape}")
    h, w = a.shape
    data = np.ascontiguousarray(a, dtype="<f4").tobytes()
    return None, data, int(w), int(h), header


def detect(image: Any, **options: Any) -> DetectReport:
    """Find stars in an image: a FITS path, a 2-D array (row 0 = FITS row 1) or an Astropy HDU.

    Options: ``fwhm`` (px; default measured), ``detect_threshold`` (σ of a peak; default: the
    level that lets one noise peak through per frame, about 5), ``saturation`` and ``gain``
    (default: the header's ``SATURATE``, ``GAIN``/``EGAIN``), ``max_sources`` (default 2000),
    ``header`` (dict or cards, for an array). Returns the report ``fayde detect`` prints:
    settings, background, noise, ``sources`` (the stars, with FITS 1-based ``x``, ``y``,
    ``sigma_px``, ``flux``, ``p_star`` and ``class``), ``streaks`` (ends, angle, width, flux)
    and ``rejected`` (hot pixels, cosmic-ray hits, anything on a streak).
    """
    path, data, w, h, header = _image_args(image)
    if header is not None and "header" not in options:
        options = dict(options, header=header)
    opts = _options(options)
    text = _native.detect(
        path, data, w, h, json.dumps(opts, allow_nan=False, default=_json_scalar)
    )
    return cast(DetectReport, json.loads(text))


def _options(options: Dict[str, Any]) -> Dict[str, Any]:
    opts: Dict[str, Any] = {}
    for k, v in options.items():
        if v is None:
            continue
        if k in ("header", "wcs"):
            opts[k] = _header(v)
        elif k == "instrument":
            if isinstance(v, dict):
                opts["instrument"] = v
            else:
                opts["instrument_path"] = os.fspath(v)
        else:
            opts[k] = _json_scalar(v)
    return opts


class Solver:
    """A UCAC4 catalogue and pattern indexes, opened once and reused.

    ``ucac4`` is the UCAC4 directory (holding ``u4b/``, optionally ``u4i/``); ``indexes`` the
    ``.fidx`` files used by :meth:`solve`. Thread-safe: the GIL is released while solving.
    """

    def __init__(self, ucac4: "str | os.PathLike[str]", indexes: Iterable["str | os.PathLike[str]"] = ()):
        self._native = _native.Solver(os.fspath(ucac4), [os.fspath(p) for p in indexes])

    def solve(self, x: Sequence[float], y: Sequence[float], **options: Any) -> Report:
        """Solve a source list. Returns the report; ``report["solved"]`` says whether it solved
        (an abstention has ``reason``, ``best_log_odds``, ``hypotheses_tested``, ``search``)."""
        return self._run("solve", x, y, options)

    def solve_image(self, image: Any, **options: Any) -> Report:
        """Find the stars in an image and solve: ``image`` is a FITS path, a 2-D array or an
        Astropy HDU. The image's header gives the priors unless ``header`` is passed. Takes
        :meth:`solve`'s options and :func:`detect`'s; the report adds ``detection``."""
        path, data, w, h, header = _image_args(image)
        if header is not None and "header" not in options:
            options = dict(options, header=header)
        text, error = self._native.run_image(
            "solve",
            path,
            data,
            w,
            h,
            json.dumps(self._options(options), allow_nan=False, default=_json_scalar),
        )
        report = cast(Report, json.loads(text))
        if error:
            raise FaydeError(cast(Dict[str, Any], report))
        return report

    def verify(self, x: Sequence[float], y: Sequence[float], wcs: Any, **options: Any) -> Report:
        """Check a WCS (dict, cards or astropy Header) against the catalogue, unchanged."""
        return self._run("verify", x, y, dict(options, wcs=wcs))

    def refit(self, x: Sequence[float], y: Sequence[float], wcs: Any, **options: Any) -> Report:
        """Re-fit a WCS against the catalogue."""
        return self._run("refit", x, y, dict(options, wcs=wcs))

    def stars_at(
        self,
        wcs: Any,
        x: Any = None,
        y: Any = None,
        *,
        ra: Any = None,
        dec: Any = None,
        sigma_px: Any = None,
        sigma_arcsec: Any = None,
        **options: Any,
    ) -> AtReport:
        """Is there a catalogue star here? For pixel positions ``x``, ``y`` and/or sky positions
        ``ra``, ``dec`` (degrees) under ``wcs``, the nearest catalogue stars at the epoch, with
        separation, combined 1σ, ``p_consistent`` (Rayleigh tail: chance of this offset or more
        if the position is the star), ``p_chance`` (chance that some star lies this close to an
        arbitrary position) and ``consistent``.

        Options: ``epoch``, ``header`` (for the epoch), ``pixel_error`` (default pixel 1σ),
        ``wcs_sigma_arcsec``, ``search_radius_arcsec`` (default 30), ``max_stars`` (default 3),
        ``origin``. Positions come back in input order, pixels first.
        """
        if (x is None) != (y is None) or (ra is None) != (dec is None):
            raise ValueError("give x with y, and ra with dec")
        opts = self._options(dict(options, wcs=wcs))
        text, error = self._native.at(
            _required("x", x),
            _required("y", y),
            _floats("sigma_px", sigma_px),
            _required("ra", ra),
            _required("dec", dec),
            _floats("sigma_arcsec", sigma_arcsec),
            json.dumps(opts, allow_nan=False, default=_json_scalar),
        )
        at = cast(AtReport, json.loads(text))
        if error:
            raise FaydeError(cast(Dict[str, Any], at))
        return at

    def _options(self, options: Dict[str, Any]) -> Dict[str, Any]:
        return _options(options)

    def _run(self, mode: str, x: Any, y: Any, options: Dict[str, Any]) -> Report:
        flux = _floats("flux", options.pop("flux", None))
        sigma_px = _floats("sigma_px", options.pop("sigma_px", None))
        p_star = _floats("p_star", options.pop("p_star", None))
        truncated = _bools("truncated", options.pop("truncated", None))
        opts = self._options(options)
        text, error = self._native.run(
            mode,
            _required("x", x),
            _required("y", y),
            flux,
            sigma_px,
            p_star,
            truncated,
            json.dumps(opts, allow_nan=False, default=_json_scalar),
        )
        report = cast(Report, json.loads(text))
        if error:
            raise FaydeError(cast(Dict[str, Any], report))
        return report
