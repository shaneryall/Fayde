"""Types of the reports `fayde.Solver` returns (the CLI's JSON reports).

Every report is a plain ``dict``; these ``TypedDict`` classes describe its keys for type
checkers. ``total=False``: which keys are present depends on the outcome (solved, abstained,
error) and the command. Pixel fields are FITS 1-based.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple, TypedDict


class Stage(TypedDict):
    """One stage of the search."""

    stage: str
    hypotheses: int
    seconds: float
    outcome: str


class Pair(TypedDict):
    """A matched source and catalogue star."""

    source: int
    x: float
    y: float
    star_id: int
    ra_deg: float
    dec_deg: float
    mag: Optional[float]
    predicted_x: float
    predicted_y: float
    residual_px: Tuple[float, float]
    inlier_prob: float


class Covariance(TypedDict):
    """Fit covariance, per tangent-plane axis, rad², over monomials of (u/norm, v/norm)."""

    terms: List[Tuple[int, int]]
    matrix_rad2: List[float]
    norm_px: float
    crpix_fits: Tuple[float, float]


class FaydeIdentity(TypedDict):
    version: str
    commit: str


class CatalogIdentity(TypedDict):
    name: str
    release: str
    location: str
    reference_epoch: str
    solve_epoch_julian_year: Optional[float]


class IndexIdentity(TypedDict):
    path: str
    content_hash: str
    fov_deg: float
    stars: int
    patterns: int
    epoch_year: float


class Identity(TypedDict):
    """What produced a report."""

    fayde: FaydeIdentity
    catalog: CatalogIdentity
    indexes: List[IndexIdentity]


class Report(TypedDict, total=False):
    """`solve`, `refit` and `verify` reports."""

    solved: bool
    error: bool
    reason: str
    warnings: List[str]
    method: str
    crval_deg: Tuple[float, float]
    pixel_scale_arcsec: float
    roll_deg: float
    sip_order: int
    parity: str
    centre_sigma_arcsec: float
    worst_sigma_arcsec: float
    loo_rms_px: Optional[float]
    log_odds: float
    nfa: float
    hypotheses_tested: int
    matches: int
    rms_px: float
    error_scale: float
    index_hash: Optional[str]
    seconds: float
    prior_violated: bool
    violated_priors: List[str]
    solved_by_stage: Optional[str]
    relaxed_priors: bool
    search: List[Stage]
    wcs_cards: List[str]
    covariance: Covariance
    pairs: List[Pair]
    best_log_odds: Optional[float]
    accepted: bool
    threshold: float
    confident_matches: int
    references: int
    sources: int
    detection: Dict[str, Any]
    priors: Dict[str, Any]
    identity: Identity


class AtStar(TypedDict):
    star_id: int
    ra_deg: float
    dec_deg: float
    mag: Optional[float]
    pixel: Optional[Tuple[float, float]]
    separation_arcsec: float
    sigma_arcsec: float
    p_consistent: float
    p_chance: float


class AtPosition(TypedDict):
    index: int
    ra_deg: float
    dec_deg: float
    pixel: Optional[Tuple[float, float]]
    position_sigma_arcsec: float
    consistent: bool
    stars: List[AtStar]


class AtReport(TypedDict, total=False):
    """`stars_at` report."""

    positions: List[AtPosition]
    consistent: int
    consistent_p: float
    options: Dict[str, Any]
    identity: Identity
    error: bool
    reason: str


class DetectedSource(TypedDict):
    """A detection reported by :func:`fayde.detect` (FITS 1-based pixels).

    The detection's class is under the key ``"class"`` (a Python keyword, so not declared
    here): ``"star"``, ``"saturated"`` or ``"extended"`` for sources, ``"sharp"`` (a hot pixel
    or cosmic-ray hit) or ``"streak"`` for rejected detections.
    """

    x: float
    y: float
    sigma_px: float
    flux: float
    snr: float
    p_star: float
    fwhm: float
    sigma_major_px: Optional[float]
    sigma_minor_px: float
    theta_deg: float
    streak: Optional[int]
    truncated: bool
    saturated: bool


class DetectedStreak(TypedDict):
    """A traced streak. Ends are FITS 1-based pixels."""

    x0: float
    y0: float
    x1: float
    y1: float
    length_px: float
    angle_deg: float
    width_sigma_px: float
    flux: float
    snr: float
    coverage: float
    touches_edge: bool
    cosmic_ray: bool
    photometry_masked: bool


class DetectImage(TypedDict):
    width: int
    height: int
    hdu: int


class DetectReport(TypedDict, total=False):
    """:func:`fayde.detect` report."""

    image: DetectImage
    fwhm_px: float
    fwhm_estimated: bool
    background: float
    noise: float
    threshold_sigma: float
    saturation: Optional[float]
    gain: Optional[float]
    count: int
    classes: Dict[str, int]
    rejected_narrow: int
    bin: int
    kernel_sigma_px: float
    passes: int
    undefined_pixels: int
    seconds: float
    streaks: List[DetectedStreak]
    warnings: List[str]
    sources: List[DetectedSource]
    rejected: List[DetectedSource]
