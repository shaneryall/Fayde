"""Stub of the compiled module (crates/fayde-py/src/lib.rs). Use `fayde.Solver`, not this."""

from typing import List, Optional, Tuple

__version__: str
commit: str

def detect(
    path: Optional[str], data: Optional[bytes], width: int, height: int, options_json: str
) -> str: ...

class Solver:
    def __init__(self, ucac4: str, indexes: List[str] = ...) -> None: ...
    def run(
        self,
        mode: str,
        x: List[float],
        y: List[float],
        flux: Optional[List[float]],
        sigma_px: Optional[List[float]],
        p_star: Optional[List[float]],
        truncated: Optional[List[bool]],
        options_json: str,
    ) -> Tuple[str, bool]: ...
    def run_image(
        self,
        mode: str,
        path: Optional[str],
        data: Optional[bytes],
        width: int,
        height: int,
        options_json: str,
    ) -> Tuple[str, bool]: ...
    def at(
        self,
        x: List[float],
        y: List[float],
        sigma_px: Optional[List[float]],
        ra: List[float],
        dec: List[float],
        sigma_arcsec: Optional[List[float]],
        options_json: str,
    ) -> Tuple[str, bool]: ...
