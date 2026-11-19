"""The ``fayde`` command, installed with the package.

It runs the same Rust code as the ``fayde`` binary from a release.
"""

from __future__ import annotations

import signal
import sys

from . import _native


def main() -> None:
    """Run ``fayde`` with the command-line arguments and exit with its status."""
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    sys.exit(_native.cli_main(sys.argv[1:]))


if __name__ == "__main__":
    main()
