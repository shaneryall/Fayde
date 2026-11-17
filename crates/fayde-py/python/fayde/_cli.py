"""The ``fayde`` command, installed with the package.

It runs the same Rust code as the ``fayde`` binary from a release.
"""

from __future__ import annotations

import signal
import sys

from . import _native


def main() -> None:
    """Run ``fayde`` with the command-line arguments and exit with its status."""
    # The command runs in Rust without the interpreter lock, so Python's own Ctrl-C handler
    # would not run until it ends. Restore the default, which stops the process at once, as
    # it does the binary.
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    sys.exit(_native.cli_main(sys.argv[1:]))


if __name__ == "__main__":
    main()
