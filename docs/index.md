# Fayde

Fayde finds the stars in a FITS image, matches them against a star catalog, and returns a world
coordinate system (WCS) solution with its measured accuracy. It runs from the command line and
from Python on Linux, macOS, and Windows.

To solve your first image, follow these pages in order:

1. [Install Fayde](install.md).
2. [Choose and connect a star catalog](catalog.md).
3. [Solve a frame](usage.md#solve-a-frame).

## How-to guides

Follow the steps for a task:

- [Install Fayde](install.md)
- [Choose and connect a star catalog](catalog.md)
- [Use the command line](usage.md)
- [Use Fayde from Python](python.md)
- [Build pattern indexes](index-set.md)

## Reference

Look up fields, formats, and rules. For the options of each command, run the command with
`--help`.

- [Python API](reference/python-api.md)
- [Source list format](source-list.md)
- [Versioning](versioning.md)
- [Algorithm sources](references.md)

## Explanation

Read how the solver and the detector work, and why they work that way:

- [Plate solver design](design/plate-solver.md)
- [Star detection](design/detection.md)
- [Literature review](research/literature-review.md)
