# Source list format

A source list contains the positions of sources that you detect in one image. Use a source
list when you detect the stars yourself. When you give Fayde the image with `--image`, you
don't need a source list.

Fayde reads a source list with the `--sources` option of `fayde solve`, `fayde verify`, and
`fayde refit`. The `fayde detect` command writes one with `-o`.

## File structure

A source list is a FITS file with a binary table extension. Each row of the table describes one
source. Fayde reads the first binary table that has `X` and `Y` columns and ignores columns
that this page doesn't list. Fayde recognizes a source list by this table, not by the filename.

Name a source list after its image, with the suffix `.sources.fits`. For example, use
`frame.sources.fits` for the source list from `frame.fits`.

The `fayde detect` command sets the table's `EXTNAME` keyword to `SOURCES` and records the Fayde
version in the `CREATOR` keyword. Fayde doesn't require these keywords in a source list that
you write yourself.

## Columns

The following table lists the columns that Fayde reads:

| Column | Required | Content |
|---|---|---|
| `X` | Yes | The horizontal position of the source, in pixels |
| `Y` | Yes | The vertical position of the source, in pixels |
| `FLUX` | No | The integrated flux, in any linear unit. Fayde uses it to order the sources by brightness |
| `SIGMA` or `SIGMA_PX` | No | The 1σ centroid error, in pixels. Fayde ignores a value that isn't positive |
| `P_STAR` | No | The prior probability that the source is a real star, from 0 to 1. Fayde ignores a value outside that range |
| `TRUNCATED` | No | Whether the frame edge or a mask cuts the source. A logical or integer column, where a nonzero value means true |

Positions follow the FITS pixel convention: the center of the first pixel is at (1, 1), and `X`
increases along a row of the image.

When a source has no `SIGMA` value, Fayde uses the value of `--pixel-error`, which is 1 pixel
by default.

## Image size

Fayde needs the width and height of the image. In the table header or the primary header, set
the `IMAGEW` and `IMAGEH` keywords to the dimensions in pixels. If the source list doesn't have
these keywords, Fayde takes the size from the frame header (`--header`) or from the `--width`
and `--height` options.

## Write a source list

In this example, `x`, `y`, and `flux` are sequences with one value for each source, and the
image is 4,096 by 4,096 pixels. To write and solve a source list, follow these steps:

1. In Python, with the `astropy` package installed, write the source list with an Astropy table:

   ```python
   from astropy.table import Table

   table = Table({"X": x, "Y": y, "FLUX": flux})
   table.meta["IMAGEW"] = 4096
   table.meta["IMAGEH"] = 4096
   table.write("frame.sources.fits", format="fits")
   ```

2. In a terminal, solve the source list with `--sources`:

   ```bash
   fayde solve --sources frame.sources.fits --header FRAME --ucac4 UCAC4_DIR -o OUTPUT
   ```

   Replace the following:

   - `FRAME`: The path to the FITS image
   - `UCAC4_DIR`: The path to the directory that contains `u4b`
   - `OUTPUT`: The path to the WCS file to write

   For more information about solving a frame, see [Solve a frame](usage.md#solve-a-frame).
