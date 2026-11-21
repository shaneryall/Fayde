# fayde-cli

The `fayde-cli` crate provides the `fayde` command for the
[Fayde](https://github.com/shaneryall/Fayde) astrometric plate solver. You provide a FITS image of
a star field. The command finds the stars, matches them against the UCAC4 star catalog, and
returns a verified world coordinate system (WCS) solution with its measured accuracy.

Fayde reads the UCAC4 star catalog from a local copy. The catalog isn't included with this
crate. To download it, see
[Get the UCAC4 catalog](https://github.com/shaneryall/Fayde/blob/main/docs/catalog.md).

## Install the command

You need Rust 1.85 or later to build the command. In a terminal, build and install it with Cargo:

```bash
cargo install --locked fayde-cli
```

To install a prebuilt binary for Linux, macOS, or Windows without Rust, see
[Install the command-line tool](https://github.com/shaneryall/Fayde/blob/main/docs/install.md#install-the-command-line-tool).

## Use the command

To solve a frame, pass the FITS image and the catalog to `fayde solve`:

```bash
fayde solve --image FRAME --ucac4 UCAC4_DIR -o OUTPUT
```

Replace the following:

- `FRAME`: The path to the FITS image, for example `frame.fits`
- `UCAC4_DIR`: The path to the directory that contains `u4b`
- `OUTPUT`: The path to the WCS file to write, for example `frame.wcs`

The command prints a JSON report to standard output.

To list the commands, run `fayde --help`. To list a command's options and default values, run
the command with `--help`, for example `fayde solve --help`.

For more information, see the following pages:

- Priors and the report:
  [Use the command line](https://github.com/shaneryall/Fayde/blob/main/docs/usage.md)
- Source lists:
  [Source list format](https://github.com/shaneryall/Fayde/blob/main/docs/source-list.md)
- Pattern indexes:
  [Pattern indexes for sub-degree fields](https://github.com/shaneryall/Fayde/blob/main/docs/index-set.md)
- Versioning:
  [Versioning](https://github.com/shaneryall/Fayde/blob/main/docs/versioning.md)

## License

Fayde uses the GNU General Public License, version 3 or later. For license terms and
attribution, see the `LICENSE` and `NOTICE` files in this crate. If you use Fayde in a
publication, cite it as described in
[`CITATION.cff`](https://github.com/shaneryall/Fayde/blob/main/CITATION.cff).
