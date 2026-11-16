# fayde-fits

The `fayde-fits` crate handles FITS input and output for the
[Fayde](https://github.com/shaner306/Fayde) astrometric plate solver. It reads FITS headers
and images, reads and writes source lists, derives priors from a frame header, reads instrument
profiles, and writes `.wcs` files.

For more information about Fayde and how to install the `fayde` command or Python package, see
the [Fayde README](https://github.com/shaner306/Fayde#readme).

All Fayde crates share one release version. Each crate requires exactly that version of the
other Fayde crates.

For more information about the algorithms, see the
[references](https://github.com/shaner306/Fayde/blob/main/docs/references.md).

## License

Fayde uses the Apache License 2.0. For license terms and attribution, see the `LICENSE` and
`NOTICE` files in this crate. If you use Fayde in a publication, cite it as described in
[`CITATION.cff`](https://github.com/shaner306/Fayde/blob/main/CITATION.cff).
