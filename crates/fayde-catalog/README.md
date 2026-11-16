# fayde-catalog

The `fayde-catalog` crate provides star catalog access for the
[Fayde](https://github.com/shaner306/Fayde) astrometric plate solver. It contains the `Catalog`
trait, an in-memory catalog, a reader for the original UCAC4 files, and proper-motion
propagation with error growth.

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
