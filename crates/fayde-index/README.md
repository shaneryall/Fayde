# fayde-index

The `fayde-index` crate builds and reads pattern indexes for the
[Fayde](https://github.com/shaneryall/Fayde) astrometric plate solver. It contains uniform star
selection on a HEALPix grid, a k-d tree, the 4-star pattern hash table, and the versioned index
file with its content hash.

For more information about Fayde and how to install the `fayde` command or Python package, see
the [Fayde README](https://github.com/shaneryall/Fayde#readme).

All Fayde crates share one release version. Each crate requires exactly that version of the
other Fayde crates.

For more information about the algorithms, see the
[references](https://github.com/shaneryall/Fayde/blob/main/docs/references.md).

## License

Fayde uses the Apache License 2.0. For license terms and attribution, see the `LICENSE` and
`NOTICE` files in this crate. If you use Fayde in a publication, cite it as described in
[`CITATION.cff`](https://github.com/shaneryall/Fayde/blob/main/CITATION.cff).
