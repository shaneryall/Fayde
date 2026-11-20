# Contribute to Fayde

Bug reports, questions, and changes are welcome. Everyone who takes part follows the
[code of conduct](CODE_OF_CONDUCT.md).

## Report a bug or ask a question

To report a bug or ask a question, open an issue on GitHub. For a
failed or wrong solve, include the command, the JSON report, the Fayde version
(`fayde --version`), and the image or source list if you can share it.

To report a security problem, don't open an issue. See the [security policy](SECURITY.md).

## Build Fayde from source

Fayde is a Cargo workspace. All crates share one release version. The workspace contains
the following crates:

| Crate | Contents |
|---|---|
| `fayde-core` | Sky positions, unit vectors, pixel conventions |
| `fayde-wcs` | TAN and SIP projection, and FITS WCS header cards |
| `fayde-catalog` | The `Catalog` trait, an in-memory catalog, the UCAC4 reader, and proper-motion propagation |
| `fayde-index` | HEALPix, k-d tree, uniform star selection, the 4-star pattern hash table, and the index file |
| `fayde-solve` | Hypothesis generation, verification, and the WCS fit |
| `fayde-fits` | FITS reader, images, source lists, priors from a header, and the `.wcs` writer |
| `fayde-detect` | Star and streak detection for FITS images |
| `fayde-api` | Request building and the JSON report, shared by the command line and Python |
| `fayde-cli` | The `fayde` command |
| `fayde-py` | The Python package |

To build Fayde from a clone of the repository, you need Git and Rust 1.85 or later. Follow
these steps:

1. In a terminal, clone the repository:

   ```bash
   git clone https://github.com/shaneryall/Fayde
   ```

2. Go to the repository directory:

   ```bash
   cd Fayde
   ```

3. Build and install the `fayde` command:

   ```bash
   cargo install --locked --path crates/fayde-cli
   ```

4. Optional: To build and install the Python package, activate a virtual environment.

5. Optional: In the active virtual environment, build and install the Python package with `pip`:

   ```bash
   python -m pip install ./crates/fayde-py
   ```

## Check a change

Continuous integration runs code, documentation, packaging, and platform checks. To check a
change locally, run these commands:

```bash
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test
RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps
python3 scripts/check_release.py
python3 scripts/gen_third_party_licenses.py --check
```

After a change to `Cargo.lock`, run `python3 scripts/gen_third_party_licenses.py` without
`--check` to update the third-party license file. Continuous integration also checks the
licenses, advisories, and sources of every dependency against `deny.toml` with
[cargo-deny](https://github.com/EmbarkStudios/cargo-deny), and measures test coverage.

The `fuzz` directory holds fuzz targets for the readers of files that may come from other
people. Continuous integration runs each target every night. To run one locally, install nightly
Rust and [cargo-fuzz](https://github.com/rust-fuzz/cargo-fuzz), and then run this command in the
`fuzz` directory:

```bash
cargo +nightly fuzz run fits_image
```

The other targets are `wcs_header`, `pattern_index`, `source_list`, and
`instrument_profile`.

Continuous integration also checks the Python package against an installed wheel. To run the
same checks, follow these steps:

1. In the virtual environment where you installed the Python package, reinstall it to include
   your change:

   ```bash
   python -m pip install ./crates/fayde-py
   ```

2. Install the test tools:

   ```bash
   python -m pip install pytest numpy mypy
   ```

3. Check the package's types:

   ```bash
   (cd crates/fayde-py/python && python -m mypy --strict fayde)
   ```

4. Run the Python tests:

   ```bash
   python -m pytest -q crates/fayde-py/tests
   ```

   The tests skip the catalog cases when the oracle file is missing. For these cases, use the
   catalog and oracle setup described in the steps below.

Astropy generates the synthetic FITS fixtures in `tests/fixtures/fits/` and the TAN vectors
in `crates/fayde-wcs/tests/fixtures/`. The vector fixture ships with the WCS crate.

Some tests compare Fayde with the real catalog. These tests use `#[ignore]` because the catalog
isn't in the repository. To include them, follow these steps:

1. Put a copy of UCAC4 at `data/ucac4/`, or set the `FAYDE_UCAC4_DIR` environment variable to
   the path to a copy.

2. Put `ucac4_oracle.json` in the same directory. This file holds cone-search results from an
   independent UCAC4 reader. The `Oracle`, `Cone`, and `OracleStar` structures in
   `crates/fayde-catalog/tests/reader_oracle.rs` define the JSON fields.

3. Run every test:

   ```bash
   cargo test -- --include-ignored
   ```

## Build the documentation site

The documentation site is built from the pages in `docs/` with
[Material for MkDocs](https://squidfunk.github.io/mkdocs-material/). Its command-line reference
comes from the `--help` text of a built `fayde` command. To build the site and serve it on your
computer, follow these steps:

1. Build the `fayde` command:

   ```bash
   cargo build --release -p fayde-cli
   ```

2. In a virtual environment, install the site tools at the versions that continuous integration
   uses:

   ```bash
   python -m pip install --require-hashes -r scripts/docs/requirements.txt
   ```

   `scripts/docs/requirements.in` lists the tools. To update them, change that file and run
   `uv pip compile --universal --generate-hashes --python-version 3.13 --no-strip-markers
   scripts/docs/requirements.in -o scripts/docs/requirements.txt`.

3. Serve the site:

   ```bash
   FAYDE_BIN=target/release/fayde mkdocs serve
   ```

   The output gives the local address of the site. Continuous integration builds the site with
   `mkdocs build --strict`, so a broken link or a warning fails the build.

## License of contributions

Fayde is licensed under the GNU General Public License, version 3 or later. When you submit a
contribution (a pull request, a patch, or any other change), you agree to the following terms:

1. You wrote the contribution, or you have the right to submit it under these terms.
2. Your contribution is licensed under the same license as Fayde.
3. You grant Shane Ryall, Fayde's copyright holder, a perpetual, worldwide, non-exclusive,
   royalty-free, irrevocable license to use, modify, and distribute your contribution under any
   license, including permissive and commercial ones. You keep the copyright in your
   contribution.

In your pull request, say that you agree to these terms. Pull requests without that statement
can't be merged.

## Rules for changes

A change is easier to review and merge when it follows these rules:

- Keep results reproducible. Fayde gives the same numbers on every platform. A change that
  alters a report or the bytes of an index says so in its description and updates the pinned
  tests on purpose. A change for speed alone leaves every report unchanged.
- Write no `unsafe` code. Every crate except the Python bindings (`fayde-py`) forbids it.
- Name the source of every algorithm. Add the source to the
  [references](docs/references.md) and cite it in the module's doc comment. If an algorithm has
  no published source, say so there.
- Add a test that fails without the change.
- Write documentation in the style of the
  [Google developer documentation style guide](https://developers.google.com/style).
- Describe a user-visible change in the `Unreleased` section of the [changelog](CHANGELOG.md).
- Start each commit subject with its type: `feat`, `bug`, `docs`, `refactor`, `test`, `chore`,
  `ci`, or `perf`. For example: `bug: keep search priors out of WCS verification`.

For more information about the algorithms, see the [references](docs/references.md).
For more information about the architecture, see the
[design document](docs/design/plate-solver.md). For the versioning rules, see
[Versioning](docs/versioning.md). For release steps, see [Release Fayde](RELEASING.md).
