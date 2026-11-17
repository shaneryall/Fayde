# Develop Fayde

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

Continuous integration runs code, documentation, packaging, and platform checks. To check a
change locally, run these commands:

```bash
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test
python3 scripts/check_release.py
python3 scripts/gen_third_party_licenses.py --check
```

After a change to `Cargo.lock`, run `python3 scripts/gen_third_party_licenses.py` without
`--check` to update the third-party license file.

The test vectors in `tests/fixtures/` are synthetic. The scripts in `scripts/` generate them
with Astropy.

Some tests compare Fayde with the real catalog. These tests use `#[ignore]` because the catalog
isn't in the repository. To include them, follow these steps:

1. Put a copy of UCAC4 at `data/ucac4/`, or set the `FAYDE_UCAC4_DIR` environment variable to
   the path to a copy.

2. Put `ucac4_oracle.json` in the same directory. This file holds cone-search results from an
   independent UCAC4 reader. The top of `crates/fayde-catalog/tests/reader_oracle.rs` describes
   its format.

3. Run every test:

   ```bash
   cargo test -- --include-ignored
   ```

The Python tests in `crates/fayde-py/tests` run with `pytest` against an installed wheel. They
skip the tests that need the catalog when it's missing.

For more information about the algorithms, see the [references](references.md).
For more information about the architecture, see the [design document](design/plate-solver.md).

## Versioning and releases

Fayde follows semantic versioning. While the version is 0.x, a breaking change increases the
minor number. The following changes are breaking:

- Removing a command-line or Python option, or changing its meaning
- Removing or renaming a JSON report field
- Changing the index file format or the bytes of an index, including its `content_hash`

Adding a report field isn't a breaking change. If a release changes indexes, you must rebuild
your indexes as directed in the release notes.

A fix that leaves reports and indexes unchanged is a patch release. So is a new command,
option, or report field that changes nothing that already exists.

A tag `vX.Y.Z` that matches the workspace version starts the release pipeline. The pipeline
builds the `fayde` command and the Python wheels for Linux, macOS, and Windows, runs each one on
its own platform, and then publishes them as a GitHub release with SHA-256 checksums. Publish
the same version to PyPI and crates.io. For release instructions, see
[Release Fayde](releasing.md).
