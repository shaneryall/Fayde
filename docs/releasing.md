# Release Fayde

To release Fayde, push a Git tag `vX.Y.Z` for a commit on the `main` branch. GitHub Actions
builds, tests, and publishes the release files. For more information about versioning, see
[Versioning and releases](development.md#versioning-and-releases).

Replace `X.Y.Z` in the commands on this page with the release version, for example `0.3.0`.

## Prepare the release commit

Before you create the release tag, update and check the release files:

1. In `Cargo.toml`, set the new version in two places:

   - The `version` field under `[workspace.package]`
   - Every `=X.Y.Z` requirement under `[workspace.dependencies]`

   All Fayde crates share one release version. Each crate requires exactly that version of the
   other Fayde crates.

2. In `CITATION.cff`, set `version` and `date-released`.

3. In `CHANGELOG.md`, rename the `## Unreleased` heading to `## X.Y.Z (DATE)`.
   Replace `DATE` with the release date in the form `2026-11-16`. This section becomes the
   release notes.

4. Update `Cargo.lock`:

   ```bash
   cargo check
   ```

5. If the previous step changed a third-party crate in `Cargo.lock`, update the third-party
   license file:

   ```bash
   python3 scripts/gen_third_party_licenses.py
   ```

6. If the release changes the bytes of an index, rebuild the index set.

7. If you rebuilt the index set, update the hashes in [`docs/index-set.md`](index-set.md).

8. Check that the release files agree with each other:

   ```bash
   python3 scripts/check_release.py --tag vX.Y.Z
   ```

9. Run every test, including the tests that need the UCAC4 catalog:

   ```bash
   cargo test -- --include-ignored
   ```

   Continuous integration doesn't run the catalog tests, because the catalog isn't in the
   repository.

10. Commit the changes.

## Publish the release on GitHub

From the prepared release commit, create and push the tag:

1. Tag the release commit:

   ```bash
   git tag -a vX.Y.Z -m "Fayde X.Y.Z"
   ```

2. Push the branch and the tag:

   ```bash
   git push origin main vX.Y.Z
   ```

The tag starts the [release workflow](../.github/workflows/release.yml), which performs the
following actions:

1. Runs the same checks as continuous integration.
2. Builds the `fayde` command for Linux (x86_64 and 64-bit Arm, glibc 2.17 or later), macOS
   (Apple silicon and Intel), and Windows (x86_64). Each archive also contains the license files,
   the index guide, and the example instrument profile.
3. Builds one Python wheel for each platform and one source distribution.
4. Runs each binary and installs and tests each wheel on its own platform. It also builds a
   wheel from the source distribution and tests it.
5. Creates the GitHub release with a `SHA256SUMS` file after every required check passes.
6. Downloads the prebuilt `fayde` binary from the new release on every platform, then checks
   and installs it with the commands from the [install guide](install.md).

A workflow run that fails before step 5 publishes nothing. Don't move the tag of a release
that the workflow has already published. To retry a release that failed before publication,
follow these steps:

1. Fix the cause of the failure.

2. Commit the fix.

3. Move the tag to the fixed commit:

   ```bash
   git tag -f -a vX.Y.Z -m "Fayde X.Y.Z"
   ```

4. Push the branch:

   ```bash
   git push origin main
   ```

5. Replace the tag on GitHub with the moved tag:

   ```bash
   git push --force origin vX.Y.Z
   ```

Archive names identify the platform without a version, for example
`fayde-linux-x86_64.tar.gz`. The download URLs in the [install guide](install.md) depend on
these names. If you change a name in the workflow, update the install guide too.

## Publish to PyPI and crates.io

Publish every release to both registries. The install commands in the
[install guide](install.md) use PyPI and crates.io, so they install a version only after you
publish it.

You can't replace a version after you publish it to PyPI or crates.io. After you check the
files in the GitHub release, start the separate publishing workflow.

### Set up the registries

Before you publish to a registry for the first time, configure its publishing credentials:

1. In the repository on GitHub, go to **Settings > Environments**.

2. Create an environment named `pypi`.

3. Create an environment named `crates-io`.

   Use the names `pypi` and `crates-io` exactly as shown. The publishing workflow uses these
   environments. It can't access a secret that you store in an environment with another name,
   such as `crates.io`.

4. Optional: Add yourself as a required reviewer to both environments. GitHub then asks for
   your approval before a publishing job runs.

5. On PyPI, add a trusted publisher for the project `fayde` with this repository, the workflow
   `publish.yml`, and the environment `pypi`. Trusted publishing doesn't require you to store
   a PyPI token.

6. On crates.io, create an API token.

7. In the `crates-io` environment, add the token as a secret named `CARGO_REGISTRY_TOKEN`.

### Publish a release to the registries

Publish the files from an existing GitHub release:

1. Check the files in the GitHub release.
2. In the repository on GitHub, go to **Actions > Publish to registries**.
3. Click **Run workflow**.
4. Enter the release tag.
5. Select the registries to publish to.
6. Click **Run workflow**.

The [publishing workflow](../.github/workflows/publish.yml) publishes these files:

- PyPI receives the wheels and source distribution from the GitHub release, which the release
  workflow tested.
- The crates.io registry receives the crates from the tagged commit in dependency order. The
  Python package, `fayde-py`, goes to PyPI.

The first release publishes nine new crates. From one account, crates.io accepts five new
crates at once and then one every 10 minutes. The first run therefore takes about 40 minutes.
The workflow waits and retries when crates.io limits publishing. If a run stops, start the
workflow again. It skips crate versions that crates.io already has.

After each registry accepts the version, the workflow installs and checks that version. It uses
`pip install` on every supported platform and `cargo install` on Linux.

## Release checks

The `scripts/check_release.py` script compares the workspace version with the following:

- The version that each Fayde crate requires of the others
- The version in `CITATION.cff`
- The sections of the changelog
- The tag, when you pass `--tag`

The script also compares the copies of `LICENSE` and `NOTICE` in each crate with the originals.
Each crate has copies because a published crate or wheel can include only files from its own
directory. Continuous integration runs the script on every push.
