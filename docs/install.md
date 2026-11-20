# Install Fayde

Use `pip install` to install both the Python module and the `fayde` command. To use the command
without Python, install the command-line tool on its own. Then
[get the UCAC4 catalog](catalog.md).

## Prerequisites

To install Fayde, use a computer that runs one of these systems:

- Linux on x86_64 or 64-bit Arm, with glibc 2.17 or later
- macOS 11 or later, on Apple silicon or an Intel processor
- Windows 10 or later, on x86_64

You also need the following:

- About 9 GB of free disk space for the UCAC4 catalog
- For the Python package, CPython 3.10 or later

To solve a frame, you also need a FITS image from your own observation. The image must meet
these requirements:

- Fayde measures point stars, straight trails that share a length and direction, and some
  curved trails. A shared curved path or field rotation can use a template learned from the
  image. See [Detect stars without solving](usage.md#detect-stars-without-solving) for
  the detection modes and how to read the result.
- The image isn't compressed with `fpack`. To decompress an image, use the `funpack` program
  from CFITSIO.
- The header gives the pointing and the pixel scale. If it doesn't, you can
  [specify them yourself](usage.md#specify-the-pointing-and-scale).

## Install the command-line tool

Each release includes a prebuilt `fayde` binary. You don't need Rust to install it. After you
install the tool, continue to [Get the UCAC4 catalog](catalog.md).

### Install on Linux or macOS

To install the latest release, follow these steps:

1. In a terminal, identify your operating system and processor architecture:

   ```bash
   uname -sm
   ```

2. Set the `platform` variable to the value that matches the output:

   | Output of `uname -sm` | Value |
   |---|---|
   | `Linux x86_64` | `linux-x86_64` |
   | `Linux aarch64` | `linux-aarch64` |
   | `Darwin arm64` | `macos-arm64` |
   | `Darwin x86_64` | `macos-x86_64` |

   For example, for the output `Darwin arm64`, set the variable to `macos-arm64`:

   ```bash
   platform=macos-arm64
   ```

3. Download the archive:

   ```bash
   curl -LO "https://github.com/shaneryall/Fayde/releases/latest/download/fayde-$platform.tar.gz"
   ```

   Optional: Before you extract the archive, check it as described in
   [Verify a download](#verify-a-download).

4. Extract the archive:

   ```bash
   tar -xzf "fayde-$platform.tar.gz"
   ```

   Extraction creates a directory named `fayde-` followed by your platform value. The
   directory contains the `fayde` binary, `README.md`, `CHANGELOG.md`, the license files
   (`LICENSE`, `NOTICE`, and `THIRD-PARTY-LICENSES.md`), and `CITATION.cff`. It also contains
   the documentation pages in its `docs` directory and the example instrument profile in its
   `profiles` directory.

5. Install the binary in `/usr/local/bin`, which is on your `PATH`:

   ```bash
   sudo mkdir -p /usr/local/bin && sudo install -m 755 "fayde-$platform/fayde" /usr/local/bin/fayde
   ```

   If `sudo` prompts for your password, enter it.

6. If macOS blocks a binary that you downloaded with a web browser instead of `curl`, remove
   the quarantine attribute to allow it to run:

   ```bash
   sudo xattr -d com.apple.quarantine /usr/local/bin/fayde
   ```

7. Check the installation:

   ```bash
   fayde --version
   ```

   The output is similar to the following:

   ```text
   fayde 0.3.0 (0123456789ab)
   ```

### Install on Windows

To install the latest release, follow these steps:

1. In PowerShell, download the archive:

   ```powershell
   curl.exe -LO https://github.com/shaneryall/Fayde/releases/latest/download/fayde-windows-x86_64.zip
   ```

   Optional: Before you extract the archive, check it as described in
   [Verify a download](#verify-a-download).

2. Extract the archive:

   ```powershell
   Expand-Archive fayde-windows-x86_64.zip -DestinationPath .
   ```

   Extraction creates a directory named `fayde-windows-x86_64`. The directory contains
   `fayde.exe`, `README.md`, `CHANGELOG.md`, the license files (`LICENSE`, `NOTICE`, and
   `THIRD-PARTY-LICENSES.md`), and `CITATION.cff`. It also contains the documentation pages in
   its `docs` directory and the example instrument profile in its `profiles` directory.

3. In the current PowerShell session, add the extracted directory to `Path`:

   ```powershell
   $env:Path += ";$PWD\fayde-windows-x86_64"
   ```

   The `Path` change lasts until you close PowerShell.

4. Check the installation:

   ```powershell
   fayde --version
   ```

   The output is similar to the following:

   ```text
   fayde 0.3.0 (0123456789ab)
   ```

Optional: To keep the `Path` change after you close PowerShell, follow these steps:

1. In the Windows search box, enter `environment variables`.
2. Select **Edit environment variables for your account**.
3. Select **Path**.
4. Click **Edit**.
5. Click **New**.
6. Enter the full path of the extracted directory.
7. Click **OK**.

### Verify a download

Each release includes a `SHA256SUMS` file that lists the SHA-256 checksum of every other file in
the release.

To check an archive on Linux or macOS, follow these steps:

1. In the directory that contains the archive, download the checksum file:

   ```bash
   curl -LO https://github.com/shaneryall/Fayde/releases/latest/download/SHA256SUMS
   ```

2. Compare the checksums. On Linux, use `sha256sum`:

   ```bash
   sha256sum --check --ignore-missing SHA256SUMS
   ```

   On macOS, use `shasum`:

   ```bash
   shasum -a 256 --check --ignore-missing SHA256SUMS
   ```

   The output names the archive, followed by `OK`.

To check the archive on Windows, follow these steps:

1. In PowerShell, in the directory that contains the archive, download the checksum file:

   ```powershell
   curl.exe -LO https://github.com/shaneryall/Fayde/releases/latest/download/SHA256SUMS
   ```

2. Print the checksum of the archive:

   ```powershell
   (Get-FileHash fayde-windows-x86_64.zip -Algorithm SHA256).Hash
   ```

3. Print the expected checksum:

   ```powershell
   Select-String fayde-windows-x86_64.zip SHA256SUMS
   ```

4. Compare the checksums from the two outputs. The archive is intact when they match.
   Uppercase and lowercase letters are equivalent.

A checksum shows that the archive is intact. Releases from 0.3.3 on also have a build
attestation, which shows that a GitHub Actions workflow of the Fayde repository built the
archive. Optional: To check it, run this [GitHub CLI](https://cli.github.com/) command in the
directory that contains the archive:

```bash
gh attestation verify ARCHIVE --repo shaneryall/Fayde
```

Replace `ARCHIVE` with the archive's filename, for example `fayde-linux-x86_64.tar.gz`. The
command fails if no attestation from the Fayde repository matches the archive.

### Install with Cargo

You need [Rust](https://www.rust-lang.org/tools/install) 1.85 or later to build the `fayde`
command from crates.io. To build and install the command, follow these steps:

1. In a terminal, install the `fayde-cli` crate:

   ```bash
   cargo install --locked fayde-cli
   ```

   Cargo installs `fayde` in its binary directory, which is `~/.cargo/bin` unless you changed it.
   If your shell can't find `fayde`, add that directory to your `PATH`.

2. Check the installation:

   ```bash
   fayde --version
   ```

   A build from crates.io shows `unknown` in place of the commit. The output is similar to the
   following:

   ```text
   fayde 0.3.0 (unknown)
   ```

## Install the Python package

The `fayde` package on PyPI includes a prebuilt wheel for every supported platform. Each wheel
supports CPython 3.10 and every later version, and you don't need Rust to install it. The
package installs both the Python module and the `fayde` command.

If PyPI has no wheel for your platform, `pip` builds the package from its source distribution.
You need Rust 1.85 or later for that build.

To install the package, follow these steps:

1. On Linux or macOS, in a terminal, create a virtual environment:

   ```bash
   python3 -m venv .venv
   ```

   On Debian and Ubuntu, this command needs the `python3-venv` package.

   On Windows, in PowerShell, create the virtual environment:

   ```powershell
   py -m venv .venv
   ```

2. On Linux or macOS, activate the virtual environment:

   ```bash
   source .venv/bin/activate
   ```

   On Windows, activate the virtual environment:

   ```powershell
   .venv\Scripts\Activate.ps1
   ```

   If PowerShell reports that running scripts is disabled, follow these steps:

   1. Allow scripts for the current session:

      ```powershell
      Set-ExecutionPolicy -Scope Process -ExecutionPolicy RemoteSigned
      ```

   2. Activate the environment again.

3. Install the package:

   ```bash
   python -m pip install fayde
   ```

4. Check the installation:

   ```bash
   python -c "import fayde; print(fayde.__version__)"
   ```

   The output is the version number, for example `0.3.0`.

5. Check that you can run the `fayde` command:

   ```bash
   fayde --version
   ```

   The output is similar to the following:

   ```text
   fayde 0.3.0 (0123456789ab)
   ```

   The command is available while the virtual environment is active.


To solve a frame from Python, see [Use Fayde from Python](python.md).
