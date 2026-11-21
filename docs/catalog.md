# Get the UCAC4 catalog

Fayde reads the UCAC4 catalog (Zacharias et al. 2013) from a local copy of its original binary
files. Fayde doesn't include the catalog. The Strasbourg astronomical Data Center (CDS)
publishes it as [catalog I/322A](https://vizier.cfa.harvard.edu/ftp/cats/I/322A/).
The download is about 8.9 GB.

Use a terminal on Linux or macOS, or PowerShell on Windows. On Windows, use `curl.exe` in place
of `curl` in the download commands.

To download the catalog, follow these steps:

1. On Linux or macOS, create the catalog directories:

   ```bash
   mkdir -p UCAC4/u4b UCAC4/u4i
   ```

   On Windows, create the directories:

   ```powershell
   mkdir UCAC4\u4b, UCAC4\u4i
   ```

2. Go to the `u4b` directory:

   ```bash
   cd UCAC4/u4b
   ```

3. Download the 900 zone files, `z001` to `z900`:

   ```bash
   curl -fO "https://vizier.cfa.harvard.edu/ftp/cats/I/322A/UCAC4/u4b/z[001-900]"
   ```

   If the download stops, start it again from the file that was downloading. For example, if it
   stopped at `z412`, change the range to `z[412-900]`.

4. Go to the `u4i` directory:

   ```bash
   cd ../u4i
   ```

5. Download the two index files:

   ```bash
   curl -fO "https://vizier.cfa.harvard.edu/ftp/cats/I/322A/UCAC4/u4i/{u4index.unf,u4hpm.dat}"
   ```

   Fayde uses `u4index.unf` to seek directly to the stars that it needs and `u4hpm.dat` for the
   proper motions of fast-moving stars. If `u4index.unf` is missing, Fayde finds the same stars
   by searching each zone file.

6. Return to the directory where you started:

   ```bash
   cd ../..
   ```

7. On Linux or macOS, count the zone files:

   ```bash
   ls UCAC4/u4b | wc -l
   ```

   On Windows, count the zone files:

   ```powershell
   (Get-ChildItem UCAC4\u4b).Count
   ```

   Check that the output is `900`. If it's lower, download the missing files again.

The catalog directory has this layout:

```text
UCAC4/
├── u4b/
│   ├── z001
│   ├── ...
│   └── z900
└── u4i/
    ├── u4hpm.dat
    └── u4index.unf
```

When a Fayde command has a `--ucac4` option, give it the path to the `UCAC4` directory, which
is the directory that contains `u4b`. Download all 900 zone files, because Fayde reports an
error when a zone that it needs is missing.

To leave out `--ucac4`, set the `FAYDE_UCAC4` environment variable to the path to the `UCAC4`
directory. For more information about this variable, see
[Set the catalog directory once](usage.md#set-the-catalog-directory-once).
