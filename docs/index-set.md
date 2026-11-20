# Pattern indexes for sub-degree fields

A pattern index lets Fayde solve a frame when the pointing or pixel scale is uncertain or
unknown. If you know the pointing and have a precise pixel scale from a frame header or
instrument profile, you need only the UCAC4 catalog.

An index covers fields near its target field of view. Fayde chooses among the indexes you
provide. To cover a wider range of scales, follow these steps:

1. Build an index for each field of view.
2. Pass each index to `fayde solve` with a separate `--index` option.

## Build the recommended indexes

The two recommended indexes cover the whole sky at epoch 2025.0. You need the full UCAC4 catalog,
with all 900 zone files, and at least 32 GB of free memory.

Replace `UCAC4_DIR` with the path to the directory that contains `u4b`. To cover both field
sizes, build both indexes:

1. Build the index for fields of about 1°:

   ```bash
   fayde index build --ucac4 UCAC4_DIR --fov 1.0 --stars-per-fov 23.7 --pattern-stars-per-fov 23.7 --mag-limit 15.5 --min-edge 0.04 --epoch 2025 -o allsky-1.0.fidx
   ```

2. Build the index for fields somewhat wider than 1°:

   ```bash
   fayde index build --ucac4 UCAC4_DIR --fov 1.3 --stars-per-fov 40 --pattern-stars-per-fov 20 --stars-per-field 12 --mag-limit 15.5 --min-edge 0.05 --epoch 2025 -o allsky-1.3.fidx
   ```

The builds used the version 0.1.0 binary on a 32-core x86_64 Linux computer. Other jobs ran on
the computer during the builds, with a load of about 38. The table lists the build results:

| Index | Stars | Patterns | File size | Build time | Peak memory |
|---|---:|---:|---:|---:|---:|
| `allsky-1.0.fidx` | 1,572,824 | 259,900,180 | 17,211,606,974 bytes | 30 min 26 s | 22.7 GB |
| `allsky-1.3.fidx` | 1,572,824 | 153,951,343 | 10,218,984,414 bytes | 21 min 50 s | 13.6 GB |

## Plan for memory

Fayde loads the entire index into memory. On an otherwise idle 64-bit Arm computer,
`fayde index info` loaded each recommended index with these results:

| Index | Load time | Peak memory |
|---|---:|---:|
| `allsky-1.0.fidx` | 48 s | 18.3 GB |
| `allsky-1.3.fidx` | 27 s | 10.9 GB |

The `fayde` command loads the index each time that it runs. To solve many frames, use the
Python package and reuse one `fayde.Solver`, which loads the index once.

## Check a content hash

Each index has a content hash, a 64-bit checksum that covers its stars, patterns, and build
options. Fayde uses the hash to identify the index that produced a solution and detect file
damage. The hash isn't a cryptographic signature. Fayde records each index's hash in the
`identity.indexes[].content_hash` field of the solve report.

To compare your index with the expected hash, follow these steps:

1. Print the index metadata:

   ```bash
   fayde index info allsky-1.0.fidx
   ```

2. Compare the `content_hash` field in the output with the following table.

An index's bytes and hash don't depend on the platform or build profile. Fayde uses
`fayde_core::Portable` for all its math. Builds from the full UCAC4 catalog with the commands
on this page produce these hashes:

| Index | `content_hash` |
|---|---|
| `allsky-1.0.fidx` | `1ddf14fa50ad59dd` |
| `allsky-1.3.fidx` | `4ad74722223054e1` |

Version 0.1.0 and every later version build the same bytes. If a release changes the bytes,
the release notes describe the change, and this page lists the new hashes.

Continuous integration doesn't build the recommended indexes, because each one takes tens of
gigabytes of memory and about half an hour.

## Index build options

To list every option and its default value, run `fayde index build --help`. The table describes
the options used in the recommended commands:

| Option | Meaning |
|---|---|
| `--fov` | The field of view that the index is built for, in degrees |
| `--stars-per-fov` | The number of stars kept per field of view, selected uniformly over HEALPix cells |
| `--pattern-stars-per-fov` | The number of those stars, taken brightest first, that form patterns |
| `--stars-per-field` | The number of stars per lattice field from which patterns are formed |
| `--mag-limit` | The faintest magnitude considered |
| `--min-edge` | The shortest pattern edge, in degrees. Below it, centroid errors dominate |
| `--epoch` | The epoch that the star positions are propagated to, as a Julian year |
| `--region RA,DEC,RADIUS` | Builds a regional index, in degrees, instead of an all-sky index. A regional build takes minutes, which makes it useful for tests |

For more information about the algorithms, see the [references](references.md).
