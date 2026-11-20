# Versioning

<a id="versioning-and-releases"></a>Fayde follows semantic versioning. While the version is 0.x,
a breaking change increases the minor number, and anything else increases the patch number.

## Breaking changes

The following changes are breaking:

- Removing a command-line or Python option, or changing its meaning
- Removing or renaming a JSON report field
- Changing the index file format or the bytes of an index, including its `content_hash`

Adding a report field isn't a breaking change. If a release changes indexes, you must rebuild
your indexes as directed in the release notes.

A fix that leaves reports and indexes unchanged is a patch release. So is a new command,
option, or report field that changes nothing that already exists.

## Deprecations

Before Fayde removes or renames a command-line option, a Python argument, or a report field, one
release deprecates it:

1. The release that deprecates it still accepts the old name and behaves as before. Fayde prints
   a warning that names the replacement, and the release notes list the deprecation.
2. A later breaking release removes it, no earlier than the next minor version.

A change that can't keep the old behavior, such as a new index format, can't be deprecated
first. Its release notes say what to do, for example to rebuild your indexes.

## Releases

Each release has a section in the [changelog](https://github.com/shaneryall/Fayde/blob/main/CHANGELOG.md)
that lists its changes and marks the breaking ones. The same version is published on GitHub,
PyPI, and crates.io.
