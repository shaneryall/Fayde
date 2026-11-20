# Security policy

## Supported versions

Security fixes go into the latest release. Fayde is at 0.x, so update to the latest version to
get a fix.

## Report a vulnerability

Report a vulnerability privately, not in a public issue:

1. Open Fayde's [private reporting form](https://github.com/shaneryall/Fayde/security/advisories/new)
   on GitHub.
2. Describe the problem, the Fayde version (`fayde --version`), and how to reproduce it. If a
   file triggers it, attach the file or say how to make one.

Only you and the maintainers see the report. The maintainers aim to answer within seven days.
When a fix is released, the release notes credit you, unless you ask not to be named.

## Scope

Fayde reads files that may come from other people: FITS images and tables, source lists, WCS
files, instrument profiles, star catalogs, and pattern indexes. Any of
these that makes Fayde crash, hang, use unbounded memory or time, or write outside the output
paths you gave it is a vulnerability.

A wrong or missing solution for a valid image is a bug, not a vulnerability. Report it in a
public issue.
