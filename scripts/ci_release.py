"""Create the GitLab release for $CI_COMMIT_TAG with a link per file in dist/.

Runs in the `release` CI job. The links point at this job's artifacts, which never expire.
Uses only the standard library and the job token.
"""

import json
import os
import pathlib
import urllib.request

env = os.environ
tag = env["CI_COMMIT_TAG"]
job_url = f"{env['CI_PROJECT_URL']}/-/jobs/{env['CI_JOB_ID']}/artifacts/raw/dist"
files = sorted(p.name for p in pathlib.Path("dist").iterdir())
if not files:
    raise SystemExit("dist/ is empty")
body = {
    "tag_name": tag,
    "name": f"Fayde {tag}",
    "description": (
        "Linux binaries (glibc 2.17+, x86_64 and aarch64) and Python wheels (abi3, CPython 3.10+"
        "; manylinux2014). Windows and macOS wheels are built by hand until those CI runners "
        "exist. Index set, build commands and content hashes: `docs/index-set.md`. Versioning "
        "rules: `docs/versioning.md`. Changes: `CHANGELOG.md`. Licence: GPL-3.0-or-later (`LICENSE`, `NOTICE`); "
        "please cite as in `CITATION.cff`."
    ),
    "assets": {"links": [{"name": f, "url": f"{job_url}/{f}", "link_type": "package"} for f in files]},
}
req = urllib.request.Request(
    f"{env['CI_API_V4_URL']}/projects/{env['CI_PROJECT_ID']}/releases",
    data=json.dumps(body).encode(),
    headers={"JOB-TOKEN": env["CI_JOB_TOKEN"], "Content-Type": "application/json"},
    method="POST",
)
with urllib.request.urlopen(req) as r:
    print(r.status, json.load(r).get("_links", {}).get("self"))
