"""Publish the Fayde crates to crates.io, one at a time, in dependency order.

    python3 scripts/publish_crates.py           # publish every crate of the workspace version
    python3 scripts/publish_crates.py --plan    # print the order and what is already published

crates.io lets one account publish 5 brand-new crates at once and then 1 every 10 minutes
(https://crates.io/docs/rate-limits), and Fayde has more crates than that. `cargo publish
--workspace` stops at the first refusal and cannot be started again once some crates are up.
So this script

- skips a crate whose version is already on crates.io, which makes a run that failed half way
  safe to start again;
- waits when crates.io answers "429 Too Many Requests" and then tries the same crate again.

The first release therefore takes about 40 minutes. Later releases publish new versions of
existing crates, which crates.io allows 30 at once, and take a few minutes.

Needs `cargo` and, to publish, a token in CARGO_REGISTRY_TOKEN. Uses only the standard library.
"""

import datetime
import email.utils
import json
import os
import pathlib
import re
import subprocess
import sys
import time
import urllib.error
import urllib.request

ROOT = pathlib.Path(__file__).resolve().parent.parent
API = os.environ.get("FAYDE_CRATES_IO_API", "https://crates.io/api/v1")
USER_AGENT = "fayde-release-script (https://github.com/shaneryall/Fayde)"
DEFAULT_WAIT_S = 600
MAX_ATTEMPTS = 12


def publishable_in_order():
    """(name, version) of every workspace crate that goes to crates.io, dependencies first."""
    out = subprocess.run(
        ["cargo", "metadata", "--format-version", "1", "--locked"],
        cwd=ROOT, check=True, capture_output=True, text=True,
    ).stdout
    meta = json.loads(out)
    members = set(meta["workspace_members"])
    crates = {
        p["name"]: p for p in meta["packages"]
        if p["id"] in members and (p["publish"] is None or "crates-io" in p["publish"])
    }
    needs = {
        name: sorted({d["name"] for d in p["dependencies"] if d["name"] in crates and d["name"] != name})
        for name, p in crates.items()
    }
    order = []
    while len(order) < len(crates):
        ready = sorted(n for n in crates if n not in order and all(d in order for d in needs[n]))
        if not ready:
            raise SystemExit("the workspace crates depend on each other in a cycle")
        order.append(ready[0])
    return [(name, crates[name]["version"]) for name in order]


def is_published(name, version):
    """Whether crates.io already has this version of this crate."""
    request = urllib.request.Request(f"{API}/crates/{name}/{version}", headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            return response.status == 200
    except urllib.error.HTTPError as error:
        if error.code == 404:
            return False
        raise SystemExit(f"crates.io answered {error.code} when asked about {name} {version}")


def wait_seconds(output):
    """How long crates.io asks us to wait, from "Please try again after <date>"."""
    found = re.search(r"try again after ([A-Za-z]{3}, [^\n]*?GMT)", output)
    if found:
        try:
            at = email.utils.parsedate_to_datetime(found.group(1))
            now = datetime.datetime.now(datetime.timezone.utc)
            return max(5.0, (at - now).total_seconds() + 5.0)
        except (TypeError, ValueError):
            pass
    return float(DEFAULT_WAIT_S)


def publish(name):
    """Run `cargo publish` for one crate. Returns its exit status and everything it printed."""
    process = subprocess.Popen(
        ["cargo", "publish", "-p", name, "--locked"],
        cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
    )
    lines = []
    for line in process.stdout:
        sys.stdout.write(line)
        sys.stdout.flush()
        lines.append(line)
    return process.wait(), "".join(lines)


def main():
    plan_only = "--plan" in sys.argv[1:]
    crates = publishable_in_order()
    print("order: " + ", ".join(name for name, _ in crates))
    done = []
    for name, version in crates:
        if is_published(name, version):
            print(f"{name} {version}: already on crates.io, skipped")
            continue
        if plan_only:
            print(f"{name} {version}: to publish")
            continue
        for attempt in range(1, MAX_ATTEMPTS + 1):
            status, output = publish(name)
            if status == 0:
                done.append(name)
                break
            if "already exists" in output or "already uploaded" in output:
                print(f"{name} {version}: already on crates.io, skipped")
                break
            if "429" not in output and "too many" not in output.lower():
                raise SystemExit(f"{name} {version}: cargo publish failed (see above)")
            if attempt == MAX_ATTEMPTS:
                raise SystemExit(f"{name} {version}: still refused after {MAX_ATTEMPTS} attempts")
            wait = wait_seconds(output)
            print(f"{name} {version}: crates.io asks to wait; next attempt in {wait:.0f} s", flush=True)
            time.sleep(wait)
    if not plan_only:
        print(f"published {len(done)} crate(s): " + (", ".join(done) or "none"))


if __name__ == "__main__":
    main()
