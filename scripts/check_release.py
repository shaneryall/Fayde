"""Check that the files a release is made from agree with each other.

    python3 scripts/check_release.py                 # consistency (run in CI on every push)
    python3 scripts/check_release.py --tag v1.2.3    # also: the tag is this version, and the
                                                     # changelog has no "Unreleased" section
    python3 scripts/check_release.py --notes         # print the changelog section of this version

Checked: the workspace version against the exact versions the Fayde crates require of each
other (Cargo.toml), against CITATION.cff and against the changelog; and the copies of LICENSE
and NOTICE in each crate against the originals (a published crate or wheel can only carry files
from its own directory). Uses only the standard library (Python 3.11 or later).
"""

import pathlib
import re
import sys
import tomllib

ROOT = pathlib.Path(__file__).resolve().parent.parent


def changelog_section(version):
    """Lines under the `## <version> (date)` heading, or None when there is no such heading."""
    lines = (ROOT / "CHANGELOG.md").read_text(encoding="utf-8").splitlines()
    start = next((i for i, l in enumerate(lines) if re.match(rf"## {re.escape(version)}( |$)", l)), None)
    if start is None:
        return None
    end = next((i for i in range(start + 1, len(lines)) if lines[i].startswith("## ")), len(lines))
    return "\n".join(lines[start + 1 : end]).strip("\n")


def problems(version, tag):
    found = []
    cargo = tomllib.loads((ROOT / "Cargo.toml").read_text(encoding="utf-8"))
    pinned = {}
    for name, dep in cargo["workspace"]["dependencies"].items():
        if isinstance(dep, dict) and "path" in dep:
            pinned[name] = dep.get("version")
            if dep.get("version") != f"={version}":
                found.append(f"Cargo.toml: {name} requires {dep.get('version')!r}, expected '={version}'")
    crates = sorted(p for p in (ROOT / "crates").iterdir() if (p / "Cargo.toml").is_file())
    for crate in crates:
        manifest = tomllib.loads((crate / "Cargo.toml").read_text(encoding="utf-8"))
        # A Fayde crate named with its own `path` would be published without a version.
        for table in ("dependencies", "dev-dependencies", "build-dependencies"):
            for name, dep in manifest.get(table, {}).items():
                if isinstance(dep, dict) and "path" in dep:
                    found.append(f"{crate.name}/Cargo.toml: {name} has its own path; use `{name}.workspace = true`")
        for name in ("LICENSE", "NOTICE"):
            copy = crate / name
            if not copy.is_file() or copy.read_bytes() != (ROOT / name).read_bytes():
                found.append(f"crates/{crate.name}/{name} is missing or differs from {name}")
    cff = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    cited = re.search(r'^version:\s*"?([^"\s]+)"?\s*$', cff, re.M)
    if not cited or cited.group(1) != version:
        found.append(f"CITATION.cff: version is {cited.group(1) if cited else 'missing'}, expected {version}")
    if changelog_section(version) is None:
        found.append(f"CHANGELOG.md: no section '## {version} (date)'")
    if tag is not None:
        if tag != f"v{version}":
            found.append(f"tag {tag} does not match the workspace version {version} (expected v{version})")
        if changelog_section("Unreleased") is not None:
            found.append("CHANGELOG.md: move the entries under 'Unreleased' to the section of this version")
    return found


def main():
    args = sys.argv[1:]
    cargo = tomllib.loads((ROOT / "Cargo.toml").read_text(encoding="utf-8"))
    version = cargo["workspace"]["package"]["version"]
    if "--notes" in args:
        section = changelog_section(version)
        if section is None:
            raise SystemExit(f"CHANGELOG.md: no section '## {version} (date)'")
        print(section)
        return
    tag = None
    if "--tag" in args:
        i = args.index("--tag")
        if i + 1 >= len(args):
            raise SystemExit("--tag needs a value, for example --tag v1.2.3")
        tag = args[i + 1]
    found = problems(version, tag)
    if found:
        raise SystemExit("release check failed:\n  " + "\n  ".join(found))
    print(f"release files agree: version {version}" + (f", tag {tag}" if tag else ""))


if __name__ == "__main__":
    main()
