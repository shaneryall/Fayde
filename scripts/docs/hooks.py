"""MkDocs hooks for Fayde's documentation site.

- Links from a page to a file outside `docs/` (a source file, a profile, the changelog) point to
  that file on GitHub, at the ref in the `DOCS_REF` environment variable (default `main`).
- The changelog page and the command-line reference are generated at build time, so neither is
  written twice. The command-line reference runs the `fayde` binary named by `FAYDE_BIN`
  (default `fayde` on the path).
"""

from __future__ import annotations

import os
import posixpath
import re
import subprocess

from mkdocs.structure.files import File

REPO = "https://github.com/shaneryall/Fayde"
LINK = re.compile(r"(\]\()([^)\s]+)(\))")
COMMANDS = [
    [], ["solve"], ["detect"], ["verify"], ["refit"], ["calibrate"],
    ["index", "build"], ["index", "info"],
    ["catalog", "download"], ["catalog", "check"], ["catalog", "import"],
    ["compare"], ["bench"],
]


def _ref() -> str:
    return os.environ.get("DOCS_REF", "main")


def _outside(page_src: str, target: str, root_relative: bool) -> str | None:
    """The repository path of `target` when it lies outside docs/, else None."""
    if re.match(r"^[a-z]+:|^#|^/", target):
        return None
    path, _, frag = target.partition("#")
    base = "" if root_relative else posixpath.join("docs", posixpath.dirname(page_src))
    repo_path = posixpath.normpath(posixpath.join(base, path))
    if repo_path.startswith("docs/") and not root_relative:
        return None
    if root_relative and repo_path.startswith("docs/") and repo_path.endswith(".md"):
        return None
    return repo_path + (f"#{frag}" if frag else "")


def on_page_markdown(markdown: str, page, config, files) -> str:
    src = page.file.src_uri
    root_relative = src == "reference/changelog.md"

    def fix(m: re.Match) -> str:
        target = m.group(2)
        repo_path = _outside(src, target, root_relative)
        if repo_path is None:
            if root_relative and target.startswith("docs/"):
                return f"{m.group(1)}../{target[len('docs/'):]}{m.group(3)}"
            return m.group(0)
        kind = "tree" if repo_path.endswith("/") else "blob"
        return f"{m.group(1)}{REPO}/{kind}/{_ref()}/{repo_path}{m.group(3)}"

    return LINK.sub(fix, markdown)


def _cli_reference() -> str:
    fayde = os.environ.get("FAYDE_BIN", "fayde")
    version = subprocess.run([fayde, "--version"], capture_output=True, text=True, check=True)
    out = [
        "# Command-line reference",
        "",
        "Each command's options and defaults, as printed by its `--help` option in "
        f"`{version.stdout.strip()}`. For tasks and examples, see "
        "[Use the command line](../usage.md).",
        "",
    ]
    for args in COMMANDS:
        run = subprocess.run([fayde, *args, "--help"], capture_output=True, text=True)
        if run.returncode != 0:
            continue  # a command this version doesn't have
        text = (run.stdout or run.stderr).rstrip()
        out += [f"## fayde {' '.join(args)}".rstrip(), "", "```text", text, "```", ""]
    return "\n".join(out)


def on_files(files, config):
    root = os.path.dirname(config.config_file_path)
    with open(os.path.join(root, "CHANGELOG.md"), encoding="utf-8") as f:
        changelog = f.read()
    files.append(File.generated(config, "reference/changelog.md", content=changelog))
    files.append(File.generated(config, "reference/cli.md", content=_cli_reference()))
    return files
