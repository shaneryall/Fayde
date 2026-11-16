//! Records the source commit for `fayde --version` and the identity block of JSON reports.
//!
//! In order: `FAYDE_COMMIT` from the build environment; git, when building in a checkout
//! (`-dirty` when tracked files changed); `.git_archival.txt`, which `git archive` (and so a
//! forge's source downloads) fills in through `export-subst` (`.gitattributes`); else
//! "unknown".

use std::process::Command;

fn main() {
    println!("cargo:rerun-if-env-changed=FAYDE_COMMIT");
    println!("cargo:rerun-if-changed=../../.git/HEAD");
    println!("cargo:rerun-if-changed=../../.git/index");
    println!("cargo:rerun-if-changed=../../.git_archival.txt");
    let commit = std::env::var("FAYDE_COMMIT").ok().or_else(|| {
        let git = |args: &[&str]| {
            Command::new("git")
                .args(args)
                .output()
                .ok()
                .filter(|o| o.status.success())
                .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
        };
        let head = git(&["rev-parse", "--short=12", "HEAD"])?;
        let dirty =
            git(&["status", "--porcelain", "--untracked-files=no"]).is_some_and(|s| !s.is_empty());
        Some(if dirty { format!("{head}-dirty") } else { head })
    });
    // A `git archive` export: the placeholder was replaced by the short hash.
    let commit = commit.or_else(|| {
        let manifest = std::env::var("CARGO_MANIFEST_DIR").ok()?;
        let text = std::fs::read_to_string(
            std::path::Path::new(&manifest).join("../../.git_archival.txt"),
        )
        .ok()?;
        let c = text.lines().find_map(|l| l.strip_prefix("commit:"))?.trim();
        (!c.is_empty() && !c.contains("$Format")).then(|| c.to_string())
    });
    println!(
        "cargo:rustc-env=FAYDE_COMMIT={}",
        commit.unwrap_or_else(|| "unknown".into())
    );
}
