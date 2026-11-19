use std::{path::Path, process::Command};

fn checkout_commit(root: &Path) -> Option<String> {
    if !root.join(".git").exists() {
        return None;
    }
    let git = |args: &[&str]| {
        let output = Command::new("git")
            .current_dir(root)
            .arg("--no-optional-locks")
            .args(args)
            .output()
            .ok()?;
        output
            .status
            .success()
            .then(|| String::from_utf8_lossy(&output.stdout).into_owned())
    };
    let head = git(&["rev-parse", "--short=12", "HEAD"])?;
    let watch_git = |name: &str| {
        if let Some(path) = git(&["rev-parse", "--git-path", name]) {
            let path = root.join(path.trim_end());
            if let Some(existing) = path.ancestors().find(|p| p.exists()) {
                println!("cargo:rerun-if-changed={}", existing.display());
            }
        }
    };
    watch_git("HEAD");
    watch_git("index");
    if let Some(branch) = git(&["symbolic-ref", "-q", "HEAD"]) {
        watch_git(branch.trim_end());
        if git(&["rev-parse", "--git-path", "packed-refs"])
            .is_some_and(|p| root.join(p.trim_end()).is_file())
        {
            watch_git("packed-refs");
        }
    }
    for file in git(&["ls-files", "-z"])?.split_terminator('\0') {
        if !file.contains(['\r', '\n']) {
            println!("cargo:rerun-if-changed={}", root.join(file).display());
        }
    }
    let dirty = !git(&["status", "--porcelain", "--untracked-files=no"])?.is_empty();
    let head = head.trim_end();
    Some(if dirty {
        format!("{head}-dirty")
    } else {
        head.to_owned()
    })
}

fn main() {
    println!("cargo:rerun-if-env-changed=FAYDE_COMMIT");
    let root = Path::new("../..");
    let commit = std::env::var("FAYDE_COMMIT")
        .ok()
        .or_else(|| checkout_commit(root))
        .or_else(|| {
            let path = root.join(".git_archival.txt");
            let text = std::fs::read_to_string(&path).ok()?;
            println!("cargo:rerun-if-changed={}", path.display());
            let c = text.lines().find_map(|l| l.strip_prefix("commit:"))?.trim();
            (!c.is_empty() && !c.contains("$Format")).then(|| c.to_string())
        });
    println!(
        "cargo:rustc-env=FAYDE_COMMIT={}",
        commit.as_deref().unwrap_or("unknown")
    );
}
