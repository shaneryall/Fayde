use std::path::PathBuf;
use std::process::Command;

fn scratch(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("fayde-cli-o-{tag}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn fayde(args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_fayde"))
        .args(args)
        .output()
        .unwrap()
}

#[test]
fn an_error_before_the_solve_removes_a_stale_output() {
    let d = scratch("stale");
    let out = d.join("frame.wcs");
    std::fs::write(&out, b"an earlier run's WCS").unwrap();
    let r = fayde(&[
        "solve",
        "--xyls",
        d.join("missing.xyls").to_str().unwrap(),
        "--ucac4",
        d.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ]);
    assert!(!r.status.success());
    assert!(
        !out.exists(),
        "a stale output must not survive a failed run"
    );
    assert!(String::from_utf8_lossy(&r.stderr).contains("removed the existing"));
    let _ = std::fs::remove_dir_all(&d);
}

#[test]
fn the_output_cannot_be_an_input() {
    let d = scratch("input");
    let wcs = d.join("given.wcs");
    std::fs::write(&wcs, b"the caller's WCS").unwrap();
    let r = fayde(&[
        "refit",
        "--xyls",
        d.join("missing.xyls").to_str().unwrap(),
        "--ucac4",
        d.to_str().unwrap(),
        "--wcs",
        wcs.to_str().unwrap(),
        "-o",
        wcs.to_str().unwrap(),
    ]);
    assert!(!r.status.success());
    assert!(String::from_utf8_lossy(&r.stderr).contains("is the --wcs input"));
    assert_eq!(std::fs::read(&wcs).unwrap(), b"the caller's WCS");
    let _ = std::fs::remove_dir_all(&d);
}

#[test]
fn detect_cannot_overwrite_its_instrument_profile() {
    let d = scratch("detect-profile");
    let profile = d.join("profile.json");
    std::fs::write(&profile, "{}").unwrap();
    let out = d.join("./profile.json");
    let r = fayde(&[
        "detect",
        "--image",
        "missing.fits",
        "--instrument",
        profile.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ]);
    assert!(!r.status.success());
    assert!(String::from_utf8_lossy(&r.stderr).contains("is the --instrument input"));
    assert_eq!(std::fs::read_to_string(profile).unwrap(), "{}");
    let _ = std::fs::remove_dir_all(&d);
}

#[test]
fn calibrate_cannot_overwrite_its_list_or_the_files_it_names() {
    let d = scratch("calibrate-inputs");
    let frame = d.join("frame.fits");
    let wcs = d.join("frame.wcs");
    let list = d.join("solved.tsv");
    let list_text = format!("{}\t{}\n", frame.display(), wcs.display());
    std::fs::write(&frame, "frame").unwrap();
    std::fs::write(&wcs, "wcs").unwrap();
    std::fs::write(&list, &list_text).unwrap();
    for (path, input, contents) in [
        (&list, "the --solved input", list_text.as_str()),
        (&frame, "a frame in the --solved list", "frame"),
        (&wcs, "a WCS in the --solved list", "wcs"),
    ] {
        let r = fayde(&[
            "calibrate",
            "--solved",
            list.to_str().unwrap(),
            "-o",
            path.to_str().unwrap(),
        ]);
        assert!(!r.status.success());
        let stderr = String::from_utf8_lossy(&r.stderr);
        assert!(stderr.contains(input), "{stderr}");
        assert_eq!(std::fs::read_to_string(path).unwrap(), contents);
    }
    let _ = std::fs::remove_dir_all(&d);
}

#[test]
fn verify_never_touches_the_path_given_as_output() {
    let d = scratch("verify");
    let keep = d.join("notes.txt");
    std::fs::write(&keep, b"not a WCS of this run").unwrap();
    let r = fayde(&[
        "verify",
        "--xyls",
        d.join("missing.xyls").to_str().unwrap(),
        "--ucac4",
        d.to_str().unwrap(),
        "--wcs",
        d.join("missing.wcs").to_str().unwrap(),
        "-o",
        keep.to_str().unwrap(),
    ]);
    assert!(!r.status.success());
    let stderr = String::from_utf8_lossy(&r.stderr);
    assert!(
        stderr.contains("-o is not used by fayde verify"),
        "{stderr}"
    );
    assert_eq!(std::fs::read(&keep).unwrap(), b"not a WCS of this run");
    let _ = std::fs::remove_dir_all(&d);
}

#[test]
fn a_command_line_mistake_keeps_the_output() {
    let d = scratch("cmdline");
    let out = d.join("frame.wcs");
    std::fs::write(&out, b"an earlier run's WCS").unwrap();
    let r = fayde(&[
        "solve",
        "--ucac4",
        d.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ]);
    assert!(!r.status.success());
    assert!(String::from_utf8_lossy(&r.stderr).contains("is required"));
    assert!(
        out.exists(),
        "a command-line error must not remove the output"
    );
    let _ = std::fs::remove_dir_all(&d);
}
