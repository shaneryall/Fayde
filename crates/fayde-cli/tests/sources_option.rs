use std::path::{Path, PathBuf};

use fayde_core::Portable;
use std::process::Command;

fn fayde(args: &[&str]) -> (std::process::ExitStatus, String, String) {
    let o = Command::new(env!("CARGO_BIN_EXE_fayde"))
        .args(args)
        .output()
        .unwrap();
    (
        o.status,
        String::from_utf8_lossy(&o.stdout).into(),
        String::from_utf8_lossy(&o.stderr).into(),
    )
}

fn scratch(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("fayde-sources-{tag}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn image(path: &Path) {
    let n = 128usize;
    let mut v = vec![100f32; n * n];
    for y in 0..n {
        for x in 0..n {
            let (dx, dy) = (x as f32 - 60.3, y as f32 - 70.6);
            v[y * n + x] += 5000.0 * (-(dx * dx + dy * dy) / 3.2).pexp();
            v[y * n + x] += ((x * 7 + y * 13) % 11) as f32;
        }
    }
    let cards = [
        "SIMPLE  =                    T",
        "BITPIX  =                  -32",
        "NAXIS   =                    2",
        "NAXIS1  =                  128",
        "NAXIS2  =                  128",
        "END",
    ];
    let mut out = Vec::new();
    for c in cards {
        out.extend(format!("{c:<80}").bytes());
    }
    out.resize(out.len().div_ceil(2880) * 2880, b' ');
    for p in &v {
        out.extend(p.to_be_bytes());
    }
    out.resize(out.len().div_ceil(2880) * 2880, 0);
    std::fs::write(path, out).unwrap();
}

#[test]
fn source_lists_and_images_are_told_apart() {
    let d = scratch("kinds");
    let img = d.join("frame.fits");
    image(&img);
    let list = d.join("frame.sources.fits");
    let (ok, out, err) = fayde(&[
        "detect",
        "--image",
        img.to_str().unwrap(),
        "-o",
        list.to_str().unwrap(),
    ]);
    assert!(ok.success(), "{err}");
    assert!(out.contains("\"count\""));
    let ucac4 = d.to_str().unwrap();

    let (ok, _, err) = fayde(&[
        "solve",
        "--sources",
        img.to_str().unwrap(),
        "--ucac4",
        ucac4,
    ]);
    assert!(
        !ok.success() && err.contains("is an image, not a source list"),
        "{err}"
    );
    let (ok, _, err) = fayde(&["solve", "--image", list.to_str().unwrap(), "--ucac4", ucac4]);
    assert!(
        !ok.success() && err.contains("is a source list, not an image"),
        "{err}"
    );
    let (ok, _, err) = fayde(&["detect", "--image", list.to_str().unwrap()]);
    assert!(
        !ok.success() && err.contains("is a source list, not an image"),
        "{err}"
    );

    let (_, _, err) = fayde(&["solve", "--xyls", img.to_str().unwrap(), "--ucac4", ucac4]);
    assert!(err.contains("is an image, not a source list"), "{err}");
    let l = list.to_str().unwrap();
    let (ok, _, err) = fayde(&["solve", "--xyls", l, "--sources", l, "--ucac4", ucac4]);
    assert!(
        !ok.success() && err.contains("another name for --sources"),
        "{err}"
    );
    let _ = std::fs::remove_dir_all(&d);
}
