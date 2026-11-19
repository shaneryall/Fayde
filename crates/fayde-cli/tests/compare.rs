use std::process::Command;

#[test]
fn hinted_only_does_not_count_a_blind_failure() {
    let dir = std::env::temp_dir().join(format!("fayde-compare-{}", std::process::id()));
    let zone = dir.join("u4b");
    std::fs::create_dir_all(&zone).unwrap();
    std::fs::write(zone.join("z460"), []).unwrap();
    let path = dir.join("frames.json");
    let frames = serde_json::json!({"frames": [{
        "width": 128,
        "height": 128,
        "sources_fits": [],
        "solve_field_wcs_cards": [
            "CTYPE1  = 'RA---TAN'", "CTYPE2  = 'DEC--TAN'", "CRVAL1  = 150.0",
            "CRVAL2  = 1.9", "CRPIX1  = 64.0", "CRPIX2  = 64.0",
        ],
    }]});
    std::fs::write(&path, frames.to_string()).unwrap();
    let result = Command::new(env!("CARGO_BIN_EXE_fayde"))
        .args(["compare", "--hinted-only", "--ucac4"])
        .arg(&dir)
        .arg("--frames")
        .arg(&path)
        .output()
        .unwrap();
    let _ = std::fs::remove_dir_all(&dir);
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(result.status.success(), "{stderr}");
    assert!(stderr.contains("skipped (--hinted-only)"));
    let report: serde_json::Value = serde_json::from_slice(&result.stdout).unwrap();
    assert_eq!(report["solve_field_solved"], 1);
    assert_eq!(report["header_priors_failed_where_solve_field_solved"], 1);
    assert_eq!(report["blind_failed_where_solve_field_solved"], 0);
    assert!(report["blind_seconds_median"].is_null());
}
