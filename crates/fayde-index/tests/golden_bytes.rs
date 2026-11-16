//! The index a build produces, byte for byte.
//!
//! The expected values were re-captured when all maths moved to `fayde_core::Portable`
//! (pure-Rust libm), and are the same on macOS aarch64, Linux x86_64 and in debug and release
//! builds. Before that, platform maths libraries gave different bytes. They pin the content
//! hash, the file length and
//! an FNV-1a 64 of the whole file (trailing hash included), so any change to selection, pattern
//! order, table layout or serialization shows up here.
//!
//! The metadata block holds the build parameters and `fayde_index::BUILDER`, so an intended
//! change to either (a new `BuildParams` field, a new builder) changes these values: re-capture
//! them then, from a build known to be good. A release that does not touch the builder keeps
//! `BUILDER`, and so these values (crate version bumps do not change index bytes).

use fayde_core::Portable;
use std::path::PathBuf;

use fayde_catalog::{ucac4::Ucac4Catalog, Catalog, CatalogId, Star, StarQuery};
use fayde_core::{Epoch, SkyPos};
use fayde_index::{build, BuildParams, Index};

/// What is pinned: content hash, file length, FNV-1a 64 of every byte of the file.
#[derive(Debug, PartialEq, Eq)]
struct Fingerprint {
    content_hash: u64,
    file_len: u64,
    file_fnv: u64,
}

fn fnv(bytes: &[u8]) -> u64 {
    bytes.iter().fold(0xcbf2_9ce4_8422_2325, |h, &b| {
        (h ^ u64::from(b)).wrapping_mul(0x0100_0000_01b3)
    })
}

fn fingerprint(index: &Index, tag: &str) -> Fingerprint {
    let path = std::env::temp_dir().join(format!("fayde-golden-{tag}-{}.fidx", std::process::id()));
    index.write(&path).unwrap();
    let bytes = std::fs::read(&path).unwrap();
    // The file must read back as the same index, from the file and from its bytes.
    assert_eq!(&Index::read(&path).unwrap(), index);
    assert_eq!(&Index::from_bytes(&bytes).unwrap(), index);
    std::fs::remove_file(&path).unwrap();
    Fingerprint {
        content_hash: index.content_hash(),
        file_len: bytes.len() as u64,
        file_fnv: fnv(&bytes),
    }
}

fn synthetic_id() -> CatalogId {
    CatalogId {
        name: "synthetic",
        release: "0",
        location: String::new(),
    }
}

/// Deterministic stars in a cap (linear congruential generator; no external data).
fn synthetic_stars(n: u64, centre: (f64, f64), radius_deg: f64) -> Vec<Star> {
    let mut s = 0x5DEE_CE66_D1CE_4E5Bu64;
    let mut rnd = || {
        s = s
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (s >> 11) as f64 / (1u64 << 53) as f64
    };
    (0..n)
        .map(|i| {
            let r = radius_deg * rnd().sqrt();
            let phi = 2.0 * std::f64::consts::PI * rnd();
            let dec = centre.1 + r * phi.psin();
            let ra = centre.0 + r * phi.pcos() / dec.to_radians().pcos();
            Star {
                id: i + 1,
                pos: SkyPos::from_degrees(ra, dec).unwrap(),
                mag: Some(6.0 + 9.0 * rnd() as f32),
                pos_err_mas: [20.0, 20.0],
                pm_mas_yr: [Some((40.0 * rnd() - 20.0) as f32), Some(0.0)],
            }
        })
        .collect()
}

#[test]
fn synthetic_regional_index_bytes_are_pinned() {
    // About 1,700 lattice fields and 60,000 patterns: every stage of the builder is exercised,
    // including duplicate patterns across overlapping fields.
    let stars = synthetic_stars(6000, (45.0, 10.0), 2.2);
    let params = BuildParams {
        fov_deg: 1.0,
        min_edge_deg: 0.04,
        stars_per_fov: 40.0,
        region: Some([45.0, 10.0, 2.0]),
        ..BuildParams::default()
    };
    let index = build(&stars, synthetic_id(), &params).unwrap();
    assert!(index.meta().pattern_count > 20_000, "{:?}", index.meta());
    assert_eq!(
        fingerprint(&index, "regional"),
        Fingerprint {
            content_hash: 3_196_662_882_412_501_652,
            file_len: 5_677_375,
            file_fnv: 9_766_863_859_479_986_696,
        }
    );
}

#[test]
fn synthetic_all_sky_index_bytes_are_pinned() {
    // No region: the whole lattice is walked (thinly, to keep the test fast).
    let stars = synthetic_stars(3000, (200.0, -35.0), 6.0);
    let params = BuildParams {
        fov_deg: 2.0,
        field_oversampling: 3.0,
        patterns_per_field: 20,
        ..BuildParams::default()
    };
    let index = build(&stars, synthetic_id(), &params).unwrap();
    assert!(index.meta().pattern_count > 100, "{:?}", index.meta());
    assert_eq!(
        fingerprint(&index, "allsky"),
        Fingerprint {
            content_hash: 12_283_015_016_262_690_935,
            file_len: 196_673,
            file_fnv: 18_274_905_964_828_468_422,
        }
    );
}

#[test]
#[ignore = "needs UCAC4 data: set FAYDE_UCAC4_DIR, run with --include-ignored"]
fn regional_ucac4_index_bytes_are_pinned() {
    // Same index as `regional_ucac4.rs`. Needs UCAC4 zones z458–z463; ignored by default (`cargo test -- --include-ignored`).
    let dir = std::env::var_os("FAYDE_UCAC4_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data/ucac4"));
    if !(458..=463).all(|z| dir.join(format!("u4b/z{z}")).is_file()) {
        panic!("UCAC4 zones z458–z463 not at {}", dir.display());
    }
    let cat = Ucac4Catalog::open(dir).unwrap();
    let params = BuildParams {
        fov_deg: 0.4,
        min_edge_deg: 0.02,
        stars_per_fov: 20.0,
        mag_limit: Some(15.0),
        epoch_year: 2025.0,
        region: Some([150.0, 2.0, 0.38]),
        ..BuildParams::default()
    };
    let query = StarQuery {
        epoch: Epoch::julian_year(params.epoch_year),
        ..StarQuery::default()
    };
    let centre = SkyPos::from_degrees(150.0, 2.0).unwrap();
    let stars = cat.cone(centre, 0.38f64.to_radians(), &query).unwrap();
    let index = build(&stars, cat.id(), &params).unwrap();
    assert_eq!(
        fingerprint(&index, "ucac4"),
        Fingerprint {
            content_hash: 16_340_545_770_047_993_857,
            file_len: 625_156,
            file_fnv: 13_819_180_066_941_641_570,
        }
    );
}
