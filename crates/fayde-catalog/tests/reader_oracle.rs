//! Cross-check against an independent UCAC4 reader on real zone files.
//!
//! Needs local UCAC4 data and an oracle file, neither of which is in the repository:
//! - `$FAYDE_UCAC4_DIR` (default `<workspace>/data/ucac4`) with `u4b/` and `u4i/`;
//! - `ucac4_oracle.json` in that directory: cone searches made with another UCAC4 reader.
//!
//! Oracle format: `{"cones": [{"name", "ra_deg", "dec_deg", "radius_arcsec", "epoch", "stars":
//! [{"id", "sep_arcsec", "ra_deg", "dec_deg", "mag"}]}]}`. `epoch` is a Julian year, or null for
//! catalogue-epoch (J2000) positions; `id` is the UCAC4 designation `UCAC4-zzz-nnnnnn`; `mag` is
//! the model magnitude, or the aperture magnitude where the model one is the 20.000 placeholder.
//!
//! When the data is absent the tests print a notice and pass without checking anything.

use std::path::PathBuf;

use fayde_catalog::{
    ucac4::{Ucac4Catalog, Ucac4Hit},
    StarQuery,
};
use fayde_core::{arcsec, to_arcsec, Epoch, SkyPos};
use serde::Deserialize;

#[derive(Deserialize)]
struct Oracle {
    cones: Vec<Cone>,
}

#[derive(Deserialize)]
struct Cone {
    name: String,
    ra_deg: f64,
    dec_deg: f64,
    radius_arcsec: f64,
    epoch: Option<f64>,
    stars: Vec<OracleStar>,
}

#[derive(Deserialize)]
struct OracleStar {
    id: String,
    sep_arcsec: f64,
    ra_deg: f64,
    dec_deg: f64,
    mag: f64,
}

fn data_dir() -> Option<PathBuf> {
    let dir = std::env::var_os("FAYDE_UCAC4_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data/ucac4"));
    let ok = dir.join("u4b").is_dir() && dir.join("ucac4_oracle.json").is_file();
    assert!(ok, "no UCAC4 data/oracle at {}", dir.display());
    Some(dir)
}

fn designation(h: &Ucac4Hit) -> String {
    format!("UCAC4-{:03}-{:06}", h.zone, h.record)
}

/// The oracle's reader measures separation with a flat-sky approximation; stars within this
/// margin of the cone edge may legitimately fall on different sides.
const EDGE_MARGIN_ARCSEC: f64 = 1.0;
/// Both readers propagate the same proper motions; they differ only at second order.
const POS_TOL_MAS: f64 = 1.0;

#[test]
#[ignore = "needs UCAC4 data: set FAYDE_UCAC4_DIR, run with --include-ignored"]
fn matches_independent_reader() {
    let Some(dir) = data_dir() else { return };
    let oracle: Oracle =
        serde_json::from_str(&std::fs::read_to_string(dir.join("ucac4_oracle.json")).unwrap())
            .unwrap();
    let cat = Ucac4Catalog::open(&dir).unwrap();
    assert!(cat.has_index());
    assert!(cat.high_pm_count() > 0);

    let mut compared = 0usize;
    let mut max_d_mas = 0.0f64;
    let mut high_pm_seen = 0usize;
    for cone in &oracle.cones {
        let centre = SkyPos::from_degrees(cone.ra_deg, cone.dec_deg).unwrap();
        let query = StarQuery {
            epoch: cone.epoch.and_then(Epoch::julian_year),
            ..StarQuery::default()
        };
        let hits = cat
            .cone_hits(centre, arcsec(cone.radius_arcsec), &query)
            .unwrap_or_else(|e| panic!("{}: {e}", cone.name));
        let inner = cone.radius_arcsec - EDGE_MARGIN_ARCSEC;

        let ours: std::collections::HashMap<String, &Ucac4Hit> =
            hits.iter().map(|h| (designation(h), h)).collect();
        for s in cone.stars.iter().filter(|s| s.sep_arcsec < inner) {
            let h = ours
                .get(&s.id)
                .unwrap_or_else(|| panic!("{}: {} missing from Fayde", cone.name, s.id));
            if h.raw.is_high_pm() {
                // The oracle's reader does not load u4hpm.dat and keeps these at J2000; Fayde
                // moves them.
                high_pm_seen += 1;
                let j2000 = h.raw.position_j2000();
                let theirs = SkyPos::from_degrees(s.ra_deg, s.dec_deg).unwrap();
                assert!(to_arcsec(j2000.separation(theirs)) * 1e3 < POS_TOL_MAS);
                if cone.epoch.is_some() {
                    assert!(to_arcsec(h.star.pos.separation(j2000)) > 10.0, "{}", s.id);
                }
                continue;
            }
            let theirs = SkyPos::from_degrees(s.ra_deg, s.dec_deg).unwrap();
            let d_mas = to_arcsec(h.star.pos.separation(theirs)) * 1e3;
            assert!(d_mas < POS_TOL_MAS, "{} {}: {d_mas} mas", cone.name, s.id);
            max_d_mas = max_d_mas.max(d_mas);
            if let Some(m) = h.star.mag {
                assert!(
                    (f64::from(m) - s.mag).abs() < 1e-4,
                    "{} {} mag",
                    cone.name,
                    s.id
                );
            }
            compared += 1;
        }
        let theirs: std::collections::HashSet<&str> =
            cone.stars.iter().map(|s| s.id.as_str()).collect();
        for h in &hits {
            let sep = to_arcsec(centre.separation(h.star.pos));
            if sep < inner && !h.raw.is_high_pm() {
                assert!(
                    theirs.contains(designation(h).as_str()),
                    "{}: {} found by Fayde only",
                    cone.name,
                    designation(h)
                );
            }
        }
        eprintln!(
            "{}: {} stars (oracle {})",
            cone.name,
            hits.len(),
            cone.stars.len()
        );
    }
    eprintln!("compared {compared} stars, max position difference {max_d_mas:.2e} mas, {high_pm_seen} high-PM");
    assert!(compared > 500, "only {compared} stars compared");
    assert!(
        high_pm_seen > 0,
        "oracle should include a high-proper-motion star"
    );
}

#[test]
#[ignore = "needs UCAC4 data: set FAYDE_UCAC4_DIR, run with --include-ignored"]
fn index_seek_equals_full_scan() {
    let Some(dir) = data_dir() else { return };
    let cat = Ucac4Catalog::open(&dir).unwrap();
    for (ra, dec, r) in [
        (150.0, 2.0, 1080.0),
        (0.05, 2.1, 720.0),
        (359.9, 1.9, 300.0),
    ] {
        let centre = SkyPos::from_degrees(ra, dec).unwrap();
        let q = StarQuery {
            epoch: Epoch::julian_year(2026.0),
            ..StarQuery::default()
        };
        let a = cat.cone_hits(centre, arcsec(r), &q).unwrap();
        let b = cat.cone_hits_full_scan(centre, arcsec(r), &q).unwrap();
        assert_eq!(a, b, "({ra}, {dec}, {r}\")");
        assert!(!a.is_empty());
    }
}

#[test]
#[ignore = "needs UCAC4 data: set FAYDE_UCAC4_DIR, run with --include-ignored"]
fn strict_query_reports_missing_zone() {
    let Some(dir) = data_dir() else { return };
    let cat = Ucac4Catalog::open(&dir).unwrap();
    let centre = SkyPos::from_degrees(10.0, -45.0).unwrap(); // zone 226: not copied locally
    if cat.zone_path(226).exists() {
        return;
    }
    assert!(cat
        .cone(centre, arcsec(60.0), &StarQuery::default())
        .is_err());
    let lax = StarQuery {
        strict: false,
        ..StarQuery::default()
    };
    assert!(cat.cone(centre, arcsec(60.0), &lax).unwrap().is_empty());
}

use fayde_catalog::Catalog;
