//! Read files written by Astropy (`scripts/gen_fits_fixtures.py`) and round-trip WCS files.

use fayde_core::Portable;
use std::path::PathBuf;

use fayde_core::{Epoch, PixelPos, SkyPos};
use fayde_fits::{header_priors, read_sources, write_wcs_file, FitsFile};
use fayde_wcs::{TanWcs, Wcs};

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/fixtures/fits")
        .join(name)
}

#[test]
fn reads_a_source_list_and_converts_to_zero_based() {
    let x = read_sources(fixture("frame.sources.fits")).unwrap();
    assert_eq!((x.width, x.height), (Some(768), Some(768)));
    assert_eq!(x.sources.len(), 4);
    let s = x.sources[1];
    assert_eq!((s.x, s.y), (199.25, 299.75)); // 200.25, 300.75 in FITS 1-based
    assert_eq!(s.flux, Some(500.5));
    assert_eq!((x.sources[3].x, x.sources[3].y), (0.0, 767.0));
}

#[test]
fn header_priors_cd_sexagesimal_iso_date() {
    let f = FitsFile::open(fixture("frame.fits")).unwrap();
    let p = header_priors(&f.primary().header);
    let (ptg, origin) = p.pointing.unwrap();
    assert_eq!(origin, "CRVAL1/CRVAL2");
    assert!(ptg.separation(SkyPos::from_degrees(80.9212084123, 0.121863492816).unwrap()) < 1e-12);
    let (scale, origin) = p.scale_arcsec.unwrap();
    assert_eq!(origin, "header linear WCS");
    let det: f64 = 0.000799901108658f64.ppowi(2) + 0.000225f64.ppowi(2);
    assert!((scale - det.sqrt() * 3600.0).abs() < 1e-9);
    // DATE-OBS 2020-11-21T16:46:40.074 is the start; mid-exposure adds 3.001 s / 2.
    let (e, _) = p.epoch.unwrap();
    let want = Epoch::from_mjd(59_174.699_075_000_0 + 1.5005 / 86_400.0).unwrap();
    assert!(
        (e.year() - want.year()).abs() * 365.25 * 86_400.0 < 0.01,
        "{e:?} vs {want:?}"
    );
    assert_eq!(p.size, Some((6, 4)));
    assert!(p.warnings.is_empty());
}

#[test]
fn header_priors_pc_cdelt_and_mjd() {
    let f = FitsFile::open(fixture("frame_pc.fits")).unwrap();
    let p = header_priors(&f.primary().header);
    assert!((p.scale_arcsec.unwrap().0 - 2.88).abs() < 1e-9);
    assert!(
        (p.epoch.unwrap().0.year() - Epoch::from_mjd(59174.69907).unwrap().year()).abs() < 1e-12
    );
}

#[test]
fn mixed_cd_pc_gives_no_scale() {
    let f = FitsFile::open(fixture("frame_mixed.fits")).unwrap();
    let p = header_priors(&f.primary().header);
    assert!(p.scale_arcsec.is_none());
    assert_eq!(p.warnings.len(), 1);
    // Pointing is still available.
    assert!(p.pointing.is_some());
}

#[test]
fn wcs_file_round_trip() {
    let w: Wcs = TanWcs::new(
        SkyPos::from_degrees(150.0, 2.0).unwrap(),
        PixelPos::new(383.5, 383.5),
        [[-8.3e-4, 2e-5], [2e-5, 8.3e-4]],
    )
    .unwrap()
    .into();
    let path = std::env::temp_dir().join(format!("fayde-wcs-{}.wcs", std::process::id()));
    let mut cards = w.to_cards();
    cards.extend(fayde_fits::comment_cards("Solved by Fayde test."));
    write_wcs_file(&path, &cards).unwrap();
    let bytes = std::fs::read(&path).unwrap();
    assert_eq!(bytes.len() % 2880, 0);
    let f = FitsFile::open(&path).unwrap();
    let back = Wcs::from_header(&f.primary().header).unwrap();
    for (x, y) in [(0.0, 0.0), (767.0, 767.0)] {
        let p = PixelPos::new(x, y);
        assert!(back.pixel_to_sky(p).separation(w.pixel_to_sky(p)) < 1e-15);
    }
    std::fs::remove_file(path).unwrap();
}
