use fayde_fits::header_priors;
use fayde_wcs::Header;

fn header(cards: &[(&str, &str)]) -> Header {
    let cards: Vec<String> = cards
        .iter()
        .map(|(k, v)| format!("{k:<8}= {v:>20}"))
        .collect();
    Header::parse_cards(cards.iter().map(String::as_str))
}

#[test]
fn implausible_epochs_are_reported_and_skipped() {
    let p = header_priors(&header(&[("MJD-OBS", "2459175.2")]));
    assert!(p.epoch.is_none());
    assert!(
        p.warnings.iter().any(|w| w.contains("MJD-OBS")),
        "{:?}",
        p.warnings
    );
    let p = header_priors(&header(&[
        ("MJD-OBS", "2459175.2"),
        ("DATE-OBS", "'2020-11-21T16:46:40'"),
    ]));
    let (e, origin) = p.epoch.unwrap();
    assert!((e.year() - 2020.89).abs() < 0.01 && origin.starts_with("DATE-OBS"));
    for date in ["'0001-01-01T00:00:00'", "'99999999999999999-01-01'"] {
        assert!(
            header_priors(&header(&[("DATE-OBS", date)]))
                .epoch
                .is_none(),
            "{date}"
        );
    }
    for date in [
        "'0001-01-01T00:00:00'",
        "'99999999999999999-01-01'",
        "'yesterday'",
    ] {
        let p = header_priors(&header(&[("DATE-OBS", date)]));
        assert!(
            p.warnings.iter().any(|w| w.contains("DATE-OBS")),
            "{date}: {:?}",
            p.warnings
        );
    }
}

#[test]
fn coordinates_in_another_frame_are_not_used_as_icrs() {
    let radec = [("RA", "'11 30 00.0'"), ("DEC", "'+02 00 00'")];
    assert!(header_priors(&header(&radec)).pointing.is_some());
    for frame in [
        &[("RADESYS", "'FK4'"), ("EQUINOX", "1950.0")][..],
        &[("EQUINOX", "2026.0")][..],
        &[("RADESYS", "'GAPPT'")][..],
    ] {
        let p = header_priors(&header(&[&radec[..], frame].concat()));
        assert!(p.pointing.is_none(), "{frame:?}");
        assert!(
            p.warnings.iter().any(|w| w.contains("not ICRS/J2000")),
            "{:?}",
            p.warnings
        );
    }
    for frame in [
        &[("RADESYS", "'ICRS'")][..],
        &[("RADESYS", "'FK5'"), ("EQUINOX", "2000.0")][..],
        &[("EPOCH", "2020.5")][..],
    ] {
        let p = header_priors(&header(&[&radec[..], frame].concat()));
        assert!(p.pointing.is_some(), "{frame:?}");
    }
}

#[test]
fn out_of_range_positions_are_rejected_not_wrapped() {
    for (ra, dec) in [
        ("'25 00 00'", "'+02 00 00'"),
        ("400.0", "2.0"),
        ("-10.0", "2.0"),
    ] {
        let p = header_priors(&header(&[("RA", ra), ("DEC", dec)]));
        assert!(p.pointing.is_none(), "{ra}");
        assert!(
            p.warnings
                .iter()
                .any(|w| w.contains("RA/DEC is not a valid position")),
            "{:?}",
            p.warnings
        );
    }
    let p = header_priors(&header(&[
        ("RA", "'25 00 00'"),
        ("DEC", "'+02 00 00'"),
        ("OBJCTRA", "'11 30 00'"),
        ("OBJCTDEC", "'+02 00 00'"),
    ]));
    assert_eq!(p.pointing.unwrap().1, "OBJCTRA/OBJCTDEC");
}

#[test]
fn a_non_finite_header_number_is_not_a_prior() {
    let p = header_priors(&header(&[
        ("CTYPE1", "'RA---TAN'"),
        ("CTYPE2", "'DEC--TAN'"),
        ("CRVAL1", "NaN"),
        ("CRVAL2", "2.0"),
        ("PIXSCALE", "inf"),
    ]));
    assert!(p.pointing.is_none() && p.scale_arcsec.is_none());
}
