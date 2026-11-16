//! Headers and pixel values that must give an error or a defined answer, never a panic or a
//! silently different WCS.

use fayde_core::{PixelPos, Portable, SkyPos};
use fayde_wcs::{Header, TanWcs, Wcs};

fn card(key: &str, value: &str) -> String {
    format!("{key:<8}= {value:>20}")
}

fn header(cards: &[(&str, &str)]) -> Header {
    let cards: Vec<String> = cards.iter().map(|(k, v)| card(k, v)).collect();
    Header::parse_cards(cards.iter().map(String::as_str))
}

const TAN: [(&str, &str); 6] = [
    ("CTYPE1", "'RA---TAN'"),
    ("CTYPE2", "'DEC--TAN'"),
    ("CRVAL1", "80.0"),
    ("CRVAL2", "10.0"),
    ("CRPIX1", "384.5"),
    ("CRPIX2", "384.5"),
];

fn with(extra: &[(&'static str, &'static str)]) -> Header {
    header(&[&TAN[..], extra].concat())
}

#[test]
fn pixel_to_sky_never_panics() {
    let w = TanWcs::new(
        SkyPos::from_degrees(80.0, 10.0).unwrap(),
        PixelPos::new(100.0, 100.0),
        [[-8e-4, 0.0], [0.0, 8e-4]],
    )
    .unwrap();
    // Offsets whose squares overflow tend to a point 90° from the reference point.
    for p in [
        PixelPos::new(1e308, 1e308),
        PixelPos::new(-1e300, 5.0),
        PixelPos::new(f64::INFINITY, 0.0),
    ] {
        let sep = w.pixel_to_sky(p).separation(w.crval()).to_degrees();
        assert!((sep - 90.0).abs() < 1e-6, "{p:?}: {sep}");
    }
    // NaN has no direction: the reference point.
    let s = w.pixel_to_sky(PixelPos::new(f64::NAN, 1.0));
    assert!(s.separation(w.crval()) < 1e-12);
    // An ordinary pixel is untouched by the guard: one pixel east of the reference point.
    let s = w.pixel_to_sky(PixelPos::new(101.0, 100.0));
    assert!((s.separation(w.crval()).to_degrees() - 8e-4).abs() < 1e-12);
}

#[test]
fn extreme_header_values_do_not_panic() {
    let cd = [("CD1_1", "-8.0E-4"), ("CD2_2", "8.0E-4")];
    // A reference pixel at the edge of the f64 range, and a pixel at the other edge.
    let h = header(
        &[
            &TAN[..4],
            &[("CRPIX1", "384.5"), ("CRPIX2", "-1e308")],
            &cd[..],
        ]
        .concat(),
    );
    let w = TanWcs::from_header(&h).unwrap();
    let _ = w.pixel_to_sky(PixelPos::new(100.0, 1e308));
    // A matrix element at the edge of the range.
    if let Ok(w) = TanWcs::from_header(&with(&[("CD1_1", "1e308"), ("CD2_2", "1e-308")])) {
        let _ = w.pixel_to_sky(PixelPos::new(1e6, 1e6));
    }
}

#[test]
fn non_finite_numbers_are_errors_by_name() {
    for bad in ["NaN", "inf", "-inf", "1e999"] {
        let e = TanWcs::from_header(&with(&[("CD1_1", bad), ("CD2_2", "8.0E-4")])).unwrap_err();
        assert!(e.to_string().contains("CD1_1"), "{bad}: {e}");
    }
    // A SIP coefficient that is not a number is an error, not an absent term.
    let sip = [
        ("CTYPE1", "'RA---TAN-SIP'"),
        ("CTYPE2", "'DEC--TAN-SIP'"),
        ("CRVAL1", "80.0"),
        ("CRVAL2", "10.0"),
        ("CRPIX1", "384.5"),
        ("CRPIX2", "384.5"),
        ("CD1_1", "-8.0E-4"),
        ("CD2_2", "8.0E-4"),
        ("A_ORDER", "2"),
        ("B_ORDER", "2"),
        ("A_2_0", "1.5E-06"),
    ];
    assert!(Wcs::from_header(&header(&[&sip[..], &[("B_0_2", "2.0E-06")]].concat())).is_ok());
    let e = Wcs::from_header(&header(&[&sip[..], &[("B_0_2", "NaN")]].concat())).unwrap_err();
    assert!(e.to_string().contains("B_0_2"), "{e}");
}

#[test]
fn the_first_of_a_duplicated_keyword_is_used() {
    let h = with(&[
        ("CD1_1", "-8.0E-4"),
        ("CD2_2", "8.0E-4"),
        ("CRVAL1", "10.0"),
    ]);
    let w = TanWcs::from_header(&h).unwrap();
    assert!((w.crval().to_degrees().0 - 80.0).abs() < 1e-12);
    // `set` replaces the value that lookups return.
    let mut h = h;
    h.set("CRVAL1", fayde_wcs::Value::Float(33.0));
    assert!((TanWcs::from_header(&h).unwrap().crval().to_degrees().0 - 33.0).abs() < 1e-12);
}

#[test]
fn crota2_rotates_the_cdelt_axes() {
    // Calabretta & Greisen 2002, §6.1: CD = [[c1·cos ρ, −c2·sin ρ], [c1·sin ρ, c2·cos ρ]].
    let (c1, c2) = (-8.3e-4f64, 8.3e-4f64);
    let (sin, cos) = 37.0f64.to_radians().psin_cos();
    let want = [[c1 * cos, -c2 * sin], [c1 * sin, c2 * cos]];
    let cdelt = [("CDELT1", "-8.3E-4"), ("CDELT2", "8.3E-4")];
    let got = TanWcs::from_header(&with(&[cdelt[0], cdelt[1], ("CROTA2", "37.0")]))
        .unwrap()
        .cd();
    for i in 0..2 {
        for j in 0..2 {
            assert!((got[i][j] - want[i][j]).abs() < 1e-15, "{got:?} {want:?}");
        }
    }
    // Fayde's roll, atan2(CD1_2, CD2_2), is −ρ for these CDELT signs (east to the left).
    assert!((got[0][1].patan2(got[1][1]).to_degrees() + 37.0).abs() < 1e-9);
    // Without CROTA2 the axes are unrotated; CROTA1 is used when it is the only one given.
    let plain = TanWcs::from_header(&with(&cdelt)).unwrap().cd();
    assert_eq!(plain, [[c1, 0.0], [0.0, c2]]);
    let via1 = TanWcs::from_header(&with(&[cdelt[0], cdelt[1], ("CROTA1", "37.0")]))
        .unwrap()
        .cd();
    assert_eq!(via1, got);
    // A PC matrix takes precedence over CROTA2, as in wcslib.
    let pc = TanWcs::from_header(&with(&[
        cdelt[0],
        cdelt[1],
        ("PC1_1", "1.0"),
        ("PC2_2", "1.0"),
        ("CROTA2", "37.0"),
    ]))
    .unwrap()
    .cd();
    assert_eq!(pc, plain);
}
