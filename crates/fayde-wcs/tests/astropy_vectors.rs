use fayde_core::{to_arcsec, PixelPos, SkyPos};
use fayde_wcs::{Header, Parity, TanWcs};
use serde::Deserialize;

#[derive(Deserialize)]
struct Fixture {
    pixel_origin: u8,
    cases: Vec<Case>,
}

#[derive(Deserialize)]
struct Case {
    name: String,
    crval_deg: [f64; 2],
    crpix_fits: [f64; 2],
    cd: [[f64; 2]; 2],
    pixels_zero_based: Vec<[f64; 2]>,
    world_deg: Vec<[f64; 2]>,
    pixel_scale_arcsec: f64,
    header_cards: Vec<String>,
}

fn fixture() -> Fixture {
    let f: Fixture = serde_json::from_str(include_str!("fixtures/tan_astropy.json")).unwrap();
    assert_eq!(f.pixel_origin, 0);
    assert_eq!(f.cases.len(), 6);
    f
}

const SKY_TOL_ARCSEC: f64 = 1e-6;
const PIX_TOL: f64 = 1e-8;

fn check(case: &Case, w: &TanWcs, source: &str) {
    for (p, expect) in case.pixels_zero_based.iter().zip(&case.world_deg) {
        let got = w.pixel_to_sky(PixelPos::new(p[0], p[1]));
        let want = SkyPos::from_degrees(expect[0], expect[1]).unwrap();
        let err = to_arcsec(got.separation(want));
        assert!(
            err < SKY_TOL_ARCSEC,
            "{} ({source}) pixel {p:?}: {err} arcsec",
            case.name
        );

        let back = w.sky_to_pixel(want).expect("in front of tangent plane");
        assert!(
            (back.x - p[0]).abs() < PIX_TOL && (back.y - p[1]).abs() < PIX_TOL,
            "{} ({source}) inverse {p:?} -> {back:?}",
            case.name
        );
    }
    assert!((w.pixel_scale_arcsec() - case.pixel_scale_arcsec).abs() < 1e-9);
}

#[test]
fn matches_astropy_from_parameters() {
    for case in fixture().cases {
        let w = TanWcs::new(
            SkyPos::from_degrees(case.crval_deg[0], case.crval_deg[1]).unwrap(),
            PixelPos::new(case.crpix_fits[0] - 1.0, case.crpix_fits[1] - 1.0),
            case.cd,
        )
        .unwrap();
        check(&case, &w, "parameters");
        let expected_parity = if case.name == "flipped_parity" {
            Parity::Positive
        } else {
            Parity::Negative
        };
        assert_eq!(w.parity(), expected_parity, "{}", case.name);
    }
}

#[test]
fn matches_astropy_from_astropy_header() {
    for case in fixture().cases {
        let h = Header::parse_cards(case.header_cards.iter().map(String::as_str));
        let w = TanWcs::from_header(&h).unwrap();
        check(&case, &w, "astropy header");
    }
}

#[test]
fn own_cards_round_trip_through_astropy_cases() {
    for case in fixture().cases {
        let h = Header::parse_cards(case.header_cards.iter().map(String::as_str));
        let w = TanWcs::from_header(&h).unwrap();
        let cards = w.to_cards();
        let w2 =
            TanWcs::from_header(&Header::parse_cards(cards.iter().map(String::as_str))).unwrap();
        check(&case, &w2, "fayde cards");
    }
}

#[test]
fn pole_longitudes_match_wcslib() {
    for (dec, pole, ra_x, ra_y) in [
        (90_i32, None, 132., 42.),
        (90, Some(0), 132., 42.),
        (90, Some(180), 312., 222.),
        (-90, None, 312., 42.),
        (-90, Some(0), 132., 222.),
        (-90, Some(180), 312., 42.),
    ] {
        let dec_card = format!("CRVAL2  = {dec}");
        let pole_card = pole.map(|v| format!("LONPOLE = {v}"));
        let cards = [
            "CTYPE1  = 'RA---TAN'",
            "CTYPE2  = 'DEC--TAN'",
            "CRVAL1  = 42.",
            &dec_card,
            "CRPIX1  = 1.",
            "CRPIX2  = 1.",
            "CDELT1  = -0.01",
            "CDELT2  = 0.01",
        ];
        let w = TanWcs::from_header(&Header::parse_cards(
            cards.into_iter().chain(pole_card.as_deref()),
        ))
        .unwrap();
        for (x, y, ra) in [(10., 0., ra_x), (0., 10., ra_y)] {
            let sky = w.pixel_to_sky(PixelPos::new(x, y));
            let want = SkyPos::from_degrees(ra, f64::from(dec.signum()) * 89.900000101539).unwrap();
            assert!(to_arcsec(sky.separation(want)) < 1e-6);
        }
    }
}

#[test]
fn high_order_sip_headers_round_trip_across_the_frame() {
    use fayde_wcs::Wcs;
    for order in [8, 9] {
        let cards = [
            "CTYPE1  = 'RA---TAN-SIP'",
            "CTYPE2  = 'DEC--TAN-SIP'",
            "CRVAL1  = 42.",
            "CRVAL2  = 10.",
            "CRPIX1  = 128.",
            "CRPIX2  = 128.",
            "CDELT1  = -0.01",
            "CDELT2  = 0.01",
            "A_ORDER = 2",
            "A_2_0   = 1e-5",
            "B_0_2   = -2e-5",
        ];
        let terms = [
            format!("B_ORDER = {order}"),
            format!("B_0_{order}   = 1e-20"),
        ];
        let mut w = Wcs::from_header(&Header::parse_cards(
            cards.into_iter().chain(terms.iter().map(String::as_str)),
        ))
        .unwrap();
        assert_eq!(w.sip_order(), order);
        w.fit_inverse(256., 256., w.sip_order() + 2);
        assert_eq!(w.sip.as_ref().unwrap().ap.order, 9);
        let restored = Wcs::from_header(&Header::parse_cards(
            w.to_cards().iter().map(String::as_str),
        ))
        .unwrap();
        for x in [0., 64., 127.5, 255.] {
            for y in [0., 127.5, 255.] {
                let back = restored
                    .sky_to_pixel(w.pixel_to_sky(PixelPos::new(x, y)))
                    .unwrap();
                assert!((back.x - x).abs().max((back.y - y).abs()) < 1e-6);
            }
        }
    }
}
