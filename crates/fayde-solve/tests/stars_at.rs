mod common;

use common::field;
use fayde_core::{PixelPos, SkyPos};
use fayde_solve::{stars_at, AtOptions, AtPosition};
use fayde_wcs::Parity;

#[test]
fn a_source_on_a_star_is_consistent_and_an_empty_spot_is_not() {
    let (cat, f) = field((150.0, 2.0), 37.0, Parity::Negative, 45);
    let (i, _) = f
        .source_star
        .iter()
        .enumerate()
        .find(|(_, s)| s.is_some())
        .expect("a real source");
    let src = &f.request.sources[i];
    let on_star = AtPosition::Pixel(PixelPos::new(src.x, src.y), Some(0.5));
    let near = f.truth.pixel_to_sky(PixelPos::new(src.x, src.y));
    let (ra, dec) = near.to_degrees();
    let off = AtPosition::Sky(SkyPos::from_degrees(ra, dec + 25.0 / 3600.0).unwrap(), None);
    let r = stars_at(
        &cat,
        &f.truth,
        &[on_star, off],
        &AtOptions {
            search_radius_arcsec: 60.0,
            ..AtOptions::default()
        },
    )
    .unwrap();
    assert!(r[0].consistent, "{:?}", r[0]);
    assert_eq!(Some(r[0].stars[0].id), f.source_star[i]);
    assert!(r[0].stars[0].separation_arcsec < 3.0 * r[0].stars[0].sigma_arcsec);
    assert!(r[0].pixel.is_some());
    let s = &r[1].stars[0];
    assert!(s.separation_arcsec > 20.0, "{s:?}");
    assert!(!r[1].consistent);
    assert!(s.p_chance > 0.0 && s.p_chance < 1.0);
    assert!(r[1].stars.len() <= 3);
    assert!(r[1]
        .stars
        .windows(2)
        .all(|w| w[0].separation_arcsec <= w[1].separation_arcsec));
    let without_stars = stars_at(
        &cat,
        &f.truth,
        &[on_star, off],
        &AtOptions {
            max_stars: 0,
            search_radius_arcsec: 60.0,
            ..AtOptions::default()
        },
    )
    .unwrap();
    for (result, expected) in without_stars.iter().zip(&r) {
        assert!(result.stars.is_empty());
        assert_eq!(result.consistent, expected.consistent);
    }
}

#[test]
fn bad_options_and_positions_are_errors() {
    let (cat, f) = field((150.0, 2.0), 37.0, Parity::Negative, 46);
    let p = [AtPosition::Pixel(PixelPos::new(10.0, 10.0), None)];
    for o in [
        AtOptions {
            search_radius_arcsec: 0.0,
            ..AtOptions::default()
        },
        AtOptions {
            default_sigma_px: f64::NAN,
            ..AtOptions::default()
        },
        AtOptions {
            wcs_sigma_arcsec: -1.0,
            ..AtOptions::default()
        },
    ] {
        assert!(stars_at(&cat, &f.truth, &p, &o).is_err());
    }
    let nan = [AtPosition::Pixel(PixelPos::new(f64::NAN, 1.0), None)];
    assert!(stars_at(&cat, &f.truth, &nan, &AtOptions::default()).is_err());
    assert!(stars_at(&cat, &f.truth, &[], &AtOptions::default())
        .unwrap()
        .is_empty());
}

#[test]
fn spread_positions_tiny_radii_and_radial_wcs_sigma() {
    let (cat, f) = field((150.0, 2.0), 37.0, Parity::Negative, 47);
    let (i, _) = f
        .source_star
        .iter()
        .enumerate()
        .find(|(_, s)| s.is_some())
        .unwrap();
    let src = &f.request.sources[i];
    let here = f.truth.pixel_to_sky(PixelPos::new(src.x, src.y));
    let (ra, dec) = here.to_degrees();
    let antipode = SkyPos::from_degrees((ra + 180.0) % 360.0, -dec).unwrap();
    let r = stars_at(
        &cat,
        &f.truth,
        &[
            AtPosition::Sky(here, Some(0.3)),
            AtPosition::Sky(antipode, None),
        ],
        &AtOptions::default(),
    )
    .unwrap();
    assert!(r[0].consistent && r[1].stars.is_empty());

    let tiny = stars_at(
        &cat,
        &f.truth,
        &[AtPosition::Sky(r[0].stars[0].pos, Some(0.001))],
        &AtOptions {
            search_radius_arcsec: 0.01,
            ..AtOptions::default()
        },
    )
    .unwrap();
    let s = &tiny[0].stars[0];
    assert!(s.p_chance < 1.0, "{s:?}");

    let exact = AtPosition::Sky(r[0].stars[0].pos, Some(0.0));
    let w = stars_at(
        &cat,
        &f.truth,
        &[exact],
        &AtOptions {
            wcs_sigma_arcsec: 2.0,
            ..AtOptions::default()
        },
    )
    .unwrap();
    let s = &w[0].stars[0];
    let cat_sigma = s.sigma_arcsec * s.sigma_arcsec - 2.0;
    assert!((-1e-9..0.01).contains(&cat_sigma), "{s:?}");
}
