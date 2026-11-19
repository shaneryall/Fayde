use fayde_core::Portable;
use std::path::PathBuf;

use fayde_catalog::{ucac4::Ucac4Catalog, Catalog, StarQuery};
use fayde_core::{arcsec, Epoch, PixelPos, SkyPos};
use fayde_index::{build, image_shape, BuildParams, Index};
use fayde_wcs::TanWcs;

fn catalog() -> Option<Ucac4Catalog> {
    let dir = std::env::var_os("FAYDE_UCAC4_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data/ucac4"));
    let ok = (458..=463).all(|z| dir.join(format!("u4b/z{z}")).is_file());
    if !ok {
        panic!("UCAC4 zones z458–z463 not at {}", dir.display());
    }
    Some(Ucac4Catalog::open(dir).unwrap())
}

const CENTRE: (f64, f64) = (150.0, 2.0);
const REGION_DEG: f64 = 0.38;

fn regional_index(cat: &Ucac4Catalog) -> Index {
    let params = BuildParams {
        fov_deg: 0.4,
        min_edge_deg: 0.02,
        stars_per_fov: 20.0,
        mag_limit: Some(15.0),
        epoch_year: 2025.0,
        region: Some([CENTRE.0, CENTRE.1, REGION_DEG]),
        ..BuildParams::default()
    };
    let query = StarQuery {
        epoch: Epoch::julian_year(params.epoch_year),
        ..StarQuery::default()
    };
    let centre = SkyPos::from_degrees(CENTRE.0, CENTRE.1).unwrap();
    let stars = cat.cone(centre, REGION_DEG.to_radians(), &query).unwrap();
    build(&stars, cat.id(), &params).unwrap()
}

fn noise(seed: &mut u64, a: f64) -> f64 {
    *seed = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    ((*seed >> 11) as f64 / (1u64 << 53) as f64 * 2.0 - 1.0) * a
}

#[test]
#[ignore = "needs UCAC4 data: set FAYDE_UCAC4_DIR, run with --include-ignored"]
fn image_patterns_find_their_catalogue_stars() {
    let Some(cat) = catalog() else { return };
    let index = regional_index(&cat);
    let m = index.meta();
    eprintln!(
        "index: {} stars, {} patterns, {} slots, nside {}",
        m.star_count, m.pattern_count, m.table_size, m.healpix_nside
    );
    assert!(m.star_count > 50 && m.pattern_count > 100);

    let scale = 3.0 / 3600.0;
    let (s, c) = 37f64.to_radians().psin_cos();
    let wcs = TanWcs::new(
        SkyPos::from_degrees(CENTRE.0 + 0.03, CENTRE.1 - 0.02).unwrap(),
        PixelPos::new(149.5, 149.5),
        [[-scale * c, scale * s], [scale * s, scale * c]],
    )
    .unwrap();
    let mut seed = 7u64;
    let mut image: Vec<(PixelPos, u32)> = index
        .stars()
        .iter()
        .enumerate()
        .filter_map(|(i, st)| {
            let sky = fayde_core::UnitVector(st.vec).to_sky();
            let p = wcs.sky_to_pixel(sky)?;
            ((0.0..300.0).contains(&p.x) && (0.0..300.0).contains(&p.y)).then(|| {
                (
                    PixelPos::new(p.x + noise(&mut seed, 0.3), p.y + noise(&mut seed, 0.3)),
                    i as u32,
                )
            })
        })
        .collect();
    image.sort_by(|a, b| {
        index.stars()[a.1 as usize]
            .mag
            .total_cmp(&index.stars()[b.1 as usize].mag)
    });
    image.truncate(8);
    assert!(
        image.len() >= 6,
        "too few index stars in the synthetic frame"
    );

    let px_rad = arcsec(3.0);
    let (mut tried, mut found, mut wrong_best) = (0, 0, 0);
    let n = image.len();
    for a in 0..n {
        for b in a + 1..n {
            for c2 in b + 1..n {
                for d in c2 + 1..n {
                    let pick = [image[a], image[b], image[c2], image[d]];
                    let pts = pick.map(|(p, _)| [p.x * px_rad, p.y * px_rad]);
                    let (shape, order) = image_shape(&pts);
                    let truth = order.map(|k| pick[k].1);
                    let range = Some((shape.largest * 0.9, shape.largest * 1.1));
                    let matches = index.lookup(&shape, 0.01, range);
                    tried += 1;
                    if matches.iter().any(|m| m.stars == truth) {
                        found += 1;
                        if matches[0].stars != truth {
                            wrong_best += 1;
                        }
                    }
                }
            }
        }
    }
    eprintln!("{found}/{tried} image patterns found in the index ({wrong_best} not ranked first)");
    assert!(found >= 3, "found {found}/{tried}");
    assert_eq!(wrong_best, 0);
}

#[test]
#[ignore = "needs UCAC4 data: set FAYDE_UCAC4_DIR, run with --include-ignored"]
fn file_round_trip_and_corruption() {
    let Some(cat) = catalog() else { return };
    let index = regional_index(&cat);
    let path = std::env::temp_dir().join(format!("fayde-test-{}.fidx", std::process::id()));
    index.write(&path).unwrap();
    let back = Index::read(&path).unwrap();
    assert_eq!(back, index);
    assert_eq!(back.content_hash(), index.content_hash());

    let mut bytes = std::fs::read(&path).unwrap();
    let mid = bytes.len() / 2;
    bytes[mid] ^= 0x01;
    assert!(Index::from_bytes(&bytes).unwrap_err().contains("hash"));
    bytes.truncate(bytes.len() - 100);
    assert!(Index::from_bytes(&bytes).is_err());
    std::fs::remove_file(&path).unwrap();
}

#[test]
#[ignore = "needs UCAC4 data: set FAYDE_UCAC4_DIR, run with --include-ignored"]
fn build_is_deterministic() {
    let Some(cat) = catalog() else { return };
    assert_eq!(
        regional_index(&cat).content_hash(),
        regional_index(&cat).content_hash()
    );
}
