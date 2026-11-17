use super::*;
use crate::request::{Frame, Policy, Priors, Source};
use fayde_catalog::memory::MemoryCatalog;
use fayde_core::Portable;

#[test]
fn crowded_local_refit_resists_offset_false_sources() {
    let truth = TanWcs::new(
        SkyPos::from_degrees(150.0, 2.0).unwrap(),
        PixelPos::new(511.5, 511.5),
        [[-3.0 / 3600.0, 0.0], [0.0, 3.0 / 3600.0]],
    )
    .unwrap();
    let mut stars = Vec::new();
    let mut sources = Vec::new();
    for i in 0..42 {
        let x = 360.0 + 40.0 * (i % 7) as f64;
        let y = 400.0 + 40.0 * (i / 7) as f64;
        stars.push(Star {
            id: i + 1,
            pos: truth.pixel_to_sky(PixelPos::new(x, y)),
            mag: Some(10.0),
            pos_err_mas: [60.0, 60.0],
            pm_mas_yr: [Some(0.0), Some(0.0)],
        });
        sources.push(Source {
            x,
            y,
            flux: Some(1.0),
            sigma_px: Some(0.2),
            p_star: None,
            truncated: false,
        });
        // Brighter false detections steal two thirds of the initial reference claims.
        // They lie 8–11 px from the truth, inside the first refit radius.
        if i % 3 != 0 {
            sources.push(Source {
                x: x + 8.0 + (i % 4) as f64,
                y,
                flux: Some(2.0),
                sigma_px: Some(0.2),
                p_star: None,
                truncated: false,
            });
        }
    }
    let cat = MemoryCatalog::new(stars);
    let req = SolveRequest {
        frame: Frame {
            width: 1024,
            height: 1024,
            trail: None,
        },
        sources,
        priors: Priors::default(),
        policy: Policy::default(),
    };
    let ctx = new_ctx(&req, &cat).unwrap();
    let set = RefSet::from_cone(&cat, truth.crval(), 1f64.to_radians(), &ctx.query).unwrap();
    let fitted = ctx
        .local_refit(&set, truth)
        .expect("refit must retain real matches");
    let centre = fitted.sky_to_pixel(truth.crval()).unwrap();
    let error = (centre.x - 511.5).phypot(centre.y - 511.5);
    eprintln!("crowded local refit centre error: {error:.6} px");
    assert!(error < 0.1, "false detections moved the fit by {error} px");
}

#[test]
fn pair_order_alternates_and_exhausts_both_sides() {
    let pairs = [
        (1.0, 0, 0),
        (2.0, 1, 1),
        (3.0, 2, 2),
        (4.0, 3, 3),
        (5.0, 4, 4),
    ];
    let order = |mean| {
        pairs_from_mean(&pairs, mean)
            .map(|p| p.1)
            .collect::<Vec<_>>()
    };
    assert_eq!(order(3.1), vec![2, 3, 1, 4, 0]);
    assert_eq!(order(2.9), vec![2, 1, 3, 0, 4]);
    assert_eq!(order(0.0), vec![0, 1, 2, 3, 4]);
    assert_eq!(order(6.0), vec![4, 3, 2, 1, 0]);
    assert_eq!(pairs_from_mean(&[], 3.0).count(), 0);
}

#[test]
fn similarity_recovers_both_parities_and_rejects_shear() {
    let pts = [
        [120.0, 210.0],
        [810.0, 240.0],
        [670.0, 730.0],
        [310.0, 610.0],
    ];
    let crpix = PixelPos::new(511.5, 511.5);
    for sign in [-1.0, 1.0] {
        let truth = TanWcs::new(
            SkyPos::from_degrees(150.0, 2.0).unwrap(),
            crpix,
            [[sign * 0.0008, 0.0002], [-sign * 0.0002, 0.0008]],
        )
        .unwrap();
        let skies = pts.map(|p| truth.pixel_to_sky(PixelPos::new(p[0], p[1])));
        let seeds = pattern_seeds(&pts, &skies, 0.002, 0.3, crpix);
        assert_eq!(
            seeds.len(),
            1,
            "asymmetric pattern has a unique correspondence"
        );
        let (rms, seed) = seeds[0];
        assert!(rms < 0.01);
        assert_eq!(seed.parity(), truth.parity());
        assert!((seed.pixel_scale_arcsec() / truth.pixel_scale_arcsec() - 1.0).abs() < 1e-5);
        let sheared = pts.map(|p| [p[0] + 0.01 * (p[1] - crpix.y), p[1]]);
        // The loose edge tolerance admits this shape, but the four-parameter fit
        // cannot absorb shear as the former affine seed could.
        assert!(pattern_seeds(&sheared, &skies, 0.02, 0.1, crpix).is_empty());
    }
    let sky = SkyPos::from_degrees(150.0, 2.0).unwrap();
    assert!(pattern_seeds(&[[0.0; 2]; 4], &[sky; 4], 0.02, 0.3, crpix).is_empty());
}

#[test]
fn similarity_seed_is_fitted_at_its_own_crval_near_the_pole() {
    // Two stars in one corner of a frame at Dec 88.5°: their centroid is 0.4° from the frame
    // centre. A transform fitted in the centroid's tangent plane and reused at CRVAL is rotated
    // by (east offset × tan δ), tens of pixels across the frame here.
    let crpix = PixelPos::new(511.5, 511.5);
    for sign in [-1.0, 1.0] {
        let truth = TanWcs::new(
            SkyPos::from_degrees(40.0, 88.5).unwrap(),
            crpix,
            [[sign * 0.0008, 0.0002], [-sign * 0.0002, 0.0008]],
        )
        .unwrap();
        let pts = [[60.0, 90.0], [240.0, 150.0]];
        let pairs = pts.map(|p| (p, truth.pixel_to_sky(PixelPos::new(p[0], p[1]))));
        let (seed, rms) = similarity_seed(&pairs, truth.parity(), crpix).unwrap();
        assert!(rms < 1e-6, "two points are fitted exactly: {rms}");
        for corner in [[0.0, 0.0], [1023.0, 1023.0], [1023.0, 0.0], [0.0, 1023.0]] {
            let p = PixelPos::new(corner[0], corner[1]);
            let err_px =
                fayde_core::to_arcsec(seed.pixel_to_sky(p).separation(truth.pixel_to_sky(p)))
                    / truth.pixel_scale_arcsec();
            assert!(err_px < 0.01, "corner {corner:?} off by {err_px} px");
        }
    }
}

/// A [`PairTable`] lookup returns what the sorted list of every pair returns: the pairs a k-d
/// tree search finds within reach, by separation and then by star numbers, cut to the range.
#[test]
fn pair_table_lookup_equals_the_sorted_list() {
    let mut state = 4242u64;
    let mut rnd = || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (state >> 11) as f64 / (1u64 << 53) as f64
    };
    // Stars over a 5° cap, some of them twice (zero separation) and some in a tight clump.
    let mut vecs: Vec<[f64; 3]> = (0..700)
        .map(|i| {
            let spread = if i % 9 == 0 { 0.02 } else { 5.0 };
            SkyPos::from_degrees(40.0 + spread * rnd(), 88.0 - spread * rnd())
                .unwrap()
                .to_unit_vector()
                .0
        })
        .collect();
    vecs.extend_from_within(..5);
    for reach_deg in [0.9f64, 0.05, 20.0] {
        let reach = chord(reach_deg.to_radians());
        // The list as it was built before: every star's neighbours from a k-d tree.
        let tree = KdTree::build(vecs.clone());
        let mut all: Vec<(f64, u32, u32)> = Vec::new();
        for i in 0..vecs.len() {
            for j in tree.within(vecs[i], reach) {
                if j > i {
                    all.push((angle(&vecs[i], &vecs[j]), i as u32, j as u32));
                }
            }
        }
        all.sort_by(|a, b| a.0.total_cmp(&b.0).then((a.1, a.2).cmp(&(b.1, b.2))));
        let mut table = PairTable::new(&vecs, reach);
        assert_eq!(table.pairs.len(), all.len(), "reach {reach_deg}°");
        let widest = all.last().map_or(0.0, |p| p.0);
        for round in 0..400 {
            let mid = widest * rnd();
            let half = match round % 5 {
                0 => 0.0,
                1 => widest,
                _ => widest * 0.02 * rnd(),
            };
            let (lo, hi) = match round % 11 {
                // Ends that fall exactly on a pair.
                0 if !all.is_empty() => {
                    let k = (rnd() * all.len() as f64) as usize;
                    (all[k].0, all[(k + 30).min(all.len() - 1)].0)
                }
                1 => (-1.0, mid),
                2 => (mid, 9.0),
                _ => (mid - half, mid + half),
            };
            let start = all.partition_point(|c| c.0 < lo);
            let end = all.partition_point(|c| c.0 <= hi);
            assert_eq!(
                table.between(lo, hi),
                all[start..end],
                "reach {reach_deg}°, {lo} to {hi}"
            );
        }
    }
    // No stars, one star, no reach.
    assert!(PairTable::new(&[], 0.1).between(0.0, 1.0).is_empty());
    assert!(PairTable::new(&vecs[..1], 0.1).between(0.0, 1.0).is_empty());
    let mut none = PairTable::new(&vecs[..40], 0.0);
    assert!(none.between(1e-6, 1.0).is_empty());
}
