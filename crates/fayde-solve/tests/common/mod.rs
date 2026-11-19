#![allow(dead_code)]

use fayde_catalog::{memory::MemoryCatalog, Star};
use fayde_core::Portable;
use fayde_core::{PixelPos, SkyPos};
use fayde_solve::{Frame, PointingPrior, Policy, Priors, ScalePrior, SolveRequest, Source};
use fayde_wcs::{Parity, Poly, Sip, TanWcs, Wcs};

pub struct Rng(pub u64);
impl Rng {
    pub fn uniform(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64 + 0.5) / (1u64 << 53) as f64
    }
    pub fn normal(&mut self) -> f64 {
        let (u, v) = (self.uniform(), self.uniform());
        (-2.0 * u.pln()).sqrt() * (2.0 * std::f64::consts::PI * v).pcos()
    }
}

pub fn sky(centre: SkyPos, radius_deg: f64, rng: &mut Rng) -> Vec<Star> {
    let r = radius_deg.to_radians();
    let area_deg2 =
        2.0 * std::f64::consts::PI * (1.0 - r.pcos()) * (180.0 / std::f64::consts::PI).ppowi(2);
    let n = (400.0 * area_deg2) as usize;
    let c = centre.to_unit_vector().0;
    let e = {
        let (s, co) = centre.ra.psin_cos();
        [-s, co, 0.0]
    };
    let nn = [
        c[1] * e[2] - c[2] * e[1],
        c[2] * e[0] - c[0] * e[2],
        c[0] * e[1] - c[1] * e[0],
    ];
    (0..n)
        .map(|i| {
            let cos_t = 1.0 - rng.uniform() * (1.0 - r.pcos());
            let sin_t = (1.0 - cos_t * cos_t).sqrt();
            let phi = 2.0 * std::f64::consts::PI * rng.uniform();
            let v =
                [0, 1, 2].map(|k| c[k] * cos_t + sin_t * (phi.pcos() * e[k] + phi.psin() * nn[k]));
            let u = rng.uniform();
            let lo = 10f64.ppowf(0.3 * 6.0);
            let hi = 10f64.ppowf(0.3 * 16.0);
            let mag = ((lo + u * (hi - lo)).plog10() / 0.3) as f32;
            Star {
                id: i as u64 + 1,
                pos: fayde_core::UnitVector(v).to_sky(),
                mag: Some(mag),
                pos_err_mas: [20.0, 20.0],
                pm_mas_yr: [Some(0.0), Some(0.0)],
            }
        })
        .collect()
}

pub struct Field {
    pub request: SolveRequest,
    pub truth: Wcs,
    pub source_star: Vec<Option<u64>>,
}

pub fn field(
    centre: (f64, f64),
    roll_deg: f64,
    parity: Parity,
    seed: u64,
) -> (MemoryCatalog, Field) {
    field_with_noise(centre, roll_deg, parity, seed, 0.3)
}

pub fn field_with_noise(
    centre: (f64, f64),
    roll_deg: f64,
    parity: Parity,
    seed: u64,
    noise_px: f64,
) -> (MemoryCatalog, Field) {
    field_full(centre, roll_deg, parity, seed, noise_px, 0.0)
}

pub fn field_scaled(
    centre: (f64, f64),
    roll_deg: f64,
    parity: Parity,
    seed: u64,
    scale: f64,
) -> (MemoryCatalog, Field) {
    field_any(centre, roll_deg, parity, seed, 0.3, 0.0, scale)
}

pub fn field_full(
    centre: (f64, f64),
    roll_deg: f64,
    parity: Parity,
    seed: u64,
    noise_px: f64,
    barrel: f64,
) -> (MemoryCatalog, Field) {
    field_any(centre, roll_deg, parity, seed, noise_px, barrel, 3.0)
}

fn field_any(
    centre: (f64, f64),
    roll_deg: f64,
    parity: Parity,
    seed: u64,
    noise_px: f64,
    barrel: f64,
    scale_arcsec: f64,
) -> (MemoryCatalog, Field) {
    let mut rng = Rng(seed);
    let c = SkyPos::from_degrees(centre.0, centre.1).unwrap();
    let stars = sky(c, 1.6, &mut rng);
    let (w, h) = (1024u32, 1024u32);
    let s = scale_arcsec / 3600.0;
    let (sr, cr) = roll_deg.to_radians().psin_cos();
    let flip = if parity == Parity::Negative {
        -1.0
    } else {
        1.0
    };
    let cd = [[flip * s * cr, s * sr], [-flip * s * sr, s * cr]];
    let tan = TanWcs::new(c, PixelPos::new(511.5, 511.5), cd).unwrap();
    let truth = if barrel == 0.0 {
        Wcs::from(tan)
    } else {
        let empty = Poly {
            order: 0,
            terms: vec![],
        };
        Wcs {
            tan,
            sip: Some(Sip {
                a: Poly {
                    order: 3,
                    terms: vec![(3, 0, barrel), (1, 2, barrel)],
                },
                b: Poly {
                    order: 3,
                    terms: vec![(2, 1, barrel), (0, 3, barrel)],
                },
                ap: empty.clone(),
                bp: empty,
            }),
        }
    };
    assert_eq!(truth.parity(), parity);

    let mut items: Vec<(f64, Source, Option<u64>)> = Vec::new();
    for st in &stars {
        let m = st.mag.unwrap();
        if m > 14.5 {
            continue;
        }
        let Some(p) = truth.sky_to_pixel(st.pos) else {
            continue;
        };
        if !(p.x >= 0.0 && p.y >= 0.0 && p.x < 1023.0 && p.y < 1023.0) || rng.uniform() < 0.2 {
            continue;
        }
        let flux = 10f64.ppowf(-0.4 * (f64::from(m) - 20.0)) * (1.0 + 0.1 * rng.normal());
        let src = Source {
            x: p.x + noise_px * rng.normal(),
            y: p.y + noise_px * rng.normal(),
            flux: Some(flux),
            sigma_px: Some(noise_px.max(0.05)),
            p_star: None,
            truncated: false,
        };
        items.push((flux, src, Some(st.id)));
    }
    let n_false = (items.len() as f64 * 0.3) as usize;
    for _ in 0..n_false {
        let pick = items[(rng.uniform() * items.len() as f64) as usize % items.len()].0;
        let src = Source {
            x: rng.uniform() * 1023.0,
            y: rng.uniform() * 1023.0,
            flux: Some(pick * (0.5 + rng.uniform())),
            sigma_px: Some(0.3),
            p_star: None,
            truncated: false,
        };
        items.push((src.flux.unwrap(), src, None));
    }
    items.sort_by_key(|it| ((it.1.x * 1e6) as i64).wrapping_mul(2654435761) % 1_000_003);
    let source_star = items.iter().map(|it| it.2).collect();
    let sources = items.into_iter().map(|it| it.1).collect();

    let off = SkyPos::from_degrees(
        centre.0 + 0.2 / centre.1.to_radians().pcos(),
        centre.1 - 0.05,
    )
    .unwrap();
    let request = SolveRequest {
        frame: Frame {
            width: w,
            height: h,
            trail: None,
        },
        sources,
        priors: Priors {
            pointing: Some(PointingPrior {
                centre: off,
                sigma: 0.2f64.to_radians(),
                max_radius: 0.5f64.to_radians(),
            }),
            scale: Some(ScalePrior {
                arcsec_per_px: scale_arcsec,
                sigma: 0.03,
            }),
            roll: None,
            parity: None,
        },
        policy: Policy::default(),
    };
    (
        MemoryCatalog::new(stars),
        Field {
            request,
            truth,
            source_star,
        },
    )
}
