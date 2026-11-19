use fayde_core::Portable;
use std::time::Instant;

use fayde_catalog::{ucac4::Ucac4Catalog, Catalog, StarQuery};
use fayde_core::{to_arcsec, Epoch, PixelPos, SkyPos, UnitVector};
use fayde_index::Index;
use fayde_solve::{solve_with_index, Frame, Outcome, Policy, Priors, SolveRequest, Source};
use fayde_wcs::TanWcs;

use crate::{num, options, seconds};

struct Rng(u64);
impl Rng {
    fn uniform(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64 + 0.5) / (1u64 << 53) as f64
    }
    fn normal(&mut self) -> f64 {
        let (u, v) = (self.uniform(), self.uniform());
        (-2.0 * u.pln()).sqrt() * (2.0 * std::f64::consts::PI * v).pcos()
    }
}

pub(crate) fn run(args: &[&str]) -> Result<(), String> {
    let o = options(
        args,
        &[
            "ucac4",
            "index",
            "frames",
            "seed",
            "size",
            "scale",
            "image-mag",
            "deadline",
            "pattern-sources",
        ],
        &[],
    )?;
    let deadline = Some(seconds(&o, "deadline")?.unwrap_or(std::time::Duration::from_secs(60)));
    let pattern_sources: usize =
        num(&o, "pattern-sources")?.unwrap_or(Policy::default().pattern_sources);
    let cat = Ucac4Catalog::open(o.get("ucac4").ok_or("--ucac4 is required")?)
        .map_err(|e| e.to_string())?;
    let index =
        Index::read(o.get("index").ok_or("--index is required")?).map_err(|e| e.to_string())?;
    let frames: usize = num(&o, "frames")?.unwrap_or(100);
    let mut rng = Rng(num(&o, "seed")?.unwrap_or(1));
    let size: u32 = num(&o, "size")?.unwrap_or(1024);
    let scale: f64 = num(&o, "scale")?.unwrap_or(3.0);
    let image_mag: f32 = num(&o, "image-mag")?.unwrap_or(15.0);
    let epoch = Epoch::julian_year(index.meta().params.epoch_year);
    let n = f64::from(size);

    let (mut real, mut solved, mut correct, mut wrong, mut too_few) = (0, 0, 0, 0, 0);
    let (mut random, mut random_solved) = (0, 0);
    let mut times = Vec::new();
    let mut hyps = Vec::new();
    for f in 0..frames {
        let is_random = f % 2 == 1;
        let z = 2.0 * rng.uniform() - 1.0;
        let phi = 2.0 * std::f64::consts::PI * rng.uniform();
        let s = (1.0 - z * z).sqrt();
        let centre = UnitVector([s * phi.pcos(), s * phi.psin(), z]).to_sky();
        let roll = 2.0 * std::f64::consts::PI * rng.uniform();
        let flip = if rng.uniform() < 0.5 { -1.0 } else { 1.0 };
        let sd = scale / 3600.0;
        let (sr, cr) = roll.psin_cos();
        let truth = TanWcs::new(
            centre,
            PixelPos::new((n - 1.0) / 2.0, (n - 1.0) / 2.0),
            [[flip * sd * cr, sd * sr], [-flip * sd * sr, sd * cr]],
        )
        .map_err(|e| e.to_string())?;
        let q = StarQuery {
            epoch,
            mag_max: Some(image_mag),
            ..StarQuery::default()
        };
        let stars = cat
            .cone(centre, (0.75 * n * sd).to_radians(), &q)
            .map_err(|e| e.to_string())?;
        let mut sources = Vec::new();
        let mut ids: Vec<Option<u64>> = Vec::new();
        for st in &stars {
            let Some(p) = truth.sky_to_pixel(st.pos) else {
                continue;
            };
            if !(p.x >= 0.0 && p.y >= 0.0 && p.x < n - 1.0 && p.y < n - 1.0) || rng.uniform() < 0.2
            {
                continue;
            }
            sources.push(Source {
                x: p.x + 0.3 * rng.normal(),
                y: p.y + 0.3 * rng.normal(),
                flux: Some(
                    10f64.ppowf(-0.4 * (f64::from(st.mag.unwrap_or(16.0)) - 20.0))
                        * (1.0 + 0.1 * rng.normal()).max(0.1),
                ),
                sigma_px: Some(0.3),
                p_star: None,
                truncated: false,
            });
            ids.push(Some(st.id));
        }
        let n_real = sources.len();
        for k in 0..(n_real * 3 / 10).max(3) {
            let flux = sources
                .get(k % n_real.max(1))
                .and_then(|s| s.flux)
                .unwrap_or(1.0)
                * (0.5 + rng.uniform());
            sources.push(Source {
                x: rng.uniform() * (n - 1.0),
                y: rng.uniform() * (n - 1.0),
                flux: Some(flux),
                sigma_px: Some(0.3),
                p_star: None,
                truncated: false,
            });
            ids.push(None);
        }
        if is_random {
            for s in &mut sources {
                s.x = rng.uniform() * (n - 1.0);
                s.y = rng.uniform() * (n - 1.0);
            }
        }
        let req = SolveRequest {
            frame: Frame {
                width: size,
                height: size,
                trail: None,
            },
            sources,
            priors: Priors::default(),
            policy: Policy {
                epoch,
                deadline,
                pattern_sources,
                ..Policy::default()
            },
        };
        let t0 = Instant::now();
        let out = solve_with_index(&req, &cat, &index).map_err(|e| e.to_string())?;
        let dt = t0.elapsed().as_secs_f64();
        let (status, detail) = match (&out, is_random) {
            (Outcome::Solved(sol), true) => {
                random += 1;
                random_solved += 1;
                (
                    "RANDOM-SOLVED",
                    format!("log-odds {:.1} NFA {:.1e}", sol.log_odds, sol.nfa),
                )
            }
            (
                Outcome::Abstained {
                    hypotheses_tested, ..
                },
                true,
            ) => {
                random += 1;
                hyps.push(*hypotheses_tested);
                (
                    "random-abstained",
                    format!("{hypotheses_tested} hypotheses"),
                )
            }
            (Outcome::Solved(sol), false) => {
                real += 1;
                solved += 1;
                let err = to_arcsec(
                    sol.wcs
                        .crval()
                        .separation(truth.pixel_to_sky(sol.wcs.crpix())),
                );
                hyps.push(sol.hypotheses_tested);
                times.push(dt);
                if err < 5.0 {
                    correct += 1;
                    (
                        "solved",
                        format!("{err:.2}″ {:?} {} hyp", sol.method, sol.hypotheses_tested),
                    )
                } else {
                    wrong += 1;
                    (
                        "WRONG",
                        format!("{err:.1}″ off, log-odds {:.1}", sol.log_odds),
                    )
                }
            }
            (
                Outcome::Abstained {
                    reason,
                    hypotheses_tested,
                    ..
                },
                false,
            ) => {
                real += 1;
                if n_real < 8 {
                    too_few += 1;
                }
                hyps.push(*hypotheses_tested);
                (
                    "abstained",
                    format!(
                        "{reason}; {n_real} real sources, {hypotheses_tested} hyp; {}",
                        diagnose(&req, &ids, &index)
                    ),
                )
            }
        };
        let (ra, dec) = centre.to_degrees();
        eprintln!("frame {f:4} ({ra:7.2}, {dec:+6.2}) {status:<16} {dt:7.3}s  {detail}");
    }
    times.sort_by(f64::total_cmp);
    let pct = |p: f64| {
        times
            .get(((times.len() as f64 - 1.0) * p).round() as usize)
            .copied()
    };
    let summary = serde_json::json!({
        "frames_real": real,
        "solved": solved,
        "solved_correct": correct,
        "solved_wrong": wrong,
        "abstained": real - solved,
        "abstained_with_fewer_than_8_real_sources": too_few,
        "frames_random": random,
        "random_solved": random_solved,
        "solve_seconds_p50": pct(0.5),
        "solve_seconds_p90": pct(0.9),
        "solve_seconds_max": pct(1.0),
        "max_hypotheses": hyps.iter().max(),
        "index_patterns": index.meta().pattern_count,
        "index_stars": index.meta().star_count,
        "frame_px": size,
        "scale_arcsec_px": scale,
        "image_mag_limit": image_mag,
    });
    crate::out(&serde_json::to_string_pretty(&summary).expect("json"))?;
    let _ = SkyPos::from_degrees(0.0, 0.0);
    Ok(())
}

fn diagnose(req: &SolveRequest, ids: &[Option<u64>], index: &Index) -> String {
    use fayde_index::pattern::{angle, shape_of};
    let by_id: std::collections::HashMap<u64, u32> = index
        .stars()
        .iter()
        .enumerate()
        .map(|(i, s)| (s.id, i as u32))
        .collect();
    let mut order: Vec<usize> = (0..req.sources.len()).collect();
    order.sort_by(|&a, &b| {
        req.sources[b]
            .flux
            .unwrap_or(0.0)
            .total_cmp(&req.sources[a].flux.unwrap_or(0.0))
            .then(a.cmp(&b))
    });
    order.truncate(req.policy.pattern_sources);
    let real = order.iter().filter(|&&i| ids[i].is_some()).count();
    let in_index: Vec<u32> = order
        .iter()
        .filter_map(|&i| ids[i].and_then(|id| by_id.get(&id).copied()))
        .collect();
    let n_index = in_index.len();
    let in_index: Vec<u32> = in_index
        .into_iter()
        .filter(|&s| index.stars()[s as usize].pattern)
        .collect();
    let (mut present, mut total) = (0, 0);
    let k = in_index.len();
    for a in 0..k {
        for b in a + 1..k {
            for c in b + 1..k {
                for d in c + 1..k {
                    let set = [in_index[a], in_index[b], in_index[c], in_index[d]];
                    let pts = set.map(|s| index.stars()[s as usize].vec);
                    let (shape, _) = shape_of(&pts, angle);
                    total += 1;
                    let mut sorted = set;
                    sorted.sort_unstable();
                    if index.lookup(&shape, 1e-6, None).iter().any(|m| {
                        let mut s = m.stars;
                        s.sort_unstable();
                        s == sorted
                    }) {
                        present += 1;
                    }
                }
            }
        }
    }
    format!(
        "diag: top {} = {real} real, {n_index} index stars, {k} pattern stars; true patterns {present}/{total} in table",
        order.len()
    )
}
