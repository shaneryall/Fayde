use std::{
    error::Error,
    fmt,
    time::{Duration, Instant},
};

use fayde_detect::FrameMode;
use fayde_fits::Image;
use fayde_solve::{solve_auto_after, Outcome};
use serde_json::{json, Value};

use crate::{
    build_request, detect_image, use_trail, Built, DetectOptions, Operation, Report,
    RequestOptions, Solver,
};

#[derive(Debug)]
pub enum ImageRequestError {
    Detect(String),
    Request(String),
}

impl fmt::Display for ImageRequestError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Detect(message) | Self::Request(message) => f.write_str(message),
        }
    }
}

impl Error for ImageRequestError {}

pub struct PreparedImage {
    image: Image,
    request: RequestOptions,
    detect: DetectOptions,
    first: ImageAttempt,
}

struct ImageAttempt {
    built: Built,
    mode: FrameMode,
    detection: Value,
}

impl PreparedImage {
    pub fn new(
        image: Image,
        mut request: RequestOptions,
        detect: DetectOptions,
    ) -> Result<Self, ImageRequestError> {
        if request.header.is_none() {
            request.header = Some(image.header.clone());
        }
        let first = ImageAttempt::new(&image, &request, &detect)?;
        Ok(Self {
            image,
            request,
            detect,
            first,
        })
    }

    pub fn warnings(&self) -> &[String] {
        &self.first.built.warnings
    }

    pub fn run(
        self,
        solver: &Solver,
        operation: &Operation,
    ) -> Result<(Built, Report), ImageRequestError> {
        let first = self.first;
        if !matches!(operation, Operation::Solve)
            || !FrameMode::POINT_RETRY_MODES.contains(&first.mode)
        {
            let report = solver.run(&first.built, operation);
            return Ok(first.with_report(report));
        }
        let solve = |attempt: &ImageAttempt, tested| {
            solve_auto_after(
                &attempt.built.request,
                solver.catalog(),
                &solver.indexes,
                tested,
            )
        };
        let started = Instant::now();
        let mut outcome = solve(&first, 0);
        let Ok(Outcome::Abstained {
            hypotheses_tested,
            search,
            best_log_odds,
            ..
        }) = &mut outcome
        else {
            return Ok(first.solve_report(solver, outcome));
        };

        let remaining = || {
            self.request
                .deadline
                .map(|limit| limit.saturating_sub(started.elapsed()))
        };
        if *hypotheses_tested >= first.built.request.policy.max_hypotheses
            || remaining() == Some(Duration::ZERO)
        {
            return Ok(first.solve_report(solver, outcome));
        }
        let points = DetectOptions {
            points_only: true,
            ..self.detect
        };
        let mut second = ImageAttempt::new(&self.image, &self.request, &points)?;
        second.built.request.policy.deadline = remaining();
        if second.built.request.policy.deadline == Some(Duration::ZERO) {
            return Ok(first.solve_report(solver, outcome));
        }

        let retry = solve(&second, *hypotheses_tested);
        let note = json!({"solved": matches!(retry, Ok(Outcome::Solved(_))), "trails_census": first.detection["census"]});
        let selected = match retry {
            Ok(Outcome::Solved(mut solution)) => {
                for stage in search.iter_mut() {
                    stage.stage = format!("{} detection: {}", first.mode.name(), stage.stage);
                }
                search.append(&mut solution.search);
                solution.search = std::mem::take(search);
                solution.elapsed = started.elapsed();
                outcome = Ok(Outcome::Solved(solution));
                second
            }
            Ok(Outcome::Abstained {
                hypotheses_tested: total,
                search: stages,
                best_log_odds: best,
                ..
            }) => {
                *hypotheses_tested = total;
                *best_log_odds = best
                    .into_iter()
                    .chain(*best_log_odds)
                    .max_by(f64::total_cmp);
                search.extend(stages.into_iter().map(|mut stage| {
                    stage.stage = format!("point retry: {}", stage.stage);
                    stage
                }));
                first
            }
            Err(error) => {
                let (built, mut report) = second.solve_report(solver, Err(error));
                report.json["retried_as_points"] = note;
                report.json["retried_as_points"]["first_search"] = crate::trace_json(search);
                report.json["retried_as_points"]["first_hypotheses_tested"] =
                    json!(*hypotheses_tested);
                return Ok((built, report));
            }
        };
        let (built, mut report) = selected.solve_report(solver, outcome);
        report.json["retried_as_points"] = note;
        Ok((built, report))
    }
}

impl ImageAttempt {
    fn solve_report(
        self,
        solver: &Solver,
        outcome: Result<Outcome, fayde_solve::SolveError>,
    ) -> (Built, Report) {
        let report = solver.outcome(outcome, &self.built, true);
        self.with_report(report)
    }

    fn with_report(self, mut report: Report) -> (Built, Report) {
        report.json["detection"] = self.detection;
        (self.built, report)
    }

    fn new(
        image: &Image,
        request: &RequestOptions,
        detect: &DetectOptions,
    ) -> Result<Self, ImageRequestError> {
        let detection = detect_image(image, detect, crate::DetectionReport::Summary)
            .map_err(ImageRequestError::Detect)?;
        let mut opts = request.clone();
        use_trail(&mut opts, detection.trail);
        let size = (
            u32::try_from(image.width).ok(),
            u32::try_from(image.height).ok(),
        );
        let built =
            build_request(detection.sources, size, &opts).map_err(ImageRequestError::Request)?;
        Ok(Self {
            built,
            mode: detection.mode,
            detection: detection.report,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fayde_detect::synth::{Scene, Star};

    #[test]
    fn point_retry_catalogue_failure_reaches_the_report() {
        let dir = std::env::temp_dir().join(format!("fayde-retry-error-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join("z001"), []).unwrap();
        let solver = Solver::open(dir.to_str().unwrap(), &[]).unwrap();
        let mut scene = Scene::new(128, 128);
        for (x, y) in [
            (16., 25.),
            (90., 40.),
            (55., 70.),
            (32., 108.),
            (100., 104.),
        ] {
            scene.stars.push(Star { x, y, flux: 20000. });
        }
        let image = Image {
            width: 128,
            height: 128,
            header: Default::default(),
            data: scene.render(42),
            hdu: 0,
        };
        let request = RequestOptions {
            ra: Some(150.),
            dec: Some(2.),
            scale: Some(2.),
            radius: Some(0.1),
            ..Default::default()
        };
        let mut prep = PreparedImage::new(image, request, DetectOptions::default()).unwrap();
        assert!(prep.first.built.request.sources.len() >= 4);
        prep.first.mode = FrameMode::Trails;
        prep.first.built.request.sources.truncate(3);
        let (_, report) = prep.run(&solver, &Operation::Solve).unwrap();
        assert!(report.error, "{}", report.json);
        assert!(report.json["reason"].as_str().unwrap().contains("zone"));
        assert_eq!(report.json["detection"]["mode"], "points");
        assert_eq!(
            report.json["retried_as_points"]["first_hypotheses_tested"],
            0
        );
        std::fs::remove_dir_all(dir).unwrap();
    }
}
