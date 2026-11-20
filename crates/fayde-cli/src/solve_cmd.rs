use std::collections::HashMap;

use fayde_api::{
    build_request, AtCheck, AtPosition, ImageRequestError, Operation, PreparedImage,
    RequestOptions, Solver,
};
use fayde_core::{PixelConvention, PixelPos, SkyPos};
use fayde_fits::{comment_cards, read_sources, write_wcs_file, FitsFile, InstrumentProfile};
use fayde_solve::{Solution, SolveRequest};
use fayde_wcs::Wcs;

use crate::{index_paths, num, options, seconds};

const KNOWN: &[&str] = &[
    "sources",
    "xyls",
    "ucac4",
    "index",
    "header",
    "instrument",
    "wcs",
    "width",
    "height",
    "ra",
    "dec",
    "radius",
    "pointing-sigma",
    "scale",
    "scale-sigma",
    "roll",
    "roll-sigma",
    "roll-max",
    "parity",
    "trail-length",
    "trail-angle",
    "epoch",
    "deadline",
    "max-sip-order",
    "pixel-error",
    "accept-log-odds",
    "max-hypotheses",
    "image",
    "fwhm",
    "detect-threshold",
    "saturation",
    "gain",
    "max-sources",
    "at",
    "at-sky",
    "at-sources",
    "at-xyls",
    "wcs-sigma",
    "search-radius",
    "max-stars",
    "o",
    "output",
];
const FLAGS: &[&str] = &[
    "no-priors",
    "search-outside-pointing",
    "relax-priors-on-failure",
    "no-pointing",
    "no-header-wcs",
    "no-trails",
];

#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum Mode {
    Solve,
    Verify,
    Refit,
}

pub(crate) fn run(args: &[&str], mode: Mode) -> Result<(), String> {
    let mut o = options(args, KNOWN, FLAGS)?;
    canonical_names(&mut o)?;
    let at_mode = ["at", "at-sky", "at-sources"]
        .iter()
        .any(|k| o.contains_key(*k));
    if at_mode {
        if mode != Mode::Verify {
            return Err("--at, --at-sky and --at-sources go with fayde verify".into());
        }
        unused_options(&o, mode);
        return run_at(&o);
    }
    unused_options(&o, mode);
    if o.contains_key("o") && o.contains_key("output") {
        return Err("-o and --output are the same option: give one".into());
    }
    match (o.contains_key("sources"), o.contains_key("image")) {
        (true, true) => return Err("give --sources or --image, not both".into()),
        (false, false) => {
            return Err("--sources (a source list) or --image (a FITS image) is required".into())
        }
        _ => {}
    }
    let Some(path) = output_path(&o, mode) else {
        return solve_once(&o, mode, None);
    };
    crate::check_output_path(path)?;
    for key in ["wcs", "sources", "image", "header", "instrument", "ucac4"] {
        if let Some(input) = o.get(key) {
            crate::check_output_against_input(path, input, &format!("the --{key} input"))?;
        }
    }
    for input in index_paths(&o) {
        crate::check_output_against_input(path, &input, "an --index input")?;
    }
    let result = solve_once(&o, mode, Some(path));
    if result.is_err() {
        match std::fs::remove_file(path) {
            Ok(()) => eprintln!("fayde: failed: removed the existing {path}"),
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
            Err(e) => eprintln!("fayde: failed, and could not remove the existing {path}: {e}"),
        }
    }
    result
}

fn output_path(o: &HashMap<String, String>, mode: Mode) -> Option<&str> {
    match mode {
        Mode::Verify => None,
        Mode::Solve | Mode::Refit => o.get("o").or_else(|| o.get("output")).map(String::as_str),
    }
}

fn solve_once(
    o: &HashMap<String, String>,
    mode: Mode,
    out_path: Option<&str>,
) -> Result<(), String> {
    let mut opts = request_options(o)?;
    let start = |warnings: &[String]| -> Result<(Solver, Operation), String> {
        for w in warnings {
            eprintln!("fayde: warning: {w}");
        }
        let indexes = match mode {
            Mode::Solve => index_paths(o),
            Mode::Verify | Mode::Refit => Vec::new(),
        };
        let solver = Solver::open(&crate::ucac4_dir(o)?, &indexes)?;
        let operation = match mode {
            Mode::Solve => Operation::Solve,
            Mode::Verify | Mode::Refit => {
                let path = o.get("wcs").ok_or("--wcs is required")?;
                let f = FitsFile::open(path).map_err(|e| format!("{path}: {e}"))?;
                let wcs =
                    Wcs::from_header(&f.primary().header).map_err(|e| format!("{path}: {e}"))?;
                if mode == Mode::Verify {
                    Operation::Verify(wcs)
                } else {
                    Operation::Refit(wcs)
                }
            }
        };
        Ok((solver, operation))
    };
    let (built, report) = if let Some(path) = o.get("image") {
        let image = fayde_fits::read_image(path).map_err(|e| wrong_kind(path, "--image", e))?;
        let mut detect = detect_options(o)?;
        detect.points_only = o.contains_key("no-trails");
        let input_error = |e| match e {
            ImageRequestError::Detect(e) => format!("{path}: {e}"),
            ImageRequestError::Request(e) => e,
        };
        let prepared = PreparedImage::new(image, opts, detect).map_err(input_error)?;
        let (solver, operation) = start(prepared.warnings())?;
        prepared.run(&solver, &operation).map_err(input_error)?
    } else {
        for k in DETECT_OPTIONS.iter().filter(|k| o.contains_key(**k)) {
            eprintln!("fayde: warning: --{k} applies to --image: ignored");
        }
        let path = &o["sources"];
        let list = read_sources(path).map_err(|e| wrong_kind(path, "--sources", e))?;
        fayde_api::use_trail(&mut opts, list.trail);
        let built = build_request(list.sources, (list.width, list.height), &opts)?;
        let (solver, operation) = start(&built.warnings)?;
        let report = solver.run(&built, &operation);
        (built, report)
    };
    crate::out(&serde_json::to_string_pretty(&report.json).expect("json"))?;
    if let Some(path) = out_path {
        match &report.solution {
            Some(sol) => write_solution(path, sol, &built.request)?,
            None => match std::fs::remove_file(path) {
                Ok(()) => eprintln!("fayde: no solution: removed the existing {path}"),
                Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
                Err(e) => return Err(format!("{path}: {e}")),
            },
        }
    }
    if report.error {
        return Err(report.json["reason"]
            .as_str()
            .unwrap_or("error")
            .to_string());
    }
    Ok(())
}

fn run_at(o: &HashMap<String, String>) -> Result<(), String> {
    let pair = |v: &str, what: &str| -> Result<(f64, f64), String> {
        let p: Vec<f64> = v
            .split(',')
            .map(|s| s.trim().parse::<f64>())
            .collect::<Result<_, _>>()
            .map_err(|_| format!("{what}: cannot parse '{v}'"))?;
        match p[..] {
            [a, b] if a.is_finite() && b.is_finite() => Ok((a, b)),
            _ => Err(format!("{what} needs two finite numbers, got '{v}'")),
        }
    };
    let mut positions = Vec::new();
    for v in o
        .get("at")
        .into_iter()
        .flat_map(|v| v.split(crate::LIST_SEP))
    {
        let (x, y) = pair(v, "--at")?;
        let p = PixelConvention::Fits.to_internal(x, y);
        positions.push(AtPosition::Pixel(p, None));
    }
    for v in o
        .get("at-sky")
        .into_iter()
        .flat_map(|v| v.split(crate::LIST_SEP))
    {
        let (ra, dec) = pair(v, "--at-sky")?;
        let sky = SkyPos::from_degrees(ra, dec).ok_or_else(|| format!("--at-sky: bad '{v}'"))?;
        positions.push(AtPosition::Sky(sky, None));
    }
    for path in o
        .get("at-sources")
        .into_iter()
        .flat_map(|v| v.split(crate::LIST_SEP))
    {
        let xy = read_sources(path).map_err(|e| wrong_kind(path, "--at-sources", e))?;
        for s in xy.sources {
            positions.push(AtPosition::Pixel(PixelPos::new(s.x, s.y), s.sigma_px));
        }
    }
    let opts = request_options(o)?;
    let solver = Solver::open(&crate::ucac4_dir(o)?, &[])?;
    let path = o.get("wcs").ok_or("--wcs is required")?;
    let f = FitsFile::open(path).map_err(|e| format!("{path}: {e}"))?;
    let wcs = Wcs::from_header(&f.primary().header).map_err(|e| format!("{path}: {e}"))?;
    let at = AtCheck {
        wcs_sigma_arcsec: num(o, "wcs-sigma")?,
        search_radius_arcsec: num(o, "search-radius")?,
        max_stars: num(o, "max-stars")?,
    };
    let report = solver.check_at(&wcs, &positions, &opts, &at);
    crate::out(&serde_json::to_string_pretty(&report.json).expect("json"))?;
    if report.error {
        return Err(report.json["reason"]
            .as_str()
            .unwrap_or("error")
            .to_string());
    }
    Ok(())
}

pub(crate) fn wrong_kind(path: &str, option: &str, e: fayde_fits::FitsError) -> String {
    use fayde_fits::FileKind;
    let kind = fayde_fits::FitsFile::open(path)
        .map(|f| fayde_fits::file_kind(&f))
        .ok();
    match (option, kind) {
        ("--sources" | "--at-sources", Some(FileKind::Image)) => {
            format!("{path} is an image, not a source list: give it with --image to find its stars")
        }
        ("--image", Some(FileKind::SourceList)) => {
            format!("{path} is a source list, not an image: give it with --sources")
        }
        _ => format!("{path}: {e}"),
    }
}

fn canonical_names(o: &mut HashMap<String, String>) -> Result<(), String> {
    for (old, new) in [("xyls", "sources"), ("at-xyls", "at-sources")] {
        if let Some(v) = o.remove(old) {
            match o.get_mut(new) {
                Some(cur) if new == "at-sources" => {
                    cur.push(crate::LIST_SEP);
                    cur.push_str(&v);
                }
                Some(_) => return Err(format!("--{old} is another name for --{new}: give one")),
                None => {
                    o.insert(new.to_string(), v);
                }
            }
        }
    }
    Ok(())
}

pub(crate) const DETECT_OPTIONS: &[&str] = &[
    "fwhm",
    "detect-threshold",
    "saturation",
    "gain",
    "max-sources",
];

pub(crate) fn detect_options(
    o: &HashMap<String, String>,
) -> Result<fayde_api::DetectOptions, String> {
    Ok(fayde_api::DetectOptions {
        fwhm: num(o, "fwhm")?,
        threshold: num(o, "detect-threshold")?,
        saturation: num(o, "saturation")?,
        gain: num(o, "gain")?,
        max_sources: num(o, "max-sources")?,
        trail_angle: match (
            num::<f64>(o, "trail-length")?,
            num::<f64>(o, "trail-angle")?,
        ) {
            (_, Some(a)) => Some(a),
            (Some(_), None) => return Err("--trail-length needs --trail-angle".into()),
            (None, None) => None,
        },
        points_only: false,
    })
}

fn request_options(o: &HashMap<String, String>) -> Result<RequestOptions, String> {
    let header = o
        .get("header")
        .map(|p| fayde_fits::read_image_header(p).map_err(|e| format!("{p}: {e}")))
        .transpose()?;
    let instrument = match o.get("instrument") {
        Some(p) => Some(InstrumentProfile::load(p).map_err(|e| format!("{p}: {e}"))?),
        None => None,
    };
    Ok(RequestOptions {
        width: num(o, "width")?,
        height: num(o, "height")?,
        ra: num(o, "ra")?,
        dec: num(o, "dec")?,
        radius: num(o, "radius")?,
        pointing_sigma: num(o, "pointing-sigma")?,
        scale: num(o, "scale")?,
        scale_sigma: num(o, "scale-sigma")?,
        roll: num(o, "roll")?,
        roll_sigma: num(o, "roll-sigma")?,
        roll_max: num(o, "roll-max")?,
        parity: o.get("parity").cloned(),
        trail_length: num(o, "trail-length")?,
        trail_angle: num(o, "trail-angle")?,
        epoch: num(o, "epoch")?,
        deadline: seconds(o, "deadline")?,
        max_sip_order: num(o, "max-sip-order")?,
        pixel_error: num(o, "pixel-error")?,
        accept_log_odds: num(o, "accept-log-odds")?,
        max_hypotheses: num(o, "max-hypotheses")?,
        no_priors: o.contains_key("no-priors"),
        no_pointing: o.contains_key("no-pointing"),
        no_header_wcs: o.contains_key("no-header-wcs"),
        search_outside_pointing: o.contains_key("search-outside-pointing"),
        relax_priors_on_failure: o.contains_key("relax-priors-on-failure"),
        header,
        instrument,
    })
}

fn write_solution(path: &str, sol: &Solution, req: &SolveRequest) -> Result<(), String> {
    let mut cards = fayde_api::solution_cards(sol, req);
    cards.extend(comment_cards(&format!(
        "Solved by fayde {} (commit {}, {:?}). log-odds {:.1}, NFA {:.2e}, {} hypotheses, {} matches, \
         rms {:.3} px, leave-one-out rms {}, error scale {:.3}, centre 1-sigma {:.3} arcsec, \
         worst 1-sigma {:.3} arcsec, catalogue {} {}, index {}. Please cite Fayde in \
         publications (Ryall 2026; see CITATION.cff).",
        env!("CARGO_PKG_VERSION"),
        fayde_api::commit(),
        sol.method,
        sol.log_odds,
        sol.nfa,
        sol.hypotheses_tested,
        sol.matches.len(),
        sol.rms_px,
        sol.loo_rms_px
            .map_or("n/a".into(), |v| format!("{v:.3} px")),
        sol.error_scale,
        sol.covariance.radial_sigma_arcsec(sol.wcs.crpix()),
        sol.worst_sigma_arcsec,
        sol.catalog.name,
        sol.catalog.release,
        sol.index_hash
            .map_or("none".into(), |h| format!("{h:016x}")),
    )));
    write_wcs_file(path, &cards).map_err(|e| format!("{path}: {e}"))
}

fn unused_options(o: &HashMap<String, String>, mode: Mode) {
    let (name, unused): (&str, &[&str]) = match mode {
        Mode::Solve => return,
        Mode::Verify => ("verify", &["index", "o", "output"]),
        Mode::Refit => ("refit", &["index"]),
    };
    for k in unused.iter().filter(|k| o.contains_key(**k)) {
        let dashes = if k.len() == 1 { "-" } else { "--" };
        eprintln!("fayde: warning: {dashes}{k} is not used by fayde {name}: ignored");
    }
}
