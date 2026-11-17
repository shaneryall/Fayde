//! `fayde solve`, `fayde verify`, `fayde refit`: one source list at a time.
//!
//! ```text
//! fayde solve  (--sources <frame.sources.fits> | --image <frame.fits>) --ucac4 <dir> [request options]
//!              [--index <a.fidx>]... [-o <out.wcs>]
//! fayde verify (--sources <frame.sources.fits> | --image <frame.fits>) --ucac4 <dir> --wcs <in.wcs>
//! fayde refit  (--sources <frame.sources.fits> | --image <frame.fits>) --ucac4 <dir> --wcs <in.wcs>
//!              [-o <out.wcs>]
//!
//! `--xyls` and `--at-xyls` are accepted as earlier names of `--sources` and `--at-sources`.
//!
//! request options:
//!   --header <frame.fits> --instrument <profile.json> --width <px> --height <px>
//!   --ra <deg> --dec <deg> --radius <deg> --pointing-sigma <deg>
//!   --scale <arcsec/px> --scale-sigma <fraction>
//!   --roll <deg> --roll-sigma <deg> --roll-max <deg> --parity negative|positive
//!   --trail-length <px> --trail-angle <deg> --epoch <year> --deadline <s> --max-sip-order <n>
//!   --pixel-error <px> --accept-log-odds <x> --max-hypotheses <n>
//!   --no-priors --no-pointing --no-header-wcs --search-outside-pointing
//!   --relax-priors-on-failure
//!
//! `--index` may be repeated (one path each); a comma-separated list still works for paths
//! without commas. `--no-priors` drops every prior; `--no-pointing` only the pointing.
//! ```
//!
//! Solved and verify reports list the matched pairs (`pairs`): source index in the input
//! table, source position (FITS 1-based, as in the source list), catalogue star id (UCAC4: the
//! running number RNM), position at the solve epoch and magnitude, predicted pixel (FITS
//! 1-based), residual (source − predicted, px) and inlier probability (the match weight).
//!
//! Priors are taken from, in order of precedence: the command line, the instrument profile,
//! the frame header. The report (JSON) goes to stdout; a WCS, with provenance comments, to `-o`.

use std::collections::HashMap;

use fayde_api::{build_request, AtCheck, AtPosition, RequestOptions, Solver};
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
    // A command-line error, like every other one, leaves an existing -o file alone: checked
    // here, before the part of the run whose failure removes it.
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
    // Found now, not after the solve: an unwritable output must not cost the solution.
    crate::check_output_path(path)?;
    // The output is removed when the run fails (below), so it must not be one of the inputs:
    // `refit --wcs a.wcs -o a.wcs` would otherwise delete the WCS it was given.
    for key in ["wcs", "sources", "image", "header", "instrument", "ucac4"] {
        if let Some(input) = o.get(key) {
            if same_file(path, input) {
                return Err(format!("-o {path} is the --{key} input: write elsewhere"));
            }
        }
    }
    for input in index_paths(&o) {
        if same_file(path, &input) {
            return Err(format!("-o {path} is an --index input: write elsewhere"));
        }
    }
    // Any failure from here on (unreadable input, catalogue, solver error, no solution, failed
    // write) leaves no file at -o: one left by an earlier run would pass for this frame's WCS.
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

/// Where this run writes its WCS, if anywhere. `fayde verify` writes no file: its `-o` is named
/// as unused and the path is never touched, neither checked nor cleared.
fn output_path(o: &HashMap<String, String>, mode: Mode) -> Option<&str> {
    match mode {
        Mode::Verify => None,
        Mode::Solve | Mode::Refit => o.get("o").or_else(|| o.get("output")).map(String::as_str),
    }
}

/// Whether `a` and `b` name the same existing file (after resolving links and `..`).
fn same_file(a: &str, b: &str) -> bool {
    match (std::fs::canonicalize(a), std::fs::canonicalize(b)) {
        (Ok(x), Ok(y)) => x == y,
        _ => false,
    }
}

/// One solve, verify or refit; `out_path` has been checked (writable, not an input).
fn solve_once(
    o: &HashMap<String, String>,
    mode: Mode,
    out_path: Option<&str>,
) -> Result<(), String> {
    let (sources, size, image, detection) = load_sources(o)?;
    let mut opts = request_options(o)?;
    // An image carries its own header: used for the priors unless --header is given.
    if opts.header.is_none() {
        opts.header = image.map(|i| i.header);
    }
    let built = build_request(sources, size, &opts)?;
    for w in &built.warnings {
        eprintln!("fayde: warning: {w}");
    }
    let indexes = match mode {
        Mode::Solve => index_paths(o),
        Mode::Verify | Mode::Refit => Vec::new(),
    };
    let solver = Solver::open(o.get("ucac4").ok_or("--ucac4 is required")?, &indexes)?;
    let given = match mode {
        Mode::Solve => None,
        Mode::Verify | Mode::Refit => {
            let path = o.get("wcs").ok_or("--wcs is required")?;
            let f = FitsFile::open(path).map_err(|e| format!("{path}: {e}"))?;
            Some(Wcs::from_header(&f.primary().header).map_err(|e| format!("{path}: {e}"))?)
        }
    };
    let report = match (mode, &given) {
        (Mode::Verify, Some(wcs)) => solver.verify(&built, wcs),
        (Mode::Refit, Some(wcs)) => solver.refit(&built, wcs),
        _ => solver.solve(&built),
    };
    let mut report = report;
    if let Some(d) = detection {
        report.json["detection"] = d;
    }
    // Every outcome, a solver error included, is a JSON report on stdout, so callers parse
    // them the same way; an error also exits non-zero. The report goes out before the output
    // file is touched: whatever happens to the file, the caller has the result.
    crate::out(&serde_json::to_string_pretty(&report.json).expect("json"))?;
    if let Some(path) = out_path {
        match &report.solution {
            Some(sol) => write_solution(path, sol, &built.request)?,
            // No solution: a file left there by an earlier run would pass for this frame's WCS.
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

/// `fayde verify --at`: is there a catalogue star at these positions under the given WCS?
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
    let solver = Solver::open(o.get("ucac4").ok_or("--ucac4 is required")?, &[])?;
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

/// The request options from the command line; the header and profile are read from their files.
/// The sources: from `--sources` (a source list), or detected in `--image` (exactly one of
/// them). With an image, also the image (for its header) and the detection summary.
#[allow(clippy::type_complexity)]
fn load_sources(
    o: &HashMap<String, String>,
) -> Result<
    (
        Vec<fayde_solve::Source>,
        (Option<u32>, Option<u32>),
        Option<fayde_fits::Image>,
        Option<serde_json::Value>,
    ),
    String,
> {
    match (o.get("sources"), o.get("image")) {
        (Some(_), Some(_)) => Err("give --sources or --image, not both".into()),
        (None, None) => {
            Err("--sources (a source list) or --image (a FITS image) is required".into())
        }
        (Some(path), None) => {
            for k in DETECT_OPTIONS.iter().filter(|k| o.contains_key(**k)) {
                eprintln!("fayde: warning: --{k} applies to --image: ignored");
            }
            let list = read_sources(path).map_err(|e| wrong_kind(path, "--sources", e))?;
            Ok((list.sources, (list.width, list.height), None, None))
        }
        (None, Some(path)) => {
            let image = fayde_fits::read_image(path).map_err(|e| wrong_kind(path, "--image", e))?;
            let (sources, mut report) = fayde_api::detect_sources(&image, &detect_options(o)?)
                .map_err(|e| format!("{path}: {e}"))?;
            // The solve report carries the summary; the per-source list is `fayde detect`'s.
            if let Some(m) = report.as_object_mut() {
                m.remove("sources");
                m.remove("rejected");
            }
            let size = (
                u32::try_from(image.width).ok(),
                u32::try_from(image.height).ok(),
            );
            Ok((sources, size, Some(image), Some(report)))
        }
    }
}

/// The error for a file given to `option` that could not be read as such. When the file holds
/// the other kind (an image given as a source list, or the reverse) the message says so.
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

/// `--xyls` and `--at-xyls` are the earlier names of `--sources` and `--at-sources`, still
/// accepted.
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

/// Options of source detection (`--image`).
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
    })
}

fn request_options(o: &HashMap<String, String>) -> Result<RequestOptions, String> {
    let header = match o.get("header") {
        Some(p) => Some(
            FitsFile::open(p)
                .map_err(|e| format!("{p}: {e}"))?
                .primary()
                .header
                .clone(),
        ),
        None => None,
    };
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

/// An option the command does not use is named on stderr rather than dropped in silence.
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
