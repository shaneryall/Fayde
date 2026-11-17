use fayde_api::{calibrate_charge_transfer, calibrate_profile, CalibrateOptions, SolvedFrame};
use fayde_core::math::round_decimal;
use fayde_fits::{header_priors, FitsFile};
use fayde_wcs::Wcs;

use crate::options;
use crate::solve_cmd::{detect_options, DETECT_OPTIONS};

pub(crate) fn run(args: &[&str]) -> Result<(), String> {
    let mut known = vec![
        "solved",
        "roll-keyword",
        "binning",
        "name",
        "o",
        "output",
        "ucac4",
        "epoch",
    ];
    known.extend_from_slice(DETECT_OPTIONS);
    let o = options(
        args,
        &known,
        &["no-header-wcs", "charge-transfer", "periodic-pattern"],
    )?;
    let transfer = o.contains_key("charge-transfer");
    let pattern = o.contains_key("periodic-pattern");
    let mut detect = detect_options(&o, None)?;
    if o.contains_key("o") && o.contains_key("output") {
        return Err("-o and --output are the same option: give one".into());
    }
    let list_path = o.get("solved").ok_or("--solved is required")?;
    let out_path = o.get("o").or_else(|| o.get("output"));
    if let Some(path) = out_path {
        crate::check_output_path(path)?;
        crate::check_output_against_input(path, list_path, "the --solved input")?;
    }
    let solver = if transfer {
        Some(fayde_api::Solver::open(&crate::ucac4_dir(&o)?, &[])?)
    } else {
        None
    };
    let binning = o.get("binning").map(|v| {
        let mut keys = v.split(',').map(|k| k.trim().to_ascii_uppercase());
        (keys.next().unwrap_or_default(), keys.next())
    });
    if binning.as_ref().is_some_and(|(x, _)| x.is_empty()) {
        return Err("--binning needs a keyword, e.g. XBINNING,YBINNING".into());
    }
    let list = std::fs::read_to_string(list_path).map_err(|e| format!("{list_path}: {e}"))?;
    let mut frames = Vec::new();
    let mut images: Vec<(String, usize)> = Vec::new();
    for (n, line) in list.lines().enumerate() {
        let line = line.trim_end_matches('\r');
        if line.trim().is_empty() || line.trim_start().starts_with('#') {
            continue;
        }
        let at = |what: &str| format!("{list_path}, line {}: {what}", n + 1);
        let (frame, wcs) = line
            .split_once('\t')
            .map(|(a, b)| (a.trim(), b.trim()))
            .filter(|(a, b)| !a.is_empty() && !b.is_empty())
            .ok_or_else(|| at("expected the frame's FITS file, a tab, and its WCS file"))?;
        if let Some(path) = out_path {
            crate::check_output_against_input(path, frame, "a frame in the --solved list")
                .map_err(|e| at(&e))?;
            crate::check_output_against_input(path, wcs, "a WCS in the --solved list")
                .map_err(|e| at(&e))?;
        }
        let header =
            fayde_fits::read_image_header(frame).map_err(|e| at(&format!("{frame}: {e}")))?;
        let (width, height) = header_priors(&header).size.ok_or_else(|| {
            at(&format!(
                "{frame}: no image size (NAXIS1, NAXIS2) in its header"
            ))
        })?;
        let solution = FitsFile::open(wcs).map_err(|e| at(&format!("{wcs}: {e}")))?;
        let wcs =
            Wcs::from_header(&solution.primary().header).map_err(|e| at(&format!("{wcs}: {e}")))?;
        images.push((frame.to_string(), n + 1));
        frames.push(SolvedFrame {
            header,
            wcs,
            width,
            height,
        });
    }
    let image = |(frame, line): &(String, usize)| {
        fayde_fits::read_image(frame)
            .map_err(|e| crate::solve_cmd::wrong_kind(frame, "--image", e))
            .map_err(|e| format!("{list_path}, line {line}: {e}"))
    };
    let mut pattern_fit = None;
    if pattern {
        let mut finder = fayde_detect::pattern::Finder::default();
        for entry in &images {
            let im = image(entry)?;
            finder.add(&im.data, im.width, im.height);
        }
        let (found, used) = finder.finish();
        let freqs: Vec<[f64; 2]> = found
            .iter()
            .map(|f| f.frequency.map(|v| round_decimal(v, 5)))
            .collect();
        detect.periodic_pattern = if freqs.is_empty() {
            None
        } else {
            Some(fayde_fits::PeriodicPattern::new(&freqs)?)
        };
        pattern_fit = Some((freqs, found, used));
    }
    let mut transfer_frames = Vec::new();
    if let Some(solver) = &solver {
        for (entry, f) in images.iter().zip(&frames) {
            let im = image(entry)?;
            transfer_frames.push(
                fayde_api::transfer_frame(solver, &im, &f.wcs, &detect)
                    .map_err(|e| format!("{list_path}, line {}: {e}", entry.1))?,
            );
        }
    }
    let mut profile = calibrate_profile(
        &frames,
        &CalibrateOptions {
            roll_keyword: o.get("roll-keyword").map(|k| k.trim().to_ascii_uppercase()),
            binning,
            name: o.get("name").cloned(),
            no_header_wcs: o.contains_key("no-header-wcs"),
        },
    )?;
    if let Some((freqs, found, used)) = pattern_fit {
        if freqs.is_empty() {
            eprintln!("fayde: warning: no periodic pattern found in the frames' backgrounds");
        } else {
            profile["periodic_pattern"] = serde_json::json!({ "frequencies": freqs });
        }
        profile["calibration"]["periodic_pattern"] = serde_json::json!({
            "frames": used,
            "peaks": found.iter().map(|f| serde_json::json!({
                "frequency": f.frequency,
                "variance_share": f.share,
            })).collect::<Vec<_>>(),
        });
    }
    if transfer {
        let (term, stats) = calibrate_charge_transfer(&transfer_frames)?;
        match &term {
            Some(t) => profile["charge_transfer"] = t.clone(),
            None => eprintln!(
                "fayde: warning: no charge-transfer term: no geometry's shift is positive and \
                 5σ from zero (see calibration.charge_transfer)"
            ),
        }
        profile["calibration"]["charge_transfer"] = stats;
    }
    let text = serde_json::to_string_pretty(&profile).expect("json");
    crate::out(&text)?;
    if let Some(path) = out_path {
        fayde_core::fsutil::write_atomic(path.as_ref(), format!("{text}\n").as_bytes())
            .map_err(|e| format!("{path}: {e}"))?;
    }
    Ok(())
}
