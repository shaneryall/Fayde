use fayde_api::{calibrate_profile, CalibrateOptions, SolvedFrame};
use fayde_fits::{header_priors, FitsFile};
use fayde_wcs::Wcs;

use crate::options;

pub(crate) fn run(args: &[&str]) -> Result<(), String> {
    let o = options(
        args,
        &["solved", "roll-keyword", "binning", "name", "o", "output"],
        &["no-header-wcs"],
    )?;
    if o.contains_key("o") && o.contains_key("output") {
        return Err("-o and --output are the same option: give one".into());
    }
    let list_path = o.get("solved").ok_or("--solved is required")?;
    let out_path = o.get("o").or_else(|| o.get("output"));
    if let Some(path) = out_path {
        crate::check_output_path(path)?;
        crate::check_output_against_input(path, list_path, "the --solved input")?;
    }
    let binning = o.get("binning").map(|v| {
        let mut keys = v.split(',').map(|k| k.trim().to_ascii_uppercase());
        (keys.next().unwrap_or_default(), keys.next())
    });
    if binning.as_ref().is_some_and(|(x, _)| x.is_empty()) {
        return Err("--binning needs a keyword, e.g. XBINNING,YBINNING".into());
    }
    let list = std::fs::read_to_string(list_path).map_err(|e| format!("{list_path}: {e}"))?;
    let mut frames = Vec::new();
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
        frames.push(SolvedFrame {
            header,
            wcs,
            width,
            height,
        });
    }
    let profile = calibrate_profile(
        &frames,
        &CalibrateOptions {
            roll_keyword: o.get("roll-keyword").map(|k| k.trim().to_ascii_uppercase()),
            binning,
            name: o.get("name").cloned(),
            no_header_wcs: o.contains_key("no-header-wcs"),
        },
    )?;
    let text = serde_json::to_string_pretty(&profile).expect("json");
    crate::out(&text)?;
    if let Some(path) = out_path {
        fayde_core::fsutil::write_atomic(path.as_ref(), format!("{text}\n").as_bytes())
            .map_err(|e| format!("{path}: {e}"))?;
    }
    Ok(())
}
