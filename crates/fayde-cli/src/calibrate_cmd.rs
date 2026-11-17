//! `fayde calibrate`: an instrument profile from frames Fayde has solved.
//!
//! ```text
//! fayde calibrate --solved <list> [--roll-keyword <KEY>] [--binning <XKEY>[,<YKEY>]]
//!                 [--name <text>] [--no-header-wcs] [-o <profile.json>]
//! ```
//!
//! `<list>` is a text file with one solved frame per line: the frame's FITS file (for its
//! header), a tab, and the WCS file of its solution (`fayde solve -o`). Blank lines and lines
//! starting with `#` are skipped. The profile (JSON) goes to stdout, and to `-o` when given.
//!
//! A profile needs no other solver: solve some frames with their pointing and a nominal scale
//! (`--scale`), then calibrate from those solutions.

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
        let header = FitsFile::open(frame)
            .map_err(|e| at(&format!("{frame}: {e}")))?
            .primary()
            .header
            .clone();
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
        std::fs::write(path, format!("{text}\n")).map_err(|e| format!("{path}: {e}"))?;
    }
    Ok(())
}
