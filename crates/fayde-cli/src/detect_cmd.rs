//! `fayde detect`: find stars in a FITS image.
//!
//! ```text
//! fayde detect --image <frame.fits> [-o <frame.sources.fits>] [--fwhm <px>] [--detect-threshold <σ>]
//!              [--saturation <adu>] [--gain <e/adu>] [--max-sources <n>]
//! ```
//!
//! The report (JSON) goes to stdout: the settings used, background and noise, every source
//! (FITS 1-based pixels, flux, error, star probability, class), the streaks, and what was set
//! aside. `-o` also writes the sources as a source list that `fayde solve --sources` reads.

use std::collections::HashMap;

use crate::options;
use crate::solve_cmd::{detect_options, DETECT_OPTIONS};

pub(crate) fn run(args: &[&str]) -> Result<(), String> {
    let mut known: Vec<&str> = vec!["image", "o", "output"];
    known.extend_from_slice(DETECT_OPTIONS);
    let o: HashMap<String, String> = options(args, &known, &[])?;
    let image_path = o.get("image").ok_or("--image is required")?;
    if o.contains_key("o") && o.contains_key("output") {
        return Err("-o and --output are the same option: give one".into());
    }
    let out = o.get("o").or_else(|| o.get("output"));
    if let Some(path) = out {
        crate::check_output_path(path)?;
        if let (Ok(a), Ok(b)) = (
            std::fs::canonicalize(path),
            std::fs::canonicalize(image_path),
        ) {
            if a == b {
                return Err(format!("-o {path} is the --image input: write elsewhere"));
            }
        }
    }
    let image = fayde_fits::read_image(image_path)
        .map_err(|e| crate::solve_cmd::wrong_kind(image_path, "--image", e))?;
    let (sources, mut report) = fayde_api::detect_sources(&image, &detect_options(&o)?)
        .map_err(|e| format!("{image_path}: {e}"))?;
    report["identity"] = serde_json::json!({
        "fayde": {"version": env!("CARGO_PKG_VERSION"), "commit": fayde_api::commit()},
    });
    crate::out(&serde_json::to_string_pretty(&report).expect("json"))?;
    if let Some(path) = out {
        let (w, h) = (
            u32::try_from(image.width).map_err(|_| "image too wide".to_string())?,
            u32::try_from(image.height).map_err(|_| "image too tall".to_string())?,
        );
        fayde_fits::write_sources(path, &sources, w, h).map_err(|e| format!("{path}: {e}"))?;
    }
    Ok(())
}
