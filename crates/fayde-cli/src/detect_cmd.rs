use std::collections::HashMap;

use crate::options;
use crate::solve_cmd::{detect_options, DETECT_OPTIONS};

pub(crate) fn run(args: &[&str]) -> Result<(), String> {
    let mut known: Vec<&str> = vec![
        "image",
        "o",
        "output",
        "trail-length",
        "trail-angle",
        "instrument",
        "epoch",
    ];
    known.extend_from_slice(DETECT_OPTIONS);
    let o: HashMap<String, String> = options(args, &known, &["no-trails"])?;
    let image_path = o.get("image").ok_or("--image is required")?;
    if o.contains_key("o") && o.contains_key("output") {
        return Err("-o and --output are the same option: give one".into());
    }
    let out = o.get("o").or_else(|| o.get("output"));
    if let Some(path) = out {
        crate::check_output_path(path)?;
        crate::check_output_against_input(path, image_path, "the --image input")?;
        if let Some(profile) = o.get("instrument") {
            crate::check_output_against_input(path, profile, "the --instrument input")?;
        }
    }
    let image = fayde_fits::read_image(image_path)
        .map_err(|e| crate::solve_cmd::wrong_kind(image_path, "--image", e))?;
    let instrument = match o.get("instrument") {
        Some(p) => Some(fayde_fits::InstrumentProfile::load(p).map_err(|e| format!("{p}: {e}"))?),
        None => None,
    };
    let mut detect = detect_options(&o, instrument.as_ref())?;
    detect.points_only = o.contains_key("no-trails");
    let (sources, mut report) =
        fayde_api::detect_sources(&image, &detect).map_err(|e| format!("{image_path}: {e}"))?;
    report["identity"] = serde_json::json!({
        "fayde": {"version": env!("CARGO_PKG_VERSION"), "commit": fayde_api::commit()},
    });
    crate::out(&serde_json::to_string_pretty(&report).expect("json"))?;
    if let Some(path) = out {
        let (w, h) = (
            u32::try_from(image.width).map_err(|_| "image too wide".to_string())?,
            u32::try_from(image.height).map_err(|_| "image too tall".to_string())?,
        );
        let trail = fayde_api::detected_trail(&report);
        fayde_fits::write_sources_trailed(path, &sources, w, h, trail)
            .map_err(|e| format!("{path}: {e}"))?;
    }
    Ok(())
}
