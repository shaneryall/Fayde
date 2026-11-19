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
        crate::check_output_against_input(path, image_path, "the --image input")?;
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
