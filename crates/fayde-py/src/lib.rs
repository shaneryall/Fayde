//! Python bindings: `fayde._native`, wrapped by the pure-Python `fayde` package.
//!
//! The Python side passes sources as float lists and every other option as one JSON object
//! (`Opts`); this module turns them into a [`fayde_api::RequestOptions`], runs the solver with
//! the GIL released, and returns the report as JSON text: the same report `fayde solve`,
//! `fayde verify` and `fayde refit` print.

use std::time::Duration;

use fayde_api::{build_request, AtCheck, AtPosition, RequestOptions, Solver as ApiSolver};
use fayde_core::{PixelConvention, SkyPos};
use fayde_fits::InstrumentProfile;
use fayde_solve::Source;
use fayde_wcs::{Header, Wcs};
use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use serde::Deserialize;
use serde_json::{Map, Value};

/// A header given as a mapping (`{"CRVAL1": 81.07, ...}`) or as 80-character cards.
#[derive(Deserialize)]
#[serde(untagged)]
enum HeaderIn {
    Map(Map<String, Value>),
    Cards(Vec<String>),
}

impl HeaderIn {
    fn to_header(&self) -> Result<Header, String> {
        match self {
            HeaderIn::Cards(cards) => Ok(Header::parse_cards(cards.iter().map(String::as_str))),
            HeaderIn::Map(m) => fayde_api::header_from_json(m),
        }
    }
}

/// Every option the Python functions accept besides the sources. Unknown keys are errors, so a
/// misspelt option is reported instead of ignored.
#[derive(Deserialize, Default)]
#[serde(deny_unknown_fields)]
struct Opts {
    width: Option<u32>,
    height: Option<u32>,
    ra: Option<f64>,
    dec: Option<f64>,
    radius: Option<f64>,
    pointing_sigma: Option<f64>,
    scale: Option<f64>,
    scale_sigma: Option<f64>,
    roll: Option<f64>,
    roll_sigma: Option<f64>,
    roll_max: Option<f64>,
    parity: Option<String>,
    trail_length: Option<f64>,
    trail_angle: Option<f64>,
    epoch: Option<f64>,
    deadline: Option<f64>,
    max_sip_order: Option<usize>,
    pixel_error: Option<f64>,
    accept_log_odds: Option<f64>,
    max_hypotheses: Option<u64>,
    #[serde(default)]
    no_priors: bool,
    #[serde(default)]
    no_pointing: bool,
    #[serde(default)]
    search_outside_pointing: bool,
    #[serde(default)]
    relax_priors_on_failure: bool,
    header: Option<HeaderIn>,
    /// Instrument profile as an object (the profile file's JSON).
    instrument: Option<Value>,
    /// Instrument profile file.
    instrument_path: Option<String>,
    /// WCS for verify/refit.
    wcs: Option<HeaderIn>,
    /// Pixel origin of the input `x`, `y`: 0 (NumPy, default) or 1 (FITS).
    #[serde(default)]
    origin: u8,
    /// `stars_at` only: 1σ of the WCS, arcsec.
    wcs_sigma_arcsec: Option<f64>,
    /// `stars_at` only: search radius, arcsec.
    search_radius_arcsec: Option<f64>,
    /// `stars_at` only: stars reported per position.
    max_stars: Option<usize>,
    /// Images only: stellar FWHM, px (default: estimated).
    fwhm: Option<f64>,
    /// Images only: pixel threshold of the filtered image, σ.
    detect_threshold: Option<f64>,
    /// Images only: saturation level (default: the header's SATURATE).
    saturation: Option<f64>,
    /// Images only: gain, e⁻ per image unit (default: the header's GAIN/EGAIN).
    gain: Option<f64>,
    /// Images only: most sources kept.
    max_sources: Option<usize>,
}

impl Opts {
    fn detect_options(&self) -> fayde_api::DetectOptions {
        fayde_api::DetectOptions {
            fwhm: self.fwhm,
            threshold: self.detect_threshold,
            saturation: self.saturation,
            gain: self.gain,
            max_sources: self.max_sources,
        }
    }

    fn has_detect_options(&self) -> bool {
        self.fwhm.is_some()
            || self.detect_threshold.is_some()
            || self.saturation.is_some()
            || self.gain.is_some()
            || self.max_sources.is_some()
    }
}

/// An image from a FITS path, or from little-endian `f32` bytes (`width × height`, row by row,
/// first row = FITS row 1) with the header the caller gives.
fn load_image(
    path: Option<String>,
    data: Option<&[u8]>,
    width: usize,
    height: usize,
    header: Option<&HeaderIn>,
) -> Result<fayde_fits::Image, String> {
    let mut image = match (path, data) {
        (Some(p), None) => fayde_fits::read_image(&p).map_err(|e| format!("{p}: {e}"))?,
        (None, Some(bytes)) => {
            let n = width
                .checked_mul(height)
                .filter(|n| n.checked_mul(4) == Some(bytes.len()))
                .ok_or_else(|| {
                    format!(
                        "image data is {} bytes, not {width} × {height} float32 values",
                        bytes.len()
                    )
                })?;
            let data = bytes
                .chunks_exact(4)
                .map(|c| f32::from_le_bytes([c[0], c[1], c[2], c[3]]))
                .collect::<Vec<_>>();
            debug_assert_eq!(data.len(), n);
            fayde_fits::Image {
                width,
                height,
                data,
                header: Header::default(),
                hdu: 0,
            }
        }
        _ => return Err("give an image path or image data".into()),
    };
    if let Some(h) = header {
        image.header = h.to_header()?;
    }
    Ok(image)
}

/// Detect stars in an image (path or `f32` bytes); returns the detection report (JSON).
#[pyfunction]
#[pyo3(signature = (path, data, width, height, options_json))]
fn detect(
    py: Python<'_>,
    path: Option<String>,
    data: Option<&[u8]>,
    width: usize,
    height: usize,
    options_json: &str,
) -> PyResult<String> {
    // detect() takes the detection options and a header; anything else would be ignored.
    let given: serde_json::Map<String, Value> = serde_json::from_str(options_json)
        .map_err(|e| PyValueError::new_err(format!("options: {e}")))?;
    const ALLOWED: &[&str] = &[
        "header",
        "fwhm",
        "detect_threshold",
        "saturation",
        "gain",
        "max_sources",
    ];
    if let Some(k) = given.keys().find(|k| !ALLOWED.contains(&k.as_str())) {
        return Err(PyValueError::new_err(format!(
            "detect() does not take {k}: its options are {}",
            ALLOWED.join(", ")
        )));
    }
    let o: Opts = serde_json::from_value(Value::Object(given))
        .map_err(|e| PyValueError::new_err(format!("options: {e}")))?;
    let image =
        load_image(path, data, width, height, o.header.as_ref()).map_err(PyValueError::new_err)?;
    let opts = o.detect_options();
    let (_, report) = py
        .detach(|| fayde_api::detect_sources(&image, &opts))
        .map_err(PyValueError::new_err)?;
    serde_json::to_string(&report).map_err(|e| PyRuntimeError::new_err(e.to_string()))
}

/// Sources from the Python lists, in the internal (0-based) convention.
fn sources(
    x: Vec<f64>,
    y: Vec<f64>,
    flux: Option<Vec<f64>>,
    sigma_px: Option<Vec<f64>>,
    p_star: Option<Vec<f64>>,
    truncated: Option<Vec<bool>>,
    origin: u8,
) -> Result<Vec<Source>, String> {
    let n = x.len();
    if y.len() != n {
        return Err(format!("x has {n} values, y has {}", y.len()));
    }
    let check = |name: &str, len: Option<usize>| match len {
        Some(l) if l != n => Err(format!("{name} has {l} values, x has {n}")),
        _ => Ok(()),
    };
    check("flux", flux.as_ref().map(Vec::len))?;
    check("sigma_px", sigma_px.as_ref().map(Vec::len))?;
    check("p_star", p_star.as_ref().map(Vec::len))?;
    check("truncated", truncated.as_ref().map(Vec::len))?;
    let convention = match origin {
        0 => PixelConvention::ZeroBased,
        1 => PixelConvention::Fits,
        o => return Err(format!("origin must be 0 or 1, not {o}")),
    };
    // Same rules as the source-list reader (fayde_fits::sources): a non-finite value is absent (NaN
    // marks "no value"; NumPy has no missing floats), and so is a sigma that is not > 0 or a
    // star probability outside [0, 1].
    let opt = |v: &Option<Vec<f64>>, i: usize| v.as_ref().map(|v| v[i]).filter(|v| v.is_finite());
    Ok((0..n)
        .map(|i| {
            let p = convention.to_internal(x[i], y[i]);
            let mut s = Source::at(p.x, p.y);
            s.flux = opt(&flux, i);
            s.sigma_px = opt(&sigma_px, i).filter(|v| *v > 0.0);
            s.p_star = opt(&p_star, i).filter(|v| (0.0..=1.0).contains(v));
            s.truncated = truncated.as_ref().is_some_and(|t| t[i]);
            s
        })
        .collect())
}

fn request_options(o: &Opts) -> Result<RequestOptions, String> {
    let instrument = match (&o.instrument, &o.instrument_path) {
        (Some(_), Some(_)) => return Err("give instrument or instrument_path, not both".into()),
        (Some(v), None) => {
            Some(InstrumentProfile::from_json(&v.to_string()).map_err(|e| e.to_string())?)
        }
        (None, Some(p)) => Some(InstrumentProfile::load(p).map_err(|e| e.to_string())?),
        (None, None) => None,
    };
    let deadline = o
        .deadline
        .map(|s| {
            Duration::try_from_secs_f64(s)
                .map_err(|_| format!("deadline: {s} is not a valid number of seconds"))
        })
        .transpose()?;
    Ok(RequestOptions {
        width: o.width,
        height: o.height,
        ra: o.ra,
        dec: o.dec,
        radius: o.radius,
        pointing_sigma: o.pointing_sigma,
        scale: o.scale,
        scale_sigma: o.scale_sigma,
        roll: o.roll,
        roll_sigma: o.roll_sigma,
        roll_max: o.roll_max,
        parity: o.parity.clone(),
        trail_length: o.trail_length,
        trail_angle: o.trail_angle,
        epoch: o.epoch,
        deadline,
        max_sip_order: o.max_sip_order,
        pixel_error: o.pixel_error,
        accept_log_odds: o.accept_log_odds,
        max_hypotheses: o.max_hypotheses,
        no_priors: o.no_priors,
        no_pointing: o.no_pointing,
        search_outside_pointing: o.search_outside_pointing,
        relax_priors_on_failure: o.relax_priors_on_failure,
        header: o.header.as_ref().map(HeaderIn::to_header).transpose()?,
        instrument,
    })
}

/// A UCAC4 catalogue and pattern indexes, opened once. Thread-safe: solves from several Python
/// threads run in parallel (the GIL is released while solving).
#[pyclass(frozen, module = "fayde._native")]
struct Solver {
    inner: ApiSolver,
}

#[pymethods]
impl Solver {
    #[new]
    #[pyo3(signature = (ucac4, indexes = Vec::new()))]
    fn new(py: Python<'_>, ucac4: String, indexes: Vec<String>) -> PyResult<Self> {
        let inner = py
            .detach(|| ApiSolver::open(&ucac4, &indexes))
            .map_err(PyValueError::new_err)?;
        Ok(Self { inner })
    }

    /// Run `mode` ("solve", "verify" or "refit"). Returns `(report_json, error)`.
    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature = (mode, x, y, flux, sigma_px, p_star, truncated, options_json))]
    fn run(
        &self,
        py: Python<'_>,
        mode: &str,
        x: Vec<f64>,
        y: Vec<f64>,
        flux: Option<Vec<f64>>,
        sigma_px: Option<Vec<f64>>,
        p_star: Option<Vec<f64>>,
        truncated: Option<Vec<bool>>,
        options_json: &str,
    ) -> PyResult<(String, bool)> {
        let o: Opts = serde_json::from_str(options_json)
            .map_err(|e| PyValueError::new_err(format!("options: {e}")))?;
        if o.has_detect_options() {
            return Err(PyValueError::new_err(
                "fwhm, detect_threshold, saturation, gain and max_sources apply to images",
            ));
        }
        let src = sources(x, y, flux, sigma_px, p_star, truncated, o.origin)
            .map_err(PyValueError::new_err)?;
        self.finish(py, mode, src, (None, None), &o, None, None)
    }

    /// Like `run`, with the sources detected in an image (path or `f32` bytes). The image's
    /// header is used for the priors unless the options give one.
    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature = (mode, path, data, width, height, options_json))]
    fn run_image(
        &self,
        py: Python<'_>,
        mode: &str,
        path: Option<String>,
        data: Option<&[u8]>,
        width: usize,
        height: usize,
        options_json: &str,
    ) -> PyResult<(String, bool)> {
        let o: Opts = serde_json::from_str(options_json)
            .map_err(|e| PyValueError::new_err(format!("options: {e}")))?;
        let image = load_image(path, data, width, height, o.header.as_ref())
            .map_err(PyValueError::new_err)?;
        let opts = o.detect_options();
        let (src, mut detection) = py
            .detach(|| fayde_api::detect_sources(&image, &opts))
            .map_err(PyValueError::new_err)?;
        if let Some(m) = detection.as_object_mut() {
            m.remove("sources");
            m.remove("rejected");
        }
        let size = (
            u32::try_from(image.width).ok(),
            u32::try_from(image.height).ok(),
        );
        self.finish(py, mode, src, size, &o, Some(image.header), Some(detection))
    }

    /// "Is there a catalogue star here?" under the WCS in `options_json` (`fayde verify --at`).
    /// Pixel positions `x`, `y` (with optional `sigma_px`) and sky positions `ra`, `dec`
    /// (degrees, optional `sigma_arcsec`); NaN sigma means "default". Returns the report JSON
    /// and whether it is an error.
    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature = (x, y, sigma_px, ra, dec, sigma_arcsec, options_json))]
    fn at(
        &self,
        py: Python<'_>,
        x: Vec<f64>,
        y: Vec<f64>,
        sigma_px: Option<Vec<f64>>,
        ra: Vec<f64>,
        dec: Vec<f64>,
        sigma_arcsec: Option<Vec<f64>>,
        options_json: &str,
    ) -> PyResult<(String, bool)> {
        let o: Opts = serde_json::from_str(options_json)
            .map_err(|e| PyValueError::new_err(format!("options: {e}")))?;
        let err = |m: &str| PyValueError::new_err(m.to_string());
        let convention = match o.origin {
            0 => PixelConvention::ZeroBased,
            1 => PixelConvention::Fits,
            v => return Err(err(&format!("origin must be 0 or 1, not {v}"))),
        };
        if x.len() != y.len() || ra.len() != dec.len() {
            return Err(err("x/y and ra/dec must have equal lengths"));
        }
        let len_ok = |v: &Option<Vec<f64>>, n: usize| v.as_ref().is_none_or(|v| v.len() == n);
        if !len_ok(&sigma_px, x.len()) || !len_ok(&sigma_arcsec, ra.len()) {
            return Err(err("sigma_px / sigma_arcsec must match their positions"));
        }
        let sig = |v: &Option<Vec<f64>>, i: usize| v.as_ref().map(|v| v[i]).filter(|v| !v.is_nan());
        let mut positions = Vec::with_capacity(x.len() + ra.len());
        for i in 0..x.len() {
            positions.push(AtPosition::Pixel(
                convention.to_internal(x[i], y[i]),
                sig(&sigma_px, i),
            ));
        }
        for i in 0..ra.len() {
            let sky = SkyPos::from_degrees(ra[i], dec[i])
                .ok_or_else(|| err(&format!("bad sky position ({}, {})", ra[i], dec[i])))?;
            positions.push(AtPosition::Sky(sky, sig(&sigma_arcsec, i)));
        }
        let wcs = match &o.wcs {
            Some(w) => {
                let h = w.to_header().map_err(|e| err(&e))?;
                Wcs::from_header(&h).map_err(|e| err(&format!("wcs: {e}")))?
            }
            None => return Err(err("stars_at needs wcs")),
        };
        let opts = request_options(&o).map_err(|e| err(&e))?;
        let at = AtCheck {
            wcs_sigma_arcsec: o.wcs_sigma_arcsec,
            search_radius_arcsec: o.search_radius_arcsec,
            max_stars: o.max_stars,
        };
        let report = py.detach(|| self.inner.check_at(&wcs, &positions, &opts, &at));
        let text = serde_json::to_string(&report.json)
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        Ok((text, report.error))
    }
}

impl Solver {
    /// Shared tail of `run` and `run_image`.
    #[allow(clippy::too_many_arguments)]
    fn finish(
        &self,
        py: Python<'_>,
        mode: &str,
        src: Vec<Source>,
        size: (Option<u32>, Option<u32>),
        o: &Opts,
        image_header: Option<Header>,
        detection: Option<Value>,
    ) -> PyResult<(String, bool)> {
        let mut opts = request_options(o).map_err(PyValueError::new_err)?;
        if opts.header.is_none() {
            opts.header = image_header;
        }
        let wcs = match (&o.wcs, mode) {
            (Some(w), "verify" | "refit") => {
                let h = w.to_header().map_err(PyValueError::new_err)?;
                Some(Wcs::from_header(&h).map_err(|e| PyValueError::new_err(format!("wcs: {e}")))?)
            }
            (None, "verify" | "refit") => {
                return Err(PyValueError::new_err(format!("{mode} needs wcs")))
            }
            (Some(_), _) => return Err(PyValueError::new_err("wcs is for verify and refit")),
            (None, _) => None,
        };
        if !matches!(mode, "solve" | "verify" | "refit") {
            return Err(PyValueError::new_err(format!("unknown mode {mode}")));
        }
        if o.wcs_sigma_arcsec.is_some() || o.search_radius_arcsec.is_some() || o.max_stars.is_some()
        {
            return Err(PyValueError::new_err(
                "wcs_sigma_arcsec, search_radius_arcsec and max_stars are stars_at options",
            ));
        }
        let built = build_request(src, size, &opts).map_err(PyValueError::new_err)?;
        let mut report = py.detach(|| match (mode, &wcs) {
            ("verify", Some(w)) => self.inner.verify(&built, w),
            ("refit", Some(w)) => self.inner.refit(&built, w),
            _ => self.inner.solve(&built),
        });
        if let Some(d) = detection {
            report.json["detection"] = d;
        }
        let text = serde_json::to_string(&report.json)
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        Ok((text, report.error))
    }
}

/// The compiled part of the `fayde` package.
#[pymodule]
fn _native(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Solver>()?;
    m.add_function(wrap_pyfunction!(detect, m)?)?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add("commit", fayde_api::commit())?;
    Ok(())
}
