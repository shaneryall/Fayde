//! Instrument profiles (design doc §6.5): what is known about a camera independently of any one
//! frame, stored once as JSON and applied to every frame's priors.
//!
//! ```json
//! {
//!   "name": "Example camera",
//!   "scale_arcsec": 3.0,
//!   "scale_sigma_arcsec": 0.003,
//!   "parity": "positive",
//!   "roll": { "keyword": "ORIENTAT", "offset_deg": 0.0, "sign": 1,
//!             "sigma_deg": 0.1, "max_dev_deg": 1.0 },
//!   "binning": { "x_keyword": "XBINNING", "y_keyword": "YBINNING" }
//! }
//! ```
//!
//! `scale_arcsec` is the unbinned scale. With a binning rule, the scale and its σ are multiplied
//! by `√(bx·by)` read from the header; a frame whose binning keywords are missing or invalid
//! gets no scale prior from the profile (its scale is unknown, not unbinned).
//!
//! Roll follows Fayde's convention: the position angle of the image +y axis, east of north
//! (`atan2(CD1_2, CD2_2)`). The rule `roll = offset + sign × header[keyword]` maps a spacecraft
//! or mount keyword onto it; its offset and sign are calibrated from solved frames. Instrument
//! rules live here, not in the solver, so the core crates stay instrument-agnostic.

use std::path::Path;

use fayde_wcs::{Header, Parity, Value};
use serde::Deserialize;

use crate::FitsError;

/// A roll prior derived from a header keyword.
#[derive(Debug, Clone, PartialEq, Deserialize)]
pub struct RollRule {
    /// Header keyword holding the instrument's own roll angle, degrees.
    pub keyword: String,
    /// Offset added after applying `sign`, degrees.
    pub offset_deg: f64,
    /// +1 or −1.
    pub sign: f64,
    /// 1σ of the derived roll, degrees.
    pub sigma_deg: f64,
    /// Hard limit on the deviation from the derived roll, degrees.
    pub max_dev_deg: f64,
}

/// Header keywords holding the on-chip binning factors.
#[derive(Debug, Clone, PartialEq, Deserialize)]
pub struct BinningRule {
    /// Binning along x (columns), e.g. `XBINNING`.
    pub x_keyword: String,
    /// Binning along y (rows), e.g. `YBINNING`; the x factor when absent.
    #[serde(default)]
    pub y_keyword: Option<String>,
}

/// What is known about an instrument.
#[derive(Debug, Clone, PartialEq, Deserialize, Default)]
pub struct InstrumentProfile {
    /// Display name.
    #[serde(default)]
    pub name: Option<String>,
    /// Pixel scale, arcsec/px. Takes precedence over header-derived scales.
    #[serde(default)]
    pub scale_arcsec: Option<f64>,
    /// 1σ of the scale, arcsec/px.
    #[serde(default)]
    pub scale_sigma_arcsec: Option<f64>,
    /// `"negative"` (sky as seen: north up, east left) or `"positive"` (mirrored).
    #[serde(default)]
    pub parity: Option<String>,
    /// Roll from a header keyword.
    #[serde(default)]
    pub roll: Option<RollRule>,
    /// Binning from header keywords: scales `scale_arcsec` per frame.
    #[serde(default)]
    pub binning: Option<BinningRule>,
    /// 1σ of the pointing a frame's header gives, degrees: how far the header's pointing is
    /// from the true field centre, per axis. With it the search starts in a small cone around
    /// the pointing and widens to the hard limit only when that finds nothing.
    #[serde(default)]
    pub pointing_sigma_deg: Option<f64>,
    /// Hard limit on the pointing error, degrees: the search goes no farther from the header's
    /// pointing.
    #[serde(default)]
    pub pointing_radius_deg: Option<f64>,
    /// Whether a world coordinate system in a frame's header may supply the pointing and the
    /// scale (the default). `false` for an instrument whose headers can carry a plate solver's
    /// solution: the pointing then comes from the header's pointing keywords only.
    #[serde(default)]
    pub header_wcs: Option<bool>,
}

/// Priors an instrument profile gives for one frame.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct InstrumentPriors {
    /// Scale and 1σ, arcsec/px.
    pub scale: Option<(f64, f64)>,
    /// Roll, 1σ and hard limit, radians.
    pub roll: Option<(f64, f64, f64)>,
    /// Parity.
    pub parity: Option<Parity>,
}

/// Parse `"negative"`/`"positive"` (also `neg`/`pos`).
pub fn parse_parity(s: &str) -> Option<Parity> {
    match s.trim().to_ascii_lowercase().as_str() {
        "negative" | "neg" => Some(Parity::Negative),
        "positive" | "pos" => Some(Parity::Positive),
        _ => None,
    }
}

impl InstrumentProfile {
    /// Read a profile from a JSON file and validate it.
    pub fn load(path: impl AsRef<Path>) -> Result<Self, FitsError> {
        Self::from_json(&std::fs::read_to_string(path)?)
    }

    /// Parse a profile from JSON text and validate it.
    pub fn from_json(text: &str) -> Result<Self, FitsError> {
        let p: Self = serde_json::from_str(text)
            .map_err(|e| FitsError::Invalid(format!("instrument profile: {e}")))?;
        p.validate()?;
        Ok(p)
    }

    fn validate(&self) -> Result<(), FitsError> {
        let bad = |s: &str| Err(FitsError::Invalid(format!("instrument profile: {s}")));
        if self
            .scale_arcsec
            .is_some_and(|s| !(s.is_finite() && s > 0.0))
        {
            return bad("scale_arcsec must be positive");
        }
        if self
            .scale_sigma_arcsec
            .is_some_and(|s| !(s.is_finite() && s >= 0.0))
        {
            return bad("scale_sigma_arcsec must be >= 0");
        }
        if let Some(p) = &self.parity {
            if parse_parity(p).is_none() {
                return bad("parity must be \"negative\" or \"positive\"");
            }
        }
        for (name, v) in [
            ("pointing_sigma_deg", self.pointing_sigma_deg),
            ("pointing_radius_deg", self.pointing_radius_deg),
        ] {
            if v.is_some_and(|v| !(v.is_finite() && v > 0.0)) {
                return bad(&format!("{name} must be positive"));
            }
        }
        if let Some(r) = &self.roll {
            if r.sign != 1.0 && r.sign != -1.0 {
                return bad("roll.sign must be 1 or -1");
            }
            if !(r.sigma_deg > 0.0 && r.max_dev_deg > 0.0 && r.offset_deg.is_finite()) {
                return bad("roll needs finite offset_deg and positive sigma_deg, max_dev_deg");
            }
        }
        Ok(())
    }

    /// The priors this profile gives for a frame with header `h`. A roll rule whose keyword is
    /// missing or not numeric gives no roll prior.
    pub fn priors(&self, h: &Header) -> InstrumentPriors {
        let read = |k: &str| match h.get(k)? {
            Value::Str(s) => s.trim().parse::<f64>().ok(),
            v => v.as_f64(),
        };
        // Binning factor √(bx·by); None when the rule's keywords are missing or not positive.
        let bin = match &self.binning {
            None => Some(1.0),
            Some(b) => {
                let bx = read(&b.x_keyword);
                let by = b.y_keyword.as_deref().map_or(bx, read);
                match (bx, by) {
                    (Some(x), Some(y))
                        if x >= 1.0 && y >= 1.0 && x.is_finite() && y.is_finite() =>
                    {
                        Some((x * y).sqrt())
                    }
                    _ => None,
                }
            }
        };
        let scale = self.scale_arcsec.zip(bin).map(|(s, b)| {
            // Without a stated σ, 0.5% is a conservative instrument-level default.
            (s * b, self.scale_sigma_arcsec.unwrap_or(0.005 * s) * b)
        });
        let roll = self.roll.as_ref().and_then(|r| {
            let v: f64 = match h.get(&r.keyword)? {
                Value::Str(s) => s.trim().parse().ok()?,
                v => v.as_f64()?,
            };
            if !v.is_finite() {
                return None;
            }
            let deg: f64 = r.offset_deg + r.sign * v;
            Some((
                deg.to_radians(),
                r.sigma_deg.to_radians(),
                r.max_dev_deg.to_radians(),
            ))
        });
        InstrumentPriors {
            scale,
            roll,
            parity: self.parity.as_deref().and_then(parse_parity),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roll_rule_and_validation() {
        let p: InstrumentProfile = serde_json::from_str(
            r#"{"scale_arcsec": 2.99, "parity": "positive",
                "roll": {"keyword": "ORIENTAT", "offset_deg": 90.0, "sign": -1,
                         "sigma_deg": 0.1, "max_dev_deg": 1.0}}"#,
        )
        .unwrap();
        p.validate().unwrap();
        let h = Header::parse_cards(["ORIENTAT=                 30.0"]);
        let pr = p.priors(&h);
        assert!((pr.roll.unwrap().0.to_degrees() - 60.0).abs() < 1e-12);
        assert_eq!(pr.parity, Some(Parity::Positive));
        assert!((pr.scale.unwrap().1 - 0.005 * 2.99).abs() < 1e-15);
        // Missing keyword: no roll prior.
        assert!(p.priors(&Header::default()).roll.is_none());
        // Binning multiplies the scale and its σ; missing binning keywords drop the scale.
        let mut b = p.clone();
        b.binning = Some(BinningRule {
            x_keyword: "XBINNING".into(),
            y_keyword: Some("YBINNING".into()),
        });
        let h2 = Header::parse_cards([
            "XBINNING=                    2",
            "YBINNING=                    2",
        ]);
        let (s, sig) = b.priors(&h2).scale.unwrap();
        assert!((s - 5.98).abs() < 1e-12 && (sig - 2.0 * 0.005 * 2.99).abs() < 1e-12);
        let h1 = Header::parse_cards([
            "XBINNING=                    1",
            "YBINNING=                    1",
        ]);
        assert!((b.priors(&h1).scale.unwrap().0 - 2.99).abs() < 1e-12);
        assert!(b.priors(&Header::default()).scale.is_none());
        let bad: InstrumentProfile = serde_json::from_str(r#"{"parity": "sideways"}"#).unwrap();
        assert!(bad.validate().is_err());
        // Pointing uncertainty and the header-WCS switch are optional and validated.
        assert_eq!(
            (p.pointing_sigma_deg, p.pointing_radius_deg, p.header_wcs),
            (None, None, None)
        );
        let q = InstrumentProfile::from_json(
            r#"{"pointing_sigma_deg": 0.02, "pointing_radius_deg": 1.5, "header_wcs": false}"#,
        )
        .unwrap();
        assert_eq!(
            (q.pointing_sigma_deg, q.pointing_radius_deg, q.header_wcs),
            (Some(0.02), Some(1.5), Some(false))
        );
        assert!(InstrumentProfile::from_json(r#"{"pointing_sigma_deg": 0}"#).is_err());
        assert!(InstrumentProfile::from_json(r#"{"pointing_radius_deg": -1}"#).is_err());
    }
}
