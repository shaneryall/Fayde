use std::path::Path;

use fayde_wcs::{Header, Parity, Value};
use serde::Deserialize;

use crate::FitsError;

#[derive(Debug, Clone, PartialEq, Deserialize)]
pub struct RollRule {
    pub keyword: String,
    pub offset_deg: f64,
    pub sign: f64,
    pub sigma_deg: f64,
    pub max_dev_deg: f64,
}

#[derive(Debug, Clone, PartialEq, Deserialize)]
pub struct BinningRule {
    pub x_keyword: String,
    #[serde(default)]
    pub y_keyword: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Deserialize, Default)]
pub struct InstrumentProfile {
    #[serde(default)]
    pub name: Option<String>,
    #[serde(default)]
    pub scale_arcsec: Option<f64>,
    #[serde(default)]
    pub scale_sigma_arcsec: Option<f64>,
    #[serde(default)]
    pub parity: Option<String>,
    #[serde(default)]
    pub roll: Option<RollRule>,
    #[serde(default)]
    pub binning: Option<BinningRule>,
    #[serde(default)]
    pub pointing_sigma_deg: Option<f64>,
    #[serde(default)]
    pub pointing_radius_deg: Option<f64>,
    #[serde(default)]
    pub header_wcs: Option<bool>,
}

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct InstrumentPriors {
    pub scale: Option<(f64, f64)>,
    pub roll: Option<(f64, f64, f64)>,
    pub parity: Option<Parity>,
}

pub fn parse_parity(s: &str) -> Option<Parity> {
    match s.trim().to_ascii_lowercase().as_str() {
        "negative" | "neg" => Some(Parity::Negative),
        "positive" | "pos" => Some(Parity::Positive),
        _ => None,
    }
}

impl InstrumentProfile {
    pub fn load(path: impl AsRef<Path>) -> Result<Self, FitsError> {
        Self::from_json(&std::fs::read_to_string(path)?)
    }

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

    pub fn priors(&self, h: &Header) -> InstrumentPriors {
        let read = |k: &str| match h.get(k)? {
            Value::Str(s) => s.trim().parse::<f64>().ok(),
            v => v.as_f64(),
        };
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
        let scale = self
            .scale_arcsec
            .zip(bin)
            .map(|(s, b)| (s * b, self.scale_sigma_arcsec.unwrap_or(0.005 * s) * b));
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
        assert!(p.priors(&Header::default()).roll.is_none());
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
