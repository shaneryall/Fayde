use fayde_core::{PixelPos, Portable, SkyPos};

use crate::{TanWcs, WcsError};

#[derive(Debug, Clone, PartialEq)]
pub enum Value {
    Str(String),
    Bool(bool),
    Int(i64),
    Float(f64),
}

impl Value {
    pub fn as_f64(&self) -> Option<f64> {
        match *self {
            Self::Int(i) => Some(i as f64),
            Self::Float(f) => Some(f),
            _ => None,
        }
    }

    pub fn as_str(&self) -> Option<&str> {
        match self {
            Self::Str(s) => Some(s),
            _ => None,
        }
    }
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct Header {
    cards: Vec<(String, Result<Value, String>)>,
}

fn is_wcs_keyword(key: &str) -> bool {
    let digits = |s: &str| !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit());
    let numbered = |stem: &str| {
        key.strip_prefix(stem).is_some_and(|rest| {
            let rest = rest
                .strip_suffix(|c: char| c.is_ascii_uppercase())
                .unwrap_or(rest);
            digits(rest)
        })
    };
    let indexed = |stem: &str| {
        key.strip_prefix(stem).is_some_and(|rest| {
            let rest = rest
                .strip_suffix(|c: char| c.is_ascii_uppercase())
                .unwrap_or(rest);
            rest.split_once('_')
                .is_some_and(|(i, j)| digits(i) && digits(j))
        })
    };
    ["CRVAL", "CRPIX", "CDELT", "CROTA", "CTYPE", "CUNIT"]
        .iter()
        .any(|s| numbered(s))
        || ["CD", "PC", "PV", "PS", "A_", "B_", "AP_", "BP_"]
            .iter()
            .any(|s| indexed(s))
        || matches!(
            key,
            "LONPOLE"
                | "LATPOLE"
                | "WCSAXES"
                | "WCSNAME"
                | "A_ORDER"
                | "B_ORDER"
                | "AP_ORDER"
                | "BP_ORDER"
                | "A_DMAX"
                | "B_DMAX"
        )
}

impl Header {
    pub fn parse_cards<'a>(cards: impl IntoIterator<Item = &'a str>) -> Self {
        let mut out = Vec::new();
        for card in cards {
            if let Some((key, value)) = parse_card(card) {
                out.push((key, value));
            }
        }
        Self { cards: out }
    }

    pub fn parse_block(text: &str) -> Self {
        if text.contains('\n') {
            Self::parse_cards(text.lines())
        } else {
            let chars: Vec<char> = text.chars().collect();
            let cards: Vec<String> = chars.chunks(80).map(|c| c.iter().collect()).collect();
            Self::parse_cards(cards.iter().map(String::as_str))
        }
    }

    pub fn set(&mut self, key: &str, value: Value) {
        let key = key.to_ascii_uppercase();
        match self.cards.iter_mut().find(|(k, _)| *k == key) {
            Some(card) => card.1 = Ok(value),
            None => self.cards.push((key, Ok(value))),
        }
    }

    pub fn get(&self, key: &str) -> Option<&Value> {
        let key = key.to_ascii_uppercase();
        self.cards
            .iter()
            .filter(|(k, _)| *k == key)
            .find_map(|(_, v)| v.as_ref().ok())
    }

    pub fn get_checked(&self, key: &str) -> Result<Option<&Value>, WcsError> {
        let key = key.to_ascii_uppercase();
        self.cards
            .iter()
            .find(|(k, _)| *k == key)
            .map(|(_, value)| {
                value.as_ref().map_err(|raw| WcsError::InvalidValue {
                    key,
                    reason: format!("cannot parse '{raw}'"),
                })
            })
            .transpose()
    }

    pub fn contains(&self, key: &str) -> bool {
        self.get(key).is_some()
    }

    pub fn without_wcs(&self) -> Self {
        Self {
            cards: self
                .cards
                .iter()
                .filter(|(k, _)| !is_wcs_keyword(k))
                .cloned()
                .collect(),
        }
    }

    pub fn number(&self, key: &str) -> Result<Option<f64>, WcsError> {
        self.get_checked(key)?
            .map(|v| {
                v.as_f64()
                    .filter(|x| x.is_finite())
                    .ok_or_else(|| WcsError::InvalidValue {
                        key: key.to_string(),
                        reason: format!("expected a finite number, got {v:?}"),
                    })
            })
            .transpose()
    }

    fn f64_or(&self, key: &str, default: f64) -> Result<f64, WcsError> {
        Ok(self.number(key)?.unwrap_or(default))
    }

    fn f64_required(&self, key: &str) -> Result<f64, WcsError> {
        self.number(key)?
            .ok_or_else(|| WcsError::MissingKey(key.to_string()))
    }

    pub(crate) fn str_trimmed(&self, key: &str) -> Result<Option<&str>, WcsError> {
        self.get_checked(key)?
            .map(|v| {
                v.as_str()
                    .map(str::trim)
                    .ok_or_else(|| WcsError::InvalidValue {
                        key: key.to_string(),
                        reason: format!("expected a string, got {v:?}"),
                    })
            })
            .transpose()
    }
}

fn parse_card(card: &str) -> Option<(String, Result<Value, String>)> {
    let card: String = card.chars().take(80).collect();
    if card.get(8..10) != Some("= ") {
        return None;
    }
    let key = card.get(..8)?.trim_end().to_ascii_uppercase();
    if key.is_empty() {
        return None;
    }
    let field = card.get(10..)?.trim();
    if field.is_empty() || field.starts_with('/') {
        return None;
    }
    Some((key, parse_value(field).ok_or_else(|| field.to_string())))
}

fn parse_value(field: &str) -> Option<Value> {
    let s = field.trim_start();
    if let Some(rest) = s.strip_prefix('\'') {
        let mut out = String::new();
        let mut chars = rest.chars().peekable();
        while let Some(c) = chars.next() {
            if c == '\'' {
                if chars.peek() == Some(&'\'') {
                    out.push('\'');
                    chars.next();
                } else {
                    return Some(Value::Str(out.trim_end().to_string()));
                }
            } else {
                out.push(c);
            }
        }
        return None;
    }
    let token = s.split('/').next()?.trim();
    match token {
        "T" => return Some(Value::Bool(true)),
        "F" => return Some(Value::Bool(false)),
        "" => return None,
        _ => {}
    }
    if let Ok(i) = token.parse::<i64>() {
        return Some(Value::Int(i));
    }
    match token.replace(['D', 'd'], "E").parse::<f64>() {
        Ok(v) if v.is_finite() => Some(Value::Float(v)),
        Ok(_) => Some(Value::Str(token.to_string())),
        Err(_) => None,
    }
}

pub fn check_reference_frame(h: &Header) -> Result<(), WcsError> {
    let sys = match h.str_trimmed("RADESYS")? {
        Some(sys) => Some(sys),
        None => h.str_trimmed("RADECSYS")?,
    }
    .map(str::to_ascii_uppercase);
    let equinox_card = match h.get_checked("EQUINOX")? {
        Some(equinox) => Some(equinox),
        None => h.get_checked("EPOCH")?,
    };
    let equinox = match equinox_card {
        None => None,
        Some(v) => Some(
            v.as_f64()
                .or_else(|| {
                    v.as_str()
                        .and_then(|s| s.trim().trim_start_matches(['J', 'j']).parse().ok())
                })
                .filter(|v| v.is_finite())
                .ok_or_else(|| WcsError::Unsupported(format!("EQUINOX = {v:?}")))?,
        ),
    };
    let j2000 = |e: Option<f64>| e.is_none_or(|e| (e - 2000.0).abs() < 1e-6);
    match sys.as_deref() {
        Some("ICRS") => Ok(()),
        Some("FK5") if j2000(equinox) => Ok(()),
        None if equinox.is_none_or(|e| e >= 1984.0) && j2000(equinox) => Ok(()),
        Some(s) => Err(WcsError::Unsupported(format!(
            "RADESYS = '{s}', EQUINOX = {equinox:?}: only ICRS and FK5 J2000 are supported"
        ))),
        None => Err(WcsError::Unsupported(format!(
            "EQUINOX = {equinox:?} without RADESYS: only ICRS and FK5 J2000 are supported"
        ))),
    }
}

pub fn is_degree_unit(u: &str) -> bool {
    ["deg", "degree", "degrees"]
        .iter()
        .any(|d| u.trim().eq_ignore_ascii_case(d))
}

pub fn format_real(v: f64) -> String {
    let s = format!("{v:e}");
    let (mantissa, exp) = s.split_once('e').expect("{:e} always has an exponent");
    let exp: i32 = exp.parse().expect("valid exponent");
    let mantissa = if mantissa.contains('.') {
        mantissa.to_string()
    } else {
        format!("{mantissa}.0")
    };
    let sign = if exp < 0 { '-' } else { '+' };
    format!("{mantissa}E{sign}{:02}", exp.abs())
}

fn card(key: &str, value: &str, comment: &str) -> String {
    let mut s = format!("{key:<8}= {value:>20}");
    if !comment.is_empty() {
        s.push_str(" / ");
        s.push_str(comment);
    }
    let mut s: String = s.chars().take(80).collect();
    while s.len() < 80 {
        s.push(' ');
    }
    s
}

pub(crate) fn num_card(key: &str, value: &str) -> String {
    card(key, value, "")
}

pub(crate) fn str_card(key: &str, value: &str, comment: &str) -> String {
    let quoted = format!("'{:<8}'", value.replace('\'', "''"));
    card(key, &format!("{quoted:<20}"), comment)
}

impl TanWcs {
    pub fn from_header(h: &Header) -> Result<Self, WcsError> {
        for (key, want) in [("CTYPE1", "RA---TAN"), ("CTYPE2", "DEC--TAN")] {
            let got = h
                .str_trimmed(key)?
                .ok_or_else(|| WcsError::MissingKey(key.to_string()))?;
            if got != want {
                return Err(WcsError::Unsupported(format!("{key} = '{got}'")));
            }
        }
        for key in ["CUNIT1", "CUNIT2"] {
            if let Some(u) = h.str_trimmed(key)? {
                if !is_degree_unit(u) {
                    return Err(WcsError::Unsupported(format!("{key} = '{u}'")));
                }
            }
        }
        check_reference_frame(h)?;
        let dec = h.f64_required("CRVAL2")?;
        let lonpole = h.f64_or("LONPOLE", if dec == 90.0 { 0.0 } else { 180.0 })?;
        let polar_flip = dec.abs() == 90.0 && lonpole.abs() <= 1e-12;
        if !polar_flip && (lonpole - 180.0).abs() > 1e-12 {
            return Err(WcsError::Unsupported(format!("LONPOLE = {lonpole}")));
        }

        const CD: [&str; 4] = ["CD1_1", "CD1_2", "CD2_1", "CD2_2"];
        const PC: [&str; 4] = ["PC1_1", "PC1_2", "PC2_1", "PC2_2"];
        let matrix_present = |keys: [&str; 4]| -> Result<bool, WcsError> {
            let mut present = false;
            for key in keys {
                present |= h.get_checked(key)?.is_some();
            }
            Ok(present)
        };
        let has_cd = matrix_present(CD)?;
        let has_pc = matrix_present(PC)?;
        let cd = match (has_cd, has_pc) {
            (true, true) => {
                return Err(WcsError::Unsupported(
                    "both CDi_j and PCi_j present (ambiguous)".into(),
                ))
            }
            (true, false) => [
                [h.f64_or(CD[0], 0.0)?, h.f64_or(CD[1], 0.0)?],
                [h.f64_or(CD[2], 0.0)?, h.f64_or(CD[3], 0.0)?],
            ],
            (false, true) => {
                let c1 = h.f64_or("CDELT1", 1.0)?;
                let c2 = h.f64_or("CDELT2", 1.0)?;
                [
                    [c1 * h.f64_or(PC[0], 1.0)?, c1 * h.f64_or(PC[1], 0.0)?],
                    [c2 * h.f64_or(PC[2], 0.0)?, c2 * h.f64_or(PC[3], 1.0)?],
                ]
            }
            (false, false) => {
                let c1 = h.f64_or("CDELT1", 1.0)?;
                let c2 = h.f64_or("CDELT2", 1.0)?;
                let rho = match h.number("CROTA2")? {
                    Some(rho) => rho,
                    None => h.f64_or("CROTA1", 0.0)?,
                };
                let (s, c) = rho.to_radians().psin_cos();
                [[c1 * c, -c2 * s], [c1 * s, c2 * c]]
            }
        };

        let ra = h.f64_required("CRVAL1")?;
        let cd = if polar_flip {
            cd.map(|row| row.map(|v| -v))
        } else {
            cd
        };
        let crval = SkyPos::from_degrees(ra, dec).ok_or_else(|| WcsError::InvalidValue {
            key: "CRVAL1/CRVAL2".into(),
            reason: format!("({ra}, {dec}) is not a valid sky position"),
        })?;
        let crpix = PixelPos::new(
            h.f64_required("CRPIX1")? - 1.0,
            h.f64_required("CRPIX2")? - 1.0,
        );
        TanWcs::new(crval, crpix, cd)
    }

    pub fn to_cards(&self) -> Vec<String> {
        let (ra, dec) = self.crval().to_degrees();
        let p = self.crpix();
        let cd = self.cd();
        vec![
            card("WCSAXES", "2", "Number of coordinate axes"),
            str_card("CTYPE1", "RA---TAN", "Gnomonic projection"),
            str_card("CTYPE2", "DEC--TAN", "Gnomonic projection"),
            str_card("CUNIT1", "deg", ""),
            str_card("CUNIT2", "deg", ""),
            str_card("RADESYS", "ICRS", ""),
            card("CRVAL1", &format_real(ra), "[deg] RA of reference point"),
            card("CRVAL2", &format_real(dec), "[deg] Dec of reference point"),
            card(
                "CRPIX1",
                &format_real(p.x + 1.0),
                "Reference pixel (FITS 1-based)",
            ),
            card(
                "CRPIX2",
                &format_real(p.y + 1.0),
                "Reference pixel (FITS 1-based)",
            ),
            card("CD1_1", &format_real(cd[0][0]), "[deg/px]"),
            card("CD1_2", &format_real(cd[0][1]), "[deg/px]"),
            card("CD2_1", &format_real(cd[1][0]), "[deg/px]"),
            card("CD2_2", &format_real(cd[1][1]), "[deg/px]"),
            card("LONPOLE", &format_real(180.0), "[deg]"),
        ]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn without_wcs_drops_only_wcs_cards() {
        let wcs = [
            "CRVAL1", "CRVAL2", "CRPIX1", "CRPIX2A", "CDELT1", "CROTA2", "CTYPE1", "CUNIT2",
            "CD1_1", "CD2_1", "PC1_2", "PC2_2B", "PV1_0", "PS1_0", "LONPOLE", "LATPOLE", "WCSAXES",
            "WCSNAME", "A_ORDER", "B_ORDER", "AP_ORDER", "BP_ORDER", "A_2_0", "B_0_2", "AP_1_1",
            "BP_0_3", "A_DMAX", "B_DMAX",
        ];
        let other = [
            "RA", "DEC", "OBJCTRA", "OBJCTDEC", "ROLLANG", "EQUINOX", "RADESYS", "NAXIS1",
            "DATE-OBS", "EXPTIME", "XBINNING", "CDF", "CD", "PCOUNT", "CRVAL", "CTYPE", "A_",
            "B_X_1", "PA", "CDELTA", "AIRMASS", "BSCALE", "PV", "CD1_", "CD_1",
        ];
        let mut h = Header::default();
        for k in wcs.iter().chain(&other) {
            h.set(k, Value::Float(1.0));
        }
        let stripped = h.without_wcs();
        for k in wcs {
            assert!(h.contains(k) && !stripped.contains(k), "{k} should go");
        }
        for k in other {
            assert!(stripped.contains(k), "{k} should stay");
        }
    }

    #[test]
    fn non_ascii_cards_are_skipped_not_a_panic() {
        let h = Header::parse_cards(["1234567é= 1", "KEY     = 'ok'"]);
        assert!(h.get("KEY").is_some());
    }

    #[test]
    fn reference_frames_other_than_icrs_or_fk5_j2000_are_rejected() {
        let ok = |cards: &[&str]| {
            check_reference_frame(&Header::parse_cards(cards.iter().copied())).is_ok()
        };
        assert!(ok(&[]));
        assert!(ok(&["RADESYS = 'ICRS'"]));
        assert!(ok(&["RADESYS = 'FK5'", "EQUINOX =               2000.0"]));
        assert!(ok(&["RADECSYS= 'FK5'", "EPOCH   =               2000.0"]));
        assert!(ok(&["EQUINOX =               2000.0"]));
        assert!(!ok(&["RADESYS = 'FK5'", "EQUINOX =               2050.0"]));
        assert!(!ok(&["RADESYS = 'FK4'"]));
        assert!(!ok(&["EQUINOX =               1950.0"]));
        assert!(!ok(&["RADESYS = 'GAPPT'"]));
    }

    #[test]
    fn parses_value_types() {
        let h = Header::parse_cards([
            "SIMPLE  =                    T / conforms",
            "NAXIS1  =                 1024",
            "CRVAL1  =  1.5D+02 / D exponent",
            "OBJECT  = 'O''Brien  '         / escaped quote",
            "COMMENT this is not a value card",
            "BROKEN  = 'unterminated",
        ]);
        assert_eq!(h.get("simple"), Some(&Value::Bool(true)));
        assert_eq!(h.get("NAXIS1"), Some(&Value::Int(1024)));
        assert_eq!(h.get("CRVAL1").and_then(Value::as_f64), Some(150.0));
        assert_eq!(h.get("OBJECT").and_then(Value::as_str), Some("O'Brien"));
        assert!(!h.contains("COMMENT") && !h.contains("BROKEN"));
    }

    #[test]
    fn real_formatting_is_lossless_and_fits_style() {
        assert_eq!(format_real(-8.3e-4), "-8.3E-04");
        assert_eq!(format_real(150.0), "1.5E+02");
        assert_eq!(format_real(0.0), "0.0E+00");
        for v in [
            0.1,
            1.0 / 3.0,
            -1.234_567_890_123_456_7e-12,
            359.999_999_999_9,
        ] {
            let back: f64 = format_real(v).parse().unwrap();
            assert_eq!(back.to_bits(), v.to_bits(), "{v}");
        }
    }

    #[test]
    fn cards_are_80_chars_and_round_trip() {
        let w = TanWcs::new(
            SkyPos::from_degrees(83.82, -5.39).unwrap(),
            PixelPos::new(1023.5, 767.5),
            [[-7.188e-4, 4.15e-4], [4.15e-4, 7.188e-4]],
        )
        .unwrap();
        let cards = w.to_cards();
        assert!(cards.iter().all(|c| c.len() == 80));
        let back =
            TanWcs::from_header(&Header::parse_cards(cards.iter().map(String::as_str))).unwrap();
        assert_eq!(back.cd(), w.cd());
        assert_eq!(back.crpix(), w.crpix());
        assert!(back.crval().separation(w.crval()) < 1e-15);
    }

    #[test]
    fn rejects_mixed_cd_pc_and_sip() {
        let base = [
            "CTYPE1  = 'RA---TAN'",
            "CTYPE2  = 'DEC--TAN'",
            "CRVAL1  = 10.0",
            "CRVAL2  = 20.0",
            "CRPIX1  = 1.0",
            "CRPIX2  = 1.0",
        ];
        let mut mixed: Vec<&str> = base.to_vec();
        mixed.extend(["CD1_1   = 1E-3", "CD2_2   = 1E-3", "PC1_1   = 1.0"]);
        assert!(matches!(
            TanWcs::from_header(&Header::parse_cards(mixed)),
            Err(WcsError::Unsupported(_))
        ));
        let mut sip: Vec<&str> = base.to_vec();
        sip[0] = "CTYPE1  = 'RA---TAN-SIP'";
        assert!(matches!(
            TanWcs::from_header(&Header::parse_cards(sip)),
            Err(WcsError::Unsupported(_))
        ));
        let missing: Vec<&str> = base[..5].to_vec();
        assert_eq!(
            TanWcs::from_header(&Header::parse_cards(missing)),
            Err(WcsError::MissingKey("CRPIX2".into()))
        );
    }

    #[test]
    fn absent_matrix_terms_default_and_malformed_terms_fail() {
        let base = [
            "CTYPE1  = 'RA---TAN'",
            "CTYPE2  = 'DEC--TAN'",
            "CRVAL1  = 10.0",
            "CRVAL2  = 20.0",
            "CRPIX1  = 1.0",
            "CRPIX2  = 1.0",
            "CD1_1   = -1E-3",
            "CD2_2   = 1E-3",
        ];
        let read = |extra| TanWcs::from_header(&Header::parse_cards(base.into_iter().chain(extra)));
        assert_eq!(read(None).unwrap().cd(), [[-1e-3, 0.0], [0.0, 1e-3]]);
        assert!(
            matches!(read(Some("CD1_2   = broken")), Err(WcsError::InvalidValue { key, .. }) if key == "CD1_2")
        );
        assert!(
            matches!(read(Some("CUNIT1  = 123")), Err(WcsError::InvalidValue { key, .. }) if key == "CUNIT1")
        );
    }
}
