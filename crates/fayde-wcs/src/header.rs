//! Minimal FITS header cards: parse 80-character card images and write TAN WCS cards.
//!
//! Card syntax follows the FITS Standard 3.0 (Pence et al. 2010, A&A 524, A42).
//!
//! This is not a general FITS library. It reads the keyword/value cards a WCS needs and writes
//! standard-conforming cards that Astropy, wcslib and astrometry.net all read.

use fayde_core::{PixelPos, Portable, SkyPos};

use crate::{TanWcs, WcsError};

/// A parsed header value.
#[derive(Debug, Clone, PartialEq)]
pub enum Value {
    /// Quoted string, trailing blanks removed.
    Str(String),
    /// Logical `T` / `F`.
    Bool(bool),
    /// Integer.
    Int(i64),
    /// Real (`E` or `D` exponent).
    Float(f64),
}

impl Value {
    /// Numeric value as `f64` (integers are widened).
    pub fn as_f64(&self) -> Option<f64> {
        match *self {
            Self::Int(i) => Some(i as f64),
            Self::Float(f) => Some(f),
            _ => None,
        }
    }

    /// String value.
    pub fn as_str(&self) -> Option<&str> {
        match self {
            Self::Str(s) => Some(s),
            _ => None,
        }
    }
}

/// Keyword → value pairs in card order. A keyword that appears more than once is looked up by
/// its first card, as cfitsio and Astropy do.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Header {
    cards: Vec<(String, Value)>,
}

impl Header {
    /// Parse card images (one card per item; each up to 80 characters). Cards without a value
    /// indicator (`COMMENT`, `HISTORY`, blank, `END`) and cards with unparseable values are
    /// skipped; lookups of those keys then report them missing.
    pub fn parse_cards<'a>(cards: impl IntoIterator<Item = &'a str>) -> Self {
        let mut out = Vec::new();
        for card in cards {
            if let Some((key, value)) = parse_card(card) {
                out.push((key, value));
            }
        }
        Self { cards: out }
    }

    /// Parse a header block: concatenated 80-character cards with no separators, as stored in a
    /// FITS file, or newline-separated cards.
    pub fn parse_block(text: &str) -> Self {
        if text.contains('\n') {
            Self::parse_cards(text.lines())
        } else {
            let chars: Vec<char> = text.chars().collect();
            let cards: Vec<String> = chars.chunks(80).map(|c| c.iter().collect()).collect();
            Self::parse_cards(cards.iter().map(String::as_str))
        }
    }

    /// Insert or replace a value.
    pub fn set(&mut self, key: &str, value: Value) {
        let key = key.to_ascii_uppercase();
        match self.cards.iter_mut().find(|(k, _)| *k == key) {
            Some(card) => card.1 = value,
            None => self.cards.push((key, value)),
        }
    }

    /// Look up a keyword (case-insensitive).
    pub fn get(&self, key: &str) -> Option<&Value> {
        let key = key.to_ascii_uppercase();
        self.cards.iter().find(|(k, _)| *k == key).map(|(_, v)| v)
    }

    /// Whether a keyword is present.
    pub fn contains(&self, key: &str) -> bool {
        self.get(key).is_some()
    }

    fn f64_or(&self, key: &str, default: f64) -> Result<f64, WcsError> {
        match self.get(key) {
            None => Ok(default),
            Some(v) => v.as_f64().ok_or_else(|| WcsError::InvalidValue {
                key: key.to_string(),
                reason: format!("expected a number, got {v:?}"),
            }),
        }
    }

    fn f64_required(&self, key: &str) -> Result<f64, WcsError> {
        match self.get(key) {
            None => Err(WcsError::MissingKey(key.to_string())),
            Some(_) => self.f64_or(key, f64::NAN),
        }
    }

    fn str_trimmed(&self, key: &str) -> Option<String> {
        self.get(key)
            .and_then(Value::as_str)
            .map(|s| s.trim().to_string())
    }
}

fn parse_card(card: &str) -> Option<(String, Value)> {
    let card: String = card.chars().take(80).collect();
    // `get` returns None where a byte offset falls inside a multi-byte character: FITS cards
    // are ASCII, so such a card is malformed, not a reason to panic.
    if card.get(8..10) != Some("= ") {
        return None;
    }
    let key = card.get(..8)?.trim_end().to_ascii_uppercase();
    if key.is_empty() {
        return None;
    }
    parse_value(card.get(10..)?).map(|v| (key, v))
}

fn parse_value(field: &str) -> Option<Value> {
    let s = field.trim_start();
    if let Some(rest) = s.strip_prefix('\'') {
        // Quoted string; '' is an escaped quote.
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
        return None; // unterminated
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
    // A non-finite number (`NaN`, `inf`, an exponent past the f64 range) is kept as text: a
    // numeric read of that keyword then fails by name instead of reaching the arithmetic.
    match token.replace(['D', 'd'], "E").parse::<f64>() {
        Ok(v) if v.is_finite() => Some(Value::Float(v)),
        Ok(_) => Some(Value::Str(token.to_string())),
        Err(_) => None,
    }
}

/// The celestial reference frame must be ICRS or FK5 J2000 (FITS Standard 4.0 §8.1):
/// `RADESYS` (or the older `RADECSYS`) `ICRS`, or `FK5` with `EQUINOX` (or the older `EPOCH`)
/// absent or 2000. Without `RADESYS`, no equinox means ICRS, an equinox before 1984 means FK4
/// (rejected) and any other must be 2000. Other frames and equinoxes would need a precession
/// or frame transformation, which Fayde does not apply, so they are rejected rather than
/// silently read as ICRS.
pub fn check_reference_frame(h: &Header) -> Result<(), WcsError> {
    let sys = h
        .str_trimmed("RADESYS")
        .or_else(|| h.str_trimmed("RADECSYS"))
        .map(|s| s.to_ascii_uppercase());
    let equinox = match h.get("EQUINOX").or_else(|| h.get("EPOCH")) {
        None => None,
        Some(v) => Some(
            v.as_f64()
                .or_else(|| {
                    v.as_str()
                        .and_then(|s| s.trim().trim_start_matches(['J', 'j']).parse().ok())
                })
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

/// Whether a `CUNITi` value means degrees. FITS requires `deg`; `degree(s)` is common in real
/// headers and unambiguous, so it is accepted too.
pub fn is_degree_unit(u: &str) -> bool {
    ["deg", "degree", "degrees"]
        .iter()
        .any(|d| u.trim().eq_ignore_ascii_case(d))
}

/// Format a real value losslessly in FITS form: shortest round-trip mantissa, upper-case `E`,
/// signed two-digit (or longer) exponent, e.g. `-8.3E-04`.
pub fn format_real(v: f64) -> String {
    let s = format!("{v:e}"); // e.g. "-8.3e-4", "1.5e2", "0e0"
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

/// A string-valued card (crate-internal helper for other modules).
pub(crate) fn str_card_pub(key: &str, value: &str, comment: &str) -> String {
    str_card(key, value, comment)
}

/// A numeric card with no comment.
pub(crate) fn num_card(key: &str, value: &str) -> String {
    card(key, value, "")
}

fn str_card(key: &str, value: &str, comment: &str) -> String {
    // FITS: string values start at column 11 and are at least 8 characters inside the quotes.
    let quoted = format!("'{:<8}'", value.replace('\'', "''"));
    card(key, &format!("{quoted:<20}"), comment)
}

impl TanWcs {
    /// Read a TAN WCS from FITS header values.
    ///
    /// Accepts `CDi_j`, or `PCi_j` with optional `CDELTi` (both default to identity/1), or
    /// `CDELTi` with the older `CROTA2` rotation (`CROTA1` if only that is given). `PCi_j` takes
    /// precedence over `CROTAn`, as in wcslib. A header carrying both CD and PC keywords is
    /// rejected as ambiguous. `CUNITi` must be
    /// `deg` when present; `LONPOLE` must be 180 when present; `RADESYS` must be `ICRS` or
    /// `FK5` with equinox J2000 (treated as ICRS; they differ by tens of mas), following the
    /// FITS defaults when `RADESYS` is absent (see [`check_reference_frame`]). `-SIP` and other projections are rejected.
    pub fn from_header(h: &Header) -> Result<Self, WcsError> {
        for (key, want) in [("CTYPE1", "RA---TAN"), ("CTYPE2", "DEC--TAN")] {
            let got = h
                .str_trimmed(key)
                .ok_or_else(|| WcsError::MissingKey(key.to_string()))?;
            if got != want {
                return Err(WcsError::Unsupported(format!("{key} = '{got}'")));
            }
        }
        for key in ["CUNIT1", "CUNIT2"] {
            if let Some(u) = h.str_trimmed(key) {
                if !is_degree_unit(&u) {
                    return Err(WcsError::Unsupported(format!("{key} = '{u}'")));
                }
            }
        }
        check_reference_frame(h)?;
        let lonpole = h.f64_or("LONPOLE", 180.0)?;
        if (lonpole - 180.0).abs() > 1e-12 {
            return Err(WcsError::Unsupported(format!("LONPOLE = {lonpole}")));
        }

        const CD: [&str; 4] = ["CD1_1", "CD1_2", "CD2_1", "CD2_2"];
        const PC: [&str; 4] = ["PC1_1", "PC1_2", "PC2_1", "PC2_2"];
        let has_cd = CD.iter().any(|k| h.contains(k));
        let has_pc = PC.iter().any(|k| h.contains(k));
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
                // AIPS convention (Calabretta & Greisen 2002, §6.1): a rotation ρ = CROTA2 of
                // the image axes, CD = [[c1·cos ρ, −c2·sin ρ], [c1·sin ρ, c2·cos ρ]].
                let rho = h.f64_or("CROTA2", h.f64_or("CROTA1", 0.0)?)?;
                let (s, c) = rho.to_radians().psin_cos();
                [[c1 * c, -c2 * s], [c1 * s, c2 * c]]
            }
        };

        let ra = h.f64_required("CRVAL1")?;
        let dec = h.f64_required("CRVAL2")?;
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

    /// FITS card images (80 characters each) for this WCS, in `CDi_j` form with a 1-based
    /// `CRPIX`. Values are written losslessly.
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
        // Older headers: RADECSYS and EPOCH.
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
}
