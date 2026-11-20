#![forbid(unsafe_code)]

mod image;
pub mod instrument;
mod priors;
mod sources;

use std::{fmt, path::Path};

use fayde_wcs::{Header, Value};

pub use image::{image_from_file, read_image, read_image_header, Image};
pub use instrument::{parse_parity, InstrumentPriors, InstrumentProfile, RollRule};
pub use priors::{
    header_priors, header_priors_for_size, mjd_from_iso, parse_sexagesimal, HeaderPriors,
};
pub use sources::{
    file_kind, read_sources, write_sources, write_sources_trailed, FileKind, SourceList,
};

fn pixel_dimension(value: f64) -> Option<u32> {
    (value.fract() == 0.0 && (1.0..=f64::from(u32::MAX)).contains(&value)).then_some(value as u32)
}

const BLOCK: usize = 2880;
const CARD: usize = 80;

#[derive(Debug)]
pub enum FitsError {
    Io(std::io::Error),
    Invalid(String),
}

impl fmt::Display for FitsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(e) => write!(f, "{e}"),
            Self::Invalid(s) => write!(f, "invalid FITS: {s}"),
        }
    }
}

impl std::error::Error for FitsError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io(e) => Some(e),
            _ => None,
        }
    }
}

impl From<std::io::Error> for FitsError {
    fn from(e: std::io::Error) -> Self {
        Self::Io(e)
    }
}

#[derive(Debug, Clone)]
pub struct Hdu {
    pub header: Header,
    pub cards: Vec<String>,
    pub data: std::ops::Range<usize>,
}

impl Hdu {
    pub fn xtension(&self) -> Option<&str> {
        self.header.get("XTENSION").and_then(|v| v.as_str())
    }

    fn number(&self, key: &str) -> Result<Option<f64>, FitsError> {
        self.header
            .number(key)
            .map_err(|e| FitsError::Invalid(e.to_string()))
    }

    fn int(&self, key: &str) -> Result<Option<i64>, FitsError> {
        match self
            .header
            .get_checked(key)
            .map_err(|e| FitsError::Invalid(e.to_string()))?
        {
            Some(Value::Int(value)) => Ok(Some(*value)),
            None => Ok(None),
            Some(_) => Err(FitsError::Invalid(format!("{key} must be an integer"))),
        }
    }

    fn dimension(&self, key: &str, default: usize) -> Result<usize, FitsError> {
        match self.int(key)? {
            Some(value) => {
                usize::try_from(value).map_err(|_| FitsError::Invalid(format!("{key} = {value}")))
            }
            None => Ok(default),
        }
    }
}

#[derive(Debug, Clone)]
pub struct FitsFile {
    bytes: Vec<u8>,
    hdus: Vec<Hdu>,
}

impl FitsFile {
    pub fn open(path: impl AsRef<Path>) -> Result<Self, FitsError> {
        Self::from_bytes(std::fs::read(path)?)
    }

    pub fn from_bytes(bytes: Vec<u8>) -> Result<Self, FitsError> {
        let mut hdus = Vec::new();
        let mut at = 0usize;
        while at + BLOCK <= bytes.len() {
            let mut cards = Vec::new();
            let mut end = false;
            let mut pos = at;
            while !end {
                if pos + BLOCK > bytes.len() {
                    return Err(FitsError::Invalid("header without END".into()));
                }
                for c in bytes[pos..pos + BLOCK].chunks(CARD) {
                    let s: String = c
                        .iter()
                        .map(|&b| if b.is_ascii() { b as char } else { '?' })
                        .collect();
                    if s.starts_with("END") && s[3..].trim().is_empty() {
                        end = true;
                        break;
                    }
                    cards.push(s);
                }
                pos += BLOCK;
            }
            if hdus.is_empty() && !cards.first().is_some_and(|c| c.starts_with("SIMPLE")) {
                return Err(FitsError::Invalid("first card is not SIMPLE".into()));
            }
            let header = Header::parse_cards(cards.iter().map(String::as_str));
            let hdu = Hdu {
                header,
                cards,
                data: pos..pos,
            };
            let bitpix = hdu
                .int("BITPIX")?
                .ok_or_else(|| FitsError::Invalid("no BITPIX".into()))?;
            if ![8, 16, 32, 64, -32, -64].contains(&bitpix) {
                return Err(FitsError::Invalid(format!("BITPIX = {bitpix}")));
            }
            let naxis = hdu.int("NAXIS")?.unwrap_or(0);
            if !(0..=999).contains(&naxis) {
                return Err(FitsError::Invalid(format!("NAXIS = {naxis}")));
            }
            let too_big = || FitsError::Invalid("data size overflows".into());
            let mut n: usize = if naxis == 0 { 0 } else { 1 };
            for k in 1..=naxis {
                n = n
                    .checked_mul(hdu.dimension(&format!("NAXIS{k}"), 0)?)
                    .ok_or_else(too_big)?;
            }
            let pcount = hdu.dimension("PCOUNT", 0)?;
            let gcount = hdu.dimension("GCOUNT", 1)?.max(1);
            let len = (bitpix.unsigned_abs() as usize / 8)
                .checked_mul(gcount)
                .and_then(|v| v.checked_mul(pcount.checked_add(n)?))
                .ok_or_else(too_big)?;
            let end = pos.checked_add(len).ok_or_else(too_big)?;
            if end > bytes.len() {
                return Err(FitsError::Invalid("data runs past end of file".into()));
            }
            hdus.push(Hdu {
                data: pos..end,
                ..hdu
            });
            at = len
                .div_ceil(BLOCK)
                .checked_mul(BLOCK)
                .and_then(|v| v.checked_add(pos))
                .ok_or_else(too_big)?;
        }
        if hdus.is_empty() {
            return Err(FitsError::Invalid("no header".into()));
        }
        Ok(Self { bytes, hdus })
    }

    pub fn hdus(&self) -> &[Hdu] {
        &self.hdus
    }

    pub fn primary(&self) -> &Hdu {
        &self.hdus[0]
    }

    pub fn data(&self, hdu: &Hdu) -> &[u8] {
        &self.bytes[hdu.data.clone()]
    }
}

pub(crate) fn pad_card(s: &str) -> String {
    let mut c: String = s.chars().take(CARD).collect();
    while c.len() < CARD {
        c.push(' ');
    }
    c
}

pub fn write_wcs_file(path: impl AsRef<Path>, cards: &[String]) -> Result<(), FitsError> {
    let path = path.as_ref();
    if cards
        .iter()
        .any(|card| !card.bytes().all(|b| (b' '..=b'~').contains(&b)))
    {
        return Err(FitsError::Invalid(
            "header cards must contain printable ASCII".into(),
        ));
    }
    let mut text = String::new();
    for c in [
        "SIMPLE  =                    T / conforms to FITS standard",
        "BITPIX  =                    8 / no data",
        "NAXIS   =                    0 / no data",
        "EXTEND  =                    T",
    ] {
        text.push_str(&pad_card(c));
    }
    for c in cards {
        text.push_str(&pad_card(c));
    }
    text.push_str(&pad_card("END"));
    while text.len() % BLOCK != 0 {
        text.push(' ');
    }
    fayde_core::fsutil::write_atomic(path, text.as_bytes()).map_err(FitsError::from)
}

pub fn comment_cards(text: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut line = String::new();
    for word in text.split_whitespace() {
        if line.len() + word.len() + 1 > 72 && !line.is_empty() {
            out.push(format!("COMMENT {line}"));
            line.clear();
        }
        if !line.is_empty() {
            line.push(' ');
        }
        for c in word.chars() {
            if line.len() == 72 {
                out.push(format!("COMMENT {line}"));
                line.clear();
            }
            line.push(if c.is_ascii_graphic() { c } else { '?' });
        }
    }
    if !line.is_empty() {
        out.push(format!("COMMENT {line}"));
    }
    out
}
