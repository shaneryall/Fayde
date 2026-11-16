//! FITS input and output for Fayde (FITS Standard 4.0, Pence et al. 2010 for 3.0).
//!
//! - [`FitsFile`]: header-and-data units of a FITS file, read without external libraries.
//! - [`read_image`]: a 2-D image as `f32`, for source extraction.
//! - [`write_sources`]: write sources as such a table (X, Y, FLUX, SIGMA, P_STAR, TRUNCATED).
//! - [`read_sources`]: a source list from a FITS binary table (Fayde's format, `frame.sources.fits`: any table
//!   with `X`/`Y` columns), converted from FITS 1-based to Fayde's zero-based pixels.
//! - [`header_priors`]: pointing, pixel scale and observation epoch from an image header, with
//!   strict validation.
//! - [`write_wcs_file`]: a header-only FITS file holding a solved WCS, like `solve-field`'s
//!   `.wcs` output.

#![forbid(unsafe_code)]

mod image;
pub mod instrument;
mod priors;
mod sources;

use std::{fmt, path::Path};

use fayde_wcs::Header;

pub use image::{image_from_file, read_image, Image};
pub use instrument::{parse_parity, InstrumentPriors, InstrumentProfile, RollRule};
pub use priors::{header_priors, mjd_from_iso, parse_sexagesimal, HeaderPriors};
pub use sources::{file_kind, read_sources, write_sources, FileKind, SourceList};

const BLOCK: usize = 2880;
const CARD: usize = 80;

/// Errors reading or writing FITS.
#[derive(Debug)]
pub enum FitsError {
    /// I/O failure.
    Io(std::io::Error),
    /// The file is not valid FITS or lacks what was asked for.
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

impl std::error::Error for FitsError {}

impl From<std::io::Error> for FitsError {
    fn from(e: std::io::Error) -> Self {
        Self::Io(e)
    }
}

/// One header-and-data unit.
#[derive(Debug, Clone)]
pub struct Hdu {
    /// Parsed header values.
    pub header: Header,
    /// Raw card images, in order.
    pub cards: Vec<String>,
    /// Byte range of the data in the file.
    pub data: std::ops::Range<usize>,
}

impl Hdu {
    /// `XTENSION` value (`None` for the primary HDU).
    pub fn xtension(&self) -> Option<&str> {
        self.header.get("XTENSION").and_then(|v| v.as_str())
    }

    fn int(&self, key: &str) -> Option<i64> {
        self.header
            .get(key)
            .and_then(|v| v.as_f64())
            .map(|v| v as i64)
    }
}

/// A FITS file held in memory.
#[derive(Debug, Clone)]
pub struct FitsFile {
    bytes: Vec<u8>,
    hdus: Vec<Hdu>,
}

impl FitsFile {
    /// Read and index a FITS file.
    pub fn open(path: impl AsRef<Path>) -> Result<Self, FitsError> {
        Self::from_bytes(std::fs::read(path)?)
    }

    /// Index FITS bytes.
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
                .int("BITPIX")
                .ok_or_else(|| FitsError::Invalid("no BITPIX".into()))?;
            if ![8, 16, 32, 64, -32, -64].contains(&bitpix) {
                return Err(FitsError::Invalid(format!("BITPIX = {bitpix}")));
            }
            let naxis = hdu.int("NAXIS").unwrap_or(0);
            if !(0..=999).contains(&naxis) {
                return Err(FitsError::Invalid(format!("NAXIS = {naxis}")));
            }
            // Sizes come from the file: every step is checked, so a hostile header gives an
            // error rather than an overflow and an out-of-range slice.
            let too_big = || FitsError::Invalid("data size overflows".into());
            let dim = |k: &str, default: i64| -> Result<usize, FitsError> {
                let v = hdu.int(k).unwrap_or(default);
                usize::try_from(v).map_err(|_| FitsError::Invalid(format!("{k} = {v}")))
            };
            let mut n: usize = if naxis == 0 { 0 } else { 1 };
            for k in 1..=naxis {
                n = n
                    .checked_mul(dim(&format!("NAXIS{k}"), 0)?)
                    .ok_or_else(too_big)?;
            }
            let pcount = dim("PCOUNT", 0)?;
            let gcount = dim("GCOUNT", 1)?.max(1);
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

    /// All HDUs.
    pub fn hdus(&self) -> &[Hdu] {
        &self.hdus
    }

    /// Primary header.
    pub fn primary(&self) -> &Hdu {
        &self.hdus[0]
    }

    /// Data bytes of an HDU.
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

/// Write a header-only FITS file: `SIMPLE`, `BITPIX = 8`, `NAXIS = 0`, then `cards` (80-character
/// images, e.g. from [`fayde_wcs::Wcs::to_cards`]), then `END`. The file is written to a temporary
/// name in the same directory and renamed into place.
pub fn write_wcs_file(path: impl AsRef<Path>, cards: &[String]) -> Result<(), FitsError> {
    let path = path.as_ref();
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
    let (mut f, tmp) = fayde_core::fsutil::create_temp_beside(path)?;
    let write = || -> std::io::Result<()> {
        use std::io::Write;
        f.write_all(text.as_bytes())?;
        f.sync_all()?;
        drop(f);
        std::fs::rename(&tmp, path)
    };
    // A failed write must not leave the temporary file behind.
    write().map_err(|e| {
        let _ = std::fs::remove_file(&tmp);
        FitsError::from(e)
    })
}

/// A comment card (`COMMENT text`), wrapped at 72 characters.
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
        line.push_str(word);
    }
    if !line.is_empty() {
        out.push(format!("COMMENT {line}"));
    }
    out
}
