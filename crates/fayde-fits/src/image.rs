//! 2-D images from FITS files (FITS Standard 4.0 §5.2 data arrays, §4.4.2.5 `BSCALE`/`BZERO`/
//! `BLANK`).

use std::path::Path;

use fayde_wcs::Header;

use crate::{FitsError, FitsFile, Hdu};

/// A 2-D image as `f32`, row by row: `data[y * width + x]`, where `(x, y)` = `(0, 0)` is the
/// first pixel of the file (FITS pixel (1, 1)), the same zero-based convention as the rest of
/// Fayde. Blank pixels (`BLANK` for integer images, NaN for float images) are NaN.
#[derive(Debug, Clone)]
pub struct Image {
    /// Columns (`NAXIS1`).
    pub width: usize,
    /// Rows (`NAXIS2`).
    pub height: usize,
    /// Pixel values after `BSCALE`/`BZERO`.
    pub data: Vec<f32>,
    /// The image's header, then the primary header's cards when the image is in an extension
    /// (the image's own cards win: lookups use the first card of a keyword).
    pub header: Header,
    /// Index of the HDU the image came from.
    pub hdu: usize,
}

impl Image {
    /// Pixel at `(x, y)` (zero-based), or `None` outside the image.
    pub fn get(&self, x: usize, y: usize) -> Option<f32> {
        (x < self.width && y < self.height).then(|| self.data[y * self.width + x])
    }
}

/// Read the first 2-D image of a FITS file: the primary HDU if it holds one, else the first
/// `IMAGE` extension. Axes beyond the second must have length 1. Tile-compressed (fpack)
/// images are reported as unsupported.
pub fn read_image(path: impl AsRef<Path>) -> Result<Image, FitsError> {
    image_from_file(&FitsFile::open(path)?)
}

/// [`read_image`] on a file already read.
pub fn image_from_file(f: &FitsFile) -> Result<Image, FitsError> {
    let bad = |s: String| Err(FitsError::Invalid(s));
    for (i, hdu) in f.hdus().iter().enumerate() {
        let xt = hdu.xtension().map(|s| s.trim().to_ascii_uppercase());
        if xt.as_deref() == Some("BINTABLE") && truthy(hdu, "ZIMAGE") {
            return bad(format!(
                "HDU {i} is a tile-compressed (fpack) image, which is not supported: \
                 decompress it first (e.g. `funpack`)"
            ));
        }
        if !(i == 0 || xt.as_deref() == Some("IMAGE")) {
            continue;
        }
        let naxis = int(hdu, "NAXIS").unwrap_or(0);
        if naxis < 2 {
            continue;
        }
        let dim = |k: usize| int(hdu, &format!("NAXIS{k}")).unwrap_or(0);
        let (w, h) = (dim(1), dim(2));
        if w <= 0 || h <= 0 {
            continue;
        }
        if (3..=naxis as usize).any(|k| dim(k) != 1) {
            return bad(format!(
                "HDU {i} has {naxis} axes: only 2-D images (or cubes with one plane) are read"
            ));
        }
        let (width, height) = (w as usize, h as usize);
        let n = width
            .checked_mul(height)
            .ok_or_else(|| FitsError::Invalid("image size overflows".into()))?;
        let data = decode(hdu, f.data(hdu), n)?;
        let mut cards: Vec<&str> = hdu.cards.iter().map(String::as_str).collect();
        if i > 0 {
            cards.extend(f.primary().cards.iter().map(String::as_str));
        }
        return Ok(Image {
            width,
            height,
            data,
            header: Header::parse_cards(cards),
            hdu: i,
        });
    }
    bad("no 2-D image in the file".into())
}

fn int(hdu: &Hdu, key: &str) -> Option<i64> {
    hdu.header
        .get(key)
        .and_then(|v| v.as_f64())
        .filter(|v| v.is_finite())
        .map(|v| v as i64)
}

fn truthy(hdu: &Hdu, key: &str) -> bool {
    matches!(
        hdu.header.get(key),
        Some(fayde_wcs::header::Value::Bool(true))
    )
}

/// Big-endian pixels → `f32` with `BSCALE`/`BZERO`; integer `BLANK` → NaN.
fn decode(hdu: &Hdu, bytes: &[u8], n: usize) -> Result<Vec<f32>, FitsError> {
    let bitpix = int(hdu, "BITPIX").unwrap_or(0);
    let size = (bitpix.unsigned_abs() / 8) as usize;
    if size == 0 || bytes.len() < n.saturating_mul(size) {
        return Err(FitsError::Invalid(
            "image data shorter than NAXIS1 × NAXIS2".into(),
        ));
    }
    let num = |k: &str, d: f64| -> Result<f64, FitsError> {
        match hdu.header.get(k) {
            None => Ok(d),
            Some(v) => v
                .as_f64()
                .filter(|x| x.is_finite())
                .ok_or_else(|| FitsError::Invalid(format!("{k} is not a finite number"))),
        }
    };
    let (scale, zero) = (num("BSCALE", 1.0)?, num("BZERO", 0.0)?);
    // BLANK is matched against the raw integer exactly (a 64-bit value through f64 would lose
    // its low bits and mask the wrong pixels).
    let blank = match (bitpix > 0, hdu.header.get("BLANK")) {
        (true, Some(fayde_wcs::header::Value::Int(b))) => Some(*b),
        _ => None,
    };
    let px = |raw: f64, int_raw: Option<i64>| -> f32 {
        if int_raw.is_some() && int_raw == blank {
            f32::NAN
        } else {
            (zero + scale * raw) as f32
        }
    };
    let b = &bytes[..n * size];
    Ok(match bitpix {
        8 => b
            .iter()
            .map(|&v| px(f64::from(v), Some(i64::from(v))))
            .collect(),
        16 => b
            .chunks_exact(2)
            .map(|c| {
                let v = i16::from_be_bytes([c[0], c[1]]);
                px(f64::from(v), Some(i64::from(v)))
            })
            .collect(),
        32 => b
            .chunks_exact(4)
            .map(|c| {
                let v = i32::from_be_bytes([c[0], c[1], c[2], c[3]]);
                px(f64::from(v), Some(i64::from(v)))
            })
            .collect(),
        64 => b
            .chunks_exact(8)
            .map(|c| {
                let v = i64::from_be_bytes(c.try_into().expect("8 bytes"));
                px(v as f64, Some(v))
            })
            .collect(),
        -32 => b
            .chunks_exact(4)
            .map(|c| {
                px(
                    f64::from(f32::from_be_bytes([c[0], c[1], c[2], c[3]])),
                    None,
                )
            })
            .collect(),
        -64 => b
            .chunks_exact(8)
            .map(|c| px(f64::from_be_bytes(c.try_into().expect("8 bytes")), None))
            .collect(),
        _ => return Err(FitsError::Invalid(format!("BITPIX = {bitpix}"))),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A FITS file: primary header `cards`, data `bytes`, padded to 2880-byte blocks.
    pub(crate) fn fits(cards: &[&str], data: &[u8]) -> Vec<u8> {
        let mut out = Vec::new();
        for c in cards.iter().chain(std::iter::once(&"END")) {
            let mut card = c.to_string();
            card.truncate(80);
            out.extend(format!("{card:<80}").bytes());
        }
        out.resize(out.len().div_ceil(2880) * 2880, b' ');
        out.extend_from_slice(data);
        out.resize(out.len().div_ceil(2880) * 2880, 0);
        out
    }

    #[test]
    fn integer_image_with_scaling_and_blank() {
        // 3×2, BITPIX 16, BZERO 32768 (unsigned 16-bit convention), BLANK -1.
        let raw: [i16; 6] = [-32768, 0, 1, -1, 32767, 100];
        let data: Vec<u8> = raw.iter().flat_map(|v| v.to_be_bytes()).collect();
        let bytes = fits(
            &[
                "SIMPLE  =                    T",
                "BITPIX  =                   16",
                "NAXIS   =                    2",
                "NAXIS1  =                    3",
                "NAXIS2  =                    2",
                "BZERO   =                32768",
                "BLANK   =                   -1",
            ],
            &data,
        );
        let img = image_from_file(&FitsFile::from_bytes(bytes).unwrap()).unwrap();
        assert_eq!((img.width, img.height), (3, 2));
        assert_eq!(img.get(0, 0), Some(0.0));
        assert_eq!(img.get(1, 0), Some(32768.0));
        assert_eq!(img.get(2, 0), Some(32769.0));
        assert!(img.get(0, 1).unwrap().is_nan());
        assert_eq!(img.get(1, 1), Some(65535.0));
        assert_eq!(img.get(3, 0), None);
    }

    #[test]
    fn a_64_bit_blank_is_matched_exactly() {
        // 2^53 + 1 is not representable as f64: through f64 it would equal 2^53.
        let blank: i64 = 9_007_199_254_740_993;
        let raw = [blank, blank - 1];
        let data: Vec<u8> = raw.iter().flat_map(|v| v.to_be_bytes()).collect();
        let bytes = fits(
            &[
                "SIMPLE  =                    T",
                "BITPIX  =                   64",
                "NAXIS   =                    2",
                "NAXIS1  =                    2",
                "NAXIS2  =                    1",
                "BLANK   =     9007199254740993",
            ],
            &data,
        );
        let img = image_from_file(&FitsFile::from_bytes(bytes).unwrap()).unwrap();
        assert!(img.get(0, 0).unwrap().is_nan());
        assert!(img.get(1, 0).unwrap().is_finite());
    }

    #[test]
    fn float_image_and_unsupported_layouts() {
        let data: Vec<u8> = [1.5f32, f32::NAN]
            .iter()
            .flat_map(|v| v.to_be_bytes())
            .collect();
        let ok = fits(
            &[
                "SIMPLE  =                    T",
                "BITPIX  =                  -32",
                "NAXIS   =                    3",
                "NAXIS1  =                    2",
                "NAXIS2  =                    1",
                "NAXIS3  =                    1",
            ],
            &data,
        );
        let img = image_from_file(&FitsFile::from_bytes(ok).unwrap()).unwrap();
        assert_eq!(img.get(0, 0), Some(1.5));
        assert!(img.get(1, 0).unwrap().is_nan());

        let cube = fits(
            &[
                "SIMPLE  =                    T",
                "BITPIX  =                  -32",
                "NAXIS   =                    3",
                "NAXIS1  =                    1",
                "NAXIS2  =                    1",
                "NAXIS3  =                    2",
            ],
            &data,
        );
        let e = image_from_file(&FitsFile::from_bytes(cube).unwrap()).unwrap_err();
        assert!(e.to_string().contains("2-D"), "{e}");

        let empty = fits(
            &[
                "SIMPLE  =                    T",
                "BITPIX  =                    8",
                "NAXIS   =                    0",
            ],
            &[],
        );
        assert!(image_from_file(&FitsFile::from_bytes(empty).unwrap()).is_err());
    }
}
