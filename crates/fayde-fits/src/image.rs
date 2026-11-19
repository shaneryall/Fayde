use std::path::Path;

use fayde_wcs::Header;

use crate::{FitsError, FitsFile, Hdu};

#[derive(Debug, Clone)]
pub struct Image {
    pub width: usize,
    pub height: usize,
    pub data: Vec<f32>,
    pub header: Header,
    pub hdu: usize,
}

impl Image {
    pub fn get(&self, x: usize, y: usize) -> Option<f32> {
        if x >= self.width || y >= self.height {
            return None;
        }
        self.data
            .get(y.checked_mul(self.width)?.checked_add(x)?)
            .copied()
    }
}

/// Read the first 2-D image. Extension header values take precedence over primary header values.
pub fn read_image(path: impl AsRef<Path>) -> Result<Image, FitsError> {
    image_from_file(&FitsFile::open(path)?)
}

pub fn image_from_file(f: &FitsFile) -> Result<Image, FitsError> {
    let (i, width, height) =
        image_hdu(f)?.ok_or_else(|| FitsError::Invalid("no 2-D image in the file".into()))?;
    let hdu = &f.hdus()[i];
    let n = width
        .checked_mul(height)
        .ok_or_else(|| FitsError::Invalid("image size overflows".into()))?;
    Ok(Image {
        width,
        height,
        data: decode(hdu, f.data(hdu), n)?,
        header: merged_header(f, i),
        hdu: i,
    })
}

/// Read the header selected by [`read_image`], or the primary header of a header-only file.
pub fn read_image_header(path: impl AsRef<Path>) -> Result<Header, FitsError> {
    let file = FitsFile::open(path)?;
    let i = image_hdu(&file)?.map_or(0, |(i, _, _)| i);
    Ok(merged_header(&file, i))
}

fn merged_header(file: &FitsFile, i: usize) -> Header {
    if i == 0 {
        return file.primary().header.clone();
    }
    Header::parse_cards(
        file.hdus()[i]
            .cards
            .iter()
            .chain(&file.primary().cards)
            .map(String::as_str),
    )
}

fn image_hdu(f: &FitsFile) -> Result<Option<(usize, usize, usize)>, FitsError> {
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
        let naxis = hdu.dimension("NAXIS", 0)?;
        if naxis < 2 {
            continue;
        }
        let (width, height) = (hdu.dimension("NAXIS1", 0)?, hdu.dimension("NAXIS2", 0)?);
        if width == 0 || height == 0 {
            continue;
        }
        for k in 3..=naxis {
            if hdu.dimension(&format!("NAXIS{k}"), 0)? != 1 {
                return bad(format!(
                    "HDU {i} has {naxis} axes: only 2-D images (or cubes with one plane) are read"
                ));
            }
        }
        return Ok(Some((i, width, height)));
    }
    Ok(None)
}

fn truthy(hdu: &Hdu, key: &str) -> bool {
    matches!(
        hdu.header.get(key),
        Some(fayde_wcs::header::Value::Bool(true))
    )
}

fn decode(hdu: &Hdu, bytes: &[u8], n: usize) -> Result<Vec<f32>, FitsError> {
    let bitpix = hdu
        .int("BITPIX")?
        .ok_or_else(|| FitsError::Invalid("no BITPIX".into()))?;
    let size = (bitpix.unsigned_abs() / 8) as usize;
    if size == 0 || bytes.len() < n.saturating_mul(size) {
        return Err(FitsError::Invalid(
            "image data shorter than NAXIS1 × NAXIS2".into(),
        ));
    }
    let scale = hdu.number("BSCALE")?.unwrap_or(1.0);
    let zero = hdu.number("BZERO")?.unwrap_or(0.0);
    let blank = if bitpix > 0 { hdu.int("BLANK")? } else { None };
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
    fn malformed_image_scaling_and_null_cards_are_errors() {
        for card in ["BSCALE  = broken", "BZERO   = 'words'", "BLANK   = 0.5"] {
            let bytes = fits(
                &[
                    "SIMPLE  = T",
                    "BITPIX  = 16",
                    "NAXIS   = 2",
                    "NAXIS1  = 1",
                    "NAXIS2  = 1",
                    card,
                ],
                &[0, 1],
            );
            let err = image_from_file(&FitsFile::from_bytes(bytes).unwrap()).unwrap_err();
            assert!(err.to_string().contains(card[..8].trim()), "{err}");
        }
    }

    #[test]
    fn integer_image_with_scaling_and_blank() {
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
    fn fractional_bitpix_is_rejected_before_decoding() {
        let bytes = fits(&["SIMPLE  = T", "BITPIX  = -32.5", "NAXIS   = 0"], &[]);
        let error = FitsFile::from_bytes(bytes).unwrap_err();
        assert!(error.to_string().contains("BITPIX must be an integer"));
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
