//! Source lists from FITS binary tables.

use std::path::Path;

use fayde_core::PixelConvention;
use fayde_solve::Source;

use crate::{pad_card, FitsError, FitsFile, Hdu};

/// A source list read from a FITS table.
#[derive(Debug, Clone, PartialEq)]
pub struct SourceList {
    /// Sources in Fayde's zero-based pixel convention, in table order.
    pub sources: Vec<Source>,
    /// `IMAGEW`, if present.
    pub width: Option<u32>,
    /// `IMAGEH`, if present.
    pub height: Option<u32>,
}

#[derive(Debug, Clone, Copy)]
struct Column {
    offset: usize,
    kind: u8,
    /// `TSCALn`, `TZEROn`: physical = zero + scale × stored (FITS Standard 4.0 §7.3.2).
    scale: f64,
    zero: f64,
    /// `TNULLn` of an integer column: the stored value meaning "undefined".
    null: Option<i64>,
}

impl Column {
    /// The physical value of this column in `row`; `None` for undefined or unreadable values.
    fn read(&self, row: &[u8]) -> Option<f64> {
        let at = self.offset;
        let b = |n: usize| row.get(at..at + n);
        let int = |v: i64| (self.null != Some(v)).then_some(v as f64);
        let raw = match self.kind {
            b'E' => f64::from(f32::from_be_bytes(b(4)?.try_into().ok()?)),
            b'D' => f64::from_be_bytes(b(8)?.try_into().ok()?),
            b'J' => int(i64::from(i32::from_be_bytes(b(4)?.try_into().ok()?)))?,
            b'I' => int(i64::from(i16::from_be_bytes(b(2)?.try_into().ok()?)))?,
            b'K' => int(i64::from_be_bytes(b(8)?.try_into().ok()?))?,
            b'B' => int(i64::from(*row.get(at)?))?,
            // Logical: not scaled.
            b'L' => return Some(f64::from(u8::from(*row.get(at)? == b'T'))),
            _ => return None,
        };
        let v = self.zero + self.scale * raw;
        v.is_finite().then_some(v)
    }
}

fn width_of(kind: u8) -> Option<usize> {
    Some(match kind {
        b'L' | b'B' | b'A' | b'X' => 1,
        b'I' => 2,
        b'J' | b'E' => 4,
        b'K' | b'D' | b'C' => 8,
        b'M' => 16,
        b'P' => 8,
        b'Q' => 16,
        _ => return None,
    })
}

/// Optional per-source columns of a source list.
struct Optional {
    flux: Option<Column>,
    truncated: Option<Column>,
    sigma: Option<Column>,
    p_star: Option<Column>,
}

/// Columns by upper-case name (`TTYPEn`), for scalar numeric columns.
fn columns(hdu: &Hdu) -> Result<Vec<(String, Column)>, FitsError> {
    let n = hdu
        .header
        .get("TFIELDS")
        .and_then(|v| v.as_f64())
        .unwrap_or(0.0) as usize;
    let mut out = Vec::new();
    let mut offset = 0usize;
    for i in 1..=n {
        let form = hdu
            .header
            .get(&format!("TFORM{i}"))
            .and_then(|v| v.as_str())
            .ok_or_else(|| FitsError::Invalid(format!("TFORM{i} missing")))?
            .trim()
            .to_string();
        let digits: String = form.chars().take_while(char::is_ascii_digit).collect();
        let repeat: usize = if digits.is_empty() {
            1
        } else {
            digits.parse().unwrap_or(1)
        };
        let kind = *form
            .as_bytes()
            .get(digits.len())
            .ok_or_else(|| FitsError::Invalid(format!("TFORM{i} = '{form}'")))?;
        let w = width_of(kind).ok_or_else(|| FitsError::Invalid(format!("TFORM{i} = '{form}'")))?;
        let name = hdu
            .header
            .get(&format!("TTYPE{i}"))
            .and_then(|v| v.as_str())
            .unwrap_or("")
            .trim()
            .to_ascii_uppercase();
        if repeat == 1 {
            let num = |k: &str| hdu.header.get(&format!("{k}{i}")).and_then(|v| v.as_f64());
            let scale = num("TSCAL").unwrap_or(1.0);
            let zero = num("TZERO").unwrap_or(0.0);
            if !(scale.is_finite() && zero.is_finite()) {
                return Err(FitsError::Invalid(format!("TSCAL{i}/TZERO{i} not finite")));
            }
            let null = num("TNULL").filter(|v| v.fract() == 0.0).map(|v| v as i64);
            out.push((
                name,
                Column {
                    offset,
                    kind,
                    scale,
                    zero,
                    null,
                },
            ));
        }
        offset = w
            .checked_mul(if kind == b'X' {
                repeat.div_ceil(8)
            } else {
                repeat
            })
            .and_then(|n| offset.checked_add(n))
            .ok_or_else(|| FitsError::Invalid(format!("TFORM{i}: row too long")))?;
    }
    Ok(out)
}

/// Write `sources` as a source list (`frame.sources.fits`): an empty primary HDU, then a binary
/// table `EXTNAME = 'SOURCES'` with `X`, `Y`
/// (FITS 1-based, `D`), `FLUX`, `SIGMA` and `P_STAR` (`D`, NaN when absent) and `TRUNCATED` (`L`), with
/// `IMAGEW`/`IMAGEH`. [`read_sources`] reads it back unchanged. Written to a temporary file and
/// renamed into place.
pub fn write_sources(
    path: impl AsRef<Path>,
    sources: &[Source],
    width: u32,
    height: u32,
) -> Result<(), FitsError> {
    const ROW: usize = 8 * 5 + 1;
    let mut out = Vec::new();
    let block = |out: &mut Vec<u8>, fill: u8| {
        let n = out.len().div_ceil(2880) * 2880;
        out.resize(n, fill);
    };
    for c in [
        "SIMPLE  =                    T",
        "BITPIX  =                    8",
        "NAXIS   =                    0",
        "EXTEND  =                    T",
        "END",
    ] {
        out.extend(pad_card(c).bytes());
    }
    block(&mut out, b' ');
    let mut cards = vec![
        "XTENSION= 'BINTABLE'".to_string(),
        "BITPIX  =                    8".into(),
        "NAXIS   =                    2".into(),
        format!("NAXIS1  = {ROW:>20}"),
        format!("NAXIS2  = {:>20}", sources.len()),
        "PCOUNT  =                    0".into(),
        "GCOUNT  =                    1".into(),
        "TFIELDS =                    6".into(),
    ];
    for (k, (name, form)) in [
        ("X", "D"),
        ("Y", "D"),
        ("FLUX", "D"),
        ("SIGMA", "D"),
        ("P_STAR", "D"),
        ("TRUNCATED", "L"),
    ]
    .iter()
    .enumerate()
    {
        cards.push(format!("TTYPE{}  = '{name:<8}'", k + 1));
        cards.push(format!("TFORM{}  = '{form:<8}'", k + 1));
    }
    // Say what the table is, so it is recognisable in any FITS viewer even if renamed.
    cards.push("EXTNAME = 'SOURCES'".into());
    cards.push(format!(
        "CREATOR = '{:<18}'",
        concat!("fayde ", env!("CARGO_PKG_VERSION"))
    ));
    cards
        .push("COMMENT Fayde source list: X, Y in FITS pixels (first pixel centre at 1, 1)".into());
    cards.push(format!("IMAGEW  = {width:>20}"));
    cards.push(format!("IMAGEH  = {height:>20}"));
    cards.push("END".into());
    for c in &cards {
        out.extend(pad_card(c).bytes());
    }
    block(&mut out, b' ');
    for s in sources {
        let (x, y) = PixelConvention::Fits.from_internal(fayde_core::PixelPos::new(s.x, s.y));
        for v in [
            x,
            y,
            s.flux.unwrap_or(f64::NAN),
            s.sigma_px.unwrap_or(f64::NAN),
            s.p_star.unwrap_or(f64::NAN),
        ] {
            out.extend(v.to_be_bytes());
        }
        out.push(if s.truncated { b'T' } else { b'F' });
    }
    block(&mut out, 0);
    let path = path.as_ref();
    let (mut f, tmp) = fayde_core::fsutil::create_temp_beside(path)?;
    let write = || -> std::io::Result<()> {
        use std::io::Write;
        f.write_all(&out)?;
        f.sync_all()?;
        drop(f);
        std::fs::rename(&tmp, path)
    };
    write().map_err(|e| {
        let _ = std::fs::remove_file(&tmp);
        FitsError::from(e)
    })
}

/// What a FITS file holds, judged by its contents (not its name).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FileKind {
    /// A binary table with `X` and `Y` columns: a source list.
    SourceList,
    /// A 2-D image (primary HDU or an `IMAGE` extension).
    Image,
    /// Neither.
    Other,
}

/// What `f` holds: a source list if any table has `X` and `Y` columns, else an image if any
/// HDU holds a 2-D array, else neither.
pub fn file_kind(f: &FitsFile) -> FileKind {
    let has =
        |cols: &[(String, Column)], n: &str| cols.iter().any(|(c, _)| c.eq_ignore_ascii_case(n));
    let table = f.hdus().iter().any(|h| {
        h.xtension()
            .is_some_and(|x| x.trim().eq_ignore_ascii_case("BINTABLE"))
            && columns(h).is_ok_and(|c| has(&c, "X") && has(&c, "Y"))
    });
    if table {
        return FileKind::SourceList;
    }
    let image = f.hdus().iter().enumerate().any(|(i, h)| {
        let image_hdu = i == 0
            || h.xtension()
                .is_some_and(|x| x.trim().eq_ignore_ascii_case("IMAGE"));
        let n = |k: &str| h.header.get(k).and_then(|v| v.as_f64()).unwrap_or(0.0);
        image_hdu && n("NAXIS") >= 2.0 && n("NAXIS1") > 0.0 && n("NAXIS2") > 0.0
    });
    if image {
        FileKind::Image
    } else {
        FileKind::Other
    }
}

/// Read the first binary table with `X` and `Y` columns (FITS 1-based pixels). Optional
/// columns: `FLUX`, `TRUNCATED` (logical or integer), `SIGMA` or `SIGMA_PX` (1σ centroid
/// error, pixels) and `P_STAR` (prior probability of a real star, 0–1). `IMAGEW`/`IMAGEH`
/// are taken from the table header or the primary header.
pub fn read_sources(path: impl AsRef<Path>) -> Result<SourceList, FitsError> {
    let f = FitsFile::open(path)?;
    // The first binary table with scalar X and Y columns: a file may carry other tables
    // (metadata, per-frame records) before the source list.
    let mut tables = 0;
    let mut first_err = None;
    let mut found = None;
    for hdu in f
        .hdus()
        .iter()
        .filter(|h| h.xtension().map(str::trim) == Some("BINTABLE"))
    {
        tables += 1;
        match columns(hdu) {
            Ok(cols) => {
                let find = |name: &str| cols.iter().find(|(n, _)| n == name).map(|(_, c)| *c);
                if let (Some(x), Some(y)) = (find("X"), find("Y")) {
                    let opt = Optional {
                        flux: find("FLUX"),
                        truncated: find("TRUNCATED"),
                        sigma: find("SIGMA").or_else(|| find("SIGMA_PX")),
                        p_star: find("P_STAR"),
                    };
                    found = Some((hdu, x, y, opt));
                    break;
                }
            }
            Err(e) => {
                first_err.get_or_insert(e);
            }
        }
    }
    let Some((hdu, x, y, opt)) = found else {
        return Err(match (tables, first_err) {
            (0, _) => FitsError::Invalid("no BINTABLE extension".into()),
            (_, Some(e)) => e,
            _ => FitsError::Invalid("no binary table with X and Y columns".into()),
        });
    };
    let row_len = hdu
        .header
        .get("NAXIS1")
        .and_then(|v| v.as_f64())
        .unwrap_or(0.0) as usize;
    let rows = hdu
        .header
        .get("NAXIS2")
        .and_then(|v| v.as_f64())
        .unwrap_or(0.0) as usize;
    let data = f.data(hdu);
    if row_len == 0 || row_len.checked_mul(rows).is_none_or(|n| data.len() < n) {
        return Err(FitsError::Invalid(
            "table data shorter than NAXIS1 × NAXIS2".into(),
        ));
    }
    let mut sources = Vec::with_capacity(rows);
    for r in 0..rows {
        let row = &data[r * row_len..(r + 1) * row_len];
        let (Some(fx), Some(fy)) = (x.read(row), y.read(row)) else {
            return Err(FitsError::Invalid(format!("row {r}: unreadable X/Y")));
        };
        let p = PixelConvention::Fits.to_internal(fx, fy);
        let mut s = Source::at(p.x, p.y);
        s.flux = opt.flux.and_then(|c| c.read(row));
        // Per-source centroid error (pixels) and prior star probability, when given; values
        // out of range are treated as absent rather than trusted.
        s.sigma_px = opt.sigma.and_then(|c| c.read(row)).filter(|v| *v > 0.0);
        s.p_star = opt
            .p_star
            .and_then(|c| c.read(row))
            .filter(|v| (0.0..=1.0).contains(v));
        // Optional TRUNCATED column (logical or integer): the source is cut by the frame edge.
        s.truncated = opt
            .truncated
            .and_then(|c| c.read(row))
            .is_some_and(|v| v != 0.0);
        sources.push(s);
    }
    let dim = |k: &str| {
        hdu.header
            .get(k)
            .or_else(|| f.primary().header.get(k))
            .and_then(|v| v.as_f64())
            .filter(|v| *v > 0.0)
            .map(|v| v as u32)
    };
    Ok(SourceList {
        sources,
        width: dim("IMAGEW"),
        height: dim("IMAGEH"),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A FITS file: an empty primary HDU, then one BINTABLE per `(cards, rows)`.
    fn fits(tables: &[(&[&str], Vec<u8>)]) -> Vec<u8> {
        let block = |cards: &[String]| {
            let mut b: Vec<u8> = cards
                .iter()
                .chain(std::iter::once(&"END".to_string()))
                .flat_map(|c| format!("{c:<80}").into_bytes())
                .collect();
            b.resize(b.len().div_ceil(2880) * 2880, b' ');
            b
        };
        let mut out = block(&[
            "SIMPLE  =                    T".into(),
            "BITPIX  =                    8".into(),
            "NAXIS   =                    0".into(),
        ]);
        for (cards, data) in tables {
            out.extend(block(
                &cards.iter().map(|c| c.to_string()).collect::<Vec<_>>(),
            ));
            let mut d = data.clone();
            d.resize(d.len().div_ceil(2880) * 2880, 0);
            out.extend(d);
        }
        out
    }

    fn read(bytes: Vec<u8>) -> Result<SourceList, FitsError> {
        // Tests run in parallel: each file needs its own name.
        static N: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
        let n = N.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        let path = std::env::temp_dir().join(format!(
            "fayde-sources-test-{}-{n}.fits",
            std::process::id()
        ));
        std::fs::write(&path, bytes).unwrap();
        let r = read_sources(&path);
        let _ = std::fs::remove_file(&path);
        r
    }

    const XY_J: &[&str] = &[
        "XTENSION= 'BINTABLE'",
        "BITPIX  =                    8",
        "NAXIS   =                    2",
        "NAXIS1  =                    8",
        "NAXIS2  =                    1",
        "PCOUNT  =                    0",
        "GCOUNT  =                    1",
        "TFIELDS =                    2",
        "TTYPE1  = 'X       '",
        "TFORM1  = '1J      '",
        "TSCAL1  =                  0.5",
        "TZERO1  =                  1.0",
        "TTYPE2  = 'Y       '",
        "TFORM2  = '1J      '",
    ];

    fn row(x: i32, y: i32) -> Vec<u8> {
        [x.to_be_bytes(), y.to_be_bytes()].concat()
    }

    #[test]
    fn column_scaling_is_applied() {
        // Stored X = 100, TSCAL 0.5, TZERO 1 → physical 51 (FITS) → 50 (zero-based).
        let xy = read(fits(&[(XY_J, row(100, 20))])).unwrap();
        assert_eq!(xy.sources[0].x, 50.0);
        assert_eq!(xy.sources[0].y, 19.0);
    }

    #[test]
    fn written_lists_read_back_unchanged() {
        let dir = std::env::temp_dir().join(format!("fayde-xyw-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("out.sources.fits");
        let mut a = Source::at(10.25, 20.5);
        a.flux = Some(1234.5);
        a.sigma_px = Some(0.05);
        a.p_star = Some(0.2);
        let mut b = Source::at(0.0, 767.0);
        b.truncated = true;
        write_sources(&path, &[a, b], 768, 768).unwrap();
        let xy = read_sources(&path).unwrap();
        assert_eq!((xy.width, xy.height), (Some(768), Some(768)));
        let f = FitsFile::open(&path).unwrap();
        assert_eq!(file_kind(&f), FileKind::SourceList);
        assert_eq!(
            f.hdus()[1].header.get("EXTNAME").and_then(|v| v.as_str()),
            Some("SOURCES")
        );
        assert_eq!(xy.sources.len(), 2);
        assert_eq!((xy.sources[0].x, xy.sources[0].y), (10.25, 20.5));
        assert_eq!(xy.sources[0].flux, Some(1234.5));
        assert_eq!(xy.sources[0].sigma_px, Some(0.05));
        assert_eq!(xy.sources[0].p_star, Some(0.2));
        assert!(!xy.sources[0].truncated);
        assert_eq!((xy.sources[1].x, xy.sources[1].y), (0.0, 767.0));
        assert_eq!((xy.sources[1].flux, xy.sources[1].sigma_px), (None, None));
        assert!(xy.sources[1].truncated);
        write_sources(&path, &[], 10, 10).unwrap();
        assert!(read_sources(&path).unwrap().sources.is_empty());
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn sigma_and_p_star_columns_are_read() {
        let cards: &[&str] = &[
            "XTENSION= 'BINTABLE'",
            "BITPIX  =                    8",
            "NAXIS   =                    2",
            "NAXIS1  =                   16",
            "NAXIS2  =                    2",
            "PCOUNT  =                    0",
            "GCOUNT  =                    1",
            "TFIELDS =                    4",
            "TTYPE1  = 'X       '",
            "TFORM1  = '1E      '",
            "TTYPE2  = 'Y       '",
            "TFORM2  = '1E      '",
            "TTYPE3  = 'SIGMA   '",
            "TFORM3  = '1E      '",
            "TTYPE4  = 'P_STAR  '",
            "TFORM4  = '1E      '",
        ];
        let rows: Vec<u8> = [[10.0f32, 20.0, 0.25, 0.9], [30.0, 40.0, -1.0, 1.5]]
            .iter()
            .flat_map(|r| r.iter().flat_map(|v| v.to_be_bytes()))
            .collect();
        let xy = read(fits(&[(cards, rows)])).unwrap();
        assert_eq!(xy.sources[0].sigma_px, Some(0.25));
        assert!((xy.sources[0].p_star.unwrap() - 0.9).abs() < 1e-6);
        // Out-of-range values are ignored, not trusted.
        assert_eq!(xy.sources[1].sigma_px, None);
        assert_eq!(xy.sources[1].p_star, None);
    }

    #[test]
    fn the_first_table_with_x_and_y_is_used() {
        let meta: &[&str] = &[
            "XTENSION= 'BINTABLE'",
            "BITPIX  =                    8",
            "NAXIS   =                    2",
            "NAXIS1  =                    4",
            "NAXIS2  =                    1",
            "PCOUNT  =                    0",
            "GCOUNT  =                    1",
            "TFIELDS =                    1",
            "TTYPE1  = 'EXPTIME '",
            "TFORM1  = '1E      '",
        ];
        let xy = read(fits(&[
            (meta, 3.0f32.to_be_bytes().to_vec()),
            (XY_J, row(100, 20)),
        ]))
        .unwrap();
        assert_eq!(xy.sources.len(), 1);
    }

    #[test]
    fn absurd_dimensions_are_an_error_not_a_panic() {
        let huge: &[&str] = &[
            "XTENSION= 'BINTABLE'",
            "BITPIX  =                    8",
            "NAXIS   =                    2",
            "NAXIS1  =  9223372036854775807",
            "NAXIS2  =                    2",
            "PCOUNT  =                    0",
            "GCOUNT  =                    1",
            "TFIELDS =                    0",
        ];
        assert!(FitsFile::from_bytes(fits(&[(huge, Vec::new())])).is_err());
    }
}
