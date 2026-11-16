//! UCAC4 (Zacharias et al. 2013, AJ 145, 44) read in place from the original distribution.
//!
//! Layout expected (as distributed by USNO/CDS):
//!
//! ```text
//! <root>/u4b/z001 … z900    78-byte little-endian records, one file per 0.2° Dec zone, by RA
//! <root>/u4i/u4index.unf    optional; start/count per (zone, 0.25° RA bin) for direct seeks
//! <root>/u4i/u4hpm.dat      optional; true proper motions of high-proper-motion stars
//! ```
//!
//! `open` also accepts the `u4b` directory itself, in which case `u4i` is looked for next to it.
//! Field definitions follow `readme_u4.txt` (notes 6, 8 and 9 for the error and proper-motion
//! encodings).

use fayde_core::Portable;
use std::{
    collections::HashMap,
    fs::File,
    io::{Read, Seek, SeekFrom},
    path::{Path, PathBuf},
};

use fayde_core::{Epoch, SkyPos};

use crate::{
    propagate::{grow_error, propagate},
    Band, Catalog, CatalogError, CatalogId, Star, StarQuery,
};

/// Bytes per record.
pub const RECORD_BYTES: usize = 78;
/// Number of declination zones.
pub const N_ZONES: u16 = 900;
const N_RA_BINS: usize = 1440;
const ZONE_HEIGHT_DEG: f64 = 0.2;
const RA_BIN_DEG: f64 = 0.25;
const MISSING_MAG: i16 = 20_000;
const HPM_SENTINEL: i16 = 32_767;
/// Largest proper motion representable inline (mas/yr).
const MAX_INLINE_PM_MAS_YR: f64 = 3_276.6;
/// Proper-motion error used when the catalogue has none: readme note 9 sets "no data" to 500
/// (0.1 mas/yr units) for output tables.
pub const DEFAULT_SIGMA_PM_MAS_YR: f64 = 50.0;

/// Proper-motion uncertainty, mas/yr, of a star flagged as high proper motion whose supplement
/// entry is missing: its motion is only known to exceed what the 2-byte field holds (32767 ×
/// 0.1 mas/yr), so its position is left at J2000 with an error that grows at that rate.
pub const HPM_MISSING_SIGMA_PM_MAS_YR: f64 = 3276.7;

/// One decoded UCAC4 record. Field names follow `readme_u4.txt`; values are raw catalogue units.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(missing_docs)]
pub struct Ucac4Record {
    /// RA at epoch J2000, mas.
    pub ra: u32,
    /// South-pole distance at J2000, mas (Dec = spd − 90°).
    pub spd: u32,
    /// Model magnitude, mmag (20000 = missing).
    pub magm: i16,
    /// Aperture magnitude, mmag (20000 = missing).
    pub maga: i16,
    pub sigmag: u8,
    pub objt: u8,
    pub cdf: u8,
    /// Position error at central epoch, stored −128 (see [`Ucac4Record::sigma_pos_mas`]).
    pub sigra: i8,
    pub sigdc: i8,
    pub na1: u8,
    pub nu1: u8,
    pub cu1: u8,
    /// Central epoch − 1900, 0.01 yr.
    pub cepra: u16,
    pub cepdc: u16,
    /// Proper motion, 0.1 mas/yr (32767 = see u4hpm.dat).
    pub pmrac: i16,
    pub pmdc: i16,
    /// Proper-motion error code (see [`decode_sigma_pm`]).
    pub sigpmr: i8,
    pub sigpmd: i8,
    pub pts_key: u32,
    /// 2MASS J, H, Ks, mmag.
    pub jhk: [i16; 3],
    pub icqflg: [u8; 3],
    pub q2mflg: [u8; 3],
    /// APASS B, V, g, r, i, mmag (20000 = missing).
    pub apasm: [i16; 5],
    pub apase: [u8; 5],
    pub gcflg: u8,
    pub mcf: i32,
    pub leda: u8,
    pub x2m: u8,
    /// Unique star identifier.
    pub rnm: u32,
    pub zn2: u16,
    pub rn2: u32,
}

struct Cursor<'a> {
    b: &'a [u8; RECORD_BYTES],
    at: usize,
}

impl Cursor<'_> {
    fn take<const N: usize>(&mut self) -> [u8; N] {
        let mut out = [0u8; N];
        out.copy_from_slice(&self.b[self.at..self.at + N]);
        self.at += N;
        out
    }
    fn u8(&mut self) -> u8 {
        self.take::<1>()[0]
    }
    fn i8(&mut self) -> i8 {
        i8::from_le_bytes(self.take())
    }
    fn u16(&mut self) -> u16 {
        u16::from_le_bytes(self.take())
    }
    fn i16(&mut self) -> i16 {
        i16::from_le_bytes(self.take())
    }
    fn u32(&mut self) -> u32 {
        u32::from_le_bytes(self.take())
    }
    fn i32(&mut self) -> i32 {
        i32::from_le_bytes(self.take())
    }
}

impl Ucac4Record {
    /// Decode one 78-byte record.
    pub fn decode(bytes: &[u8; RECORD_BYTES]) -> Self {
        let mut c = Cursor { b: bytes, at: 0 };
        let r = Self {
            ra: c.u32(),
            spd: c.u32(),
            magm: c.i16(),
            maga: c.i16(),
            sigmag: c.u8(),
            objt: c.u8(),
            cdf: c.u8(),
            sigra: c.i8(),
            sigdc: c.i8(),
            na1: c.u8(),
            nu1: c.u8(),
            cu1: c.u8(),
            cepra: c.u16(),
            cepdc: c.u16(),
            pmrac: c.i16(),
            pmdc: c.i16(),
            sigpmr: c.i8(),
            sigpmd: c.i8(),
            pts_key: c.u32(),
            jhk: [c.i16(), c.i16(), c.i16()],
            icqflg: [c.u8(), c.u8(), c.u8()],
            q2mflg: [c.u8(), c.u8(), c.u8()],
            apasm: [c.i16(), c.i16(), c.i16(), c.i16(), c.i16()],
            apase: [c.u8(), c.u8(), c.u8(), c.u8(), c.u8()],
            gcflg: c.u8(),
            mcf: c.i32(),
            leda: c.u8(),
            x2m: c.u8(),
            rnm: c.u32(),
            zn2: c.u16(),
            rn2: c.u32(),
        };
        debug_assert_eq!(c.at, RECORD_BYTES);
        r
    }

    /// Position at epoch J2000.0 (ICRS).
    pub fn position_j2000(&self) -> SkyPos {
        let ra = (f64::from(self.ra) / 3.6e6).to_radians();
        let dec = (f64::from(self.spd) / 3.6e6 - 90.0).to_radians();
        // spd is at most 180° by format; clamp guards a corrupt record.
        SkyPos::new(
            ra,
            dec.clamp(-std::f64::consts::FRAC_PI_2, std::f64::consts::FRAC_PI_2),
        )
        .expect("finite by construction")
    }

    /// Position errors at the central epochs, mas, `[RA·cos δ, Dec]` (readme note 6).
    pub fn sigma_pos_mas(&self) -> [f64; 2] {
        [
            f64::from(i16::from(self.sigra) + 128),
            f64::from(i16::from(self.sigdc) + 128),
        ]
    }

    /// Central epochs of the mean RA and Dec.
    pub fn central_epochs(&self) -> [Epoch; 2] {
        let e = |v: u16| Epoch::julian_year(1900.0 + 0.01 * f64::from(v)).expect("finite");
        [e(self.cepra), e(self.cepdc)]
    }

    /// Whether the true proper motion is in the high-proper-motion supplement.
    pub fn is_high_pm(&self) -> bool {
        self.pmrac == HPM_SENTINEL || self.pmdc == HPM_SENTINEL
    }

    /// Proper motion `[μα*, μδ]` in mas/yr. A component is `None` when the catalogue flags it as
    /// having no data (error code 255), or when the star is in the high-proper-motion supplement
    /// and `hpm` has no entry for it.
    pub fn pm_mas_yr(&self, hpm: Option<[f64; 2]>) -> [Option<f64>; 2] {
        if self.is_high_pm() {
            return match hpm {
                Some(pm) => [Some(pm[0]), Some(pm[1])],
                None => [None, None],
            };
        }
        let comp = |pm: i16, sig: i8| decode_sigma_pm(sig).map(|_| 0.1 * f64::from(pm));
        [comp(self.pmrac, self.sigpmr), comp(self.pmdc, self.sigpmd)]
    }

    /// Proper-motion errors `[σμα*, σμδ]`, mas/yr, defaulting to
    /// [`DEFAULT_SIGMA_PM_MAS_YR`] where the catalogue has none.
    pub fn sigma_pm_mas_yr(&self) -> [f64; 2] {
        [
            decode_sigma_pm(self.sigpmr).unwrap_or(DEFAULT_SIGMA_PM_MAS_YR),
            decode_sigma_pm(self.sigpmd).unwrap_or(DEFAULT_SIGMA_PM_MAS_YR),
        ]
    }

    /// Magnitude in `band`, or `None` if the catalogue has no value.
    pub fn mag(&self, band: Band) -> Option<f32> {
        let m = |raw: i16| (raw != MISSING_MAG).then(|| f32::from(raw) / 1000.0);
        match band {
            // Model magnitude unless it is the 20.000 placeholder, then aperture.
            Band::Default => {
                if self.magm.unsigned_abs() < 19_999 {
                    m(self.magm)
                } else {
                    m(self.maga)
                }
            }
            Band::UcacModel => m(self.magm),
            Band::UcacAperture => m(self.maga),
            Band::ApassB => m(self.apasm[0]),
            Band::ApassV => m(self.apasm[1]),
            Band::ApassG => m(self.apasm[2]),
            Band::ApassR => m(self.apasm[3]),
            Band::ApassI => m(self.apasm[4]),
            // 2MASS uses 0 / 20000-style placeholders inconsistently; treat both as missing.
            Band::TwoMassJ => m(self.jhk[0]).filter(|v| *v > 0.0),
            Band::TwoMassH => m(self.jhk[1]).filter(|v| *v > 0.0),
            Band::TwoMassK => m(self.jhk[2]).filter(|v| *v > 0.0),
        }
    }

    /// Build a [`Star`] at `epoch` (or J2000 if `None`).
    pub fn to_star(&self, hpm: Option<[f64; 2]>, epoch: Option<Epoch>, band: Band) -> Star {
        let pm = self.pm_mas_yr(hpm);
        // Without its supplement entry, a high-PM star's catalogue errors describe the sentinel,
        // not its motion: the error must cover the (unknown, large) motion instead.
        let sig_pm = if self.is_high_pm() && hpm.is_none() {
            [HPM_MISSING_SIGMA_PM_MAS_YR; 2]
        } else {
            self.sigma_pm_mas_yr()
        };
        let target = epoch.unwrap_or(Epoch::J2000);
        let pos = propagate(
            self.position_j2000(),
            [pm[0].unwrap_or(0.0), pm[1].unwrap_or(0.0)],
            Epoch::J2000,
            target,
        );
        let sig = self.sigma_pos_mas();
        let cep = self.central_epochs();
        Star {
            id: u64::from(self.rnm),
            pos,
            mag: self.mag(band),
            pos_err_mas: [
                grow_error(sig[0], cep[0], sig_pm[0], target) as f32,
                grow_error(sig[1], cep[1], sig_pm[1], target) as f32,
            ],
            pm_mas_yr: [pm[0].map(|v| v as f32), pm[1].map(|v| v as f32)],
        }
    }
}

/// Decode a proper-motion error code (readme note 9) to mas/yr. `None` means "no data".
pub fn decode_sigma_pm(raw: i8) -> Option<f64> {
    match i16::from(raw) + 128 {
        v @ 0..=250 => Some(0.1 * f64::from(v)),
        251 => Some(27.5),
        252 => Some(32.5),
        253 => Some(37.5),
        254 => Some(45.0),
        _ => None,
    }
}

/// 1-based zone containing a J2000 declination (degrees).
pub fn zone_number(dec_deg: f64) -> u16 {
    let z = ((dec_deg + 90.0) / ZONE_HEIGHT_DEG).floor() as i64 + 1;
    z.clamp(1, i64::from(N_ZONES)) as u16
}

/// Start/count table from `u4index.unf`.
struct ZoneIndex {
    n0: Vec<i32>,
    nn: Vec<i32>,
}

impl ZoneIndex {
    fn load(path: &Path) -> Result<Self, CatalogError> {
        let bytes = std::fs::read(path).map_err(|source| CatalogError::Io {
            path: path.into(),
            source,
        })?;
        let n = usize::from(N_ZONES) * N_RA_BINS;
        if bytes.len() != 2 * n * 4 {
            return Err(CatalogError::Malformed {
                path: path.into(),
                reason: format!("expected {} bytes, found {}", 2 * n * 4, bytes.len()),
            });
        }
        let ints: Vec<i32> = bytes
            .chunks_exact(4)
            .map(|c| i32::from_le_bytes([c[0], c[1], c[2], c[3]]))
            .collect();
        Ok(Self {
            n0: ints[..n].to_vec(),
            nn: ints[n..].to_vec(),
        })
    }

    /// Record range `[start, end)` (0-based) covering RA bins `lo..=hi` of `zone`.
    fn range(&self, zone: u16, lo: usize, hi: usize) -> (u64, u64) {
        let i = |bin: usize| usize::from(zone - 1) + usize::from(N_ZONES) * bin;
        let start = self.n0[i(lo)].max(0) as u64;
        // In i64: the two counts of a damaged index can sum past i32.
        let end = (i64::from(self.n0[i(hi)].max(0)) + i64::from(self.nn[i(hi)].max(0))) as u64;
        (start, end.max(start))
    }
}

/// A star found by a cone search, with its location in the catalogue.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ucac4Hit {
    /// Zone, 1–900.
    pub zone: u16,
    /// Record number within the zone, 1-based (as in UCAC4 designations `UCAC4-zzz-nnnnnn`).
    pub record: u32,
    /// Raw record.
    pub raw: Ucac4Record,
    /// Star at the query epoch.
    pub star: Star,
}

/// UCAC4 read in place.
pub struct Ucac4Catalog {
    zone_dir: PathBuf,
    index: Option<ZoneIndex>,
    hpm: HashMap<u32, [f64; 2]>,
    max_pm_mas_yr: f64,
}

impl Ucac4Catalog {
    /// Open a UCAC4 distribution: either its root (holding `u4b/` and optionally `u4i/`) or the
    /// `u4b/` zone directory itself.
    pub fn open(path: impl AsRef<Path>) -> Result<Self, CatalogError> {
        let path = path.as_ref();
        let (zone_dir, aux_dir) = if path.join("u4b").is_dir() {
            (path.join("u4b"), path.join("u4i"))
        } else if path.is_dir() {
            let parent = path.parent().unwrap_or(path);
            (path.to_path_buf(), parent.join("u4i"))
        } else {
            return Err(CatalogError::NotFound(format!(
                "{} is not a UCAC4 directory",
                path.display()
            )));
        };
        let has_any_zone = (1..=N_ZONES).any(|z| zone_dir.join(zone_file(z)).is_file());
        if !has_any_zone {
            return Err(CatalogError::NotFound(format!(
                "no zone files z001…z900 in {}",
                zone_dir.display()
            )));
        }
        let index_path = aux_dir.join("u4index.unf");
        let index = if index_path.is_file() {
            Some(ZoneIndex::load(&index_path)?)
        } else {
            None
        };
        let hpm_path = aux_dir.join("u4hpm.dat");
        let hpm = if hpm_path.is_file() {
            load_hpm(&hpm_path)?
        } else {
            HashMap::new()
        };
        let max_pm_mas_yr = hpm
            .values()
            .map(|p| p[0].phypot(p[1]))
            .fold(MAX_INLINE_PM_MAS_YR * std::f64::consts::SQRT_2, f64::max);
        Ok(Self {
            zone_dir,
            index,
            hpm,
            max_pm_mas_yr,
        })
    }

    /// Whether `u4index.unf` was found (cone searches seek instead of reading whole zones).
    pub fn has_index(&self) -> bool {
        self.index.is_some()
    }

    /// Number of high-proper-motion supplement entries loaded.
    pub fn high_pm_count(&self) -> usize {
        self.hpm.len()
    }

    /// Path of a zone file.
    pub fn zone_path(&self, zone: u16) -> PathBuf {
        self.zone_dir.join(zone_file(zone))
    }

    /// Read records `[start, end)` (0-based) of a zone. `Ok(None)` if the zone file is absent.
    fn read_records(
        &self,
        zone: u16,
        range: Option<(u64, u64)>,
    ) -> Result<Option<Vec<Ucac4Record>>, CatalogError> {
        let path = self.zone_path(zone);
        let mut f = match File::open(&path) {
            Ok(f) => f,
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => return Ok(None),
            Err(source) => return Err(CatalogError::Io { path, source }),
        };
        let io = |source| CatalogError::Io {
            path: path.clone(),
            source,
        };
        let mut buf = Vec::new();
        match range {
            None => {
                f.read_to_end(&mut buf).map_err(io)?;
                if buf.len() % RECORD_BYTES != 0 {
                    return Err(CatalogError::Malformed {
                        path: path.clone(),
                        reason: format!("size {} is not a multiple of {RECORD_BYTES}", buf.len()),
                    });
                }
            }
            Some((start, end)) => {
                if end == start {
                    return Ok(Some(Vec::new()));
                }
                // Check the range against the file before sizing the buffer from it: the index
                // is a separate file and may be damaged.
                let available = f.metadata().map_err(io)?.len() / RECORD_BYTES as u64;
                if end > available {
                    return Err(CatalogError::Malformed {
                        path: path.clone(),
                        reason: format!(
                            "index asks for records {start}..{end}, zone has {available}; index and zone disagree"
                        ),
                    });
                }
                f.seek(SeekFrom::Start(start * RECORD_BYTES as u64))
                    .map_err(io)?;
                buf.resize(((end - start) as usize) * RECORD_BYTES, 0);
                f.read_exact(&mut buf)
                    .map_err(|e| CatalogError::Malformed {
                        path: path.clone(),
                        reason: format!(
                            "records {start}..{end} not readable ({e}); index and zone disagree"
                        ),
                    })?;
            }
        }
        Ok(Some(
            buf.chunks_exact(RECORD_BYTES)
                .map(|c| Ucac4Record::decode(c.try_into().expect("exact chunk")))
                .collect(),
        ))
    }

    /// Every record of a zone, in catalogue order. `Ok(None)` if the zone file is absent.
    pub fn read_zone(&self, zone: u16) -> Result<Option<Vec<Ucac4Record>>, CatalogError> {
        self.read_records(zone, None)
    }

    /// Cone search returning catalogue locations and raw records as well as stars.
    pub fn cone_hits(
        &self,
        centre: SkyPos,
        radius: f64,
        query: &StarQuery,
    ) -> Result<Vec<Ucac4Hit>, CatalogError> {
        self.search(centre, radius, query, true)
    }

    /// Same as [`Ucac4Catalog::cone_hits`] but reads whole zones instead of seeking with
    /// `u4index.unf`. Slower; used to check the index path.
    pub fn cone_hits_full_scan(
        &self,
        centre: SkyPos,
        radius: f64,
        query: &StarQuery,
    ) -> Result<Vec<Ucac4Hit>, CatalogError> {
        self.search(centre, radius, query, false)
    }

    fn search(
        &self,
        centre: SkyPos,
        radius: f64,
        query: &StarQuery,
        use_index: bool,
    ) -> Result<Vec<Ucac4Hit>, CatalogError> {
        if !radius.is_finite() || !(0.0..std::f64::consts::PI).contains(&radius) {
            return Err(CatalogError::InvalidQuery(format!("radius {radius} rad")));
        }
        // Zones and RA bins are by J2000 position: widen the search by the largest motion any
        // star can have made since J2000.
        let pad = query
            .epoch
            .map(|e| e.years_since(Epoch::J2000).abs() * self.max_pm_mas_yr * crate::propagate::MAS)
            .unwrap_or(0.0);
        let r = radius + pad + 1e-9;
        let (ra0, dec0) = (centre.ra.to_degrees(), centre.dec.to_degrees());
        let r_deg = r.to_degrees();
        let z_lo = zone_number(dec0 - r_deg);
        let z_hi = zone_number(dec0 + r_deg);
        let bins = ra_bins(ra0, dec0, r);

        let mut out = Vec::new();
        for zone in z_lo..=z_hi {
            let chunks: Vec<(Option<(u64, u64)>, u64)> = match (&self.index, use_index) {
                (Some(idx), true) => bins
                    .iter()
                    .map(|&(lo, hi)| {
                        let range = idx.range(zone, lo, hi);
                        (Some(range), range.0)
                    })
                    .collect(),
                _ => vec![(None, 0)],
            };
            for (range, first) in chunks {
                let Some(records) = self.read_records(zone, range)? else {
                    if query.strict {
                        return Err(CatalogError::NotFound(format!(
                            "zone file {} (needed for this cone)",
                            self.zone_path(zone).display()
                        )));
                    }
                    break;
                };
                for (k, raw) in records.iter().enumerate() {
                    let hpm = raw
                        .is_high_pm()
                        .then(|| self.hpm.get(&raw.rnm).copied())
                        .flatten();
                    let star = raw.to_star(hpm, query.epoch, query.band);
                    if centre.separation(star.pos) > radius {
                        continue;
                    }
                    if let Some(limit) = query.mag_max {
                        if !star.mag.is_some_and(|m| m <= limit) {
                            continue;
                        }
                    }
                    if query.require_proper_motion
                        && !(star.pm_mas_yr[0].is_some() && star.pm_mas_yr[1].is_some())
                    {
                        continue;
                    }
                    out.push(Ucac4Hit {
                        zone,
                        record: (first + k as u64 + 1) as u32,
                        raw: *raw,
                        star,
                    });
                }
            }
        }
        out.sort_by_key(|h| (h.zone, h.record));
        out.dedup_by_key(|h| (h.zone, h.record));
        Ok(out)
    }
}

impl Catalog for Ucac4Catalog {
    fn id(&self) -> CatalogId {
        CatalogId {
            name: "UCAC4",
            release: "2012",
            location: self.zone_dir.display().to_string(),
        }
    }

    fn cone(
        &self,
        centre: SkyPos,
        radius: f64,
        query: &StarQuery,
    ) -> Result<Vec<Star>, CatalogError> {
        Ok(self
            .cone_hits(centre, radius, query)?
            .into_iter()
            .map(|h| h.star)
            .collect())
    }

    fn scan(&self, query: &StarQuery, visit: &mut dyn FnMut(Star)) -> Result<(), CatalogError> {
        for zone in 1..=N_ZONES {
            let Some(records) = self.read_zone(zone)? else {
                if query.strict {
                    return Err(CatalogError::NotFound(format!(
                        "zone file {}",
                        self.zone_path(zone).display()
                    )));
                }
                continue;
            };
            for raw in &records {
                if let Some(limit) = query.mag_max {
                    // Cheap pre-check on the raw magnitude before propagating.
                    if !raw.mag(query.band).is_some_and(|m| m <= limit) {
                        continue;
                    }
                }
                let hpm = raw
                    .is_high_pm()
                    .then(|| self.hpm.get(&raw.rnm).copied())
                    .flatten();
                let star = raw.to_star(hpm, query.epoch, query.band);
                if query.require_proper_motion
                    && !(star.pm_mas_yr[0].is_some() && star.pm_mas_yr[1].is_some())
                {
                    continue;
                }
                visit(star);
            }
        }
        Ok(())
    }
}

fn zone_file(zone: u16) -> String {
    format!("z{zone:03}")
}

/// Inclusive RA-bin ranges covering a cone of radius `r` (rad) at (`ra0`, `dec0`) in degrees.
fn ra_bins(ra0: f64, dec0: f64, r: f64) -> Vec<(usize, usize)> {
    let all = vec![(0, N_RA_BINS - 1)];
    let dec0r = dec0.to_radians();
    if dec0r.abs() + r >= std::f64::consts::FRAC_PI_2 {
        return all;
    }
    let s = r.psin() / dec0r.pcos();
    if s >= 1.0 {
        return all;
    }
    let half = s.pasin().to_degrees();
    let bin = |ra: f64| ((ra.rem_euclid(360.0) / RA_BIN_DEG).floor() as usize).min(N_RA_BINS - 1);
    let (lo, hi) = (ra0 - half, ra0 + half);
    if hi - lo >= 360.0 {
        return all;
    }
    let (blo, bhi) = (bin(lo), bin(hi));
    if blo <= bhi && lo >= 0.0 && hi < 360.0 {
        vec![(blo, bhi)]
    } else {
        vec![(blo, N_RA_BINS - 1), (0, bhi)]
    }
}

fn load_hpm(path: &Path) -> Result<HashMap<u32, [f64; 2]>, CatalogError> {
    let text = std::fs::read_to_string(path).map_err(|source| CatalogError::Io {
        path: path.into(),
        source,
    })?;
    let mut out = HashMap::new();
    for (lineno, line) in text.lines().enumerate() {
        let cols: Vec<&str> = line.split_whitespace().collect();
        if cols.is_empty() {
            continue;
        }
        let bad = || CatalogError::Malformed {
            path: path.into(),
            reason: format!("line {}: expected rnm zn rnz pmrc pmd …", lineno + 1),
        };
        if cols.len() < 5 {
            return Err(bad());
        }
        let rnm: u32 = cols[0].parse().map_err(|_| bad())?;
        let pmrc: f64 = cols[3].parse().map_err(|_| bad())?;
        let pmd: f64 = cols[4].parse().map_err(|_| bad())?;
        if !(pmrc.is_finite() && pmd.is_finite()) {
            return Err(bad());
        }
        out.insert(rnm, [0.1 * pmrc, 0.1 * pmd]);
    }
    Ok(out)
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    /// Encode a record (inverse of `decode`) for tests.
    pub(crate) fn encode(r: &Ucac4Record) -> [u8; RECORD_BYTES] {
        let mut v = Vec::with_capacity(RECORD_BYTES);
        v.extend(r.ra.to_le_bytes());
        v.extend(r.spd.to_le_bytes());
        v.extend(r.magm.to_le_bytes());
        v.extend(r.maga.to_le_bytes());
        v.extend([r.sigmag, r.objt, r.cdf]);
        v.extend(r.sigra.to_le_bytes());
        v.extend(r.sigdc.to_le_bytes());
        v.extend([r.na1, r.nu1, r.cu1]);
        v.extend(r.cepra.to_le_bytes());
        v.extend(r.cepdc.to_le_bytes());
        v.extend(r.pmrac.to_le_bytes());
        v.extend(r.pmdc.to_le_bytes());
        v.extend(r.sigpmr.to_le_bytes());
        v.extend(r.sigpmd.to_le_bytes());
        v.extend(r.pts_key.to_le_bytes());
        for m in r.jhk {
            v.extend(m.to_le_bytes());
        }
        v.extend(r.icqflg);
        v.extend(r.q2mflg);
        for m in r.apasm {
            v.extend(m.to_le_bytes());
        }
        v.extend(r.apase);
        v.push(r.gcflg);
        v.extend(r.mcf.to_le_bytes());
        v.extend([r.leda, r.x2m]);
        v.extend(r.rnm.to_le_bytes());
        v.extend(r.zn2.to_le_bytes());
        v.extend(r.rn2.to_le_bytes());
        v.try_into().expect("78 bytes")
    }

    pub(crate) fn sample() -> Ucac4Record {
        Ucac4Record {
            ra: 540_000_000,  // 150°
            spd: 331_200_000, // 92° → Dec +2°
            magm: 12_345,
            maga: 12_400,
            sigmag: 3,
            objt: 0,
            cdf: 0,
            sigra: -113, // 15 mas
            sigdc: -108, // 20 mas
            na1: 10,
            nu1: 9,
            cu1: 4,
            cepra: 9_850, // 1998.50
            cepdc: 9_900, // 1999.00
            pmrac: 1_234, // 123.4 mas/yr
            pmdc: -567,   // -56.7 mas/yr
            sigpmr: -108, // 20 → 2.0 mas/yr
            sigpmd: 125,  // 253 → 37.5 mas/yr
            pts_key: 7,
            jhk: [11_000, 10_800, 10_700],
            icqflg: [5, 5, 5],
            q2mflg: [0, 0, 0],
            apasm: [13_000, 12_500, MISSING_MAG, 12_200, 12_000],
            apase: [10, 10, 0, 10, 10],
            gcflg: 0,
            mcf: 100_000_000,
            leda: 0,
            x2m: 0,
            rnm: 123_456_789,
            zn2: 0,
            rn2: 0,
        }
    }

    #[test]
    fn damaged_index_is_an_error_not_an_allocation() {
        // One zone of three records, and an index whose every count is i32::MAX: the range it
        // gives must be checked against the zone file before a buffer is sized from it.
        let dir = std::env::temp_dir().join(format!("fayde-bad-u4index-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(dir.join("u4b")).unwrap();
        std::fs::create_dir_all(dir.join("u4i")).unwrap();
        let rec = encode(&sample());
        std::fs::write(dir.join("u4b").join("z461"), [rec, rec, rec].concat()).unwrap();
        let n = usize::from(N_ZONES) * N_RA_BINS;
        let mut index = vec![0u8; n * 4];
        index.extend(std::iter::repeat_n(i32::MAX.to_le_bytes(), n).flatten());
        std::fs::write(dir.join("u4i").join("u4index.unf"), index).unwrap();
        let cat = Ucac4Catalog::open(&dir).unwrap();
        // Mid-zone 461 (Dec +2.0° to +2.2°), so the cone needs no other zone file.
        let centre = SkyPos::from_degrees(150.0, 2.1).unwrap();
        let err = cat
            .cone_hits(centre, 0.01f64.to_radians(), &StarQuery::default())
            .unwrap_err();
        assert!(matches!(err, CatalogError::Malformed { .. }), "{err}");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn decode_round_trip_and_fields() {
        let r = sample();
        assert_eq!(Ucac4Record::decode(&encode(&r)), r);
        let (ra, dec) = r.position_j2000().to_degrees();
        assert!((ra - 150.0).abs() < 1e-12 && (dec - 2.0).abs() < 1e-12);
        assert_eq!(r.sigma_pos_mas(), [15.0, 20.0]);
        let [e1, e2] = r.central_epochs();
        assert!((e1.year() - 1998.5).abs() < 1e-9 && (e2.year() - 1999.0).abs() < 1e-9);
        assert_eq!(
            r.pm_mas_yr(None)[0].map(|v| (v * 10.0).round()),
            Some(1234.0)
        );
        assert_eq!(r.sigma_pm_mas_yr(), [2.0, 37.5]);
    }

    #[test]
    fn high_pm_star_without_supplement_gets_an_error_covering_its_motion() {
        let mut r = sample();
        r.pmrac = HPM_SENTINEL;
        let at = |hpm| {
            r.to_star(
                hpm,
                Some(Epoch::julian_year(2025.0).unwrap()),
                Band::Default,
            )
            .pos_err_mas
        };
        // With the supplement entry: catalogue errors, small.
        assert!(at(Some([4000.0, -100.0]))[0] < 200.0);
        // Without it: the position stays at J2000 and the error grows at ≥ 3.3″/yr.
        let e = at(None);
        assert!(e[0] > 25.0 * 3000.0 && e[1] > 25.0 * 3000.0, "{e:?}");
    }

    #[test]
    fn magnitudes_and_fallbacks() {
        let mut r = sample();
        assert_eq!(r.mag(Band::Default), Some(12.345));
        assert_eq!(r.mag(Band::ApassG), None);
        assert_eq!(r.mag(Band::ApassR), Some(12.2));
        r.magm = MISSING_MAG;
        assert_eq!(r.mag(Band::Default), Some(12.4));
        r.maga = MISSING_MAG;
        assert_eq!(r.mag(Band::Default), None);
    }

    #[test]
    fn sigma_pm_codes() {
        let code = |v: i16| (v - 128) as i8;
        assert_eq!(decode_sigma_pm(code(0)), Some(0.0));
        assert_eq!(decode_sigma_pm(code(250)), Some(25.0));
        assert_eq!(decode_sigma_pm(code(251)), Some(27.5));
        assert_eq!(decode_sigma_pm(code(254)), Some(45.0));
        assert_eq!(decode_sigma_pm(code(255)), None);
    }

    #[test]
    fn missing_pm_component_and_high_pm() {
        let mut r = sample();
        r.sigpmd = 127; // 255: no data
        assert!(r.pm_mas_yr(None)[1].is_none());
        assert_eq!(r.sigma_pm_mas_yr()[1], DEFAULT_SIGMA_PM_MAS_YR);
        r.pmrac = HPM_SENTINEL;
        assert_eq!(r.pm_mas_yr(None), [None, None]);
        assert_eq!(
            r.pm_mas_yr(Some([4108.7, 3141.3])),
            [Some(4108.7), Some(3141.3)]
        );
    }

    #[test]
    fn star_at_epoch_moves_and_error_grows() {
        let r = sample();
        let s2000 = r.to_star(None, None, Band::Default);
        let s2026 = r.to_star(None, Epoch::julian_year(2026.0), Band::Default);
        let moved = fayde_core::to_arcsec(s2000.pos.separation(s2026.pos)) * 1000.0;
        let expect = 26.0 * (123.4f64.phypot(56.7));
        assert!((moved - expect).abs() < 0.01, "{moved} vs {expect} mas");
        assert!(
            s2026.pos_err_mas[1] > 900.0,
            "37.5 mas/yr × 27 yr dominates"
        );
    }

    #[test]
    fn zone_numbers_and_ra_bins() {
        assert_eq!(zone_number(-90.0), 1);
        assert_eq!(zone_number(90.0), 900);
        assert_eq!(zone_number(2.0), 461);
        assert_eq!(ra_bins(10.0, 0.0, 0.01), vec![(37, 42)]);
        // 0.01 rad = 0.573°: 0.1 − 0.573 = 359.527° → bin 1438.
        assert_eq!(ra_bins(0.1, 0.0, 0.01), vec![(1438, 1439), (0, 2)]);
        assert_eq!(ra_bins(10.0, 89.9, 0.01), vec![(0, 1439)]);
    }
}
