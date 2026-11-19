use fayde_core::Portable;
use std::{
    collections::HashMap,
    fs::File,
    io::{Read, Seek, SeekFrom},
    path::{Path, PathBuf},
};

use fayde_core::{Epoch, SkyPos};

use crate::{
    propagate::{grow_error, propagate, propagate_in, LocalFrame},
    Band, Catalog, CatalogError, CatalogId, Star, StarQuery,
};

pub const RECORD_BYTES: usize = 78;
pub const N_ZONES: u16 = 900;
const N_RA_BINS: usize = 1440;
const ZONE_HEIGHT_DEG: f64 = 0.2;
const RA_BIN_DEG: f64 = 0.25;
const RA_BIN_MAS: u64 = 900_000;
const MISSING_MAG: i16 = 20_000;
const HPM_SENTINEL: i16 = 32_767;
const MAX_INLINE_PM_MAS_YR: f64 = 3_276.6;
pub const DEFAULT_SIGMA_PM_MAS_YR: f64 = 50.0;

pub const HPM_MISSING_SIGMA_PM_MAS_YR: f64 = 3276.7;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Ucac4Record {
    pub ra: u32,
    pub spd: u32,
    pub magm: i16,
    pub maga: i16,
    pub sigmag: u8,
    pub objt: u8,
    pub cdf: u8,
    pub sigra: i8,
    pub sigdc: i8,
    pub na1: u8,
    pub nu1: u8,
    pub cu1: u8,
    pub cepra: u16,
    pub cepdc: u16,
    pub pmrac: i16,
    pub pmdc: i16,
    pub sigpmr: i8,
    pub sigpmd: i8,
    pub pts_key: u32,
    pub jhk: [i16; 3],
    pub icqflg: [u8; 3],
    pub q2mflg: [u8; 3],
    pub apasm: [i16; 5],
    pub apase: [u8; 5],
    pub gcflg: u8,
    pub mcf: i32,
    pub leda: u8,
    pub x2m: u8,
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

    pub fn position_j2000(&self) -> SkyPos {
        let ra = (f64::from(self.ra) / 3.6e6).to_radians();
        let dec = (f64::from(self.spd) / 3.6e6 - 90.0).to_radians();
        SkyPos::new(
            ra,
            dec.clamp(-std::f64::consts::FRAC_PI_2, std::f64::consts::FRAC_PI_2),
        )
        .expect("finite by construction")
    }

    pub fn sigma_pos_mas(&self) -> [f64; 2] {
        [
            f64::from(i16::from(self.sigra) + 128),
            f64::from(i16::from(self.sigdc) + 128),
        ]
    }

    pub fn central_epochs(&self) -> [Epoch; 2] {
        let e = |v: u16| Epoch::julian_year(1900.0 + 0.01 * f64::from(v)).expect("finite");
        [e(self.cepra), e(self.cepdc)]
    }

    pub fn is_high_pm(&self) -> bool {
        self.pmrac == HPM_SENTINEL || self.pmdc == HPM_SENTINEL
    }

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

    pub fn sigma_pm_mas_yr(&self) -> [f64; 2] {
        [
            decode_sigma_pm(self.sigpmr).unwrap_or(DEFAULT_SIGMA_PM_MAS_YR),
            decode_sigma_pm(self.sigpmd).unwrap_or(DEFAULT_SIGMA_PM_MAS_YR),
        ]
    }

    pub fn mag(&self, band: Band) -> Option<f32> {
        let m = |raw: i16| (raw != MISSING_MAG).then(|| f32::from(raw) / 1000.0);
        match band {
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
            Band::TwoMassJ => m(self.jhk[0]).filter(|v| *v > 0.0),
            Band::TwoMassH => m(self.jhk[1]).filter(|v| *v > 0.0),
            Band::TwoMassK => m(self.jhk[2]).filter(|v| *v > 0.0),
        }
    }

    pub fn to_star(&self, hpm: Option<[f64; 2]>, epoch: Option<Epoch>, band: Band) -> Star {
        let pm = self.pm_mas_yr(hpm);
        let target = epoch.unwrap_or(Epoch::J2000);
        let pos = propagate(
            self.position_j2000(),
            [pm[0].unwrap_or(0.0), pm[1].unwrap_or(0.0)],
            Epoch::J2000,
            target,
        );
        self.star_at(pos, pm, hpm, target, band)
    }

    fn star_at(
        &self,
        pos: SkyPos,
        pm: [Option<f64>; 2],
        hpm: Option<[f64; 2]>,
        target: Epoch,
        band: Band,
    ) -> Star {
        let sig_pm = if self.is_high_pm() && hpm.is_none() {
            [HPM_MISSING_SIGMA_PM_MAS_YR; 2]
        } else {
            self.sigma_pm_mas_yr()
        };
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

pub fn zone_number(dec_deg: f64) -> u16 {
    let z = ((dec_deg + 90.0) / ZONE_HEIGHT_DEG).floor() as i64 + 1;
    z.clamp(1, i64::from(N_ZONES)) as u16
}

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

    fn range(&self, zone: u16, lo: usize, hi: usize) -> (u64, u64) {
        let i = |bin: usize| usize::from(zone - 1) + usize::from(N_ZONES) * bin;
        let start = self.n0[i(lo)].max(0) as u64;
        let end = (i64::from(self.n0[i(hi)].max(0)) + i64::from(self.nn[i(hi)].max(0))) as u64;
        (start, end.max(start))
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ucac4Hit {
    pub zone: u16,
    pub record: u32,
    pub raw: Ucac4Record,
    pub star: Star,
}

pub struct Ucac4Catalog {
    zone_dir: PathBuf,
    index: Option<ZoneIndex>,
    hpm: HashMap<u32, [f64; 2]>,
    max_pm_mas_yr: f64,
}

impl Ucac4Catalog {
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

    pub fn has_index(&self) -> bool {
        self.index.is_some()
    }

    pub fn high_pm_count(&self) -> usize {
        self.hpm.len()
    }

    pub fn zone_path(&self, zone: u16) -> PathBuf {
        self.zone_dir.join(zone_file(zone))
    }

    fn open_zone(&self, zone: u16) -> Result<Option<(File, PathBuf)>, CatalogError> {
        let path = self.zone_path(zone);
        match File::open(&path) {
            Ok(f) => Ok(Some((f, path))),
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(None),
            Err(source) => Err(CatalogError::Io { path, source }),
        }
    }

    fn read_records(
        &self,
        zone: u16,
        range: Option<(u64, u64)>,
    ) -> Result<Option<Vec<Ucac4Record>>, CatalogError> {
        let Some((mut f, path)) = self.open_zone(zone)? else {
            return Ok(None);
        };
        Ok(Some(
            read_range(&mut f, &path, range)?
                .chunks_exact(RECORD_BYTES)
                .map(|c| Ucac4Record::decode(c.try_into().expect("exact chunk")))
                .collect(),
        ))
    }

    pub fn read_zone(&self, zone: u16) -> Result<Option<Vec<Ucac4Record>>, CatalogError> {
        self.read_records(zone, None)
    }

    fn lookup(&self) -> Lookup {
        if self.index.is_some() {
            Lookup::Index
        } else {
            Lookup::RaSearch
        }
    }

    pub fn cone_hits(
        &self,
        centre: SkyPos,
        radius: f64,
        query: &StarQuery,
    ) -> Result<Vec<Ucac4Hit>, CatalogError> {
        self.hits(centre, radius, query, self.lookup())
    }

    pub fn cone_hits_full_scan(
        &self,
        centre: SkyPos,
        radius: f64,
        query: &StarQuery,
    ) -> Result<Vec<Ucac4Hit>, CatalogError> {
        self.hits(centre, radius, query, Lookup::FullScan)
    }

    pub fn cone_hits_ra_search(
        &self,
        centre: SkyPos,
        radius: f64,
        query: &StarQuery,
    ) -> Result<Vec<Ucac4Hit>, CatalogError> {
        self.hits(centre, radius, query, Lookup::RaSearch)
    }

    fn hits(
        &self,
        centre: SkyPos,
        radius: f64,
        query: &StarQuery,
        lookup: Lookup,
    ) -> Result<Vec<Ucac4Hit>, CatalogError> {
        let test = ConeTest::new(self, centre, radius, query);
        let mut out = Vec::new();
        self.chunks(centre, radius, query, lookup, &mut |zone, first, bytes| {
            for (k, c) in bytes.chunks_exact(RECORD_BYTES).enumerate() {
                let raw = Ucac4Record::decode(c.try_into().expect("exact chunk"));
                if let Some(star) = test.star(&raw) {
                    out.push(Ucac4Hit {
                        zone,
                        record: (first + k as u64 + 1) as u32,
                        raw,
                        star,
                    });
                }
            }
        })?;
        out.sort_by_key(|h| (h.zone, h.record));
        out.dedup_by_key(|h| (h.zone, h.record));
        Ok(out)
    }

    fn chunks(
        &self,
        centre: SkyPos,
        radius: f64,
        query: &StarQuery,
        lookup: Lookup,
        visit: &mut dyn FnMut(u16, u64, Vec<u8>),
    ) -> Result<(), CatalogError> {
        if !radius.is_finite() || !(0.0..std::f64::consts::PI).contains(&radius) {
            return Err(CatalogError::InvalidQuery(format!("radius {radius} rad")));
        }
        let r = radius + self.motion_pad(query) + 1e-9;
        let (ra0, dec0) = (centre.ra.to_degrees(), centre.dec.to_degrees());
        let r_deg = r.to_degrees();
        let z_lo = zone_number(dec0 - r_deg);
        let z_hi = zone_number(dec0 + r_deg);
        let bins = ra_bins(ra0, dec0, r);

        for zone in z_lo..=z_hi {
            let Some((mut f, path)) = self.open_zone(zone)? else {
                if query.strict {
                    return Err(CatalogError::NotFound(format!(
                        "zone file {} (needed for this cone)",
                        self.zone_path(zone).display()
                    )));
                }
                continue;
            };
            let (ranges, searched): (Vec<Option<(u64, u64)>>, bool) = match (lookup, &self.index) {
                (Lookup::Index, Some(idx)) => (
                    bins.iter()
                        .map(|&(lo, hi)| Some(idx.range(zone, lo, hi)))
                        .collect(),
                    false,
                ),
                (Lookup::FullScan, _) => (vec![None], false),
                (Lookup::RaSearch | Lookup::Index, _) => {
                    let records = record_count(&f, &path)?;
                    let mut ranges = Vec::with_capacity(bins.len());
                    for &(lo, hi) in &bins {
                        ranges.push(Some(ra_range(&mut f, &path, records, lo, hi)?));
                    }
                    (ranges, true)
                }
            };
            for range in ranges {
                let bytes = read_range(&mut f, &path, range)?;
                if searched {
                    check_ra_order(&bytes, &path)?;
                }
                visit(zone, range.map_or(0, |r| r.0), bytes);
            }
        }
        Ok(())
    }

    fn motion_pad(&self, query: &StarQuery) -> f64 {
        query
            .epoch
            .map(|e| e.years_since(Epoch::J2000).abs() * self.max_pm_mas_yr * crate::propagate::MAS)
            .unwrap_or(0.0)
    }

    fn brightest(
        &self,
        centre: SkyPos,
        radius: f64,
        query: &StarQuery,
        n: usize,
    ) -> Result<Vec<Star>, CatalogError> {
        struct Candidate {
            mag: f32,
            id: u32,
            zone: u16,
            record: u32,
            chunk: u32,
            at: u32,
        }
        let order = |a: &Candidate, b: &Candidate| {
            a.mag
                .total_cmp(&b.mag)
                .then(a.id.cmp(&b.id))
                .then((a.zone, a.record).cmp(&(b.zone, b.record)))
        };
        let mut chunks: Vec<Vec<u8>> = Vec::new();
        let mut candidates: Vec<Candidate> = Vec::new();
        self.chunks(
            centre,
            radius,
            query,
            self.lookup(),
            &mut |zone, first, bytes| {
                for (k, c) in bytes.chunks_exact(RECORD_BYTES).enumerate() {
                    let raw = Ucac4Record::decode(c.try_into().expect("exact chunk"));
                    let Some(mag) = raw.mag(query.band).filter(|m| m.is_finite()) else {
                        continue;
                    };
                    if query.mag_max.is_some_and(|limit| mag > limit) {
                        continue;
                    }
                    candidates.push(Candidate {
                        mag,
                        id: raw.rnm,
                        zone,
                        record: (first + k as u64 + 1) as u32,
                        chunk: chunks.len() as u32,
                        at: k as u32,
                    });
                }
                chunks.push(bytes);
            },
        )?;
        let test = ConeTest::new(self, centre, radius, query);
        let mut out = Vec::new();
        let mut last = None;
        let mut start = 0;
        let mut batch = n.saturating_add(n / 2).saturating_add(64);
        while start < candidates.len() && out.len() < n {
            let end = start.saturating_add(batch).min(candidates.len());
            if end < candidates.len() {
                candidates[start..].select_nth_unstable_by(end - start, order);
            }
            candidates[start..end].sort_unstable_by(order);
            for c in &candidates[start..end] {
                if last == Some((c.zone, c.record)) {
                    continue;
                }
                last = Some((c.zone, c.record));
                let at = c.at as usize * RECORD_BYTES;
                let bytes = &chunks[c.chunk as usize][at..at + RECORD_BYTES];
                let raw = Ucac4Record::decode(bytes.try_into().expect("one record"));
                if let Some(star) = test.star(&raw) {
                    out.push(star);
                    if out.len() == n {
                        break;
                    }
                }
            }
            start = end;
            batch = batch.saturating_mul(2);
        }
        Ok(out)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Lookup {
    Index,
    RaSearch,
    FullScan,
}

struct ConeTest<'a> {
    catalog: &'a Ucac4Catalog,
    query: &'a StarQuery,
    centre: SkyPos,
    centre_vec: [f64; 3],
    radius: f64,
    reject_below: f64,
    sure: Option<(f64, f64)>,
}

const EDGE_BAND: f64 = 1e-9;
const CHORD_TEST_MAX_RADIUS: f64 = 1.0;

impl<'a> ConeTest<'a> {
    fn new(catalog: &'a Ucac4Catalog, centre: SkyPos, radius: f64, query: &'a StarQuery) -> Self {
        let reach = radius + catalog.motion_pad(query) + EDGE_BAND;
        let chord2 = |angle: f64| (2.0 * (angle / 2.0).psin()).ppowi(2);
        Self {
            catalog,
            query,
            centre,
            centre_vec: centre.to_unit_vector().0,
            radius,
            reject_below: if reach < std::f64::consts::PI {
                reach.pcos() - 1e-12
            } else {
                f64::NEG_INFINITY
            },
            sure: (radius <= CHORD_TEST_MAX_RADIUS).then(|| {
                (
                    chord2((radius - EDGE_BAND).max(0.0)),
                    chord2(radius + EDGE_BAND),
                )
            }),
        }
    }

    fn star(&self, raw: &Ucac4Record) -> Option<Star> {
        let q = self.query;
        if let Some(limit) = q.mag_max {
            if !raw.mag(q.band).is_some_and(|m| m <= limit) {
                return None;
            }
        }
        let hpm = raw
            .is_high_pm()
            .then(|| self.catalog.hpm.get(&raw.rnm).copied())
            .flatten();
        let pm = raw.pm_mas_yr(hpm);
        if q.require_proper_motion && !(pm[0].is_some() && pm[1].is_some()) {
            return None;
        }
        let j2000 = raw.position_j2000();
        let frame = LocalFrame::at(j2000);
        let c = self.centre_vec;
        if frame.u0[0] * c[0] + frame.u0[1] * c[1] + frame.u0[2] * c[2] < self.reject_below {
            return None;
        }
        let target = q.epoch.unwrap_or(Epoch::J2000);
        let (pos, v) = propagate_in(
            &frame,
            j2000,
            [pm[0].unwrap_or(0.0), pm[1].unwrap_or(0.0)],
            Epoch::J2000,
            target,
        );
        let exact = || {
            let outside = self.centre.separation(pos) > self.radius;
            !outside
        };
        let inside = match self.sure {
            Some((inside, outside)) => {
                let d2 = (v[0] - c[0]).ppowi(2) + (v[1] - c[1]).ppowi(2) + (v[2] - c[2]).ppowi(2);
                d2 < inside || (d2 <= outside && exact())
            }
            None => exact(),
        };
        inside.then(|| raw.star_at(pos, pm, hpm, target, q.band))
    }
}

fn record_count(f: &File, path: &Path) -> Result<u64, CatalogError> {
    let len = f
        .metadata()
        .map_err(|source| CatalogError::Io {
            path: path.into(),
            source,
        })?
        .len();
    if len % RECORD_BYTES as u64 != 0 {
        return Err(CatalogError::Malformed {
            path: path.into(),
            reason: format!("size {len} is not a multiple of {RECORD_BYTES}"),
        });
    }
    Ok(len / RECORD_BYTES as u64)
}

fn ra_range(
    f: &mut File,
    path: &Path,
    records: u64,
    lo: usize,
    hi: usize,
) -> Result<(u64, u64), CatalogError> {
    let start = first_at_or_after(f, path, records, lo as u64 * RA_BIN_MAS)?;
    let end = if hi + 1 >= N_RA_BINS {
        records
    } else {
        first_at_or_after(f, path, records, (hi as u64 + 1) * RA_BIN_MAS)?
    };
    Ok((start, end.max(start)))
}

fn first_at_or_after(
    f: &mut File,
    path: &Path,
    records: u64,
    ra_mas: u64,
) -> Result<u64, CatalogError> {
    const BLOCK: u64 = 64;
    let io = |source| CatalogError::Io {
        path: path.into(),
        source,
    };
    let ra = |c: &[u8]| u64::from(u32::from_le_bytes([c[0], c[1], c[2], c[3]]));
    let (mut lo, mut hi) = (0u64, records);
    let mut one = [0u8; 4];
    while hi - lo > BLOCK {
        let mid = lo + (hi - lo) / 2;
        read_at(f, mid * RECORD_BYTES as u64, &mut one).map_err(io)?;
        if ra(&one) < ra_mas {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    let mut block = vec![0u8; (hi - lo) as usize * RECORD_BYTES];
    read_at(f, lo * RECORD_BYTES as u64, &mut block).map_err(io)?;
    let before = block
        .chunks_exact(RECORD_BYTES)
        .position(|c| ra(c) >= ra_mas)
        .unwrap_or((hi - lo) as usize);
    Ok(lo + before as u64)
}

fn read_at(f: &mut File, offset: u64, buf: &mut [u8]) -> std::io::Result<()> {
    #[cfg(unix)]
    {
        std::os::unix::fs::FileExt::read_exact_at(f, buf, offset)
    }
    #[cfg(not(unix))]
    {
        f.seek(SeekFrom::Start(offset))?;
        f.read_exact(buf)
    }
}

fn check_ra_order(bytes: &[u8], path: &Path) -> Result<(), CatalogError> {
    let ra = |c: &[u8]| u32::from_le_bytes([c[0], c[1], c[2], c[3]]);
    let mut records = bytes.chunks_exact(RECORD_BYTES).map(ra);
    let Some(mut before) = records.next() else {
        return Ok(());
    };
    for next in records {
        if next < before {
            return Err(CatalogError::Malformed {
                path: path.into(),
                reason: "records are not in RA order, which a search without u4index.unf needs"
                    .into(),
            });
        }
        before = next;
    }
    Ok(())
}

fn read_range(
    f: &mut File,
    path: &Path,
    range: Option<(u64, u64)>,
) -> Result<Vec<u8>, CatalogError> {
    let io = |source| CatalogError::Io {
        path: path.into(),
        source,
    };
    let mut buf = Vec::new();
    match range {
        None => {
            f.seek(SeekFrom::Start(0)).map_err(io)?;
            f.read_to_end(&mut buf).map_err(io)?;
            if buf.len() % RECORD_BYTES != 0 {
                return Err(CatalogError::Malformed {
                    path: path.into(),
                    reason: format!("size {} is not a multiple of {RECORD_BYTES}", buf.len()),
                });
            }
        }
        Some((start, end)) => {
            if end == start {
                return Ok(buf);
            }
            let available = f.metadata().map_err(io)?.len() / RECORD_BYTES as u64;
            if end > available {
                return Err(CatalogError::Malformed {
                    path: path.into(),
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
                    path: path.into(),
                    reason: format!(
                        "records {start}..{end} not readable ({e}); index and zone disagree"
                    ),
                })?;
        }
    }
    Ok(buf)
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
        let test = ConeTest::new(self, centre, radius, query);
        let mut out: Vec<(u16, u64, Star)> = Vec::new();
        self.chunks(
            centre,
            radius,
            query,
            self.lookup(),
            &mut |zone, first, bytes| {
                for (k, c) in bytes.chunks_exact(RECORD_BYTES).enumerate() {
                    let raw = Ucac4Record::decode(c.try_into().expect("exact chunk"));
                    if let Some(star) = test.star(&raw) {
                        out.push((zone, first + k as u64, star));
                    }
                }
            },
        )?;
        out.sort_by_key(|&(zone, record, _)| (zone, record));
        out.dedup_by_key(|&mut (zone, record, _)| (zone, record));
        Ok(out.into_iter().map(|(_, _, star)| star).collect())
    }

    fn cone_brightest(
        &self,
        centre: SkyPos,
        radius: f64,
        query: &StarQuery,
        n: usize,
    ) -> Result<Vec<Star>, CatalogError> {
        self.brightest(centre, radius, query, n)
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
            ra: 540_000_000,
            spd: 331_200_000,
            magm: 12_345,
            maga: 12_400,
            sigmag: 3,
            objt: 0,
            cdf: 0,
            sigra: -113,
            sigdc: -108,
            na1: 10,
            nu1: 9,
            cu1: 4,
            cepra: 9_850,
            cepdc: 9_900,
            pmrac: 1_234,
            pmdc: -567,
            sigpmr: -108,
            sigpmd: 125,
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

    fn synthetic(name: &str, spoil: bool) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(format!("fayde-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(dir.join("u4b")).unwrap();
        let mut rnm = 0u32;
        for zone in 460u16..=462 {
            let spd0 = (u32::from(zone) - 1) * 720_000;
            let mut ras: Vec<u32> = (0..200u32)
                .map(|i| 532_800_000 + i * 72_000 + (i * 7919 + u32::from(zone)) % 50_000)
                .chain((0..40u32).map(|i| i * 90_000 + 13))
                .chain((0..40u32).map(|i| 1_295_999_000 - i * 90_000))
                .collect();
            ras.sort_unstable();
            if spoil {
                ras.swap(100, 180);
            }
            let mut bytes = Vec::new();
            for (k, ra) in ras.into_iter().enumerate() {
                rnm += 1;
                let k = k as u32;
                bytes.extend(encode(&Ucac4Record {
                    ra,
                    spd: spd0 + (k * 104_729) % 720_000,
                    magm: if k % 11 == 0 {
                        MISSING_MAG
                    } else {
                        8_000 + ((k * 37) % 80) as i16 * 100
                    },
                    maga: MISSING_MAG,
                    pmrac: ((k % 21) as i16 - 10) * 400,
                    pmdc: ((k % 13) as i16 - 6) * 700,
                    rnm,
                    ..sample()
                }));
            }
            std::fs::write(dir.join("u4b").join(zone_file(zone)), bytes).unwrap();
        }
        dir
    }

    #[test]
    fn search_without_index_equals_full_scan() {
        let dir = synthetic("ra-search", false);
        let cat = Ucac4Catalog::open(&dir).unwrap();
        assert!(!cat.has_index());
        let mut found = 0;
        for (i, (ra, dec, r)) in [
            (150.0, 2.1, 0.5),
            (150.0, 2.1, 0.03),
            (148.01, 2.3, 0.2),
            (151.99, 1.9, 0.2),
            (0.01, 2.1, 0.4),
            (359.98, 2.1, 0.4),
            (0.0, 2.1, 0.0005),
            (200.0, 2.1, 0.3),
            (150.0, 2.1, 2.9),
        ]
        .into_iter()
        .enumerate()
        {
            let centre = SkyPos::from_degrees(ra, dec).unwrap();
            let q = StarQuery {
                epoch: [None, Epoch::julian_year(2021.3), Epoch::julian_year(1960.0)][i % 3],
                mag_max: (i % 4 == 3).then_some(12.0),
                require_proper_motion: i % 5 == 4,
                strict: false,
                ..StarQuery::default()
            };
            let r = f64::to_radians(r);
            let by_search = cat.cone_hits(centre, r, &q).unwrap();
            assert_eq!(
                by_search,
                cat.cone_hits_full_scan(centre, r, &q).unwrap(),
                "({ra}, {dec}, {r})"
            );
            let mut plain = Vec::new();
            for zone in 460..=462 {
                for raw in cat.read_zone(zone).unwrap().unwrap() {
                    let star = raw.to_star(None, q.epoch, q.band);
                    let outside = centre.separation(star.pos) > r;
                    let keep = !outside
                        && q.mag_max.is_none_or(|m| star.mag.is_some_and(|v| v <= m))
                        && (!q.require_proper_motion || star.pm_mas_yr.iter().all(Option::is_some));
                    if keep {
                        plain.push(star);
                    }
                }
            }
            let stars = cat.cone(centre, r, &q).unwrap();
            assert_eq!(stars, plain, "({ra}, {dec}, {r})");
            assert_eq!(stars, by_search.iter().map(|h| h.star).collect::<Vec<_>>());
            for n in [0, 1, 5, 60, usize::MAX] {
                assert_eq!(
                    cat.cone_brightest(centre, r, &q, n).unwrap(),
                    crate::brightest_first(stars.clone(), n),
                    "({ra}, {dec}, {r}) brightest {n}"
                );
            }
            found += stars.len();
        }
        assert!(found > 300, "only {found} stars in the test cones");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn zone_out_of_ra_order_is_an_error_without_index() {
        let dir = synthetic("ra-order", true);
        let cat = Ucac4Catalog::open(&dir).unwrap();
        let centre = SkyPos::from_degrees(150.0, 2.1).unwrap();
        let lax = StarQuery {
            strict: false,
            ..StarQuery::default()
        };
        let err = cat.cone(centre, 2f64.to_radians(), &lax).unwrap_err();
        assert!(
            matches!(&err, CatalogError::Malformed { reason, .. } if reason.contains("RA order")),
            "{err}"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn damaged_index_is_an_error_not_an_allocation() {
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
        assert!(at(Some([4000.0, -100.0]))[0] < 200.0);
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
        r.sigpmd = 127;
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
        assert_eq!(ra_bins(0.1, 0.0, 0.01), vec![(1438, 1439), (0, 2)]);
        assert_eq!(ra_bins(10.0, 89.9, 0.01), vec![(0, 1439)]);
    }
}
