//! Index file format, version 1. All integers little-endian.
//!
//! ```text
//! magic        8  b"FAYDEIDX"
//! version      u32
//! meta_len     u32
//! meta         meta_len bytes, JSON (IndexMeta)
//! star_count   u64
//! stars        star_count × { x f64, y f64, z f64, mag f32, id u64, flags u8 } (37 bytes;
//!              flags bit 0 = pattern tier). Version 1 files lack `flags`; all their stars
//!              count as pattern stars.
//! slot_count   u64
//! slots        slot_count × { s0 u32, s1 u32, s2 u32, s3 u32, check u16, largest f32 } (22 bytes)
//! hash         u64   FNV-1a 64 of every preceding byte (Fowler–Noll–Vo; IETF draft-eastlake-fnv)
//! ```

use std::{
    fmt,
    io::{BufReader, BufWriter, Read, Seek, SeekFrom, Write},
    path::{Path, PathBuf},
};

use crate::{
    pattern::{PatternTable, Slot},
    Index, IndexMeta, IndexStar, FORMAT_VERSION,
};

const MAGIC: &[u8; 8] = b"FAYDEIDX";
/// Bytes of one slot record.
const SLOT_BYTES: usize = 22;
/// Chunk size for streamed hashing and serialization.
const CHUNK: usize = 1 << 16;
/// Buffer size of file readers and writers.
const FILE_BUFFER: usize = 1 << 20;

/// Errors reading or writing an index file.
#[derive(Debug)]
pub enum IndexFileError {
    /// I/O failure.
    Io(PathBuf, std::io::Error),
    /// Not a Fayde index, wrong version, truncated, or corrupted.
    Invalid(PathBuf, String),
}

impl fmt::Display for IndexFileError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(p, e) => write!(f, "{}: {e}", p.display()),
            Self::Invalid(p, s) => write!(f, "{}: {s}", p.display()),
        }
    }
}

impl std::error::Error for IndexFileError {}

struct Fnv(u64);

impl Fnv {
    fn new() -> Self {
        Self(0xcbf2_9ce4_8422_2325)
    }
    fn update(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.0 ^= u64::from(b);
            self.0 = self.0.wrapping_mul(0x0100_0000_01b3);
        }
    }
}

/// Feed the serialized index to `sink` in chunks: every byte the trailing hash covers, in file
/// order. Nothing is materialized. The slot table of an all-sky index is tens of gigabytes, and
/// a copy of it used to be the largest single allocation of a build.
fn emit<E>(index: &Index, mut sink: impl FnMut(&[u8]) -> Result<(), E>) -> Result<(), E> {
    let meta = serde_json::to_vec(&index.meta).expect("metadata serializes");
    let mut b: Vec<u8> = Vec::with_capacity(CHUNK + 64);
    b.extend_from_slice(MAGIC);
    b.extend(FORMAT_VERSION.to_le_bytes());
    b.extend((meta.len() as u32).to_le_bytes());
    b.extend(&meta);
    b.extend((index.stars.len() as u64).to_le_bytes());
    for s in &index.stars {
        for c in s.vec {
            b.extend(c.to_le_bytes());
        }
        b.extend(s.mag.to_le_bytes());
        b.extend(s.id.to_le_bytes());
        b.push(u8::from(s.pattern));
        if b.len() >= CHUNK {
            sink(&b)?;
            b.clear();
        }
    }
    b.extend((index.table.slots().len() as u64).to_le_bytes());
    for s in index.table.slots() {
        for id in s.stars {
            b.extend(id.to_le_bytes());
        }
        b.extend(s.check.to_le_bytes());
        b.extend(s.largest.to_le_bytes());
        if b.len() >= CHUNK {
            sink(&b)?;
            b.clear();
        }
    }
    sink(&b)
}

pub(crate) fn content_hash(index: &Index) -> u64 {
    let mut h = Fnv::new();
    let hashed: Result<(), std::convert::Infallible> = emit(index, |bytes| {
        h.update(bytes);
        Ok(())
    });
    match hashed {
        Ok(()) => h.0,
        Err(never) => match never {},
    }
}

/// Why parsing stopped.
enum ParseError {
    Io(std::io::Error),
    Invalid(String),
}

/// A reader over the hashed content of an index file, counting what is left of it.
struct Input<R> {
    r: R,
    left: u64,
}

impl<R: Read> Input<R> {
    fn fill(&mut self, out: &mut [u8]) -> Result<(), ParseError> {
        if self.left < out.len() as u64 {
            return Err(ParseError::Invalid("truncated".into()));
        }
        self.r.read_exact(out).map_err(ParseError::Io)?;
        self.left -= out.len() as u64;
        Ok(())
    }
    fn arr<const N: usize>(&mut self) -> Result<[u8; N], ParseError> {
        let mut b = [0u8; N];
        self.fill(&mut b)?;
        Ok(b)
    }
    fn u32(&mut self) -> Result<u32, ParseError> {
        self.arr().map(u32::from_le_bytes)
    }
    fn u64(&mut self) -> Result<u64, ParseError> {
        self.arr().map(u64::from_le_bytes)
    }
    /// How many `size`-byte records the remaining content can hold.
    fn room_for(&self, size: usize) -> usize {
        usize::try_from(self.left / size as u64).unwrap_or(usize::MAX)
    }
}

/// Parse the content of an index file (everything before the trailing hash, which the caller
/// has verified), reading it once, front to back.
fn parse<R: Read>(r: &mut Input<R>) -> Result<Index, ParseError> {
    let invalid = ParseError::Invalid;
    r.arr::<8>()?; // magic, checked by the caller
    let version = r.u32()?;
    if version != FORMAT_VERSION && version != 1 {
        return Err(invalid(format!(
            "format version {version}, this build reads {FORMAT_VERSION}"
        )));
    }
    let meta_len = r.u32()? as usize;
    if r.left < meta_len as u64 {
        return Err(invalid("truncated".into()));
    }
    let mut meta = vec![0u8; meta_len];
    r.fill(&mut meta)?;
    let meta: IndexMeta =
        serde_json::from_slice(&meta).map_err(|e| invalid(format!("metadata: {e}")))?;
    let n = r.u64()? as usize;
    let mut stars = Vec::with_capacity(n.min(r.room_for(36)));
    for _ in 0..n {
        let x = f64::from_le_bytes(r.arr()?);
        let y = f64::from_le_bytes(r.arr()?);
        let z = f64::from_le_bytes(r.arr()?);
        let mag = f32::from_le_bytes(r.arr()?);
        let id = r.u64()?;
        let pattern = if version >= 2 {
            r.arr::<1>()?[0] & 1 == 1
        } else {
            true
        };
        stars.push(IndexStar {
            vec: [x, y, z],
            mag,
            id,
            pattern,
        });
    }
    let m = r.u64()? as usize;
    let mut slots = Vec::with_capacity(m.min(r.room_for(SLOT_BYTES)));
    for _ in 0..m {
        let b = r.arr::<SLOT_BYTES>()?;
        let word = |i: usize| u32::from_le_bytes([b[i], b[i + 1], b[i + 2], b[i + 3]]);
        let ids = [word(0), word(4), word(8), word(12)];
        if ids[0] != u32::MAX && ids.iter().any(|&i| i as usize >= n) {
            return Err(invalid(
                "pattern refers to a star outside the star table".into(),
            ));
        }
        slots.push(Slot {
            stars: ids,
            check: u16::from_le_bytes([b[16], b[17]]),
            largest: f32::from_le_bytes([b[18], b[19], b[20], b[21]]),
        });
    }
    if r.left != 0 {
        return Err(invalid("trailing bytes after slot table".into()));
    }
    let table = PatternTable::from_slots(slots);
    if table.len() != meta.pattern_count || n != meta.star_count {
        return Err(invalid("counts disagree with metadata".into()));
    }
    // A matching checksum proves the bytes are as written, not that they are usable: check the
    // invariants lookups rely on.
    if !(2..=crate::pattern::MAX_BINS).contains(&meta.bins) {
        return Err(invalid(format!(
            "bins = {} outside 2..={}",
            meta.bins,
            crate::pattern::MAX_BINS
        )));
    }
    Ok(Index::new(meta, stars, table))
}

impl Index {
    /// Write to a file (atomically: a temporary file in the same directory, then rename). The
    /// bytes are streamed; memory use does not grow with the index.
    pub fn write(&self, path: impl AsRef<Path>) -> Result<(), IndexFileError> {
        let path = path.as_ref();
        let (file, tmp) = fayde_core::fsutil::create_temp_beside(path)
            .map_err(|e| IndexFileError::Io(path.into(), e))?;
        let write = || -> std::io::Result<()> {
            let mut w = BufWriter::with_capacity(FILE_BUFFER, file);
            let mut h = Fnv::new();
            emit(self, |bytes| {
                h.update(bytes);
                w.write_all(bytes)
            })?;
            w.write_all(&h.0.to_le_bytes())?;
            let file = w.into_inner().map_err(|e| e.into_error())?;
            file.sync_all()?;
            std::fs::rename(&tmp, path)
        };
        // A failed write must not leave the temporary file behind: it can be tens of GB.
        write().map_err(|e| {
            let _ = std::fs::remove_file(&tmp);
            IndexFileError::Io(path.into(), e)
        })
    }

    /// Read and verify a file. Two streamed passes: the hash first, so no count in the file is
    /// trusted before the file is known to be intact, then the content. The file is never held
    /// in memory next to the index built from it.
    pub fn read(path: impl AsRef<Path>) -> Result<Self, IndexFileError> {
        let path = path.as_ref();
        let io = |e| IndexFileError::Io(path.into(), e);
        let invalid = |s: &str| IndexFileError::Invalid(path.into(), s.into());
        let mut file = std::fs::File::open(path).map_err(io)?;
        let len = file.metadata().map_err(io)?.len();
        if len < 8 + 8 {
            return Err(invalid("not a Fayde index file"));
        }
        let content = len - 8;
        let mut h = Fnv::new();
        let mut tail = [0u8; 8];
        {
            let mut r = BufReader::with_capacity(FILE_BUFFER, &mut file);
            let mut buf = vec![0u8; CHUNK];
            let mut left = content;
            while left > 0 {
                let n = left.min(CHUNK as u64) as usize;
                r.read_exact(&mut buf[..n]).map_err(io)?;
                if left == content && &buf[..8] != MAGIC {
                    return Err(invalid("not a Fayde index file"));
                }
                h.update(&buf[..n]);
                left -= n as u64;
            }
            r.read_exact(&mut tail).map_err(io)?;
        }
        if h.0 != u64::from_le_bytes(tail) {
            return Err(invalid("content hash mismatch (corrupted or truncated)"));
        }
        file.seek(SeekFrom::Start(0)).map_err(io)?;
        let mut input = Input {
            r: BufReader::with_capacity(FILE_BUFFER, file),
            left: content,
        };
        parse(&mut input).map_err(|e| match e {
            ParseError::Io(e) => io(e),
            ParseError::Invalid(s) => invalid(&s),
        })
    }

    /// Parse and verify file bytes.
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, String> {
        if bytes.len() < 8 + 8 || &bytes[..8] != MAGIC {
            return Err("not a Fayde index file".into());
        }
        let (content, tail) = bytes.split_at(bytes.len() - 8);
        let stored = u64::from_le_bytes(tail.try_into().expect("8 bytes"));
        let mut h = Fnv::new();
        h.update(content);
        if h.0 != stored {
            return Err("content hash mismatch (corrupted or truncated)".into());
        }
        let mut input = Input {
            r: content,
            left: content.len() as u64,
        };
        parse(&mut input).map_err(|e| match e {
            ParseError::Io(e) => e.to_string(),
            ParseError::Invalid(s) => s,
        })
    }
}
