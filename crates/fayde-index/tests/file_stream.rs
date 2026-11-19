use fayde_catalog::{CatalogId, Star};
use fayde_core::SkyPos;
use fayde_index::{build, BuildParams, Index, IndexFileError};

fn index() -> Index {
    let mut s = 0x0123_4567_89AB_CDEFu64;
    let mut rnd = || {
        s = s
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (s >> 11) as f64 / (1u64 << 53) as f64
    };
    let stars: Vec<Star> = (0..1500u64)
        .map(|i| Star {
            id: i + 1,
            pos: SkyPos::from_degrees(80.0 + 3.0 * rnd(), -20.0 + 3.0 * rnd()).unwrap(),
            mag: Some(7.0 + 8.0 * rnd() as f32),
            pos_err_mas: [20.0, 20.0],
            pm_mas_yr: [Some(0.0), Some(0.0)],
        })
        .collect();
    let params = BuildParams {
        region: Some([81.5, -18.5, 1.2]),
        ..BuildParams::default()
    };
    let id = CatalogId {
        name: "synthetic",
        release: "0",
        location: String::new(),
    };
    build(&stars, id, &params).unwrap()
}

fn fnv(bytes: &[u8]) -> u64 {
    bytes.iter().fold(0xcbf2_9ce4_8422_2325, |h, &b| {
        (h ^ u64::from(b)).wrapping_mul(0x0100_0000_01b3)
    })
}

fn rehash(bytes: &mut Vec<u8>) {
    let content = bytes.len() - 8;
    let h = fnv(&bytes[..content]);
    bytes.truncate(content);
    bytes.extend(h.to_le_bytes());
}

fn star_count_at(bytes: &[u8]) -> usize {
    16 + u32::from_le_bytes(bytes[12..16].try_into().unwrap()) as usize
}

struct Scratch(std::path::PathBuf);

impl Scratch {
    fn new(tag: &str, bytes: &[u8]) -> Self {
        let path =
            std::env::temp_dir().join(format!("fayde-stream-{tag}-{}.fidx", std::process::id()));
        std::fs::write(&path, bytes).unwrap();
        Self(path)
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.0);
    }
}

fn rejection(tag: &str, bytes: &[u8]) -> String {
    let from_memory = Index::from_bytes(bytes).unwrap_err();
    let file = Scratch::new(tag, bytes);
    match Index::read(&file.0).unwrap_err() {
        IndexFileError::Invalid(_, from_file) => assert_eq!(from_file, from_memory, "{tag}"),
        IndexFileError::Io(_, e) => panic!("{tag}: I/O error {e} instead of {from_memory}"),
    }
    from_memory
}

#[test]
fn round_trip_through_a_file_and_through_bytes() {
    let index = index();
    assert!(index.meta().pattern_count > 1000, "{:?}", index.meta());
    let file = Scratch::new("round", &[]);
    index.write(&file.0).unwrap();
    let bytes = std::fs::read(&file.0).unwrap();
    assert_eq!(
        u64::from_le_bytes(bytes[bytes.len() - 8..].try_into().unwrap()),
        index.content_hash()
    );
    assert_eq!(fnv(&bytes[..bytes.len() - 8]), index.content_hash());
    for back in [
        Index::read(&file.0).unwrap(),
        Index::from_bytes(&bytes).unwrap(),
    ] {
        assert_eq!(back, index);
        assert_eq!(back.content_hash(), index.content_hash());
    }
    index.write(&file.0).unwrap();
    assert_eq!(std::fs::read(&file.0).unwrap(), bytes);
}

#[test]
fn damaged_files_are_rejected_before_parsing() {
    let index = index();
    let file = Scratch::new("good", &[]);
    index.write(&file.0).unwrap();
    let good = std::fs::read(&file.0).unwrap();

    let stars_at = star_count_at(&good);
    for at in [
        9,
        20,
        stars_at + 3,
        stars_at + 100,
        good.len() / 2,
        good.len() - 9,
        good.len() - 1,
    ] {
        let mut bad = good.clone();
        bad[at] ^= 0x01;
        assert!(rejection("flip", &bad).contains("hash"), "byte {at}");
    }
    for keep in [good.len() - 1, good.len() - 100, good.len() / 2, 40, 16] {
        assert!(rejection("short", &good[..keep]).contains("hash"), "{keep}");
    }
    for keep in [15, 8, 0] {
        assert_eq!(rejection("tiny", &good[..keep]), "not a Fayde index file");
    }
    let mut longer = good.clone();
    longer.extend([0u8; 5]);
    assert!(rejection("long", &longer).contains("hash"));
    let mut other = good.clone();
    other[0] = b'X';
    assert_eq!(rejection("magic", &other), "not a Fayde index file");
    assert!(matches!(
        Index::read(file.0.with_extension("missing")),
        Err(IndexFileError::Io(..))
    ));
}

#[test]
fn intact_hash_over_a_broken_structure_is_still_rejected() {
    let index = index();
    let file = Scratch::new("structure", &[]);
    index.write(&file.0).unwrap();
    let good = std::fs::read(&file.0).unwrap();
    let stars_at = star_count_at(&good);
    let n = u64::from_le_bytes(good[stars_at..stars_at + 8].try_into().unwrap()) as usize;
    let slots_at = stars_at + 8 + 37 * n;

    let edit = |f: &dyn Fn(&mut Vec<u8>)| {
        let mut b = good.clone();
        f(&mut b);
        rehash(&mut b);
        b
    };
    let b = edit(&|b| b[stars_at..stars_at + 8].copy_from_slice(&u64::MAX.to_le_bytes()));
    assert_eq!(rejection("stars", &b), "truncated");
    let m = u64::from_le_bytes(good[slots_at..slots_at + 8].try_into().unwrap());
    for count in [m + 1, u64::MAX / 2] {
        let b = edit(&|b| b[slots_at..slots_at + 8].copy_from_slice(&count.to_le_bytes()));
        assert_eq!(rejection("slots", &b), "truncated");
    }
    let b = edit(&|b| b[slots_at..slots_at + 8].copy_from_slice(&(m - 1).to_le_bytes()));
    assert_eq!(rejection("leftover", &b), "trailing bytes after slot table");
    let first_used = (0..m as usize)
        .map(|k| slots_at + 8 + 22 * k)
        .find(|&at| good[at..at + 4] != u32::MAX.to_le_bytes())
        .unwrap();
    let b = edit(&|b| b[first_used + 4..first_used + 8].copy_from_slice(&(n as u32).to_le_bytes()));
    assert_eq!(
        rejection("star", &b),
        "pattern refers to a star outside the star table"
    );
    let b = edit(&|b| b[8..12].copy_from_slice(&9u32.to_le_bytes()));
    assert!(rejection("version", &b).starts_with("format version 9"));
    let b = edit(&|b| b[16] = b'!');
    assert!(rejection("meta", &b).starts_with("metadata:"));
}

#[test]
fn version_1_files_still_read() {
    let index = index();
    let file = Scratch::new("v2", &[]);
    index.write(&file.0).unwrap();
    let v2 = std::fs::read(&file.0).unwrap();
    let stars_at = star_count_at(&v2);
    let n = u64::from_le_bytes(v2[stars_at..stars_at + 8].try_into().unwrap()) as usize;

    let mut v1 = v2[..stars_at + 8].to_vec();
    v1[8..12].copy_from_slice(&1u32.to_le_bytes());
    for k in 0..n {
        let at = stars_at + 8 + 37 * k;
        v1.extend_from_slice(&v2[at..at + 36]);
    }
    v1.extend_from_slice(&v2[stars_at + 8 + 37 * n..]);
    rehash(&mut v1);
    assert_eq!(v1.len(), v2.len() - n);

    let old = Scratch::new("v1", &v1);
    for back in [
        Index::read(&old.0).unwrap(),
        Index::from_bytes(&v1).unwrap(),
    ] {
        assert_eq!(
            back.content_hash(),
            u64::from_le_bytes(v1[v1.len() - 8..].try_into().unwrap())
        );
        assert_eq!(back.stars().len(), n);
        assert!(back.stars().iter().all(|s| s.pattern));
        for (a, b) in back.stars().iter().zip(index.stars()) {
            assert_eq!((a.vec, a.mag, a.id), (b.vec, b.mag, b.id));
        }
        assert_eq!(back.table(), index.table());
        assert_eq!(back.meta(), index.meta());
    }
}
