use std::path::{Path, PathBuf};
use std::process::Command;

use fayde_catalog::{ucac4::Ucac4Catalog, Catalog, StarQuery};
use fayde_core::{Epoch, PixelPos, Portable, SkyPos};
use fayde_wcs::TanWcs;
use serde_json::Value;

const W: usize = 512;
const SCALE: f64 = 3.0;
const ROLL: f64 = 23.0;
const CENTRE: (f64, f64) = (150.05, 1.98);

fn ucac4_dir() -> PathBuf {
    let dir = std::env::var_os("FAYDE_UCAC4_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data/ucac4"));
    assert!(
        (458..=463).all(|z| dir.join(format!("u4b/z{z}")).is_file()),
        "UCAC4 zones z458–z463 not at {}",
        dir.display()
    );
    dir
}

fn truth() -> TanWcs {
    let s = SCALE / 3600.0;
    let (sr, cr) = ROLL.to_radians().psin_cos();
    TanWcs::new(
        SkyPos::from_degrees(CENTRE.0, CENTRE.1).unwrap(),
        PixelPos::new((W as f64 - 1.0) / 2.0, (W as f64 - 1.0) / 2.0),
        [[-s * cr, s * sr], [s * sr, s * cr]],
    )
    .unwrap()
}

fn render(path: &Path, header: &[String]) {
    let cat = Ucac4Catalog::open(ucac4_dir()).unwrap();
    let wcs = truth();
    let q = StarQuery {
        epoch: Epoch::julian_year(2025.0),
        mag_max: Some(15.5),
        ..StarQuery::default()
    };
    let stars = cat.cone(wcs.crval(), 0.45f64.to_radians(), &q).unwrap();
    let mut seed = 7u64;
    let mut rnd = || {
        seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((seed >> 11) as f64 + 0.5) / (1u64 << 53) as f64
    };
    let mut img = vec![0f64; W * W];
    for v in img.iter_mut() {
        let (u, t) = (rnd(), rnd());
        *v = 500.0 + 5.0 * (-2.0 * u.pln()).sqrt() * (2.0 * std::f64::consts::PI * t).pcos();
    }
    let sig = 3.0 / 2.354_820_045;
    for st in &stars {
        let Some(p) = wcs.sky_to_pixel(st.pos) else {
            continue;
        };
        let flux = 5.0 * 10f64.ppowf(-0.4 * (f64::from(st.mag.unwrap_or(15.0)) - 20.0));
        let (cx, cy) = (p.x.round() as i64, p.y.round() as i64);
        for y in cy - 8..=cy + 8 {
            for x in cx - 8..=cx + 8 {
                if x < 0 || y < 0 || x >= W as i64 || y >= W as i64 {
                    continue;
                }
                let (dx, dy) = (x as f64 - p.x, y as f64 - p.y);
                img[y as usize * W + x as usize] += flux / (2.0 * std::f64::consts::PI * sig * sig)
                    * (-(dx * dx + dy * dy) / (2.0 * sig * sig)).pexp();
            }
        }
    }
    for _ in 0..30 {
        let i = (rnd() * (W * W) as f64) as usize;
        img[i] += 4000.0;
    }
    let mut cards = vec![
        "SIMPLE  =                    T".to_string(),
        "BITPIX  =                  -32".into(),
        "NAXIS   =                    2".into(),
        format!("NAXIS1  = {W:>20}"),
        format!("NAXIS2  = {W:>20}"),
        "DATE-OBS= '2025-01-01T00:00:00'".into(),
    ];
    cards.extend(header.iter().cloned());
    cards.push("END".into());
    let mut out = Vec::new();
    for c in &cards {
        out.extend(format!("{c:<80}").bytes());
    }
    out.resize(out.len().div_ceil(2880) * 2880, b' ');
    for v in &img {
        out.extend((*v as f32).to_be_bytes());
    }
    out.resize(out.len().div_ceil(2880) * 2880, 0);
    std::fs::write(path, out).unwrap();
}

fn fayde(args: &[&str]) -> (Value, std::process::ExitStatus, String) {
    let o = Command::new(env!("CARGO_BIN_EXE_fayde"))
        .args(args)
        .output()
        .unwrap();
    let report = serde_json::from_slice(&o.stdout).unwrap_or(Value::Null);
    (report, o.status, String::from_utf8_lossy(&o.stderr).into())
}

fn check_solution(r: &Value) {
    assert_eq!(r["solved"], true, "{r}");
    let (ra, dec) = (
        r["crval_deg"][0].as_f64().unwrap(),
        r["crval_deg"][1].as_f64().unwrap(),
    );
    let sep = SkyPos::from_degrees(ra, dec)
        .unwrap()
        .separation(truth().crval())
        .to_degrees()
        * 3600.0;
    assert!(sep < 0.5, "centre off by {sep:.2}″");
    let scale = r["pixel_scale_arcsec"].as_f64().unwrap();
    assert!((scale / SCALE - 1.0).abs() < 1e-3, "scale {scale}");
    let roll = r["roll_deg"].as_f64().unwrap();
    assert!((roll.abs() - ROLL).abs() < 0.05, "roll {roll}");
}

fn scratch(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("fayde-image-{tag}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

#[test]
#[ignore = "needs UCAC4 data: set FAYDE_UCAC4_DIR, run with --include-ignored"]
fn solves_a_raw_image_with_command_line_priors() {
    let d = scratch("cli");
    let img = d.join("frame.fits");
    render(&img, &[]);
    let u = ucac4_dir();
    let (r, status, err) = fayde(&[
        "solve",
        "--image",
        img.to_str().unwrap(),
        "--ucac4",
        u.to_str().unwrap(),
        "--ra",
        "150.0",
        "--dec",
        "2.0",
        "--radius",
        "0.2",
        "--scale",
        "3.0",
        "--epoch",
        "2025.0",
    ]);
    assert!(status.success(), "{err}");
    check_solution(&r);
    let det = &r["detection"];
    assert!(det["count"].as_u64().unwrap() >= 30, "{det}");
    assert!(det["classes"]["sharp"].as_u64().unwrap() >= 25, "{det}");
    assert_eq!(det["rejected_narrow"], det["classes"]["sharp"]);
    assert!(
        (det["fwhm_px"].as_f64().unwrap() - 3.0).abs() < 0.4,
        "{det}"
    );
    let _ = std::fs::remove_dir_all(&d);
}

#[test]
#[ignore = "needs UCAC4 data: set FAYDE_UCAC4_DIR, run with --include-ignored"]
fn solves_a_raw_image_from_its_own_header_and_detect_writes_a_usable_list() {
    let d = scratch("hdr");
    let img = d.join("frame.fits");
    let s = SCALE * 1.005 / 3600.0;
    render(
        &img,
        &[
            "CTYPE1  = 'RA---TAN'".into(),
            "CTYPE2  = 'DEC--TAN'".into(),
            format!("CRVAL1  = {:>20}", CENTRE.0 + 0.4),
            format!("CRVAL2  = {:>20}", CENTRE.1 + 0.03),
            format!("CRPIX1  = {:>20}", 256.5 - 0.4 / s),
            format!("CRPIX2  = {:>20}", 256.5),
            format!("CDELT1  = {:>20}", -s),
            format!("CDELT2  = {:>20}", s),
        ],
    );
    let original = std::fs::read(&img).unwrap();
    let block = |cards: &[&str]| {
        let mut bytes: Vec<u8> = cards
            .iter()
            .flat_map(|card| format!("{card:<80}").into_bytes())
            .collect();
        bytes.resize(2880, b' ');
        bytes
    };
    let mut file = block(&[
        "SIMPLE  = T",
        "BITPIX  = 8",
        "NAXIS   = 0",
        "EXTEND  = T",
        "DATE-OBS= '2025-01-01T00:00:00'",
        "CRVAL1  = 12.",
        "END",
    ]);
    let mut cards: Vec<_> = original[..2880]
        .chunks_exact(80)
        .map(|c| std::str::from_utf8(c).unwrap().trim_end())
        .take_while(|&c| c != "END")
        .filter(|c| !c.starts_with("DATE-OBS"))
        .collect();
    cards[0] = "XTENSION= 'IMAGE   '";
    cards.splice(5..5, ["PCOUNT  = 0", "GCOUNT  = 1"]);
    cards.push("END");
    file.extend(block(&cards));
    file.extend_from_slice(&original[2880..]);
    std::fs::write(&img, file).unwrap();
    let u = ucac4_dir();
    let (r, status, err) = fayde(&[
        "solve",
        "--image",
        img.to_str().unwrap(),
        "--ucac4",
        u.to_str().unwrap(),
        "--radius",
        "0.1",
        "--epoch",
        "2025.0",
    ]);
    assert!(status.success(), "{err}");
    check_solution(&r);
    assert_eq!(r["priors"]["pointing"]["origin"], "WCS image centre");
    assert_eq!(r["detection"]["image"]["hdu"], 1);

    let xy = d.join("frame.sources.fits");
    let (det, status, err) = fayde(&[
        "detect",
        "--image",
        img.to_str().unwrap(),
        "-o",
        xy.to_str().unwrap(),
    ]);
    assert!(status.success(), "{err}");
    let n = det["count"].as_u64().unwrap();
    assert_eq!(det["sources"].as_array().unwrap().len() as u64, n);
    let (r2, status, err) = fayde(&[
        "solve",
        "--sources",
        xy.to_str().unwrap(),
        "--ucac4",
        u.to_str().unwrap(),
        "--header",
        img.to_str().unwrap(),
        "--radius",
        "0.1",
        "--epoch",
        "2025.0",
    ]);
    assert!(status.success(), "{err}");
    check_solution(&r2);
    assert_eq!(r2["crval_deg"], r["crval_deg"]);
    let solution = d.join("reference.wcs");
    fayde_fits::write_wcs_file(&solution, &truth().to_cards()).unwrap();
    let frames = d.join("frames.tsv");
    std::fs::write(
        &frames,
        format!("{}\t{}\n", img.display(), solution.display()).repeat(5),
    )
    .unwrap();
    let (profile, status, err) = fayde(&["calibrate", "--solved", frames.to_str().unwrap()]);
    assert!(status.success(), "{err}");
    assert_eq!(profile["calibration"]["frames"], 5);
    let _ = std::fs::remove_dir_all(&d);
}

#[test]
fn image_option_errors() {
    let d = scratch("err");
    let (_, status, err) = fayde(&[
        "solve",
        "--sources",
        "a.sources.fits",
        "--image",
        "b.fits",
        "--ucac4",
        ".",
    ]);
    assert!(!status.success() && err.contains("not both"), "{err}");
    let (_, status, err) = fayde(&[
        "detect",
        "--image",
        d.join("missing.fits").to_str().unwrap(),
    ]);
    assert!(!status.success() && err.contains("missing.fits"), "{err}");
    let _ = std::fs::remove_dir_all(&d);
}
