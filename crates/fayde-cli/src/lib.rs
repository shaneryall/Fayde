//! The `fayde` command-line tool, as a library: the `fayde` binary and the `fayde` command that
//! the Python package installs both call [`run`], so the two are the same program.
//!
//! ```text
//! fayde index build --ucac4 <dir> --fov <deg> -o <file> [--epoch <year>] [--mag-limit <m>]
//!                   [--stars-per-fov <n>] [--oversampling <x>] [--patterns-per-field <n>]
//!                   [--stars-per-field <n>] [--min-edge <deg>] [--region <ra,dec,r>]
//! fayde index info <file>
//! ```

use std::{collections::HashMap, time::Instant};

use fayde_catalog::ucac4::Ucac4Catalog;
use fayde_index::{build_from_catalog_with_progress, BuildParams, Index};

mod bench;
mod calibrate_cmd;
mod compare;
mod detect_cmd;
mod solve_cmd;

/// Run the `fayde` command with `args`, the arguments after the program name. Prints the
/// command's output and errors, and returns the process exit status: 0 when the command
/// completes, 1 when it fails.
pub fn run(args: &[String]) -> u8 {
    if args.iter().any(|a| a == "--help" || a == "-h") || args.is_empty() {
        let help = match args.first().map(String::as_str) {
            Some("solve" | "verify" | "refit") => REQUEST_HELP,
            Some("index") if args.get(1).is_some_and(|a| a == "build") => INDEX_BUILD_HELP,
            Some("detect") => DETECT_HELP,
            Some("calibrate") => CALIBRATE_HELP,
            _ => USAGE,
        };
        return finish(out(help));
    }
    if args.iter().any(|a| a == "--version" || a == "-V") {
        return finish(out(&format!(
            "fayde {} ({})",
            env!("CARGO_PKG_VERSION"),
            fayde_api::commit()
        )));
    }
    let result = match args
        .iter()
        .map(String::as_str)
        .collect::<Vec<_>>()
        .as_slice()
    {
        ["index", "build", rest @ ..] => index_build(rest),
        ["index", "info", file] => index_info(file),
        ["bench", rest @ ..] => bench::run(rest),
        ["compare", rest @ ..] => compare::run(rest),
        ["detect", rest @ ..] => detect_cmd::run(rest),
        ["calibrate", rest @ ..] => calibrate_cmd::run(rest),
        ["solve", rest @ ..] => solve_cmd::run(rest, solve_cmd::Mode::Solve),
        ["verify", rest @ ..] => solve_cmd::run(rest, solve_cmd::Mode::Verify),
        ["refit", rest @ ..] => solve_cmd::run(rest, solve_cmd::Mode::Refit),
        _ => Err(USAGE.to_string()),
    };
    finish(result)
}

fn finish(result: Result<(), String>) -> u8 {
    match result {
        Ok(()) => 0,
        Err(e) => {
            eprintln!("fayde: {e}");
            1
        }
    }
}

/// Print a line to stdout. A reader that has gone away (`fayde solve ... | head -1`) is not an
/// error, and `println!` would panic on it; any other write failure is reported.
pub(crate) fn out(text: &str) -> Result<(), String> {
    use std::io::Write;
    let mut stdout = std::io::stdout().lock();
    match writeln!(stdout, "{text}").and_then(|()| stdout.flush()) {
        Ok(()) => Ok(()),
        Err(e) if e.kind() == std::io::ErrorKind::BrokenPipe => Ok(()),
        Err(e) => Err(format!("stdout: {e}")),
    }
}

/// Check that `path` can be written, before the work whose result goes there: an output in a
/// missing or read-only directory must not cost a solve or an index build.
pub(crate) fn check_output_path(path: &str) -> Result<(), String> {
    let p = std::path::Path::new(path);
    if p.is_dir() {
        return Err(format!("{path}: is a directory"));
    }
    let dir = match p.parent() {
        Some(d) if !d.as_os_str().is_empty() => d,
        _ => std::path::Path::new("."),
    };
    // A new temporary file beside the output (O_CREAT|O_EXCL: never follows or truncates an
    // existing file or symlink), removed again at once.
    let (_, probe) = fayde_core::fsutil::create_temp_beside(p)
        .map_err(|e| format!("{path}: cannot write in {}: {e}", dir.display()))?;
    let _ = std::fs::remove_file(&probe);
    Ok(())
}

const USAGE: &str = "usage:
  fayde index build --ucac4 <dir> --fov <deg> -o <file> [--epoch <year>] [--mag-limit <m>]
                    [--stars-per-fov <n>] [--oversampling <x>] [--patterns-per-field <n>]
                    [--stars-per-field <n>] [--pattern-stars-per-fov <n>] [--min-edge <deg>]
                    [--region <ra,dec,r>]
  fayde index info <file>
  fayde detect --image <frame.fits> [-o <frame.sources.fits>] [detection options]
  fayde solve  (--sources <frame.sources.fits> | --image <frame.fits>) --ucac4 <dir> [request options]
               [--index <a.fidx>]... [-o <out.wcs>]
  fayde verify (--sources <frame.sources.fits> | --image <frame.fits>) --ucac4 <dir> --wcs <in.wcs> [request options]
  fayde verify --ucac4 <dir> --wcs <in.wcs> (--at <x,y>... | --at-sky <ra,dec>... | --at-sources <file>)
               [--header <frame.fits>] [--epoch <year>] [--pixel-error <px>] [--wcs-sigma <arcsec>]
               [--search-radius <arcsec>] [--max-stars <n>]
  fayde refit  (--sources <frame.sources.fits> | --image <frame.fits>) --ucac4 <dir> --wcs <in.wcs> [request options]
               [-o <out.wcs>]
     detection options (with --image): --fwhm <px> --detect-threshold <sigma> --saturation <adu>
       --gain <e/adu> --max-sources <n>
     request options: --header <frame.fits> --instrument <profile.json> --width <px> --height <px>
       --ra <deg> --dec <deg> --radius <deg> --pointing-sigma <deg>
       --scale <arcsec/px> --scale-sigma <fraction>
       --roll <deg> --roll-sigma <deg> --roll-max <deg> --parity negative|positive
       --trail-length <px> --trail-angle <deg> --epoch <year> --deadline <s> --max-sip-order <n>
       --pixel-error <px> --accept-log-odds <x> --max-hypotheses <n>
       --no-priors --no-pointing --no-header-wcs --search-outside-pointing
       --relax-priors-on-failure
  fayde calibrate --solved <list> [--roll-keyword <KEY>] [--binning <XKEY>[,<YKEY>]]
                  [--name <text>] [--no-header-wcs] [-o <profile.json>]
  fayde compare --ucac4 <dir> --index <a.fidx>[,<b.fidx>] --frames <frames.json> [--limit <n>]
                [--instrument <profile.json>] [--no-header-wcs] [--hinted-only] [--deadline <s>]
  fayde bench --ucac4 <dir> --index <file> [--frames <n>] [--seed <n>] [--size <px>] [--scale <arcsec/px>]
              [--image-mag <m>] [--deadline <s>] [--pattern-sources <n>]";

/// `fayde solve|verify|refit --help`: every request option with its default.
const REQUEST_HELP: &str =
    "fayde solve  (--sources <frame.sources.fits> | --image <frame.fits>) --ucac4 <dir> [options] [--index <a.fidx>]...
             [-o <out.wcs>]
fayde verify (--sources <frame.sources.fits> | --image <frame.fits>) --ucac4 <dir> --wcs <in.wcs> [options]
fayde verify --ucac4 <dir> --wcs <in.wcs> (--at <x,y> | --at-sky <ra,dec> | --at-sources <file>)...
fayde refit  (--sources <frame.sources.fits> | --image <frame.fits>) --ucac4 <dir> --wcs <in.wcs> [options]
             [-o <out.wcs>]

Once the inputs are read, the report (JSON) goes to stdout for every outcome: a solution, an
abstention, or a solver error (with \"error\": true). A mistake on the command line or an
unreadable input is reported on stderr only. The exit status is 1 for any error, 0 otherwise.
Priors come from, in order: these options, the instrument profile, the frame header.

Inputs
  --image <frame.fits>     a FITS image: Fayde finds the stars (see fayde detect --help for the
                           detection options); its header is used unless --header is given
  --sources <frame.sources.fits>
                           a source list from your own detector: a FITS binary table with
                           X, Y (FITS pixels, first pixel centre at 1, 1); optional FLUX,
                           SIGMA or SIGMA_PX, P_STAR, TRUNCATED; --xyls is an older name
  --ucac4 <dir>            UCAC4 directory (holding u4b/, optionally u4i/)
  --index <a.fidx>         pattern index; repeat for several (solve only; needed without a
                           pointing and tight scale)
  --header <frame.fits>    frame header: pointing, scale, epoch, size, binning
  --instrument <json>      instrument profile (scale, parity, roll rule, binning rule,
                           pointing uncertainty); fayde calibrate writes one
  --width, --height <px>   image size              default: source list IMAGEW/IMAGEH, else header
Pointing
  --ra, --dec <deg>        pointing centre         default: header
  --radius <deg>           hard limit              default: profile, else 2
  --pointing-sigma <deg>   1-sigma                 default: profile, else radius/4
                           (a 1-sigma of a sixth of the limit or less: the search starts
                           within 3 sigma and widens to the limit if that finds nothing)
  --no-pointing            drop the pointing only (keep scale, roll, parity)
  --no-header-wcs          do not read the pointing or scale from a WCS in the header (it
                           may be another solver's solution): RA/DEC or OBJCTRA/OBJCTDEC only
  --search-outside-pointing  search beyond the limit if nothing is found inside
Scale, roll, parity
  --scale <arcsec/px>      pixel scale             default: profile, else header
  --scale-sigma <fraction> 1-sigma                 default 0.02 (profile: its own)
  --roll <deg>             roll                    default: profile rule (1-sigma 0.5, limit 3)
  --roll-sigma, --roll-max <deg>  override either source's 1-sigma and hard limit
  --parity negative|positive                       default: profile
  --no-priors              drop every prior
  --relax-priors-on-failure  retry once with wider hard limits
Frame
  --trail-length <px> --trail-angle <deg>   trail of every source (truncated sources)
  --epoch <year>           Julian year, 1800-2200  default: header (mid-exposure), else 2000
Search and fit
  --pixel-error <px>       centroid 1-sigma of sources without SIGMA     default 1.0
  --accept-log-odds <x>    acceptance threshold    default ln(1e9) = 20.7
  --max-hypotheses <n>     hypothesis budget       default 5000000
  --deadline <s>           wall-clock budget from the start            default none
  --max-sip-order <n>      highest SIP order tried (BIC chooses), up to 9   default 3
Output
  -o <out.wcs>             solved WCS with provenance comments (solve, refit). A run that
                           writes none (no solution, or an error once the options are
                           accepted) removes a file already at that path. It may not name
                           an input of the run.
verify --at
  --wcs-sigma <arcsec>     radial 1-sigma of the WCS (a solve's worst_sigma_arcsec) default 0
  --search-radius <arcsec> stars reported within   default 30
  --max-stars <n>          stars per position      default 3";

/// `fayde calibrate --help`.
const CALIBRATE_HELP: &str = "fayde calibrate --solved <list> [options] [-o <profile.json>]

Writes an instrument profile (for --instrument) from frames Fayde has solved: the pixel scale,
the parity, the rule that gives the roll from a header keyword, and how far the header's
pointing is from the solved frame centre. Values are medians of the solutions, so a few poor
solutions do not move them; at least 5 solved frames are needed, and more are better.

A profile needs no other solver. Solve some frames with the pointing from their headers and a
nominal scale (fayde solve --scale <arcsec/px> ... -o frame.wcs), then calibrate from those.

  --solved <list>          text file, one solved frame per line: the frame's FITS file (for
                           its header), a tab, the WCS file of its solution; # starts a comment
  --roll-keyword <KEY>     header keyword with the instrument's own roll, degrees: the profile
                           gets the rule roll = offset + sign x KEY       default: no roll rule
  --binning <XKEY>[,<YKEY>]  header keywords of the binning factors: the profile's scale is
                           then the unbinned one                           default: none
  --name <text>            the profile's name
  --no-header-wcs          measure the header pointing from RA/DEC or OBJCTRA/OBJCTDEC, not a
                           WCS in the header, and make the profile tell the solver the same
  -o <profile.json>        also write the profile to a file

The profile goes to stdout. Its \"calibration\" member holds the numbers behind it: frames
used, spreads, worst deviations.";

/// `fayde detect --help`, and the detection options of `--image`.
const DETECT_HELP: &str = "fayde detect --image <frame.fits> [-o <frame.sources.fits>] [options]

Finds stars (point sources: sidereal tracking or stare) in a FITS image and prints a JSON
report: the settings used, the stars' width, background and noise, every source with its
position (FITS 1-based), error, flux, star probability and class, the streaks traced (ends,
angle, width, flux), and what was set aside (hot pixels, cosmic-ray hits, anything on a
streak). -o also writes the sources as a source list for --sources.
The same detection runs inside `fayde solve|verify|refit --image <frame.fits>`, which also
reads the frame's header from the image unless --header is given.

  --image <frame.fits>        FITS image: BITPIX 8/16/32/64/-32/-64 with BSCALE/BZERO/BLANK,
                              primary HDU or first IMAGE extension (not fpack-compressed)
  --fwhm <px>                 stellar FWHM                     default: measured from the image
  --detect-threshold <sigma>  peak threshold, filtered image   default: one noise peak per
                              frame (about 5 for 4 megapixels)
  --saturation <adu>          saturation level                 default: SATURATE, else the
                              level of flat-topped stars
  --gain <e/adu>              gain, for the stars' own noise   default: GAIN or EGAIN, else none
  --max-sources <n>           brightest sources kept           default 2000
  -o <frame.sources.fits>     also write the sources as a source list (Fayde's FITS table:
                              X, Y, FLUX, SIGMA, P_STAR, TRUNCATED); name it after the frame";

/// `fayde index build --help`.
const INDEX_BUILD_HELP: &str = "fayde index build --ucac4 <dir> --fov <deg> -o <file> [options]

  --fov <deg>                     field of view the index serves (required)
  --epoch <year>                  star positions propagated to       default 2025
  --mag-limit <m>                 faintest magnitude                 default none
  --stars-per-fov <n>             stars kept per field of view       default 20
  --pattern-stars-per-fov <n>     brightest of those used in patterns default 20
  --stars-per-field <n>           stars per lattice field            default 12
  --patterns-per-field <n>        patterns per lattice field         default 50
  --oversampling <x>              lattice fields covering each sky point  default 100
  --min-edge <deg>                shortest pattern edge              default 0.05
  --region <ra,dec,r>             regional index (degrees)           default all sky

Recommended set, sizes and memory: docs/index-set.md. `fayde index info <file>` prints the
content hash.";

/// `--key value` pairs.
/// Parse `--key value` options and `--flag` switches. `known` lists the keys that take a value,
/// `flags` the switches (stored with value "1"). Unknown, duplicated or valueless keys are errors,
/// so a typo such as `--raduis` is reported instead of silently ignored.
pub(crate) fn options(
    args: &[&str],
    known: &[&str],
    flags: &[&str],
) -> Result<HashMap<String, String>, String> {
    let mut out = HashMap::new();
    let mut it = args.iter();
    while let Some(k) = it.next() {
        let key = k
            .strip_prefix("--")
            .or_else(|| k.strip_prefix('-'))
            .ok_or_else(|| format!("unexpected argument '{k}' (see fayde --help)"))?;
        let value = if flags.contains(&key) {
            "1".to_string()
        } else if known.contains(&key) {
            it.next()
                .ok_or_else(|| format!("missing value for {k}"))?
                .to_string()
        } else {
            return Err(format!("unknown option {k} (see fayde --help)"));
        };
        if REPEATABLE.contains(&key) {
            out.entry(key.to_string())
                .and_modify(|v: &mut String| {
                    v.push(LIST_SEP);
                    v.push_str(&value);
                })
                .or_insert(value);
        } else if out.insert(key.to_string(), value).is_some() {
            return Err(format!("{k} given more than once"));
        }
    }
    Ok(out)
}

/// Options that may be given more than once; their values are joined with [`LIST_SEP`].
const REPEATABLE: &[&str] = &["index", "at", "at-sky", "at-sources", "at-xyls"];
pub(crate) const LIST_SEP: char = '\u{1f}';

/// The paths given with `--index`: one per `--index`. For compatibility, a value that is not an
/// existing file is also split on commas (`--index a.fidx,b.fidx`).
pub(crate) fn index_paths(o: &HashMap<String, String>) -> Vec<String> {
    o.get("index")
        .map(|v| {
            v.split(LIST_SEP)
                .flat_map(|p| {
                    if std::path::Path::new(p).exists() {
                        vec![p.to_string()]
                    } else {
                        p.split(',').map(|q| q.trim().to_string()).collect()
                    }
                })
                .filter(|p| !p.is_empty())
                .collect()
        })
        .unwrap_or_default()
}

pub(crate) fn num<T: std::str::FromStr>(
    o: &HashMap<String, String>,
    k: &str,
) -> Result<Option<T>, String> {
    o.get(k)
        .map(|v| v.parse().map_err(|_| format!("--{k}: cannot parse '{v}'")))
        .transpose()
}

/// `--<k>` as a duration in seconds: finite, non-negative and representable.
pub(crate) fn seconds(
    o: &HashMap<String, String>,
    k: &str,
) -> Result<Option<std::time::Duration>, String> {
    num::<f64>(o, k)?
        .map(|v| {
            std::time::Duration::try_from_secs_f64(v)
                .map_err(|_| format!("--{k}: '{v}' is not a valid number of seconds"))
        })
        .transpose()
}

fn index_build(args: &[&str]) -> Result<(), String> {
    let o = options(
        args,
        &[
            "ucac4",
            "o",
            "output",
            "fov",
            "epoch",
            "mag-limit",
            "stars-per-fov",
            "oversampling",
            "patterns-per-field",
            "stars-per-field",
            "pattern-stars-per-fov",
            "min-edge",
            "region",
        ],
        &[],
    )?;
    let dir = o.get("ucac4").ok_or("--ucac4 is required")?;
    if o.contains_key("o") && o.contains_key("output") {
        return Err("-o and --output are the same option: give one".into());
    }
    let out = o
        .get("o")
        .or_else(|| o.get("output"))
        .ok_or("-o is required")?;
    let mut p = BuildParams {
        fov_deg: num(&o, "fov")?.ok_or("--fov is required")?,
        ..BuildParams::default()
    };
    if let Some(v) = num(&o, "epoch")? {
        p.epoch_year = v;
    }
    p.mag_limit = num(&o, "mag-limit")?;
    if let Some(v) = num(&o, "stars-per-fov")? {
        p.stars_per_fov = v;
    }
    if let Some(v) = num(&o, "oversampling")? {
        p.field_oversampling = v;
    }
    if let Some(v) = num(&o, "patterns-per-field")? {
        p.patterns_per_field = v;
    }
    if let Some(v) = num(&o, "stars-per-field")? {
        p.stars_per_field = v;
    }
    if let Some(v) = num(&o, "pattern-stars-per-fov")? {
        p.pattern_stars_per_fov = v;
    }
    if let Some(v) = num(&o, "min-edge")? {
        p.min_edge_deg = v;
    }
    if let Some(r) = o.get("region") {
        let v: Vec<f64> = r
            .split(',')
            .map(|s| s.trim().parse().map_err(|_| format!("--region: '{r}'")))
            .collect::<Result<_, _>>()?;
        let [ra, dec, rad] = v[..] else {
            return Err("--region needs ra,dec,radius".into());
        };
        p.region = Some([ra, dec, rad]);
    }
    // Both checked before the catalogue scan, which takes minutes on the full sky.
    p.validate().map_err(|e| e.to_string())?;
    check_output_path(out)?;
    let cat = Ucac4Catalog::open(dir).map_err(|e| e.to_string())?;
    let t0 = Instant::now();
    let index = build_from_catalog_with_progress(&cat, &p, &mut |m| {
        eprintln!("[{:>7.1}s] {m}", t0.elapsed().as_secs_f64());
    })
    .map_err(|e| e.to_string())?;
    eprintln!("[{:>7.1}s] writing {out}", t0.elapsed().as_secs_f64());
    index.write(out).map_err(|e| e.to_string())?;
    print_info(&index, Some(t0.elapsed().as_secs_f64()))
}

fn index_info(file: &str) -> Result<(), String> {
    let index = Index::read(file).map_err(|e| e.to_string())?;
    print_info(&index, None)
}

fn print_info(index: &Index, build_seconds: Option<f64>) -> Result<(), String> {
    let m = index.meta();
    let info = serde_json::json!({
        "meta": m,
        "content_hash": format!("{:016x}", index.content_hash()),
        "file_bytes_estimate": 32 + m.star_count * 36 + m.table_size * 22,
        "build_seconds": build_seconds,
    });
    out(&serde_json::to_string_pretty(&info).expect("json"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn index_may_be_repeated_and_paths_may_contain_commas() {
        let dir = std::env::temp_dir().join(format!("fayde-cli-test-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let comma = dir.join("a,b.fidx");
        std::fs::write(&comma, b"").unwrap();
        let comma = comma.to_str().unwrap();
        let o = options(
            &["--index", comma, "--index", "x.fidx,y.fidx"],
            &["index"],
            &[],
        )
        .unwrap();
        assert_eq!(index_paths(&o), vec![comma, "x.fidx", "y.fidx"]);
        let _ = std::fs::remove_dir_all(&dir);
        // Other options still refuse repeats.
        assert!(options(&["--ucac4", "a", "--ucac4", "b"], &["ucac4"], &[]).is_err());
    }

    #[test]
    fn output_paths_are_checked_before_the_work() {
        let dir = std::env::temp_dir().join(format!("fayde-cli-out-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let s = |p: std::path::PathBuf| p.to_str().unwrap().to_string();
        assert!(check_output_path(&s(dir.join("a.wcs"))).is_ok());
        // Nothing is left behind by the check.
        assert_eq!(std::fs::read_dir(&dir).unwrap().count(), 0);
        let e = check_output_path(&s(dir.clone())).unwrap_err();
        assert!(e.contains("is a directory"), "{e}");
        let e = check_output_path(&s(dir.join("missing").join("a.wcs"))).unwrap_err();
        assert!(e.contains("cannot write in"), "{e}");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn deadlines_must_be_valid_durations() {
        let o = |v: &str| HashMap::from([("deadline".to_string(), v.to_string())]);
        for bad in ["-1", "NaN", "inf", "1e300"] {
            assert!(seconds(&o(bad), "deadline").is_err(), "{bad}");
        }
        assert_eq!(
            seconds(&o("2.5"), "deadline").unwrap(),
            Some(std::time::Duration::from_millis(2500))
        );
        assert_eq!(seconds(&HashMap::new(), "deadline").unwrap(), None);
    }
}
