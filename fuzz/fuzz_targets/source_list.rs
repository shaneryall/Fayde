//! Any bytes given as a source-list file: the reader returns sources or an error. The reader
//! takes a path, so each input is written to a file of this process first.
#![no_main]

use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let path = std::env::temp_dir().join(format!("fayde-fuzz-{}.fits", std::process::id()));
    if std::fs::write(&path, data).is_ok() {
        let _ = fayde_fits::read_sources(&path);
    }
});
