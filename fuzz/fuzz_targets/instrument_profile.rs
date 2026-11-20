//! Any text given as an instrument profile: the reader returns a profile or an error.
#![no_main]

use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    if let Ok(text) = std::str::from_utf8(data) {
        let _ = fayde_fits::InstrumentProfile::from_json(text);
    }
});
