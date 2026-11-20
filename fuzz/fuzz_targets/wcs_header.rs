//! Any text given as FITS header cards: the WCS readers return a WCS or an error.
#![no_main]

use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    if let Ok(text) = std::str::from_utf8(data) {
        let header = fayde_wcs::Header::parse_block(text);
        let _ = fayde_wcs::TanWcs::from_header(&header);
        let _ = fayde_wcs::Wcs::from_header(&header);
    }
});
