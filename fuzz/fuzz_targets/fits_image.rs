//! Any bytes given as a FITS file: the reader returns an image or an error, within bounded time
//! and memory.
#![no_main]

use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    if let Ok(file) = fayde_fits::FitsFile::from_bytes(data.to_vec()) {
        let _ = fayde_fits::image_from_file(&file);
    }
});
