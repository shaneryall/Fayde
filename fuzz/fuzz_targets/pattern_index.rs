//! Any bytes given as a pattern index file: the reader returns an index or an error.
#![no_main]

use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let _ = fayde_index::Index::from_bytes(data);
});
