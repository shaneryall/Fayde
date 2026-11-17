//! A fast hash for keys that are already well-spread numbers.
//!
//! The standard `HashMap` hashes with SipHash-1-3, which is built to resist keys chosen to
//! collide. The maps here are keyed by numbers the solver made itself (grid cells, star
//! identifiers and coordinates) and are looked up thousands of times per solve, where SipHash
//! was a measurable part of the time. This is the multiply-rotate hash of Firefox and the Rust
//! compiler ("Fx hash"; the `rustc-hash` crate). It has no random seed: the same keys always
//! hash the same, and nothing here depends on a map's iteration order.

use std::hash::{BuildHasherDefault, Hasher};

/// The hasher. See the module documentation.
#[derive(Default, Clone, Copy)]
pub(crate) struct Mix(u64);

impl Mix {
    fn add(&mut self, word: u64) {
        self.0 = (self.0.rotate_left(5) ^ word).wrapping_mul(0x517c_c1b7_2722_0a95);
    }
}

impl Hasher for Mix {
    fn finish(&self) -> u64 {
        self.0
    }
    fn write(&mut self, bytes: &[u8]) {
        for chunk in bytes.chunks(8) {
            let mut word = [0u8; 8];
            word[..chunk.len()].copy_from_slice(chunk);
            self.add(u64::from_le_bytes(word));
        }
    }
    fn write_u32(&mut self, v: u32) {
        self.add(u64::from(v));
    }
    fn write_u64(&mut self, v: u64) {
        self.add(v);
    }
    fn write_i64(&mut self, v: i64) {
        self.add(v as u64);
    }
    fn write_usize(&mut self, v: usize) {
        self.add(v as u64);
    }
}

/// A `HashMap` hashed with [`Mix`].
pub(crate) type Map<K, V> = std::collections::HashMap<K, V, BuildHasherDefault<Mix>>;
