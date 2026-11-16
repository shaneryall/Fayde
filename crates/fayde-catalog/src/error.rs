use std::{fmt, path::PathBuf};

/// Errors reading a catalogue.
#[derive(Debug)]
pub enum CatalogError {
    /// A file could not be read.
    Io {
        /// File.
        path: PathBuf,
        /// Underlying error.
        source: std::io::Error,
    },
    /// A file exists but its content is not what the catalogue format requires.
    Malformed {
        /// File.
        path: PathBuf,
        /// What was wrong.
        reason: String,
    },
    /// The directory does not look like the expected catalogue layout.
    NotFound(String),
    /// A query argument was invalid (non-finite, negative radius, ...).
    InvalidQuery(String),
}

impl fmt::Display for CatalogError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io { path, source } => write!(f, "{}: {source}", path.display()),
            Self::Malformed { path, reason } => write!(f, "{}: {reason}", path.display()),
            Self::NotFound(what) => write!(f, "catalogue not found: {what}"),
            Self::InvalidQuery(what) => write!(f, "invalid catalogue query: {what}"),
        }
    }
}

impl std::error::Error for CatalogError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io { source, .. } => Some(source),
            _ => None,
        }
    }
}
