use std::{fmt, path::PathBuf};

#[derive(Debug)]
pub enum CatalogError {
    Io {
        path: PathBuf,
        source: std::io::Error,
    },
    Malformed {
        path: PathBuf,
        reason: String,
    },
    NotFound(String),
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
