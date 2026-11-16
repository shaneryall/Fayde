use std::fmt;

/// Errors building a WCS or reading one from a FITS header.
#[derive(Debug, Clone, PartialEq)]
pub enum WcsError {
    /// A required value was NaN or infinite.
    NonFinite(&'static str),
    /// The CD matrix is singular, so pixels cannot be recovered from sky positions.
    SingularMatrix,
    /// A required header keyword is missing.
    MissingKey(String),
    /// A header keyword has the wrong type or an unparseable value.
    InvalidValue {
        /// Keyword.
        key: String,
        /// What was wrong.
        reason: String,
    },
    /// The header describes something this crate does not support (yet).
    Unsupported(String),
}

impl fmt::Display for WcsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonFinite(what) => write!(f, "non-finite {what}"),
            Self::SingularMatrix => write!(f, "singular CD matrix"),
            Self::MissingKey(k) => write!(f, "missing header keyword {k}"),
            Self::InvalidValue { key, reason } => write!(f, "invalid {key}: {reason}"),
            Self::Unsupported(what) => write!(f, "unsupported WCS: {what}"),
        }
    }
}

impl std::error::Error for WcsError {}
