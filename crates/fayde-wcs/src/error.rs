use std::fmt;

#[derive(Debug, Clone, PartialEq)]
pub enum WcsError {
    NonFinite(&'static str),
    SingularMatrix,
    MissingKey(String),
    InvalidValue { key: String, reason: String },
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
