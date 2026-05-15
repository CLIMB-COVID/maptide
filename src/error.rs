use noodles::core::region;
use std::error::Error;
use std::fmt::{self, Display};
use std::io;

#[derive(Debug)]
pub enum MapTideError {
    KeyNotFound,
    IndexNotFound,
    InvalidBase,
    IntegerOverflow,
    AlignmentStartNotFound,
    AlignmentEndNotFound,
    MappingQualityNotFound,
    QualityScoreNotFound,
    ReferenceSequenceIDNotFound,
    IOError(io::Error),
    ParseError(region::ParseError),
    Error(String),
}

impl From<io::Error> for MapTideError {
    fn from(e: io::Error) -> Self {
        MapTideError::IOError(e)
    }
}

impl From<region::ParseError> for MapTideError {
    fn from(e: region::ParseError) -> Self {
        MapTideError::ParseError(e)
    }
}

impl Display for MapTideError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            MapTideError::KeyNotFound => f.write_str("KeyNotFound"),
            MapTideError::IndexNotFound => f.write_str("IndexNotFound"),
            MapTideError::InvalidBase => f.write_str("InvalidBase"),
            MapTideError::IntegerOverflow => f.write_str("IntegerOverflow"),
            MapTideError::AlignmentStartNotFound => f.write_str("AlignmentStartNotFound"),
            MapTideError::AlignmentEndNotFound => f.write_str("AlignmentEndNotFound"),
            MapTideError::MappingQualityNotFound => f.write_str("MappingQualityNotFound"),
            MapTideError::QualityScoreNotFound => f.write_str("QualityScoreNotFound"),
            MapTideError::ReferenceSequenceIDNotFound => f.write_str("ReferenceSequenceIDNotFound"),
            MapTideError::IOError(ref _e) => f.write_str("IOError"),
            MapTideError::ParseError(ref _e) => f.write_str("ParseError"),
            MapTideError::Error(ref msg) => f.write_str(msg),
        }
    }
}

impl Error for MapTideError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match *self {
            MapTideError::IOError(ref e) => Some(e),
            MapTideError::ParseError(ref e) => Some(e),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::error::Error;

    #[test]
    fn test_error_variant_display() {
        assert_eq!(
            format!("{}", MapTideError::Error("something went wrong".to_string())),
            "something went wrong"
        );
    }

    #[test]
    fn test_from_io_error() {
        let io_err = io::Error::new(io::ErrorKind::NotFound, "file not found");
        let e = MapTideError::from(io_err);
        assert!(matches!(e, MapTideError::IOError(_)));
    }

    #[test]
    fn test_from_parse_error() {
        let parse_err = region::ParseError::Invalid;
        let e = MapTideError::from(parse_err);
        assert!(matches!(e, MapTideError::ParseError(_)));
    }

    #[test]
    fn test_source_wrapped_variants() {
        let e = MapTideError::IOError(io::Error::new(io::ErrorKind::Other, "inner"));
        assert!(e.source().is_some());

        let e = MapTideError::ParseError(region::ParseError::Invalid);
        assert!(e.source().is_some());
    }

    #[test]
    fn test_source_unit_variants() {
        assert!(MapTideError::KeyNotFound.source().is_none());
        assert!(MapTideError::InvalidBase.source().is_none());
        assert!(MapTideError::IntegerOverflow.source().is_none());
    }
}
