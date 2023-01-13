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
            MapTideError::IntegerOverflow => f.write_str("IntegerOverlow"),
            MapTideError::AlignmentStartNotFound => f.write_str("AlignmentStartNotFound"),
            MapTideError::AlignmentEndNotFound => f.write_str("AlignmentEndNotFound"),
            MapTideError::MappingQualityNotFound => f.write_str("MappingQualityNotFound"),
            MapTideError::QualityScoreNotFound => f.write_str("QualityScoreNotFound"),
            MapTideError::ReferenceSequenceIDNotFound => f.write_str("ReferenceSequenceIDNotFound"),
            MapTideError::IOError(ref _e) => f.write_str("IOError"),
            MapTideError::ParseError(ref _e) => f.write_str("ParseError"),
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
