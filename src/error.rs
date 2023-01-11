use noodles::core::region;
use std::error::Error;
use std::fmt::{self, Display};
use std::io;

#[derive(Debug)]
pub enum BaseMapError {
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

impl From<io::Error> for BaseMapError {
    fn from(e: io::Error) -> Self {
        BaseMapError::IOError(e)
    }
}

impl From<region::ParseError> for BaseMapError {
    fn from(e: region::ParseError) -> Self {
        BaseMapError::ParseError(e)
    }
}

impl Display for BaseMapError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            BaseMapError::KeyNotFound => f.write_str("KeyNotFound"),
            BaseMapError::IndexNotFound => f.write_str("IndexNotFound"),
            BaseMapError::InvalidBase => f.write_str("InvalidBase"),
            BaseMapError::IntegerOverflow => f.write_str("IntegerOverlow"),
            BaseMapError::AlignmentStartNotFound => f.write_str("AlignmentStartNotFound"),
            BaseMapError::AlignmentEndNotFound => f.write_str("AlignmentEndNotFound"),
            BaseMapError::MappingQualityNotFound => f.write_str("MappingQualityNotFound"),
            BaseMapError::QualityScoreNotFound => f.write_str("QualityScoreNotFound"),
            BaseMapError::ReferenceSequenceIDNotFound => f.write_str("ReferenceSequenceIDNotFound"),
            BaseMapError::IOError(ref _e) => f.write_str("IOError"),
            BaseMapError::ParseError(ref _e) => f.write_str("ParseError"),
        }
    }
}

impl Error for BaseMapError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match *self {
            BaseMapError::IOError(ref e) => Some(e),
            BaseMapError::ParseError(ref e) => Some(e),
            _ => None,
        }
    }
}
