use noodles::core::region;
use std::error::Error;
use std::fmt::{self, Display};
use std::io;

#[derive(Debug)]
pub enum BaseCountError {
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

impl From<io::Error> for BaseCountError {
    fn from(e: io::Error) -> Self {
        BaseCountError::IOError(e)
    }
}

impl From<region::ParseError> for BaseCountError {
    fn from(e: region::ParseError) -> Self {
        BaseCountError::ParseError(e)
    }
}

impl Display for BaseCountError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            BaseCountError::KeyNotFound => f.write_str("KeyNotFound"),
            BaseCountError::IndexNotFound => f.write_str("IndexNotFound"),
            BaseCountError::InvalidBase => f.write_str("InvalidBase"),
            BaseCountError::IntegerOverflow => f.write_str("IntegerOverlow"),
            BaseCountError::AlignmentStartNotFound => f.write_str("AlignmentStartNotFound"),
            BaseCountError::AlignmentEndNotFound => f.write_str("AlignmentEndNotFound"),
            BaseCountError::MappingQualityNotFound => f.write_str("MappingQualityNotFound"),
            BaseCountError::QualityScoreNotFound => f.write_str("QualityScoreNotFound"),
            BaseCountError::ReferenceSequenceIDNotFound => {
                f.write_str("ReferenceSequenceIDNotFound")
            }
            BaseCountError::IOError(ref _e) => f.write_str("IOError"),
            BaseCountError::ParseError(ref _e) => f.write_str("ParseError"),
        }
    }
}

impl Error for BaseCountError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match *self {
            BaseCountError::IOError(ref e) => Some(e),
            BaseCountError::ParseError(ref e) => Some(e),
            _ => None,
        }
    }
}
