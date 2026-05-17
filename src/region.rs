use crate::error::MapTideError;
use noodles::core::region::Interval;
use noodles::core::Region;
use noodles::sam::alignment::Record;
use regex::Regex;

/// Validate that the `region` is one of the formats: `CHROM`, `CHROM:START` `CHROM:START-END`.
pub fn validate_region(region: &str) -> Result<(), MapTideError> {
    let re = Regex::new(r"^[^:]+(:\d+(-\d+)?)?$").unwrap();
    if !re.is_match(region) {
        Err(MapTideError::Error(format!(
            "invalid region: '{}' is not one of the supported formats: CHROM, CHROM:START, CHROM:START-END",
            region
        )))
    } else {
        Ok(())
    }
}

/// Check the interval defined by the alignment of `record` intersects the interval defined in `region`.
pub fn intersects_region(record: &Record, region: &Region) -> Result<bool, MapTideError> {
    let seq_start = record
        .alignment_start()
        .ok_or_else(|| MapTideError::AlignmentStartNotFound)?;

    let seq_end = record
        .alignment_end()
        .ok_or_else(|| MapTideError::AlignmentEndNotFound)?;

    let seq_interval = Interval::from(seq_start..=seq_end);

    if region.interval().intersects(seq_interval) {
        Ok(true)
    } else {
        Ok(false)
    }
}

#[cfg(test)]
mod tests {
    use super::{intersects_region, validate_region};
    use noodles::core::Position;
    use noodles::core::Region;
    use noodles::sam::alignment::Record;
    use noodles::sam::record::cigar::{op::Kind, Op};
    use noodles::sam::record::Cigar;

    fn build_record(start: usize, cigar: Cigar) -> Record {
        Record::builder()
            .set_alignment_start(Position::new(start).unwrap())
            .set_cigar(cigar)
            .build()
    }

    #[test]
    fn test_validate_region_valid() {
        assert!(validate_region("chr1").is_ok());
        assert!(validate_region("chr1:100").is_ok());
        assert!(validate_region("chr1:100-200").is_ok());
    }

    #[test]
    fn test_validate_region_invalid() {
        assert!(validate_region("chr1:").is_err());
        assert!(validate_region("chr1:100-").is_err());
        assert!(validate_region("chr1:-100").is_err());
        assert!(validate_region("chr1:100-200:300").is_err());
        assert!(validate_region("chr1:100-200-300").is_err());
        assert!(validate_region("chr1:abc").is_err());
        assert!(validate_region("chr1:100-abc").is_err());
        assert!(validate_region(":100-200").is_err());
    }

    #[test]
    fn test_intersects_region_overlapping() {
        // read spans 50-150, region is 100-200
        let record = build_record(
            50,
            Cigar::try_from(vec![Op::new(Kind::Match, 101)]).unwrap(),
        );
        let region: Region = "chr1:100-200".parse().unwrap();
        assert!(intersects_region(&record, &region).unwrap());
    }

    #[test]
    fn test_intersects_region_contained() {
        // read spans 110-120, region is 100-200
        let record = build_record(
            110,
            Cigar::try_from(vec![Op::new(Kind::Match, 11)]).unwrap(),
        );
        let region: Region = "chr1:100-200".parse().unwrap();
        assert!(intersects_region(&record, &region).unwrap());
    }

    #[test]
    fn test_intersects_region_non_overlapping() {
        // read spans 1-50, region is 100-200
        let record = build_record(1, Cigar::try_from(vec![Op::new(Kind::Match, 50)]).unwrap());
        let region: Region = "chr1:100-200".parse().unwrap();
        assert!(!intersects_region(&record, &region).unwrap());
    }
}
