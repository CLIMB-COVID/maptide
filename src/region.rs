use crate::error::MapTideError;
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

#[cfg(test)]
mod tests {
    use super::validate_region;

    #[test]
    fn test_validate_region() {
        // Valid regions
        assert!(validate_region("chr1").is_ok());
        assert!(validate_region("chr1:100").is_ok());
        assert!(validate_region("chr1:100-200").is_ok());

        // Invalid regions
        assert!(validate_region("chr1:").is_err());
        assert!(validate_region("chr1:100-").is_err());
        assert!(validate_region("chr1:-100").is_err());
        assert!(validate_region("chr1:100-200:300").is_err());
        assert!(validate_region("chr1:100-200-300").is_err());
        assert!(validate_region("chr1:abc").is_err());
        assert!(validate_region("chr1:100-abc").is_err());
        assert!(validate_region(":100-200").is_err());
    }
}
