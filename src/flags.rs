use crate::error::MapTideError;
use noodles::sam::record::Flags;

/// Get flags for filtering records.
pub fn get_filter_flags(
    keep_supplementary: bool,
    keep_secondary: bool,
    keep_qc_fail: bool,
    keep_duplicate: bool,
) -> Result<Flags, MapTideError> {
    // Unmapped flag is always set
    let mut flags = Flags::UNMAPPED;

    if !keep_supplementary {
        flags |= Flags::SUPPLEMENTARY;
    }
    if !keep_secondary {
        flags |= Flags::SECONDARY;
    }
    if !keep_qc_fail {
        flags |= Flags::QC_FAIL;
    }
    if !keep_duplicate {
        flags |= Flags::DUPLICATE;
    }
    Ok(flags)
}

#[cfg(test)]
mod tests {
    use noodles::sam::record::Flags;

    use super::get_filter_flags;

    #[test]
    fn test_get_filter_flags() {
        // Test with only filtering unmapped reads
        let flags = get_filter_flags(true, true, true, true).unwrap();
        assert_eq!(flags.bits(), Flags::UNMAPPED.bits());

        // Test with filtering all flags
        let flags = get_filter_flags(false, false, false, false).unwrap();
        assert_eq!(
            flags.bits(),
            (Flags::UNMAPPED
                | Flags::SUPPLEMENTARY
                | Flags::SECONDARY
                | Flags::QC_FAIL
                | Flags::DUPLICATE)
                .bits()
        );

        // Test with keeping supplementary and secondary reads
        let flags = get_filter_flags(true, true, false, false).unwrap();
        assert_eq!(
            flags.bits(),
            (Flags::UNMAPPED | Flags::QC_FAIL | Flags::DUPLICATE).bits()
        );

        // Test with keeping qc fail and duplicate reads
        let flags = get_filter_flags(false, false, true, true).unwrap();
        assert_eq!(
            flags.bits(),
            (Flags::UNMAPPED | Flags::SUPPLEMENTARY | Flags::SECONDARY).bits()
        );
    }
}
