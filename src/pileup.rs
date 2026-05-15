use crate::error::MapTideError;
use crate::types::{CoordinateMap, RefArr};
use noodles::core::Position;
use noodles::sam::alignment::Record;
use noodles::sam::record::cigar::op::Kind;
use noodles::sam::record::sequence::{Base, Sequence};
use noodles::sam::record::QualityScores;

/// Check the quality score for the base at `seq_pos` is greater than or equal to `base_quality`.
fn min_base_quality(
    quals: &QualityScores,
    seq_pos: Position,
    base_quality: usize,
) -> Result<bool, MapTideError> {
    let base_qual = usize::from(
        quals
            .get(seq_pos)
            .ok_or_else(|| MapTideError::QualityScoreNotFound)?
            .get(),
    );

    if base_qual >= base_quality {
        Ok(true)
    } else {
        Ok(false)
    }
}

/// Check the mapping score for `record` is greater than or equal to `mapping_quality`.
pub fn min_mapping_quality(record: &Record, mapping_quality: usize) -> Result<bool, MapTideError> {
    let map_qual = usize::from(
        record
            .mapping_quality()
            .ok_or_else(|| MapTideError::MappingQualityNotFound)?
            .get(),
    );

    if map_qual >= mapping_quality {
        Ok(true)
    } else {
        Ok(false)
    }
}

/// Add the base from `seq` at `seq_pos` to `ref_arr`.
fn count_ref_base(
    ref_arr: &mut RefArr,
    seq: &Sequence,
    ref_pos: usize,
    offset: usize,
    seq_pos: Position,
) -> Result<(), MapTideError> {
    // Match the base at the given seq_pos, and update the CoordinateMap
    match seq.get(seq_pos) {
        Some(&Base::A) => {
            ref_arr[ref_pos - offset - 1][0] += 1;
            Ok(())
        }
        Some(&Base::C) => {
            ref_arr[ref_pos - offset - 1][1] += 1;
            Ok(())
        }
        Some(&Base::G) => {
            ref_arr[ref_pos - offset - 1][2] += 1;
            Ok(())
        }
        Some(&Base::T) => {
            ref_arr[ref_pos - offset - 1][3] += 1;
            Ok(())
        }
        Some(&Base::N) => {
            ref_arr[ref_pos - offset - 1][5] += 1;
            Ok(())
        }
        Some(_) => Err(MapTideError::InvalidBase),
        None => Err(MapTideError::KeyNotFound),
    }
}

/// Add the base from `seq` at `(seq_pos, ins_pos)` to `ins_map`.
fn count_ins_base(
    ins_map: &mut CoordinateMap,
    seq: &Sequence,
    ref_pos: usize,
    seq_pos: Position,
    ins_pos: usize,
) -> Result<(), MapTideError> {
    // Match the base at the given seq_pos, and update the CoordinateMap
    match seq.get(seq_pos) {
        Some(&Base::A) => {
            ins_map.entry((ref_pos, ins_pos)).or_insert_with(|| [0; 6])[0] += 1;
            Ok(())
        }
        Some(&Base::C) => {
            ins_map.entry((ref_pos, ins_pos)).or_insert_with(|| [0; 6])[1] += 1;
            Ok(())
        }
        Some(&Base::G) => {
            ins_map.entry((ref_pos, ins_pos)).or_insert_with(|| [0; 6])[2] += 1;
            Ok(())
        }
        Some(&Base::T) => {
            ins_map.entry((ref_pos, ins_pos)).or_insert_with(|| [0; 6])[3] += 1;
            Ok(())
        }
        Some(&Base::N) => {
            ins_map.entry((ref_pos, ins_pos)).or_insert_with(|| [0; 6])[5] += 1;
            Ok(())
        }
        Some(_) => Err(MapTideError::InvalidBase),
        None => Err(MapTideError::KeyNotFound),
    }
}

/// Use the CIGAR information of `record` to count each base in its sequence, and add them to `ref_arr`, or `ins_map`.
///
/// Bases are ignored if their quality score is less than `base_quality`.
pub fn count_record(
    ref_arr: &mut RefArr,
    offset: usize,
    ins_map: &mut CoordinateMap,
    record: &Record,
    base_quality: usize,
    region_start: usize,
    region_end: usize,
) -> Result<(), MapTideError> {
    // Positions are 1-based
    // This is the start position of the read in the reference
    let mut ref_pos = record
        .alignment_start()
        .ok_or_else(|| MapTideError::AlignmentStartNotFound)?
        .get();

    // This is the position locally along the sequence (minimum is 1)
    let mut seq_pos = Position::MIN;

    // The read sequence
    let seq = record.sequence();

    // The read sequence quality scores
    let quals = record.quality_scores();

    // Iterate through CIGAR information
    for cig in record.cigar().iter() {
        match cig.kind() {
            // Match/mismatch consumes both the reference and sequence
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                for _ in 1..=cig.len() {
                    if ref_pos >= region_start
                        && ref_pos <= region_end
                        && min_base_quality(quals, seq_pos, base_quality)?
                    {
                        count_ref_base(ref_arr, seq, ref_pos, offset, seq_pos)?;
                    }

                    ref_pos += 1;
                    seq_pos = seq_pos
                        .checked_add(1)
                        .ok_or_else(|| MapTideError::IntegerOverflow)?;
                }
            }

            // Insertion consumes the sequence only
            Kind::Insertion => {
                for i in 1..=cig.len() {
                    if ref_pos >= region_start
                        && ref_pos <= region_end
                        && min_base_quality(quals, seq_pos, base_quality)?
                    {
                        // ref_pos was already incremented for the match/mismatch before the insertion
                        // so here we use ref_pos - 1
                        count_ins_base(ins_map, seq, ref_pos - 1, seq_pos, i)?;
                    }

                    seq_pos = seq_pos
                        .checked_add(1)
                        .ok_or_else(|| MapTideError::IntegerOverflow)?;
                }
            }

            // Deletion/skip consumes the reference only
            Kind::Deletion | Kind::Skip => {
                for _ in 1..=cig.len() {
                    if ref_pos >= region_start && ref_pos <= region_end {
                        ref_arr[ref_pos - offset - 1][4] += 1;
                    }

                    ref_pos += 1;
                }
            }

            // Softclip consumes the sequence only
            Kind::SoftClip => {
                seq_pos = seq_pos
                    .checked_add(cig.len())
                    .ok_or_else(|| MapTideError::IntegerOverflow)?;
            }

            // Hardclip and padding don't consume the reference or the sequence
            Kind::HardClip | Kind::Pad => {}
        };
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{count_record, CoordinateMap};
    use noodles::core::Position;
    use noodles::sam::alignment::Record;
    use noodles::sam::record::cigar::{op::Kind, Op};
    use noodles::sam::record::sequence::Base;
    use noodles::sam::record::{Cigar, MappingQuality, QualityScores, Sequence};

    fn build_record(
        start: usize,
        bases: Vec<Base>,
        cigar: Cigar,
        mapq: u8,
        baseqs: Vec<u8>,
    ) -> Record {
        Record::builder()
            .set_alignment_start(Position::new(start).unwrap())
            .set_sequence(Sequence::from(bases))
            .set_cigar(cigar)
            .set_mapping_quality(MappingQuality::new(mapq).unwrap())
            .set_quality_scores(QualityScores::try_from(baseqs).unwrap())
            .build()
    }

    #[test]
    fn test_count_record_matches() {
        let bases = vec![Base::A, Base::C, Base::G, Base::T, Base::N];
        let num_bases = bases.len();
        let record = build_record(
            1,
            bases,
            Cigar::try_from(vec![
                Op::new(Kind::Match, 1),
                Op::new(Kind::SequenceMatch, 1),
                Op::new(Kind::SequenceMismatch, 1),
                Op::new(Kind::Match, 2),
            ])
            .unwrap(),
            60,
            vec![40; num_bases],
        );
        let mut ref_arr = vec![[0usize; 6]; num_bases];
        let mut ins_map = CoordinateMap::new();
        count_record(&mut ref_arr, 0, &mut ins_map, &record, 0, 1, num_bases).unwrap();
        assert_eq!(ref_arr[0], [1, 0, 0, 0, 0, 0]); // A
        assert_eq!(ref_arr[1], [0, 1, 0, 0, 0, 0]); // C
        assert_eq!(ref_arr[2], [0, 0, 1, 0, 0, 0]); // G
        assert_eq!(ref_arr[3], [0, 0, 0, 1, 0, 0]); // T
        assert_eq!(ref_arr[4], [0, 0, 0, 0, 0, 1]); // N
    }

    #[test]
    fn test_count_record_insertions() {
        let bases = vec![Base::A, Base::G, Base::T];
        let num_bases = bases.len();
        let record = build_record(
            1,
            bases,
            Cigar::try_from(vec![
                Op::new(Kind::Match, 1),
                Op::new(Kind::Insertion, 1),
                Op::new(Kind::Match, 1),
            ])
            .unwrap(),
            60,
            vec![40; num_bases],
        );
        let mut ref_arr = vec![[0usize; 6]; num_bases];
        let mut ins_map = CoordinateMap::new();
        count_record(&mut ref_arr, 0, &mut ins_map, &record, 0, 1, num_bases).unwrap();
        assert_eq!(ref_arr[0], [1, 0, 0, 0, 0, 0]); // A at pos 1
        assert_eq!(ins_map[&(1, 1)], [0, 0, 1, 0, 0, 0]); // G ins at (pos 1, ins 1)
        assert_eq!(ref_arr[1], [0, 0, 0, 1, 0, 0]); // T at pos 2
    }

    #[test]
    fn test_count_record_deletions() {
        let bases = vec![Base::A, Base::C];
        let num_bases = bases.len();
        let num_dels = 2;
        let record = build_record(
            1,
            bases,
            Cigar::try_from(vec![
                Op::new(Kind::Match, 1),
                Op::new(Kind::Deletion, 2),
                Op::new(Kind::Match, 1),
            ])
            .unwrap(),
            60,
            vec![40; num_bases],
        );
        let mut ref_arr = vec![[0usize; 6]; num_bases + num_dels];
        let mut ins_map = CoordinateMap::new();
        count_record(
            &mut ref_arr,
            0,
            &mut ins_map,
            &record,
            0,
            1,
            num_bases + num_dels,
        )
        .unwrap();
        assert_eq!(ref_arr[0], [1, 0, 0, 0, 0, 0]); // A at pos 1
        assert_eq!(ref_arr[1], [0, 0, 0, 0, 1, 0]); // del at pos 2
        assert_eq!(ref_arr[2], [0, 0, 0, 0, 1, 0]); // del at pos 3
        assert_eq!(ref_arr[3], [0, 1, 0, 0, 0, 0]); // C at pos 4
    }

    #[test]
    fn test_count_record_min_base_quality() {
        let bases = vec![Base::A, Base::C];
        let num_bases = bases.len();
        let record = build_record(
            1,
            bases,
            Cigar::try_from(vec![Op::new(Kind::Match, num_bases)]).unwrap(),
            60,
            vec![40, 5],
        );
        let mut ref_arr = vec![[0usize; 6]; num_bases];
        let mut ins_map = CoordinateMap::new();
        count_record(&mut ref_arr, 0, &mut ins_map, &record, 20, 1, num_bases).unwrap();
        assert_eq!(ref_arr[0], [1, 0, 0, 0, 0, 0]); // A passes
        assert_eq!(ref_arr[1], [0, 0, 0, 0, 0, 0]); // C excluded
    }

    #[test]
    fn test_count_record_soft_clips() {
        let record = build_record(
            1,
            vec![Base::T, Base::A, Base::C, Base::T],
            Cigar::try_from(vec![
                Op::new(Kind::SoftClip, 1),
                Op::new(Kind::Match, 2),
                Op::new(Kind::SoftClip, 1),
            ])
            .unwrap(),
            60,
            vec![40, 40, 40, 40],
        );
        let mut ref_arr = vec![[0usize; 6]; 2];
        let mut ins_map = CoordinateMap::new();
        count_record(&mut ref_arr, 0, &mut ins_map, &record, 0, 1, 2).unwrap();
        assert_eq!(ref_arr[0], [1, 0, 0, 0, 0, 0]); // A at pos 1
        assert_eq!(ref_arr[1], [0, 1, 0, 0, 0, 0]); // C at pos 2
    }

    #[test]
    fn test_count_record_hard_clips() {
        let record = build_record(
            1,
            vec![Base::T, Base::A, Base::C],
            Cigar::try_from(vec![Op::new(Kind::HardClip, 1), Op::new(Kind::Match, 3)]).unwrap(),
            60,
            vec![40, 40, 40],
        );
        let mut ref_arr = vec![[0usize; 6]; 3];
        let mut ins_map = CoordinateMap::new();
        count_record(&mut ref_arr, 0, &mut ins_map, &record, 0, 1, 3).unwrap();
        assert_eq!(ref_arr[0], [0, 0, 0, 1, 0, 0]); // T at pos 1
        assert_eq!(ref_arr[1], [1, 0, 0, 0, 0, 0]); // A at pos 2
        assert_eq!(ref_arr[2], [0, 1, 0, 0, 0, 0]); // C at pos 3
    }

    #[test]
    fn test_count_record_region_bounds_respected() {
        let record = build_record(
            1,
            vec![Base::A, Base::C, Base::G, Base::T],
            Cigar::try_from(vec![Op::new(Kind::Match, 4)]).unwrap(),
            60,
            vec![40, 40, 40, 40],
        );
        let mut ref_arr = vec![[0usize; 6]; 4];
        let mut ins_map = CoordinateMap::new();
        count_record(&mut ref_arr, 0, &mut ins_map, &record, 0, 2, 3).unwrap();
        assert_eq!(ref_arr[0], [0, 0, 0, 0, 0, 0]); // pos 1 outside region
        assert_eq!(ref_arr[1], [0, 1, 0, 0, 0, 0]); // C at pos 2
        assert_eq!(ref_arr[2], [0, 0, 1, 0, 0, 0]); // G at pos 3
        assert_eq!(ref_arr[3], [0, 0, 0, 0, 0, 0]); // pos 4 outside region
    }
}
