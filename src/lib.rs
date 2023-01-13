use noodles::bam::bai;
use noodles::core::region::{Interval, ParseError};
use noodles::core::{Position, Region};
use noodles::sam::alignment::Record;
use noodles::sam::record::cigar::op::Kind;
use noodles::sam::record::sequence::{Base, Sequence};
use noodles::sam::record::{Flags, QualityScores};
use pyo3::exceptions::{PyException, PyIOError, PyIndexError, PyKeyError, PyOverflowError};
use pyo3::prelude::*;
use std::collections::HashMap;
use std::fs::File;

mod error;
use error::MapTideError;

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Clone, Copy)]
struct Coordinate(usize, usize);

impl IntoPy<PyObject> for Coordinate {
    fn into_py(self, py: Python<'_>) -> PyObject {
        (self.0, self.1).into_py(py)
    }
}

impl From<MapTideError> for PyErr {
    fn from(e: MapTideError) -> Self {
        match e {
            MapTideError::KeyNotFound => PyKeyError::new_err(e.to_string()),
            MapTideError::IndexNotFound => PyIndexError::new_err(e.to_string()),
            MapTideError::IntegerOverflow => PyOverflowError::new_err(e.to_string()),
            MapTideError::IOError(e) => PyIOError::new_err(e.to_string()),
            _ => PyException::new_err(e.to_string()),
        }
    }
}

type CoordinateMap = HashMap<Coordinate, [usize; 6]>;

type RefArr = Vec<[usize; 6]>;

type RefMap = HashMap<String, (RefArr, usize)>;

type MapTide = HashMap<String, CoordinateMap>;

type RefLengths = HashMap<String, usize>;

/// Open the BAM file located at `bam_path` and return a reader.
fn get_reader(
    bam_path: String,
) -> Result<noodles::bam::Reader<noodles::bgzf::Reader<File>>, MapTideError> {
    // Open file
    let file = File::open(bam_path)?;

    // Create a reader from the file
    let mut reader = noodles::bam::Reader::new(file);

    // Read the SAM header
    reader.read_header()?;

    // Return the reader
    Ok(reader)
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
            ins_map
                .entry(Coordinate(ref_pos, ins_pos))
                .or_insert_with(|| [0; 6])[0] += 1;
            Ok(())
        }
        Some(&Base::C) => {
            ins_map
                .entry(Coordinate(ref_pos, ins_pos))
                .or_insert_with(|| [0; 6])[1] += 1;
            Ok(())
        }
        Some(&Base::G) => {
            ins_map
                .entry(Coordinate(ref_pos, ins_pos))
                .or_insert_with(|| [0; 6])[2] += 1;
            Ok(())
        }
        Some(&Base::T) => {
            ins_map
                .entry(Coordinate(ref_pos, ins_pos))
                .or_insert_with(|| [0; 6])[3] += 1;
            Ok(())
        }
        Some(&Base::N) => {
            ins_map
                .entry(Coordinate(ref_pos, ins_pos))
                .or_insert_with(|| [0; 6])[5] += 1;
            Ok(())
        }
        Some(_) => Err(MapTideError::InvalidBase),
        None => Err(MapTideError::KeyNotFound),
    }
}

/// Use the CIGAR information of `record` to count each base in its sequence, and add them to `ref_arr`, or `ins_map`.
///
/// Bases are ignored if their quality score is less than `base_quality`.
fn count_record(
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
                        count_ins_base(ins_map, seq, ref_pos, seq_pos, i)?;
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

/// Check the interval defined by the alignment of `record` intersects the interval defined in `region`.
fn intersects(record: &Record, region: &Region) -> Result<bool, MapTideError> {
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
fn min_mapping_quality(record: &Record, mapping_quality: usize) -> Result<bool, MapTideError> {
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

fn init_maps() -> (RefMap, MapTide, RefLengths) {
    // Map of reference names to vector for storing base counts
    let ref_arrs: RefMap = HashMap::new();

    // Map of reference names to CoordinateMap, for storing insertion data
    let ins_maps: MapTide = MapTide::new();

    // Map of reference names to reference lengths
    let ref_lengths: RefLengths = HashMap::new();

    (ref_arrs, ins_maps, ref_lengths)
}

/// Initialise arrays in `ref_arrs`, and CoordinateMaps in `ins_maps`.
///
/// If `region` is `None`, initialises arrays for all positions across all references.
///
/// If `region` is `Some`, initialises array over the region specified.
fn init_coordinates(
    ref_arrs: &mut RefMap,
    ins_maps: &mut MapTide,
    ref_lengths: &HashMap<String, usize>,
    region: Option<&Region>,
) -> Result<(), MapTideError> {
    if let Some(reg) = region {
        let region_name = reg.name();
        let interval = reg.interval();

        // Get length of the region name's sequence
        let ref_length = ref_lengths
            .get(region_name)
            .ok_or_else(|| MapTideError::KeyNotFound)?;

        // Handle unbounded region start
        let region_start = match interval.start() {
            Some(x) => x.get(),
            None => 1,
        };

        // Handle unbounded region end
        let region_end = match interval.end() {
            Some(x) => x.get(),
            None => *ref_length,
        };

        // Add reference to ref_arrs and ins_maps
        ref_arrs.entry(region_name.to_owned()).or_insert_with(|| {
            (
                vec![[0; 6]; region_end - region_start + 1],
                region_start - 1,
            )
        });

        ins_maps
            .entry(region_name.to_owned())
            .or_insert_with(|| CoordinateMap::new());
    } else {
        // Add every reference to ref_arrs and ins_maps
        for (ref_name, ref_length) in ref_lengths.iter() {
            ref_arrs
                .entry(ref_name.to_owned())
                .or_insert_with(|| (vec![[0; 6]; *ref_length], 0));

            ins_maps
                .entry(ref_name.to_owned())
                .or_insert_with(|| CoordinateMap::new());
        }
    }

    Ok(())
}

/// Merge `ref_arrs` into `ins_maps` to have a single `MapTide` containing all coordinates and counts.
fn merge_into_base_map(
    ref_arrs: &RefMap,
    mut ins_maps: MapTide,
    ref_lengths: &HashMap<String, usize>,
) -> Result<MapTide, MapTideError> {
    for (ref_name, _) in ref_lengths.iter() {
        let (ref_arr, offset) = ref_arrs
            .get(ref_name)
            .ok_or_else(|| MapTideError::KeyNotFound)?;

        let ins_map = ins_maps
            .get_mut(ref_name)
            .ok_or_else(|| MapTideError::KeyNotFound)?;

        for (i, row) in ref_arr.iter().enumerate() {
            ins_map.entry(Coordinate(i + offset + 1, 0)).or_insert(*row);
        }
    }
    Ok(ins_maps)
}

#[pyfunction]
fn all_(bam_path: String, mapping_quality: usize, base_quality: usize) -> PyResult<MapTide> {
    // Create initial maps
    let (mut ref_arrs, mut ins_maps, mut ref_lengths) = init_maps();

    // Reader for iterating through records
    let mut reader = get_reader(bam_path)?;

    // Reference sequence information
    let ref_seqs = reader.read_reference_sequences()?;

    // Add reference sequence information to HashMaps
    for reff in ref_seqs.iter() {
        ref_lengths.insert(reff.0.to_owned(), reff.1.length().get());
    }

    // Initialise coordinates
    init_coordinates(&mut ref_arrs, &mut ins_maps, &ref_lengths, None)?;

    // Define flags for filtering records
    let flags = Flags::from(
        Flags::UNMAPPED.bits()
            + Flags::SUPPLEMENTARY.bits()
            + Flags::SECONDARY.bits()
            + Flags::QC_FAIL.bits()
            + Flags::DUPLICATE.bits(),
    );

    for result in reader.records() {
        let record = result?;

        if record.flags().intersects(flags) || !min_mapping_quality(&record, mapping_quality)? {
            continue;
        }

        let ref_seq_id = record
            .reference_sequence_id()
            .ok_or_else(|| MapTideError::ReferenceSequenceIDNotFound)?;

        let ref_name = ref_seqs
            .get_index(ref_seq_id)
            .ok_or_else(|| MapTideError::KeyNotFound)?
            .0;

        let ref_length = ref_lengths
            .get(ref_name)
            .ok_or_else(|| MapTideError::KeyNotFound)?;

        let (ref_arr, offset) = ref_arrs
            .get_mut(ref_name)
            .ok_or_else(|| MapTideError::KeyNotFound)?;

        let ins_map = ins_maps
            .get_mut(ref_name)
            .ok_or_else(|| MapTideError::KeyNotFound)?;

        count_record(
            ref_arr,
            *offset,
            ins_map,
            &record,
            base_quality,
            1,
            *ref_length,
        )?;
    }

    let base_map = merge_into_base_map(&ref_arrs, ins_maps, &ref_lengths)?;
    Ok(base_map)
}

#[pyfunction]
fn query_(
    bam_path: String,
    bai_path: Option<String>,
    region: String,
    mapping_quality: usize,
    base_quality: usize,
) -> PyResult<MapTide> {
    // Create initial maps
    let (mut ref_arrs, mut ins_maps, mut ref_lengths) = init_maps();

    // Reader for iterating through records
    let mut reader = get_reader(bam_path)?;

    // Reference sequence information
    let ref_seqs = reader.read_reference_sequences()?;

    // Add reference sequence information to HashMaps
    for reff in ref_seqs.iter() {
        ref_lengths.insert(reff.0.to_owned(), reff.1.length().get());
    }

    // Parse region
    let region: Region = region
        .parse()
        .map_err(|x: ParseError| PyException::new_err(x.to_string()))?;
    let region_name = region.name();

    // Handle unbounded region start
    let region_start = match region.interval().start() {
        Some(x) => x.get(),
        None => 1,
    };

    // Handle unbounded region end
    let region_end = match region.interval().end() {
        Some(x) => x.get(),
        None => *ref_lengths
            .get(region_name)
            .ok_or_else(|| MapTideError::KeyNotFound)?,
    };

    // Initialise coordinates
    init_coordinates(&mut ref_arrs, &mut ins_maps, &ref_lengths, Some(&region))?;

    // Define flags for filtering records
    let flags = Flags::from(
        Flags::UNMAPPED.bits()
            + Flags::SUPPLEMENTARY.bits()
            + Flags::SECONDARY.bits()
            + Flags::QC_FAIL.bits()
            + Flags::DUPLICATE.bits(),
    );

    let (ref_arr, offset) = ref_arrs
        .get_mut(region_name)
        .ok_or_else(|| MapTideError::KeyNotFound)?;

    let ins_map = ins_maps
        .get_mut(region_name)
        .ok_or_else(|| MapTideError::KeyNotFound)?;

    if let Some(b_path) = bai_path {
        // Read the index file
        let index = bai::read(b_path)?;

        // Create query iterator over reads intersecting the region
        let query = reader.query(&ref_seqs, &index, &region)?;

        for result in query {
            let record = result?;
            if record.flags().intersects(flags) || !min_mapping_quality(&record, mapping_quality)? {
                continue;
            }

            count_record(
                ref_arr,
                *offset,
                ins_map,
                &record,
                base_quality,
                region_start,
                region_end,
            )?;
        }
    } else {
        for result in reader.records() {
            let record = result?;
            let record_ref_name = ref_seqs
                .get_index(
                    record
                        .reference_sequence_id()
                        .ok_or_else(|| MapTideError::ReferenceSequenceIDNotFound)?,
                )
                .ok_or_else(|| MapTideError::IndexNotFound)?
                .0;

            if record.flags().intersects(flags)
                || record_ref_name != region.name()
                || !intersects(&record, &region)?
                || !min_mapping_quality(&record, mapping_quality)?
            {
                continue;
            }

            count_record(
                ref_arr,
                *offset,
                ins_map,
                &record,
                base_quality,
                region_start,
                region_end,
            )?;
        }
    }

    let base_map = merge_into_base_map(&ref_arrs, ins_maps, &ref_lengths)?;
    Ok(base_map)
}

#[pyfunction]
fn parse_region_(region: String) -> PyResult<(String, Option<usize>, Option<usize>)> {
    let region: Region = region
        .parse()
        .map_err(|x: ParseError| PyException::new_err(x.to_string()))?;
    let interval = region.interval();
    let start = match interval.start() {
        Some(x) => Some(x.get()),
        None => None,
    };
    let end = match interval.end() {
        Some(x) => Some(x.get()),
        None => None,
    };

    Ok((region.name().to_string(), start, end))
}

/// A Python module implemented in Rust.
#[pymodule]
fn maptide(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(all_, m)?)?;
    m.add_function(wrap_pyfunction!(query_, m)?)?;
    m.add_function(wrap_pyfunction!(parse_region_, m)?)?;

    Ok(())
}
