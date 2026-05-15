mod error;
mod flags;
mod pileup;
mod reader;
mod region;
mod types;

use error::MapTideError;
use flags::get_filter_flags;
use noodles::bam::bai;
use noodles::core::region::ParseError;
use noodles::core::Region;
use pileup::{count_record, min_mapping_quality};
use pyo3::exceptions::{PyException, PyIOError, PyIndexError, PyKeyError, PyOverflowError};
use pyo3::prelude::*;
use reader::get_reader;
use region::{intersects_region, validate_region};
use types::{CoordinateMap, MapTide, RefLengths, RefMap};

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

fn init_maps() -> (RefMap, MapTide, RefLengths) {
    // Map of reference names to vector for storing base counts
    let ref_arrs: RefMap = RefMap::new();

    // Map of reference names to CoordinateMap, for storing insertion data
    let ins_maps: MapTide = MapTide::new();

    // Map of reference names to reference lengths
    let ref_lengths: RefLengths = RefLengths::new();

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
    ref_lengths: &RefLengths,
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
    ref_lengths: &RefLengths,
) -> Result<MapTide, MapTideError> {
    for (ref_name, _) in ref_lengths.iter() {
        let (ref_arr, offset) = ref_arrs
            .get(ref_name)
            .ok_or_else(|| MapTideError::KeyNotFound)?;

        let ins_map = ins_maps
            .get_mut(ref_name)
            .ok_or_else(|| MapTideError::KeyNotFound)?;

        for (i, row) in ref_arr.iter().enumerate() {
            ins_map.entry((i + offset + 1, 0)).or_insert(*row);
        }
    }
    Ok(ins_maps)
}

#[pyfunction(signature = (bam_path, mapping_quality, base_quality, supplementary, secondary, qc_fail, duplicate))]
fn all(
    bam_path: String,
    mapping_quality: usize,
    base_quality: usize,
    supplementary: bool,
    secondary: bool,
    qc_fail: bool,
    duplicate: bool,
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

    // Initialise coordinates
    init_coordinates(&mut ref_arrs, &mut ins_maps, &ref_lengths, None)?;

    // Define flags for filtering records
    let flags = get_filter_flags(supplementary, secondary, qc_fail, duplicate)?;

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

#[pyfunction(signature = (bam_path, bai_path, region, mapping_quality, base_quality, supplementary, secondary, qc_fail, duplicate))]
fn query(
    bam_path: String,
    bai_path: Option<String>,
    region: String,
    mapping_quality: usize,
    base_quality: usize,
    supplementary: bool,
    secondary: bool,
    qc_fail: bool,
    duplicate: bool,
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
    validate_region(&region)?;
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
    let flags = get_filter_flags(supplementary, secondary, qc_fail, duplicate)?;

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
                || !intersects_region(&record, &region)?
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

#[pyfunction(signature = (region))]
fn parse_region(region: String) -> PyResult<(String, Option<usize>, Option<usize>)> {
    validate_region(&region)?;
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
fn maptide(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(all, m)?)?;
    m.add_function(wrap_pyfunction!(query, m)?)?;
    m.add_function(wrap_pyfunction!(parse_region, m)?)?;

    Ok(())
}
