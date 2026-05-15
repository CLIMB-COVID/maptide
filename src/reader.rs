use crate::error::MapTideError;
use std::fs::File;

/// Open the BAM file located at `bam_path` and return a reader.
pub fn get_reader(
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
