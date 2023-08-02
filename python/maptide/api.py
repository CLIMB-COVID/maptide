import os
from typing import Dict, Tuple, List
from . import maptide  #  type: ignore


def query(
    bam: str,
    region: str | None = None,
    bai: str | None = None,
    mapping_quality: int = 0,
    base_quality: int = 0,
) -> Dict[str, Dict[Tuple[int, int], List[int]]]:
    """Performs a pileup over a region, obtaining per-position base frequencies for the provided BAM file.

    Parameters
    ----------
    bam : str
        Path to the BAM file.
    region : str, optional
        Region to query, in the form `CHROM:START-END` (default: all positions)
    bai : str, optional
        Path to index file (default: same path as the BAM file, but with .bai appended)
    mapping_quality : int, optional
        Minimum mapping quality for a read to be included in the pileup (default: 0)
    base_quality : int, optional
        Minimum base quality for a base within a read to be included in the pileup (default: 0)

    Returns
    -------
    dict
        Mapping: reference -> (reference position, insert position) -> [base frequencies].
    """

    if region:
        if not bai and os.path.isfile(bam + ".bai"):
            bai = bam + ".bai"
        return maptide.query(bam, bai, region, mapping_quality, base_quality)
    else:
        return maptide.all(bam, mapping_quality, base_quality)


def parse_region(region: str) -> Tuple[str, int, int]:
    """Parses a region of the form `CHROM:START-END`, returning the tuple `(CHROM, START, END)`.

    Parameters
    ----------
    region : str
        Region to parse.

    Returns
    -------
    tuple
        Parsed region tuple.
    """
    return maptide.parse_region(region)
