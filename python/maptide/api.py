import os
from typing import Dict, Tuple, Optional, Any
from . import maptide  #  type: ignore


BASES = ["A", "C", "G", "T", "DS", "N"]


def query(
    bam: str,
    region: Optional[str] = None,
    bai: Optional[str] = None,
    mapping_quality: int = 0,
    base_quality: int = 0,
    supplementary: bool = False,
    secondary: bool = False,
    qc_fail: bool = False,
    duplicate: bool = False,
    annotated: bool = False,
) -> Dict[str, Dict[Tuple[int, int], Any]]:
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
    supplementary : bool, optional
        Include supplementary alignments in the pileup (default: False)
    secondary : bool, optional
        Include secondary alignments in the pileup (default: False)
    qc_fail : bool, optional
        Include reads failing quality control in the pileup (default: False)
    duplicate : bool, optional
        Include duplicate reads in the pileup (default: False)
    annotated : bool, optional
        Return frequencies annotated with their bases, as a `dict[str, int]`. Default is to return frequencies only, as a `list[int]` (default: False)

    Returns
    -------
    dict
        Mapping: reference -> (reference position, insert position) -> [base frequencies].
    """

    if region:
        if not bai and os.path.isfile(bam + ".bai"):
            bai = bam + ".bai"
        data = maptide.query(
            bam,
            bai,
            region,
            mapping_quality,
            base_quality,
            supplementary,
            secondary,
            qc_fail,
            duplicate,
        )
    else:
        data = maptide.all(
            bam,
            mapping_quality,
            base_quality,
            supplementary,
            secondary,
            qc_fail,
            duplicate,
        )

    if annotated:
        for _, positions in data.items():
            for position, frequencies in positions.items():
                positions[position] = dict(zip(BASES, frequencies))

    return data


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
