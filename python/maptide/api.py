import os
from typing import Dict, Tuple, List
from . import maptide  #  type: ignore


def query(
    bam: str,
    region: str | None = None,
    bai: str | None = None,
    mapping_quality: int = 0,
    base_quality: int = 0,
) -> Dict[Tuple[int, int], List[int]]:
    """Obtain base frequencies from the provided BAM file.

    Required arguments:
        `bam`: Path to the BAM file.
    Optional arguments:
        `region`
        `bai`
        `mapping_quality`
        `base_quality`
    Returns:
        A `dict` mapping tuples of the form `(int, int)` to base frequencies.
    """
    if region:
        if not bai and os.path.isfile(bam + ".bai"):
            bai = bam + ".bai"
        return maptide.query(bam, bai, region, mapping_quality, base_quality)
    else:
        return maptide.all(bam, mapping_quality, base_quality)


def parse_region(region: str):
    return maptide.parse_region(region)
