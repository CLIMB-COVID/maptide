from .maptide import all_, query_, parse_region_
from typing import Optional


def query(
    bam: str,
    region: Optional[str] = None,
    bai: Optional[str] = None,
    mapping_quality: int = 0,
    base_quality: int = 0,
    indexed=True,
):
    """Obtain base frequencies from the provided BAM file.

    Required arguments:
        `bam`: Path to the BAM file.
    Optional arguments:
        `region`
        `bai`
        `mapping_quality`
        `base_quality`
        `indexed`
    Returns:
        A `dict` mapping tuples of the form `(int, int)` to base frequencies.     
    """

    if not indexed and bai:
        raise Exception("Cannot set indexed=False while also providing a BAI file")

    if region:
        if indexed:
            if not bai:
                bai = bam + ".bai"
            return query_(bam, bai, region, mapping_quality, base_quality)
        else:
            return query_(bam, None, region, mapping_quality, base_quality)
    else:
        return all_(bam, mapping_quality, base_quality)


def parse_region(region: str):
    return parse_region_(region)
