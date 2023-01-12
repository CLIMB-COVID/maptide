from .basemap import all_, query_, parse_region_
from typing import Optional


def query(
    bam: str,
    region: Optional[str] = None,
    bai: Optional[str] = None,
    mapping_quality: int = 0,
    base_quality: int = 0,
    indexed=True,
):
    if not region:
        return all_(bam, mapping_quality, base_quality)

    if not indexed:
        return query_(bam, None, region, mapping_quality, base_quality)

    if not bai:
        bai = bam + ".bai"

    return query_(bam, bai, region, mapping_quality, base_quality)


def parse_region(region: str):
    return parse_region_(region)
