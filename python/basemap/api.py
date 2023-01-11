from .basemap import all_, query_, parse_region_
from typing import Optional


def all(bam_path: str, mapping_quality: int = 0, base_quality: int = 0):
    return all_(bam_path, mapping_quality, base_quality)


def query(bam_path: str, region: str, mapping_quality: int = 0, base_quality: int = 0):
    return query_(bam_path, None, region, mapping_quality, base_quality)


def iquery(
    bam_path: str,
    region: str,
    bai_path: Optional[str] = None,
    mapping_quality: int = 0,
    base_quality: int = 0,
):
    if not bai_path:
        bai_path = bam_path + ".bai"

    return query_(bam_path, bai_path, region, mapping_quality, base_quality)


def parse_region(region: str):
    return parse_region_(region)
