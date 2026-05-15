import argparse
import math
import sys
import csv
from importlib.metadata import version
from . import api


def entropy(probabilities, normalised=False):
    ent = sum([-(x * math.log2(x)) if x != 0 else 0 for x in probabilities])

    if normalised:
        return ent / math.log2(len(probabilities))
    else:
        return ent


def get_stats(counts, decimals=3):
    coverage = sum(counts)
    probabilities = [count / coverage if coverage > 0 else 0.0 for count in counts]
    percentages = [100 * probability for probability in probabilities]
    ent = entropy(probabilities, normalised=True)
    secondary_count = list(counts)
    secondary_count.pop(counts.index(max(counts)))
    secondary_coverage = sum(secondary_count)
    secondary_probabilities = [
        count / secondary_coverage if secondary_coverage > 0 else 0.0
        for count in secondary_count
    ]
    secondary_ent = entropy(secondary_probabilities, normalised=True)
    return (
        [coverage]
        + counts
        + [round(x, decimals) for x in percentages + [ent, secondary_ent]]
    )


def iterate(data, region=None, stats=False, decimals=3):
    if region:
        chrom, start, end = api.parse_region(region)
        for (pos, ins_pos), row in sorted(data[chrom].items()):
            if (not start or pos >= start) and (not end or pos <= end):
                if stats:
                    yield [chrom, pos, ins_pos] + get_stats(row, decimals=decimals)
                else:
                    yield [chrom, pos, ins_pos, sum(row)] + row
    else:
        for chrom, chrom_data in data.items():
            for (pos, ins_pos), row in sorted(chrom_data.items()):
                if stats:
                    yield [chrom, pos, ins_pos] + get_stats(row, decimals=decimals)
                else:
                    yield [chrom, pos, ins_pos, sum(row)] + row


def run():
    parser = argparse.ArgumentParser()
    parser.add_argument("bam", help="Path to BAM file")
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=version("maptide"),
    )
    parser.add_argument(
        "-r",
        "--region",
        help="Region to view, specified in the form CHROM:START-END (default: everything)",
    )
    parser.add_argument(
        "-i",
        "--index",
        help="Path to index (BAI) file (default: </path/to/bam>.bai)",
    )
    parser.add_argument(
        "-m",
        "--mapping-quality",
        type=int,
        default=0,
        help="Minimum mapping quality (default: %(default)s)",
    )
    parser.add_argument(
        "-b",
        "--base-quality",
        type=int,
        default=0,
        help="Minimum base quality (default: %(default)s)",
    )
    parser.add_argument(
        "--supplementary",
        action="store_true",
        default=False,
        help="Include supplementary alignments (default: %(default)s)",
    )
    parser.add_argument(
        "--secondary",
        action="store_true",
        default=False,
        help="Include secondary alignments (default: %(default)s)",
    )
    parser.add_argument(
        "--qc-fail",
        action="store_true",
        default=False,
        help="Include reads failing quality control (default: %(default)s)",
    )
    parser.add_argument(
        "--duplicate",
        action="store_true",
        default=False,
        help="Include duplicate reads (default: %(default)s)",
    )
    parser.add_argument(
        "-s",
        "--stats",
        action="store_true",
        default=False,
        help="Output additional per-position statistics (default: %(default)s)",
    )
    parser.add_argument(
        "-d",
        "--decimals",
        type=int,
        default=3,
        help="Number of decimal places to display (default: %(default)s)",
    )

    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    columns = [
        "chrom",
        "pos",
        "ins",
        "cov",
    ] + [base.lower() for base in api.BASES]

    if args.stats:
        columns.extend(
            [f"pc_{base.lower()}" for base in api.BASES]
            + [
                "entropy",
                "secondary_entropy",
            ]
        )

    writer = csv.writer(sys.stdout, delimiter="\t")
    writer.writerow(columns)

    data = api.query(
        bam=args.bam,
        region=args.region,
        bai=args.index,
        mapping_quality=args.mapping_quality,
        base_quality=args.base_quality,
    )

    for row in iterate(
        data, region=args.region, stats=args.stats, decimals=args.decimals
    ):
        writer.writerow(row)
