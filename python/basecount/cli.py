from .api import query, parse_region
import argparse
import math
import sys
import csv


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
        chrom, start, end = parse_region(region)
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
    index_group = parser.add_mutually_exclusive_group()
    parser.add_argument("bam", help="Path to BAM file")
    parser.add_argument(
        "--region",
        help="Region to view, specified in the form CHROM:START-END (default: everything)",
    )
    index_group.add_argument(
        "--index",
        default=None,
        help="Path to index (BAI) file (default: </path/to/bam>.bai)",
    )
    parser.add_argument(
        "--mapq",
        type=int,
        default=0,
        help="Minimum mapping quality (default: %(default)s)",
    )
    parser.add_argument(
        "--baseq",
        type=int,
        default=0,
        help="Minimum base quality (default: %(default)s)",
    )
    index_group.add_argument(
        "--noindex",
        action="store_true",
        default=False,
        help="Do not use an index file when querying the BAM file (default: %(default)s)",
    )
    parser.add_argument(
        "--stats",
        action="store_true",
        default=False,
        help="Output additional per-position statistics (default: %(default)s)",
    )
    parser.add_argument(
        "--decimals",
        type=int,
        default=3,
        help="Number of decimal places to display (default: %(default)s)",
    )

    args = parser.parse_args()

    columns = [
        "chrom",
        "pos",
        "ins",
        "cov",
        "a",
        "c",
        "g",
        "t",
        "ds",
        "n",
    ]

    if args.stats:
        columns.extend(
            [
                "pc_a",
                "pc_c",
                "pc_g",
                "pc_t",
                "pc_ds",
                "pc_n",
                "entropy",
                "secondary_entropy",
            ]
        )

    writer = csv.writer(sys.stdout, delimiter="\t")
    writer.writerow(columns)

    data = query(
        bam=args.bam,
        region=args.region,
        bai=args.index,
        mapping_quality=args.mapq,
        base_quality=args.baseq,
        indexed=not args.noindex,
    )

    for row in iterate(
        data, region=args.region, stats=args.stats, decimals=args.decimals
    ):
        writer.writerow(row)
