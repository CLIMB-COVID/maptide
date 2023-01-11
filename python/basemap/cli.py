from .api import all, query, iquery, parse_region
import argparse
import math
import sys
import csv


def add_bam(parser):
    parser.add_argument("bam", help="Path to BAM file")


def add_region(parser):
    parser.add_argument(
        "region", help="Region to view, specified in the form CHROM:START-END"
    )


def add_index(parser):
    parser.add_argument(
        "--index",
        default=None,
        help="Path to index (BAI) file (default: </path/to/bam>.bai)",
    )


def add_mapq(parser):
    parser.add_argument(
        "--mapq",
        type=int,
        default=0,
        help="Minimum mapping quality (default: %(default)s)",
    )


def add_baseq(parser):
    parser.add_argument(
        "--baseq",
        type=int,
        default=0,
        help="Minimum base quality (default: %(default)s)",
    )


def add_stats(parser):
    parser.add_argument(
        "--stats",
        action="store_true",
        default=False,
        help="Output additional per-position statistics",
    )


def add_decimals(parser):
    parser.add_argument(
        "--decimals",
        type=int,
        default=3,
        help="Number of decimal places to display (default: %(default)s)",
    )


def get_entropy(probabilities, normalised=False):
    entropy = sum([-(x * math.log2(x)) if x != 0 else 0 for x in probabilities])

    if normalised:
        return entropy / math.log2(len(probabilities))
    else:
        return entropy


def get_stats(counts, decimals=3):
    coverage = sum(counts)
    probabilities = [count / coverage if coverage > 0 else 0.0 for count in counts]
    percentages = [100 * probability for probability in probabilities]
    entropy = get_entropy(probabilities, normalised=True)
    secondary_count = list(counts)
    secondary_count.pop(counts.index(max(counts)))
    secondary_coverage = sum(secondary_count)
    secondary_probabilities = [
        count / secondary_coverage if secondary_coverage > 0 else 0.0
        for count in secondary_count
    ]
    secondary_entropy = get_entropy(secondary_probabilities, normalised=True)
    return (
        [coverage]
        + counts
        + [round(x, decimals) for x in percentages + [entropy, secondary_entropy]]
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

    command = parser.add_subparsers(dest="command", required=True)

    all_parser = command.add_parser("all", help="Count bases at all positions.")
    add_bam(all_parser)
    add_mapq(all_parser)
    add_baseq(all_parser)
    add_stats(all_parser)
    add_decimals(all_parser)

    query_parser = command.add_parser(
        "query", help="Count bases in specific region (without an index file)."
    )
    add_bam(query_parser)
    add_region(query_parser)
    add_mapq(query_parser)
    add_baseq(query_parser)
    add_stats(query_parser)
    add_decimals(query_parser)

    iquery_parser = command.add_parser(
        "iquery", help="Count bases in specific region (with an index file)."
    )
    add_bam(iquery_parser)
    add_region(iquery_parser)
    add_index(iquery_parser)
    add_mapq(iquery_parser)
    add_baseq(iquery_parser)
    add_stats(iquery_parser)
    add_decimals(iquery_parser)

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

    if args.command == "all":
        data = all(args.bam, args.mapq, args.baseq)

        for row in iterate(data, stats=args.stats, decimals=args.decimals):
            writer.writerow(row)

    elif args.command == "query":
        data = query(args.bam, args.region, args.mapq, args.baseq)

        for row in iterate(
            data, region=args.region, stats=args.stats, decimals=args.decimals
        ):
            writer.writerow(row)

    elif args.command == "iquery":
        data = iquery(args.bam, args.region, args.index, args.mapq, args.baseq)

        for row in iterate(
            data, region=args.region, stats=args.stats, decimals=args.decimals
        ):
            writer.writerow(row)
