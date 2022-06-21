#!/usr/bin/env python
"""collate_results.py

Produce the output CSV file (melted) containing lineage and abundance/prop, and
some of the amplicon stats for each of those samples

 rg --no-filename --no-heading -g '*ampliconstats.txt' "^(FDEPTH|FVDEPTH)" > amp.stats
"""
import sys
from ast import literal_eval
import argparse
import logging
from pathlib import Path, PurePath
import re
import statistics
import csv
from contextlib import contextmanager


META_PREFIX = "meta_"
DEMIX_PATTERN = re.compile(r"^summarized\s+(.*)$", re.MULTILINE)
AMPLICON_STATS_PATTERN = re.compile(rf"^(FDEPTH.*)$", re.MULTILINE) 


def get_parser():
    p = argparse.ArgumentParser(Path(__file__).name)
    p.add_argument("demix", help="Output from `freyja demix` command", type=Path)
    p.add_argument("ampliconstats", help="Output from `samtools ampliconstats` command", type=Path)
    p.add_argument("mosdepth", help="Mosdepth summary output `{prefix}.mosdepth.summary.txt`", type=Path)
    p.add_argument("-o", "--output", type=Path)
    p.add_argument("--meta-sample-id")
    p.add_argument("--meta-sample-site-code")
    p.add_argument("--meta-timestamp-sample-collected")
    p.add_argument("--meta-sequencing-lab-code")
    p.add_argument("--meta-original-sample-id")
    p.add_argument("--meta-sequencing-sample-id")
    p.add_argument("--meta-sequencing-run-id")
    p.add_argument("--coverage-threshold", default=20, type=int)
    p.add_argument("--overwrite", action="store_true", help="Overwrite the output file, if output already exists it won't be overwritten")
    p.add_argument(
        "--verbose",
        "-v",
        action="count",
        default=0,
        help="increase verbosity (-v, -vv, ...)",
    )
    p.add_argument("--version", "-V", action="version", version="%(prog)s 0.0.1")
    return p


LOG_LEVELS = {
    0: logging.CRITICAL,
    1: logging.ERROR,
    2: logging.WARN,
    3: logging.INFO,
    4: logging.DEBUG,
}


def catcher(logger, *exceptions):
    def decorator(func):
        def inner(*args, **kwargs):
            try:
                return func(*args, **kwargs)
            except exceptions as err:
                msg = getattr(err, "message", str(err))
                logger.error(msg)
                sys.exit(msg if msg else 1)
            except:
                raise
        return inner
    return decorator


@catcher(logging.getLogger(__name__), ValueError, FileNotFoundError)
def search_file(filepath, pattern):
    """Search a file for a pattern, returning the first group.
    """
    if not isinstance(filepath, PurePath):
        filepath = Path(filepath)

    if filepath.is_file():
        if match := pattern.search(filepath.open().read()):
            return match.group(1)
        else:
            raise ValueError(f"Pattern {pattern.pattern!r} not found in {filepath}")
    else:
        raise FileNotFoundError(f"File {filepath} not found")


def main():
    p = get_parser()
    args = p.parse_args()
    logging.basicConfig(
        level=LOG_LEVELS.get(args.verbose, logging.DEBUG),
        format="[%(asctime)s - %(levelname)s] %(message)s",
    )
    logger = logging.getLogger(__name__)
    meta_vals = {}
    for flag, value in vars(args).items():
        # Log our command line params
        logger.info(f"{flag}={value}")
        # Get the meta-data
        if flag.startswith(META_PREFIX):
            meta_vals[flag.removeprefix(META_PREFIX)] = value

    if not args.overwrite and args.output.is_file():
        msg = f"Output file {args.output} already exists! Use `--overwrite` to overwrite it"
        logger.error(msg)
        sys.exit(msg)

    logger.info(f"Parsing {args.demix}")
    lineages = literal_eval(search_file(args.demix, DEMIX_PATTERN))
    lineages = dict(lineages)

    logger.info(f"Parsing {args.ampliconstats}")
    fdepth_line = search_file(args.ampliconstats, AMPLICON_STATS_PATTERN)

    _, bam_file, stats = fdepth_line.split(maxsplit=2)
    amplicons = literal_eval(re.sub(r"\s", ",", stats))
    non_zero_amplicons = list(filter(lambda v: v > 0, amplicons))
    if len(non_zero_amplicons) <= 1:
        logger.warning("Zero or one amplicons has coverage")

    amplicon_column_name = "amplicon_{:>03}_mean_depth"
    amplicon_depths = {amplicon_column_name.format(i): d for i, d in enumerate(amplicons, start=1)}

    record = {}
    logger.info(f"Parsing {args.mosdepth}")
    with open(args.mosdepth, newline="") as mosdepth_summary:
        # Assume that the first record is the one that we are interested in
        #   as we are only using a single reference genome with a single contig
        # We could check that the reference name matches a the contig names in 
        #   the ampliconstats file
        record = next(csv.DictReader(mosdepth_summary, delimiter="\t"))

    logger.info(f"Calculate amplicon stats statistics")
    amplicon_features = {
        "amplicon_mean": statistics.fmean(amplicons),
        "non_zero_amplicon_mean": statistics.fmean(non_zero_amplicons) if non_zero_amplicons else None,
        "amplicon_median": statistics.median(amplicons),
        "non_zero_amplicon_median": statistics.median(non_zero_amplicons) if non_zero_amplicons else None,
        f"count_gte_{args.coverage_threshold}": sum(v >= args.coverage_threshold for v in amplicons),
        f"count_lt_{args.coverage_threshold}": sum(v < args.coverage_threshold for v in amplicons),
        "stdev": statistics.stdev(amplicons),
        "non_zero_stdev": statistics.stdev(non_zero_amplicons) if len(non_zero_amplicons) > 1 else None,
        "lineage": None,
        "abundance": None,
        "mean_genome_coverage": record.get("mean", "")
    }
    amplicon_features.update(meta_vals)
    amplicon_features.update(amplicon_depths)

    logger.info(f"Writing CSV {args.output}")
    with open(args.output, mode="w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=amplicon_features.keys())
        writer.writeheader()
        if not lineages:
            writer.writerow(amplicon_features)
        else:
            for variant, abundance in lineages.items():
                amplicon_features["lineage"] = variant
                amplicon_features["abundance"] = abundance
                writer.writerow(amplicon_features)
    logger.info(f"Written CSV, Goodbye!")


if __name__ == "__main__":
    main()
