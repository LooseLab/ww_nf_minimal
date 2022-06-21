#!/usr/bin/env python
"""collate_lineages.py

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

"""
results/liverpool/220107/freyja/6_seqloc_G4_2136511_G4_2109758_32.demix.txt
3:lineages      ['BA.1.1' 'B.1.1.529' 'BA.1' 'B.1.1.161']
4-abundances    [0.27652267 0.27652267 0.27652267 0.170432  ]
"""

META_PREFIX = "meta_"
LINEAGE_PATTERN = re.compile(r"^lineages\s+(.*)abundances", re.MULTILINE | re.DOTALL)
ABUNDANCE_PATTERN = re.compile(r"^abundances\s+(.*)resid", re.MULTILINE | re.DOTALL)


def get_parser():
    p = argparse.ArgumentParser(Path(__file__).name)
    p.add_argument("demix", help="Output from `freyja demix` command", type=Path)
    p.add_argument("-o", "--output", type=Path)
    p.add_argument("--meta-sample-id")
    p.add_argument("--meta-sample-site-code")
    p.add_argument("--meta-timestamp-sample-collected")
    p.add_argument("--meta-sequencing-lab-code")
    p.add_argument("--meta-original-sample-id")
    p.add_argument("--meta-sequencing-sample-id")
    p.add_argument("--meta-sequencing-run-id")
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
    lineages_str = search_file(args.demix, LINEAGE_PATTERN)
    lineages_str = lineages_str.strip().replace('"', "")
    lineages = literal_eval(re.sub("\s+", ",", lineages_str))
    abundances_str = search_file(args.demix, ABUNDANCE_PATTERN)
    abundances_str = abundances_str.strip().replace('"', "")
    abundances = literal_eval(re.sub("\s+", ",", abundances_str))

    lineages = dict(zip(lineages, abundances))

    logger.info(f"Calculate amplicon stats statistics")
    record = {
        "lineage": None,
        "abundance": None,
    }
    record.update(meta_vals)

    logger.info(f"Writing CSV {args.output}")
    with open(args.output, mode="w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=record.keys())
        writer.writeheader()
        if not lineages:
            writer.writerow(record)
        else:
            for variant, abundance in lineages.items():
                record["lineage"] = variant
                record["abundance"] = abundance
                writer.writerow(record)
    logger.info(f"Written CSV, Goodbye!")


if __name__ == "__main__":
    main()
