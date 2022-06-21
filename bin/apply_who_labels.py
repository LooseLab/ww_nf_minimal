#!/usr/bin/env python
"""collate_lineages.py

"""
import sys
import argparse
import logging
from pathlib import Path, PurePath

import pandas as pd


def get_parser():
    p = argparse.ArgumentParser(Path(__file__).name)
    p.add_argument("lineages", help="all_lineages.csv", type=Path)
    p.add_argument("who_map", help="VOC lineage map", type=Path)
    p.add_argument("sub_map", help="Sub-lineage map", type=Path)
    # p.add_argument("-o", "--output", type=argparse.FileType("w"), default=sys.stdout)
    p.add_argument("-o", "--output", type=Path)
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


def main():
    p = get_parser()
    args = p.parse_args()
    logging.basicConfig(
        level=LOG_LEVELS.get(args.verbose, logging.DEBUG),
        format="[%(asctime)s - %(levelname)s] %(message)s",
    )
    logger = logging.getLogger(__name__)

    for k, v in vars(args).items():
        logger.info(f"{k}={v}")

    logger.info(f"Loading {args.lineages}")
    df = pd.read_csv(args.lineages)
    logger.info(f"Loaded {args.lineages}")

    logger.info(f"Loading {args.who_map}")
    with open(args.who_map, "r") as fh:
        d = dict(l.strip().split(",") for l in fh)
    logger.info(f"Loaded {args.who_map}")

    logger.info(f"Applying WHO map")
    df["WHO_lineage"] = df["lineage"].map(d)
    logger.info(f"Applied WHO map")

    logger.info(f"Loading {args.sub_map}")
    with open(args.sub_map, "r") as fh:
        d = dict(l.strip().split(",") for l in fh)
    logger.info(f"Loaded {args.sub_map}")

    logger.info(f"Applying sub-lineage map")
    df["parent_WHO_lineage"] = df["lineage"].map(d)
    logger.info(f"Applied sub-lineage map")

    logger.info(f"Writing CSV")
    df.to_csv(args.output, index=False)
    logger.info(f"Written CSV, Goodbye!")


if __name__ == "__main__":
    main()
