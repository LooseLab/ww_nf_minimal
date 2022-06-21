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
    p.add_argument("variant_table", help="Variant table CSV", type=Path)
    p.add_argument("coverage", help="Per-base coverage tsv", type=Path)
    p.add_argument("acgt", help="_acgt.tsv file", type=Path)
    p.add_argument("indel", help="_varscan.indel.tsv", type=Path)
    p.add_argument("snps", help="SNP list, `|` delimited file", type=Path)
    p.add_argument("-o", "--output", type=Path)
    p.add_argument("--sample", type=Path)
    p.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite the output file, if output already exists it won't be overwritten",
    )
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


def get_cols(row, snp_list):
    if row["VarPos"] in snp_list:
        if len(row["VarPattern"]) == 1:
            return row[row["RefPattern"]], row[row["VarPattern"]], row["VarName"]
        else:
            return row[row["RefPattern"]], 0, "other2"
    else:
        if all(pd.isna(row[x]) for x in ["Cons", "Reads1", "Reads2"]):
            # No indel at this position
            return row["coverage"], 0, row["VarName"]
        else:
            # Might be an indel
            if row["VarPattern"] in row["Cons"]:
                return row["Reads1"], row["Reads2"], row["VarName"]
            else:
                return row["Reads1"], 0, "other"
    return 0, 0, row["VarName"]


def main():
    p = get_parser()
    args = p.parse_args()
    logging.basicConfig(
        level=LOG_LEVELS.get(args.verbose, logging.DEBUG),
        format="[%(asctime)s - %(levelname)s] %(message)s",
    )
    logger = logging.getLogger(__name__)
    for flag, value in vars(args).items():
        # Log our command line params
        logger.info(f"{flag}={value}")
        # Get the meta-data

    if not args.overwrite and args.output.is_file():
        msg = f"Output file {args.output} already exists! Use `--overwrite` to overwrite it"
        logger.error(msg)
        sys.exit(msg)

    col_map = {
        "VarPos": "MutPos",
        "VarID": "MutID",
        "RefPattern": "Refbase",
        "VarPattern": "VarPatt",
        "VarName": "MutName",
        "coverage": "TotReads",
    }
    out_order = [
        "Sample",
        "MutPos",
        "MutID",
        "TotReads",
        "Refbase",
        "VarPatt",
        "ReadsRef",
        "ReadsAlt",
        "MutName",
    ]

    df = pd.read_csv(args.variant_table)
    acgt = pd.read_csv(
        args.acgt, sep="\t", usecols=["n_base", "A", "C", "G", "T"]
    ).rename(columns={"n_base": "VarPos"})
    cov = pd.read_csv(
        args.coverage, sep="\t", usecols=["position", "coverage"]
    ).rename(columns={"position": "VarPos"})
    indel = pd.read_csv(
        args.indel,
        sep="\t",
        usecols=["Position", "Cons", "Reads1", "Reads2"],
    ).rename(columns={"Position": "VarPos"})

    with open(args.snps) as fh:
        snp_list = set(map(int, fh.read().split("|")))

    merge_params = {"how": "left", "on": "VarPos"}
    df = (
        df.merge(acgt, **merge_params)
        .merge(cov, **merge_params)
        .merge(indel, **merge_params)
    )

    df["Sample"] = args.sample
    df[["ReadsRef", "ReadsAlt", "VarName"]] = df.apply(get_cols, 1, result_type="expand", snp_list=snp_list)

    df = df.rename(columns=col_map)
    df[out_order].to_csv(args.output, sep="\t", index=False)
    logger.info(f"Written CSV, Goodbye!")


if __name__ == "__main__":
    main()
