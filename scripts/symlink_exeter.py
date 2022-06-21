#!/usr/bin/env python3

import os
import sys
import csv
import re
import argparse
from pathlib import Path

from tqdm import tqdm

LABS = {
    "NOTTINGHAM_UNI": "nottingham",
    "LIVERPOOL_UNI": "liverpool",
    "EXETER_UNI": "exeter",
    "EA": "nottingham",
}

def get_parser():
    p = argparse.ArgumentParser()
    p.add_argument("csv", help="Input CSV file")
    p.add_argument("src", help="Source folder of FASTQ files", type=Path)
    p.add_argument("dest", help="Destination folder for created symlinks", type=Path)
    return p

def main():
    p = get_parser()
    args = p.parse_args()
    print(args)

    args.dest.mkdir(parents=True, exist_ok=True)
    files = set(args.src.rglob("*.fastq.gz"))
    print(len(files))

    with open(args.csv, newline="") as f:
        reader = csv.DictReader(f)
        for row in tqdm(reader):
            # Get run ID
            run_id = row["sequencing_run_id"]
            # Set and create destination directory
            dst = args.dest / run_id
            dst.mkdir(parents=True, exist_ok=True)
            # Create pattern for finding matching FASTQ file
            sample = row["sequencing_sample_id"]
            loc = f"^{args.src}/.*[\d]{{8}}/{run_id}_{sample}.*R[12].*\.fastq\.gz$"
            # loc = f"^{args.src}/.*[\d]{{8}}/{run_id}_{sample}.*R[12].*\.fastq\.gz$"
            fn_pat = re.compile(loc)
            # Create symlink pattern
            sym = f"{dst}/{{}}"
            found = []

            # Iterate all found FASTQ saving matching paths
            for f in files:
                if fn_pat.match(str(f)):
                    found.append((f, f.name.removeprefix(f"{run_id}_")))

            # Create the symlinks and remove from found FASTQ
            for s, d in found:
                d = Path(sym.format(d))
                if not d.is_file():
                    d.symlink_to(s)
                # discard (no raise) seen files
                files.discard(s)

if __name__ == "__main__":
    main()
