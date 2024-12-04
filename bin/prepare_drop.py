import argparse
from pathlib import Path
import csv
import sys

"""
Takes a single FRASER or OUTRIDER file
Filters out rows below statistical threshold (p or FDR)
Splits it per sample into out_dir
"""


def main(
    in_path: Path,
    out_dir: Path,
    stat_col: str,
    stat_cutoff: float,
    hgnc_symbol_col: str,
    sample_col: str,
):
    rows_all = []
    rows_to_keep = []
    print(f"Initial number of entries, counting file: {in_path}")
    with open(in_path) as in_fh:
        reader = csv.DictReader(in_fh, delimiter="\t")
        for row in reader:
            rows_all.append(row)
            if float(row[stat_col]) < stat_cutoff:
                rows_to_keep.append(row)
    print(f"All rows: {len(rows_all)}")
    print(f"Number of rows left: {len(rows_to_keep)}")


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--in_path", required=True)
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--stat_col", required=True)
    parser.add_argument("--stat_cutoff", required=True, type=float)
    parser.add_argument("--hgnc_symbol_col", required=True)
    parser.add_argument("--sample_col", required=True)
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_arguments()
    main(
        args.in_path,
        args.out_dir,
        args.stat_col,
        args.stat_cutoff,
        args.hgnc_symbol_col,
        args.sample_col,
    )
