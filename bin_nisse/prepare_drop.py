#!/usr/bin/env python3

import argparse
from pathlib import Path
import csv

"""
Takes a single FRASER or OUTRIDER file
Filters out rows below statistical threshold (p or FDR)
Splits it per sample into out_dir
"""


def main(
    in_path: Path,
    out_path: Path,
    sample_col: str,
    sample: str,
    stat_col: str,
    stat_cutoff: float,
    hgnc_symbol_col: str,
    hgnc_symbol_id_map: str,
    verbose: bool,
):
    
    # Add HGNC ID
    hgnc_symbol_id_dict = {}
    with open(hgnc_symbol_id_map, "r") as in_fh:
        for line in in_fh:
            line = line.rstrip()
            (hgnc_symbol, hgnc_id) = line.split("\t")
            hgnc_symbol_id_dict[hgnc_symbol] = hgnc_id

    stats = {"all": 0, "stat": 0, "hgnc_symbol": 0, "hgnc_id": 0}

    output_rows = []
    print(f"Initial number of entries, counting file: {in_path}")
    headers = None
    with open(in_path) as in_fh:
        reader = csv.DictReader(in_fh, delimiter="\t")
        for row in reader:
            if headers == None:
                headers = list(row.keys()) + ["hgncId"]

            if row[sample_col] != sample:
                continue

            pass_stat_cutoff = float(row[stat_col]) < stat_cutoff
            has_hgnc_symbol = not row[hgnc_symbol_col].startswith("ENS")
            hgnc_symbol = row[hgnc_symbol_col]
            hgnc_id = hgnc_symbol_id_dict.get(hgnc_symbol)
            has_hgnc_id = hgnc_id is not None

            stats['all'] += 1
            if pass_stat_cutoff:
               stats['stat'] += 1
            if has_hgnc_symbol:
                stats['hgnc_symbol'] += 1
            if has_hgnc_id:
                stats['hgnc_id'] += 1

            if pass_stat_cutoff and has_hgnc_symbol and has_hgnc_id:
                row['hgncId'] = hgnc_id
                output_rows.append(row)
    
    if not headers:
        raise ValueError("No headers found. Is the file empty?")
    
    print(f"Writing {len(output_rows)} rows to {out_path}")
    with open(out_path, "w", newline='') as out_fh:
        writer = csv.DictWriter(out_fh, fieldnames=headers, delimiter="\t")

        writer.writeheader()
        writer.writerows(output_rows)

    if verbose:
        print("Number passing each filter")
        for key, count in stats.items():
            print(f"{key}: {count}")
    

def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()

    parser.add_argument("--in_path", required=True)
    parser.add_argument("--out_path", required=True)

    parser.add_argument("--sample_col", required=True, help="Only return hits for the specific sample")
    parser.add_argument("--sample", required=True)

    parser.add_argument("--stat_col", required=True)
    parser.add_argument("--stat_cutoff", required=True, type=float)
    parser.add_argument("--hgnc_symbol_col", required=True)

    parser.add_argument("--hgnc_symbol_id_map", required=True)
    parser.add_argument("--verbose", action="store_true")

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_arguments()
    main(
        args.in_path,
        args.out_path,
        args.sample_col,
        args.sample,
        args.stat_col,
        args.stat_cutoff,
        args.hgnc_symbol_col,
        args.hgnc_symbol_id_map,
        args.verbose
    )
