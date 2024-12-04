import argparse
from pathlib import Path
import csv
import sys
from collections import defaultdict

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
    hgnc_symbol_id_map: str,
):
    
    # Add HGNC ID
    hgnc_symbol_id_dict = {}
    with open(hgnc_symbol_id_map, "r") as in_fh:
        for line in in_fh:
            line = line.rstrip()
            (hgnc_symbol, hgnc_id) = line.split("\t")
            hgnc_symbol_id_dict[hgnc_symbol] = hgnc_id

    stats = defaultdict(int)

    output_rows = []
    print(f"Initial number of entries, counting file: {in_path}")
    with open(in_path) as in_fh:
        reader = csv.DictReader(in_fh, delimiter="\t")
        for row in reader:
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
                row['hgncID'] = hgnc_id
                output_rows.append(row)
        
        print(stats)
            # print(hgnc_id)
            # print(row)
    


    # print(f"All rows: {len(rows_all)}")
    # print(f"Number of rows left: {len(rows_passing_stat)}")


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--in_path", required=True)
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--stat_col", required=True)
    parser.add_argument("--stat_cutoff", required=True, type=float)
    parser.add_argument("--hgnc_symbol_col", required=True)
    parser.add_argument("--sample_col", required=True)

    parser.add_argument("--hgnc_symbol_id_map", required=True)

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
        args.hgnc_symbol_id_map
    )
