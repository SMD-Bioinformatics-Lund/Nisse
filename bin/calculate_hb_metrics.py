#!/usr/bin/env python3

import argparse
from pathlib import Path
from typing import Optional, Set
import csv
from collections import defaultdict


def main(hbs_path: Path, gene_counts: Path, strandedness: Optional[str]):

    hb_genes = get_hb_genes(hbs_path)

    if strandedness not in ["forward", "reverse"] and strandedness is not None:
        raise ValueError(
            f"strandedness option should be forward, reverse or not assinged. Found: {strandedness}"
        )

    nbr_rows = 0
    counts = defaultdict(int)
    with open(gene_counts, "r") as in_fh:
        for line in in_fh:
            (gene_id_raw, both, fw, rv) = line.rstrip().split("\t")
            gene_id = gene_id_raw.split(".")[0]

            if gene_id.startswith("N_"):
                print(f"Skipping N_ prefixed row: {line.rstrip()}")
                continue

            if strandedness == "reverse":
                value = int(rv)
            elif strandedness == "forward":
                value = int(fw)
            else:
                value = int(both)

            counts["all"] += value
            if gene_id in hb_genes:
                counts["hb"] += value
                counts[gene_id] += value

            nbr_rows += 1

    # FIXME: Wrap the output here. Write to file.
    print(f"Nbr rows processed: {nbr_rows}")
    hb_perc = 100 * (counts["hb"] / counts["all"])
    print(counts)
    print(f"Perc hb: {round(hb_perc, 3)}%")


def get_hb_genes(hb_map_path: Path) -> Set[str]:
    hb_genes = set()
    with hb_map_path.open("r") as hb_in:
        hb_genes_reader = csv.DictReader(hb_in, delimiter="\t")
        for row in hb_genes_reader:
            ensembl_val = row.get("ensembl")
            if ensembl_val is None:
                raise ValueError(f"Expected ensembl value, found None for row: {row}")
            hb_genes.add(ensembl_val)
    return hb_genes


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--hb_map", required=True, help="TSV file with one header 'ensembl'")
    parser.add_argument("--gene_counts", required=True, help="htseq-lib like output")
    parser.add_argument(
        "--strandedness", required=True, help="forward, reverse or none", default=None
    )
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_arguments()
    main(Path(args.hb_map), Path(args.gene_counts), args.strandedness)
