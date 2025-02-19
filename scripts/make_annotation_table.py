#!/usr/bin/env python3

import pandas as pd
import argparse
from pathlib import Path

DESCRIPTION = """Prepare annotation table for DROP database"""

HEADERS = {
    "RNA_ID": None,
    "RNA_BAM_FILE": "NA",
    "DNA_VCF_FILE": "NA",
    "DNA_ID": "NA",
    "DROP_GROUP": "outrider,fraser",
    "PAIRED_END": "True",
    "COUNT_MODE": "IntersectionStrict",
    "COUNT_OVERLAPS": "True",
    "SPLICE_COUNTS_DIR": None,
    "STRAND": None,
    "HPO_TERMS": "NA",
    "GENE_ANNOTATION": None,
    "GENOME": "NA",
    "SEX": None,
    "GENE_COUNTS_FILE": None,
}


def main(
    source_tsv: str,
    rna_id_col: str,
    sex_col: str,
    splice_counts_dir: str,
    gene_counts_file: str,
    strand: str,
    gencode_annotation: str,
    output_path: str,
):
    # Initial columns
    annotation_df = pd.read_csv(source_tsv, sep="\t")[[rna_id_col, sex_col]]

    # Get the names right
    annotation_df = annotation_df.rename({rna_id_col: "RNA_ID", sex_col: "SEX"})

    # Fill in known constants
    for header, fill_value in HEADERS.items():
        if fill_value is not None:
            annotation_df.loc[:, header] = fill_value

    # Assign values from arguments
    annotation_df.loc[:, "SPLICE_COUNTS_DIR"] = splice_counts_dir
    annotation_df.loc[:, "GENE_COUNTS_FILE"] = gene_counts_file
    annotation_df.loc[:, "STRAND"] = strand
    annotation_df.loc[:, "GENE_ANNOTATION"] = gencode_annotation

    annotation_df.to_csv(output_path, sep="\t", index=False)

def parse_arguments():
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument(
        "--source_tsv", required=True, help="Table containing sample IDs and sex information"
    )
    parser.add_argument("--rna_id_col", required=True, help="Used for field RNA_ID")
    parser.add_argument("--sex_col", required=True, help="Used for field SEX")
    parser.add_argument(
        "--splice_counts_dir", required=True, help="Used for field SPLICE_COUNTS_DIR"
    )
    parser.add_argument("--gene_counts_file", required=True, help="Used for field GENE_COUNTS_FILE")
    parser.add_argument("--strand", required=True, help="forward or reverse, used for field STRAND")
    parser.add_argument(
        "--gencode_annotation",
        required=True,
        help="Needs to match base name of used annotation file, i.e. params.gtf=gencode.v46.annotation.gtf.gz means gencode.v46.annotation. Used for field GENE_ANNOTATION",
    )
    parser.add_argument("--output", required=True, help="Output path")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_arguments()

    if not Path(args.source_tsv).exists():
        raise ValueError("--source_tsv must exist")

    if not Path(args.splice_counts_dir).exists():
        raise ValueError("--splice_counts_dir must exist")

    if not Path(args.gene_counts_file).exists():
        raise ValueError("--gene_counts_file must exist")


    main(
        args.source_tsv,
        args.rna_id_col,
        args.sex_col,
        args.splice_counts_dir,
        args.gene_counts_file,
        args.strand,
        args.gencode_annotation,
        args.output
    )
