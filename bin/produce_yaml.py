#!/usr/bin/env python3

import argparse
from pathlib import Path

class Sample:
    def __init__(
        self,
        analysis_type: str,
        sample_id: str,
        sample_name: str,
        sex: str,
        phenotype: str,
        tissue_type: str,
        bam_path: Path,
        rna_alignment_path: Path,
        rna_coverage_bigwig: Path,
        splice_junctions_bed: Path
    ):
        self.analysis_type = analysis_type
        self.sample_id = sample_id
        self.sample_name = sample_name
        self.sex = sex
        self.phenotype = phenotype
        self.tissue_type = tissue_type
        self.bam_path = bam_path
        self.rna_alignment_path = rna_alignment_path
        self.rna_coverage_bigwig = rna_coverage_bigwig
        self.splice_junctions_bed = splice_junctions_bed
    
    def __str__(self) -> str:
        fields = [
            f"- analysis_type: {self.analysis_type}",
            f"  sample_id: {self.sample_id}",
            f"  sample_name: {self.sample_name}",
            f"  sex: {self.sex}"
            f"  phenotype: {self.phenotype}"
            f"  tissue_type: {self.tissue_type}"
            f"  bam_path: {self.bam_path}"
            f"  rna_alignment_path: {self.rna_alignment_path}"
            f"  rna_coverage_bigwig: {self.rna_coverage_bigwig}"
            f"  splice_junctions_bed: {self.splice_junctions_bed}"
        ]
        return "\n".join(fields)
        

def main(
    sample_id:str,
    sex: str,
    phenotype: str,
    tissue: str,
    bam_path: Path,
    rna_bigwig: Path,
    splice_junctions: Path,
    vcf: Path,
    fraser: Path,
    outrider: Path,
):
    yaml_dict = {
        'owner': 'rnaseq',
        'family': sample_id,
        'family_name': sample_id,
        'synopsis': ['First batch of Tomte samples'],
        'samples': [
            str(Sample(
                analysis_type="wgs",
                sample_id=sample_id,
                sample_name=sample_id,
                sex=sex,
                phenotype=phenotype,
                tissue_type=tissue,
                bam_path=bam_path,
                rna_alignment_path=bam_path,
                rna_coverage_bigwig=rna_bigwig,
                splice_junctions_bed=splice_junctions,
            ))
        ],
        'vcf_snv': vcf,
        'omics_files': [
            f'fraser: ${fraser}',
            f'outrider: ${outrider}',
        ],
        "default_gene_panels": "[]",
        "gene_panels": "[]"
    }


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample_id")
    parser.add_argument("--tomte_results")
    parser.add_argument("--fraser_path")
    parser.add_argument("--outrider_path")
    parser.add_argument("--vcf_path")
    parser.add_argument("--sex")
    parser.add_argument("--phenotype")
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_arguments()
    main(
        args.sample_id,
        args.sex,
        args.phenotype,
        args.tomte_results,
        args.fraser_path,
        args.outrider_path,
        args.vcf_path
    )
