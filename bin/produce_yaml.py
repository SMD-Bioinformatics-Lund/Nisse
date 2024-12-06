#!/usr/bin/env python3

import argparse
from pathlib import Path
from typing import List

def get_space(level: int) -> str:
    return " " * 4 * level

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
    
    def get_lines(self) -> List[str]:

        keys = self.__dict__.keys()
        fields = [f"{key}: {self.__dict__[key]}" for key in keys]
        return fields
        

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
            Sample(
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
            )
        ],
        'vcf_snv': vcf,
        'omics_files': [
            f'fraser: ${fraser}',
            f'outrider: ${outrider}',
        ],
        "default_gene_panels": "[]",
        "gene_panels": "[]",
        "human_genome_build": "'38'",
        "rna_human_genome_build": "'38'"
    }

    for (key, value) in yaml_dict.items():
        if isinstance(value, list):
            value_it = value
            print(f"{key}:")
            for val in value_it:
                if isinstance(val, Sample):
                    sample_rows = val.get_lines()
                    is_first = True
                    for row in sample_rows:
                        prefix = "  "
                        if is_first:
                            prefix = "- "
                            is_first = False
                        print(get_space(1) + prefix + row)
                else:
                    print(get_space(1) + val)
        elif isinstance(value, Sample):
            my_sample = value
            print(f"{key}:")
            for val in my_sample.get_lines():
                print(get_space(1) + val)
        else:
            print(f"{key}: {value}")


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample_id", required=True)
    parser.add_argument("--tomte_results", required=True)
    parser.add_argument("--fraser_path", required=True)
    parser.add_argument("--outrider_path", required=True)
    parser.add_argument("--vcf_path", required=True)
    parser.add_argument("--sex", required=True)
    parser.add_argument("--phenotype", required=True)
    parser.add_argument("--tissue", required=True)
    parser.add_argument("--fraser", required=True)
    parser.add_argument("--outrider", required=True)
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_arguments()
    main(
        args.sample_id,
        args.sex,
        args.phenotype,
        args.tissue,
        args.tomte_results,
        args.fraser_path,
        args.outrider_path,
        args.vcf_path,
        args.fraser,
        args.outrider,
    )
