#!/usr/bin/env python3

import argparse
from pathlib import Path
from typing import List

DESCRIPTION = """Generate YAML config used to load sample into Scout"""

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
    sample_id: str,
    sample_description: str,
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
    not_found: list[Path] = []
    if not bam_path.exists():
        print(f"File in {bam_path} not found (bam_path)")
        not_found.append(bam_path)
    if not rna_bigwig.exists():
        print(f"File in {rna_bigwig} not found (rna_bigwig)")
        not_found.append(rna_bigwig)
    if not splice_junctions.exists():
        print(f"File in {splice_junctions} not found (splice_junctions)")
        not_found.append(splice_junctions)
    if not vcf.exists():
        print(f"File in {vcf} not found (vcf)")
        not_found.append(vcf)
    if not fraser.exists():
        print(f"File in {fraser} not found (fraser)")
        not_found.append(fraser)
    if not outrider.exists():
        print(f"File in {outrider} not found (outrider)")
        not_found.append(outrider)

    if len(not_found) > 0:
        not_found_str = '\n'.join([str(p) for p in not_found])
        raise ValueError(f"All required paths not present. Missing:\n{not_found_str}")

    yaml_dict = {
        'owner': 'rnaseq',
        'family': sample_id,
        'family_name': sample_id,
        'synopsis': [sample_description],
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
            f'fraser: {fraser}',
            f'outrider: {outrider}',
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
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument("--sample_id", required=True)
    parser.add_argument("--description", required=True, help="Synopsis field in YAML")
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--sex", required=True)
    parser.add_argument("--phenotype", required=True)
    parser.add_argument("--tissue", required=True)
    parser.add_argument("--fraser", required=True)
    parser.add_argument("--outrider", required=True)
    parser.add_argument("--bam_path", required=True)
    parser.add_argument("--splice_junctions", required=True)
    parser.add_argument("--rna_bigwig", required=True)
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_arguments()
    main(
        sample_id=args.sample_id,
        sample_description=args.description,
        vcf=Path(args.vcf),
        sex=args.sex,
        phenotype=args.phenotype,
        tissue=args.tissue,
        fraser=Path(args.fraser),
        outrider=Path(args.outrider),
        bam_path=Path(args.bam_path),
        splice_junctions=Path(args.splice_junctions),
        rna_bigwig=Path(args.rna_bigwig),
    )
