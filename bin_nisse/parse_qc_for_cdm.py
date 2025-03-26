#!/usr/bin/env python3

import argparse
import json
import numpy as np
from ast import literal_eval
from pathlib import Path
from collections import defaultdict
from typing import Any, Dict, List, Mapping, Set, Tuple, Union, Optional

VERSION = "1.1.0"
DESCRIPTION = """
Parse QC-output from Tomte and Nisse into a JSON format
compatible with CDM.
"""


class QCEntry:
    def __init__(
        self, results: Dict[str, Any], label: Optional[str], software: str, version: str, url: str
    ):
        self.results = results

        self.label = label
        self.software = software
        self.version = version
        self.url = url

    def get_result_dict(self) -> Dict[str, Any]:

        entry_dict: Dict[str, Any] = {}

        if self.label:
            entry_dict["label"] = self.label

        if self.software:
            entry_dict["software"] = self.software

        if self.version:
            entry_dict["version"] = self.version

        if self.url:
            entry_dict["url"] = self.url

        entry_dict["results"] = self.results

        return entry_dict


def main(
    multiqc_general_stats: str,
    picard_rna_coverage: Optional[str],
    multiqc_star: Optional[str],
    merged_hb_estimate: Optional[str],
    hetcalls_vcfs: Optional[List[str]],
    output_file: str,
    sample_id: Optional[str],
    debug: bool,
):
    samples_qc_dict = {}

    sample_qc_tups = parse_multiqc_general_stats(multiqc_general_stats)
    for sample, qc_dict in sample_qc_tups:
        qc_entry = QCEntry(
            qc_dict, "MultiQC general stats", "MultiQC", "No version", "https://seqera.io/multiqc/"
        )
        samples_qc_dict[sample] = [qc_entry]

    if merged_hb_estimate:
        # Add HB estimate data
        hb_data = process_hb_estimate_data(merged_hb_estimate)

        for sample, hb_values in hb_data.items():
            if sample in samples_qc_dict:
                entry = QCEntry(
                    hb_values,
                    "Hemoglobin fraction",
                    "calculate_perc_mapping.py",
                    "No version",
                    "https://github.com/genomic-medicine-sweden/tomte/blob/dev/bin/calculate_perc_mapping.py",
                )
                samples_qc_dict[sample].append(entry)
            else:
                continue

    # Calculate slope for Rna gene body coverage
    if picard_rna_coverage:
        coverage_data = load_coverage_data(picard_rna_coverage)

        for sample, coverage_values in coverage_data.items():
            slope = calculate_slope(coverage_values)

            cov_data = {
                "genebody_cov_slope": round(slope * 1000, 4),
                "genebody_cov": list(coverage_values),
            }

            entry = QCEntry(
                cov_data,
                "Genebody coverage slope",
                "Picard",
                "No version",
                "https://broadinstitute.github.io/picard/",
            )
            samples_qc_dict[sample].append(entry)

    # Process multiqc star stats
    if multiqc_star:
        star_data = process_multiqc_star_stats(multiqc_star)
        for sample, star_stats in star_data.items():
            entry = QCEntry(
                star_stats,
                "Multiqc STAR stats",
                "MultiQC",
                "No version",
                "https://seqera.io/multiqc/",
            )
            samples_qc_dict[sample].append(entry)

    # Heterogenicity calls
    # FIXME: We'll need multiple sample paths for this one
    if hetcalls_vcfs:
        for hetcalls_vcf in hetcalls_vcfs:
            sample, het_qc = calculate_het_qcs(hetcalls_vcf)
            # for sample, het_qc in het_qcs.items():
            entry = QCEntry(
                het_qc,
                "Heterozygosity fraction",
                "BCFTools",
                "No version",
                "https://samtools.github.io/bcftools/bcftools.html",
            )
            samples_qc_dict[sample].append(entry)

    if output_file:
        write_results(samples_qc_dict, output_file, sample_id, debug)


def parse_multiqc_general_stats(multiqc_general_stats_fp: str) -> List[Tuple[str, Dict[str, str]]]:
    header = []
    row_triplets = []
    with open(multiqc_general_stats_fp, "r") as fh:
        curr_triplet = []
        for line in fh:
            line = line.strip()
            fields = line.split("\t")

            if len(header) == 0:
                header = fields
                continue

            fields = line.split("\t")
            curr_triplet.append(fields)
            if len(curr_triplet) == 3:
                row_triplets.append(curr_triplet)
                curr_triplet = []

    # samples_qcs = {}

    sample_qc_tups: List[Tuple[str, Dict[str, str]]] = []

    for row_triplet in row_triplets:
        sample, qc_dict = parse_multiqc_sample(
            header, row_triplet[0], row_triplet[1], row_triplet[2]
        )
        sample_qc_tups.append((sample, qc_dict))

    return sample_qc_tups


def calculate_het_qcs(het_call_vcf) -> Tuple[str, Dict[str, Any]]:
    calls = {}

    sample = None

    with open(het_call_vcf, "r") as fh:
        for line in fh:
            line = line.rstrip()

            if line.startswith("##"):
                continue

            if line.startswith("#"):
                header_fields = line.split("\t")
                if len(header_fields) != 10:
                    raise ValueError(f"Expected a proband-only VCF, found header: {header_fields}")
                sample = header_fields[-1]
                continue

            fields = line.split("\t")

            id_col = 2
            sample_col = 9

            rsid = fields[id_col]
            call = fields[sample_col].split(":")[0]
            calls[rsid] = call

    nbr_calls = 0
    non_calls = 0
    nbr_het_calls = 0
    for rsid, call in calls.items():
        if call == "./.":
            non_calls += 1
            continue

        nbr_calls += 1
        ref, alt = call.split("/")
        if ref == alt:
            nbr_het_calls += 1

    qc_dict = {
        "nbr_calls": nbr_calls,
        "non_calls": non_calls,
        "nbr_het_calls": nbr_het_calls,
        "het_calls": calls,
        "pct_het_calls": round((nbr_het_calls / nbr_calls) * 100,1),
    }

    if not sample:
        raise ValueError("No header line (prefixed by a single #) was found in the VCF, aborting")

    return sample, qc_dict


def build_dict(headers: List[str], fields: List[str]) -> Dict[str, str]:
    assert len(headers) == len(
        fields
    ), f"Headers and fields nbrs differs {len(headers)} {len(fields)}. Headers: {headers}, fields: {fields}"

    value_dict = {}
    for i, header in enumerate(headers):
        value = fields[i]
        value_dict[header] = value
    return value_dict


def parse_multiqc_sample(
    headers: List[str],
    shared_fields: List[str],
    fw_fields: List[str],
    rv_fields: List[str],
) -> Tuple[str, Dict[str, str]]:
    sample = shared_fields[0]

    nbr_shared = len(shared_fields) - 1

    parsed_headers = [header.split("_generalstats_")[-1] for header in headers]

    shared_headers = parsed_headers[1 : nbr_shared + 1]
    single_headers = parsed_headers[nbr_shared + 1 :]

    shared_dict = build_dict(shared_headers, shared_fields[1:])
    fw_dict = build_dict(single_headers, fw_fields[nbr_shared + 1 :])
    rv_dict = build_dict(single_headers, rv_fields[nbr_shared + 1 :])

    # Shared keys (fw+rv)
    shared_keys = {}
    shared_keys["snvs"] = "bcftools_stats-number_of_SNPs"
    shared_keys["indels"] = "bcftools_stats-number_of_indels"
    shared_keys["ts_tv"] = "bcftools_stats-tstv"
    shared_keys["median_insert_size"] = "picard_insertsizemetrics-summed_median"
    shared_keys["pct_rrna"] = "picard_rnaseqmetrics-PCT_RIBOSOMAL_BASES"
    shared_keys["pct_mrna"] = "picard_rnaseqmetrics-PCT_MRNA_BASES"
    shared_keys["n_reads"] = "star-total_reads"
    shared_keys["m_aligned"] = "star-uniquely_mapped"
    shared_keys["pct_aligned"] = "star-uniquely_mapped_percent"
    shared_keys["pct_dup"] = "fastp-pct_duplication"
    shared_keys["m_reads_after_filtering"] = "fastp-filtering_result_passed_filter_reads"
    shared_keys["gc"] = "fastp-after_filtering_gc_content"
    shared_keys["pct_pass_filter"] = "fastp-after_filtering_q30_rate"
    shared_keys["pct_adapter"] = "fastp-pct_adapter"

    # Single keys (fw or rv)
    single_keys = {}
    single_keys["dups"] = "fastqc-percent_duplicates"
    single_keys["gc"] = "fastqc-percent_gc"
    single_keys["median_read_length"] = "fastqc-median_sequence_length"
    single_keys["m_seqs"] = "fastqc-total_sequences"

    shared_out_dict = {}
    for new_key, old_key in shared_keys.items():
        value = shared_dict.get(old_key)
        if value is None:
            valid_keys = list(shared_dict.keys())
            raise ValueError(
                f"Tried getting key: {old_key}, valid keys are: {', '.join(valid_keys)})"
            )
        shared_out_dict[new_key] = value

    fw_out_dict = {}
    for new_key, old_key in single_keys.items():
        value = fw_dict.get(old_key)
        if value is None:
            valid_keys = list(fw_dict.keys())
            raise ValueError(
                f"Tried getting key: {old_key}, valid keys are: {', '.join(valid_keys)})"
            )
        fw_out_dict[f"fw_{new_key}"] = value

    rv_out_dict = {}
    for new_key, old_key in single_keys.items():
        value = rv_dict[old_key]
        if value is None:
            valid_keys = list(rv_dict.keys())
            raise ValueError(
                f"Tried getting key: {old_key}, valid keys are: {', '.join(valid_keys)})"
            )
        rv_out_dict[f"rv_{new_key}"] = value

    shared_out_dict.update(fw_out_dict)
    shared_out_dict.update(rv_out_dict)

    return sample, shared_out_dict


def load_coverage_data(file_path: str) -> Dict[str, List[float]]:
    """
    Load RNA coverage data from a file.

    Args:
        file_path (str): Path to the input file.

    Returns:
        dict: A dictionary where keys are sample names and values are numpy arrays of coverage values.

    Example:
        Example data file content:
        ```
        Sample\t0\t1\t2\t3\t4\t5
        Sample-01\t(0, 0.301104\t(1, 0.358652)\t(2, 0.451889)\t(3, 0.540978)\t(4, 0.610814)\t...
        Sample-02\t(0, 0.401104)\t(1, 0.458652)\t(2, 0.551889)\t(3, 0.640978)\t(4, 0.710814)\t...
        ```
    """
    data = {}
    with open(file_path) as file:
        lines = file.readlines()[1:]
        for line in lines:
            col = line.strip().split("\t")
            key = col[0]
            values = []
            for v in col[1:]:
                value = literal_eval(v)
                values.append(value[1])
            data[key] = values
    return data


def calculate_slope(coverage_values: List[float]) -> float:
    """
    Calculate the slope of RNA coverage values using linear regression.
    Considering only the 20-80% percentiles of the gene body.

    Args:
        coverage_values (numpy.array): Array of RNA coverage values.

    Returns:
        float: The slope of the linear regression line.
    """
    coverage_values_np = np.array(coverage_values[19:80], dtype=float)
    percentages = np.array(range(len(coverage_values_np)))
    slope, _intercept = np.polyfit(percentages, coverage_values_np, 1)
    return slope


def process_multiqc_star_stats(
    multiqc_star: str,
) -> Dict[str, Dict[str, Union[int, float]]]:
    """
    Process multiqc star stats file.

    Args:
        multiqc_star (str): Path to the input file containing stats for star multiqc.

    Returns:
        dict: A dictionary where keys are sample names and values are dictionaries with stats.

    Columns_Extracted:
        - total_reads
        - num_noncanonical_splices
        - num_splices
        - multimapped_percent
        - unmapped_mismatches_percent
        - mapped_percent

    Example:
        Example data file content:
        ```
        Sample\ttotal_reads\tavg_input_read_length\tuniquely_mapped\t....
        Sample-01\t100000\t80\t10000\t...
        Sample-02\t200000\t90\t20000\t...
        ```

    """
    columns_needed: Dict[str, type] = {
        "total_reads": int,
        "num_noncanonical_splices": int,
        "num_splices": int,
        "multimapped_percent": float,
        "unmapped_mismatches_percent": float,
        "mapped_percent": float,
    }

    new_column_names: Dict[str, str] = {
        "total_reads": "nbr_read_pairs",
        "num_noncanonical_splices": "non_canon_splice",
        "num_splices": "canon_splice",
        "multimapped_percent": "multimap_pct",
        "unmapped_mismatches_percent": "mismatch_pct",
        "mapped_percent": "mapped_pct",
    }

    with open(multiqc_star) as file:
        lines = file.readlines()
        header = lines[0].strip().split("\t")
        data = lines[1:]

    sliced_data = {}
    indices = {}
    for col in columns_needed.keys():
        if col not in header:
            raise ValueError(f"Column {col} not found in {multiqc_star}")
        indices[col] = header.index(col)

    for line in data:
        col = line.strip().split("\t")
        key = col[0]
        values = {}
        for col_name, col_type in columns_needed.items():
            value = col[indices[col_name]]
            try:
                if col_type is int and "." in value:
                    values[new_column_names[col_name]] = int(float(value))
                else:
                    values[new_column_names[col_name]] = col_type(value)
            except ValueError:
                values[new_column_names[col_name]] = float(value)
        sliced_data[key] = values

        # Calculate splice_ratio
        sliced_data[key]["splice_ratio"] = int(
            round(
                (sliced_data[key]["canon_splice"] / sliced_data[key]["non_canon_splice"]),
                0,
            )
        )

    return sliced_data


def process_hb_estimate_data(
    hb_estimate: str,
) -> Dict[str, Dict[str, Any]]:
    """
    Process HB data for a sample.
    """
    hb_json_per_sample = {}
    hb_sample_id = hb_estimate.replace("_perc_mapping.json", "")
    with open(hb_estimate) as file:
        for line in file:
            json_data = line.strip()
            hb_json_per_sample[hb_sample_id] = json.loads(json_data)
            hb_json_per_sample[hb_sample_id]["pct_reads_mapped_to_hb_genes"] = round(
                hb_json_per_sample[hb_sample_id].pop("target_genes_frac") * 100, 1
            )

    return hb_json_per_sample


def write_results(
    qc_results: Dict[str, List[QCEntry]],
    output_path: str,
    stdout_sample: Optional[str],
    debug: bool,
) -> None:
    """
    Write json blob with general statistics and rna cov values to a file.

    Args:
        data (dict): Dictionary with sample names as keys and coverage arrays as values.
        output_path (str): Path to the output file.
    """

    if stdout_sample:
        if stdout_sample not in qc_results:
            all_valid_ids = list(qc_results.keys())
            raise ValueError(
                f"Tried getting sample_id {stdout_sample}. Valid IDs are: {', '.join(all_valid_ids)}"
            )
        qc_result = qc_results[stdout_sample]
        result_dicts = [qc.get_result_dict() for qc in qc_result]
        out_json = json.dumps(result_dicts)
        if debug:
            print(out_json)
        else:
            with open(output_path, "w") as out_fh:
                out_fh.write(out_json)
    else:
        with open(output_path, "w") as output_file:
            for sample_name, qc_result in qc_results.items():
                result_dicts = [qc.get_result_dict() for qc in qc_result]
                result_json = json.dumps(result_dicts)
                output_file.write(f"{sample_name}\t{json.dumps(result_json)}\n")
        print(f"Results for {len(qc_results)} samples written to {output_path}.")


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument(
        "--multiqc_general_stats", help="MultiQC file: multiqc_general_stats.txt", required=True
    )
    parser.add_argument(
        "--sample_id",
        default=None,
        help="Limit output to one sample",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Run in debug mode, printing to STDOUT instead of to file, if limiting output to a single sample",
    )
    parser.add_argument(
        "--picard_rna_coverage", help="Path to the input file containing RNA coverage data."
    )
    parser.add_argument(
        "--multiqc_star",
        help="Path to the input file containing stats for star multiqc.",
    )
    parser.add_argument(
        "--hb_estimate",
        help="Path to the input file containing merged HB estimate data.",
    )
    parser.add_argument(
        "--hetcalls_vcfs", help="VCF containing heterozygosity calls for selected SNPs", nargs="*"
    )
    parser.add_argument("--output_file", required=True)
    parser.add_argument(
        "--version",
        action="version",
        version=VERSION,
        help="Show the program's version and exit",
    )

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_arguments()

    main(
        args.multiqc_general_stats,
        args.picard_rna_coverage,
        args.multiqc_star,
        args.hb_estimate,
        args.hetcalls_vcfs,
        args.output_file,
        args.sample_id,
        args.debug,
    )
