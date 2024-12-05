#!/usr/bin/env python3

import argparse
import json
from typing import Any, Optional, Dict, List, Tuple
import numpy as np
from ast import literal_eval


def main(
    multiqc_general_stats: str, picard_rna_coverage: str, output_file: str, sample_id: Optional[str]
):

    in_fp = multiqc_general_stats

    header = []
    batches = []
    with open(in_fp, "r") as in_fh:
        curr_batch = []
        for line in in_fh:
            line = line.strip()
            fields = line.split("\t")

            if len(header) == 0:
                header = fields
                continue

            fields = line.split("\t")
            curr_batch.append(fields)
            if len(curr_batch) == 3:
                batches.append(curr_batch)
                curr_batch = []

    for batch in batches:
        sample, combined_dict = make_combined_dict(header, batch[0], batch[1], batch[2])

    results = {}

    for batch in batches:
        sample, combined_dict = make_combined_dict(header, batch[0], batch[1], batch[2])
        results[sample] = combined_dict

    # Calculate slope for Rna gene body coverage
    if picard_rna_coverage:
        coverage_data = load_coverage_data(picard_rna_coverage)

        for sample, coverage_values in coverage_data.items():
            slope = calculate_slope(coverage_values)
            results[sample]["genebody_cov_slope"] = slope
            results[sample]["genebody_cov"] = list(coverage_values)

    if output_file:
        write_results(results, output_file, sample_id)


def write_results(data: Dict[str, Any], output_path: str, sample_id: Optional[str]) -> None:
    """
    Write json blob with general statistics and rna cov values to a file.

    Args:
        data (dict): Dictionary with sample names as keys and coverage arrays as values.
        output_path (str): Path to the output file.
    """
    with open(output_path, "w") as output_file:
        if sample_id is not None:
            print(f"Extracting only sample: {sample_id}")
            sample_only_data = data.get(sample_id)
            if sample_only_data is None:
                all_valid_ids = list(data.keys())
                raise ValueError(f"Tried getting sample_id {sample_id}. Valid IDs are: {', '.join(all_valid_ids)}")
            output_file.write(json.dumps(sample_only_data))
        else:
            for sample_name, values in data.items():
                output_file.write(f"{sample_name}\t{json.dumps(values)}\n")
    print(f"Results written to {output_path}.")


def build_dict(headers: List[str], fields: List[str]) -> Dict[str, str]:
    assert len(headers) == len(
        fields
    ), f"Headers and fields nbrs differs {len(headers)} {len(fields)}"

    value_dict = {}
    for i, header in enumerate(headers):
        value = fields[i]
        value_dict[header] = value
    return value_dict


def make_combined_dict(
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
    shared_keys["insert_size"] = "picard_insertsizemetrics-summed_median"
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
            raise ValueError(f"Tried getting key: {old_key}, valid keys are: {', '.join(valid_keys)})")
        shared_out_dict[new_key] = value

    fw_out_dict = {}
    for new_key, old_key in single_keys.items():
        value = fw_dict.get(old_key)
        if value is None:
            valid_keys = list(fw_dict.keys())
            raise ValueError(f"Tried getting key: {old_key}, valid keys are: {', '.join(valid_keys)})")
        fw_out_dict[f"fw_{new_key}"] = value

    rv_out_dict = {}
    for new_key, old_key in single_keys.items():
        value = rv_dict[old_key]
        if value is None:
            valid_keys = list(rv_dict.keys())
            raise ValueError(f"Tried getting key: {old_key}, valid keys are: {', '.join(valid_keys)})")
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


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--multiqc_general_stats", required=True)
    parser.add_argument(
        "--sample_id", default=None, help="Limit output to one sample and print to STDOUT"
    )
    parser.add_argument(
        "--picard_rna_coverage",
        help="Path to the input file containing RNA coverage data.",
    )
    parser.add_argument("--output_file", required=True)

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_arguments()
    main(args.multiqc_general_stats, args.picard_rna_coverage, args.output_file, args.sample_id)