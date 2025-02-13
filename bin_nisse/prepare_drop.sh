#!/bin/bash

# What is this script even doing
# Should it be Python

if [[ $# -ne 6 ]]; then
    echo "Usage: $0 <in path> <out dir> <fdr col> <fdr cutoff> <hgncSymbol col> <sample col>"
    exit 1
fi

in_path=$1
out_dir=$2
fdr_col=$3
fdr_cutoff=$4
hgnc_symbol_col=$5
sample_col=$6

mkdir -p "${out_dir}"

echo "Initial number of entries, counting file: \"${in_path}\""
tail -n +2 "${in_path}" | wc -l

after_cutoff="${out_dir}/cutoff.tsv"

awk -v fdr_col="${fdr_col}" -v fdr_cutoff="${fdr_cutoff}" 'NR == 1 || $fdr_col <= fdr_cutoff {print $0}' "${in_path}" > "${after_cutoff}"
echo "$(cat "${after_cutoff}" | wc -l) lines written to ${after_cutoff}"

after_ens="${out_dir}/after_ens.tsv"
cat "${after_cutoff}" |
    awk -v hgnc_col="${hgnc_symbol_col}" 'NR == 1 || $hgnc_col !~ /^ENS/ {print $0}' > "${after_ens}"
echo "$(wc -l "${after_ens}") lines written to ${after_ens}"

after_hgnc_id="${out_dir}/after_hgnc_id.tsv"
python3 add_col.py \
    --input "${after_ens}" \
    --map active_hgncid_map.tsv \
    --target hgncSymbol \
    --new_header hgnc_id \
    --output "${after_hgnc_id}"
echo "$(cat "${after_hgnc_id}" | wc -l) lines written to ${after_hgnc_id}"

# Split on sample
samples_out_dir="${out_dir}/samples"
mkdir -p "${samples_out_dir}"
tail -n +2 "${after_hgnc_id}" | cut -f "${sample_col}" | sort | uniq | while read -r sample; do
    sample_out_fp="${samples_out_dir}/${sample}.tsv"
    #echo "Writing entries to ${sample_out_fp}"
    awk \
        -v sample_col="${sample_col}" \
        -v sample="${sample}" \
        'NR == 1 || $sample_col == sample' "${after_hgnc_id}" > "${sample_out_fp}"
done
echo "$(find "${samples_out_dir}" | wc -l) sample TSVs written to ${samples_out_dir}"













