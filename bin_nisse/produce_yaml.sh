#!/bin/bash

if [[ $# -ne 7 ]]; then
    echo "Usage: $0 <template> <sample ID> <results dir> <parsed DROP> <parsed VCF> <out fp>"
    exit 1
fi

template=$1
sample_id=$2
results_dir=$3
fraser_path=$4
outrider_path=$5
vcf_path=$6
out_fp=$7

cat "${template}" > "${out_fp}"
sed -i "s|<FAMILY_ID>|${sample_id}|" "${out_fp}"
sed -i "s|<FAMILY_NAME>|${sample_id}|" "${out_fp}"
sed -i "s|<SAMPLE_ID>|${sample_id}|" "${out_fp}"
sed -i "s|<SAMPLE_NAME>|${sample_id}|" "${out_fp}"

cram_path="${results_dir}/alignment/${sample_id}.cram"
if [ -f "${cram_path}" ]; then
    sed -i "s|<BAM>|${cram_path}|" "${out_fp}"
    sed -i "s|<RNA_BAM>|${cram_path}|" "${out_fp}"
else
    echo "No CRAM file found at path: ${cram_path}"
    exit 1
    #sed -i "/<BAM>/d" ${out_fp}
fi

bigwig_path="${results_dir}/ucsc/${sample_id}.bw"
[ -f "${bigwig_path}" ] || { echo "Error: ${bigwig_path} is not a valid file"; exit 1; }
sed -i "s|<RNA_BIGWIG>|${bigwig_path}|" "${out_fp}"

splice_path="${results_dir}/junction/${sample_id}_junction.bed.gz"
[ -f "${splice_path}" ] || { echo "Error: ${splice_path} is not a valid file"; exit 1; }
sed -i "s|<RNA_SPLICE>|${splice_path}|" "${out_fp}"

# variants_path="${results_dir}/call_variants/${sample_id}.filtered.vcf.gz"
variants_path="${vcf_path}"
[ -f "${variants_path}" ] || { echo "Error: ${variants_path} is not a valid file"; exit 1; }
sed -i "s|<VCF>|${variants_path}|" "${out_fp}"

# fraser_path="${parsed_drop}/fraser/samples/${sample_id}.tsv"
fraser_exists=true
if [ -f "${fraser_path}" ]; then
    sed -i "s|<FRASER>|${fraser_path}|" "${out_fp}"
else
    sed -i "/<FRASER>/d" "${out_fp}"
    fraser_exists=false
    echo "FRASER not found ${fraser_path}, skipping"
fi

# outrider_path="${parsed_drop}/outrider/samples/${sample_id}.tsv"
outrider_exists=true
if [ -f "${outrider_path}" ]; then
    sed -i "s|<OUTRIDER>|${outrider_path}|" "${out_fp}"
else
    sed -i "/<OUTRIDER>/d" "${out_fp}"
    outrider_exists=false
    echo "OUTRIDER not found ${outrider_path}, skipping"
fi

if [ "${outrider_exists}" = false ] && [ "${fraser_exists}" = false ]; then
    sed -i "/omics_files:/d" "${out_fp}"
    echo "hit"
fi



