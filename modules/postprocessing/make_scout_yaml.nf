process MAKE_SCOUT_YAML {

    tag "${meta.sample}"
    label "process_low"
    container "${params.containers.base}"

    input:
    tuple val(meta), path(nisse_parsed_fraser), path(nisse_parsed_outrider), path(nisse_parsed_vcf), path(vcf_tbi)
    val tomte_results_dir
    val nisse_output_dir

    output:
    path "${meta.sample}_scout.yaml", emit: yaml

    script:
    // FIXME: Pull these out from the process
    def phenotype = "affected"
    def tissue = "blood"
    def sex = "${meta.sex}" == "M" ? "male" : "${meta.sex}" == "F" ? "female" : "unknown"
    """
    produce_yaml.py \
        --sample_id "${meta.sample}" \
        --fraser "${nisse_output_dir}/drop/${nisse_parsed_fraser}" \
        --outrider "${nisse_output_dir}/drop/${nisse_parsed_outrider}" \
        --vcf "${nisse_output_dir}/vcf/${nisse_parsed_vcf}" \
        --sex "${sex}" \
        --phenotype "${phenotype}" \
        --tissue "${tissue}" \
        --bam_path "${tomte_results_dir}/alignment/${meta.sample}.cram" \
        --splice_junctions "${nisse_output_dir}/junction/${meta.sample}_junction.bed.gz" \
        --rna_bigwig "${tomte_results_dir}/ucsc/${meta.sample}.bw" > ${meta.sample}_scout.yaml
    """

    stub:
    """
    touch "${meta.sample}_scout.yaml"
    """
}
