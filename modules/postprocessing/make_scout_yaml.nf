process MAKE_SCOUT_YAML {

    tag "${meta.sample}"
    label "process_low"
    container "${params.containers.base}"

    input:
    tuple val(meta), path(nisse_parsed_fraser), path(nisse_parsed_outrider), path(nisse_parsed_vcf), path(nisse_junctions), path(nisse_junctions_tbi), path(vcf_tbi), path(cram), path(cram_crai), path(bigwig)
    val tomte_results_dir
    val nisse_output_dir
    val phenotype
    val tissue

    output:
    path "${meta.sample}_scout.yaml", emit: yaml

    script:
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
        --bam_path "${tomte_results_dir}/alignment/${cram}" \
        --splice_junctions "${nisse_output_dir}/junction/${nisse_junctions}" \
        --rna_bigwig "${tomte_results_dir}/ucsc/${bigwig}" > ${meta.sample}_scout.yaml
    """

    stub:
    """
    touch "${meta.sample}_scout.yaml"
    """
}
