process MAKE_SCOUT_YAML {

    tag "${meta.sample}"
    label "process_low"
    container "${params.containers.base}"

    input:
    tuple val(meta), path(fraser), path(outrider), path(vcf), path(vcf_tbi), path(junctions), path(junctions_tbi), path(cram), path(cram_crai), path(bigwig), path(peddy)
    val tomte_results_dir
    val nisse_output_dir
    val phenotype
    val tissue

    output:
    path "${meta.sample}_scout.yaml", emit: yaml

    script:
    def sex = "${meta.sex}" == "M" ? "male" : "${meta.sex}" == "F" ? "female" : "unknown"
    // FIXME: produce_yaml.py should take the Peddy sex check
    """
    produce_yaml.py \
        --sample_id "${meta.sample}" \
        --fraser "${nisse_output_dir}/drop/${fraser}" \
        --outrider "${nisse_output_dir}/drop/${outrider}" \
        --vcf "${nisse_output_dir}/vcf/${vcf}" \
        --sex "${sex}" \
        --phenotype "${phenotype}" \
        --tissue "${tissue}" \
        --bam_path "${tomte_results_dir}/alignment/${cram}" \
        --splice_junctions "${nisse_output_dir}/junction/${junctions}" \
        --rna_bigwig "${tomte_results_dir}/ucsc/${bigwig}" > ${meta.sample}_scout.yaml
    """

    stub:
    """
    touch "${meta.sample}_scout.yaml"
    """
}
