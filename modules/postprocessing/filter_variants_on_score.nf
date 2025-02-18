process FILTER_VARIANTS_ON_SCORE {

    tag "${meta.sample}"
    label "process_low"
    container "${params.containers.base}"

    input:
    tuple val(meta), path(scored_vcf), path(scored_vcf_tbi)
    val score_threshold

    output:
    tuple val(meta), path("${meta.sample}_filtered.vcf"), emit: vcf

    script:
    """
    filter_score.py "${scored_vcf}" "${score_threshold}" > "${meta.sample}_filtered.vcf"
    """

    stub:
    """
    touch "${meta.sample}_filtered.vcf"
    """
}
