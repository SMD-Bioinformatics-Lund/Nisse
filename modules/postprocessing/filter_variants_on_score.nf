process FILTER_VARIANTS_ON_SCORE {

    tag "${meta.sample}"
	label "process_small"
	container "${params.containers.base}"

    input:
    tuple val(meta), path(scored_vcf), path(scored_vcf_tbi)
    val(threshold)

    output:
    path "${meta.sample}_out.txt"

    script:
    """
    bash hello.sh ${scored_vcf} > "${meta.sample}_out.txt"
    """

    stub:
    """
    touch ${meta.sample}_out.txt
    """
}