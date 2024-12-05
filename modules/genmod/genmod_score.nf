process GENMOD_SCORE {

    tag "${meta.sample}"
	label "process_low"
	container "${params.containers.genmod}"

    input:
        tuple val(meta), path(vcf)
        tuple val(meta2), path(ped)
        path(score_config)

    output:
        tuple val(meta), path("*_score.vcf"), emit: vcf

    script:
    """
    genmod \\
        score \\
        --score_config ${score_config} \\
        --family_file ${ped} \\
        --rank_results \\
        --outfile ${meta.sample}_score.vcf \\
        ${vcf}
    """

    stub:
    """
    touch ${meta.sample}_score.vcf
    """
}