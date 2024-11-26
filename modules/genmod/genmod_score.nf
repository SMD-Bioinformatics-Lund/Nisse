process GENMOD_SCORE {

    input:
        tuple val(prefix), path(vcf)
        path(score_config)
        path(ped)

    output:
        tuple val(prefix), path("*_score.vcf"), emit: vcf

    script:
    """
    genmod \\
        score \\
        --score_config ${score_config} \\
        --family_file ${ped} \\
        --rank_results \\
        --outfile ${prefix}_score.vcf \\
        ${vcf}
    """

    stub:
    """
    touch ${prefix}_score.vcf
    """
}