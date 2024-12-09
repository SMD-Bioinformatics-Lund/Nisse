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
        tuple val(meta), path("*_versions.yml"), emit: versions

    script:
    """
    genmod \\
        score \\
        --score_config ${score_config} \\
        --family_file ${ped} \\
        --rank_results \\
        --outfile ${meta.sample}_score.vcf \\
        ${vcf}
    ${genmodscore_version(task)}
    """

    stub:
    """
    touch ${meta.sample}_score.vcf
    ${genmodscore_version(task)}
    """
}
def genmodscore_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    genmod: \$(echo \$(genmod --version 2>&1) | sed -e "s/^.*genmod version: //")
	END_VERSIONS	
	"""
}