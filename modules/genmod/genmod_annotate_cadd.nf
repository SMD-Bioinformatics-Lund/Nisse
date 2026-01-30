process GENMOD_ANNOTATE_CADD {

    tag "${meta.sample}"
    label "process_low"
    container "${params.containers.genmod}"

    input:
    tuple val(meta), path(vcf)
    tuple path(cadd), val(cadd_tbi)

    output:
    tuple val(meta), path("*_cadd_annotate.vcf"), emit: vcf
    path("*_versions.yml"), emit: versions

    script:
    """
    genmod \\
        annotate \\
        --cadd-file "${cadd}" \\
        --outfile ${meta.sample}_cadd_annotate.vcf \\
        ${vcf}

    ${genmodscore_version(task)}
    """

    stub:
    """
    touch ${meta.sample}_cadd_annotate.vcf
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
