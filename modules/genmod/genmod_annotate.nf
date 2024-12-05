process GENMOD_ANNOTATE {

    tag "${meta.sample}"
	label "process_low"
	container "${params.containers.genmod}"

    input:
        tuple val(meta), path(vcf)

    output:
        tuple val(meta), path("*_annotate.vcf"), emit: vcf
        path("*_versions.yml"), emit: versions

    script:
    """
    genmod \\
        annotate \\
        --annotate_regions \\
        --genome-build 38 \\
        --outfile ${meta.sample}_annotate.vcf \\
        ${vcf}

    ${genmodscore_version(task)}
    """

    stub:
    """
    touch ${meta.sample}_annotate.vcf
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