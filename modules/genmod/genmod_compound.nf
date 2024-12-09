process GENMOD_COMPOUND {

    tag "${meta.sample}"
	label "process_medium"
	container "${params.containers.genmod}"

    input:
        tuple val(meta), path(vcf)

    output:
        tuple val(meta), path("*_compound.vcf"), emit: vcf
        tuple val(meta), path("*_versions.yml"), emit: versions

    script:
    """
    genmod \\
        compound \\
        --processes ${task.cpus} \\
        --outfile ${meta.sample}_compound.vcf \\
        --penalty 0 \\
        ${vcf}
    ${genmodscore_version(task)}
    """

    stub:
    """
    touch ${meta.sample}_compound.vcf
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