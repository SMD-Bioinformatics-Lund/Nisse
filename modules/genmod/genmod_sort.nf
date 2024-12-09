process GENMOD_SORT {

    tag "${meta.sample}"
	label "process_low"
	container "${params.containers.genmod}"

    input:
        tuple val(meta), path(vcf)

    output:
        tuple val(meta), path("${meta.sample}_scored.vcf"), emit: vcf
        tuple val(meta), path("*_versions.yml"), emit: versions

    script:
    """
    genmod \\
        sort \\
        -p \\
        -f "${meta.sample}" \\
        -o "${meta.sample}_scored.vcf" \\
        "${vcf}"
    
    ${genmodscore_version(task)}
    """

    stub:
    """
    touch "${meta.sample}_scored.vcf"
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