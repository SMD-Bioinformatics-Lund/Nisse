process GENMOD_MODELS {

    tag "${meta.sample}"
    label "process_low"
    container "${params.containers.genmod}"

    input:
    tuple val(meta), path(vcf)
    path(family_all_cases)

    output:
    tuple val(meta), path("${meta.sample}_models.vcf"), emit: vcf
    tuple val(meta), path("*_versions.yml"), emit: versions

    script:
    """
    grep -E "^#|^${meta.case}\\s" "${family_all_cases}" > "${meta.case}.ped"
    genmod \\
        models \\
        --whole_gene \\
        --processes "${task.cpus}" \\
        --family_file "${meta.case}.ped" \\
        --outfile "${meta.sample}_models.vcf" \\
        "${vcf}"
    
    ${genmodscore_version(task)}
    """

    stub:
    """
    touch "${meta.sample}_models.vcf"
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
