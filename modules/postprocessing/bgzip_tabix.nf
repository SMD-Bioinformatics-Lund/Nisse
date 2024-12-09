process BGZIP_TABIX {

    tag "${meta.sample}"
    label "process_low"
    container "${params.containers.bcftools}"

    input:
    tuple val(meta), path(junction_bed)

    output:
    tuple val(meta), path("${junction_bed}.gz"), path("${junction_bed}.gz.tbi"), emit: gz
    tuple val(meta), path("*_versions.yml"), emit: versions

    script:
    """
    bgzip "${junction_bed}"
    tabix "${junction_bed}.gz"

    ${versions(task)}
    """

    stub:
    """
    touch "${junction_bed}.bed.gz"
    touch "${junction_bed}.bed.gz.tbi"

    ${versions(task)}
    """
}
def versions(task) {
    """
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    tabix: \$(echo \$(tabix --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
	END_VERSIONS	
	"""
}
