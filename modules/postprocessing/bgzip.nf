process BGZIP {

    tag "${meta.sample}"
    label "process_low"
    container "${params.containers.bcftools}"

    input:
    tuple val(meta), path(bed_or_vcf)

    output:
    tuple val(meta), path("${bed_or_vcf}.gz"), emit: vcf, optional: true
    tuple val(meta), path("${bed_or_vcf}.gz"), emit: bed, optional: true
    path("*_versions.yml"), emit: versions

    script:
    """
    bgzip "${bed_or_vcf}"

    ${versions(task)}
    """

    stub:
    """
    touch "${bed_or_vcf}.gz"

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
