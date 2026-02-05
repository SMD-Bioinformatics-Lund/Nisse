process BGZIP_TABIX {

    tag "${meta.sample}"
    label "process_low"
    container "${params.containers.bcftools}"

    input:
    tuple val(meta), path(bed_or_vcf)

    output:
    tuple val(meta), path("${bed_or_vcf}.gz"), path("${bed_or_vcf}.gz.tbi"), emit: vcf_tbi, optional: true
    tuple val(meta), path("${bed_or_vcf}.gz"), path("${bed_or_vcf}.gz.tbi"), emit: bed_tbi, optional: true
    path("*_versions.yml"), emit: versions

    script:
    """
    bgzip "${bed_or_vcf}"
    tabix "${bed_or_vcf}.gz"

    ${versions(task)}
    """

    stub:
    """
    touch "${bed_or_vcf}.gz"
    touch "${bed_or_vcf}.gz.tbi"

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
