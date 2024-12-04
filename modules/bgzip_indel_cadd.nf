process BGZIP_INDEL_CADD {

    tag "${meta.sample}"
	label "process_medium"
	container "${params.container_bcftools}"

	input:
		tuple val(meta), file(cadd_scores)
	
	output:
		tuple val(meta), file("${meta.sample}.cadd.gz"), file("${meta.sample}.cadd.gz.tbi"), emit: cadd
		tuple val(meta), file("*versions.yml"), emit: versions
	
	script:
		"""
		gunzip -c ${cadd_scores} > "${meta.sample}.cadd"
		bgzip -@ ${task.cpus} "${meta.sample}.cadd"
		tabix -p vcf "${meta.sample}.cadd.gz"
		${bgzip_indel_cadd_version(task)}
		"""
	
	stub:
		"""
		touch "${meta.sample}.cadd.gz"
		touch "${meta.sample}.cadd.gz.tbi"
		${bgzip_indel_cadd_version(task)}
		"""
}
def bgzip_indel_cadd_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    tabix: \$(echo \$(tabix --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
	    bcftools: \$(echo \$(bcftools --version 2>&1) | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
	END_VERSIONS
	"""
}
