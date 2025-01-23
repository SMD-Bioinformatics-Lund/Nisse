process VCF_COMPLETION {

    tag "${meta.sample}"
	label "process_low"
	container "${params.containers.base}"

	input:
		tuple val(meta), path(vcf)

	output:
		tuple val(meta), path("${meta.sample}_scored.vcf.gz"), path("${meta.sample}_scored.vcf.gz.tbi"), emit: vcf
		tuple val(meta), path("*versions.yml"), emit: versions

	script:
		"""
		sed 's/^MT/M/' -i $vcf
		sed 's/ID=MT,length/ID=M,length/' -i $vcf
		bgzip -@ ${task.cpus} $vcf -f
		tabix ${vcf}.gz -f

		${vcf_completion_version(task)}
		"""

	stub:
		"""
		touch "${meta.sample}_scored.vcf.gz"
		touch "${meta.sample}_scored.vcf.gz.tbi"

		${vcf_completion_version(task)}
		"""
}
def vcf_completion_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    tabix: \$(echo \$(tabix --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
	END_VERSIONS	
	"""
}