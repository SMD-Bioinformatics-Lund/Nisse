process VCF_COMPLETION {

    tag "${meta.sample}"
	label "process_small"
	container "${params.base}"

	input:
		tuple val(meta), file(vcf)

	output:
		tuple val(meta), file("${meta.sample}.scored.vcf.gz"), file("${meta.sample}.scored.vcf.gz.tbi"), emit: vcf
		tuple val(meta), file("*versions.yml"), emit: versions

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
		touch "${meta.sample}.scored.vcf.gz"
		touch "${meta.sample}.scored.vcf.gz.tbi"

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