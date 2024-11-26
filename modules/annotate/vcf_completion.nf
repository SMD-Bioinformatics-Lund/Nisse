process VCF_COMPLETION {
	cpus 16
	publishDir "${OUTDIR}/vcf", mode: 'copy', overwrite: 'true', pattern: '*.vcf.gz*'
	tag "$meta.id"
	time '1h'
	memory '5 GB'

	input:
		set meta, type, file(vcf) from scored_vcf

	output:
		set meta, type, file("${meta.id}.scored.vcf.gz"), file("${meta.id}.scored.vcf.gz.tbi")
		set meta, file("*versions.yml")

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
		touch "${meta.id}.scored.vcf.gz"
		touch "${meta.id}.scored.vcf.gz.tbi"

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