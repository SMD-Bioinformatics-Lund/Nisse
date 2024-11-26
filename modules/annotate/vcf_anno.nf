process VCFANNO {
	memory '1GB'
	time '1h'
	errorStrategy 'retry'
	maxErrors 5
	cpus 2

	input:
		set meta, file(vcf) from vep

	output:
		set meta, file("${meta.id}.clinvar.loqusdb.gene.vcf") into vcfanno_vcf
		set meta, file("*versions.yml") into ch_vcfanno_versions

	script:
		"""
		vcfanno_linux64 -lua $params.VCFANNO_LUA $params.vcfanno $vcf > ${meta.id}.clinvar.loqusdb.gene.vcf
		${vcfanno_version(task)}
		"""

	stub:
		"""
		touch "${meta.id}.clinvar.loqusdb.gene.vcf"
		${vcfanno_version(task)}
		"""
}
def vcfanno_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    vcfanno: \$(echo \$(vcfanno_linux64 2>&1 | grep version | cut -f3 -d' ')  )
	END_VERSIONS	
	"""
}