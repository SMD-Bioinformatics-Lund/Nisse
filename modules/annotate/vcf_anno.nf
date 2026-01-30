process VCF_ANNO {

	tag "${meta.sample}"
    label 'process_low'
	errorStrategy 'retry'
	maxErrors 5
	container "${params.containers.base}"

	input:
		tuple val(meta), path(vcf)
		val(vep_params)

	output:
		tuple val(meta), path("${meta.sample}.clinvar.loqusdb.gene.vcf"), emit: vcf
		path("*versions.yml"), emit: versions

	script:
		"""
		vcfanno_linux64 -lua "${vep_params.VCFANNO_LUA}" "${vep_params.VCFANNO}" "${vcf}" > ${meta.sample}.clinvar.loqusdb.gene.vcf
		${vcfanno_version(task)}
		"""

	stub:
		"""
		touch "${meta.sample}.clinvar.loqusdb.gene.vcf"
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