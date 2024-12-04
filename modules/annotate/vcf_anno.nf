process VCF_ANNO {

	tag "${meta.id}"
	container "${params.containers.base}"
    label 'process_small'
	errorStrategy 'retry'
	maxErrors 5

	input:
		tuple val(meta), file(vcf)

	output:
		tuple val(meta), file("${meta.id}.clinvar.loqusdb.gene.vcf"), emit: vcf
		tuple val(meta), file("*versions.yml"), emit: versions

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