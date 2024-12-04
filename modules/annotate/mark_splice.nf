process MARK_SPLICE {

	tag "${meta.sample}"
    label 'process_small'
	container "${params.containers.base}"


	input:
		tuple val(meta), file(vcf)

	output:
		tuple val(meta), file("${meta.sample}.marksplice.vcf"), emit: vcf

	script:
		"""
		/opt/bin/mark_spliceindels.pl $vcf > ${meta.sample}.marksplice.vcf
		"""

	stub:
		"""
		touch "${meta.sample}.marksplice.vcf"
		"""
}