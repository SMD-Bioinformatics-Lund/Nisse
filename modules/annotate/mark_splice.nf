process MARK_SPLICE {

	tag "${meta.sample}"
	// container "${params.containers.base}"
    label 'process_small'


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