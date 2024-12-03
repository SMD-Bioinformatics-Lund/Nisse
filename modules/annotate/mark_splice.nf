process MARK_SPLICE {

	tag "${meta.id}"
	container "${params.containers.base}"
    label 'process_small'


	input:
		tuple val(meta), file(vcf)

	output:
		tuple val(meta), file("${meta.id}.marksplice.vcf")

	script:
		"""
		/opt/bin/mark_spliceindels.pl $vcf > ${meta.id}.marksplice.vcf
		"""

	stub:
		"""
		touch "${meta.id}.marksplice.vcf"
		"""
}