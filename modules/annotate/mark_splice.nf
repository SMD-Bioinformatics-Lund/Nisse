process MARK_SPLICE {
	cpus 2
	tag "$meta.id"
	memory '1 GB'
	time '1h'

	input:
		set meta, file(vcf)

	output:
		set meta, file("${meta.id}.marksplice.vcf")

	script:
		"""
		/opt/bin/mark_spliceindels.pl $vcf > ${meta.id}.marksplice.vcf
		"""

	stub:
		"""
		touch "${meta.id}.marksplice.vcf"
		"""
}