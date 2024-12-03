process MODIFY_VCF {
	
	tag "${meta.id}"
    label 'process_small'

	input:
		tuple val(meta), file(vcf)

	output:
		tuple val(meta), file("${meta.id}.mod.vcf")

	script:
		"""
		modify_vcf_scout.pl $vcf > ${meta.id}.mod.vcf
		"""

	stub:
		"""
		touch "${meta.id}.mod.vcf"
		"""
} 