process MODIFY_VCF {
	
	tag "${meta.sample}"
    label 'process_small'
	container "${params.containers.base}"

	input:
		tuple val(meta), file(vcf)

	output:
		tuple val(meta), file("${meta.sample}.mod.vcf"), emit: vcf

	script:
		"""
		modify_vcf_scout.pl $vcf > ${meta.sample}.mod.vcf
		"""

	stub:
		"""
		touch "${meta.sample}.mod.vcf"
		"""
} 