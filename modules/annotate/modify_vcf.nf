process MODIFY_VCF {
	cpus 2
	tag "$meta.id"
	memory '1 GB'
	time '1h'

	input:
		set meta, file(vcf) from vcfanno_vcf

	output:
		set meta, file("${meta.id}.mod.vcf") into mod_vcf

	script:
		"""
		modify_vcf_scout.pl $vcf > ${meta.id}.mod.vcf
		"""

	stub:
		"""
		touch "${meta.id}.mod.vcf"
		"""
} 