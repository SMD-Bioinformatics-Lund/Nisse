process extract_indels_for_cadd {
	cpus 2
	tag "$group"
	memory '1 GB'
	time '1h'

	input:
		set group, id, file(vcf)
	
	output:
		set group, file("${group}.only_indels.vcf")
		set group, file("*versions.yml")

	script:
		"""
		bcftools view $vcf -V snps -o ${group}.only_indels.vcf
		${extract_indels_for_cadd_version(task)}
		"""

	stub:
		"""
		touch "${group}.only_indels.vcf"
		${extract_indels_for_cadd_version(task)}
		""" 
}
def extract_indels_for_cadd_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    bcftools: \$(echo \$(bcftools --version 2>&1) | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
	END_VERSIONS	
	"""
}