process EXTRACT_INDELS_FOR_CADD {
	label "process_low"

	input:
		tuple val(meta), file(vcf), file(vcf_tbi)
	
	output:
		tuple val(meta), file("${meta.case}.only_indels.vcf"), emit: vcf
		tuple val(meta), file("*versions.yml"), emit: versions

	script:
		"""
		bcftools view $vcf -V snps -o ${meta.case}.only_indels.vcf
		${extract_indels_for_cadd_version(task)}
		"""

	stub:
		"""
		touch "${meta.case}.only_indels.vcf"
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