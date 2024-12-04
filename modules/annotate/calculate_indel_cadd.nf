process CALCULATE_INDEL_CADD {

    tag "${meta.sample}"
	label "process_medium"
	// container "${params.cadd}"

	input:
		tuple val(meta), file(vcf)

	output:
		tuple val(meta), file("${meta.sample}.indel_cadd.gz"), emit: vcf
		tuple val(meta), file("*versions.yml"), emit: versions

	script:
		"""
		/CADD-scripts/CADD.sh -c ${task.cpus} -g GRCh38 -o ${meta.sample}.indel_cadd.gz $vcf
		${calculate_indel_cadd_version(task)}
		"""

	stub:
		"""
		touch "${meta.sample}.indel_cadd.gz"
		${calculate_indel_cadd_version(task)}
		"""
}
def calculate_indel_cadd_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    cadd: \$(echo \$(/CADD-scripts/CADD.sh -v 2>&1) | sed -e "s/^.*CADD-v// ; s/ (c).*//")
	END_VERSIONS	
	"""
}