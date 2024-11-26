process CALCULATE_INDEL_CADD {
	cpus 2
	container = '/fs1/resources/containers/cadd_v1.6.sif'
	tag "${meta.id}"
	memory '20 GB'
	time '3h'

	input:
		set meta, file(vcf)

	output:
		set meta, file("${meta.id}.indel_cadd.gz")
		set meta, file("*versions.yml")

	script:
		"""
		/CADD-scripts/CADD.sh -c ${task.cpus} -g GRCh38 -o ${meta.id}.indel_cadd.gz $vcf
		${calculate_indel_cadd_version(task)}
		"""

	stub:
		"""
		touch "${meta.id}.indel_cadd.gz"
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