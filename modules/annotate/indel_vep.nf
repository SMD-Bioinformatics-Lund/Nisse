process INDEL_VEP {
	// container "${params.container_vep}"
	label "process_medium"

	input:
		tuple val(meta), file(vcf)

	output:
		tuple val(meta), file("${meta.id}.only_indels.vep.filtered.vcf"), emit: vcf
		tuple val(meta), file("*versions.yml"), emit: versions

	script:
		"""
		vep \\
			-i "${vcf}" \\
			-o "${meta.group}.only_indels.vep.vcf" \\
			--offline \\
			--cache \\
			--merged \\
			--vcf \\
			--synonyms "${params.vep.VEP_SYNONYMS}" \\
			--fasta "${params.vep.VEP_FASTA}" \\
			-custom "${params.vep.GNOMAD_GENOMES},gnomADg,vcf,exact,0,AF" \\
			-custom "${params.vep.GNOMAD_MT},gnomAD_mt,vcf,exact,0,AF_hom,AF_het" \\
			--dir_cache "${params.vep.VEP_CACHE}" \\
			--force_overwrite \\
			--no_stats \\
			--fork "${task.cpus}"
		filter_indels.pl "${meta.id}.only_indels.vep.vcf" > "${meta.id}.only_indels.vep.filtered.vcf"
		${indel_vep_version(task)}
		"""

	stub:
		"""
		touch "${meta.id}.only_indels.vep.filtered.vcf"
		${indel_vep_version(task)}
		"""
}
def indel_vep_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    vep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : // ; s/ .*\$//')
	END_VERSIONS	
	"""
}