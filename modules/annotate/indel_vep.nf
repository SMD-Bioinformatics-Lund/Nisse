process INDEL_VEP {

    tag "${meta.sample}"
	label "process_medium"
	container "${params.containers.vep}"

	input:
		tuple val(meta), path(vcf)
		val(vep_params)

	output:
		tuple val(meta), path("${meta.sample}.only_indels.vep.filtered.vcf"), emit: vcf
		path("*versions.yml"), emit: versions

	script:
		"""
		vep \\
			-i "${vcf}" \\
			-o "${meta.sample}.only_indels.vep.vcf" \\
			--offline \\
			--cache \\
			--merged \\
			--vcf \\
			--synonyms "${vep_params.VEP_SYNONYMS}" \\
			--fasta "${vep_params.VEP_FASTA}" \\
			-custom "${vep_params.GNOMAD_GENOMES},gnomADg,vcf,exact,0,AF" \\
			-custom "${vep_params.GNOMAD_MT},gnomAD_mt,vcf,exact,0,AF_hom,AF_het" \\
			--dir_cache "${vep_params.VEP_CACHE}" \\
			--force_overwrite \\
			--no_stats \\
			--fork "${task.cpus}" \\
			--format vcf
		filter_indels.pl "${meta.sample}.only_indels.vep.vcf" > "${meta.sample}.only_indels.vep.filtered.vcf"
		${indel_vep_version(task)}
		"""

	stub:
		"""
		touch "${meta.sample}.only_indels.vep.filtered.vcf"
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