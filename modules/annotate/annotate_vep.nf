process ANNOTATE_VEP {

	tag "${meta.sample}"
	label 'process_high'
	container "${params.containers.vep}"

	input:
	tuple val(meta), path(vcf), path(vcf_tbi)
	val vep_params

	output:
	tuple val(meta), path("${meta.sample}.vep.vcf"), emit: vcf
	tuple val(meta), path("*versions.yml"), emit: versions

	script:
	def group = "${meta.sample}"
	"""
		vep \\
			-i "${vcf}" \\
			-o "${group}.vep.vcf" \\
			--offline \\
			--sift b --polyphen b --ccds --hgvs --symbol --numbers --domains --regulatory --canonical --protein --biotype --af --af_1kg --max_af --pubmed --uniprot --mane --tsl --appris --variant_class --gene_phenotype --mirna \\
			--merged \\
			--vcf \\
			--no_stats \\
			--synonyms "${vep_params.VEP_SYNONYMS}" \\
			--fork "${task.cpus}" \\
			--force_overwrite \\
			--fasta "${vep_params.VEP_FASTA}" \\
			--dir_cache "${vep_params.VEP_CACHE}" \\
			--dir_plugins "${vep_params.VEP_PLUGINS}" \\
			--distance "${vep_params.VEP_TRANSCRIPT_DISTANCE}" \\
			-cache \\
			--format vcf \\
			--plugin "CADD,${vep_params.CADD}" \\
			--plugin "LoFtool" \\
			--plugin "MaxEntScan,${vep_params.MAXENTSCAN},SWA,NCSS" \\
			--plugin "dbNSFP,${vep_params.DBNSFP},transcript_match=1,REVEL_score,REVEL_rankscore" \\
			-custom "${vep_params.GNOMAD_EXOMES},gnomADe,vcf,exact,0,AF_popmax,AF,popmax" \\
			-custom "${vep_params.GNOMAD_GENOMES},gnomADg,vcf,exact,0,AF_popmax,AF,popmax" \\
			-custom "${vep_params.GNOMAD_MT},gnomAD_mt,vcf,exact,0,AF_hom,AF_het" \\
			-custom "${vep_params.PHYLOP},phyloP100way,bigwig" \\
			-custom "${vep_params.PHASTCONS},phastCons,bigwig"

		${annotate_vep_version(task)}
		"""

	stub:
	"""
		touch "${meta.sample}.vep.vcf"
		${annotate_vep_version(task)}
		"""
}
def annotate_vep_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    vep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : // ; s/ .*\$//')
	END_VERSIONS
	"""
}
