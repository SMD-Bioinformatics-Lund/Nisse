process ANNOTATE_VEP {

	tag "${meta.sample}"
    label 'process_large'
	container "${params.containers.vep}"

	input:
		tuple val(meta), file(vcf)

	output:
		tuple val(meta), file("${meta.sample}.vep.vcf"), emit: vcf
		tuple val(meta), file("*versions.yml"), emit: versions

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
			--synonyms "${params.vep.VEP_SYNONYMS}" \\
			--fork "${task.cpus}" \\
			--force_overwrite \\
			--fasta "${params.vep.VEP_FASTA}" \\
			--dir_cache "${params.vep.VEP_CACHE}" \\
			--dir_plugins "${params.vep.VEP_PLUGINS}" \\
			--distance "${params.vep.VEP_TRANSCRIPT_DISTANCE}" \\
			-cache \\
			--plugin "CADD,${params.vep.CADD}" \\
			--plugin "LoFtool" \\
			--plugin "MaxEntScan,${params.vep.MAXENTSCAN},SWA,NCSS" \\
			--plugin "dbNSFP,${params.vep.DBNSFP},transcript_match=1,REVEL_score,REVEL_rankscore" \\
			-custom "${params.vep.GNOMAD_EXOMES},gnomADe,vcf,exact,0,AF_popmax,AF,popmax" \\
			-custom "${params.vep.GNOMAD_GENOMES},gnomADg,vcf,exact,0,AF_popmax,AF,popmax" \\
			-custom "${params.vep.GNOMAD_MT},gnomAD_mt,vcf,exact,0,AF_hom,AF_het" \\
			-custom "${params.vep.PHYLOP},phyloP100way,bigwig" \\
			-custom "${params.vep.PHASTCONS},phastCons,bigwig"

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