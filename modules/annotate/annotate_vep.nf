// FIXME: Can I pass in the params as a value blob? Instead of doing global access

process ANNOTATE_VEP {
	container = "${params.container_vep}"
	cpus 30
	tag "${meta.id}"
	memory '50 GB'
	time '5h'

	input:
		set meta, file(vcf)

	output:
		set meta, file("${group}.vep.vcf")
		set meta, file("*versions.yml")

	script:
		"""
		vep \\
			-i ${vcf} \\
			-o ${group}.vep.vcf \\
			--offline \\
			--sift b --polyphen b --ccds --hgvs --symbol --numbers --domains --regulatory --canonical --protein --biotype --af --af_1kg --max_af --pubmed --uniprot --mane --tsl --appris --variant_class --gene_phenotype --mirna \\
			--merged \\
			--vcf \\
			--no_stats \\
			--synonyms $params.VEP_SYNONYMS \\
			--fork ${task.cpus} \\
			--force_overwrite \\
			--fasta $params.VEP_FASTA \\
			--dir_cache $params.VEP_CACHE \\
			--dir_plugins $params.VEP_PLUGINS \\
			--distance $params.VEP_TRANSCRIPT_DISTANCE \\
			-cache \\
			--plugin CADD,$params.CADD \\
			--plugin LoFtool \\
			--plugin MaxEntScan,$params.MAXENTSCAN,SWA,NCSS \\
			--plugin dbNSFP,$params.DBNSFP,transcript_match=1,REVEL_score,REVEL_rankscore \\
			-custom $params.GNOMAD_EXOMES,gnomADe,vcf,exact,0,AF_popmax,AF,popmax \\
			-custom $params.GNOMAD_GENOMES,gnomADg,vcf,exact,0,AF_popmax,AF,popmax \\
			-custom $params.GNOMAD_MT,gnomAD_mt,vcf,exact,0,AF_hom,AF_het \\
			-custom $params.PHYLOP,phyloP100way,bigwig \\
			-custom $params.PHASTCONS,phastCons,bigwig

		${annotate_vep_version(task)}
		"""

	stub:
		"""
		touch "${group}.vep.vcf"
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