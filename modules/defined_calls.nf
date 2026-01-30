def bcftools_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    bcftools: \$(echo \$(bcftools --version 2>&1) | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
	END_VERSIONS
	"""
}

process IDSNP_CALL {
    label 'process_single'
    tag "${meta.id}"
    container "${params.containers.bcftools}"

    input:
        tuple val(meta), path(bam), path(bai)
        val idsnp_params

    output:
        tuple val(meta), path("*final.vcf"), emit: vcf
        path("*versions.yml"), emit: versions

    script:
        def prefix  = "${meta.sample}"
        def max_depth = 1000
        def min_map_qual = 10
        """
        bcftools mpileup \\
            -Ou \\
            -R "${idsnp_params.idsnp_bed}" \\
            -f "${idsnp_params.genome}" \\
            -d "${max_depth}" \\
            -q "${min_map_qual}" \\
            "${bam}" | \\
        bcftools call \\
            -A \\
            -C alleles \\
            -T "${idsnp_params.idsnp_tsv}" \\
            -m \\
            -Ov \\
            > "${prefix}.raw.vcf"
        
        bcftools annotate \\
            -a "${idsnp_params.variant_rsids_bed}" \\
            -c "CHROM,FROM,TO,ID" \\
            -h "${idsnp_params.header}" \\
            -o "${prefix}.final.vcf" "${prefix}.raw.vcf"

        ${bcftools_version(task)}
        """

    stub:
        def prefix = "${meta.sample}"
        """
        touch ${prefix}.final.vcf
        touch ${prefix}.genotypes.json

        ${bcftools_version(task)}
        """
}

process IDSNP_VCF_TO_JSON {
    label 'process_single'
    tag "${meta.id}"
    container "${params.containers.python}"

    input:
        tuple val(meta), path(vcf)
    
    output:
        tuple val(meta), path("*.json"), emit: json
        path("*versions.yml"), emit: versions
    
    script:
    def prefix = "${meta.sample}"
    """
    genotype_to_json.py "${vcf}" "${prefix}.genotypes.json"

    ${python_version(task)}
    """

    stub:
    def prefix = "${meta.sample}"
    """
    touch "${prefix}.genotypes.json"

    ${python_version(task)}
    """
}
def python_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    python: \$(echo \$(python --version | cut -f2 -d" "))
	END_VERSIONS
	"""
}

process PERC_HETEROZYGOTES {
    tag "${meta.sample}"
    label "process_low"
    container "${params.containers.bcftools}"

    input:
        tuple val(meta), path(bam), path(bai)
        val hetcall_params
    
    output:
        tuple val(meta), path("*_heterozygosity_calls.vcf"), emit: vcf
        path("*versions.yml"), emit: versions

    script:
    def prefix = "${meta.sample}"
    """
    bcftools mpileup \\
        -Ou \\
        -R "${hetcall_params.targets_bed}" \\
        -f "${hetcall_params.genome}" \\
        "${bam}" | \\
    bcftools call \\
        -A \\
        -C alleles \\
        -T "${hetcall_params.targets_tsv}" \\
        -m -Ov - \\
        > "${prefix}.raw.vcf"
    
    bcftools annotate \\
        -a "${hetcall_params.variant_rsids_bed}" \\
        -c "CHROM,FROM,TO,ID" \\
        -h "${hetcall_params.header}" \\
        -o "${prefix}_heterozygosity_calls.vcf" "${prefix}.raw.vcf"
    
    ${bcftools_version(task)}
    """

    stub:
    def prefix = "${meta.sample}"
    """
    touch "${prefix}_heterozygosity_calls.vcf"

    ${bcftools_version(task)}
    """
}

