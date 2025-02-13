def bcftools_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    bcftools: \$(echo \$(bcftools --version 2>&1) | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
	END_VERSIONS
	"""
}

// FIXME: publishDir
process IDSNP_CALL {
    label 'process_single'
    tag "${meta.id}"

    // Let's start with all beds
    input:
        tuple val(meta), path(bam), path(bai)
        val idsnp_params
        // path idsnp_bed
        // path idSnp_bed_gz
        // path idSnp_std_bed_gz
        // path header
        // path genome

    output:
        tuple val(meta), path("*final.vcf"), path("*genotypes.json"),       emit:   sample_id_genotypes
        path "versions.yml",                                                            emit:   versions

    script:
        def prefix  = "${meta.sample}"
        // Document flags
        // -d 1000
        // -q 10
        // FIXME: Check the params.bed arguments
        // FIXME: Is there also a variant argument needed?
        """
        bcftools mpileup \\
            -Ou \\
            -R "${idsnp_params.idsnp_bed}" \\
            -f "${idsnp_params.genome}" \\
            -d 1000 \\
            -q 10 \\
            "${bam}" | \\
        bcftools call \\
            -A \\
            -C alleles \\
            -T "${idsnp_params.idSnp_bed_gz}" \\
            -m \\
            -Ov \\
            > "${prefix}.raw.vcf"
        
        bcftools annotate \\
            -a "${idsnp_params.idSnp_std_bed_gz}" \\
            -c "CHROM,FROM,TO,ID" \\
            -h "${idsnp_params.header}" \\
            -o "${meta.sample}.final.vcf" "${prefix}.raw.vcf"
        """

    stub:
        def prefix = "${meta.sample}"
        """
        touch ${prefix}.final.vcf
        touch ${prefix}.genotypes.json

        ${bcftools_version(task)}
        """
}

process PERC_HETEROZYGOTES {
    tag "${meta.sample}"
    label "process_low"
    container "${params.containers.bcftools}"

    input:
        // What is the targets_bed vs targets_tsv?
        tuple val(meta), path(bam), path(targets_bed), path(targets_tsv), path(genome)
    
    output:
        tuple val(meta)

    script:
    """
    bcftools mpileup \\
        -Ou \\
        -R "${targets_bed}" \\
        -f genome.fa \\
        input.bam | \\
    bcftools call \\
        -A \\
        -C alleles \\
        -T targets-file.tsv.gz \\
        -m -Ov - \\
        > "${meta.sample}_heterozygosity_calls.vcf"
    """

    stub:
    """
    touch "${meta.sample}_heterozygosity_calls.vcf"
    """
}

