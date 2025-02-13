// NEXT: What are the input / outputs?
// Start drafting the processes

// /fs1/paul/idsnp/README
// /fs1/jakob/src/util/data_utils/idsnp/idsnp.sh
// bcftools only? And some script?
// Check the somatic pipeline as well

def bcftools_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    bcftools: \$(echo \$(bcftools --version 2>&1) | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
	END_VERSIONS
	"""
}

process ALLELE_CALL {
    label 'process_single'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(bam), file(bai)

    output:
        tuple val(group), val(meta), file("*final.vcf"), file("*genotypes.json"),       emit:   sample_id_genotypes
        path "versions.yml",                                                            emit:   versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix  = task.ext.prefix ?: "${meta.id}"
        def args    = task.ext.args  ?: ""
        def args2   = task.ext.args2 ?: ""
        def args3   = task.ext.args3 ?: ""
        """
        bcftools mpileup $args $bam | bcftools call $args2   > ${prefix}.raw.vcf
        bcftools annotate $args3 -o ${prefix}.final.vcf ${prefix}.raw.vcf
        bcftools query -f '%ID\\t[%GT]\\n' ${prefix}.final.vcf > ${prefix}.genotypes
        genotype2json.py ${prefix}.genotypes ${prefix}.genotypes.json

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$(echo \$(bcftools --version 2>&1) | sed 's/bcftools //; s/ .*//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        touch ${prefix}.final.vcf
        touch ${prefix}.genotypes.json

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$(echo \$(bcftools --version 2>&1) | sed 's/bcftools //; s/ .*//')
        END_VERSIONS
        """
}

process IDSNPS {
    tag "${meta.sample}"
    label "process_low"
    container "${params.containers.bcftools}"

    input:
        tuple val(meta), path(bam), path(bai), path(fasta), path(id_bed), path(id_tsv)
    
    output:
        val(meta), path("${meta.sample}_idsnps.vcf")
    
    script:
    """
    bcftools mpileup -Ou -R "${id_bed}" -f "${fasta}" "${bam}" | \\
        bcftools call -A -C alleles -T "${id_tsv}" -m -Ov - > "${meta.sample}_idsnps.vcf"

    ${bcftools_version(task)}
    """

    stub:
    """
    touch "${meta.sample}_idsnps.vcf"

    ${bcftools_version(task)}
    """
}



// BCFTOOLS here as Paul did it?
// Notes from Paul:

// its all in the `/fs2/paul/rnaseq_scanb_genotyping` . The README is describing how to get calls from a target-file using bcftools. The files used are the following:  
// `genome-fa` -> `/fs1/resources/ref/hg38/fasta/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna`  
// `target.bed` -> `/fs2/paul/rnaseq_scanb_genotyping/genotyping-213-snp_feb2018_50nt-flanks.bed`  
// `targets-file.tsv.gz`  -> `/fs2/paul/rnaseq_scanb_genotyping/genotyping-213-snp_feb2018.chr.tsv.gz`  
// In the bam-folder I have collected bam-files run on tomte without downsampling.
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

