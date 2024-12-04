process GENMOD_MODELS {

    tag "${meta.sample}"
	label "process_small"
	container "${params.containers.base}"

    input:
        tuple val(meta), path(vcf)
        tuple val(meta2), path(family)

    output:
        tuple val(meta), path("${meta.sample}_models.vcf"), emit: vcf

    script:
    """
    genmod \\
        models \\
        --whole_gene \\
        --processes ${task.cpus} \\
        --family_file ${family} \\
        --outfile ${meta.sample}_models.vcf \\
        ${vcf}
    """

    stub:
    """
    touch ${meta.sample}_models.vcf
    """
}