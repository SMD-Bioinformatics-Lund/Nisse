process GENMOD_MODELS {

    cpus 40

    input:
        tuple val(meta), val(prefix), path(vcf)
        path(family)

    output:
        tuple val(meta), val(prefix), path("*_models.vcf"), emit: vcf

    script:
    """
    genmod \\
        models \\
        --whole_gene \\
        --processes ${task.cpus} \\
        --family_file ${family} \\
        --outfile ${prefix}_models.vcf \\
        ${vcf}
    """

    stub:
    """
    touch ${prefix}_models.vcf
    """
}