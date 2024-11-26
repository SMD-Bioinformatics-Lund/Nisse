process GENMOD_COMPOUND {

    cpus 40

    input:
        tuple val(prefix), path(vcf)

    output:
        tuple val(prefix), path("*_compound.vcf"), emit: vcf

    script:
    """
    genmod \\
        compound \\
        --processes ${task.cpus} \\
        --outfile ${prefix}_compound.vcf \\
        ${vcf}
    """

    stub:
    """
    touch ${prefix}_compound.vcf
    """
}