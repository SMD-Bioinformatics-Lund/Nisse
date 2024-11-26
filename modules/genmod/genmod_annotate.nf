process GENMOD_ANNOTATE {
    input:
        tuple val(prefix), path(vcf)

    output:
        tuple val(prefix), path("*_annotate.vcf"), emit: vcf

    script:
    """
    genmod \\
        annotate \\
        --annotate_regions \\
        --genome-build 38 \\
        --outfile ${prefix}_annotate.vcf \\
        ${vcf}
    """

    stub:
    """
    touch ${prefix}_annotate.vcf
    """
}