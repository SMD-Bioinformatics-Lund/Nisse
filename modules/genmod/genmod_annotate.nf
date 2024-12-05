process GENMOD_ANNOTATE {

    tag "${meta.sample}"
	label "process_low"
	container "${params.containers.genmod}"

    input:
        tuple val(meta), path(vcf)

    output:
        tuple val(meta), path("*_annotate.vcf"), emit: vcf

    script:
    """
    genmod \\
        annotate \\
        --annotate_regions \\
        --genome-build 38 \\
        --outfile ${meta.sample}_annotate.vcf \\
        ${vcf}
    """

    stub:
    """
    touch ${meta.sample}_annotate.vcf
    """
}