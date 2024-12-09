process BGZIP_TABIX {

    tag "${meta.sample}"
    label "process_low"
    container "${params.containers.bcftools}"

    input:
    tuple val(meta), path(junction_bed)

    output:
    tuple val(meta), path("${meta.sample}.bed.gz"), path("${meta.sample}.bed.gz.tbi"), emit: bed

    script:
    """
    bgzip "${junction_bed}"
    tabix "${junction_bed}.gz"
    """

    stub:
    """
    touch "${junction_bed}.bed.gz"
    touch "${junction_bed}.bed.gz.tbi"
    """
}
