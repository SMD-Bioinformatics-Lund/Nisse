process PARSE_TOMTE_QC {

    tag "${meta.sample}"
	label "process_low"
	container "${params.containers.base}"

    input:
    tuple val(meta), path(multiqc_general_stats), path(picard_rna_coverage)

    output:
    path("${meta.sample}_out.json"), emit: json

    script:
    """
    bash hello.sh "placeholder" > "${meta.sample}_out.json"
    """

    stub:
    """
    touch "${meta.sample}_out.txt"
    """
}
