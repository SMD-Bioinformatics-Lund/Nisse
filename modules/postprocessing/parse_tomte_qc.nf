process PARSE_TOMTE_QC {

    tag "${meta.sample}"
	label "process_low"
    // FIXME: Does this container work?
	container "${params.containers.base}"

    input:
    tuple val(meta), path(multiqc_general_stats), path(picard_rna_coverage)

    output:
    path("${meta.sample}_out.json"), emit: json

    script:
    """
    parse_tomte_qc.py \\
        --multiqc_general_stats "${multiqc_general_stats}" \\
        --sample_id "${meta.sample}" \\
        --picard_rna_coverage "${picard_rna_coverage}" \\
        --output_file "${meta.sample}_out.json"
    """

    stub:
    """
    touch "${meta.sample}_out.txt"
    """
}
