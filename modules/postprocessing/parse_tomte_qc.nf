process PARSE_TOMTE_QC {

    tag "${meta.sample}"
	label "process_low"
	container "${params.containers.base}"

    input:
    tuple val(meta), path(multiqc_general_stats), path(picard_rna_coverage), path(hb_estimates)

    output:
    path("${meta.sample}_out.json"), emit: json
    path "versions.yml", emit: versions

    script:
    """

    parse_tomte_qc.py \\
        --multiqc_general_stats "${multiqc_general_stats}" \\
        --sample_id "${meta.sample}" \\
        --picard_rna_coverage "${picard_rna_coverage}" \\
        --hb_estimate "${hb_estimates}" \\
        --output_file "${meta.sample}_out.json"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parse_tomte_qc: \$(parse_tomte_qc.py --version)
    END_VERSIONS
    """

    stub:
    """
    touch "${meta.sample}_out.json"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parse_tomte_qc: \$(parse_tomte_qc.py --version)
    END_VERSIONS
    """
}
