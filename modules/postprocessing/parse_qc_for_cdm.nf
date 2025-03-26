process PARSE_QC_FOR_CDM {

    tag "${meta.sample}"
	label "process_low"
	container "${params.containers.base}"

    input:
    tuple val(meta), path(multiqc_general_stats), path(multiqc_star_summary), path(multiqc_picard_rna_coverage), path(hb_estimates), path(hetcalls_vcf)

    output:
    path("${meta.sample}_out.json"), emit: json
    tuple val(meta), path("versions.yml"), emit: versions

    script:
    """
    parse_qc_for_cdm.py \\
        --multiqc_general_stats "${multiqc_general_stats}" \\
        --multiqc_star "${multiqc_star_summary}" \\
        --sample_id "${meta.sample}" \\
        --picard_rna_coverage "${multiqc_picard_rna_coverage}" \\
        --hb_estimate "${hb_estimates}" \\
        --hetcalls_vcf "${hetcalls_vcf}" \\
        --output_file "${meta.sample}_out.json"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parse_qc_for_cdm: \$(parse_qc_for_cdm.py --version)
    END_VERSIONS
    """

    stub:
    """
    touch "${meta.sample}_out.json"

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        parse_qc_for_cdm: \$(parse_qc_for_cdm.py --version)
    END_VERSIONS
    """
}
