process PREPARE_DROP {

    tag "${meta.sample}"
    label "process_low"
    container "${params.containers.base}"

    input:
    val(label)
    tuple val(meta), path(drop_full)
    path(hgnc_map)
    val(stat_col)
    val(stat_cutoff)

    output:
    tuple val(meta), path("${meta.sample}_${label}_parsed.tsv"), emit: drop

    script:
def hgnc_symbol_col = "hgncSymbol"
    """
    prepare_drop.py \\
        --in_path "${drop_full}" \\
        --out_path "${meta.sample}_${label}_parsed.tsv" \\
        --sample_col "sampleID" \\
        --sample "${meta.sample}" \\
        --stat_col "${stat_col}" \\
        --stat_cutoff "${stat_cutoff}" \\
        --hgnc_symbol_col "${hgnc_symbol_col}" \\
        --hgnc_symbol_id_map "${hgnc_map}"
    """

    stub:
    """
    touch "${meta.sample}_${label}_parsed.tsv"
    """
}
