process PREPARE_DROP {

    tag "${meta.sample}"
    label "process_low"
    container "${params.containers.base}"
    publishDir "${params.outdir}/drop/", mode: 'copy', pattern: "*.tsv", overwrite: true

    input:
    val(label)
    tuple val(meta), path(drop_full)
    path(hgnc_map)
    val(stat_col)
    val(stat_cutoff)

    output:
    path "${meta.sample}_drop_parsed.tsv", emit: drop

    script:
    def hgnc_symbol_col = "hgncSymbol"
    """
    prepare_drop.py \\
        --in_path "${drop_full}" \\
        --out_path "${meta.sample}_drop_parsed.tsv" \\
        --stat_col "${stat_col}" \\
        --stat_cutoff "${stat_cutoff}" \\
        --hgnc_symbol_col "${hgnc_symbol_col}" \\
        --hgnc_symbol_id_map "${hgnc_map}"
    """

    stub:
    """
    touch "${meta.sample}_${label}_drop.tsv"
    """
}
