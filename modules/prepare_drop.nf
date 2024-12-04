process PREPARE_DROP {

    tag "${meta.sample}"
    label "process_low"
    container "${params.containers.base}"
    publishDir "${params.outdir}/drop/", mode: 'copy', pattern: "*.tsv", overwrite: true

    input:
    val(label)
    tuple val(meta), path(drop_full)
    path(hgnc_map)

    output:
    path "${meta.sample}_drop_parsed.tsv"

    script:
    // FIXME: What to use here? What did we use during validation?
    def stat_col = "pValue" // pValue or padjust
    def stat_cutoff = 0.5
    def hgnc_symbol_col = "hgncSymbol"
    """
    prepare_drop.py \\
        --in_path "${drop_full}" \\
        --out_path "${meta.sample}_drop_parsed.tsv" \\
        --stat_col "${stat_col}" \\
        --stat_cutoff "${stat_cutoff}" \\
        --hgnc_symbol_col "${hgnc_symbol_col}"
    """

    stub:
    """
    touch "${meta.sample}_${label}_drop.tsv"
    """
}
