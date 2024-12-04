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
    def in_path = ""
    def out_path = ""
    def fdr_col = "pValue" // pValue or padjust
    def fdr_cutoff = ""
    def hgnc_symbol_col = ""
    def sample_col = ""
    """
    prepare_drop.py \\
        --in_path ${drop_full} \\
        --out_path ${meta.sample}_drop_parsed.tsv \\
        --stat_col "pValue" \\
        --stat_cutoff 0.5 \\
        --hgnc_symbol_col "hgncSymbol" \
        --sample_col "sampleID"
    """

    stub:
    """
    touch "${meta.sample}_${label}_drop.tsv"
    """
}
