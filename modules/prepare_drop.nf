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
    path "${meta.sample}_${label}_drop.tsv"

    script:
    def in_path = ""
    def out_dir = ""
    def fdr_col = "pValue/padjust"
    def fdr_cutoff = ""
    def hgnc_symbol_col = ""
    def sample_col = ""
    """
    prepare_drop.sh ${hgnc_map} > "${meta.sample}_${label}_drop.tsv"
    """

    stub:
    """
    touch "${meta.sample}_${label}_drop.tsv"
    """
}
