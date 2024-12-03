// What exactly is this part performing?

// HGNC ID map
// Split 

process PREPARE_DROP {
    input:
        val(meta)
        val(label)
        path(drop_full)
        path(hgn_map)

    output:
        path "${meta.sample}_${label}_drop.tsv"

    script:
    """
    prepare_drop.sh ${hgn_map} > "${meta.sample}_${label}_drop.tsv"
    """

    stub:
    """
    touch "${meta.sample}_${label}_drop.tsv"
    """
}
