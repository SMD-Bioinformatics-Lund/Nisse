process PREPARE_DROP {

    tag "${meta.sample}"
	label "process_small"
	container "${params.containers.base}"

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
