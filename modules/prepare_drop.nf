// What exactly is this part performing?

// HGNC ID map

process PREPARE_DROP {
    input:
    val(meta)
    path(drop_full)
    path(hgn_map)

    output:
    path 'hello_out.txt'

    script:
    """
    prepare_drop.sh ${hgn_map} > "hello_out.txt"
    """

    stub:
    """
    touch hello_out.txt
    """
}
