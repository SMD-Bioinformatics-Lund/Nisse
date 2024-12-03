process PARSE_TOMTE_QC {
    input:
    tuple path(multiqc_general_stats), path(picard_rna_coverage)

    output:
    path('hello_out.txt')

    script:
    """
    bash hello.sh "placeholder" > "hello_out.txt"
    """

    stub:
    """
    touch hello_out.txt
    """
}
