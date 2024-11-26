process parse_tomte_qc {
    input:
    path multiqc_general_stats, path picard_rna_coverage

    output:
    path 'hello_out.txt'

    script:
    """
    bash hello.sh ${input_file} > "hello_out.txt"
    """

    stub:
    """
    touch hello_out.txt
    """
}
