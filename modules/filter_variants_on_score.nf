process filter_variants_on_score {
    input:
    path scored_vcf, path scored_vcf_tbi
    val threshold

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